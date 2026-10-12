/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstdlib>
#include <cstring>

using namespace std;
using namespace kraken;

string Input_DB_filename, Output_DB_filename, Index_filename;
string Tmp_dir;  // from -T, else $TMPDIR, else /tmp (resolved in parse_command_line)
double Mem_budget_gib = 1.0;  // -m in GiB; default 1.0 GiB
uint8_t Bin_key_nt = 15;
int Num_threads = 1;
bool Zero_vals = false;
bool Operate_in_RAM = false;
// Global until I can find a way to pass this to the sorting function
size_t Key_len = 8;

static int pair_cmp(const void *a, const void *b);
static void parse_command_line(int argc, char **argv);
static void bin_and_sort_data_ram(KrakenDB &kdb, char *data, KrakenDBIndex &idx);
static void bin_and_sort_data_disk(KrakenDB &kdb, KrakenDBIndex &idx,
                                   ofstream &output_file);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  cerr << "db_sort: Getting database into memory ...";
  QuickFile input_db_file(Input_DB_filename);
  KrakenDB *input_db = new KrakenDB(input_db_file.ptr());
  Key_len = input_db->get_key_len();
  uint64_t val_len = input_db->get_val_len();
  uint64_t key_ct = input_db->get_key_ct();
  input_db->make_index(Index_filename, Bin_key_nt);
  QuickFile index_file(Index_filename);
  KrakenDBIndex db_index(index_file.ptr());

  uint64_t skip_len = input_db->header_size();
  char *header = new char[ skip_len ];
  memcpy(header, input_db_file.ptr(), skip_len);

  delete input_db;
  // No longer legit for traversal, but allows access to KDB functions
  input_db = new KrakenDB(header);
  input_db_file.close_file();  // Stop using memory-mapped file

  cerr << "db_sort: Sorting ...";
  if (Operate_in_RAM) {
    // original behavior: whole database loaded in memory
    char *data = new char[ key_ct * (Key_len + val_len) ];
    // Populate data w/ pairs from DB and sort bins in parallel
    bin_and_sort_data_ram(*input_db, data, db_index);

    cerr << "db_sort: Sorting complete - writing database to disk ..." << endl;
    ofstream output_file(Output_DB_filename.c_str(), std::ofstream::binary);
    if (!output_file.is_open())
      errx(EX_CANTCREAT, "can't open output file %s", Output_DB_filename.c_str());
    output_file.write(header, skip_len);
    if (!output_file)
      err(EX_IOERR, "error writing output file %s", Output_DB_filename.c_str());
    output_file.write(data, key_ct * (Key_len + val_len));
    output_file.close();
    if (!output_file)
      err(EX_IOERR, "error writing output file %s", Output_DB_filename.c_str());
    delete [] data;
  }
  else {
    // disk-backed: bins are staged in a temp file and sorted one chunk at a time
    ofstream output_file(Output_DB_filename.c_str(), std::ofstream::binary);
    if (!output_file.is_open())
      errx(EX_CANTCREAT, "can't open output file %s", Output_DB_filename.c_str());
    output_file.write(header, skip_len);
    if (!output_file)
      err(EX_IOERR, "error writing output file %s", Output_DB_filename.c_str());
    bin_and_sort_data_disk(*input_db, db_index, output_file);
    output_file.close();
    if (!output_file)
      errx(EX_IOERR, "error writing output file %s", Output_DB_filename.c_str());
    cerr << "db_sort: Sorting complete - database written to disk" << endl;
  }

  return 0;
}

static void bin_and_sort_data_ram(KrakenDB &kdb, char *data, KrakenDBIndex &idx) {
  uint8_t nt = idx.indexed_nt();
  uint64_t *offsets = idx.get_array();
  uint64_t entries = 1ull << (nt * 2);
  uint64_t key_len = kdb.get_key_len();
  uint64_t val_len = kdb.get_val_len();
  uint64_t pair_size = key_len + val_len;
  // Use VLA-like small buffer for reads; pair_size is small and known
  vector<char> pair(pair_size);

  ifstream input_file(Input_DB_filename.c_str(), std::ifstream::binary);
  input_file.seekg(kdb.header_size(), ios_base::beg);

  // Create a copy of the offsets array for use as insertion positions
  vector<uint64_t> pos(offsets, offsets + entries);
  for (uint64_t i = 0; i < kdb.get_key_ct(); i++) {
    input_file.read(&pair[0], pair_size);
    if ((uint64_t) input_file.gcount() != pair_size)
      errx(EX_DATAERR, "unexpected end of input database");
    uint64_t kmer = 0;
    memcpy(&kmer, &pair[0], key_len);
    uint64_t b_key = kdb.bin_key(kmer, nt);
    char *pair_pos = data + pair_size * pos[b_key]++;
    // Copy pair into correct bin (but not final position)
    memcpy(pair_pos, &pair[0], pair_size);
    if (Zero_vals)
      memset(pair_pos + key_len, 0, val_len);
  }
  input_file.close();

  // Sort all bins
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (uint64_t i = 0; i < entries; i++) {
    qsort(data + offsets[i] * pair_size,
          offsets[i+1] - offsets[i], pair_size,
          pair_cmp);
  }
}

// Write exactly len bytes at offset, retrying on short writes
static void full_pwrite(int fd, const char *buf, size_t len, off_t off) {
  while (len > 0) {
    ssize_t n = pwrite(fd, buf, len, off);
    if (n < 0) {
      if (errno == EINTR) continue;
      err(EX_IOERR, "write to temporary file failed");
    }
    if (n == 0)
      errx(EX_IOERR, "short write to temporary file");
    buf += n;
    off += n;
    len -= (size_t) n;
  }
}

// Read exactly len bytes at offset, retrying on short reads
static void full_pread(int fd, char *buf, size_t len, off_t off) {
  while (len > 0) {
    ssize_t n = pread(fd, buf, len, off);
    if (n < 0) {
      if (errno == EINTR) continue;
      err(EX_IOERR, "read from temporary file failed");
    }
    if (n == 0)
      errx(EX_IOERR, "unexpected end of temporary file");
    buf += n;
    off += n;
    len -= (size_t) n;
  }
}

// Disk-backed implementation: only a chunk of bins is held in memory at once.
// Final output is identical to the in-memory path.
static void bin_and_sort_data_disk(KrakenDB &kdb, KrakenDBIndex &idx,
                                   ofstream &output_file) {
  uint8_t nt = idx.indexed_nt();
  uint64_t *offsets = idx.get_array();
  uint64_t entries = 1ull << (nt * 2);
  uint64_t key_len = kdb.get_key_len();
  uint64_t val_len = kdb.get_val_len();
  uint64_t pair_size = key_len + val_len;
  uint64_t key_ct = kdb.get_key_ct();
  vector<char> pair(pair_size);

  // Convert GiB budget to bytes (GiB == 2^30 bytes)
  if (!(Mem_budget_gib > 0.0))
    errx(EX_USAGE, "memory budget must be a positive number of GiB");
  const long double max_gib = (long double) ( (long double) UINT64_MAX / (long double)(1ULL << 30) );
  if ((long double) Mem_budget_gib > max_gib)
    errx(EX_USAGE, "memory budget too large");
  const uint64_t budget = (uint64_t) (Mem_budget_gib * (1ULL << 30));

  // Temp file, unlinked immediately so it is cleaned up on any exit
  string tmpl = Tmp_dir + "/db_sort_XXXXXX";
  vector<char> tmpl_buf(tmpl.begin(), tmpl.end());
  tmpl_buf.push_back('\0');
  int tmpfd = mkstemp(&tmpl_buf[0]);
  if (tmpfd < 0)
    err(EX_CANTCREAT, "can't create temporary file in %s", Tmp_dir.c_str());
  unlink(&tmpl_buf[0]);

  if (ftruncate(tmpfd, (off_t) (key_ct * pair_size)) != 0)
    err(EX_IOERR, "can't size temporary file in %s (disk full?)", Tmp_dir.c_str());

  ifstream input_file(Input_DB_filename.c_str(), std::ifstream::binary);
  input_file.seekg(kdb.header_size(), ios_base::beg);

  // Copy of offsets used as insertion positions
  vector<uint64_t> pos(offsets, offsets + entries);
  for (uint64_t i = 0; i < key_ct; i++) {
    input_file.read(&pair[0], pair_size);
    if ((uint64_t) input_file.gcount() != pair_size)
      errx(EX_DATAERR, "unexpected end of input database");
    uint64_t kmer = 0;
    memcpy(&kmer, &pair[0], key_len);
    uint64_t b_key = kdb.bin_key(kmer, nt);
    if (Zero_vals)
      memset(&pair[0] + key_len, 0, val_len);
    full_pwrite(tmpfd, &pair[0], pair_size, (off_t) (pos[b_key]++ * pair_size));
  }
  input_file.close();

  // Group consecutive bins into chunks of at most `budget` bytes. Bins are
  // contiguous in the temp file, so each chunk is one sequential read. Bins
  // within a chunk are sorted in parallel; chunks are written in order.
  // A single bin larger than the budget forms its own chunk.
  char *buf = NULL;
  size_t cap = 0;
  uint64_t i = 0;
  while (i < entries) {
    uint64_t first_off = offsets[i];
    uint64_t j = i;
    // grow j while total bytes from first_off to offsets[j+1] fit in budget
    while (j < entries && (offsets[j+1] - first_off) * pair_size <= budget)
      j++;
    if (j == i)  // bin i alone exceeds the budget
      j = i + 1;
    uint64_t bytes64 = (offsets[j] - first_off) * pair_size;
    if (bytes64 == 0) {  // only empty bins
      i = j;
      continue;
    }
    size_t bytes = (size_t) bytes64;
    if (bytes > cap) {
      free(buf);
      buf = (char *) malloc(bytes);
      if (buf == NULL)
        errx(EX_OSERR, "can't allocate %llu bytes for bin(s) starting at %llu; "
             "try a larger -n value or a larger -m",
             (unsigned long long) bytes, (unsigned long long) i);
      cap = bytes;
    }
    full_pread(tmpfd, buf, bytes, (off_t) (first_off * pair_size));

    // Parallel sort of the bins in this chunk
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (uint64_t b = i; b < j; b++) {
      qsort(buf + (offsets[b] - first_off) * pair_size,
            offsets[b+1] - offsets[b], pair_size, pair_cmp);
    }

    output_file.write(buf, bytes);
    if (!output_file)
      err(EX_IOERR, "error writing output file %s", Output_DB_filename.c_str());
    i = j;
  }
  free(buf);
  if (close(tmpfd) != 0)
    err(EX_IOERR, "error closing temporary file");
}

static int pair_cmp(const void *a, const void *b) {
  uint64_t aval = 0, bval = 0;
  memcpy(&aval, a, Key_len);
  memcpy(&bval, b, Key_len);
  if (aval < bval)
    return -1;
  else if (aval == bval)
    return 0;
  else
    return 1;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "n:d:o:i:t:T:m:zM")) != -1) {
    switch (opt) {
      case 'n' :
        sig = atoll(optarg);
        if (sig < 1 || sig > 31)
          errx(EX_USAGE, "bin key length out of range");
        Bin_key_nt = (uint8_t) sig;
        break;
      case 'd' :
        Input_DB_filename = optarg;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      case 'T' :
        Tmp_dir = optarg;
        break;
      case 'm' :
      {
        char *endptr = NULL;
        double v = strtod(optarg, &endptr);
        if (endptr == optarg || v <= 0.0)
          errx(EX_USAGE, "memory budget must be a positive number (GiB)");
        Mem_budget_gib = v;
        break;
      }
      case 't' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive thread count");
        #ifdef _OPENMP
        if (sig > omp_get_num_procs())
          errx(EX_USAGE, "thread count exceeds number of processors");
        Num_threads = sig;
        omp_set_num_threads(Num_threads);
        #endif
        break;
      case 'z' :
        Zero_vals = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (Input_DB_filename.empty() || Output_DB_filename.empty()
      || Index_filename.empty())
    usage();

  // Resolve temp directory: -T, then $TMPDIR, then /tmp
  if (Tmp_dir.empty()) {
    const char *env = getenv("TMPDIR");
    Tmp_dir = (env != NULL && env[0] != '\0') ? env : "/tmp";
  }
  while (Tmp_dir.size() > 1 && Tmp_dir[Tmp_dir.size() - 1] == '/')
    Tmp_dir.erase(Tmp_dir.size() - 1);
  if (!Operate_in_RAM) {
    struct stat st;
    if (stat(Tmp_dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode))
      errx(EX_USAGE, "temporary directory %s does not exist or is not a directory",
           Tmp_dir.c_str());
  }
}

void usage(int exit_code) {
  cerr << "Usage: db_sort [-z] [-M] [-t threads] [-T tmpdir] [-m GiB] <-d input db> <-o output db> <-i output idx>\n"
       << "  -M  sort entirely in RAM (original behavior); otherwise sort on disk\n"
       << "  -T  directory for temporary file (default: $TMPDIR, else /tmp)\n"
       << "  -m  memory budget in GiB for disk-mode sorting (default: 1.0)\n";
  exit(exit_code);
}
