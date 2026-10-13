/*
 * Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
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
#include "krakenutil.hpp"
#include "seqreader.hpp"
#include "taxdb.hpp"
#include "uid_mapping.hpp"
#include <unordered_map>
#include <map>
#include <algorithm>
#include <cmath>
#include <sys/file.h>

#define SKIP_LEN 50000

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_files();
void process_single_file();
void process_file(string filename, uint32_t taxid);
void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish, bool is_contaminant_taxid = false);
void add_sequence(uint32_t taxid, string &seq, bool is_contaminant_taxid);
void flush_kmers();

int Num_threads = 1;
string DB_filename, Index_filename,
  Output_DB_filename, TaxDB_filename,
  Kmer_count_filename,
  File_to_taxon_map_filename,
  ID_to_taxon_map_filename, Multi_fasta_filename;
bool Force_contaminant_taxid = false;
bool Reset_taxid = false;
uint32_t New_taxid_start = 1000000000;

bool Allow_extra_kmers = false;
bool verbose = false;
bool Operate_in_RAM = false;
bool One_FASTA_file = false;
bool Add_taxIds_for_Assembly = false;
bool Add_taxIds_for_Sequences = false;
bool Use_uids_instead_of_taxids = false;
bool Output_UID_map_to_STDOUT = false;
bool Pretend = false;
uint32_t Minimum_sequence_size = 0;

// Memory budget (GiB) for the batched k-mer update buffer
double Memory_budget_gib = 1.0;

string UID_map_filename;
ofstream UID_map_file;

uint32_t current_uid = 0;
unordered_map<uint32_t, uint32_t> Parent_map;
vector< const TaxidSet*  > UID_to_taxids_vec;
map< TaxidSet, uint32_t> Taxids_to_UID_map;

unordered_map<string, uint32_t> ID_to_taxon_map;
unordered_map<uint32_t, bool> SeqId_added;
KrakenDB Database;
TaxonomyDB<uint32_t> taxdb;

const string prefix = "kraken:taxid|";

// do not add sequence taxIDs for host sequences (currently only human and mouse)
const uint32_t TID_HUMAN = 9606;
const uint32_t TID_MOUSE = 10090;

// k-mers appearing in contaminant sequences will keep the contaminant
//  sequence taxid, even if they also appear in a genome
const uint32_t TID_CONTAMINANT1 = 32630; // 'synthetic construct'
const uint32_t TID_CONTAMINANT2 = 81077; // 'artificial sequences'

// ---------------------------------------------------------------------
// Batched update machinery (non-UID mode)
// ---------------------------------------------------------------------
struct KmerRec {
  uint64_t kmer;   // canonical k-mer
  uint32_t bin;    // DB index bin of the k-mer
  uint32_t seq;    // index into Seq_info (preserves original update order)
};

struct SeqInfo {
  uint32_t taxid;
  bool is_contaminant;
};

static inline bool kmer_rec_less(const KmerRec &a, const KmerRec &b) {
  if (a.bin != b.bin) return a.bin < b.bin;
  if (a.kmer != b.kmer) return a.kmer < b.kmer;
  return a.seq < b.seq;
}

vector<KmerRec> Kmer_buffer;
vector<SeqInfo> Seq_info;
size_t Max_records = 0;
size_t Group_chunks = 1;

// Cache of lca() results (lca is a pure function of its arguments)
unordered_map<uint64_t, uint32_t> Lca_cache;
const size_t LCA_CACHE_MAX = 4u << 20;

static inline uint32_t cached_lca(uint32_t a, uint32_t b) {
  const uint64_t key = ((uint64_t)a << 32) | b;
  auto it = Lca_cache.find(key);
  if (it != Lca_cache.end())
    return it->second;
  uint32_t r = lca(Parent_map, a, b);
  if (Lca_cache.size() >= LCA_CACHE_MAX)
    Lca_cache.clear();
  Lca_cache.emplace(key, r);
  return r;
}

// Apply one update to a k-mer value (same semantics as the original code)
static inline void apply_update(uint32_t &val, uint32_t taxid, bool is_contaminant_taxid) {
  if (Reset_taxid) {
    val = 0;
  } else if (!Force_contaminant_taxid) {
    val = cached_lca(taxid, val);
  } else {
    if (val == TID_CONTAMINANT1 || val == TID_CONTAMINANT2) {
      // keep value
    } else if (is_contaminant_taxid) {
      // assign the taxid of the (last) sequence to k-mers
      val = taxid;
    } else {
      val = cached_lca(taxid, val);
    }
  }
}

void copy_database_file(const string &source, const string &destination) {
  struct stat source_stat, destination_stat;
  if (stat(source.c_str(), &source_stat) < 0)
    err(EX_NOINPUT, "can't stat database %s", source.c_str());
  if (stat(destination.c_str(), &destination_stat) == 0) {
    if (source_stat.st_dev == destination_stat.st_dev &&
        source_stat.st_ino == destination_stat.st_ino)
      return;
  } else if (errno != ENOENT) {
    err(EX_CANTCREAT, "can't stat output database %s", destination.c_str());
  }

  ifstream input(source.c_str(), ios::in | ios::binary);
  if (!input.is_open())
    err(EX_NOINPUT, "can't open database %s", source.c_str());
  ofstream output(destination.c_str(), ios::out | ios::binary | ios::trunc);
  if (!output.is_open())
    err(EX_CANTCREAT, "can't open output database %s", destination.c_str());

  char buffer[1024 * 1024];
  while (input) {
    input.read(buffer, sizeof(buffer));
    streamsize bytes_read = input.gcount();
    if (bytes_read > 0)
      output.write(buffer, bytes_read);
    if (!output)
      err(EX_IOERR, "error writing output database %s", destination.c_str());
  }
  if (!input.eof())
    err(EX_IOERR, "error reading database %s", source.c_str());
  output.close();
  if (!output)
    err(EX_IOERR, "error closing output database %s", destination.c_str());
}

void validate_database_layout(KrakenDB &database, size_t db_file_size,
                              KrakenDBIndex &index, size_t index_file_size) {
  if (index_file_size < 8)
    errx(EX_DATAERR, "database index is truncated");
  const uint8_t indexed_nt = index.indexed_nt();
  if (indexed_nt == 0 || indexed_nt > 15 || indexed_nt > database.get_k())
    errx(EX_DATAERR, "database index has an invalid minimizer length");

  const uint64_t bin_count = 1ull << (indexed_nt * 2);
  if (bin_count > (SIZE_MAX - 16) / sizeof(uint64_t))
    errx(EX_DATAERR, "database index is too large");
  const size_t required_index_size =
      8 + (size_t)(bin_count + 1) * sizeof(uint64_t);
  if (index_file_size < required_index_size)
    errx(EX_DATAERR, "database index is truncated");

  if (db_file_size < 72 || database.get_key_bits() == 0 ||
      database.get_key_bits() >= 64)
    errx(EX_DATAERR, "database has an invalid header");
  const size_t header_size = database.header_size();
  const uint64_t pair_size = database.pair_size();
  if (header_size > db_file_size || pair_size == 0 ||
      database.get_key_ct() > (db_file_size - header_size) / pair_size)
    errx(EX_DATAERR, "database is truncated");

  if (index.at(0) != 0)
    errx(EX_DATAERR, "database index does not start at the first key");
  uint64_t previous_offset = 0;
  const uint64_t key_mask = (1ull << database.get_key_bits()) - 1;
  for (uint64_t bin = 0; bin < bin_count; ++bin) {
    const uint64_t begin = index.at(bin);
    const uint64_t end = index.at(bin + 1);
    if (begin != previous_offset || end < begin ||
        end > database.get_key_ct())
      errx(EX_DATAERR, "database index offsets are invalid");

    uint64_t previous_key = 0;
    bool have_previous_key = false;
    char *pair = database.get_pair_ptr() + begin * pair_size;
    for (uint64_t position = begin; position < end; ++position) {
      uint64_t key = 0;
      memcpy(&key, pair, database.get_key_len());
      key &= key_mask;
      if (database.bin_key(key) != bin)
        errx(EX_DATAERR, "database index does not match database keys");
      if (have_previous_key && key < previous_key)
        errx(EX_DATAERR, "database keys are not sorted within an index bin");
      previous_key = key;
      have_previous_key = true;
      pair += pair_size;
    }
    previous_offset = end;
  }
  if (previous_offset != database.get_key_ct())
    errx(EX_DATAERR, "database index does not cover all keys");
}


int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  if (!TaxDB_filename.empty()) {
    taxdb = TaxonomyDB<uint32_t>(TaxDB_filename);
    Parent_map = taxdb.getParentMap();
  } else {
    cerr << "TaxDB argument is required!" << endl;
    return 1;
  }

  if (Use_uids_instead_of_taxids) {
    UID_map_file.open(UID_map_filename, ios_base::out | ios_base::binary);

    if (!UID_map_file.is_open()) {
      cerr << "Something went wrong while creating the file." << endl;
      exit(1);
    }
  }

  const bool write_to_output = !Output_DB_filename.empty() && !Pretend;
  const string active_db_filename =
      write_to_output ? Output_DB_filename : DB_filename;
  int db_lock_fd = open(active_db_filename.c_str(),
                        O_RDWR | (write_to_output ? O_CREAT : 0), 0666);
  if (db_lock_fd < 0)
    err(EX_OSERR, "unable to open database %s", active_db_filename.c_str());
  if (flock(db_lock_fd, LOCK_EX | LOCK_NB) < 0)
    err(EX_OSERR, "unable to lock database %s", active_db_filename.c_str());
  if (write_to_output) {
    copy_database_file(DB_filename, Output_DB_filename);
    DB_filename = Output_DB_filename;
  }

  QuickFile db_file(DB_filename, "rw");
  size_t db_file_size = db_file.size();
  if (db_file_size < 72)
    errx(EX_DATAERR, "database is truncated");
  vector<char> dat;
  if (Operate_in_RAM) {
    db_file.close_file();
    dat = slurp_file(DB_filename, db_file_size);
    Database = KrakenDB(dat.data(), db_file_size);
  } else {
    Database = KrakenDB(db_file.ptr(), db_file_size);
  }

  KmerScanner::set_k(Database.get_k());

  QuickFile idx_file(Index_filename);
  if (idx_file.size() < 8)
    errx(EX_DATAERR, "database index is truncated");
  KrakenDBIndex db_index(idx_file.ptr());
  Database.set_index(&db_index);
  validate_database_layout(Database, db_file_size, db_index, idx_file.size());

  // Set up the batch buffer (not used in UID mode)
  if (!Use_uids_instead_of_taxids) {
    double bytes = Memory_budget_gib * 1024.0 * 1024.0 * 1024.0;
    double max_rec_d = bytes / (double)sizeof(KmerRec);
    const double cap = (double)(SIZE_MAX / sizeof(KmerRec) / 2);
    if (max_rec_d > cap) max_rec_d = cap;
    Max_records = (size_t)max_rec_d;
    if (Max_records < 2 * (size_t)SKIP_LEN)
      Max_records = 2 * (size_t)SKIP_LEN;
    Group_chunks = (Max_records / 2) / SKIP_LEN;
    if (Group_chunks < 1) Group_chunks = 1;
    // Reserve once so the buffer never reallocates (keeps memory bounded)
    Kmer_buffer.reserve(Max_records);
  }

  if (One_FASTA_file)
    process_single_file();
  else
    process_files();

  // Apply any remaining pending updates before syncing / counting
  flush_kmers();
  { vector<KmerRec>().swap(Kmer_buffer); }

  if (!Operate_in_RAM && !Pretend &&
      msync(db_file.ptr(), db_file_size, MS_SYNC) < 0)
    err(EX_IOERR, "unable to sync database %s", DB_filename.c_str());

  if (!Kmer_count_filename.empty()) {
    ofstream ofs(Kmer_count_filename.c_str());
    cerr << "Writing kmer counts to " << Kmer_count_filename << "..." << endl;
    auto counts = Database.count_taxons();
    for (auto it = counts.begin(); it != counts.end(); ++it) {
      ofs << it->first << '\t' << it->second << '\n';
    }
    ofs.close();
  }

  if (Operate_in_RAM && !Pretend) {
    if (Output_DB_filename.size() > 0) {
      DB_filename = Output_DB_filename;
    }
    cerr << "Writing database from RAM back to " << DB_filename << " ..." << endl;
    ofstream ofs(DB_filename.c_str(), ofstream::binary);
    ofs.write(dat.data(), db_file_size);
    ofs.close();
    dat.clear();
  }

  UID_map_file.close();
  close(db_lock_fd);

  // Write new TaxDB file if new taxids were added
  if ((Add_taxIds_for_Sequences || Add_taxIds_for_Assembly) && !TaxDB_filename.empty() && !Pretend) {
    cerr << "Writing new TaxDB ..." << endl;
    ofstream ofs(TaxDB_filename.c_str());
    taxdb.writeTaxonomyIndex(ofs);
    ofs.close();
  }

  return 0;
}

inline 
uint32_t get_new_taxid(
    unordered_map<string, uint32_t>& name_to_taxid_map, 
    unordered_map<uint32_t,uint32_t>& Parent_map,
    string name, uint32_t parent_taxid, const string & rank_name) {

  auto it = name_to_taxid_map.find(name);
  if (it == name_to_taxid_map.end()) {
    uint32_t new_taxid = ++New_taxid_start;
    bool insert_res = taxdb.insert(new_taxid, parent_taxid, rank_name, name);
    if (!insert_res) {
      return 0;
    }
    Parent_map[new_taxid] = parent_taxid;
    name_to_taxid_map[name] = new_taxid;
    return new_taxid;
   } else {
    return it->second;
   }
}

unordered_map<string,uint32_t> read_seqid_to_taxid_map(string ID_to_taxon_map_filename, 
    TaxonomyDB<uint32_t>& taxdb, unordered_map<uint32_t,uint32_t>& Parent_map, 
    bool Add_taxIds_for_Assembly, bool Add_taxIds_for_Sequences) {

  cerr << "Reading sequence ID to taxonomy ID mapping ... ";

  unordered_map<string, uint32_t> ID_to_taxon_map;
  ifstream map_file(ID_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", ID_to_taxon_map_filename.c_str());
  }
  string line, seq_id, name;
  uint32_t taxid;

  if (Add_taxIds_for_Assembly || Add_taxIds_for_Sequences) {
    for (auto it = taxdb.entries.begin(); it != taxdb.entries.end(); ++it) {
      if (it->first >= New_taxid_start) {
        New_taxid_start = it->first+100;
      } 
    }
    cerr << "[starting new taxonomy IDs with " << (New_taxid_start+1) << ']';
  }

  // Used when adding new taxids for assembly or sequence
  unordered_map<string, uint32_t> name_to_taxid_map;

  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    istringstream iss(line);
    iss >> seq_id >> taxid;

    auto it = ID_to_taxon_map.find(seq_id);
    if (it != ID_to_taxon_map.end()) {
      // The sequence ID has been seen before, ignore
      continue;
    }

    uint32_t orig_taxid = taxid;

    if (Add_taxIds_for_Assembly && iss.good()) {
      iss.get();
      getline(iss, name);
      if (!name.empty())
        taxid = get_new_taxid(name_to_taxid_map, Parent_map, name, taxid, "assembly");
    }

    if (Add_taxIds_for_Sequences && orig_taxid != TID_HUMAN && orig_taxid != TID_MOUSE) {
      taxid = get_new_taxid(name_to_taxid_map, Parent_map, seq_id, taxid, "sequence");
    }
    if (Add_taxIds_for_Assembly || Add_taxIds_for_Sequences) {
      cout << seq_id << '\t' << taxid << '\n';
    }
    ID_to_taxon_map[seq_id] = taxid;
  }
  if (ID_to_taxon_map.size() == 0) {
    cerr << "Error: No ID mappings present!!" << endl;
  }
  cerr << " got " << ID_to_taxon_map.size() << " mappings." << endl;
  return ID_to_taxon_map;
}

void process_single_file() {
  ID_to_taxon_map = read_seqid_to_taxid_map(ID_to_taxon_map_filename, taxdb, Parent_map, Add_taxIds_for_Assembly, Add_taxIds_for_Sequences);

  FastaReader reader(Multi_fasta_filename);
  DNASequence dna;
  uint32_t seqs_processed = 0;
  uint32_t seqs_skipped = 0;
  uint32_t seqs_no_taxid = 0;

  while (reader.is_valid()) {
    dna = reader.next_sequence();
    if (! reader.is_valid())
      break;

    if ( dna.seq.empty() ) {
      ++seqs_skipped;
      continue;
    }

    // Get the taxid. If the header specifies kraken:taxid, use that
    uint32_t taxid = 0;
    auto it = ID_to_taxon_map.find(dna.id);
    if (it != ID_to_taxon_map.end()) {
      taxid = it->second;
    } else {
      // Check if it can be mapped by removing .\d suffix
      size_t pos = dna.id.find_last_of('.');
      bool is_num_suffix = pos != string::npos;
      if (is_num_suffix) {
        for (size_t i = pos + 1; i < dna.id.size(); ++i) {
          if (!isdigit(dna.id[i])) {
            is_num_suffix = false;
            break;
          }
        }
      }
      if (is_num_suffix) {
        it = ID_to_taxon_map.find(dna.id.substr(0, pos));
        if (it != ID_to_taxon_map.end()) {
          taxid = it->second;
        }
      }
    }
    
    if (taxid == 0 && dna.id.size() >= prefix.size() && dna.id.substr(0,prefix.size()) == prefix) {
      // if the AC is not in the map, check if the fasta entry starts with '>kraken:taxid'
        taxid = std::stol(dna.id.substr(prefix.size()));
        if (taxid == 0) {
          cerr << "Error: taxonomy ID is zero for sequence '" << dna.id << "'?!" << endl;
        }
        const auto strBegin = dna.header_line.find_first_not_of("\t ");
        if (strBegin != std::string::npos)
            dna.header_line = dna.header_line.substr(strBegin);
    } 
    
    if (taxid == 0) {
        cerr << "Error! Didn't find taxonomy ID mapping for sequence " <<  dna.id << "!!" << endl;
        ++seqs_skipped;
        continue;
    }

    if (Minimum_sequence_size > 0 && dna.seq.size() < Minimum_sequence_size) {
      cerr << "Skipping sequence " << dna.id << " as it's too short (" << dna.seq.size() << ")" << endl;
      ++ seqs_skipped;
      continue;
    }

    auto it_p = Parent_map.find(taxid);
    if (it_p == Parent_map.end()) {
      cerr << "Skipping sequence " << dna.id << " since taxonomy ID " << taxid << " is not in taxonomy database!" << endl;
      ++ seqs_skipped;
      continue;
    }
    
    bool is_contaminant_taxid = taxid == TID_CONTAMINANT1 || taxid == TID_CONTAMINANT2;
    if (Add_taxIds_for_Sequences && taxid != TID_HUMAN && it_p->second != TID_HUMAN && taxid != TID_MOUSE && it_p->second != TID_MOUSE) {
      // Update entry based on header line
      auto entryIt = taxdb.entries.find(taxid);
      if (entryIt == taxdb.entries.end()) {
        cerr << "Error! Didn't find taxid " << taxid << " in TaxonomyDB - can't update it!! ["<<dna.header_line<<"]" << endl;
      } else {
        entryIt->second.scientificName = dna.header_line;
      }
    }

    if (taxid) {
      if (Parent_map.find(taxid) == Parent_map.end() || taxdb.entries.find(taxid) == taxdb.entries.end()) {
        cerr << "Ignoring sequence for taxID " << taxid << " - not in taxDB\n";
      } else {
        add_sequence(taxid, dna.seq, is_contaminant_taxid);
        ++seqs_processed;
      }
    } else {
      if (verbose) 
        cerr << "Skipping sequence with header [" << dna.header_line << "] - no taxid" << endl;
      ++seqs_no_taxid;
    }

    cerr << "\rProcessed " << seqs_processed << " sequences";
  }
  cerr << "\r                                                                            ";
  cerr << "\rFinished processing " << seqs_processed << " sequences (skipping "<< seqs_skipped <<" empty sequences, and " << seqs_no_taxid<<" sequences with no taxonomy mapping)" << endl;
}

void process_files() {
  cerr << "Processing files in " << File_to_taxon_map_filename.c_str() << endl;
  ifstream map_file(File_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", File_to_taxon_map_filename.c_str());
  }
  string line;
  uint32_t seqs_processed = 0;

  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    string filename;
    uint32_t taxid;
    istringstream iss(line);
    iss >> filename;
    iss >> taxid;
    // TODO: Support a mapping file with only file names, not taxids
    process_file(filename, taxid);
    cerr << "\rProcessed " << ++seqs_processed << " sequences";
  }
  cerr << "\r                                                       ";
  cerr << "\rFinished processing " << seqs_processed << " sequences" << endl;
}

void process_file(string filename, uint32_t taxid) {
  FastaReader reader(filename);
  DNASequence dna;
  
  // For the purposes of this program, we assume these files are
  // single-fasta files.
  dna = reader.next_sequence();

  add_sequence(taxid, dna.seq, false);
}

// ---------------------------------------------------------------------
// Batched path
// ---------------------------------------------------------------------

// Scan a region and append pending updates to the shared buffer.
static void collect_kmers(uint32_t seq_idx, string &seq, size_t start, size_t finish) {
  KmerScanner scanner(seq, start, finish);
  uint64_t *kmer_ptr;
  vector<KmerRec> local;
  local.reserve(SKIP_LEN + 64);

  while ((kmer_ptr = scanner.next_kmer()) != NULL) {
    if (scanner.ambig_kmer())
      continue;
    KmerRec r;
    r.kmer = Database.canonical_representation(*kmer_ptr);
    r.bin = (uint32_t)Database.bin_key(r.kmer);
    r.seq = seq_idx;
    local.push_back(r);
  }
  if (local.empty())
    return;
#ifdef _OPENMP
  #pragma omp critical(kmer_buffer_append)
#endif
  Kmer_buffer.insert(Kmer_buffer.end(), local.begin(), local.end());
}

// Sort pending updates into DB order and apply them in a single pass.
// Updates to the same k-mer are applied in original (sequence) order.
void flush_kmers() {
  if (Kmer_buffer.empty())
    return;
  sort(Kmer_buffer.begin(), Kmer_buffer.end(), kmer_rec_less);

  KrakenDBIndex *index = Database.get_index();
  char *pairs = Database.get_pair_ptr();
  const uint64_t pair_sz = Database.pair_size();
  const uint64_t key_len = Database.get_key_len();
  const uint64_t key_mask = (1ull << Database.get_key_bits()) - 1;

  const size_t n = Kmer_buffer.size();
  bool have_bin = false;
  uint32_t cur_bin = 0;
  uint64_t bin_end = 0;
  uint64_t cursor = 0;  // keys before cursor in this bin are < current k-mer

  size_t i = 0;
  while (i < n) {
    const uint64_t kmer = Kmer_buffer[i].kmer;
    const uint32_t bin = Kmer_buffer[i].bin;
    size_t j = i + 1;
    while (j < n && Kmer_buffer[j].kmer == kmer)
      ++j;

    if (!have_bin || bin != cur_bin) {
      have_bin = true;
      cur_bin = bin;
      cursor = index->at(bin);
      bin_end = index->at((uint64_t)bin + 1);
    }

    // lower_bound for kmer within [cursor, bin_end)
    uint64_t lo = cursor, hi = bin_end;
    while (lo < hi) {
      const uint64_t mid = lo + (hi - lo) / 2;
      uint64_t key = 0;
      memcpy(&key, pairs + pair_sz * mid, key_len);
      key &= key_mask;
      if (key < kmer)
        lo = mid + 1;
      else
        hi = mid;
    }
    bool found = false;
    if (lo < bin_end) {
      uint64_t key = 0;
      memcpy(&key, pairs + pair_sz * lo, key_len);
      key &= key_mask;
      found = (key == kmer);
    }

    if (!found) {
      cursor = lo;
      if (!Allow_extra_kmers) {
        errx(EX_DATAERR, "kmer found in sequence that is not in database");
      } else if (verbose) {
        for (size_t x = i; x < j; ++x)
          cerr << "kmer found in sequence w/ taxid " << Seq_info[Kmer_buffer[x].seq].taxid
               << " that is not in database" << endl;
      }
    } else {
      char *val_loc = pairs + pair_sz * lo + key_len;
      uint32_t original;
      memcpy(&original, val_loc, sizeof(original));
      uint32_t val = original;

      bool have_prev = false;
      uint32_t prev_taxid = 0;
      bool prev_contam = false;
      for (size_t x = i; x < j; ++x) {
        const SeqInfo &si = Seq_info[Kmer_buffer[x].seq];
        // Consecutive identical updates are idempotent - skip them
        if (have_prev && si.taxid == prev_taxid && si.is_contaminant == prev_contam)
          continue;
        apply_update(val, si.taxid, si.is_contaminant);
        have_prev = true;
        prev_taxid = si.taxid;
        prev_contam = si.is_contaminant;
      }
      // Only write when changed, to avoid dirtying mmapped pages needlessly
      if (val != original)
        memcpy(val_loc, &val, sizeof(val));
      cursor = lo + 1;
    }
    i = j;
  }
  Kmer_buffer.clear();
}

// UID mode applies updates immediately (UID numbering depends on encounter
// order); otherwise updates are batched.
void add_sequence(uint32_t taxid, string &seq, bool is_contaminant_taxid) {
  if (Use_uids_instead_of_taxids) {
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < seq.size(); i += SKIP_LEN)
      set_lcas(taxid, seq, i, i + SKIP_LEN + Database.get_k() - 1, is_contaminant_taxid);
    return;
  }

  // No pending records means no outstanding references into Seq_info
  if (Kmer_buffer.empty())
    Seq_info.clear();

  const uint32_t seq_idx = (uint32_t)Seq_info.size();
  SeqInfo si;
  si.taxid = taxid;
  si.is_contaminant = is_contaminant_taxid;
  Seq_info.push_back(si);

  const size_t group_bases = Group_chunks * (size_t)SKIP_LEN;
  for (size_t g = 0; g < seq.size(); g += group_bases) {
    const size_t gend = min(seq.size(), g + group_bases);
    // Upper bound on k-mers produced by this group
    const size_t need = gend - g;
    if (Kmer_buffer.size() + need > Max_records)
      flush_kmers();
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = g; i < gend; i += SKIP_LEN)
      collect_kmers(seq_idx, seq, i, i + SKIP_LEN + Database.get_k() - 1);
  }
}

// Original immediate-update implementation; only used for the UID path.
void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish, bool is_contaminant_taxid) {
  KmerScanner scanner(seq, start, finish);
  uint64_t *kmer_ptr;
  uint32_t *val_ptr;

  while ((kmer_ptr = scanner.next_kmer()) != NULL) {
    if (scanner.ambig_kmer())
      continue;
    val_ptr = Database.kmer_query(
                Database.canonical_representation(*kmer_ptr)
    );
    if (val_ptr == NULL) {
      if (! Allow_extra_kmers) {
        errx(EX_DATAERR, "kmer found in sequence that is not in database");
      } 
      else if (verbose) {
        cerr << "kmer found in sequence w/ taxid " << taxid << " that is not in database" << endl;
      }
      continue;
    }

#ifdef _OPENMP
    #pragma omp critical(new_uid)
#endif
    *val_ptr = uid_mapping(Taxids_to_UID_map, UID_to_taxids_vec, taxid, *val_ptr, current_uid, UID_map_file);
  }
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "f:d:i:t:n:m:F:xMTRvb:aApI:o:Sc:E:B:")) != -1) {
    switch (opt) {
      case 'f' :
        File_to_taxon_map_filename = optarg;
        break;
      case 'I' :
        Use_uids_instead_of_taxids = true;
        UID_map_filename = optarg;
        break;
      case 'd' :
        DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 'F' :
        Multi_fasta_filename = optarg;
        break;
      case 'm' :
        ID_to_taxon_map_filename = optarg;
        break;
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
      case 'T' :
        Force_contaminant_taxid = true;
        break;
      case 'R' :
        Reset_taxid = true;
        break;
      case 'v' :
        verbose = true;
        break;
      case 'x' :
        Allow_extra_kmers = true;
        break;
      case 'a' :
        Add_taxIds_for_Sequences = true;
        break;
      case 'A' :
        Add_taxIds_for_Assembly = true;
        break;
      case 'b' :
        TaxDB_filename = optarg;
        break;
      case 'c' :
        Kmer_count_filename = optarg;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'E' :
        Minimum_sequence_size = atoi(optarg);
        break;
      case 'B' : {
        char *endp = NULL;
        errno = 0;
        double gib = strtod(optarg, &endp);
        if (errno || endp == optarg || *endp != '\0' ||
            !std::isfinite(gib) || gib <= 0)
          errx(EX_USAGE, "memory budget must be a positive number of GiB");
        Memory_budget_gib = gib;
        break;
      }
      case 'p' :
        Pretend = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filename.empty() || Index_filename.empty() ||
      TaxDB_filename.empty())
    usage();
  if (File_to_taxon_map_filename.empty() &&
      (Multi_fasta_filename.empty() || ID_to_taxon_map_filename.empty()))
    usage();

  if (! File_to_taxon_map_filename.empty())
    One_FASTA_file = false;
  else
    One_FASTA_file = true;
}

void usage(int exit_code) {
  cerr << "Usage: set_lcas [options]" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "* -b filename      Taxonomy DB file" << endl
       << "  -t #             Number of threads" << endl
       << "  -M               Copy DB to RAM during operation" << endl
       << "  -B #             Memory budget in GiB for the sorted k-mer update buffer (default 1; fractions allowed)" << endl
       << "  -o filename      Output database to filename, instead of overwriting the input database" << endl
       << "  -x               K-mers not found in DB do not cause errors" << endl
       << "  -f filename      File to taxon map" << endl
       << "  -F filename      Multi-FASTA file with sequence data" << endl
       << "  -m filename      Sequence ID to taxon map" << endl
       << "  -a               Add taxonomy IDs (starting with "<<(New_taxid_start+1)<<") for assemblies (third column in seqid2taxid.map) to Taxonomy DB" << endl
       << "  -A               Add taxonomy IDs (starting with "<<(New_taxid_start+1)<<") for sequences to Taxonomy DB" << endl
       << "  -T               When a k-mer appears in a 'synthetic construct' sequence, force the taxID to be the 'synthetic construct' taxID, instead of the LCA." << endl
       << "  -E #             Exclude sequences that are shorter than the threshold." << endl
       << "  -I filename      Write UIDs into database, and output (binary) UID-to-taxid map to filename" << endl
       << "  -p               Pretend - do not write database back to disk (when working in RAM)" << endl
       << "  -v               Verbose output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "-F and -m must be specified together.  If -f is given, "
       << "-F/-m are ignored." << endl;
  exit(exit_code);
}
