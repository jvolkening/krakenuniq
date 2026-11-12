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

#define __STDC_FORMAT_MACROS 1 // for PRIu64, etc.

#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"
#include "readcounts.hpp"
#include "taxdb.hpp"
#include "gzstream.h"
#include "uid_mapping.hpp"
#include <sstream>
#include <sys/mman.h>
#if defined(_OPENMP) && defined(__GLIBCXX__)
#include <parallel/algorithm>
#endif
#include <inttypes.h>
#include <cassert>
#include <cstdio>
#include <inttypes.h>

const size_t DEF_WORK_UNIT_SIZE = 500000;
const size_t DEF_BATCH_SIZE = 16000000;
int New_taxid_start = 1000000000;

using namespace std;
using namespace kraken;

#define USE_KHSET_FOR_EXACT_COUNTING

#ifdef EXACT_COUNTING
  #ifdef USE_KHSET_FOR_EXACT_COUNTING
    #include "khset.h"
    using READCOUNTS = ReadCounts< kh::khset64_t >;
  #else
    #include <unordered_set>
    using READCOUNTS = ReadCounts< unordered_set<uint64_t> >;
  #endif
#else
  using READCOUNTS = ReadCounts<HyperLogLogPlusMinus<uint64_t> >;
#endif


#ifndef _OPENMP
  int omp_get_thread_num() { return 0; }
#endif

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename);
void process_file_with_db_chunk(char *filename);
void classify_sequence_with_db_chunk(std::pair<DNASequence, uint32_t> & seq, std::fstream & fp, const uint32_t db_chunk_id, const uint32_t db_id);
bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, READCOUNTS>&);
bool finish_classification(DNASequence &dna, const vector<uint32_t> &taxa,
                           const vector<char> &ambig_list,
                           unordered_map<uint32_t, uint32_t> &hit_counts,
                           uint32_t hits, uint32_t taxon, ostringstream &koss,
                           ostringstream &coss, ostringstream &uoss,
                           unordered_map<uint32_t, READCOUNTS>& my_taxon_counts);
void process_file_batched(char *filename);
inline void print_sequence(ostream* oss_ptr, const DNASequence& dna);
string hitlist_string(const vector<uint32_t> &taxa, const vector<char>& ambig_list);


set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);
double get_seconds(struct timeval time1, struct timeval time2);
unordered_map<uint32_t, READCOUNTS> taxon_counts; // stats per taxon

int Num_threads = 1;
vector<string> DB_filenames;
vector<string> Index_filenames;
bool Quick_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Print_kraken_report = false;
bool Populate_memory = false;
uint64_t Populate_memory_size = 0;
bool Only_classified_kraken_output = false;
bool Print_sequence = false;
bool Print_Progress = true;
bool full_report = false;

bool Map_UIDs = false;
string UID_to_TaxID_map_filename;
map<uint32_t, vector<uint32_t> > UID_to_taxids_map;
QuickFile UID_to_TaxID_map_file;

uint32_t Minimum_hit_count = 1;
unordered_map<uint32_t, uint32_t> Parent_map;
unordered_map<uint32_t, vector<uint32_t> > Uid_dict;
string Classified_output_file, Unclassified_output_file, Kraken_output_file, Report_output_file, TaxDB_file;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Kraken_output;
ostream *Report_output;
vector<ofstream*> Open_fstreams;
vector<ogzstream*> Open_gzstreams;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;
size_t Batch_size = DEF_BATCH_SIZE;  // bp per batch of reads; 0 = classify read by read
bool Page_hints = true;
// Mapped DB/index files (address, size), filled in by main()
static vector<std::pair<const char *, size_t> > Db_maps, Idx_maps;
TaxonomyDB<uint32_t> taxdb;
static vector<KrakenDB*> KrakenDatabases (DB_filenames.size());

struct db_status {
  db_status() : current_bin_key(0), current_min_pos(1), current_max_pos(0) {}
  uint64_t current_bin_key;
  int64_t current_min_pos;
  int64_t current_max_pos;
};

unsigned long long total_classified = 0;
unsigned long long total_sequences = 0;
unsigned long long total_bases = 0;
uint32_t ambig_taxon = -1;

inline bool ends_with(std::string const & value, std::string const & ending)
{
        if (ending.size() > value.size()) return false;
            return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

ostream* cout_or_file(string file, bool append = false) {
    if (file == "-")
      return &cout;

    if (ends_with(file, ".gz")) {
      ogzstream* ogzs = new ogzstream(file.c_str());
      ogzs->exceptions( ifstream::failbit | ifstream::badbit );
      Open_gzstreams.push_back(ogzs);
      return ogzs;
    } else {
      ofstream* ofs = append? new ofstream(file.c_str(), std::ofstream::app) : new ofstream(file.c_str());
      ofs->exceptions( ifstream::failbit | ifstream::badbit );
      Open_fstreams.push_back(ofs);
      return ofs;
    }
}

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);
  
  if (Map_UIDs) {
    if (DB_filenames.size() > 1) {
      cerr << "Cannot use more than one database with UID mapping!" << endl;
      return 1;
    }

    cerr << "Reading UID mapping file " << UID_to_TaxID_map_filename << endl;
    UID_to_TaxID_map_file.open_file(UID_to_TaxID_map_filename);

    // Always Populate memory
    //if (Populate_memory) {
    UID_to_TaxID_map_file.load_file();
    //}
  }

  if (Populate_memory && Populate_memory_size == 0)
    cerr << "Loading database(s)... " << endl;

  static vector<QuickFile> idx_files (DB_filenames.size());
  static vector<QuickFile> db_files (DB_filenames.size());
  static vector<KrakenDBIndex> db_indices (DB_filenames.size());


  const bool batch_mode = Batch_size > 0 && !Quick_mode && !(Populate_memory && Populate_memory_size > 0);
  if (batch_mode)
    QuickFile::willneed_on_open = false;  // don't read whole (possibly huge) files in

  // TODO: Check DB_filenames and Index_filesnames have the same length
  for (size_t i=0; i < DB_filenames.size(); ++i) {
    cerr << " Database " << DB_filenames[i] << endl;
    db_files[i].open_file(DB_filenames[i]);

    KrakenDatabases.push_back(new KrakenDB(db_files[i].ptr()));
    idx_files[i].open_file(Index_filenames[i]);
    // NOTE: we switched the order, i.e., we are creating the objects before loading everything into main memory
    db_indices[i] = KrakenDBIndex(idx_files[i].ptr());
    KrakenDatabases[i]->set_index(&db_indices[i]);
    Db_maps.push_back(std::make_pair((const char *) db_files[i].ptr(), db_files[i].size()));
    Idx_maps.push_back(std::make_pair((const char *) idx_files[i].ptr(), idx_files[i].size()));

    if (Populate_memory && Populate_memory_size == 0) // only when no chunk size is passed!
    {
      db_files[i].load_file();
      idx_files[i].load_file();
    }
    else if (Populate_memory && Populate_memory_size > 0)
    {
      KrakenDatabases[i]->prepare_chunking(Populate_memory_size);
    }
  }

  // Check all databases have the same k
  uint8_t kmer_size = KrakenDatabases[0]->get_k();
  for (size_t i = 1; i < KrakenDatabases.size(); ++i) {
    uint8_t kmer_size_i = KrakenDatabases[i]->get_k();
    if (kmer_size_i != kmer_size) {
      fprintf(stderr, "Different k-mer sizes in databases 1 and %lu: %i vs %i!\n", i+1, (int)kmer_size, (int)kmer_size_i);
      exit(1);
    }
  };
  KmerScanner::set_k(kmer_size);

  if (Populate_memory && Populate_memory_size == 0)
    cerr << "\ncomplete." << endl;


  if (!TaxDB_file.empty()) {
      taxdb = TaxonomyDB<uint32_t>(TaxDB_file, false);
      Parent_map = taxdb.getParentMap();
  } else {
      cerr << "TaxDB argument is required!" << endl;
      return 1;
  }

  if (Print_classified) {
    Classified_output = cout_or_file(Classified_output_file);
  }

  if (Print_unclassified) {
    Unclassified_output = cout_or_file(Unclassified_output_file);
  }

  if (! Kraken_output_file.empty()) {
    if (Kraken_output_file == "off" || Kraken_output_file == "-") {
      Print_kraken = false;
    //else if (Kraken_output_file == "-") {
    //  Kraken_output = &cout;
    } else {
      cerr << "Writing Kraken output to " << Kraken_output_file << endl;
      Kraken_output = cout_or_file(Kraken_output_file);
    }
  } else {
    Kraken_output = &cout;
  }

  //cerr << "Print_kraken: " << Print_kraken << "; Print_kraken_report: " << Print_kraken_report << "; k: " << uint32_t(KrakenDatabases[0]->get_k()) << endl;

  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (int i = optind; i < argc; i++) {
    if (Populate_memory && Populate_memory_size > 0)
      process_file_with_db_chunk(argv[i]);
    else if (Batch_size > 0 && !Quick_mode)
      process_file_batched(argv[i]);
    else
      process_file(argv[i]);
  }
  gettimeofday(&tv2, NULL);

  report_stats(tv1, tv2);

  if (!Report_output_file.empty() && Report_output_file != "off") {
    gettimeofday(&tv1, NULL);
    std::cerr << "Writing report file to " << Report_output_file <<"  ..\n";
    for (size_t i = 0; i < DB_filenames.size(); ++i) {
      const auto fname = DB_filenames[i] + ".counts";
      ifstream ifs(fname);
      bool counts_file_gd = false;
      if (ifs.good()) {
        if (ifs.peek() == std::ifstream::traits_type::eof()) {
          cerr << "Kmer counts file is empty - trying to regenerate ..." << endl;
        } else {
          ifs.close();
          counts_file_gd = true;
        }
      } 
      if (!counts_file_gd) {
        ofstream ofs(fname);
        cerr << "Writing kmer counts to " << fname << "... [only once for this database, may take a while] " << endl;
        auto counts = KrakenDatabases[i]->count_taxons();
        for (auto it = counts.begin(); it != counts.end(); ++it) {
          ofs << it->first << '\t' << it->second << '\n';
        }
        ofs.close();
      }
      taxdb.readGenomeSizes(fname);
    }
     Report_output = cout_or_file(Report_output_file, true);
  
    TaxReport<uint32_t,READCOUNTS> rep = TaxReport<uint32_t, READCOUNTS>(*Report_output, taxdb, taxon_counts, false);
    if (HLL_PRECISION > 0) {
      if (full_report) {
        rep.setReportCols(vector<string> { 
          "%",
          "reads", 
          "taxReads",
          "kmers",
          "taxKmers",
          "kmersDB",
          "taxKmersDB",
          "dup",
          "cov", 
          "taxID", 
          "rank", 
          "taxName"});
      } else {
        rep.setReportCols(vector<string> { 
          "%",
          "reads", 
          "taxReads",
          "kmers",
          "dup",
          "cov", 
          "taxID", 
          "rank", 
          "taxName"});
      }
    } else {
      rep.setReportCols(vector<string> { 
        "%",
        "reads", 
        "taxReads",
        "taxID", 
        "rank", 
        "taxName"});
    }
    rep.printReport("kraken");
    gettimeofday(&tv2, NULL);
    fprintf(stderr, "Report finished in %.3f seconds.\n", get_seconds(tv1,tv2));
  }
  cerr << "Finishing up ..." << endl;

  for (size_t i = 0; i < Open_fstreams.size(); ++i) {
    ofstream* ofs = Open_fstreams[i];
    ofs->close();
  }

  for (size_t i = 0; i < Open_gzstreams.size(); ++i) {
    ogzstream* ogzs = Open_gzstreams[i];
    ogzs->close();
  }

  for (size_t i=0; i < KrakenDatabases.size(); ++i) {
    delete KrakenDatabases[i];
  }

  return 0;
}

double get_seconds(struct timeval time1, struct timeval time2) {
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;
  return(seconds);
}

void report_stats(struct timeval time1, struct timeval time2) {
  double seconds = get_seconds(time1, time2);

  cerr << "\r";
  fprintf(stderr, 
          "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
          (unsigned long long) total_sequences, total_bases / 1.0e6, seconds,
          total_sequences / 1.0e3 / (seconds / 60),
          total_bases / 1.0e6 / (seconds / 60) );
  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) (total_sequences - total_classified),
          (total_sequences - total_classified) * 100.0 / total_sequences);
}

bool determine_input_file_type(char* filename)
{
  bxz::ifstream file;
  file.open(filename);
  if (file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", filename);
  }
  std::istreambuf_iterator<char> file_iter(file);
  bool ret = *file_iter == '@';
  file.close();
  return ret;
};

void merge_intermediate_results_by_workers(const bool first_intermediate_output, const std::string & tmp_file_name) {
  const std::string filename_merged_summary = tmp_file_name;

  FILE *fp_prev_merged_summary = NULL;
  std::string filename_prev_merged_summary;
  if (!first_intermediate_output)
  {
    filename_prev_merged_summary = filename_merged_summary + ".prev";
    rename(filename_merged_summary.c_str(), filename_prev_merged_summary.c_str());
    fp_prev_merged_summary = fopen(filename_prev_merged_summary.c_str(), "rb");
  }

  std::fstream fp_merged_summary(filename_merged_summary, std::fstream::out | std::fstream::binary | std::fstream::trunc);
  fp_merged_summary.exceptions(std::fstream::badbit);

  std::vector<FILE*> worker_files(Num_threads);
  std::vector<uint32_t> next_read_id(Num_threads);
  for (int worker_id = 0; worker_id < Num_threads; ++worker_id)
  {
    const std::string worker_filename = filename_merged_summary + "." + std::to_string(worker_id);
    worker_files[worker_id] = fopen(worker_filename.c_str(), "rb");

    if (fread(&next_read_id[worker_id], sizeof(next_read_id[worker_id]), 1, worker_files[worker_id]) != 1)
      next_read_id[worker_id] = 0; // 0 indicates that there is no read left in this file (reads are indexed starting at 1)
  }

  std::vector<uint32_t> taxa, taxa_prev_summary;
  uint32_t seq_idx = 1;
  while (true)
  {
    // determine worker to get this read info from
    int worker_to_retrieve_from = -1;
    for (int worker_id = 0; worker_id < Num_threads; ++worker_id) {
      if (seq_idx == next_read_id[worker_id]) {
        worker_to_retrieve_from = worker_id;
        break;
      }
    }
    if (worker_to_retrieve_from == -1)
      break; // finished processing all files

    // retrieve from that worker and load next sequence name
    uint32_t taxa_size;
    fread(&taxa_size, sizeof(uint32_t), 1, worker_files[worker_to_retrieve_from]);
    taxa.resize(taxa_size);
    fread(&taxa[0], sizeof(uint32_t), taxa_size, worker_files[worker_to_retrieve_from]);

    if (!first_intermediate_output)
    {
      uint32_t taxa_prev_size;
      fread(&taxa_prev_size, sizeof(uint32_t), 1, fp_prev_merged_summary);
      assert(taxa_size == taxa_prev_size);
      taxa_prev_summary.resize(taxa_prev_size);
      fread(&taxa_prev_summary[0], sizeof(uint32_t), taxa_prev_size, fp_prev_merged_summary);

      for (uint32_t i = 0; i < taxa_prev_summary.size(); ++i)
      {
        assert(!(taxa[i] && taxa_prev_summary[i])); // a kmer cannot be found in more than one database chunk
        if (taxa_prev_summary[i])
        {
          taxa[i] = taxa_prev_summary[i]; // if previous chunk/database got a match, keep the tax info
        }
      }
    }

    try {
      fp_merged_summary.write((char*) &taxa_size, 1 * sizeof(uint32_t));
      fp_merged_summary.write((char*) &taxa[0], taxa_size * sizeof(uint32_t));
    }
    catch (const std::fstream::failure& e) {
      printf("ERROR: Could not write to temporary file: %s\n", e.what());
      exit(1);
    }

    // load the next read id from that file
    if (fread(&next_read_id[worker_to_retrieve_from], sizeof(next_read_id[worker_to_retrieve_from]), 1, worker_files[worker_to_retrieve_from]) != 1)
      next_read_id[worker_to_retrieve_from] = 0; // 0 indicates that there is no read left in this file (reads are indexed starting at 1)

    ++seq_idx;
  }

  for (int worker_id = 0; worker_id < Num_threads; ++worker_id)
  {
    fclose(worker_files[worker_id]);
    const std::string worker_filename = filename_merged_summary + "." + std::to_string(worker_id);
    unlink(worker_filename.c_str());
  }

  if (!first_intermediate_output)
  {
    fclose(fp_prev_merged_summary);
    unlink(filename_prev_merged_summary.c_str());
  }

  fp_merged_summary.close();
}

void process_file(char *filename) {
  string file_str(filename);
  DNASequenceReader *reader;
  DNASequence dna;

  Fastq_input = determine_input_file_type(filename);

  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    vector<DNASequence> work_unit;
    ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;

    while (reader->is_valid()) {
      work_unit.clear();
      size_t total_nt = 0;

#ifdef _OPENMP
      #pragma omp critical(get_input)
#endif
      {
        while (total_nt < Work_unit_size) {
          dna = reader->next_sequence();
          if (! reader->is_valid())
            break;
          work_unit.push_back(dna);
          total_nt += dna.seq.size();
        }
      }
      if (total_nt == 0)
        break;
      
      unordered_map<uint32_t, READCOUNTS> my_taxon_counts;
      uint64_t my_total_classified = 0;
      kraken_output_ss.str("");
      classified_output_ss.str("");
      unclassified_output_ss.str("");
      for (size_t j = 0; j < work_unit.size(); j++) {
        my_total_classified += 
            classify_sequence( work_unit[j], kraken_output_ss,
                           classified_output_ss, unclassified_output_ss,
                           my_taxon_counts);
      }
 
#ifdef _OPENMP
      #pragma omp critical(write_output)
#endif
      {
        total_classified += my_total_classified;
        for (auto it = my_taxon_counts.begin(); it != my_taxon_counts.end(); ++it) {
          taxon_counts[it->first] += std::move(it->second);
        }

        if (Print_kraken)
          (*Kraken_output) << kraken_output_ss.str();
        if (Print_classified)
          (*Classified_output) << classified_output_ss.str();
        if (Print_unclassified)
          (*Unclassified_output) << unclassified_output_ss.str();
        total_sequences += work_unit.size();
        total_bases += total_nt;
        //if (Print_Progress && total_sequences % 100000 < work_unit.size()) 
        if (Print_Progress) {  
          fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
                          total_sequences, total_classified * 100.0 / total_sequences);
        }
      }
    }
  }  // end parallel section

  delete reader;
}

// ---------------------------------------------------------------------------
// Batched classification.
//
// Instead of looking up each k-mer in every database as it is scanned (random
// access into databases that may be larger than RAM), a batch of reads is
// processed in phases: (1) extract all canonical k-mers, (2) for each database
// in turn, sort the still-unresolved k-mers by minimizer bin and look them up
// in that order (so the DB file is read front to back, duplicate k-mers are
// looked up once, and cache-line prefetches / paging hints are issued ahead of
// the lookups), (3) assemble per-read results exactly as classify_sequence does.
// The output is identical to the unbatched code (modulo read order with >1 thread).
// ---------------------------------------------------------------------------

namespace {

struct BatchQuery {
  uint64_t bin;
  uint64_t kmer;
  uint32_t idx;   // position in the batch's flat k-mer array
};

inline bool batch_query_less(const BatchQuery &a, const BatchQuery &b) {
  return a.bin != b.bin ? a.bin < b.bin : a.kmer < b.kmer;
}

enum : uint8_t { KM_PENDING = 0, KM_AMBIG = 1, KM_RESOLVED = 2 };

const size_t PF_DIST_INDEX = 32;     // queries ahead: prefetch index cache line
const size_t PF_DIST_DATA = 16;      // queries ahead: prefetch DB cache lines
const size_t HINT_DIST_INDEX = 1024; // queries ahead: page-in hint for index
const size_t HINT_DIST_DATA = 256;   // queries ahead: page-in hint for DB
const size_t QUERY_CHUNK = 2048;     // queries per work item
const size_t MAX_HINT_BYTES = 1 << 20;

// Fraction of (sampled) pages of a mapping currently in the page cache
double resident_fraction(const char *base, size_t size) {
  const size_t page = getpagesize();
  const size_t n_pages = size / page;
  if (n_pages == 0)
    return 1.0;
  const size_t samples = std::min<size_t>(n_pages, 512);
  size_t resident = 0, tested = 0;
  for (size_t s = 0; s < samples; ++s) {
    size_t p = (n_pages / samples) * s + (s * 7919) % std::max<size_t>(n_pages / samples, 1);
    unsigned char vec = 0;
    if (mincore((void *) (base + p * page), 1, &vec) == 0) {
      ++tested;
      resident += vec & 1;
    }
  }
  return tested ? (double) resident / tested : 1.0;
}

inline void hint_range(const char *lo, const char *hi, const char **last_end) {
  if (hi <= lo || (size_t)(hi - lo) > MAX_HINT_BYTES)
    return;
  const uintptr_t page = getpagesize();
  uintptr_t a = (uintptr_t) lo & ~(page - 1);
  uintptr_t b = ((uintptr_t) hi + page - 1) & ~(page - 1);
  if ((uintptr_t) *last_end > a)  // queries are sorted: skip pages already hinted
    a = (uintptr_t) *last_end;
  if (a >= b)
    return;
  madvise((void *) a, b - a, MADV_WILLNEED);
  *last_end = (const char *) b;
}

void lookup_sorted_queries(size_t db_id, const std::vector<BatchQuery> &qs,
                           std::vector<uint32_t> &taxa, std::vector<uint8_t> &state,
                           bool hints) {
  KrakenDB *db = KrakenDatabases[db_id];
  const size_t n = qs.size();
  const size_t n_chunks = (n + QUERY_CHUNK - 1) / QUERY_CHUNK;

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t c = 0; c < n_chunks; ++c) {
    const size_t a = c * QUERY_CHUNK, b = std::min(n, a + QUERY_CHUNK);
    const char *hinted_idx_end = NULL, *hinted_db_end = NULL;
    uint64_t prev_kmer = ~0ull;
    uint32_t prev_taxon = 0;
    bool prev_found = false;

    for (size_t j = a; j < b; ++j) {
      if (hints && j + HINT_DIST_INDEX < b) {
        const char *p = db->get_index()->prefetch_bin(qs[j + HINT_DIST_INDEX].bin);
        hint_range(p, p + 16, &hinted_idx_end);
      }
      if (hints && j + HINT_DIST_DATA < b) {
        const char *lo, *hi;
        db->bin_bytes(qs[j + HINT_DIST_DATA].bin, &lo, &hi);
        hint_range(lo, hi, &hinted_db_end);
      }
      if (j + PF_DIST_INDEX < b)
        db->get_index()->prefetch_bin(qs[j + PF_DIST_INDEX].bin);
      if (j + PF_DIST_DATA < b) {
        const char *lo, *hi;
        db->bin_bytes(qs[j + PF_DIST_DATA].bin, &lo, &hi);
        __builtin_prefetch(lo);
        __builtin_prefetch(lo + (hi - lo) / 2);
        __builtin_prefetch(hi > lo ? hi - 1 : hi);
      }

      const BatchQuery &q = qs[j];
      uint32_t taxon;
      bool found;
      if (j > a && q.kmer == prev_kmer) {
        taxon = prev_taxon;
        found = prev_found;
      } else {
        uint32_t *val_ptr = db->kmer_query_in_bin(q.kmer, q.bin);
        found = val_ptr != NULL;
        taxon = found ? *val_ptr : 0;
        prev_kmer = q.kmer;
        prev_taxon = taxon;
        prev_found = found;
      }
      if (found) {
        taxa[q.idx] = taxon;
        state[q.idx] = KM_RESOLVED;
      }
    }
  }
}

} // namespace

void process_file_batched(char *filename) {
  string file_str(filename);
  DNASequenceReader *reader;
  Fastq_input = determine_input_file_type(filename);
  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);

  const size_t k = KrakenDatabases[0]->get_k();
  const size_t n_dbs = KrakenDatabases.size();
  vector<size_t> unit_start;

  vector<DNASequence> reads;
  vector<uint64_t> off, cnt;
  vector<uint64_t> kmers;
  vector<uint32_t> taxa;
  vector<uint8_t> state;
  #ifdef _OPENMP
  const int n_threads = omp_get_max_threads();
  #else
  const int n_threads = 1;
  #endif
  vector<BatchQuery> qs;
  vector<vector<BatchQuery> > local(n_threads);  // kept across batches to avoid page-faulting fresh memory

  while (reader->is_valid()) {
    reads.clear();
    unit_start.clear();
    size_t batch_nt = 0;
    // A batch is made of whole work units, delimited like in process_file(),
    // so that per-unit statistics are merged in the same way
    while (batch_nt < Batch_size) {
      const size_t first = reads.size();
      size_t unit_nt = 0;
      while (unit_nt < Work_unit_size) {
        DNASequence dna = reader->next_sequence();
        if (!reader->is_valid())
          break;
        unit_nt += dna.seq.size();
        reads.push_back(std::move(dna));
      }
      if (reads.size() == first)
        break;
      unit_start.push_back(first);
      batch_nt += unit_nt;
    }
    if (reads.empty())
      break;
    unit_start.push_back(reads.size());

    const size_t n_reads = reads.size();
    off.assign(n_reads + 1, 0);
    for (size_t r = 0; r < n_reads; ++r)
      off[r + 1] = off[r] + (reads[r].seq.size() >= k ? reads[r].seq.size() - k + 1 : 0);
    const size_t total = off[n_reads];
    if (total >= (1ull << 32))
      errx(EX_USAGE, "batch too large (reduce -B)");
    cnt.assign(n_reads, 0);
    kmers.assign(total, 0);
    taxa.assign(total, 0);
    state.assign(total, KM_PENDING);

    // Phase 1: k-mer extraction
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (size_t r = 0; r < n_reads; ++r) {
      if (reads[r].seq.size() < k)
        continue;
      KmerScanner scanner(reads[r].seq);
      uint64_t *kmer_ptr;
      size_t i = 0;
      while ((kmer_ptr = scanner.next_kmer()) != NULL) {
        if (scanner.ambig_kmer())
          state[off[r] + i] = KM_AMBIG;
        else
          kmers[off[r] + i] = KrakenDatabases[0]->canonical_representation(*kmer_ptr);
        ++i;
      }
      cnt[r] = i;
    }
    // Slots not produced by the scanner (should not happen) never get queried
    for (size_t r = 0; r < n_reads; ++r)
      for (size_t i = cnt[r]; i < off[r + 1] - off[r]; ++i)
        state[off[r] + i] = KM_AMBIG;

    // Phase 2: look up pending k-mers, one database at a time
    for (size_t d = 0; d < n_dbs; ++d) {
      for (int t = 0; t < n_threads; ++t)
        local[t].clear();
      #ifdef _OPENMP
      #pragma omp parallel
      #endif
      {
        #ifdef _OPENMP
        vector<BatchQuery> &mine = local[omp_get_thread_num()];
        #pragma omp for schedule(static)
        #else
        vector<BatchQuery> &mine = local[0];
        #endif
        for (size_t i = 0; i < total; ++i) {
          if (state[i] == KM_PENDING) {
            BatchQuery q;
            q.kmer = kmers[i];
            q.idx = i;
            q.bin = KrakenDatabases[d]->bin_of(q.kmer);
            mine.push_back(q);
          }
        }
      }
      size_t n_q = 0;
      vector<size_t> starts(n_threads);
      for (int t = 0; t < n_threads; ++t) {
        starts[t] = n_q;
        n_q += local[t].size();
      }
      if (n_q == 0)
        break;
      qs.resize(n_q);
      #ifdef _OPENMP
      #pragma omp parallel for schedule(static, 1)
      #endif
      for (int t = 0; t < n_threads; ++t)
        std::copy(local[t].begin(), local[t].end(), qs.begin() + starts[t]);

      // KRAKEN_BATCH_FORCE_HINTS (testing): take the out-of-cache path even if DBs are resident
      const bool hints = Page_hints &&
          (getenv("KRAKEN_BATCH_FORCE_HINTS") != NULL ||
           resident_fraction(Db_maps[d].first, Db_maps[d].second) < 0.98 ||
           resident_fraction(Idx_maps[d].first, Idx_maps[d].second) < 0.98);
      // Sorting only pays off when the DB has to be paged in; for resident DBs
      // the prefetching in lookup_sorted_queries() is enough
      if (hints) {
        #if defined(_OPENMP) && defined(__GLIBCXX__)
        __gnu_parallel::sort(qs.begin(), qs.end(), batch_query_less);
        #else
        std::sort(qs.begin(), qs.end(), batch_query_less);
        #endif
      }
      lookup_sorted_queries(d, qs, taxa, state, hints);
    }

    // Phase 3: per-read classification, output in input order
    const size_t n_slices = unit_start.size() - 1;
    vector<ostringstream> kos(n_slices), cos_(n_slices), uos(n_slices);
    uint64_t batch_classified = 0;

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
      uint64_t my_classified = 0;
      vector<uint32_t> read_taxa;
      vector<char> ambig_list;
      unordered_map<uint32_t, uint32_t> hit_counts;

      #ifdef _OPENMP
      #pragma omp for schedule(dynamic)
      #endif
      for (size_t s = 0; s < n_slices; ++s) {
        unordered_map<uint32_t, READCOUNTS> unit_counts;
        for (size_t r = unit_start[s]; r < unit_start[s + 1]; ++r) {
          read_taxa.clear();
          ambig_list.clear();
          hit_counts.clear();
          uint32_t taxon = 0;
          for (size_t i = 0; i < cnt[r]; ++i) {
            const size_t x = off[r] + i;
            taxon = 0;
            if (state[x] == KM_AMBIG) {
              ambig_list.push_back(1);
            } else {
              ambig_list.push_back(0);
              taxon = taxa[x];
              unit_counts[taxon].add_kmer(kmers[x]);
              if (taxon)
                hit_counts[taxon]++;
            }
            read_taxa.push_back(taxon);
          }
          my_classified += finish_classification(reads[r], read_taxa, ambig_list, hit_counts,
                                                 0, taxon, kos[s], cos_[s], uos[s],
                                                 unit_counts);
        }
        #ifdef _OPENMP
        #pragma omp critical(write_output)
        #endif
        {
          for (auto it = unit_counts.begin(); it != unit_counts.end(); ++it)
            taxon_counts[it->first] += std::move(it->second);
        }
      }

      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      batch_classified += my_classified;
    }

    for (size_t s = 0; s < n_slices; ++s) {
      if (Print_kraken)
        (*Kraken_output) << kos[s].str();
      if (Print_classified)
        (*Classified_output) << cos_[s].str();
      if (Print_unclassified)
        (*Unclassified_output) << uos[s].str();
    }
    total_classified += batch_classified;
    total_sequences += n_reads;
    total_bases += batch_nt;
    if (Print_Progress)
      fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
              total_sequences, total_classified * 100.0 / total_sequences);
  }

  delete reader;
}

void process_file_with_db_chunk(char *filename) {
  string file_str(filename);
  DNASequence dna;

  Fastq_input = determine_input_file_type(filename);

  std::string dir_for_tmp_file = "./";
  if (!Kraken_output_file.empty())
    dir_for_tmp_file = get_directory(Kraken_output_file);
  else if (!Classified_output_file.empty())
    dir_for_tmp_file = get_directory(Classified_output_file);
  else if (!Unclassified_output_file.empty())
    dir_for_tmp_file = get_directory(Unclassified_output_file);
  else if (!Report_output_file.empty())
    dir_for_tmp_file = get_directory(Report_output_file);
  const std::string tmp_file_name = tempnam(dir_for_tmp_file.c_str(), "tmp");

  // iterate over databases
  bool first_intermediate_output = true;
  for (size_t i=0; i<KrakenDatabases.size(); ++i)
  {
    for (uint32_t db_chunk_id = 0; db_chunk_id < KrakenDatabases[0]->chunks(); ++db_chunk_id)
    {
      total_sequences = 0;
      total_bases = 0;
      KrakenDatabases[i]->load_chunk(db_chunk_id);

      DNASequenceReader *reader;
      if (Fastq_input)
        reader = new FastqReader(file_str);
      else
        reader = new FastaReader(file_str);
      uint32_t seq_idx = 1;

#ifdef _OPENMP
      #pragma omp parallel
#endif
      {
        vector <std::pair<DNASequence, uint32_t> > work_unit;
        const int worker_id = omp_get_thread_num();
        const std::string worker_filename = tmp_file_name + "." + std::to_string(worker_id);

        std::fstream fp(worker_filename, std::fstream::out | std::fstream::binary | std::fstream::trunc);
        fp.exceptions(std::fstream::badbit);

        while (reader->is_valid()) {
          work_unit.clear();
          size_t total_nt = 0;

#ifdef _OPENMP
          #pragma omp critical(get_input)
#endif
          {
            while (total_nt < Work_unit_size) {
              dna = reader->next_sequence();
              if (!reader->is_valid())
                break;
              work_unit.emplace_back(dna, seq_idx);
              total_nt += dna.seq.size();
              ++seq_idx;
            }
          }
          if (total_nt == 0)
            break;

          for (size_t j = 0; j < work_unit.size(); j++) {
            classify_sequence_with_db_chunk(work_unit[j], fp, db_chunk_id, 0);
          }

#ifdef _OPENMP
          #pragma omp critical(progress)
#endif
          {
            total_sequences += work_unit.size();
            total_bases += total_nt;
            if (Print_Progress) {
              fprintf(stderr, "\r Processed %llu sequences (database chunk %" PRIu32" of %" PRIu32 ")",
                      total_sequences, db_chunk_id + 1, KrakenDatabases[i]->chunks());
            }
          }
        }
        fp.close();
      }  // end parallel section

      delete reader;
      merge_intermediate_results_by_workers(first_intermediate_output, tmp_file_name);
      first_intermediate_output = false;
    }
  }

  fprintf(stderr, "\r Processed %llu sequences\n", total_sequences);

  // merge intermediate results of all chunks and classify
  // TODO: parallelize this (need to buffer reads), we also do not care about the final output order
  total_sequences = 0;
  total_classified = 0;

  DNASequenceReader *reader;
  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);

  unordered_map<uint32_t, uint32_t> hit_counts;
  vector<uint32_t> taxa;
  vector<char> ambig_list;

  const std::string taxa_summary_filename = tmp_file_name;
  FILE* fp_taxa_summary = fopen(taxa_summary_filename.c_str(), "rb");

  while (reader->is_valid()) {
    hit_counts.clear();
    taxa.clear();
    ambig_list.clear();

    dna = reader->next_sequence();
    if (!reader->is_valid())
      break;

    // merge results from chunks into 'taxa'
    uint32_t hits = 0;

    // get number of elements (taxa.size())
    uint32_t taxa_size;
    fread(&taxa_size, sizeof(taxa_size), 1, fp_taxa_summary);

    taxa.resize(taxa_size);
    fread(&taxa[0], sizeof(uint32_t), taxa_size, fp_taxa_summary);

    for (uint32_t i = 0; i < taxa.size(); ++i)
    {
      if (taxa[i])
      {
        ++hit_counts[taxa[i]];
        // TODO: stop querying this read with other databases
        if (Quick_mode && ++hits >= Minimum_hit_count)
          goto quick_mode_call;
      }
    }

quick_mode_call:
    uint64_t *kmer_ptr;
    uint32_t taxon = 0;
    if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
      KmerScanner scanner(dna.seq);
      uint32_t taxa_idx = 0;
      while ((kmer_ptr = scanner.next_kmer()) != NULL) {
        if (scanner.ambig_kmer()) {
          ambig_list.push_back(1);
        } else {
          ambig_list.push_back(0);
          uint64_t cannonical_kmer = KrakenDatabases[0]->canonical_representation(*kmer_ptr);
          taxon = taxa[taxa_idx];
          taxon_counts[taxon].add_kmer(cannonical_kmer);
        }
        ++taxa_idx;
      }
    }

    uint32_t call = 0;
    if (Map_UIDs) {
      if (Quick_mode) {
        cerr << "Quick mode not available when mapping UIDs" << endl;
        exit(1);
      } else {
        call = resolve_uids3(hit_counts, Parent_map, Uid_dict,
                             UID_to_TaxID_map_file.ptr(), UID_to_TaxID_map_file.size());
      }
    } else {
      if (Quick_mode)
        call = hits >= Minimum_hit_count ? taxon : 0;
      else
        call = resolve_tree(hit_counts, Parent_map);
    }

    total_classified += (call != 0);
    taxon_counts[call].incrementReadCount();

    if (Print_unclassified && !call)
      print_sequence(Unclassified_output, dna);

    if (Print_classified && call)
      print_sequence(Classified_output, dna);

    if (!Print_kraken) {
      goto next_read; // return call;
    }

    if (call) {
      (*Kraken_output) << "C\t";
    }
    else {
      if (Only_classified_kraken_output) {
        goto next_read; // return false;
      }
      (*Kraken_output) << "U\t";
    }
    (*Kraken_output) << dna.id << '\t' << call << '\t' << dna.seq.size() << '\t';

    if (Quick_mode) {
      (*Kraken_output) << "Q:" << hits;
    }
    else {
      if (taxa.empty())
        (*Kraken_output) << "0:0";
      else
        (*Kraken_output) << hitlist_string(taxa, ambig_list);
    }

    if (Print_sequence)
      (*Kraken_output) << "\t" << dna.seq;

    (*Kraken_output) << "\n";
    // return call;
next_read:
    ++total_sequences;
    fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
            total_sequences, total_classified * 100.0 / total_sequences);
    (void)1;
  }

  fclose(fp_taxa_summary);
  unlink(taxa_summary_filename.c_str());

  delete reader;
}


inline void print_sequence(ostream* oss_ptr, const DNASequence& dna) {
      if (Fastq_input) {
        (*oss_ptr) << "@" << dna.header_line << endl
            << dna.seq << endl
            << "+" << endl
            << dna.quals << endl;
      }
      else {
        (*oss_ptr) << ">" << dna.header_line << endl
            << dna.seq << endl;
      }
}

/*
inline
void append_hitlist_string(string& hitlist_string, uint32_t& last_taxon, uint32_t& last_counter, uint32_t current_taxon) {
  if (last_taxon == current_taxon) {
    ++last_counter;
  } else {
    if (last_counter > 0) {
      if (last_taxon == ambig_taxon) {
        hitlist_string += "A:" + std::to_string(last_counter) + ' ';
      } else {
        hitlist_string += std::to_string(last_taxon) + ':' + std::to_string(last_counter) + ' ';
      }
    }
    last_counter = 1;
    last_taxon = current_taxon;
  }
}
*/

string hitlist_string(const vector<uint32_t> &taxa, const vector<char> &ambig)
{
  int64_t last_code;
  int code_count = 1;
  ostringstream hitlist;

  if (ambig[0])   { last_code = -1; }
  else            { last_code = taxa[0]; }

  for (size_t i = 1; i < taxa.size(); i++) {
    int64_t code;
    if (ambig[i]) { code = -1; }
    else          { code = taxa[i]; }

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code >= 0) {
    hitlist << last_code << ":" << code_count;
  }
  else {
    hitlist << "A:" << code_count;
  }
  return hitlist.str();
}

/*
string hitlist_string_depr(const vector<uint32_t> &taxa)
{
  uint32_t last_code = taxa[0];
  int code_count = 1;
  ostringstream hitlist;

  for (size_t i = 1; i < taxa.size(); i++) {
    uint32_t code = taxa[i];

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code == -1) {
    hitlist << "A:" << code_count;
  }
  else {
    hitlist << last_code << ":" << code_count;
  }
  return hitlist.str();
}
*/

bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, READCOUNTS>& my_taxon_counts) {
  vector<uint32_t> taxa;
  vector<char> ambig_list;
  unordered_map<uint32_t, uint32_t> hit_counts;
  uint64_t *kmer_ptr;
  uint32_t taxon = 0;
  uint32_t hits = 0;  // only maintained if in quick mode

  //string hitlist_string;
  //uint32_t last_taxon;
  //uint32_t last_counter;

  vector<db_status> db_statuses(KrakenDatabases.size());

  if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
    size_t n_kmers = dna.seq.size()-KrakenDatabases[0]->get_k()+1;
    taxa.reserve(n_kmers);
    ambig_list.reserve(n_kmers);
    KmerScanner scanner(dna.seq);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (scanner.ambig_kmer()) {
        //append_hitlist_string(hitlist_string, last_taxon, last_counter, ambig_taxon);
        ambig_list.push_back(1);
      }
      else {
        uint64_t cannonical_kmer = KrakenDatabases[0]->canonical_representation(*kmer_ptr);
        ambig_list.push_back(0);
        // go through multiple databases to map k-mer
        for (size_t i=0; i<KrakenDatabases.size(); ++i) {
          uint32_t* val_ptr = KrakenDatabases[i]->kmer_query(
            cannonical_kmer, &db_statuses[i].current_bin_key,
            &db_statuses[i].current_min_pos, &db_statuses[i].current_max_pos);
          if (val_ptr) {
            taxon = *val_ptr;
            break;
          }
        }

        // cerr << "taxon for " << *kmer_ptr << " is " << taxon << endl;
        my_taxon_counts[taxon].add_kmer(cannonical_kmer);

        if (taxon) {
          hit_counts[taxon]++;
          if (Quick_mode && ++hits >= Minimum_hit_count)
            break;
        }
      }
      taxa.push_back(taxon);
      //append_hitlist_string(hitlist_string, last_taxon, last_counter, taxon);
    }
  }

  return finish_classification(dna, taxa, ambig_list, hit_counts, hits, taxon,
                               koss, coss, uoss, my_taxon_counts);
}

// Everything that follows the k-mer lookups: calls the read and writes output
bool finish_classification(DNASequence &dna, const vector<uint32_t> &taxa,
                           const vector<char> &ambig_list,
                           unordered_map<uint32_t, uint32_t> &hit_counts,
                           uint32_t hits, uint32_t taxon, ostringstream &koss,
                           ostringstream &coss, ostringstream &uoss,
                           unordered_map<uint32_t, READCOUNTS>& my_taxon_counts) {
  uint32_t call = 0;
  if (Map_UIDs) {
    if (Quick_mode) {
      cerr << "Quick mode not available when mapping UIDs" << endl;
      exit(1);
    } else {
      call = resolve_uids3(hit_counts, Parent_map, Uid_dict,
        UID_to_TaxID_map_file.ptr(), UID_to_TaxID_map_file.size());
    }
  } else {
    if (Quick_mode)
      call = hits >= Minimum_hit_count ? taxon : 0;
    else
      call = resolve_tree(hit_counts, Parent_map);
  }

  my_taxon_counts[call].incrementReadCount();

  if (Print_unclassified && !call) 
    print_sequence(&uoss, dna);

  if (Print_classified && call)
    print_sequence(&coss, dna);


  if (! Print_kraken)
    return call;

  if (call) {
    koss << "C\t";
  }
  else {
    if (Only_classified_kraken_output)
      return false;
    koss << "U\t";
  }
  koss << dna.id << '\t' << call << '\t' << dna.seq.size() << '\t';

  if (Quick_mode) {
    koss << "Q:" << hits;
  }
  else {
    if (taxa.empty())
      koss << "0:0";
    else
      koss << hitlist_string(taxa, ambig_list);
    //if (hitlist_string.empty() && last_counter == 0)
    //  koss << "0:0";
    //else {
    //  koss << hitlist_string
    //       << (last_taxon == ambig_taxon? "A" :  std::to_string(last_taxon))
    //       << ':' << std::to_string(last_counter);
    //}
  }

  if (Print_sequence)
      koss << "\t" << dna.seq;

  koss << "\n";
  return call;
}

void classify_sequence_with_db_chunk(std::pair<DNASequence, uint32_t> & seq, std::fstream & fp, const uint32_t db_chunk_id, const uint32_t db_id) {
  vector<uint32_t> taxa;
  uint64_t *kmer_ptr;
  uint32_t taxon;

  auto & dna = seq.first;
  const uint32_t & seq_idx = seq.second;

  vector<db_status> db_statuses(KrakenDatabases.size());

  if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
    size_t n_kmers = dna.seq.size()-KrakenDatabases[0]->get_k()+1;
    taxa.reserve(n_kmers);
    KmerScanner scanner(dna.seq);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (!scanner.ambig_kmer()) {
        uint64_t cannonical_kmer = KrakenDatabases[db_id]->canonical_representation(*kmer_ptr);
        const uint64_t minimizer = KrakenDatabases[db_id]->bin_key(cannonical_kmer); // TODO: inefficient because minimizer will also be computed in kmer_query()

        if (KrakenDatabases[db_id]->is_minimizer_in_chunk(minimizer, db_chunk_id)) {
          uint32_t* val_ptr = KrakenDatabases[db_id]->kmer_query_with_db_chunks(
                  cannonical_kmer, &db_statuses[db_id].current_bin_key,
                  &db_statuses[db_id].current_min_pos, &db_statuses[db_id].current_max_pos);
          if (val_ptr)
            taxon = *val_ptr;
        }
      }
      taxa.push_back(taxon);
    }
  }

  try {
    const uint32_t taxa_size = taxa.size();
    fp.write((char*) &seq_idx, 1 * sizeof(uint32_t)); // seq_idx
    fp.write((char*) &taxa_size, 1 * sizeof(uint32_t)); // number of elements
    fp.write((char*) &taxa[0], taxa_size * sizeof(uint32_t)); // elements
  }
  catch (const std::fstream::failure& e) {
    printf("ERROR: Could not write to temporary file: %s\n", e.what());
    exit(1);
  }
}

set<uint32_t> get_ancestry(uint32_t taxon) {
  set<uint32_t> path;

  while (taxon > 0) {
    path.insert(taxon);
    taxon = Parent_map[taxon];
  }
  return path;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qcC:U:Ma:r:sI:p:x:B:H")) != -1) {
    switch (opt) {
      case 'd' :
        DB_filenames.push_back(optarg);
        break;
      case 'i' :
        Index_filenames.push_back(optarg);
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
      case 'p' :
        HLL_PRECISION = stoi(optarg);
        break;
      case 'q' :
        Quick_mode = true;
        break;
      case 'm' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive minimum hit count");
        Minimum_hit_count = sig;
        break;
      case 'c' :
        Only_classified_kraken_output = true;
        break;
      case 'C' :
        Print_classified = true;
        Classified_output_file = optarg;
        break;
      case 'U' :
        Print_unclassified = true;
        Unclassified_output_file = optarg;
        break;
      case 'o' :
        Kraken_output_file = optarg;
        break;
      case 'r' :
        Report_output_file = optarg;
        break;
      case 's' :
        Print_sequence = true;
        break;
      case 'a' :
        TaxDB_file = optarg;
        break;
      case 'u' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive work unit size");
        Work_unit_size = sig;
        break;
      case 'M' :
        Populate_memory = true;
        break;
      case 'x' :
        Populate_memory = true;
        Populate_memory_size = parse_human_readable_size(optarg); // strtoull(optarg, NULL, 0);
        break;
      case 'B' :
        {
        const string size_arg(optarg);
        const size_t digits = size_arg.find_first_not_of("0123456789");
        if (size_arg.empty() || digits == 0 ||
            (digits != string::npos &&
             (digits + 1 != size_arg.size() ||
              string("kKmMgG").find(size_arg[digits]) == string::npos)))
          errx(EX_USAGE, "invalid batch size: %s", optarg);
        Batch_size = parse_human_readable_size(optarg);
        if (!Batch_size && size_arg.find_first_of("123456789") != string::npos)
          errx(EX_USAGE, "invalid or overflowing batch size: %s", optarg);
        }
        break;
      case 'H' :
        Page_hints = false;
        break;
      case 'I' :
        UID_to_TaxID_map_filename = optarg;
        Map_UIDs = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filenames.empty()) {
    cerr << "Missing mandatory option -d" << endl;
    usage();
  }
  if (Index_filenames.empty()) {
    cerr << "Missing mandatory option -i" << endl;
    usage();
  }
  if (optind == argc && !Populate_memory) {
    cerr << "No sequence data files specified" << endl;
  }
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -r filename      Output file for Kraken report output" << endl
       << "  -a filename      TaxDB" << endl
       << "  -I filename      UID to TaxId map" << endl
       << "  -p #             Precision for unique k-mer counting, between 10 and 18" << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -x size          Preload database files using x amount of RAM (e.g. 10G)" << endl
       << "  -B size          Look up k-mers in batches of reads totalling this many bp (default 16M, 0 = read by read)" << endl
       << "  -H               Don't sort k-mers and issue paging hints for out-of-cache databases (batch mode)" << endl
       << "  -s               Print read sequence in Kraken output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}
