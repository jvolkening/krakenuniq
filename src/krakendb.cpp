/*
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
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
#include "krakendb.hpp"
#include "quickfile.hpp"
#include <unordered_map>

using std::string;
using std::vector;

namespace kraken {

// File type code for Jellyfish/Kraken DBs
static const char * DATABASE_FILE_TYPE = "JFLISTDN";

// File type code on Kraken DB index
// Next byte determines # of indexed nt
static const char * KRAKEN_INDEX_STRING = "KRAKIDX";

// File type code for Kraken DB index (v2)
// v2 index corresponds to DB sorted by scrambled order
// Next byte determines # of indexed nt
static const char * KRAKEN_INDEX2_STRING = "KRAKIX2";

// (INDEX2_XOR_MASK, the XOR mask for minimizer bin keys, is in compactdb.hpp)

// Basic constructor
KrakenDB::KrakenDB() {
  fptr = NULL;
  index_ptr = NULL;
  key_ct = 0;
  val_len = 0;
  key_len = 0;
  key_bits = 0;
  k = 0;
  _filesize = 0;
  _chunks = 0;
  data = NULL;
  data_size = 0;
  data_offset = 0;
  compact = false;
  rec_bits = 0;
  dict_bits = 0;
  dict_size = 0;
  dict = NULL;
  cdata = NULL;
  cdata_bytes = 0;
}

// Assumes ptr points to start of a readable mmap'ed file
KrakenDB::KrakenDB(char *ptr, size_t filesize) {
  _filesize = filesize;
  index_ptr = NULL;
  fptr = ptr;
  _chunks = 0;
  data = NULL;
  data_size = 0;
  data_offset = 0;
  compact = false;
  rec_bits = 0;
  dict_bits = 0;
  dict_size = 0;
  dict = NULL;
  cdata = NULL;
  cdata_bytes = 0;
  if (ptr == NULL) {
    errx(EX_DATAERR, "pointer is NULL");
  }
  if (!strncmp(ptr, COMPACT_DB_MAGIC, strlen(COMPACT_DB_MAGIC))) {
    uint64_t hdr[11];
    memcpy(hdr, ptr, sizeof(hdr));
    compact = true;
    k = hdr[1];
    uint64_t nt = hdr[2];
    key_ct = hdr[3];
    dict_size = hdr[4];
    dict_bits = hdr[5];
    uint64_t pos_bits = hdr[6];
    rec_bits = hdr[7];
    uint64_t max_taxid_bits = hdr[8];
    cdata_bytes = hdr[9];
    if (k < 1 || k > 31 || nt < 1 || nt > k)
      errx(EX_DATAERR, "compact database has invalid k (%llu) or minimizer length (%llu)",
           (unsigned long long) k, (unsigned long long) nt);
    if (max_taxid_bits > COMPACT_MAX_TAXID_BITS)
      errx(EX_DATAERR, "compact database allows %llu-bit taxids, reader supports at most %u",
           (unsigned long long) max_taxid_bits, COMPACT_MAX_TAXID_BITS);
    key_bits = 2 * k;
    val_len = 4;
    key_len = key_bits / 8 + !! (key_bits % 8);
    codec.init(k, nt, INDEX2_XOR_MASK);
    if (codec.pos_bits != pos_bits || codec.t_bits + dict_bits != rec_bits || codec.t_bits > 64
        || dict_bits > 32 || dict_size == 0)
      errx(EX_DATAERR, "compact database header is inconsistent");
    uint64_t dict_bytes = (dict_size * 4 + 7) & ~7ull;
    if (filesize != 0 && COMPACT_DB_HEADER_BYTES + dict_bytes + cdata_bytes + COMPACT_PAD_BYTES > filesize)
      errx(EX_DATAERR, "compact database is truncated");
    if (cdata_bytes != (key_ct * rec_bits + 7) / 8)
      errx(EX_DATAERR, "compact database size does not match its key count");
    dict = (const uint32_t *) (ptr + COMPACT_DB_HEADER_BYTES);
    cdata = ptr + COMPACT_DB_HEADER_BYTES + dict_bytes;
    for (uint64_t i = 0; i < dict_size; i++) {
      if (dict[i] >> max_taxid_bits)
        errx(EX_DATAERR, "compact database contains taxid %u which needs more than %llu bits",
             dict[i], (unsigned long long) max_taxid_bits);
    }
    std::cerr << "Parsed compact database header with " << key_ct << " keys with k of " << (size_t)k
              << " [minimizer length " << nt << ", " << dict_size << " distinct taxids, "
              << rec_bits << " bits per key]." << std::endl;
    return;
  }
  if (strncmp(ptr, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE))) {
    errx(EX_DATAERR,"database in improper format - found %s", string(ptr, strlen(DATABASE_FILE_TYPE)).c_str());
  }
  memcpy(&key_bits, ptr + 8, 8);
  memcpy(&val_len, ptr + 16, 8);
  memcpy(&key_ct, ptr + 48, 8);
  if (val_len != 4)
    errx(EX_DATAERR, "can only handle 4 byte DB values");
  k = key_bits / 2;
  key_len = key_bits / 8 + !! (key_bits % 8);
  std::cerr << "Parsed database header with " << key_ct << " keys with k of " << (size_t)k << " [val_len " << val_len << ", key_len " << key_len << "]." << std::endl;
}

// destructor
KrakenDB::~KrakenDB() {
  if (data != NULL && data_size > 0)
    munmap(data, data_size);
}

size_t KrakenDB::filesize() const {
  return _filesize;
}

//using std::map to have the keys sorted
std::map<uint32_t,uint64_t> KrakenDB::count_taxons() {
  if (compact)
    return count_taxons_compact();
  char *ptr = get_pair_ptr();
  size_t pair_sz = pair_size();

  std::map<uint32_t, uint64_t> taxon_counts;
  if (ptr == NULL) { 
    std::cerr << "Kraken database pointer is NULL [pair_sz: " << pair_sz << ", key_ct: "<<key_ct<<", key_len: "<< key_len<<"]!" << std::endl;
    exit(1);
  }
  for (uint64_t i = 0; i < key_ct; i++) {
    if (i % 10000000 == 1) {
      fprintf(stderr, "\r %.2f %%", 100*(double)i/(double)key_ct );
    }
    uint32_t* taxon = (uint32_t *) (ptr + pair_sz * i + key_len);
    if (taxon == NULL) {
        std::cerr << "taxon is NULL (i is " << i << " and key_ct is " << key_ct << ")" << std::endl;
    } else {
        uint32_t taxon_i = *taxon;
        ++taxon_counts[taxon_i];
    }
  }
  fprintf(stderr, "\r");
  return taxon_counts;
}


std::map<uint32_t,uint64_t> KrakenDB::count_taxons_compact() {
  std::vector<uint64_t> counts(dict_size, 0);
  for (uint64_t i = 0; i < key_ct; i++) {
    if (i % 10000000 == 1)
      fprintf(stderr, "\r %.2f %%", 100*(double)i/(double)key_ct );
    uint64_t di = compact_read_bits(cdata, i * rec_bits + codec.t_bits, dict_bits);
    if (di >= dict_size)
      errx(EX_DATAERR, "compact database is corrupt (taxid index out of range)");
    counts[di]++;
  }
  fprintf(stderr, "\r");
  std::map<uint32_t, uint64_t> taxon_counts;
  for (uint64_t i = 0; i < dict_size; i++)
    if (counts[i])
      taxon_counts[dict[i]] = counts[i];
  return taxon_counts;
}

// Creates an index, indicating starting positions of each bin
// Bins contain k-mer/taxon pairs with k-mers that share a bin key
void KrakenDB::make_index(string index_filename, uint8_t nt) {
  uint64_t entries = 1ull << (nt * 2);
  vector<uint64_t> bin_counts(entries);
  char *ptr = get_pair_ptr();

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic,400)
#endif
  for (uint64_t i = 0; i < key_ct; i++) {
    uint64_t kmer = 0;
    memcpy(&kmer, ptr + i * pair_size(), key_len);
    uint64_t b_key = bin_key(kmer, nt);
#ifdef _OPENMP
    #pragma omp atomic
#endif
    bin_counts[b_key]++;
  }

  uint64_t *bin_offsets = new uint64_t[ entries + 1 ];
  bin_offsets[0] = 0;
  for (uint64_t i = 1; i <= entries; i++)
    bin_offsets[i] = bin_offsets[i-1] + bin_counts[i-1];

  QuickFile idx_file(index_filename, "w",
    strlen(KRAKEN_INDEX2_STRING) + 1 + sizeof(*bin_offsets) * (entries + 1));
  char *idx_ptr = idx_file.ptr();
  memcpy(idx_ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING));
  idx_ptr += strlen(KRAKEN_INDEX2_STRING);
  memcpy(idx_ptr++, &nt, 1);
  memcpy(idx_ptr, bin_offsets, sizeof(*bin_offsets) * (entries + 1));
}

// Simple accessor
char *KrakenDB::get_ptr() {
  return fptr;
}

// Returns start of k-mer/taxon pair array (skips header)
char *KrakenDB::get_pair_ptr() {
  return fptr == NULL ? NULL : fptr + header_size();
}

// Simple accessor
KrakenDBIndex *KrakenDB::get_index() {
  return index_ptr;
}

// Associates the index with this database
void KrakenDB::set_index(KrakenDBIndex *i_ptr) {
  if (i_ptr != NULL && compact != i_ptr->is_compact())
    errx(EX_DATAERR, "compact databases need compact indexes (and vice versa)");
  if (i_ptr != NULL && compact && i_ptr->indexed_nt() != codec.nt)
    errx(EX_DATAERR, "index minimizer length does not match compact database");
  if (i_ptr != NULL && compact && i_ptr->ef_max != key_ct)
    errx(EX_DATAERR, "compact index does not match compact database (key count differs)");
  index_ptr = i_ptr;
}

// Simple accessors/convenience methods
uint8_t KrakenDB::get_k() { return k; }
uint64_t KrakenDB::get_key_bits() { return key_bits; }
uint64_t KrakenDB::get_key_len() { return key_len; }
uint64_t KrakenDB::get_val_len() { return val_len; }
uint64_t KrakenDB::get_key_ct() { return key_ct; }
uint64_t KrakenDB::pair_size() { return key_len + val_len; }
size_t KrakenDB::header_size() { return 72 + 2 * (4 + 8 * key_bits); }

// Bin key: each k-mer is made of several overlapping m-mers, m < k
// The bin key is the m-mer whose canonical representation is "smallest"
// ("smallest" can refer to lexico. ordering or some other ordering)
uint64_t KrakenDB::bin_key(uint64_t kmer, uint64_t idx_nt) {
  return min_bin_key(kmer, idx_nt, INDEX2_XOR_MASK);
}

// The reverse complement of each m-mer window is a window of the reverse
// complement of the whole k-mer, so it is computed only once per k-mer.
uint64_t KrakenDB::min_bin_key(uint64_t kmer, uint8_t nt, uint64_t xor_mask) {
  const unsigned kk = key_bits / 2;
  const uint64_t mask = (1ull << (nt * 2)) - 1;
  xor_mask &= mask;
  const uint64_t revcom = reverse_complement(kmer, kk);
  uint64_t min_bin_key = ~0;
  for (unsigned i = 0; i < kk - nt + 1; i++) {
    const uint64_t fwd = (kmer >> (2 * i)) & mask;
    const uint64_t rev = (revcom >> (2 * (kk - nt - i))) & mask;
    const uint64_t temp_bin_key = xor_mask ^ (fwd < rev ? fwd : rev);
    if (temp_bin_key < min_bin_key)
      min_bin_key = temp_bin_key;
  }
  return min_bin_key;
}

// Separate functions to avoid a conditional in the function
// This probably isn't necessary...
uint64_t KrakenDB::bin_key(uint64_t kmer) {
  uint8_t nt = index_ptr->indexed_nt();
  uint8_t idx_type = index_ptr->index_type();
  return min_bin_key(kmer, nt, idx_type == 1 ? 0 : INDEX2_XOR_MASK);
}

// Code mostly from Jellyfish 1.6 source
uint64_t KrakenDB::reverse_complement(uint64_t kmer, uint8_t n) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

// Code mostly from Jellyfish 1.6 source
uint64_t KrakenDB::reverse_complement(uint64_t kmer) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));
}

// Lexicographically smallest of k-mer and reverse comp. of k-mer
uint64_t KrakenDB::canonical_representation(uint64_t kmer, uint8_t n) {
  uint64_t revcom = reverse_complement(kmer, n);
  return kmer < revcom ? kmer : revcom;
}

uint64_t KrakenDB::canonical_representation(uint64_t kmer) {
  uint64_t revcom = reverse_complement(kmer, k);
  return kmer < revcom ? kmer : revcom;
}

// perform search over last range to speed up queries
// NOTE: retry_on_failure implies all pointer params are non-NULL
uint32_t *KrakenDB::kmer_query(uint64_t kmer, uint64_t *last_bin_key,
                               int64_t *min_pos, int64_t *max_pos,
                               bool retry_on_failure)
{
  if (compact)
    return compact_query(kmer, last_bin_key, min_pos, max_pos, retry_on_failure, false);
  int64_t min, max, mid;
  uint64_t comp_kmer;
  uint64_t b_key;
  char *ptr = get_pair_ptr();
  size_t pair_sz = pair_size();

  // Use provided values if they exist and are valid
  if (retry_on_failure && *min_pos <= *max_pos) {
    b_key = *last_bin_key;
    min = *min_pos;
    max = *max_pos;
  }
  else {
    b_key = bin_key(kmer);
    min = index_ptr->at(b_key);
    max = index_ptr->at(b_key + 1) - 1;
    // Invalid min/max values + retry_on_failure means min/max need to be
    // initialized and set in caller
    if (retry_on_failure) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }

  // Binary search with large window
  while (min + 15 <= max) {
    mid = min + (max - min) / 2;
    comp_kmer = 0;
    memcpy(&comp_kmer, ptr + pair_sz * mid, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer > comp_kmer)
      min = mid + 1;
    else if (kmer < comp_kmer)
      max = mid - 1;
    else
      return (uint32_t *) (ptr + pair_sz * mid + key_len);
  }
  // Linear search once window shrinks
  for (mid = min; mid <= max; mid++) {
    comp_kmer = 0;
    memcpy(&comp_kmer, ptr + pair_sz * mid, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer == comp_kmer)
      return (uint32_t *) (ptr + pair_sz * mid + key_len);
  }

  uint32_t *answer = NULL;
  // ROF implies the provided values might be out of date
  // If they are, we'll update them and search again
  if (retry_on_failure) {
    b_key = bin_key(kmer);
    // If bin key hasn't changed, search fails
    if (b_key == *last_bin_key)
      return NULL;
    min = index_ptr->at(b_key);
    max = index_ptr->at(b_key + 1) - 1;
    // Recursive call w/ adjusted search params and w/o retry
    answer = kmer_query(kmer, &b_key, &min, &max, false);
    // Update caller's search params due to bin key change
    if (last_bin_key != NULL) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }
  return answer;
}

// Lookup in a compact database. The record key is t (see CompactKmerCodec),
// which is derived from the query k-mer; the bin's records are sorted by t.
uint32_t *KrakenDB::compact_query(uint64_t kmer, uint64_t *last_bin_key,
                                  int64_t *min_pos, int64_t *max_pos,
                                  bool retry_on_failure, bool use_chunks)
{
  uint64_t b_key, t;
  codec.encode(kmer, &b_key, &t);

  int64_t min, max;
  if (retry_on_failure && *last_bin_key == b_key && *min_pos <= *max_pos) {
    min = *min_pos;
    max = *max_pos;
  }
  else {
    uint64_t start, end;
    index_ptr->bin_range(b_key, &start, &end);
    min = start;
    max = (int64_t) end - 1;
    if (retry_on_failure) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }

  const char *base = use_chunks ? data : cdata;
  const uint64_t bit_shift = use_chunks ? data_offset * 8 : 0;
  const unsigned t_bits = codec.t_bits;
  const uint64_t R = rec_bits;

  while (min + 7 <= max) {
    int64_t mid = min + (max - min) / 2;
    uint64_t comp = compact_read_bits(base, mid * R - bit_shift, t_bits);
    if (t > comp)
      min = mid + 1;
    else if (t < comp)
      max = mid - 1;
    else {
      min = max = mid;
      break;
    }
  }
  for (int64_t mid = min; mid <= max; mid++) {
    uint64_t pos = mid * R - bit_shift;
    if (compact_read_bits(base, pos, t_bits) == t) {
      uint64_t di = compact_read_bits(base, pos + t_bits, dict_bits);
      return di < dict_size ? (uint32_t *) &dict[di] : NULL;
    }
  }
  return NULL;
}

uint64_t KrakenDB::bin_of(uint64_t kmer) {
  if (compact) {
    uint64_t b_key, t;
    codec.encode(kmer, &b_key, &t);
    return b_key;
  }
  return bin_key(kmer);
}

void KrakenDB::bin_bytes(uint64_t bin, const char **lo, const char **hi) {
  uint64_t start, end;
  index_ptr->bin_range(bin, &start, &end);
  if (compact) {
    *lo = cdata + ((start * rec_bits) >> 3);
    *hi = cdata + ((end * rec_bits + 7) >> 3) + 8;
  } else {
    *lo = get_pair_ptr() + start * pair_size();
    *hi = get_pair_ptr() + end * pair_size();
  }
}

uint32_t *KrakenDB::kmer_query_in_bin(uint64_t kmer, uint64_t bin) {
  if (compact) {
    uint64_t b_key = bin;
    int64_t min = 0, max = -1;
    return compact_query(kmer, &b_key, &min, &max, false, false);
  }
  int64_t min = index_ptr->at(bin);
  int64_t max = index_ptr->at(bin + 1) - 1;
  const char *ptr = get_pair_ptr();
  const size_t pair_sz = pair_size();
  const uint64_t key_mask = (1ull << key_bits) - 1;
  uint64_t comp_kmer;
  while (min + 15 <= max) {
    int64_t mid = min + (max - min) / 2;
    comp_kmer = 0;
    memcpy(&comp_kmer, ptr + pair_sz * mid, key_len);
    comp_kmer &= key_mask;
    if (kmer > comp_kmer)
      min = mid + 1;
    else if (kmer < comp_kmer)
      max = mid - 1;
    else
      return (uint32_t *) (ptr + pair_sz * mid + key_len);
  }
  for (int64_t mid = min; mid <= max; mid++) {
    comp_kmer = 0;
    memcpy(&comp_kmer, ptr + pair_sz * mid, key_len);
    comp_kmer &= key_mask;
    if (kmer == comp_kmer)
      return (uint32_t *) (ptr + pair_sz * mid + key_len);
  }
  return NULL;
}

// Binary search w/in the k-mer's bin
uint32_t *KrakenDB::kmer_query(uint64_t kmer) {
  return kmer_query(kmer, NULL, NULL, NULL, false);
}

// perform search over last range to speed up queries
// NOTE: retry_on_failure implies all pointer params are non-NULL
uint32_t *KrakenDB::kmer_query_with_db_chunks(uint64_t kmer, uint64_t *last_bin_key,
                               int64_t *min_pos, int64_t *max_pos,
                               bool retry_on_failure)
{
  if (compact)
    return compact_query(kmer, last_bin_key, min_pos, max_pos, retry_on_failure, true);
  int64_t min, max, mid;
  uint64_t comp_kmer;
  uint64_t b_key;
  size_t pair_sz = pair_size();

  // Use provided values if they exist and are valid
  if (retry_on_failure && *min_pos <= *max_pos) {
    b_key = *last_bin_key;
    min = *min_pos;
    max = *max_pos;
  }
  else {
    b_key = bin_key(kmer);
    min = index_ptr->at_with_db_chunks(b_key);
    max = index_ptr->at_with_db_chunks(b_key + 1) - 1;
    // Invalid min/max values + retry_on_failure means min/max need to be
    // initialized and set in caller
    if (retry_on_failure) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }

  // Binary search with large window
  while (min + 15 <= max) {
    mid = min + (max - min) / 2;
    comp_kmer = 0;
    memcpy(&comp_kmer, data + pair_sz * mid - data_offset, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer > comp_kmer)
      min = mid + 1;
    else if (kmer < comp_kmer)
      max = mid - 1;
    else
      return (uint32_t *) (data + pair_sz * mid - data_offset + key_len);
  }
  // Linear search once window shrinks
  for (mid = min; mid <= max; mid++) {
    comp_kmer = 0;
    memcpy(&comp_kmer, data + pair_sz * mid - data_offset, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer == comp_kmer)
      return (uint32_t *) (data + pair_sz * mid - data_offset + key_len);
  }

  uint32_t *answer = NULL;
  // ROF implies the provided values might be out of date
  // If they are, we'll update them and search again
  if (retry_on_failure) {
    b_key = bin_key(kmer);
    // If bin key hasn't changed, search fails
    if (b_key == *last_bin_key)
      return NULL;
    min = index_ptr->at_with_db_chunks(b_key);
    max = index_ptr->at_with_db_chunks(b_key + 1) - 1;
    // Recursive call w/ adjusted search params and w/o retry
    answer = kmer_query_with_db_chunks(kmer, &b_key, &min, &max, false);
    // Update caller's search params due to bin key change
    if (last_bin_key != NULL) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }
  return answer;
}

// Binary search w/in the k-mer's bin
uint32_t *KrakenDB::kmer_query_with_db_chunks(uint64_t kmer) {
  return kmer_query_with_db_chunks(kmer, NULL, NULL, NULL, false);
}

uint32_t KrakenDB::chunks() const {
  return this->_chunks;
}

void KrakenDB::load_chunk(const uint32_t db_chunk_id) {
  if (compact) {
    // Only the bit-packed records are chunked; the compact index stays mmap'ed
    const uint64_t start_byte = (dbx_chunk_bounds[db_chunk_id] * rec_bits) >> 3;
    const uint64_t end_byte = (dbx_chunk_bounds[db_chunk_id + 1] * rec_bits + 7) >> 3;
    memcpy(data, cdata + start_byte, end_byte - start_byte);
    memset(data + (end_byte - start_byte), 0, COMPACT_PAD_BYTES);
    data_offset = start_byte;
    return;
  }
  char* db_chunk_start = get_pair_ptr() + dbx_chunk_bounds[db_chunk_id] * pair_size();
  data_offset = dbx_chunk_bounds[db_chunk_id] * pair_size();
  const size_t db_chunk_len = (dbx_chunk_bounds[db_chunk_id + 1] - dbx_chunk_bounds[db_chunk_id]) * pair_size();
  memcpy(data, db_chunk_start, db_chunk_len);

  char* id_chunk_start = index_ptr->fptr + strlen(KRAKEN_INDEX_STRING) + 1 + (idx_chunk_bounds[db_chunk_id] * 8);
  index_ptr->data_offset = idx_chunk_bounds[db_chunk_id] * 8;
  size_t id_chunk_len = (idx_chunk_bounds[db_chunk_id + 1] - idx_chunk_bounds[db_chunk_id]) * 8;
  if (db_chunk_id < dbx_chunk_bounds.size() - 1)
    id_chunk_len += 8;

  index_ptr->data = data + db_chunk_len;
  memcpy(index_ptr->data, id_chunk_start, id_chunk_len);
}

// returns first element LARGER
// half-open interval. first is inclusive, last is exclusive
// the return value has to be treated as inclusive (closed interval)
uint64_t KrakenDB::upper_bound(uint64_t first, const uint64_t last)
{
  uint64_t it;
  uint64_t count, step;
  count = last - first;

  const uint64_t orig_first = first;
  const uint64_t data_orig_first = index_ptr->mmap_at(orig_first);

  // printf("upper_bound(%llu, %llu)\n", orig_first, last);
  // float max_size = 1.0f * data_size / (1024 * 1024 * 1024);
  // printf("- max_size = %f\n", max_size);

  while (count > 0) {
    it = first;
    step = count / 2;
    it += step;

    uint64_t size_index = (it + 1 - orig_first) * 8;
    uint64_t size_data = (index_ptr->mmap_at(it + 1) - data_orig_first) * pair_size();
    // float total_size_gb = 1.0f * (size_index + size_data) / (1024 * 1024 * 1024);

    // printf("- it = %llu (size: %f)\n", it, total_size_gb);
    if (size_index + size_data + 8 <= data_size) { // if (!(data_size < size_index + size_data)) // if (!comp(value, *it))
      first = ++it;
      count -= step + 1;
    }
    else
      count = step;
  }
  return first;
}

// Chunk boundaries are whole bins; a chunk holds the bytes covering its records
void KrakenDB::prepare_chunking_compact(const uint64_t max_bytes_for_db) {
  const uint64_t n_bins = 1ull << (2 * index_ptr->indexed_nt());
  uint64_t idx_pos = 0;
  idx_chunk_bounds.push_back(0);
  dbx_chunk_bounds.push_back(0);
  while (idx_pos < n_bins) {
    const uint64_t rec_start = dbx_chunk_bounds.back();
    const uint64_t start_byte = (rec_start * rec_bits) >> 3;
    auto fits = [&](uint64_t bin_end) {
      uint64_t end_byte = (index_ptr->ef_at(bin_end) * rec_bits + 7) >> 3;
      return end_byte - start_byte + COMPACT_PAD_BYTES <= max_bytes_for_db;
    };
    // largest bin_end in [idx_pos, n_bins] such that bins [idx_pos, bin_end) fit
    uint64_t lo = idx_pos, hi = n_bins;
    while (lo < hi) {
      uint64_t mid = lo + (hi - lo + 1) / 2;
      if (fits(mid))
        lo = mid;
      else
        hi = mid - 1;
    }
    if (lo == idx_pos)
      errx(EX_USAGE, "database chunk size too small for a single minimizer bin");
    idx_pos = lo;
    const uint64_t rec_end = index_ptr->ef_at(idx_pos);
    // trailing bins without k-mers are folded into the next chunk
    if (rec_end == rec_start)
      continue;
    idx_chunk_bounds.push_back(idx_pos);
    dbx_chunk_bounds.push_back(rec_end);
    this->_chunks++;
  }
  // the last chunk must cover all remaining (empty) bins
  idx_chunk_bounds.back() = n_bins;
}

void KrakenDB::prepare_chunking(const uint64_t max_bytes_for_db) {
  data = (char*) mmap(0, max_bytes_for_db, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (data == (void *) -1)
  {
    printf("ERROR: Could not allocated memory for loading database chunks (error code %i: %s)\n", errno, strerror(errno));
    exit(-1);
  }
  data_size = max_bytes_for_db;
  mlock(data, data_size);

  this->_chunks = 0;

  if (compact) {
    prepare_chunking_compact(max_bytes_for_db);
    return;
  }

  const uint64_t max_idx_pos = 1ull << (2 * index_ptr->indexed_nt()); // excl bound (half-open interval)
  // const uint64_t max_dbx_pos = index_ptr->mmap_at(max_idx_pos); // excl bound (half-open interval)
  uint64_t idx_pos = 0;
  uint64_t dbx_pos = 0; // == index_ptr->mmap_at(0);
  idx_chunk_bounds.push_back(idx_pos);
  dbx_chunk_bounds.push_back(dbx_pos);

  while (idx_pos < max_idx_pos)
  {
    // upper_bound returns element idx_pos that is the first that "is greater" (i.e., that does not fit into "data")
    // in other terms, idx_pos closes our half-open interval (excl. upper bound), hence we don't need to decrement it
    idx_pos = upper_bound(idx_pos, max_idx_pos);
    // get starting position for next larger minimizer (that does not fit into "data")
    // in other terms, dbx_pos closes our half-open interval (excl. upper bound), hence we don't need to decrement it
    dbx_pos = index_ptr->mmap_at(idx_pos);

    // uint64_t idx_size = (idx_pos - idx_chunk_bounds.back()) * 8;
    // uint64_t dbx_size = (dbx_pos - dbx_chunk_bounds.back()) * pair_size();
    // printf("- %" PRIu64 ": %s\n", idx_size + dbx_size, (idx_size + dbx_size <= data_size ? "ok" : "NOT OK, does not fit into data"));

    // especially the last chunk can have 0 k-mers but a huge part of minimizers in it (that did not occur in any k-mer in the database)
    // in that case we can skip that chunk!
    if (dbx_pos == dbx_chunk_bounds.back())
      continue;

    idx_chunk_bounds.push_back(idx_pos);
    dbx_chunk_bounds.push_back(dbx_pos);

    this->_chunks++;
  }

  // printf("Chunk bounds\n");
  // for (uint32_t i = 0; i < this->_chunks + 1; ++i)
  //   printf("%llu\t%llu\n", idx_chunk_bounds[i], dbx_chunk_bounds[i]);

  // printf("Chunk sizes\n");
  // for (uint32_t i = 1; i < this->_chunks + 1; ++i) {
    // uint64_t idx_size = (idx_chunk_bounds[i] - idx_chunk_bounds[i - 1]) * 8;
    // uint64_t dbx_size = (dbx_chunk_bounds[i] - dbx_chunk_bounds[i - 1]) * pair_size();
    // float idx_size_gb = 1.0f * idx_size / (1024 * 1024 * 1024);
    // float dbx_size_gb = 1.0f * dbx_size / (1024 * 1024 * 1024);
    // printf("Chunk %u: %.2f GB\t%.2f GB (upper bounds: %llu\t%llu)\n", i - 1, idx_size_gb, dbx_size_gb, idx_chunk_bounds[i], dbx_chunk_bounds[i]);
    // if (idx_size + dbx_size > data_size) {
    //   printf("ERROR: %llu > %llu\n", idx_size + dbx_size, data_size);
    //   exit(1);
    // }
  // }
}

bool KrakenDB::is_minimizer_in_chunk(const uint64_t minimizer, const uint32_t db_chunk_id) const {
    return idx_chunk_bounds[db_chunk_id] <= minimizer && minimizer < idx_chunk_bounds[db_chunk_id + 1];
}

KrakenDBIndex::KrakenDBIndex() {
  fptr = NULL;
  idx_type = 1;
  nt = 0;
  data = NULL;
  data_offset = 0;
  ef_n = ef_max = 0;
  ef_l = 0;
  ef_lower = NULL;
  ef_upper = ef_samples = NULL;
}

KrakenDBIndex::KrakenDBIndex(char *ptr) {
  fptr = ptr;
  idx_type = 1;
  data = NULL;
  data_offset = 0;
  ef_n = ef_max = 0;
  ef_l = 0;
  ef_lower = NULL;
  ef_upper = ef_samples = NULL;
  if (strncmp(ptr, KRAKEN_INDEX_STRING, strlen(KRAKEN_INDEX_STRING))) {
    idx_type = 2;
    if (strncmp(ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING))) {
      idx_type = 3;
      if (strncmp(ptr, COMPACT_INDEX_MAGIC, strlen(COMPACT_INDEX_MAGIC)))
        errx(EX_DATAERR, "illegal Kraken DB index format");
    }
  }
  ptr += strlen(KRAKEN_INDEX_STRING);
  memcpy(&nt, ptr, 1);
  if (idx_type == 3) {
    uint64_t hdr[6];
    memcpy(hdr, fptr + 8, sizeof(hdr));
    ef_n = hdr[0];
    ef_max = hdr[1];
    ef_l = hdr[2];
    uint64_t upper_bits = hdr[3];
    uint64_t n_samples = hdr[4];
    if (nt < 1 || nt > 15 || ef_n != (1ull << (2 * nt)) + 1 || ef_l > 56
        || n_samples != (ef_n + (1ull << COMPACT_SELECT_SHIFT) - 1) >> COMPACT_SELECT_SHIFT)
      errx(EX_DATAERR, "compact index header is inconsistent");
    uint64_t lower_words = (ef_n * ef_l + 63) / 64 + 1;
    uint64_t upper_words = (upper_bits + 63) / 64 + 1;
    ef_lower = fptr + COMPACT_INDEX_HEADER_BYTES;
    ef_upper = (const uint64_t *) (ef_lower + lower_words * 8);
    ef_samples = ef_upper + upper_words;
  }
}

// Position of the i-th set bit in the upper bit vector
uint64_t KrakenDBIndex::ef_select(uint64_t i) const {
  uint64_t s = ef_samples[i >> COMPACT_SELECT_SHIFT];
  uint64_t rem = i & ((1ull << COMPACT_SELECT_SHIFT) - 1);
  uint64_t w = s >> 6;
  uint64_t word = ef_upper[w] & (~0ull << (s & 63));
  for (;;) {
    uint64_t c = __builtin_popcountll(word);
    if (rem < c)
      break;
    rem -= c;
    word = ef_upper[++w];
  }
  while (rem--)
    word &= word - 1;
  return w * 64 + __builtin_ctzll(word);
}

uint64_t KrakenDBIndex::ef_value(uint64_t i, uint64_t upper_pos) const {
  return ((upper_pos - i) << ef_l) | compact_read_bits(ef_lower, i * ef_l, ef_l);
}

uint64_t KrakenDBIndex::ef_at(uint64_t i) const {
  return ef_value(i, ef_select(i));
}

// Records [start, end) belong to the bin
void KrakenDBIndex::bin_range(uint64_t bin, uint64_t *start, uint64_t *end) {
  if (idx_type != 3) {
    *start = at(bin);
    *end = at(bin + 1);
    return;
  }
  uint64_t pos = ef_select(bin);
  *start = ef_value(bin, pos);
  // next set bit after pos
  uint64_t w = pos >> 6;
  uint64_t word = ef_upper[w] & ~((2ull << (pos & 63)) - 1);
  while (word == 0)
    word = ef_upper[++w];
  *end = ef_value(bin + 1, w * 64 + __builtin_ctzll(word));
}

const char *KrakenDBIndex::prefetch_bin(uint64_t bin) {
  const char *p;
  if (idx_type == 3) {
    __builtin_prefetch(ef_samples + (bin >> COMPACT_SELECT_SHIFT));
    p = ef_lower + ((bin * ef_l) >> 3);
  } else {
    p = (const char *) (get_array() + bin);
  }
  __builtin_prefetch(p);
  return p;
}

// Index version (v2 uses different minimizer sort order)
uint8_t KrakenDBIndex::index_type() {
  return idx_type;
}

// How long are bin keys (i.e., what is minimizer length in bp?)
uint8_t KrakenDBIndex::indexed_nt() {
  return nt;
}

// Convenience method, allows for testing guard
uint64_t KrakenDBIndex::mmap_at(uint64_t idx) {
  if (idx_type == 3)
    return ef_at(idx);
  uint64_t *array = (uint64_t*) (fptr + strlen(KRAKEN_INDEX_STRING) + 1);
  #ifdef TESTING
  if (idx > 1 + (1ull << (nt * 2)))
    errx(EX_SOFTWARE, "KrakenDBIndex::at() called with illegal index");
  #endif
  return array[idx];
}

// Return start of index array (skips header)
uint64_t *KrakenDBIndex::get_array_with_db_chunks() { // TODO: remove wrapper
  return (uint64_t *) data;
}

uint64_t KrakenDBIndex::at_with_db_chunks(uint64_t idx) {
  if (idx_type == 3)
    return ef_at(idx);
  uint64_t *array = get_array_with_db_chunks();
#ifdef TESTING
  if ((idx - (data_offset / 8)) > 1 + (1ull << (nt * 2)))
errx(EX_SOFTWARE, "KrakenDBIndex::at() called with illegal index");
#endif
  return array[idx - (data_offset / 8)];
}

// Return start of index array (skips header)
uint64_t *KrakenDBIndex::get_array() {
  return (uint64_t *) (fptr + strlen(KRAKEN_INDEX_STRING) + 1);
}

// Convenience method, allows for testing guard
uint64_t KrakenDBIndex::at(uint64_t idx) {
  if (idx_type == 3)
    return ef_at(idx);
  uint64_t *array = get_array();
  #ifdef TESTING
  if (idx > 1 + (1ull << (nt * 2)))
    errx(EX_SOFTWARE, "KrakenDBIndex::at() called with illegal index");
  #endif
  return array[idx];
}

} // namespace
