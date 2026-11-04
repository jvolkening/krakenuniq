/*
 * compact_db: convert a final (sorted, indexed) KrakenUniq database into the
 * compact format described in compactdb.hpp. The input files are not modified.
 *
 * The conversion is lossless: every bin is decoded again and compared with the
 * source while writing, and -v additionally looks up every source k-mer in the
 * finished compact database.
 */

#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"
#include "compactdb.hpp"
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace kraken;

static string In_db, In_idx, Out_db, Out_idx;
static bool Verify = false;

static void usage(int exit_code = EX_USAGE) {
  cerr << "Usage: compact_db -d db -i index -o compact_db -x compact_index [-v]" << endl
       << endl
       << "  -d filename   Input database (after db_sort/set_lcas)" << endl
       << "  -i filename   Input index (v2)" << endl
       << "  -o filename   Output compact database" << endl
       << "  -x filename   Output compact index" << endl
       << "  -v            Afterwards look up every k-mer in the compact database" << endl
       << "  -h            Print this message" << endl;
  exit(exit_code);
}

static uint64_t file_size(const string &fn) {
  struct stat sb;
  if (stat(fn.c_str(), &sb) != 0)
    err(EX_OSERR, "unable to stat %s", fn.c_str());
  return sb.st_size;
}

static void write_index(const uint64_t *offsets, uint8_t nt, const string &fn) {
  const uint64_t n = (1ull << (2 * nt)) + 1;
  const uint64_t max_value = offsets[n - 1];
  unsigned l = 0;
  if (max_value > n)
    while ((n << (l + 1)) <= max_value)
      l++;
  const uint64_t upper_bits = n + (max_value >> l) + 1;
  const uint64_t n_samples = (n + (1ull << COMPACT_SELECT_SHIFT) - 1) >> COMPACT_SELECT_SHIFT;
  const uint64_t lower_words = (n * l + 63) / 64 + 1;
  const uint64_t upper_words = (upper_bits + 63) / 64 + 1;
  const uint64_t total = COMPACT_INDEX_HEADER_BYTES + (lower_words + upper_words + n_samples) * 8;

  QuickFile f(fn, "w", total);
  char *base = f.ptr();
  memcpy(base, COMPACT_INDEX_MAGIC, strlen(COMPACT_INDEX_MAGIC));
  base[7] = (char) nt;
  uint64_t hdr[6] = { n, max_value, l, upper_bits, n_samples, 0 };
  memcpy(base + 8, hdr, sizeof(hdr));
  char *lower = base + COMPACT_INDEX_HEADER_BYTES;
  uint64_t *upper = (uint64_t *) (lower + lower_words * 8);
  uint64_t *samples = upper + upper_words;
  for (uint64_t i = 0; i < n; i++) {
    uint64_t v = offsets[i];
    if (i && v < offsets[i - 1])
      errx(EX_DATAERR, "index offsets are not monotonic at bin %llu", (unsigned long long) i);
    uint64_t pos = (v >> l) + i;
    upper[pos >> 6] |= 1ull << (pos & 63);
    if ((i & ((1ull << COMPACT_SELECT_SHIFT) - 1)) == 0)
      samples[i >> COMPACT_SELECT_SHIFT] = pos;
    compact_write_bits(lower, i * l, l, v & ((1ull << l) - 1));
  }
  cerr << "compact_db: wrote index (" << total << " bytes, " << l << " low bits)" << endl;
}

int main(int argc, char **argv) {
  int opt;
  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:i:o:x:vh")) != -1) {
    switch (opt) {
      case 'd': In_db = optarg; break;
      case 'i': In_idx = optarg; break;
      case 'o': Out_db = optarg; break;
      case 'x': Out_idx = optarg; break;
      case 'v': Verify = true; break;
      case 'h': usage(0); break;
      default: usage();
    }
  }
  if (In_db.empty() || In_idx.empty() || Out_db.empty() || Out_idx.empty())
    usage();

  QuickFile db_file(In_db);
  QuickFile idx_file(In_idx);
  KrakenDB db(db_file.ptr(), db_file.size());
  KrakenDBIndex index(idx_file.ptr());
  if (db.is_compact() || index.is_compact())
    errx(EX_USAGE, "input is already compact");
  if (index.index_type() != 2)
    errx(EX_DATAERR, "only v2 indexes are supported");
  db.set_index(&index);

  const uint8_t nt = index.indexed_nt();
  if (nt < 1 || nt > 15)
    errx(EX_DATAERR, "minimizer length %u is not supported", (unsigned) nt);
  const uint64_t n_bins = 1ull << (2 * nt);
  const uint64_t key_ct = db.get_key_ct();
  const uint64_t key_len = db.get_key_len();
  const uint64_t key_bits = db.get_key_bits();
  const uint64_t pair_sz = db.pair_size();
  const unsigned k = db.get_k();
  if (key_ct == 0)
    errx(EX_DATAERR, "database is empty");
  if (k > 31 || nt > k)
    errx(EX_DATAERR, "unsupported k (%u) / minimizer length (%u)", k, (unsigned) nt);
  if (db.header_size() + key_ct * pair_sz > db_file.size())
    errx(EX_DATAERR, "database file is truncated");
  if (idx_file.size() < 8 + 8 * (n_bins + 1))
    errx(EX_DATAERR, "index file is truncated");

  const uint64_t *offsets = index.get_array();
  if (offsets[0] != 0 || offsets[n_bins] != key_ct)
    errx(EX_DATAERR, "index does not match database (offsets %llu..%llu, %llu keys)",
         (unsigned long long) offsets[0], (unsigned long long) offsets[n_bins],
         (unsigned long long) key_ct);

  const char *pairs = db.get_pair_ptr();
  const uint64_t key_mask = key_bits >= 64 ? ~0ull : (1ull << key_bits) - 1;

  // Pass 1: dictionary of taxids
  cerr << "compact_db: collecting taxids ..." << endl;
  unordered_set<uint32_t> taxid_set;
  uint32_t last_taxid = 0;
  bool have_last = false;
  for (uint64_t i = 0; i < key_ct; i++) {
    uint32_t taxid;
    memcpy(&taxid, pairs + i * pair_sz + key_len, 4);
    if (have_last && taxid == last_taxid)
      continue;
    if (taxid >> COMPACT_MAX_TAXID_BITS)
      errx(EX_DATAERR, "taxid %u (k-mer number %llu) needs more than %u bits; "
           "cannot create compact database", taxid, (unsigned long long) i,
           COMPACT_MAX_TAXID_BITS);
    taxid_set.insert(taxid);
    last_taxid = taxid;
    have_last = true;
  }
  vector<uint32_t> dict(taxid_set.begin(), taxid_set.end());
  sort(dict.begin(), dict.end());
  unordered_map<uint32_t, uint32_t> dict_index;
  for (size_t i = 0; i < dict.size(); i++)
    dict_index[dict[i]] = i;
  const unsigned dict_bits = dict.size() > 1 ? compact_ceil_log2(dict.size()) : 0;

  CompactKmerCodec codec;
  codec.init(k, nt, INDEX2_XOR_MASK);
  if (codec.t_bits > 64)
    errx(EX_DATAERR, "k=%u with minimizer length %u does not fit the compact format (needs %u bits)",
         k, (unsigned) nt, codec.t_bits);
  const uint64_t rec_bits = codec.t_bits + dict_bits;
  const uint64_t data_bytes = (key_ct * rec_bits + 7) / 8;
  const uint64_t dict_bytes = (dict.size() * 4 + 7) & ~7ull;
  const uint64_t total = COMPACT_DB_HEADER_BYTES + dict_bytes + data_bytes + COMPACT_PAD_BYTES;
  cerr << "compact_db: " << dict.size() << " distinct taxids, " << rec_bits
       << " bits per k-mer (" << codec.t_bits << " key + " << dict_bits << " taxid)" << endl;

  QuickFile out(Out_db, "w", total);
  char *base = out.ptr();
  memcpy(base, COMPACT_DB_MAGIC, 8);
  uint64_t hdr[11] = { 0, k, nt, key_ct, dict.size(), dict_bits, codec.pos_bits,
                       rec_bits, COMPACT_MAX_TAXID_BITS, data_bytes, 0 };
  memcpy(base + 8, hdr + 1, sizeof(hdr) - 8);
  memcpy(base + COMPACT_DB_HEADER_BYTES, dict.data(), dict.size() * 4);
  char *rec_base = base + COMPACT_DB_HEADER_BYTES + dict_bytes;

  // Pass 2: re-encode bin by bin
  cerr << "compact_db: encoding k-mers ..." << endl;
  vector<pair<uint64_t, uint32_t> > bin;
  for (uint64_t b = 0; b < n_bins; b++) {
    const uint64_t start = offsets[b], end = offsets[b + 1];
    if (end < start || end > key_ct)
      errx(EX_DATAERR, "invalid index offsets for bin %llu", (unsigned long long) b);
    if (end == start)
      continue;
    bin.clear();
    for (uint64_t r = start; r < end; r++) {
      const char *p = pairs + r * pair_sz;
      uint64_t kmer = 0;
      memcpy(&kmer, p, key_len);
      kmer &= key_mask;
      uint32_t taxid;
      memcpy(&taxid, p + key_len, 4);
      uint64_t bin_key, t;
      codec.encode(kmer, &bin_key, &t);
      if (bin_key != b)
        errx(EX_DATAERR, "k-mer %llu is stored in bin %llu but belongs to bin %llu "
             "(is the database sorted with this index?)", (unsigned long long) r,
             (unsigned long long) b, (unsigned long long) bin_key);
      if (codec.decode(bin_key, t) != kmer)
        errx(EX_SOFTWARE, "internal error: k-mer %llu does not survive re-encoding",
             (unsigned long long) r);
      bin.push_back(make_pair(t, dict_index[taxid]));
    }
    sort(bin.begin(), bin.end());
    for (size_t j = 0; j < bin.size(); j++) {
      if (j && bin[j].first == bin[j - 1].first)
        errx(EX_DATAERR, "duplicate k-mer in bin %llu", (unsigned long long) b);
      uint64_t pos = (start + j) * rec_bits;
      compact_write_bits(rec_base, pos, codec.t_bits, bin[j].first);
      compact_write_bits(rec_base, pos + codec.t_bits, dict_bits, bin[j].second);
    }
    if ((b & 0xFFFFF) == 0xFFFFF && isatty(fileno(stderr)))
      fprintf(stderr, "\r  %.1f %%", 100.0 * b / n_bins);
  }
  out.close_file();

  write_index(offsets, nt, Out_idx);

  const uint64_t old_size = file_size(In_db) + file_size(In_idx);
  const uint64_t new_size = file_size(Out_db) + file_size(Out_idx);
  cerr << "compact_db: database " << file_size(In_db) << " -> " << file_size(Out_db)
       << " bytes, index " << file_size(In_idx) << " -> " << file_size(Out_idx)
       << " bytes, total " << old_size << " -> " << new_size << " bytes ("
       << (100.0 * new_size / old_size) << " % of original)" << endl;

  if (Verify) {
    cerr << "compact_db: verifying ..." << endl;
    QuickFile cdb_file(Out_db);
    QuickFile cidx_file(Out_idx);
    KrakenDB cdb(cdb_file.ptr(), cdb_file.size());
    KrakenDBIndex cidx(cidx_file.ptr());
    cdb.set_index(&cidx);
    if (cdb.get_key_ct() != key_ct || cdb.get_k() != k)
      errx(EX_SOFTWARE, "verification failed: header mismatch");
    const uint64_t step = n_bins > 4096 ? n_bins >> 12 : 1;
    for (uint64_t b = 0; b < n_bins; b += step) {
      uint64_t s, e;
      cidx.bin_range(b, &s, &e);
      if (s != offsets[b] || e != offsets[b + 1])
        errx(EX_SOFTWARE, "verification failed: bin %llu has range [%llu,%llu), expected [%llu,%llu)",
             (unsigned long long) b, (unsigned long long) s, (unsigned long long) e,
             (unsigned long long) offsets[b], (unsigned long long) offsets[b + 1]);
    }
    for (uint64_t i = 0; i < key_ct; i++) {
      const char *p = pairs + i * pair_sz;
      uint64_t kmer = 0;
      memcpy(&kmer, p, key_len);
      kmer &= key_mask;
      uint32_t taxid;
      memcpy(&taxid, p + key_len, 4);
      uint32_t *v = cdb.kmer_query(kmer);
      if (v == NULL || *v != taxid)
        errx(EX_SOFTWARE, "verification failed for k-mer %llu", (unsigned long long) i);
    }
    if (cdb.count_taxons() != db.count_taxons())
      errx(EX_SOFTWARE, "verification failed: taxon counts differ");
    cerr << "compact_db: verification passed" << endl;
  }
  return 0;
}
