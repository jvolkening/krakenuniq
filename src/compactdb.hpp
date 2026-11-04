/*
 * Compact on-disk representation of a (sorted, indexed) KrakenUniq database.
 *
 * Compact DB file (all header fields little-endian uint64):
 *    0  magic "KRAKCDB1"
 *    8  k
 *   16  minimizer length (nt)
 *   24  number of k-mers
 *   32  dictionary size D (number of distinct taxids)
 *   40  dictionary index width in bits
 *   48  bits used for the minimizer position within a k-mer
 *   56  record width R in bits (t_bits + dictionary index bits)
 *   64  max bits allowed for a taxid (validated when the file is opened)
 *   72  size of the bit-packed record array in bytes
 *   80  reserved
 *   88  reserved
 *   96  dictionary: D * uint32 (sorted taxids), padded to a multiple of 8 bytes
 *  ...  bit-packed records, followed by 16 bytes of zero padding
 *
 * Each record is a fixed-width bit field: t (t_bits) then dictionary index.
 * Records of a bin (all k-mers sharing a minimizer) are stored consecutively,
 * sorted by t, in the same bin order as the original database, so record
 * number r of the original DB's index maps to bit offset r * R.  The k-mer is
 * not stored directly; t holds (minimizer position, orientation, remaining
 * bases) and the minimizer itself is implied by the bin.
 *
 * Compact index file ("KRAKCIX" + nt byte, then little-endian uint64 header):
 *    8  number of values N (= 4^nt + 1)
 *   16  largest value (= number of k-mers)
 *   24  number of low bits l
 *   32  length of the upper bit vector in bits
 *   40  number of select samples
 *   48,56 reserved
 *   64  lower bits (N * l bits, padded), upper bit vector, select samples
 * The values are the cumulative bin offsets of a standard KrakenUniq v2 index,
 * stored as an Elias-Fano sequence.
 */

#ifndef COMPACTDB_HPP
#define COMPACTDB_HPP

#include <cstdint>
#include <cstring>

namespace kraken {

static const char *const COMPACT_DB_MAGIC = "KRAKCDB1";        // 8 bytes
static const char *const COMPACT_INDEX_MAGIC = "KRAKCIX";      // 7 bytes + nt byte
static const uint64_t COMPACT_DB_HEADER_BYTES = 96;
static const uint64_t COMPACT_INDEX_HEADER_BYTES = 64;
static const uint64_t COMPACT_PAD_BYTES = 16;
static const unsigned COMPACT_MAX_TAXID_BITS = 26;
// XOR mask for minimizer bin keys of v2 indexes (scrambles minimizer sort order)
static const uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;
static const unsigned COMPACT_SELECT_SHIFT = 7;          // one sample per 128 ones

// Read n (<= 64) bits starting at bit position pos. Needs 9 readable bytes
// after the byte holding pos (guaranteed by the padding).
static inline uint64_t compact_read_bits(const char *base, uint64_t pos, unsigned n) {
  if (n == 0)
    return 0;
  const char *p = base + (pos >> 3);
  unsigned sh = pos & 7;
  uint64_t lo;
  memcpy(&lo, p, 8);
  uint64_t v = lo >> sh;
  if (sh + n > 64)
    v |= (uint64_t)(unsigned char)p[8] << (64 - sh);
  return n == 64 ? v : (v & ((1ull << n) - 1));
}

// Bits must be zero beforehand (freshly created files are zero-filled).
static inline void compact_write_bits(char *base, uint64_t pos, unsigned n, uint64_t v) {
  if (n == 0)
    return;
  if (n < 64)
    v &= (1ull << n) - 1;
  char *p = base + (pos >> 3);
  unsigned sh = pos & 7;
  uint64_t lo;
  memcpy(&lo, p, 8);
  lo |= v << sh;
  memcpy(p, &lo, 8);
  if (sh + n > 64)
    p[8] |= (char)(v >> (64 - sh));
}

static inline uint64_t compact_revcomp(uint64_t kmer, unsigned n) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

static inline unsigned compact_ceil_log2(uint64_t x) {
  unsigned b = 0;
  while (b < 64 && (1ull << b) < x)
    b++;
  return b;
}

// Reversible mapping between a canonical k-mer and (bin key, t).
// xor_mask must be the (already masked) index v2 XOR mask.
struct CompactKmerCodec {
  unsigned k, nt, pos_bits, rest_bits, t_bits;
  uint64_t xor_mask, win_mask;

  CompactKmerCodec() : k(0), nt(0), pos_bits(0), rest_bits(0), t_bits(0),
                       xor_mask(0), win_mask(0) {}

  void init(unsigned k_, unsigned nt_, uint64_t full_xor_mask) {
    k = k_;
    nt = nt_;
    pos_bits = compact_ceil_log2(k - nt + 1);
    rest_bits = 2 * (k - nt);
    t_bits = pos_bits + 1 + rest_bits;
    win_mask = (1ull << (2 * nt)) - 1;
    xor_mask = full_xor_mask & win_mask;
  }

  // Leftmost (lowest-order) window attaining the minimum is the minimizer
  // occurrence; this matches KrakenDB::bin_key().
  inline void encode(uint64_t kmer, uint64_t *bin_key, uint64_t *t) const {
    uint64_t best = ~0ull;
    unsigned best_pos = 0, best_orient = 0;
    uint64_t km = kmer;
    for (unsigned i = 0; i <= k - nt; i++) {
      uint64_t w = km & win_mask;
      uint64_t r = compact_revcomp(w, nt);
      uint64_t c = w < r ? w : r;
      uint64_t v = xor_mask ^ c;
      if (v < best) {
        best = v;
        best_pos = i;
        best_orient = (w != c);
      }
      km >>= 2;
    }
    uint64_t low = kmer & ((1ull << (2 * best_pos)) - 1);
    uint64_t high = kmer >> (2 * (best_pos + nt));
    uint64_t rest = (high << (2 * best_pos)) | low;
    *bin_key = best;
    *t = ((uint64_t)((best_pos << 1) | best_orient) << rest_bits) | rest;
  }

  inline uint64_t decode(uint64_t bin_key, uint64_t t) const {
    uint64_t rest = t & ((1ull << rest_bits) - 1);
    uint64_t po = t >> rest_bits;
    unsigned orient = po & 1;
    unsigned pos = po >> 1;
    uint64_t c = bin_key ^ xor_mask;
    uint64_t w = orient ? compact_revcomp(c, nt) : c;
    uint64_t low = rest & ((1ull << (2 * pos)) - 1);
    uint64_t high = rest >> (2 * pos);
    return (high << (2 * (pos + nt))) | (w << (2 * pos)) | low;
  }
};

} // namespace kraken

#endif
