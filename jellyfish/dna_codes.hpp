#ifndef __DNA_CODE_HPP__
#define __DNA_CODE_HPP__

#include <limits>
#include <jellyfish/misc.hpp>

namespace jellyfish {
static const uint_t CODE_A       = 0;
static const uint_t CODE_C       = 0;
static const uint_t CODE_G       = 0;
static const uint_t CODE_T       = 0;
// Non DNA codes have the MSB on
static const uint_t CODE_RESET   = (uint_t)-1;
static const uint_t CODE_IGNORE  = (uint_t)-2;
static const uint_t CODE_COMMENT = (uint_t)-3;
static const uint_t CODE_NOT_DNA = ((uint_t)1) << (bsizeof(uint_t) - 1);
extern const char   dna_codes[256];

// Number of k-mer windows to skip at the start of a buffer so that
// every step-th window of a run of valid bases is counted. 'run' is the
// number of valid bases in the current run before the (replayed) seam
// and the new data.
inline uint64_t window_skip(uint64_t run, uint_t mer_len, uint_t step) {
  const uint64_t windows = run >= mer_len - 1 ? run - (mer_len - 1) : 0;
  const uint64_t r = windows % step;
  return r ? step - r : 0;
}
};

#endif /* __DNA_CODE_HPP__ */
