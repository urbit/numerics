/* Exhaustive verification of +posit-unit at posit8: encode every one of
 * the 2^32 possible k-bit numerators u (k=4n=32, n=8) via a C mirror of
 * /lib/unum's +bit (the same encode() logic already validated in
 * posit_unit_check.py), tally per-pattern counts, and print them.
 *
 * This is a from-scratch re-derivation, not a transliteration copy-paste
 * of posit_unit_check.py's Python -- deliberately, so a translation bug
 * in one isn't invisible to the other. Cross-validated against the
 * Python encode() for a large random sample (see cross_check.py in this
 * directory) before trusting the full 2^32 run.
 *
 * Since the k-bit dyadic construction only ever encodes nonnegative
 * values (u/2^k is always in [0, 1 - 2^-k]), neg=0 always here -- see
 * posit_unit_check.py's reachable_patterns() for why negative-signed
 * patterns are unreachable by this construction.
 *
 * Build:  cc -O3 -o posit8_exhaustive posit8_exhaustive.c
 * Run:    ./posit8_exhaustive > posit8_exhaustive_counts.txt
 */
#include <stdio.h>
#include <stdint.h>

#define N 8
#define MAXPOS ((1 << (N - 1)) - 1)   /* 63 */
#define MSK    ((1 << N) - 1)         /* 0xff */

/* Floor division/modulus for signed operands -- explicit, not relying on
 * C's truncating "/" or implementation-defined behavior of ">>" on
 * negative values, even though in practice most compilers do arithmetic
 * shift. Matches Python's ">>"/"%" semantics on negative operands exactly. */
static int64_t fdiv(int64_t a, int64_t b) {
  int64_t q = a / b, r = a % b;
  if (r != 0 && ((r < 0) != (b < 0))) q--;
  return q;
}

static int bitlen64(uint64_t a) {
  int n = 0;
  while (a) { n++; a >>= 1; }
  return n;
}

/* Mirrors posit_unit_check.py's encode(neg=0, e, a, n=8). */
static unsigned encode8(int64_t e, uint64_t a) {
  if (a == 0) return 0;
  int lead = bitlen64(a) - 1;
  int64_t x = e + (int64_t)lead;
  uint64_t frac = a & ((lead < 64) ? ((1ULL << lead) - 1) : ~0ULL);
  int64_t r = fdiv(x, 4);
  int64_t elo = x - 4 * r;              /* always in [0,4), matches Python's x - 4*r */

  if (r >= N - 2) return (unsigned)MAXPOS;
  if (r <= -(N - 1)) return 1u;

  int64_t regval, regwid;
  if (r >= 0) { regval = ((1LL << (r + 1)) - 1) << 1; regwid = r + 2; }
  else        { regval = 1;                          regwid = -r + 1; }

  int64_t totw = regwid + 2 + lead;
  uint64_t pay = ((uint64_t)regval << (2 + lead)) | ((uint64_t)elo << lead) | frac;
  int64_t pw = N - 1;
  uint64_t mag;

  if (totw <= pw) {
    mag = pay << (pw - totw);
  } else {
    int64_t sh = totw - pw;
    uint64_t keep   = pay >> sh;
    uint64_t guard  = (pay >> (sh - 1)) & 1;
    uint64_t sticky = (sh >= 2) ? ((pay & ((1ULL << (sh - 1)) - 1)) != 0) : 0;
    uint64_t lsbit  = keep & 1;
    if (guard && (sticky || lsbit)) keep++;
    if (keep > (uint64_t)MAXPOS) keep = (uint64_t)MAXPOS;
    mag = keep;
  }
  return (unsigned)mag;
}

int main(void) {
  uint64_t hist[1 << N] = {0};
  const int64_t e = -32; /* k = 4n = 32 for posit8 */
  /* u ranges over exactly [0, 2^32) -- a 32-bit counter wrapping to 0
   * after UINT32_MAX covers this with a do/while, avoiding a 64-bit
   * comparison against 1ULL<<32 on every iteration of a 4-billion-
   * iteration hot loop. */
  uint32_t u = 0;
  do {
    hist[encode8(e, u)]++;
    u++;
  } while (u != 0);

  for (unsigned p = 0; p < (1u << N); p++)
    if (hist[p]) printf("%u %llu\n", p, (unsigned long long)hist[p]);
  return 0;
}
