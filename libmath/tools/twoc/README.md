# `/lib/twoc` jet validation harnesses

Standalone C programs that check the arithmetic in the `/lib/twoc` jet
(`libmath/vere/noun/jets/i/twoc.c`) against an independent two's-complement
reference, BEFORE any ship build.  Each copies the jet's arithmetic VERBATIM and
compares to a from-scratch reference built on hardware `__int128` / `mpz`.

- `twoc_native_test.c` — the native paths (`c3_d` for wid<=64, `__int128` for
  wid<=128).  Widths 8/16/17/32/64/100, all sign edges (min, -1, max, 0), 12
  ops.  **8280 checks.**
- `twoc_gmp_test.c` — the GNU MP path, cross-checked vs the same hardware
  reference at widths 8..100, plus hand-verified 200-bit cases GMP-only.
  **5692 checks.**

## Run

```
cc -O2 -Wall -Wno-unused-function twoc_native_test.c -o /tmp/t && /tmp/t
cc -O2 -Wall -Wno-unused-function -I/opt/homebrew/include -L/opt/homebrew/lib \
   twoc_gmp_test.c -lgmp -o /tmp/tg && /tmp/tg
```

Both print `checks=<N> fails=0` on success.
