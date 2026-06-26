# libmath C jets (vere mirror)

C jet implementations for the `math.hoon` transcendental library, **vendored
from the vere runtime** (`pkg/noun`). The jets do not build in this repo — these
directories are a maintained mirror so the C and the Hoon evolve together. There
is no automatic cross-repo sync; apply changes to vere by hand.

Two copies, one per runtime word size (mirroring `/lib/lagoon`):

- `libmath/vere/`   — the 32-bit runtime (`urbit/vere`)
- `libmath/vere64/` — the 64-bit runtime (`urbit/vere-ml64`, built `-Dvere64=true`)

Each holds `noun/jets/{i/math.c, 135/tree.c, q.h, w.h}` — the full vendored files,
ready to diff/apply against the corresponding `pkg/noun/jets/…`.

## JET ATTACHMENT — all four doors fire (RESOLVED 2026-06)
All four width doors (`@rd`/`@rs`/`@rh`/`@rq`) and every transcendental arm attach
and fire on-ship, bit-exact to the Hoon: `-test %/tests/lib` is **247/247** on a
fresh hoon-version-136 fakezod, on both the 32- and 64-bit runtimes. Live per-call
cost ≈ 1 µs (`@rd`/`@rs`/`@rh`) and ≈ 1.8 µs (`@rq` f128), vs ≈ 250 µs interpreted.

**Root cause of the earlier "only @rd fires" symptom was mundane registration, not
structure:** the test ships boot **hoon-version 136**, the runtime consults exactly
one kelvin tree, and only the **135** `tree.c` carried the full four-door block —
136/137 had been given an `@rd`-only block. Porting the `@rs`/`@rh`/`@rq` cores into
136/137 (and listing all four in `_13{6,7}_non__math_d`) made every door fire.

The nesting-**depth** theory (`non → math → rd → exp` being "too deep") was a **red
herring**: equally-deep kernel jets fire (`hex → aes → ecbX → arm`,
`pen → ut → nest → nest-in`), `~/` always compiles parent axis 7 (`hoon.hoon`
`[%sgfs *] -> [%sgcn p [%$ 7] ~ q]`), and the C registration is byte-identical across
the four doors. The fix is upstream in **urbit/vere #1044** (32-bit) and **#1045**
(ml64); the in-tree `135/136/137 tree.c` here carry the full block.

## The 32-bit and 64-bit `math.c` differ by exactly 2 lines
Marshalling is **chub-based** (`u3r_chub`/`u3i_chubs`), so it is word-size-agnostic
— the divergence that broke the `@rq sub` jet cannot recur, and all 15 algorithm
cores plus the 12 single-arg wrappers are byte-identical across the two copies.
The **only** difference is the two-arg sample extraction (atan2, pow/pow-n),
because the `u3r_mean` macro's API diverged between the runtimes:

```
64-bit (vere64):  u3r_mean(cor, {u3x_sam_2, &y}, {u3x_sam_3, &x})   // brace pairs
32-bit (vere):    u3r_mean(cor,  u3x_sam_2, &y,   u3x_sam_3, &x, 0)  // flat varargs
```

Keep the two copies in sync **except** those two lines.

## What's covered
All 15 transcendentals for **all four doors — `@rd` (double), `@rs` (single),
`@rh` (half), `@rq` (quad)**, bit-exact to their Hoon door (same reduction,
coefficients, Horner order — not merely faithful): `exp log sin cos tan atan
atan2 asin acos sqt cbt pow pow-n log-2 log-10`. The `@rs`/`@rh` cores are
single-/half-precision twins of the `@rd` ones (SoftFloat `f32`/`f16`, `uint32_t`/
`uint16_t` bit pattern); `@rh` is fully native f16 (no widen-to-f32). The `@rq`
cores are native `float128_t` (SoftFloat `f128M_*` via by-value wrappers), with
the same algorithms at higher minimax degree. All use the same chub I/O — `@rq`
reads/writes **two** chubs (low 64, high 64).

`@rq` exp uses the fdlibm rational reconstruction (`1 - ((lo - r·c/(2-c)) - hi)`)
rather than a flat Horner — the flat form is only ~1.1 ULP (the dominant `1+r`
gets rounded through the whole chain), which finely-sampled MPFR sweeps expose;
the same latent issue exists in the `@rd`/`@rs`/`@rh` exp arms and is a deferred
retrofit. All other `@rq` kernels are faithful (≤1 ULP) as written.

### Rounding modes (composite arms honor the door's `r`)
The math doors carry `r=?(%n %u %d %z)` (bunt `%z`). The transcendental KERNELS
are correctly-rounded and hardcode `~(mul ^rd %n)` (no rounding axis). The
*composite* arms (`pow`/`atan2`/`tan`/`pow-n`) assemble their result with bare
door ops that honor `r` — so the jet must too. A file-static `_math_rnd` is set
from the door's `r` (gate **axis 60**: sample `[r rtol]` → door-axis 12 →
`peg(7,12)`; base hoon `^rd` with a bare `r` is axis 30) via `'n'/'u'/'d'/'z'` ->
`near_even`/`max`/`min`/`minMag`, and brackets only the composite bare ops;
kernel calls run at `near_even`.

#### Why the kernels hardcode round-nearest-even (no rounding axis)
A faithful transcendental promises that its output is the *true* value
`exp(x)`/`log(x)`/… rounded to ≤1 ULP (correctly rounded = ≤½ ULP). That promise
is defined **with respect to round-to-nearest**: the minimax coefficients, the
argument reduction, and the whole error budget are derived so the polynomial lands
within that bound *under nearest rounding*. Round-nearest-even is therefore not a
free parameter — it is part of the function's definition.

You cannot get directed rounding by simply flipping the internal ops to, say,
round-toward-zero. That yields the *approximation's* accumulated directed-rounding
error, not "the true `exp(x)` rounded toward zero" — it can be many ULPs off and
non-monotonic. Correctly directed-rounding a transcendental requires resolving the
true result to more precision than the output, just to know which way it falls at
the rounding boundary (the Table-Maker's Dilemma) — a far harder problem this
library does not attempt. So a correctly/faithfully-rounded transcendental is a
function of `x` alone: it fixes its internal arithmetic to nearest-even and takes
**no** rounding-mode axis.

The door's `r` only governs the *composite* arithmetic that assembles already-
rounded kernel outputs — `pow = exp(n·log x)` (the `n·log x` multiply), `atan2`'s
`div`/`±π`, `tan = sin/cos`'s final `div`, `pow-n`'s repeated `mul`. Those are
ordinary exactly-rounded operations where honoring `r` is well-defined and useful.
The jet matches the Hoon op-for-op: nearest in the kernels, `r` in the composites.

All three kelvin trees (`135`/`136`/`137`) now carry the **full four-door** math
block. This matters because the runtime consults exactly the tree matching the
ship's `hoon-version`: the current `brass.pill` test ships are **136**, so the
jets fire there (and on 135/137) — see the attachment note above.

## How to apply to a vere tree (`pkg/noun`)
1. Copy `noun/jets/i/math.c` → `pkg/noun/jets/i/math.c` (use the matching word-size
   copy — the two-arg `u3r_mean` form must match that runtime's macro).
2. Add the math registration block to `pkg/noun/jets/135/tree.c`: the 15
   `_135_non__math_rd_*_a[]` harms, the `_135_non__math_rd_d[]` / `_135_non__math_d[]`
   cores, and the `{ "math", 7, 0, _135_non__math_d, no_hashes }` entry in
   `_135_non_d[]` as a sibling of `lagoon`. (`noun/jets/135/tree.c` here already
   has it — diff it in.)
3. Declare the 15 `u3qi_rd_*` in `pkg/noun/jets/q.h` and `u3wi_rd_*` in `w.h`
   (after the lagoon decls). `atan2`/`pow`/`pow_n` take two atoms.
4. Add `"jets/i/math.c"` to the jet source list in `pkg/noun/build.zig`.
5. Rebuild; on a 64-bit build pass `-Dvere64=true`.

## CRITICAL: the Hoon jet structure (the gotcha that cost a day)
`math.hoon` must mirror `/lib/lagoon`'s jet structure exactly:

```
=<  math                       :: export the engine so `/+ math` + `rd:math` work
~%  %non  ..part  ~            :: file |% anchored to the shared %non chapter
|%
++  math
  ^|
  =+  [rnd=*?(%n %u %d %z)]     :: <-- REQUIRED: the engine must be a DOOR
  ~/  %math                    ::     (sample-bearing), not a bare |%
  |%
  ++  rd
    ~/  %rd
    ^|
    |_  [r=... rtol=...]
    ++  exp  ~/  %exp  |=(x=@rd ...)
```

A **sample-less `~/`-hinted core anchors its parent to the dashboard ROOT**, so
`%math` (and everything under it) misses the hot dashboard — the core SPOTs in
`_cj_mile` but `_cj_hot_mean` returns `jax=0` and no C jet ever attaches. The
`^| =+ [rnd]` sample makes `~/` resolve the parent to `%non`, exactly as
lagoon's `++ la ^| =+ [rnd] ~/ %lagoon |%` does. Verified via a `_cj_hot_mean`
trace: `math par=0 MISS` (bare |%) -> `math par=non FOUND` (door).

The chapter chain is `pen -> hex -> non -> math -> <width door> -> <fn>`, e.g.
`non/math/rd/exp`.

## Verifying the jets
The idiomatic check is the desk's own `-test` suite — run on a **jetted** binary,
`expect-eq` compares the jet's output against the hardcoded bit pattern, so a
wrong jet fails the test:

```
-test %/tests/lib            :: math-exp/log/sqrt/ainv/atan/tan/trig/derived...
```

All pass `ok=%.y` on both the 32-bit and 64-bit jetted runtimes (proving the
chub-ABI claim: identical bit-exact results at both word sizes). `test/` here
also has a standalone C harness (`build.sh` + `rd_check.c`, sharing the master
`math.c` via `-DMATH_JET_HARNESS`) and an on-ship differ.

To confirm a jet actually fires: `~>(%bout (sin:rd:math .~1))` reports ~µs/call
(jetted) vs many ms (interpreted); `sin × 200k` under `%bout` was ~430ms.

## WARNING: never `C-c` a running serf in a test loop
SIGINT mid-event makes the runtime longjmp out with the loom half-mutated (e.g.
inside a `u3i_chubs` alloc), producing a `palloc: double free` / `loom: corrupt`
that is **persisted into the snapshot** and re-detected on every later boot as a
"serf unexpectedly shut down" in `u3j_boot`. It looks exactly like a
non-deterministic jet memory bug but is purely the SIGINT. Drive the dojo with
`C-u` only (clears the line without signalling), one query at a time — or use a
generator / `-test` thread that computes everything in one event.
