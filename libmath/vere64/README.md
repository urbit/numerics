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
All 15 `@rd` transcendentals, bit-exact to the Hoon `++ rd` door (same reduction,
coefficients, Horner order — not merely faithful): `exp log sin cos tan atan
atan2 asin acos sqt cbt pow pow-n log-2 log-10`. The other width doors
(`rs`/`rh`/`rq`) are not yet jetted.

Only the **135** kelvin tree is vendored: the stock test ships (`~dev`, fresh
fakeships from `brass.pill`) are `hoon-version` 135, and that is where these jets
fire. Apply to 136/137 as well if a target ship runs those.

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
