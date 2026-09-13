# Porting the lagoon C jets to Rust (phase 5 working notes)

Read this before touching `crates/lagoon-jets/src/jets/*.rs`.

## What a jet must do

Reproduce, bit for bit, what the Hoon arm in `hoon/lib/lagoon.hoon`
(a generated copy of `../desk/lib/lagoon.hoon`) returns for `%i754` rays,
using `sdfloat`/`sdblas`. The vere C jet in
`../vere/noun/jets/i/lagoon.c` (`u3qi_la_<arm>_i754`, wrapper `u3wi_la_<arm>`)
is the reference for *how* (which BLAS routine, accumulation order, special
cases); the Hoon is the reference for *what*. When they disagree, the Hoon
wins here: NockVM's test mode compares the jet against the Nock.

Rules:

- **Punt, never bail.** Any input the jet does not handle (another `kind`,
  wrong rank, invalid ray, unknown rounding mode) returns `Err(JetErr::Punt)`
  so the Nock runs. A jet that bails where the Hoon succeeds is a bug.
- Only `%i754` (bloq 4..=7). Other kinds punt (the C also handles `%int2`
  for some arms; that is a later phase).
- Read the door's rounding mode with `ray::rounding(subject, &space)`
  (axis 30). Elementwise ops honor it; comparisons do not need it.
- Results reuse the input `meta` noun (`ray::build(context, a.meta, data)`)
  unless the Hoon builds a new one (`mmul`, `scalar-to-ray` reductions,
  `transpose`, `diag`, `range`, `linspace`): then `ray::build_meta`. Match the
  Hoon's `tail` (often `~`, i.e. `D(0)`, not the input's).
- Reductions that return a scalar ray use the Hoon's `+scalar-to-ray`: shape
  `(reap (lent shape) 1)`, same bloq/kind/tail as the input.
- NaN results come out canonical from sdfloat; the Hoon produces the same
  canonical NaN, so no extra handling.
- Width dispatch: `by_bloq!(a.bloq, worker(args))` with a generic
  `fn worker<F: BlasFloat>(..)`. Elements move as raw `u128` bit patterns
  (`ray::elems`, `ray::pack`); convert with `from_raw::<F>` / `to_raw`.
- Add `ray::trace("<arm>")` as the first line of each jet.

## Verify

```sh
cargo build --release
# does it fire?  (prints one line per call)
LAGOON_JET_TRACE=1 LAGOON_TESTS=<substr> ./target/release/lagoon-kick hoon/run-tests.jam
# is it right?  (interpreter runs jet AND nock, crashes with %jest on mismatch)
NOCK_TEST_JETS=k.138/one/two/tri/qua/pen/non/lagoon/<hint> LAGOON_TESTS=<substr> \
  ./target/release/lagoon-kick hoon/run-tests.jam
```

`<hint>` is the `~/  %name` in the Hoon (e.g. `add-rays`, `mod-scal`, `argmax`),
comma-separate several. `LAGOON_TESTS` filters rows by substring of
`<file>/<arm>` (e.g. `arithmetic/test-mul`, `compare-reduce/`, `jet-parity/`).
With no filter the whole suite (153 tests) runs in ~1.5 s. Finish with the
whole suite in test mode for your arms.

Do not run `hoonc` or edit anything under `hoon/`; the Hoon is frozen for
this pass (if `hoon/.jam-ok` is missing, a rebuild is in progress: wait).
If you find a Hoon bug, note it in your report instead.

## Known divergences (not yours to fix)

- `%mod` for `%i754` was reconciled in numerics #78 (2026-09-13): the Hoon
  is now `a - b * san(need(toi_r(a / b)))` with the quotient rounded in the
  door mode and a crash (`need ~`) on a non-finite quotient, and the vere C
  jets do the same (vere #1057 carries it upstream). The jet follows that:
  punt on a non-finite quotient so the Nock crashes.
- Hoon `++fl` overflows to infinity in directed modes where IEEE saturates
  (urbit/urbit#7426). sdfloat saturates. Any test that overflows under
  `%z`/`%d`/`%u` will mismatch until the Hoon is fixed; report, don't work
  around.

## Findings from the first full port (2026-09-12)

All 28 arms in `hot.rs` are ported; the whole suite under test mode for all
28 hints gives one failure, non-jet:

- `lagoon-compare-reduce/test-argmax`: a **NockVM bug**, not Hoon or jet. The
  built-in `find` jet (`nockvm/src/jets/list.rs`, `util::find`) compares
  elements with `raw_equals` (word/pointer identity). `+argmax` looks up a
  `@rq` produced by `cut` in a list produced by `rip`; both are fresh
  indirect atoms with equal contents, so `find` returns `~` and `+:~` bails.
  Widths 4..6 are direct atoms and pass. With the `argmax` jet on and test
  mode off the test passes, because the jet never calls `find`. Fix belongs
  upstream (structural equality in `util::find`).

Hoon-vs-vere-C divergences the port surfaced (the Rust follows the Hoon):

- `+abs` is `?:((gte b .0) b (mul b .-1))`, not a sign-bit clear: `-0` stays
  `-0`, NaN comes back canonical. `f32_abs` in lagoon.c differs on both.
- `+range`/`+linspace` iterate and compare per step; lagoon.c derives an
  element count up front. Rounding-sensitive steps can differ in length.
- `+dot` and `+cumsum` accumulate right to left from `+0`; sdblas `dot`
  (like SoftBLAS) accumulates left to right, so `dot` does not use it.
- `+ravel`'s `~/  %ravel` hint is commented out in the Hoon, so the `ravel`
  jet is registered but unreachable.
- `+argmax`/`+argmin` indexed the boxed extreme with `~[0 0]` (rank 2 only);
  fixed in the desk to `(reap (lent shape) 0)` on 2026-09-12.

Debug switches (all read from the environment; never set in production):
`LAGOON_JET_TRACE`, `LAGOON_JET_SABOTAGE` (add-rays wrong on purpose),
`LAGOON_JET_DISABLE=max,argmax,...` (named reduction jets punt).

## The float doors (phase 6)

`crates/hoon-float-jets` is a separate crate with no lagoon dependency; the
same rules apply (punt on anything odd, rounding mode at axis 30, operands
read as the low `width` bits like `sea:ff`). Verified three ways: the
lagoon suite under all 68 hints (unchanged: the same two non-jet failures),
`hoon/float-tests.hoon` under the 40 float hints (only the overflow rows
mismatch), and the per-pair diff of `hoon/float-diff.hoon` with and without
the jets (592 of 6400 pairs differ, all ±MAX vs ±inf). Do not "fix" the
jets to overflow to infinity: the Hoon is what changes (urbit/urbit#7426).
