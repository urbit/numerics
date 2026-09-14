# Spec: pure-Rust SoftFloat for NockApp jets

Status: draft, 2026-09-12. Written to be handed to an implementing agent.
Companion facts live in `~/.claude/projects/-Users-neal-urbit-numerics/memory/nockapp-lib-porting.md`.

## 1. Goal

A pure-Rust, dependency-free, allocation-free IEEE 754 software-float
library whose results are bit-identical to the Berkeley SoftFloat 3e
calls made by the numerics C jets (`lagoon/vere/noun/jets/i/lagoon.c`,
`libmath/vere/noun/jets/i/{math,unum,twoc}.c`), so that Rust jets for
NockVM produce exactly what the Hoon produces. The Hoon (`++fl` via
`++rh`/`++rs`/`++rd`/`++rq`, and `/lib/lagoon`) is the ground truth;
SoftFloat is the oracle we test against because vere's C jets already
pass against the Hoon using it.

Non-goals: exception flags, signaling-NaN semantics, `extF80`,
`softfloat_round_odd`, fused multiply-add, decimal conversion, any
transcendental (those are Chebyshev in Hoon and get their own jets
later, built on this library).

## 2. Why pure Rust, not FFI

Decided in discussion; summary so the implementer does not relitigate:

- SoftFloat keeps the rounding mode in a global (`softfloat_roundingMode`);
  math.c sets it in 93 places. NockVM runs several kernels per process.
  We want the mode as an explicit parameter, visible to the type system.
- A `cc`/bindgen dependency would be inherited by every downstream NockApp.
- Correctly rounded results are unique. Bit-exactness comes from being
  correctly rounded plus matching a handful of conventions (NaN, signed
  zero, overflow-to-int), not from copying SoftFloat's algorithms. So the
  port may use simpler integer algorithms where convenient, as long as
  every op is exactly rounded in every mode.
- The Ares-era Rust SoftFloat jets do not survive anywhere; start clean.

## 3. Scope: operations x widths

Widths: `f16` (@rh), `f32` (@rs), `f64` (@rd), `f128` (@rq). Represent
each as a newtype over the raw bits: `F16(u16)`, `F32(u32)`, `F64(u64)`,
`F128(u128)`. f128 uses native `u128` arithmetic rather than SoftFloat's
64-bit pairs.

Required per width (this is the full set the C jets call; nothing more):

| Group | Ops |
|---|---|
| Arithmetic | `add`, `sub`, `mul`, `div`, `sqrt` |
| Compare | `eq`, `lt`, `le` (IEEE quiet compares; NaN compares false) |
| Rounding | `round_to_int` (float -> float, integral, honoring mode) |
| To int | `to_i32`, `to_i64` (honoring mode; also used for `%uint`) |
| From int | `from_i32`, `from_i64`, `from_u32`, `from_u64` |
| Width | `f16<->f32`, `f32<->f64`, `f64<->f128`, and the composites the jets use (`f32->f16`, `i32->f128`, `i64->f128`) |
| Predicates | `is_nan`, `is_inf`, `is_zero`, `abs`, `neg` (pure bit ops) |

Actual C call counts, for prioritizing: mul/add/sub/lt dominate;
f128 is used only in math.c through nine `f128M_*` calls
(add sub mul div sqrt lt le eq to_i64).

## 4. Semantics that must match

1. **Rounding modes.** Hoon has five: `%n` nearest-even, `%u` toward
   +inf, `%d` toward -inf, `%z` toward zero, `%a` nearest-ties-away
   (SoftFloat `near_maxMag`). `enum Round { NearEven, Up, Down, Zero,
   NearAway }`. Every op that can round takes `Round` explicitly.
   The lagoon door sample is `?(%n %u %d %z)` but the library supports
   all five; math.c uses `%a`.
2. **NaN.** Hoon's `++ff` `bif` produces exactly one NaN per width:
   sign 0, exponent all ones, MSB of the significand set, rest zero.
   That is `0x7e00`, `0x7fc0_0000`, `0x7ff8_0000_0000_0000`,
   `0x7fff_8000_..._0000`, and it matches vere's `HALFNAN/SINGNAN/
   DOUBNAN/QUADNAN`. Rule: any NaN result is canonicalized to this value
   at the library boundary (do what vere's `_nan_unify` does). Inputs
   that are NaN with other payloads are treated as NaN but never
   propagated. No signaling/quiet distinction.
3. **Signed zero.** IEEE 754 §6.3: `x - x` and `(-x) + x` are `+0` in
   every mode except `Down`, where they are `-0`. Products and
   quotients take the XOR sign. Verify against `++fl` in phase 1, not
   just SoftFloat.
4. **Subnormals.** Full gradual underflow, no flush-to-zero. Tininess
   detection (before/after rounding) only affects the underflow *flag*,
   never the result, so it is irrelevant here.
5. **Overflow.** Directed modes saturate to the largest finite value
   in the direction that does not overflow (e.g. `Zero` and `Down`
   never round a positive result up to +inf). Nearest modes go to inf.
6. **Float -> int.** Match SoftFloat: NaN and out-of-range return the
   saturating value SoftFloat returns (`i32::MAX` for positive
   overflow and NaN, `i32::MIN` for negative). Check how each C call
   site uses the result before relying on this; some sites pre-check
   `nonfin`.
7. **Comparisons.** `eq` is IEEE equality (`-0 == +0`, NaN != NaN).
   `lt`/`le` are false when either side is NaN. No total order.
8. **sqrt.** Correctly rounded in all five modes. `sqrt(-0) = -0`,
   `sqrt(x<0) = NaN`, `sqrt(+inf) = +inf`.
9. **round_to_int.** Ties per mode (`NearEven` to even, `NearAway`
   away), preserves sign of zero, inf/NaN pass through.

## 5. API shape

```rust
#![no_std]
#![forbid(unsafe_code)]

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Round { NearEven, Up, Down, Zero, NearAway }

pub trait SoftFloat: Copy + Eq {
    type Bits: Copy;                 // u16 / u32 / u64 / u128
    const EXP_BITS: u32;
    const SIG_BITS: u32;
    const NAN: Self;                 // canonical (section 4.2)

    fn from_bits(b: Self::Bits) -> Self;
    fn to_bits(self) -> Self::Bits;

    fn add(self, o: Self, r: Round) -> Self;
    fn sub(self, o: Self, r: Round) -> Self;
    fn mul(self, o: Self, r: Round) -> Self;
    fn div(self, o: Self, r: Round) -> Self;
    fn sqrt(self, r: Round) -> Self;
    fn round_to_int(self, r: Round) -> Self;

    fn eq(self, o: Self) -> bool;
    fn lt(self, o: Self) -> bool;
    fn le(self, o: Self) -> bool;

    fn to_i32(self, r: Round) -> i32;
    fn to_i64(self, r: Round) -> i64;
    fn from_i32(v: i32, r: Round) -> Self;
    fn from_i64(v: i64, r: Round) -> Self;
    fn from_u32(v: u32, r: Round) -> Self;
    fn from_u64(v: u64, r: Round) -> Self;

    fn is_nan(self) -> bool;
    fn is_inf(self) -> bool;
    fn is_zero(self) -> bool;
    fn abs(self) -> Self;
    fn neg(self) -> Self;
}

pub struct F16(pub u16); pub struct F32(pub u32);
pub struct F64(pub u64); pub struct F128(pub u128);

// width conversions as inherent fns or a `Convert<T>` trait:
impl F32 { pub fn to_f16(self, r: Round) -> F16; pub fn to_f64(self) -> F64; }
// ... exact widenings take no Round; narrowings do.
```

Implementation notes:

- One generic core over `(EXP_BITS, SIG_BITS)` with a wide-enough
  unsigned type (`u32` for f16, `u64` for f32, `u128` for f64, and a
  hand-rolled 256-bit or a `u128`-pair path for f128 mul/div/sqrt).
  Round-pack once, in one function, used by every op. SoftFloat's
  `s_roundPackToF32.c` / `s_normRoundPackToF32.c` are the model.
- No `f32`/`f64` hardware types anywhere in the crate, including tests
  (a `#![deny]` lint or a grep in CI). This is what makes the crate a
  valid oracle-of-record for the VM.
- No `unsafe`, no allocation, no global state, no dependencies.
- `const fn` where stable Rust permits; not required.
- SoftFloat 3e sources for reference are vendored at
  `~/urbit/SoftBLAS/SoftFloat/source/` (per-op files ~100-200 lines).

## 6. Crate layout

Under `lagoon/nockapp/` in the numerics repo (repurposing the directory
PR #17 proposed for a Hoon copy; the Hoon stays single-sourced in
`lagoon/desk` with `..part` -> `..ut`):

```
lagoon/nockapp/
  Cargo.toml              # workspace
  softfloat/              # this spec (phases 1-3)
  softblas/               # axpy, dot, gemm, scal per width (phase 4)
  lagoon-jets/            # HotEntry table + jet fns (phase 5)
  README.md
```

`softfloat` and `softblas` must not depend on nockvm; only
`lagoon-jets` does. Crate name on disk `softfloat`; if it is ever
published, the name is taken on crates.io, so pick then.

## 7. Testing

Three layers, cheapest first. All must be green before a phase is done.

1. **Oracle fuzz against C SoftFloat.** Add `softfloat-sys` (crates.io,
   Berkeley SoftFloat 3 bindings) as a *dev-dependency only*. For every
   op, width, and mode: exhaustive over all f16 pairs where feasible
   (2^32 pairs for binary ops is fine for a nightly job; sample in the
   default test), and a structured sampler for f32/f64/f128 that
   overweights zeros, subnormals, exponent boundaries, exact ties,
   inf, NaN, and values whose exact result lies within an ulp of a
   rounding boundary. Compare bits after canonicalizing NaN on both
   sides. Set the C global mode before each C call.
2. **Hoon-derived vectors.** The existing numerics test suites already
   encode Hoon results (lagoon `tests/lib/*.hoon`, libmath tests).
   Generate a vector file from a fakezod for the scalar ops (a few
   thousand cases per width/mode via `++rs` etc. under `~(. rs r)`),
   check it into `softfloat/tests/vectors/`, and assert against it.
   This is the only check that catches a SoftFloat-vs-Hoon divergence.
3. **NockVM jet test mode (phase 5).** NockVM has a built-in
   differential mode: with `NOCK_TEST_JETS` set to a list of jet paths,
   the interpreter runs both the jet and the raw Nock and bails with a
   deterministic exit on any mismatch (see `parse_test_jets` in
   `crates/nockapp/src/kernel/boot.rs` for the path syntax, and
   `interpreter.rs` around the `get_jet_test_mode` call). Run the
   lagoon Hoon test suite under it. This replaces the vere
   "rename the `~%` parent" de-jetting trick.

## 8. NockVM integration facts (phase 5)

- Reference checkout: `~/zorp/nockchain-official` (org is `nockchain/`).
  Hoon kelvin 138; `hoonc` 0.2.0 on PATH compiles a lib dir with
  `lib/` and `sur/` (both `/-` and `/+` work).
- Jet signature: `fn(&mut Context, Noun) -> nockvm::jets::Result`.
  Gate sample is `slot(subject, 6, &space)`; the lagoon door's `rnd`
  sits in the door sample (the C jets read the math door's `r` at
  axis 60; confirm lagoon's axis on a compiled core before hardcoding).
- Registration: `HotEntry = (path, axis, fn)`; path for a lib with
  parent `..ut` is
  `[K_138, "one","two","tri","qua","pen", "non","lagoon", <arm>]`
  (model: `%zeke` entries in `crates/zkvm-jetpack/src/hot.rs`).
  Requires the Hoon to say `~%  %non  ..ut  ~` (not `..part`).
- The app passes the table as `boot::setup(&kernel, cli, &HOT, name, None)`.
- Atom I/O: `Atom::as_u64`, `Atom::as_u64_pair` (already exists,
  "SoftFloat-compatible ordered pair of 64-bit words", little-word
  order, use it for f128 I/O), `as_ne_bytes` for bulk ray data,
  `Atom::from_ubig` / `IndirectAtom::new_raw_mut_bytes` for results.
- Ray layout (`/sur/lagoon`): `data` is one atom, row-major, element
  `i` at bit offset `i * 2^bloq` from the LSB, with a single 1 pin bit
  immediately above the last element. Strip/restore the pin in the jet.
- Errors: return `BAIL_EXIT` (deterministic) for bad input, never
  `BAIL_FAIL`, so test mode can compare with the Nock crash.

## 9. Phases and acceptance

| Phase | Deliverable | Done when |
|---|---|---|
| 0 | Workspace, `Round`, trait, f32 `add`/`mul`, oracle harness wired | oracle fuzz green for f32 add/mul, 5 modes |
| 1 | f32 + f64 complete (section 3) | oracle fuzz + Hoon vectors green; signed-zero and NaN rules verified against `++fl` |
| 2 | f16 | exhaustive binary-op check vs oracle |
| 3 | f128 | oracle fuzz green; `as_u64_pair` round-trip test |
| 4 | `softblas`: axpy, dot, gemm, scal x 4 widths, same accumulation order as SoftBLAS | bit-equal to SoftBLAS on fuzzed inputs (SoftBLAS via FFI as dev-dep, or vectors dumped from its munit suite) |
| 5 | `lagoon-jets`: hot table for the arms lagoon.c covers today (28), Hoon `..ut` change, hoonc build of `lagoon/desk` | lagoon Hoon tests pass under `NOCK_TEST_JETS` with every jet enabled |

Do phases in order; each is independently mergeable. Keep the per-op
files small and named after the SoftFloat file they replace so review
can be done side by side.

## 10. Open questions for the implementer to settle early

- Confirm `%a` is reachable from any lagoon arm today (math.c uses it;
  lagoon.c maps it but the door type excludes it). Support it anyway.
- Confirm the exact saturation/NaN return of `to_i32`/`to_i64` at each
  C call site before choosing between SoftFloat semantics and a
  bail (section 4.6).
- Confirm the lagoon door's `rnd` axis on a hoonc-compiled core.

## 11. Errata and decisions (2026-09-12, RustFloat implementation)

Recorded during implementation in `sigilante/RustFloat`; the sections
above are left as written.

- **§4.1 `%a` is directed rounding away from zero, not `near_maxMag`.**
  In `++fl`, `rau` maps `%a` to `%ce` (magnitude ceiling) and `swr` never
  flips it; `toj` treats `%a` like `%u` on the magnitude. Verified on the
  un-jetted `++ff` core: `toi` of 2.25 under `%a` is 3, and `1 + 1e-10`
  under `%a` is `0x3f80.0001`; the checked-in Hoon vectors decide every
  discriminating case (718 of 718) for directed-away. lagoon.c and math.c
  map `%a` to `softfloat_round_near_maxMag`, which is wrong relative to the
  Hoon, but no door (`rd`/`rs`/`rq`/`rh`, math, lagoon) admits `%a`, so it
  is unreachable. The library implements `Round::Away` as directed.
- **§4.5 / Hoon overflow.** Hoon `++fl` overflows to infinity in every
  mode (`lug` ends with a mode-independent max-exponent check). The
  library follows IEEE 754 §7.4 and saturates in directed modes; the Hoon
  is to be fixed (urbit/urbit#7426, urbit/numerics#82). The Hoon-vector
  test counts rows of exactly that shape rather than skipping them.
- **§4.6 float to int on NaN.** Decision: the library returns what
  SoftFloat returns under the specialization vere links, which is
  ARM-VFPv2 (`vere/ext/softfloat/build.zig`): NaN gives 0, positive
  overflow `MAX`, negative overflow `MIN`. (8086 would give `MIN` for all
  three; the `i32::MAX` for NaN stated in §4.6 corresponds to the RISC-V
  specialization and matches no vere build.) The jet layer bails on NaN
  and out-of-range inputs instead of relying on these values; math.c
  already pre-checks NaN at every `to_i32`/`to_i64` call site.
- **§4.9 round_to_int.** Hoon has no float-to-float rounding; `toi` gives
  `(unit @s)` and `san` of it turns `-0.4` into `+0`. The library keeps
  IEEE 754 §5.9 (sign of zero preserved).
- **§5 crate naming and §6 layout.** The package is `sdfloat` (Berkeley's
  name is not ours to use); the repo is `RustFloat` with `crates/sdfloat`
  and `crates/sdfloat-oracle`. Phases 4 and 5 stay in urbit/numerics under
  `lagoon/nockapp/`.
- **§7.1 oracle.** `softfloat-sys` 0.1.4 does not build on macOS arm64
  (its `c99` dependency's CMake, then "build rules are not implemented for
  the current target_arch and target_os"). The oracle is instead a
  dev-only crate compiling the vendored SoftFloat 3e sources with `cc`,
  ARM-VFPv2 specialization, rounding mode thread-local and passed per
  call. `Away` is checked as `round_max` for non-negative results and
  `round_min` for negative ones.
- **§7.2 Hoon vectors.** `urbit eval` on the `rs`/`rd`/`rq`/`rh` doors
  runs vere's SoftFloat jets, not the Hoon. Vectors are generated by
  calling `++ff` directly (`~(. ff [[8 23 --127] %z])` and so on), which
  has no jets. See `scripts/gen-hoon-vectors.py`.
