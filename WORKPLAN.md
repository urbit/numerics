# Numerics Work Plan

Last updated: 2026-06-30.  Based on `doc/audit-2026-06-30.md` (55 findings,
11H/27M/17L).  All 53 of 55 findings addressed in-session (2026-06-30); two
remain open.

---

## Phase 0 — Code bugs  ✓ COMPLETE

- [x] **H1** `libmath` — `++rs` `+huge` wrong bit pattern fixed:
  `` `@rs`0x7f7f.ffff ``.
- [x] **H4** `lagoon` — `++range` @rh (bloq 4) now updates `shape.meta`.
  C jet was already correct; Hoon now matches.
  _Still needed: add a `lagoon-builders` test for `(range meta [.~~0 .~~1] .~~0.25)` shape._
- [x] **H5** `lagoon` — `++mod` `%cplx` now emits `~|(%cplx-mod-unsupported !!)`;
  `++mod` doccord updated with per-kind behavior.

---

## Phase 1 — High-severity doc fixes

### libmath  ✓ COMPLETE

- [x] **H2** `math.hoon` — `+round` examples fixed; "decimal places" → "significant figures."
- [x] **H3** `math.hoon` — `@rq` `+tan` pi example replaced with `.~~~0  :: TODO: verify on ship`.
- [x] **H11** `complex.hoon` — Three doors fixed: `+rs` → `+rd`/`+rh`/`+rq`.
- [x] **H12** `unum.hoon` — `++fdp` silent truncation documented.

### lagoon + saloon  ✓ COMPLETE (except H10)

- [x] **H6** `lagoon/README.md` — `%real` → `%i754`; `eq, ne` struck from to-make.
- [x] **H7** `lagoon/README.md` — Current Status rewritten to describe shipped state.
- [x] **H8** `EIG-DESIGN.md` — Header updated to "A1+A2 shipped (PR #47)."
- [x] **H9** `EIG-DESIGN.md` — §8 ±0.0 claim corrected; `+cnear` behavior documented.
- [x] **H10** `saloon/tests` — **DONE (2026-06-30).** Added
  `saloon/desk/tests/lib/saloon-rays.hoon` with 12 tests covering `++exp`,
  `++sin`, `++cos`, `++tan`, `++log`, `++log-2`, `++log-10`, `++sqrt`,
  `++cbrt`, `++pow-n`, `++pow`, and shape-preservation; all 12 pass on ship.

---

## Phase 2 — Medium doc fixes  ✓ COMPLETE (except M1)

- [x] **M1** `math.hoon` — **DONE (2026-06-30).** Chose option (b): genuine
  relative semantics.  Added `atol=_.0` to all four doors (@rs/@rd/@rh/@rq)
  after `rtol`; changed `++is-close` to `|p-r| < atol + rtol×|r|` (NumPy
  convention); fixed `++factorial` base-case guards from `(is-close x .0/1)`
  to `(lth (abs x) rtol)` / `(lth (abs (sub x .1)) rtol)`; updated all
  doccords; fixed pre-existing `.1000001` typo in `++rs` doccord example.
- [x] **M2–M8** `math.hoon` — All seven medium math fixes applied.
- [x] **M9–M12** `lagoon.hoon` — All four medium lagoon fixes applied.
- [x] **M13** `lagoon/README.md` — Arm list updated; section retitled.
- [x] **M14–M20** `saloon.hoon` / `EIG-DESIGN.md` — All seven medium saloon fixes applied.
- [x] **M21–M23** `complex.hoon` — All three medium complex fixes applied.
- [x] **M24–M25** `unum.hoon` — Both medium unum fixes applied.
- [x] **M26** `libmath/README.md` — Posit count corrected (3 standard + 2 extensions).
- [x] **M27** `libmath/tools/NEXT-STEPS.md` — Created.

---

## Phase 3 — Low-severity doc fixes  ✓ COMPLETE

- [x] **L1–L5** `math.hoon` — All five low math fixes applied.
- [x] **L6–L9** `lagoon.hoon` / `sur/lagoon.hoon` — All four low lagoon fixes applied.
- [x] **L10** `lagoon/README.md` — `@rpq` added with caveat.
- [x] **L11–L13** `saloon.hoon` — All three low saloon fixes applied.
- [x] **L14** `complex.hoon` — `++abs` return-type note added in all four doors.
- [x] **L15** `unum.hoon` — Seven quire op arms now have per-arm doccords.
- [x] **L16** `twoc.hoon` — `++mul` sign-independence rationale added.
- [x] **L17** `README.md` (top-level) — `/lib/complex` added to desk contents.

---

## Open items (0 remaining)

All 55 findings resolved.

---

## Jet compliance tracking

Jets exist for `%i754` operations in `libmath` (math doors: exp/sin/cos/tan/
asin/acos/atan/atan2/pow-n/log/log-10/log-2/pow/sqt/cbt) and `lagoon` (array
ops: range/linspace/max/argmax/min/argmin/cumsum/stack/transpose/diag/trace/
dot/dotc/mmul/abs/conj/add-scal…div-scal/mod-scal/add-rays…mod-rays/gth/gte/
lth/lte).  All other kinds (`%cplx`, `%unum`, `%fixp`, `%int2`) fall through
to Hoon — no jet drift risk for those kinds.

**Divergences found and resolved in this audit:**

| Arm | Was | Now | Jet status |
|-----|-----|-----|------------|
| `++range` bloq=4 | Hoon bug: shape not updated | Fixed | Jet was already correct (lagoon.c:3577) |
| `++mod` %cplx | Bare `!!` crash, no message | Diagnostic `~|` added | Jet returns `u3_none` → falls to Hoon (same result) |
| `+huge` @rs | `0x7f80.0000` (+Inf) | `0x7f7f.ffff` (max finite) | Not referenced by jet; Hoon-only fix |

**No open jet divergences.**

**Policy for future jetted arms:**
1. Verify the Hoon spec is correct first (against oracle/SoftFloat/SoftPosit).
2. Implement the C jet to match the Hoon exactly.
3. Add tests in both the Hoon test file and any C-level harness (`vere64/test/`)
   at every supported bloq.
4. If a kind falls through to Hoon from the jet (`u3_none`), verify the Hoon
   path is correct — the jet absence must not mask a Hoon bug.

---

## Paper readiness gate

Before submitting a USTJ article:
- [x] All Phase 0 bugs fixed.
- [x] Phase 0 bugs ship-verified on ~zod (numjet-zod pier, 2026-06-30).
- [x] All Phase 1 doc fixes applied.
- [x] **H10** resolved — Saloon ray-level test coverage (12 tests, all pass).
- [x] **M1** resolved — `is-close` now uses `|p-r| < atol + rtol×|r|`.
- [x] H2/H3 corrected — no wrong examples.
- [x] Jet compliance table current with no open divergences.
