# libmath tooling — next steps

Status as of 2026-06-30.

## Oracle harnesses (Python)

- `posit_check.py` — SoftPosit-based oracle for `/lib/unum`.  Exhaustive
  posit8 pair testing, sampled posit16/32.  Reference: SoftPosit `pX2_*` path
  (es=2 at any width).
- `complex_check.py` — NumPy `complex64`/`complex128` oracle for `/lib/complex`.
- `fixed_check.py` — Python integer arithmetic oracle for `/lib/fixed`.
- `cheb_check.py` — mpmath oracle for Chebyshev transcendental verification
  (PR #18 follow-up).

## C test harnesses

- `rd_exp_check.c` — spot-check double exp jet against SoftFloat.
- `rq_check.c` / `rq_inv_check.py` — quad-precision arithmetic verification.
- `rh_check.c` / `rh_sweep.c` — half-precision jet verification; sweep tests
  all 65536 f16 values for selected operations.

## Known remaining work

- [ ] Chebyshev rewrite of transcendentals (PR #18): replace naive Taylor
  series in `/lib/math` with range-reduced Chebyshev approximations.  The
  oracle is `cheb_check.py`.  Acceptance criteria: domain guards, iteration
  caps, domain notes in doccords.
- [ ] `@rq` `+invlog2` value verified (audit 2026-06-30): constant is correct;
  `:: TODO check` comment removed.
- [ ] Populate doccord examples for `@rh` `+log-10`, `@rh` `+log-2`, `@rq`
  `+log-10`, `@rq` `+log-2` (requires ship with hoon-136+ and jetted @rh/@rq).
- [ ] `++fdp` equal-length assertion: consider adding `?>` guard for mismatched
  list lengths (currently silently truncates).
- [ ] Posit transcendental accuracy: current exp/sin/cos/log are naive Taylor
  series accurate near 0.  Full range requires range reduction + quire
  accumulation (see libmath/README.md §2).
- [ ] Valids (`@rvb`/`@rvh`/`@rvs`): interval arithmetic — lowest priority.
