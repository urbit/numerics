# Numerics review triage — Gnome/Angel pass (libmath, lagoon, saloon)

Local working doc — NOT committed. Tracks the six-agent review findings and the
fix plan. Severity: P0 crash-on-valid-input in shipped code; P1 footguns; P2
math.hoon transcendental hazards (folded into the Chebyshev rewrite #18); P3
docs/traceability.

Verification note: agents could not build Hoon; "confirmed" = reproduced via
Python IEEE arithmetic or tight source-tracing. Ship-verify the crash paths.

---

## P0 — fix now, three PRs   [ALL OPEN: #51, #52, #53]

### PR-1  lagoon: %fixp conj (and dotc) crash   [#51 OPEN, green]
- `trans-scalar` %fixp handles only %abs; %conj falls to `!!`.  `dotc = (dot
  (conj a) b)` so dotc crashes on every %fixp ray.  conj/dotc doccords claim
  "identity on real kinds" — a contract lie.
- FIX: add `%conj |=(b=@ b)` to the %fixp branch (matches every other real
  kind).  Add a lagoon-fixp conj/dotc test.
- Decision: none (mechanical).

### PR-2  lagoon: %i754 %mod   [#52 OPEN -- Hoon spec landed; VERE jet follow-on flagged]
- `(need (toi (div a b)))` crashes on b=0 / Inf / NaN (toi -> ~, need ~ -> bail)
  instead of IEEE NaN/Inf.  %unum %mod already guards.
- Test suite is incoherent: test-mod-nearest-6r expects 5 mod 3 = -1 (IEEE
  remainder, round-nearest); test-mods-7r expects = 2 (C fmod).  Can't both hold.
- DECISION (proposed): **C fmod / truncate-toward-zero**, matching %unum (which
  the user already chose this session) — so test-mods-7r is right, fix
  test-mod-nearest-6r.  Implement trunc + guard divisor/operand non-finite ->
  return NaN (not crash).
- Ship-verify the crash and the jet-vs-fallback behavior of test-mods-7r.

### PR-3  saloon: eig assert robustness   [#53 OPEN, green; conflicts w/ #50 -> rebase]
- `symmetric()` uses exact `=`: a Gram matrix `M^T M` from lagoon mmul is
  symmetric only up to ULP -> crash on a legitimately-symmetric input.
- `hermitian()` exact-equal + `cconj` does `0 - Im`: under %d/%u rounding
  `0 - 0 = -0` (IEEE 6.3), so conj(a+0i) != a+0i bit-exactly -> crash on any
  Hermitian matrix with non-nearest rounding.
- DECISION (proposed): relax both to an **approximate** check (|m_ij -
  (conj) m_ji| <= atol + rtol*|m_ji|) so genuinely (skew-)symmetric-up-to-
  rounding input passes while a truly asymmetric matrix still crashes (preserves
  the "don't silently run the wrong algorithm" intent; the design rejected
  silent symmetrization).  Handles ±0.0 too.
- Also fold in the cheap P1 doccord surfacing for this PR: note the 60-sweep
  cap / ~& warning, the rtol-width requirement, and that bare `sa` rtol=0x1 is
  unusable (call +sake).

---

## P1 — footguns   [DONE: PR #56 (numerics) -- cplx-reduce docs, eig rtol sane default + width guard; verified on ~hex]
- bare-`sa` default rtol=0x1 (denormal) -> 60-sweep cap always fires. (touched in PR-3 docs; consider a sane default or a require-sake guard)
- rtol width must match component width -> silent wrong results. (PR-3 docs; consider a check)
- max/min/any/all on %cplx crash via fun-scalar %gth->ord. (document at those arms)

## P2 — math.hoon transcendentals -> Chebyshev rewrite (#18)   [DONE: documented in #18 comment as acceptance criteria]
Confirmed hazards reinforcing #18; #18 should add domain guards + iteration caps,
not just better polynomials:
- log(z): O(z^2) hang for large z; diverges for z<=0.
- exp(large -x): catastrophic cancellation -> garbage/NaN.
- asin/acos crash on |x|>1 not exactly +-1.
- sqt(-0.0) hard-crash; sqt(2^-149) -> +inf.
- NaN-detection guards are no-ops (dis is unsigned, gte ... 0 always true).
Saloon element-wise exp/sin/log inherit these (eig path is insulated by capped fsqt).
ACTION: append these to the #18 follow-up comment; not fixed here.

## P3 — docs/traceability   [DONE: PR #54 (round-bankers half-even, rd.log(-inf)=NaN, tiny, sur, check/change, neq/data=, READMEs)]
- round-bankers is half-up, not banker's (name+doc lie).
- `tiny` constant: docstring/example/comment mutually contradict (value = smallest subnormal).
- rd.log(-inf) returns 0 not NaN (copy-paste from exp).  [also behavioral - small]
- sur/lagoon.hoon: `tail=*` comment is CORRECT in intent (deliberate: lets per-kind specialization data, e.g. %fixp prec [a b], slot in WITHOUT rewriting jet axes each time). Angel overstated this -> AUGMENT comment to mention %fixp uses it; do NOT call it a lie. (Separately: %cplx comment lists @cs/@cd but four widths ship @ch-@cq -- that one is real staleness.)
- lagoon + saloon READMEs stale (%real old name, kinds "planned" but shipped, missing arms, no eig).
- saloon typos: "(skew-)symmetry", neq "equal to", `data=data=` in gth/gte/geq examples; EIG-DESIGN.md "design proposal" though A1/A2 shipped.
- complex abs returns complex [|z|,0] not real (buried); equ bit-equality (+0!=-0, NaN!=NaN) undocumented.
- lagoon `check` undocumented; `change` negative->uint XXX untracked.

---

## PR sequencing
PR-1, PR-2 both touch lagoon.hoon (different regions: trans-scalar conj vs the
%mod branches) -> whichever merges second rebases.  PR-3 is saloon-only,
independent.  All branch off main.
