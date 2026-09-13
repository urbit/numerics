  ::  /tests/lib/i754rand
::::
::    /lib/i754rand (the porcelain): ++uni (float: +rs/+rd/+rh/+rq/+rs-oo/
::    +rd-oo), ++alias (Vose's method), ++dist (++rd + ++rs).  ++uni's
::    float arms are exact bit constructions (+sun then a power-of-two
::    multiply), so their expected values are checked against an
::    independent Python IEEE-754 encoder, not the Hoon float printer.
::    ++dist's value spot-checks (including the alias table) involve
::    /lib/math transcendentals (log/sqrt), so their expected values are
::    read off this same ship's actual output (sanity-checked by hand
::    against the closed-form algorithm first) rather than independently
::    rederived in Python, which isn't guaranteed to match this
::    codebase's correctly-rounded kernels to the last bit; ++dist's
::    moment tests are regression checks against the THEORETICAL targets
::    (mean/variance) at a fixed seed, with tolerances wide enough to
::    absorb real sampling noise at n=50000 but tight enough to catch an
::    actual algorithm bug.  Two real bugs were caught by this on-ship
::    testing during development (alias-table's worklist assignment
::    reversed; +binomial comparing against the wrong accumulator) --
::    see NEXT-STEPS.md.
::
/+  *test,
    rand,
    i754rand,
    math
|%
::  +rs/+rd/+rh/+rq: exact bit constructions, checked against an
::  independent Python IEEE-754 encoder (not the Hoon float printer).
++  test-float-rs  ^-  tang
  %+  expect-eq
    !>(`@rs`0x3dee.6d78)
    !>(out:(rs:uni:i754rand (from-atom:seed:rand %sm64 0)))
++  test-float-rd  ^-  tang
  %+  expect-eq
    !>(`@rd`0x3f95.072f.63b9.b5e0)
    !>(out:(rd:uni:i754rand (from-atom:seed:rand %sm64 0)))
++  test-float-rh  ^-  tang
  %+  expect-eq
    !>(`@rh`0x39af)
    !>(out:(rh:uni:i754rand (from-atom:seed:rand %sm64 0)))
++  test-float-rq  ^-  tang
  %+  expect-eq
    !>(`@rq`0x3ffd.3cd5.4372.cbe9.c441.5072.f63b.9b5e)
    !>(out:(rq:uni:i754rand (from-atom:seed:rand %sm64 0)))
++  test-float-rs-oo  ^-  tang
  %+  expect-eq
    !>(`@rs`0x3dee.6d7c)
    !>(out:(rs-oo:uni:i754rand (from-atom:seed:rand %sm64 0)))
++  test-float-rd-oo  ^-  tang
  %+  expect-eq
    !>(`@rd`0x3f95.072f.63b9.b5f0)
    !>(out:(rd-oo:uni:i754rand (from-atom:seed:rand %sm64 0)))
::  ++dist: value spot-checks against actual on-ship computation (these
::  involve /lib/math transcendentals -- log/sqrt -- so the expected
::  values below are read off this same ship's output, not independently
::  rederived in Python, which isn't guaranteed to match this codebase's
::  correctly-rounded kernels to the last bit).  Each is sanity-checked
::  against the closed-form algorithm by hand (e.g. -ln(u)/lambda for
::  +expon) before being trusted as a regression constant.
++  test-normal  ^-  tang
  %+  expect-eq
    !>(`@rd`.~-0.9479938949723624)
    !>(out:(normal:rd:dist:i754rand (from-atom:seed:rand %sm64 0)))
::  determinism: same seed -> same output (spot check; +normal's rejection
::  loop is the one arm in ++dist where a threading bug would most likely
::  show up as nondeterminism).
++  test-normal-determinism  ^-  tang
  %+  expect-eq
    !>((normal:rd:dist:i754rand (from-atom:seed:rand %phil 7)))
    !>((normal:rd:dist:i754rand (from-atom:seed:rand %phil 7)))
++  test-normal-mv  ^-  tang
  %+  expect-eq
    !>(`@rd`.~8.104012210055275)
    !>(out:(normal-mv:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~10 .~2))
++  test-normal-mv-bad-sigma  ^-  tang
  %-  expect-fail
  |.((normal-mv:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0 .~-1))
::  -ln(u)/lambda, u = +rd-oo of the first SplitMix64 output (0.0205...):
::  -ln(0.0205352...)/2 ~ 1.9428 -- matches to displayed precision.
++  test-expon  ^-  tang
  %+  expect-eq
    !>(`@rd`.~1.942806871605541)
    !>(out:(expon:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~2))
++  test-expon-bad-rate  ^-  tang
  %-  expect-fail
  |.((expon:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0))
::  alpha=2 (>=1, direct Marsaglia-Tsang path).
++  test-gamma-ge1  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.7179344171650754)
    !>(out:(gamma:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~2))
::  alpha=0.5 (<1, boost path: gamma(1.5) * u^(1/0.5)).
++  test-gamma-lt1  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.05447897345403265)
    !>(out:(gamma:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0.5))
++  test-gamma-bad-shape  ^-  tang
  %-  expect-fail
  |.((gamma:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0))
++  test-beta  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.23622515961139348)
    !>(out:(beta:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~2 .~3))
++  test-chi2  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.8261342997997718)
    !>(out:(chi2:rd:dist:i754rand (from-atom:seed:rand %sm64 0) 3))
++  test-chi2-bad-df  ^-  tang
  %-  expect-fail
  |.((chi2:rd:dist:i754rand (from-atom:seed:rand %sm64 0) 0))
++  test-student-t  ^-  tang
  %+  expect-eq
    !>(`@rd`.~-0.7137613421937223)
    !>(out:(student-t:rd:dist:i754rand (from-atom:seed:rand %sm64 0) 5))
++  test-student-t-bad-df  ^-  tang
  %-  expect-fail
  |.((student-t:rd:dist:i754rand (from-atom:seed:rand %sm64 0) 0))
::  u ~ 0.0205 < 0.5 -> %.y.
++  test-bernoulli  ^-  tang
  %+  expect-eq
    !>(`?`%.y)
    !>(out:(bernoulli:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0.5))
++  test-bernoulli-bad-prob  ^-  tang
  %-  expect-fail
  |.((bernoulli:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~1.5))
++  test-geometric  ^-  tang
  %+  expect-eq
    !>(`@ud`11)
    !>(out:(geometric:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0.3))
::  p=1 edge case: returns 1 without consuming a draw (rng unchanged).
++  test-geometric-p1  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  out  r  (geometric:rd:dist:i754rand r .~1)
  ;:  weld
    %+  expect-eq  !>(`@ud`1)  !>(out)
    %+  expect-eq  !>((from-atom:seed:rand %sm64 0))  !>(r)
  ==
++  test-geometric-bad-prob  ^-  tang
  %-  expect-fail
  |.((geometric:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0))
::  Moment tests (rand-spec.md section 10 item 4): sample mean/variance at
::  50k draws, FIXED seed -- deterministic, so this is a regression test
::  against the theoretical constants, not a flaky statistical one.
::  Tolerances are generous relative to the actual sampling error at
::  n=50000 (stderr(mean) ~ sigma/sqrt(n), stderr(var) ~ sigma^2*sqrt(2/n))
::  so a real algorithm bug (wrong constant, sign error, wrong formula)
::  fails loudly while ordinary run-to-run noise from a fixed seed does not
::  (there being only one fixed seed, "run-to-run" here really means
::  "robust to this exact computation shifting by a few ULP if the
::  underlying math kernels ever change").
++  test-moments-normal  ^-  tang
  =/  mv  (moments (from-atom:seed:rand %phil 0) normal:rd:dist:i754rand 50.000)
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.02]) -.mv .~0))
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) +.mv .~1))
  ==
++  test-moments-expon  ^-  tang
  =/  mv  (moments (from-atom:seed:rand %phil 0) |=(r=rng:rand (expon:rd:dist:i754rand r .~1)) 50.000)
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) -.mv .~1))
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) +.mv .~1))
  ==
++  test-moments-gamma  ^-  tang
  =/  mv  (moments (from-atom:seed:rand %phil 0) |=(r=rng:rand (gamma:rd:dist:i754rand r .~2)) 50.000)
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) -.mv .~2))
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) +.mv .~2))
  ==
::  +moments: draw .n samples via .f from .r0, return [mean=@rd var=@rd]
::  (population variance).  Shared by the three moment tests above.
++  moments
  |=  [r0=rng:rand f=$-(rng:rand [@rd rng:rand]) n=@]
  ^-  [mean=@rd var=@rd]
  =/  mth  ~(. rd:math [%n .~1e-13 .~0])
  =/  r    r0
  =/  i    0
  =/  sum  .~0
  =/  ssq  .~0
  |-  ^-  [@rd @rd]
  ?:  =(i n)
    =/  ct    (sun:mth n)
    =/  mean  (div:mth sum ct)
    [mean (sub:mth (div:mth ssq ct) (mul:mth mean mean))]
  =^  x  r  (f r)
  %=  $
    i    +(i)
    sum  (add:mth sum x)
    ssq  (add:mth ssq (mul:mth x x))
  ==
::  +alias +build: worked example matching rand-spec.md section 6.1's own
::  doc comment, hand-verified (weights [1,1,2]/4 -> true probabilities
::  [0.25,0.25,0.5]; this table's prob/alias reproduce that -- see the
::  arm's doccord for the by-hand trace).  This test caught a real bug
::  during development: the small/large worklist assignment was reversed,
::  which +build-loop's own doc comment now explains.
++  test-alias-build  ^-  tang
  %+  expect-eq
    !>(`alias-table:i754rand`[n=3 prob=~[.~0.75 .~0.75 .~1] alias=~[2 2 2]])
    !>((build:alias:i754rand ~[.~1 .~1 .~2]))
::  +alias +build: the other worklist-draining order (large empties last).
++  test-alias-build-other-order  ^-  tang
  %+  expect-eq
    !>(`alias-table:i754rand`[n=4 prob=~[.~1 .~0.6666666666666666 .~0.6666666666666666 .~0.6666666666666666] alias=~[0 0 0 0]])
    !>((build:alias:i754rand ~[.~3 .~1 .~1 .~1]))
++  test-alias-build-empty  ^-  tang
  %-  expect-fail
  |.((build:alias:i754rand ~))
++  test-alias-build-bad-prob  ^-  tang
  %-  expect-fail
  |.((build:alias:i754rand ~[.~0 .~0]))
::  +alias +draw: value regression, from the table above.
++  test-alias-draw  ^-  tang
  =/  t  (build:alias:i754rand ~[.~1 .~1 .~2])
  %+  expect-eq
    !>(`[out=@ud rng:rand]`[2 [%sm64 s=0x3c6e.f372.fe94.f82a]])
    !>((draw:alias:i754rand t (from-atom:seed:rand %sm64 0)))
::  +categorical: thin wrapper over +draw:alias -- same table, same
::  result, by construction.
++  test-categorical  ^-  tang
  =/  t  (build:alias:i754rand ~[.~1 .~1 .~2])
  %+  expect-eq
    !>((draw:alias:i754rand t (from-atom:seed:rand %sm64 0)))
    !>((categorical:rd:dist:i754rand (from-atom:seed:rand %sm64 0) t))
::  +poisson: Knuth branch (lambda<10) and PTRS branch (lambda>=10) value
::  regressions, plus determinism-implying mean sanity already covered by
::  hand-verified means during development (2000-draw sample means came
::  out within 0.5% of lambda for both 5 and 20 -- not re-asserted here as
::  an exact regression to keep this test fast).  Crashes if lambda<=0.
++  test-poisson-knuth  ^-  tang
  %+  expect-eq
    !>(`@ud`2)
    !>(-:(poisson:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~5))
++  test-poisson-ptrs  ^-  tang
  %+  expect-eq
    !>(`@ud`14)
    !>(-:(poisson:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~20))
++  test-poisson-bad-rate  ^-  tang
  %-  expect-fail
  |.((poisson:rd:dist:i754rand (from-atom:seed:rand %sm64 0) .~0))
::  +binomial: value regression (also exercises the fix for a real bug
::  caught during development -- the accept test originally compared the
::  uniform draw against the individual pmf term instead of the
::  accumulated CDF, which either returned 0 too often or ran the
::  recurrence past n and crashed with subtract-underflow; a 2000-draw
::  sample mean check during development came out at 5.98 vs the
::  theoretical n*p=6). Crashes if n*min(p,1-p) >= 30.
++  test-binomial  ^-  tang
  %+  expect-eq
    !>(`@ud`2)
    !>(-:(binomial:rd:dist:i754rand (from-atom:seed:rand %sm64 0) 20 .~0.3))
++  test-binomial-n-too-large  ^-  tang
  %-  expect-fail
  |.((binomial:rd:dist:i754rand (from-atom:seed:rand %sm64 0) 1.000.000 .~0.5))
::  +dirichlet: value regression (sums to 1 by construction; spot-check
::  the sum explicitly too, since that's the actual mathematical
::  invariant, not just a frozen regression number).  Crashes on an
::  empty alphas list.
++  test-dirichlet  ^-  tang
  =/  out  `(list @rd)`-:(dirichlet:rd:dist:i754rand (from-atom:seed:rand %sm64 0) ~[.~1 .~2 .~3])
  =/  mth  ~(. rd:math [%n .~1e-13 .~0])
  =/  total  (add:mth (snag 0 out) (add:mth (snag 1 out) (snag 2 out)))
  ;:  weld
    %+  expect-eq
      !>(`(list @rd)`~[.~0.04688412483230968 .~0.4265186606662172 .~0.5265972145014731])
      !>(out)
    %+  expect-eq  !>(%.y)  !>((is-close:mth total .~1))
  ==
++  test-dirichlet-empty  ^-  tang
  %-  expect-fail
  |.((dirichlet:rd:dist:i754rand (from-atom:seed:rand %sm64 0) ~))
::  ++rs mirror: representative spot-checks (not exhaustive re-coverage
::  of every ++rd test) confirming the mechanical @rd->@rs
::  re-instantiation actually produces working, independently-computed
::  @rs output -- not just that it compiles.
++  test-rs-normal  ^-  tang
  %+  expect-eq
    !>(`@rs`.-0.5933847)
    !>(out:(normal:rs:dist:i754rand (from-atom:seed:rand %sm64 0)))
++  test-rs-gamma  ^-  tang
  %+  expect-eq
    !>(`@rs`.1.0119848)
    !>(out:(gamma:rs:dist:i754rand (from-atom:seed:rand %sm64 0) .2))
++  test-rs-expon  ^-  tang
  %+  expect-eq
    !>(`@rs`.1.0752765)
    !>(out:(expon:rs:dist:i754rand (from-atom:seed:rand %sm64 0) .2))
++  test-rs-poisson-ptrs  ^-  tang
  %+  expect-eq
    !>(`@ud`14)
    !>(out:(poisson:rs:dist:i754rand (from-atom:seed:rand %sm64 0) .20))
++  test-rs-binomial  ^-  tang
  %+  expect-eq
    !>(`@ud`4)
    !>(out:(binomial:rs:dist:i754rand (from-atom:seed:rand %sm64 0) 20 .0.3))
++  test-rs-dirichlet  ^-  tang
  =/  out  `(list @rs)`out:(dirichlet:rs:dist:i754rand (from-atom:seed:rand %sm64 0) ~[.1 .2 .3])
  =/  mth  ~(. rs:math [%n .1e-6 .0])
  =/  total  (add:mth (snag 0 out) (add:mth (snag 1 out) (snag 2 out)))
  %+  expect-eq  !>(%.y)  !>((is-close:mth total .1))
--
