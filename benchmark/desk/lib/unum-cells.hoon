::  unum-cells: per-(door,arm) timing cell for /lib/unum (posit) benchmarks,
::  mirroring lib/bench-cells (the /lib/math harness).  +cell precomputes the
::  posit input list (sun/div kept OUT of the hot path) and folds the arm over
::  it via (time ..), whose ~>(%bout) slogs "took ..".  Per-call cost =
::  (cell(arm) - cell(%base)) / n.  +fdp-cell times one fused dot product over
::  length-n vectors (per-element = took / n).
::
::  The posit width is selected at runtime by %*-specializing the generic ++pp
::  core on bloq (3/4/5 = posit8/16/32); the unum jets read bloq from that
::  sample, so this one body benchmarks every width and both the jetted (hints
::  on) and interpreted (hints commented) builds.
::
/+  *bench-core, unum
|%
::    +dor:  door @tas -> bloq
++  dor  |=(door=@tas ?:(=(door %rpb) 3 ?:(=(door %rph) 4 5)))
::    +cell:  time `n` folds of `arm` at posit width `door`.
++  cell
  |=  [door=@tas arm=@tas n=@ud]
  ^-  @
  =/  d  %*(. pp:unum bloq (dor door))
  ::  inputs in [0.25, 1.25]: positive, in-range for exp/log/sqt/tan/cbrt/
  ::  log-2/log-10 and all the binary ops.  asin/acos need |x|<1 (1.25 would
  ::  be an instant NaR short-circuit, skewing the timing low) -- +inp-ainv
  ::  below is their domain-safe counterpart, in (-1,1).
  =/  inp  |=(k=@ud ^-(@ (div:d (sun:d +((mod k 5))) (sun:d 4))))
  ::  asin/acos need |x|<1 -- same cycling-small-integers style as +inp, just
  ::  denominator 8 instead of 4 so all 5 values (.125/.25/.375/.5/.625) stay
  ::  safely interior (1.25 from +inp would be an instant NaR short-circuit
  ::  for these two arms, skewing the timing low).
  =/  inp-ainv  |=(k=@ud ^-(@ (div:d (sun:d +((mod k 5))) (sun:d 8))))
  =/  safe  ?:(|(=(arm %asin) =(arm %acos)) inp-ainv inp)
  =/  xs=(list [x=@ y=@])
    %+  turn  (gulf 0 (dec n))
    |=  k=@ud  ^-([@ @] [(safe k) (safe +(k))])
  =/  step
    |=  [p=[x=@ y=@] acc=@]  ^-  @
    %+  add  acc
    ?+  arm  ~|([%bad-arm arm] !!)
      %base     x.p
      %neg      (neg:d x.p)
      %abs      (abs:d x.p)
      %sqt      (sqt:d x.p)
      %exp      (exp:d x.p)
      %log      (log:d x.p)
      %log-2    (log-2:d x.p)
      %log-10   (log-10:d x.p)
      %sin      (sin:d x.p)
      %cos      (cos:d x.p)
      %tan      (tan:d x.p)
      %atan     (atan:d x.p)
      %asin     (asin:d x.p)
      %acos     (acos:d x.p)
      %cbrt     (cbrt:d x.p)
      %pow-n    (pow-n:d x.p 3)
      %add      (add:d x.p y.p)
      %sub      (sub:d x.p y.p)
      %mul      (mul:d x.p y.p)
      %div      (div:d x.p y.p)
      %fma      (fma:d x.p y.p x.p)
      %pow      (pow:d x.p y.p)
      %lth      ?:((lth:d x.p y.p) 1 0)
    ==
  (time xs step)
::    +fdp-cell:  time one fused dot product over length-n posit vectors.
++  fdp-cell
  |=  [door=@tas n=@ud]
  ^-  @
  =/  d  %*(. pp:unum bloq (dor door))
  =/  inp  |=(k=@ud ^-(@ (div:d (sun:d +((mod k 5))) (sun:d 4))))
  =/  av=(list @)  (turn (gulf 0 (dec n)) inp)
  =/  bv=(list @)  (turn (gulf 0 (dec n)) |=(k=@ud (inp +(k))))
  ~>  %bout
  (fdp:d av bv)
--
