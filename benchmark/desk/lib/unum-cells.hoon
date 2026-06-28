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
  ::  inputs in [0.25, 1.25]: positive, in-range for exp/log/sqt and all ops.
  =/  inp  |=(k=@ud ^-(@ (div:d (sun:d +((mod k 5))) (sun:d 4))))
  =/  xs=(list [x=@ y=@])
    %+  turn  (gulf 0 (dec n))
    |=  k=@ud  ^-([@ @] [(inp k) (inp +(k))])
  =/  step
    |=  [p=[x=@ y=@] acc=@]  ^-  @
    %+  add  acc
    ?+  arm  ~|([%bad-arm arm] !!)
      %base    x.p
      %neg     (neg:d x.p)
      %abs     (abs:d x.p)
      %sqt     (sqt:d x.p)
      %exp     (exp:d x.p)
      %log     (log:d x.p)
      %sin     (sin:d x.p)
      %cos     (cos:d x.p)
      %atan    (atan:d x.p)
      %add     (add:d x.p y.p)
      %sub     (sub:d x.p y.p)
      %mul     (mul:d x.p y.p)
      %div     (div:d x.p y.p)
      %fma     (fma:d x.p y.p x.p)
      %pow     (pow:d x.p y.p)
      %lth     ?:((lth:d x.p y.p) 1 0)
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
