::  bench-cells: the per-(impl,door,arm) timing cell as a library arm, so both
::  gen/bench-math (one cell) and gen/bench-grid (the whole grid in one dojo
::  invocation) call it.  +cell builds the PRECOMPUTED input list (sun:rd kept
::  out of the hot path) then folds the arm over it via (time ..), whose
::  ~>(%bout) slogs "took ..".
::
/+  *bench-core, bd=bench-domains, math, math-taylor
|%
++  cell
  |=  [impl=?(%cheb %taylor) door=@tas arm=@tas n=@ud]
=/  dm   (get:bd ?:(=(arm %base) %exp arm) ?:(=(impl %taylor) %taylor %cheb) door)
=/  tay  =(impl %taylor)
^-  @
?+  door  ~|([%bad-door door] !!)
::
  %rd
    =/  den  (sun:rd den.dm)
    =/  lof  (div:rd (san:rd lo.dm) den)
    =/  wid  (sub:rd (div:rd (san:rd hi.dm) den) lof)
    =/  spn  (sun:rd span.dm)
    =/  inp  |=(k=@ud ^-(@rd (add:rd lof (mul:rd (div:rd (sun:rd (mod k span.dm)) spn) wid))))
    =/  xs=(list [x=@ y=@])
      %+  turn  (gulf 0 (dec n))
      |=  k=@ud  ^-  [@ @]
      ?+  arm  [(inp k) 0]
        %atan2  [(inp k) (inp +(k))]
        %pow    [(inp k) .~2.5]
        %pow-n  [(inp k) (sun:rd (add 2 (mod k 7)))]
      ==
    =/  step
      |=  [p=[x=@ y=@] acc=@]  ^-  @
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x.p
        %exp     ?:(tay (~(exp rd:math-taylor [%n .~1e-10]) x.p) (~(exp rd:math [%n .~1e-10]) x.p))
        %log     ?:(tay (~(log rd:math-taylor [%n .~1e-10]) x.p) (~(log rd:math [%n .~1e-10]) x.p))
        %sin     ?:(tay (~(sin rd:math-taylor [%n .~1e-10]) x.p) (~(sin rd:math [%n .~1e-10]) x.p))
        %cos     ?:(tay (~(cos rd:math-taylor [%n .~1e-10]) x.p) (~(cos rd:math [%n .~1e-10]) x.p))
        %tan     ?:(tay (~(tan rd:math-taylor [%n .~1e-10]) x.p) (~(tan rd:math [%n .~1e-10]) x.p))
        %atan    ?:(tay (~(atan rd:math-taylor [%n .~1e-10]) x.p) (~(atan rd:math [%n .~1e-10]) x.p))
        %asin    ?:(tay (~(asin rd:math-taylor [%n .~1e-10]) x.p) (~(asin rd:math [%n .~1e-10]) x.p))
        %acos    ?:(tay (~(acos rd:math-taylor [%n .~1e-10]) x.p) (~(acos rd:math [%n .~1e-10]) x.p))
        %sqt     ?:(tay (~(sqt rd:math-taylor [%n .~1e-10]) x.p) (~(sqt rd:math [%n .~1e-10]) x.p))
        %cbt     ?:(tay (~(cbt rd:math-taylor [%n .~1e-10]) x.p) (~(cbt rd:math [%n .~1e-10]) x.p))
        %log-2   ?:(tay (~(log-2 rd:math-taylor [%n .~1e-10]) x.p) (~(log-2 rd:math [%n .~1e-10]) x.p))
        %log-10  ?:(tay (~(log-10 rd:math-taylor [%n .~1e-10]) x.p) (~(log-10 rd:math [%n .~1e-10]) x.p))
        %atan2   ?:(tay (~(atan2 rd:math-taylor [%n .~1e-10]) x.p y.p) (~(atan2 rd:math [%n .~1e-10]) x.p y.p))
        %pow     ?:(tay (~(pow rd:math-taylor [%n .~1e-10]) x.p y.p) (~(pow rd:math [%n .~1e-10]) x.p y.p))
        %pow-n   ?:(tay (~(pow-n rd:math-taylor [%n .~1e-10]) x.p y.p) (~(pow-n rd:math [%n .~1e-10]) x.p y.p))
      ==
    (time xs step)
::
  %rs
    =/  den  (sun:rs den.dm)
    =/  lof  (div:rs (san:rs lo.dm) den)
    =/  wid  (sub:rs (div:rs (san:rs hi.dm) den) lof)
    =/  spn  (sun:rs span.dm)
    =/  inp  |=(k=@ud ^-(@rs (add:rs lof (mul:rs (div:rs (sun:rs (mod k span.dm)) spn) wid))))
    =/  xs=(list [x=@ y=@])
      %+  turn  (gulf 0 (dec n))
      |=  k=@ud  ^-  [@ @]
      ?+  arm  [(inp k) 0]
        %atan2  [(inp k) (inp +(k))]
        %pow    [(inp k) .2.5]
        %pow-n  [(inp k) (sun:rs (add 2 (mod k 7)))]
      ==
    =/  step
      |=  [p=[x=@ y=@] acc=@]  ^-  @
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x.p
        %exp     ?:(tay (~(exp rs:math-taylor [%n .1e-5]) x.p) (~(exp rs:math [%n .1e-5]) x.p))
        %log     ?:(tay (~(log rs:math-taylor [%n .1e-5]) x.p) (~(log rs:math [%n .1e-5]) x.p))
        %sin     ?:(tay (~(sin rs:math-taylor [%n .1e-5]) x.p) (~(sin rs:math [%n .1e-5]) x.p))
        %cos     ?:(tay (~(cos rs:math-taylor [%n .1e-5]) x.p) (~(cos rs:math [%n .1e-5]) x.p))
        %tan     ?:(tay (~(tan rs:math-taylor [%n .1e-5]) x.p) (~(tan rs:math [%n .1e-5]) x.p))
        %atan    ?:(tay (~(atan rs:math-taylor [%n .1e-5]) x.p) (~(atan rs:math [%n .1e-5]) x.p))
        %asin    ?:(tay (~(asin rs:math-taylor [%n .1e-5]) x.p) (~(asin rs:math [%n .1e-5]) x.p))
        %acos    ?:(tay (~(acos rs:math-taylor [%n .1e-5]) x.p) (~(acos rs:math [%n .1e-5]) x.p))
        %sqt     ?:(tay (~(sqt rs:math-taylor [%n .1e-5]) x.p) (~(sqt rs:math [%n .1e-5]) x.p))
        %cbt     ?:(tay (~(cbt rs:math-taylor [%n .1e-5]) x.p) (~(cbt rs:math [%n .1e-5]) x.p))
        %log-2   ?:(tay (~(log-2 rs:math-taylor [%n .1e-5]) x.p) (~(log-2 rs:math [%n .1e-5]) x.p))
        %log-10  ?:(tay (~(log-10 rs:math-taylor [%n .1e-5]) x.p) (~(log-10 rs:math [%n .1e-5]) x.p))
        %atan2   ?:(tay (~(atan2 rs:math-taylor [%n .1e-5]) x.p y.p) (~(atan2 rs:math [%n .1e-5]) x.p y.p))
        %pow     ?:(tay (~(pow rs:math-taylor [%n .1e-5]) x.p y.p) (~(pow rs:math [%n .1e-5]) x.p y.p))
        %pow-n   ?:(tay (~(pow-n rs:math-taylor [%n .1e-5]) x.p y.p) (~(pow-n rs:math [%n .1e-5]) x.p y.p))
      ==
    (time xs step)
::
  %rh
    =/  den  (sun:rh den.dm)
    =/  lof  (div:rh (san:rh lo.dm) den)
    =/  wid  (sub:rh (div:rh (san:rh hi.dm) den) lof)
    =/  spn  (sun:rh span.dm)
    =/  inp  |=(k=@ud ^-(@rh (add:rh lof (mul:rh (div:rh (sun:rh (mod k span.dm)) spn) wid))))
    =/  xs=(list [x=@ y=@])
      %+  turn  (gulf 0 (dec n))
      |=  k=@ud  ^-  [@ @]
      ?+  arm  [(inp k) 0]
        %atan2  [(inp k) (inp +(k))]
        %pow    [(inp k) .~~2.5]
        %pow-n  [(inp k) (sun:rh (add 2 (mod k 7)))]
      ==
    =/  step
      |=  [p=[x=@ y=@] acc=@]  ^-  @
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x.p
        %exp     ?:(tay (~(exp rh:math-taylor [%n .~~1e-2]) x.p) (~(exp rh:math [%n .~~1e-2]) x.p))
        %log     ?:(tay (~(log rh:math-taylor [%n .~~1e-2]) x.p) (~(log rh:math [%n .~~1e-2]) x.p))
        %sin     ?:(tay (~(sin rh:math-taylor [%n .~~1e-2]) x.p) (~(sin rh:math [%n .~~1e-2]) x.p))
        %cos     ?:(tay (~(cos rh:math-taylor [%n .~~1e-2]) x.p) (~(cos rh:math [%n .~~1e-2]) x.p))
        %tan     ?:(tay (~(tan rh:math-taylor [%n .~~1e-2]) x.p) (~(tan rh:math [%n .~~1e-2]) x.p))
        %atan    ?:(tay (~(atan rh:math-taylor [%n .~~1e-2]) x.p) (~(atan rh:math [%n .~~1e-2]) x.p))
        %asin    ?:(tay (~(asin rh:math-taylor [%n .~~1e-2]) x.p) (~(asin rh:math [%n .~~1e-2]) x.p))
        %acos    ?:(tay (~(acos rh:math-taylor [%n .~~1e-2]) x.p) (~(acos rh:math [%n .~~1e-2]) x.p))
        %sqt     ?:(tay (~(sqt rh:math-taylor [%n .~~1e-2]) x.p) (~(sqt rh:math [%n .~~1e-2]) x.p))
        %cbt     ?:(tay (~(cbt rh:math-taylor [%n .~~1e-2]) x.p) (~(cbt rh:math [%n .~~1e-2]) x.p))
        %log-2   ?:(tay (~(log-2 rh:math-taylor [%n .~~1e-2]) x.p) (~(log-2 rh:math [%n .~~1e-2]) x.p))
        %log-10  ?:(tay (~(log-10 rh:math-taylor [%n .~~1e-2]) x.p) (~(log-10 rh:math [%n .~~1e-2]) x.p))
        %atan2   ?:(tay (~(atan2 rh:math-taylor [%n .~~1e-2]) x.p y.p) (~(atan2 rh:math [%n .~~1e-2]) x.p y.p))
        %pow     ?:(tay (~(pow rh:math-taylor [%n .~~1e-2]) x.p y.p) (~(pow rh:math [%n .~~1e-2]) x.p y.p))
        %pow-n   ?:(tay (~(pow-n rh:math-taylor [%n .~~1e-2]) x.p y.p) (~(pow-n rh:math [%n .~~1e-2]) x.p y.p))
      ==
    (time xs step)
::
  %rq
    =/  den  (sun:rq den.dm)
    =/  lof  (div:rq (san:rq lo.dm) den)
    =/  wid  (sub:rq (div:rq (san:rq hi.dm) den) lof)
    =/  spn  (sun:rq span.dm)
    =/  inp  |=(k=@ud ^-(@rq (add:rq lof (mul:rq (div:rq (sun:rq (mod k span.dm)) spn) wid))))
    =/  xs=(list [x=@ y=@])
      %+  turn  (gulf 0 (dec n))
      |=  k=@ud  ^-  [@ @]
      ?+  arm  [(inp k) 0]
        %atan2  [(inp k) (inp +(k))]
        %pow    [(inp k) .~~~2.5]
        %pow-n  [(inp k) (sun:rq (add 2 (mod k 7)))]
      ==
    =/  step
      |=  [p=[x=@ y=@] acc=@]  ^-  @
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x.p
        %exp     ?:(tay (~(exp rq:math-taylor [%n .~~~1e-10]) x.p) (~(exp rq:math [%n .~~~1e-10]) x.p))
        %log     ?:(tay (~(log rq:math-taylor [%n .~~~1e-10]) x.p) (~(log rq:math [%n .~~~1e-10]) x.p))
        %sin     ?:(tay (~(sin rq:math-taylor [%n .~~~1e-10]) x.p) (~(sin rq:math [%n .~~~1e-10]) x.p))
        %cos     ?:(tay (~(cos rq:math-taylor [%n .~~~1e-10]) x.p) (~(cos rq:math [%n .~~~1e-10]) x.p))
        %tan     ?:(tay (~(tan rq:math-taylor [%n .~~~1e-10]) x.p) (~(tan rq:math [%n .~~~1e-10]) x.p))
        %atan    ?:(tay (~(atan rq:math-taylor [%n .~~~1e-10]) x.p) (~(atan rq:math [%n .~~~1e-10]) x.p))
        %asin    ?:(tay (~(asin rq:math-taylor [%n .~~~1e-10]) x.p) (~(asin rq:math [%n .~~~1e-10]) x.p))
        %acos    ?:(tay (~(acos rq:math-taylor [%n .~~~1e-10]) x.p) (~(acos rq:math [%n .~~~1e-10]) x.p))
        %sqt     ?:(tay (~(sqt rq:math-taylor [%n .~~~1e-10]) x.p) (~(sqt rq:math [%n .~~~1e-10]) x.p))
        %cbt     ?:(tay (~(cbt rq:math-taylor [%n .~~~1e-10]) x.p) (~(cbt rq:math [%n .~~~1e-10]) x.p))
        %log-2   ?:(tay (~(log-2 rq:math-taylor [%n .~~~1e-10]) x.p) (~(log-2 rq:math [%n .~~~1e-10]) x.p))
        %log-10  ?:(tay (~(log-10 rq:math-taylor [%n .~~~1e-10]) x.p) (~(log-10 rq:math [%n .~~~1e-10]) x.p))
        %atan2   ?:(tay (~(atan2 rq:math-taylor [%n .~~~1e-10]) x.p y.p) (~(atan2 rq:math [%n .~~~1e-10]) x.p y.p))
        %pow     ?:(tay (~(pow rq:math-taylor [%n .~~~1e-10]) x.p y.p) (~(pow rq:math [%n .~~~1e-10]) x.p y.p))
        %pow-n   ?:(tay (~(pow-n rq:math-taylor [%n .~~~1e-10]) x.p y.p) (~(pow-n rq:math [%n .~~~1e-10]) x.p y.p))
      ==
    (time xs step)
==
--
