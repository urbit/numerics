::  bench-math: run ONE timing cell.  +bench-math [impl door arm n]
::    impl  ?(%cheb %taylor)   (jetted = %cheb on a jet binary)
::    door  ?(%rh %rs %rd %rq)
::    arm   @tas               (a transcendental, or %base for the baseline loop)
::    n     @ud                iterations
::  Returns the folded accumulator; ~>(%bout ..) prints elapsed (driver scrapes it).
::  Per-call = (cell(arm) - cell(%base)) / n.
::
/+  *bench-core, bd=bench-domains, math, math-taylor
:-  %say
|=  [* [impl=?(%cheb %taylor) door=?(%rh %rs %rd %rq) arm=@tas n=@ud ~] ~]
:-  %noun
=/  dm   (get:bd ?:(=(arm %base) %exp arm) ?:(=(impl %taylor) %taylor %cheb))
=/  tay  =(impl %taylor)
^-  @
?-  door
::
  %rd
    =/  inp
      |=  i=@ud  ^-  @rd
      =/  lof  (div:rd (san:rd lo.dm) (sun:rd den.dm))
      =/  hif  (div:rd (san:rd hi.dm) (sun:rd den.dm))
      =/  k    (div:rd (sun:rd (mod i span.dm)) (sun:rd span.dm))
      (add:rd lof (mul:rd k (sub:rd hif lof)))
    =/  step
      |=  [i=@ud acc=@]  ^-  @
      =/  x  (inp i)
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x
        %exp     ?:(tay (~(exp rd:math-taylor [%n .~1e-10]) x) (~(exp rd:math [%n .~1e-10]) x))
        %log     ?:(tay (~(log rd:math-taylor [%n .~1e-10]) x) (~(log rd:math [%n .~1e-10]) x))
        %sin     ?:(tay (~(sin rd:math-taylor [%n .~1e-10]) x) (~(sin rd:math [%n .~1e-10]) x))
        %cos     ?:(tay (~(cos rd:math-taylor [%n .~1e-10]) x) (~(cos rd:math [%n .~1e-10]) x))
        %tan     ?:(tay (~(tan rd:math-taylor [%n .~1e-10]) x) (~(tan rd:math [%n .~1e-10]) x))
        %atan    ?:(tay (~(atan rd:math-taylor [%n .~1e-10]) x) (~(atan rd:math [%n .~1e-10]) x))
        %asin    ?:(tay (~(asin rd:math-taylor [%n .~1e-10]) x) (~(asin rd:math [%n .~1e-10]) x))
        %acos    ?:(tay (~(acos rd:math-taylor [%n .~1e-10]) x) (~(acos rd:math [%n .~1e-10]) x))
        %sqt     ?:(tay (~(sqt rd:math-taylor [%n .~1e-10]) x) (~(sqt rd:math [%n .~1e-10]) x))
        %cbt     ?:(tay (~(cbt rd:math-taylor [%n .~1e-10]) x) (~(cbt rd:math [%n .~1e-10]) x))
        %log-2   ?:(tay (~(log-2 rd:math-taylor [%n .~1e-10]) x) (~(log-2 rd:math [%n .~1e-10]) x))
        %log-10  ?:(tay (~(log-10 rd:math-taylor [%n .~1e-10]) x) (~(log-10 rd:math [%n .~1e-10]) x))
        %atan2   =/  y  (inp +(i))
                 ?:(tay (~(atan2 rd:math-taylor [%n .~1e-10]) x y) (~(atan2 rd:math [%n .~1e-10]) x y))
        %pow     ?:(tay (~(pow rd:math-taylor [%n .~1e-10]) x .~2.5) (~(pow rd:math [%n .~1e-10]) x .~2.5))
        %pow-n   =/  ne  (sun:rd (add 2 (mod i 7)))
                 ?:(tay (~(pow-n rd:math-taylor [%n .~1e-10]) x ne) (~(pow-n rd:math [%n .~1e-10]) x ne))
      ==
    (time n step)
::
  %rs
    =/  inp
      |=  i=@ud  ^-  @rs
      =/  lof  (div:rs (san:rs lo.dm) (sun:rs den.dm))
      =/  hif  (div:rs (san:rs hi.dm) (sun:rs den.dm))
      =/  k    (div:rs (sun:rs (mod i span.dm)) (sun:rs span.dm))
      (add:rs lof (mul:rs k (sub:rs hif lof)))
    =/  step
      |=  [i=@ud acc=@]  ^-  @
      =/  x  (inp i)
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x
        %exp     ?:(tay (~(exp rs:math-taylor [%n .1e-5]) x) (~(exp rs:math [%n .1e-5]) x))
        %log     ?:(tay (~(log rs:math-taylor [%n .1e-5]) x) (~(log rs:math [%n .1e-5]) x))
        %sin     ?:(tay (~(sin rs:math-taylor [%n .1e-5]) x) (~(sin rs:math [%n .1e-5]) x))
        %cos     ?:(tay (~(cos rs:math-taylor [%n .1e-5]) x) (~(cos rs:math [%n .1e-5]) x))
        %tan     ?:(tay (~(tan rs:math-taylor [%n .1e-5]) x) (~(tan rs:math [%n .1e-5]) x))
        %atan    ?:(tay (~(atan rs:math-taylor [%n .1e-5]) x) (~(atan rs:math [%n .1e-5]) x))
        %asin    ?:(tay (~(asin rs:math-taylor [%n .1e-5]) x) (~(asin rs:math [%n .1e-5]) x))
        %acos    ?:(tay (~(acos rs:math-taylor [%n .1e-5]) x) (~(acos rs:math [%n .1e-5]) x))
        %sqt     ?:(tay (~(sqt rs:math-taylor [%n .1e-5]) x) (~(sqt rs:math [%n .1e-5]) x))
        %cbt     ?:(tay (~(cbt rs:math-taylor [%n .1e-5]) x) (~(cbt rs:math [%n .1e-5]) x))
        %log-2   ?:(tay (~(log-2 rs:math-taylor [%n .1e-5]) x) (~(log-2 rs:math [%n .1e-5]) x))
        %log-10  ?:(tay (~(log-10 rs:math-taylor [%n .1e-5]) x) (~(log-10 rs:math [%n .1e-5]) x))
        %atan2   =/  y  (inp +(i))
                 ?:(tay (~(atan2 rs:math-taylor [%n .1e-5]) x y) (~(atan2 rs:math [%n .1e-5]) x y))
        %pow     ?:(tay (~(pow rs:math-taylor [%n .1e-5]) x .2.5) (~(pow rs:math [%n .1e-5]) x .2.5))
        %pow-n   =/  ne  (sun:rs (add 2 (mod i 7)))
                 ?:(tay (~(pow-n rs:math-taylor [%n .1e-5]) x ne) (~(pow-n rs:math [%n .1e-5]) x ne))
      ==
    (time n step)
::
  %rh
    =/  inp
      |=  i=@ud  ^-  @rh
      =/  lof  (div:rh (san:rh lo.dm) (sun:rh den.dm))
      =/  hif  (div:rh (san:rh hi.dm) (sun:rh den.dm))
      =/  k    (div:rh (sun:rh (mod i span.dm)) (sun:rh span.dm))
      (add:rh lof (mul:rh k (sub:rh hif lof)))
    =/  step
      |=  [i=@ud acc=@]  ^-  @
      =/  x  (inp i)
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x
        %exp     ?:(tay (~(exp rh:math-taylor [%n .~~1e-2]) x) (~(exp rh:math [%n .~~1e-2]) x))
        %log     ?:(tay (~(log rh:math-taylor [%n .~~1e-2]) x) (~(log rh:math [%n .~~1e-2]) x))
        %sin     ?:(tay (~(sin rh:math-taylor [%n .~~1e-2]) x) (~(sin rh:math [%n .~~1e-2]) x))
        %cos     ?:(tay (~(cos rh:math-taylor [%n .~~1e-2]) x) (~(cos rh:math [%n .~~1e-2]) x))
        %tan     ?:(tay (~(tan rh:math-taylor [%n .~~1e-2]) x) (~(tan rh:math [%n .~~1e-2]) x))
        %atan    ?:(tay (~(atan rh:math-taylor [%n .~~1e-2]) x) (~(atan rh:math [%n .~~1e-2]) x))
        %asin    ?:(tay (~(asin rh:math-taylor [%n .~~1e-2]) x) (~(asin rh:math [%n .~~1e-2]) x))
        %acos    ?:(tay (~(acos rh:math-taylor [%n .~~1e-2]) x) (~(acos rh:math [%n .~~1e-2]) x))
        %sqt     ?:(tay (~(sqt rh:math-taylor [%n .~~1e-2]) x) (~(sqt rh:math [%n .~~1e-2]) x))
        %cbt     ?:(tay (~(cbt rh:math-taylor [%n .~~1e-2]) x) (~(cbt rh:math [%n .~~1e-2]) x))
        %log-2   ?:(tay (~(log-2 rh:math-taylor [%n .~~1e-2]) x) (~(log-2 rh:math [%n .~~1e-2]) x))
        %log-10  ?:(tay (~(log-10 rh:math-taylor [%n .~~1e-2]) x) (~(log-10 rh:math [%n .~~1e-2]) x))
        %atan2   =/  y  (inp +(i))
                 ?:(tay (~(atan2 rh:math-taylor [%n .~~1e-2]) x y) (~(atan2 rh:math [%n .~~1e-2]) x y))
        %pow     ?:(tay (~(pow rh:math-taylor [%n .~~1e-2]) x .~~2.5) (~(pow rh:math [%n .~~1e-2]) x .~~2.5))
        %pow-n   =/  ne  (sun:rh (add 2 (mod i 7)))
                 ?:(tay (~(pow-n rh:math-taylor [%n .~~1e-2]) x ne) (~(pow-n rh:math [%n .~~1e-2]) x ne))
      ==
    (time n step)
::
  %rq
    =/  inp
      |=  i=@ud  ^-  @rq
      =/  lof  (div:rq (san:rq lo.dm) (sun:rq den.dm))
      =/  hif  (div:rq (san:rq hi.dm) (sun:rq den.dm))
      =/  k    (div:rq (sun:rq (mod i span.dm)) (sun:rq span.dm))
      (add:rq lof (mul:rq k (sub:rq hif lof)))
    =/  step
      |=  [i=@ud acc=@]  ^-  @
      =/  x  (inp i)
      %+  add  acc
      ?+  arm  ~|([%bad-arm arm] !!)
        %base    x
        %exp     ?:(tay (~(exp rq:math-taylor [%n .~~~1e-10]) x) (~(exp rq:math [%n .~~~1e-10]) x))
        %log     ?:(tay (~(log rq:math-taylor [%n .~~~1e-10]) x) (~(log rq:math [%n .~~~1e-10]) x))
        %sin     ?:(tay (~(sin rq:math-taylor [%n .~~~1e-10]) x) (~(sin rq:math [%n .~~~1e-10]) x))
        %cos     ?:(tay (~(cos rq:math-taylor [%n .~~~1e-10]) x) (~(cos rq:math [%n .~~~1e-10]) x))
        %tan     ?:(tay (~(tan rq:math-taylor [%n .~~~1e-10]) x) (~(tan rq:math [%n .~~~1e-10]) x))
        %atan    ?:(tay (~(atan rq:math-taylor [%n .~~~1e-10]) x) (~(atan rq:math [%n .~~~1e-10]) x))
        %asin    ?:(tay (~(asin rq:math-taylor [%n .~~~1e-10]) x) (~(asin rq:math [%n .~~~1e-10]) x))
        %acos    ?:(tay (~(acos rq:math-taylor [%n .~~~1e-10]) x) (~(acos rq:math [%n .~~~1e-10]) x))
        %sqt     ?:(tay (~(sqt rq:math-taylor [%n .~~~1e-10]) x) (~(sqt rq:math [%n .~~~1e-10]) x))
        %cbt     ?:(tay (~(cbt rq:math-taylor [%n .~~~1e-10]) x) (~(cbt rq:math [%n .~~~1e-10]) x))
        %log-2   ?:(tay (~(log-2 rq:math-taylor [%n .~~~1e-10]) x) (~(log-2 rq:math [%n .~~~1e-10]) x))
        %log-10  ?:(tay (~(log-10 rq:math-taylor [%n .~~~1e-10]) x) (~(log-10 rq:math [%n .~~~1e-10]) x))
        %atan2   =/  y  (inp +(i))
                 ?:(tay (~(atan2 rq:math-taylor [%n .~~~1e-10]) x y) (~(atan2 rq:math [%n .~~~1e-10]) x y))
        %pow     ?:(tay (~(pow rq:math-taylor [%n .~~~1e-10]) x .~~~2.5) (~(pow rq:math [%n .~~~1e-10]) x .~~~2.5))
        %pow-n   =/  ne  (sun:rq (add 2 (mod i 7)))
                 ?:(tay (~(pow-n rq:math-taylor [%n .~~~1e-10]) x ne) (~(pow-n rq:math [%n .~~~1e-10]) x ne))
      ==
    (time n step)
==
