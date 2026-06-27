::  bench-acc: per-(door,arm) accuracy sweep.  +bench-acc [door arm pts]
::    door  ?(%rh %rs %rd %rq)
::    arm   @tas   (a single-argument transcendental)
::    pts   @ud    number of points swept across the arm's (cheb) domain
::  Returns [max=@ud mean=@ud] ULP distance between the Chebyshev result and the
::  Taylor result at each point -- an on-ship ALGORITHM-AGREEMENT metric (both are
::  approximations, so this bounds their mutual disagreement, not absolute error).
::  GROUND-TRUTH ULP vs a high-precision MPFR oracle is produced by tools/rq_check.c.
::  Two-argument arms (atan2/pow/pow-n) report 0 here -- TODO, see rq_check.c.
::
::  ULP distance uses the standard monotone map of the IEEE bit-pattern: positive
::  values get the sign bit set, negative values are bit-flipped, so adjacent
::  representable floats differ by 1.
::
/+  bd=bench-domains, math, math-taylor
:-  %say
|=  [* [door=?(%rh %rs %rd %rq) arm=@tas pts=@ud ~] ~]
:-  %noun
^-  [max=@ud mean=@ud]
=/  bw    ?-(door %rh 16, %rs 32, %rd 64, %rq 128)
=/  top   (bex (dec bw))
=/  alo   (dec (bex bw))
=/  mono  |=(b=@ ^-(@ ?:(=(0 (dis b top)) (con b top) (mix b alo))))
=/  ulp   |=([a=@ b=@] ^-(@ =/(x (mono a) =/(y (mono b) ?:((gth x y) (sub x y) (sub y x))))))
=/  dm    (get:bd arm %cheb door)
=/  sweep
  |=  [inp=$-(@ud @) res=$-(@ [c=@ t=@])]  ^-  [@ud @ud]
  =/  i=@ud  0
  =/  mx=@   0
  =/  sm=@   0
  |-  ^-  [@ud @ud]
  ?:  =(i pts)  [mx ?:(=(0 pts) 0 (div sm pts))]
  =/  ct  (res (inp i))
  =/  u   (ulp c.ct t.ct)
  $(i +(i), mx (max mx u), sm (add sm u))
?-  door
::
  %rd
    =/  den  (sun:rd den.dm)
    =/  lof  (div:rd (san:rd lo.dm) den)
    =/  wd   (sub:rd (div:rd (san:rd hi.dm) den) lof)
    =/  inp  |=(k=@ud ^-(@ (add:rd lof (mul:rd (div:rd (sun:rd k) (sun:rd pts)) wd))))
    %+  sweep  inp
    |=  x=@  ^-  [c=@ t=@]
    ?+  arm  [x x]
      %exp     [(~(exp rd:math [%n .~1e-10]) x) (~(exp rd:math-taylor [%n .~1e-10]) x)]
      %log     [(~(log rd:math [%n .~1e-10]) x) (~(log rd:math-taylor [%n .~1e-10]) x)]
      %sin     [(~(sin rd:math [%n .~1e-10]) x) (~(sin rd:math-taylor [%n .~1e-10]) x)]
      %cos     [(~(cos rd:math [%n .~1e-10]) x) (~(cos rd:math-taylor [%n .~1e-10]) x)]
      %tan     [(~(tan rd:math [%n .~1e-10]) x) (~(tan rd:math-taylor [%n .~1e-10]) x)]
      %atan    [(~(atan rd:math [%n .~1e-10]) x) (~(atan rd:math-taylor [%n .~1e-10]) x)]
      %asin    [(~(asin rd:math [%n .~1e-10]) x) (~(asin rd:math-taylor [%n .~1e-10]) x)]
      %acos    [(~(acos rd:math [%n .~1e-10]) x) (~(acos rd:math-taylor [%n .~1e-10]) x)]
      %sqt     [(~(sqt rd:math [%n .~1e-10]) x) (~(sqt rd:math-taylor [%n .~1e-10]) x)]
      %cbt     [(~(cbt rd:math [%n .~1e-10]) x) (~(cbt rd:math-taylor [%n .~1e-10]) x)]
      %log-2   [(~(log-2 rd:math [%n .~1e-10]) x) (~(log-2 rd:math-taylor [%n .~1e-10]) x)]
      %log-10  [(~(log-10 rd:math [%n .~1e-10]) x) (~(log-10 rd:math-taylor [%n .~1e-10]) x)]
    ==
::
  %rs
    =/  den  (sun:rs den.dm)
    =/  lof  (div:rs (san:rs lo.dm) den)
    =/  wd   (sub:rs (div:rs (san:rs hi.dm) den) lof)
    =/  inp  |=(k=@ud ^-(@ (add:rs lof (mul:rs (div:rs (sun:rs k) (sun:rs pts)) wd))))
    %+  sweep  inp
    |=  x=@  ^-  [c=@ t=@]
    ?+  arm  [x x]
      %exp     [(~(exp rs:math [%n .1e-5]) x) (~(exp rs:math-taylor [%n .1e-5]) x)]
      %log     [(~(log rs:math [%n .1e-5]) x) (~(log rs:math-taylor [%n .1e-5]) x)]
      %sin     [(~(sin rs:math [%n .1e-5]) x) (~(sin rs:math-taylor [%n .1e-5]) x)]
      %cos     [(~(cos rs:math [%n .1e-5]) x) (~(cos rs:math-taylor [%n .1e-5]) x)]
      %tan     [(~(tan rs:math [%n .1e-5]) x) (~(tan rs:math-taylor [%n .1e-5]) x)]
      %atan    [(~(atan rs:math [%n .1e-5]) x) (~(atan rs:math-taylor [%n .1e-5]) x)]
      %asin    [(~(asin rs:math [%n .1e-5]) x) (~(asin rs:math-taylor [%n .1e-5]) x)]
      %acos    [(~(acos rs:math [%n .1e-5]) x) (~(acos rs:math-taylor [%n .1e-5]) x)]
      %sqt     [(~(sqt rs:math [%n .1e-5]) x) (~(sqt rs:math-taylor [%n .1e-5]) x)]
      %cbt     [(~(cbt rs:math [%n .1e-5]) x) (~(cbt rs:math-taylor [%n .1e-5]) x)]
      %log-2   [(~(log-2 rs:math [%n .1e-5]) x) (~(log-2 rs:math-taylor [%n .1e-5]) x)]
      %log-10  [(~(log-10 rs:math [%n .1e-5]) x) (~(log-10 rs:math-taylor [%n .1e-5]) x)]
    ==
::
  %rh
    =/  den  (sun:rh den.dm)
    =/  lof  (div:rh (san:rh lo.dm) den)
    =/  wd   (sub:rh (div:rh (san:rh hi.dm) den) lof)
    =/  inp  |=(k=@ud ^-(@ (add:rh lof (mul:rh (div:rh (sun:rh k) (sun:rh pts)) wd))))
    %+  sweep  inp
    |=  x=@  ^-  [c=@ t=@]
    ?+  arm  [x x]
      %exp     [(~(exp rh:math [%n .~~1e-2]) x) (~(exp rh:math-taylor [%n .~~1e-2]) x)]
      %log     [(~(log rh:math [%n .~~1e-2]) x) (~(log rh:math-taylor [%n .~~1e-2]) x)]
      %sin     [(~(sin rh:math [%n .~~1e-2]) x) (~(sin rh:math-taylor [%n .~~1e-2]) x)]
      %cos     [(~(cos rh:math [%n .~~1e-2]) x) (~(cos rh:math-taylor [%n .~~1e-2]) x)]
      %tan     [(~(tan rh:math [%n .~~1e-2]) x) (~(tan rh:math-taylor [%n .~~1e-2]) x)]
      %atan    [(~(atan rh:math [%n .~~1e-2]) x) (~(atan rh:math-taylor [%n .~~1e-2]) x)]
      %asin    [(~(asin rh:math [%n .~~1e-2]) x) (~(asin rh:math-taylor [%n .~~1e-2]) x)]
      %acos    [(~(acos rh:math [%n .~~1e-2]) x) (~(acos rh:math-taylor [%n .~~1e-2]) x)]
      %sqt     [(~(sqt rh:math [%n .~~1e-2]) x) (~(sqt rh:math-taylor [%n .~~1e-2]) x)]
      %cbt     [(~(cbt rh:math [%n .~~1e-2]) x) (~(cbt rh:math-taylor [%n .~~1e-2]) x)]
      %log-2   [(~(log-2 rh:math [%n .~~1e-2]) x) (~(log-2 rh:math-taylor [%n .~~1e-2]) x)]
      %log-10  [(~(log-10 rh:math [%n .~~1e-2]) x) (~(log-10 rh:math-taylor [%n .~~1e-2]) x)]
    ==
::
  %rq
    =/  den  (sun:rq den.dm)
    =/  lof  (div:rq (san:rq lo.dm) den)
    =/  wd   (sub:rq (div:rq (san:rq hi.dm) den) lof)
    =/  inp  |=(k=@ud ^-(@ (add:rq lof (mul:rq (div:rq (sun:rq k) (sun:rq pts)) wd))))
    %+  sweep  inp
    |=  x=@  ^-  [c=@ t=@]
    ?+  arm  [x x]
      %exp     [(~(exp rq:math [%n .~~~1e-10]) x) (~(exp rq:math-taylor [%n .~~~1e-10]) x)]
      %log     [(~(log rq:math [%n .~~~1e-10]) x) (~(log rq:math-taylor [%n .~~~1e-10]) x)]
      %sin     [(~(sin rq:math [%n .~~~1e-10]) x) (~(sin rq:math-taylor [%n .~~~1e-10]) x)]
      %cos     [(~(cos rq:math [%n .~~~1e-10]) x) (~(cos rq:math-taylor [%n .~~~1e-10]) x)]
      %tan     [(~(tan rq:math [%n .~~~1e-10]) x) (~(tan rq:math-taylor [%n .~~~1e-10]) x)]
      %atan    [(~(atan rq:math [%n .~~~1e-10]) x) (~(atan rq:math-taylor [%n .~~~1e-10]) x)]
      %asin    [(~(asin rq:math [%n .~~~1e-10]) x) (~(asin rq:math-taylor [%n .~~~1e-10]) x)]
      %acos    [(~(acos rq:math [%n .~~~1e-10]) x) (~(acos rq:math-taylor [%n .~~~1e-10]) x)]
      %sqt     [(~(sqt rq:math [%n .~~~1e-10]) x) (~(sqt rq:math-taylor [%n .~~~1e-10]) x)]
      %cbt     [(~(cbt rq:math [%n .~~~1e-10]) x) (~(cbt rq:math-taylor [%n .~~~1e-10]) x)]
      %log-2   [(~(log-2 rq:math [%n .~~~1e-10]) x) (~(log-2 rq:math-taylor [%n .~~~1e-10]) x)]
      %log-10  [(~(log-10 rq:math [%n .~~~1e-10]) x) (~(log-10 rq:math-taylor [%n .~~~1e-10]) x)]
    ==
==
