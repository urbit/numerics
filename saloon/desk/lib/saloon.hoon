  ::
::::  Saloon:  Scientific ALgorithms in hOON
::
::  Transcendental functions library for Urbit.
::
::  Pure Hoon implementations are generally naive formally correct algorithms,
::  awaiting efficient jetting.
::
/-  ls=lagoon
/+  *lagoon,
    math
::                                                    ::
::::                    ++sa                          ::  (2v) vector/matrix ops
~%  %saloon  ..part  ~
|%
::  +sake: set +sa params
::
::    rnd: rounding mode
::    rtol: relative tolerance, use the correct bit width @r
::
++  sake
  |=  [inrnd=rounding-mode inrtol=@r]
  %*(. sa rnd inrnd, rtol inrtol)
::
++  sa
  ^|
  =+  [rnd=*rounding-mode rtol=`@r`0x1]
  ~/  %sa-core
  |%
  ++  exp
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %exp))
  ++  log-2
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %log-2))
  ++  sin
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %sin))
  ++  sqrt
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %sqrt))
  ++  neg
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %neg))
  ++  trans-scalar
    |=  [=bloq =kind fun=@ta]
    ^-  $-(@ @)
    ?+    kind  ~|(kind !!)
        %real
      ?+    bloq  !!
          %7
        ?+  fun  !!
          %exp  ~(exp rq:math [rnd rtol])
          %log-2  ~(log-2 rq:math [rnd rtol])
          %neg  ~(neg rq:math [rnd rtol])
          %sin  ~(sin rq:math [rnd rtol])
          %sqrt  ~(sqrt rq:math [rnd rtol])
        ==
          %6
        ?+  fun  !!
          %exp  ~(exp rd:math [rnd rtol])
          %log-2  ~(log-2 rd:math [rnd rtol])
          %neg  ~(neg rd:math [rnd rtol])
          %sin  ~(sin rd:math [rnd rtol])
          %sqrt  ~(sqrt rd:math [rnd rtol])
        ==
          %5
        ?+  fun  !!
          %exp  ~(exp rs:math [rnd rtol])
          %log-2  ~(log-2 rs:math [rnd rtol])
          %neg  ~(neg rs:math [rnd rtol])
          %sin  ~(sin rs:math [rnd rtol])
          %sqrt  ~(sqrt rs:math [rnd rtol])
        ==
          %4
        ?+  fun  !!
          %exp  ~(exp rh:math [rnd rtol])
          %log-2  ~(log-2 rh:math [rnd rtol])
          %neg  ~(neg rh:math [rnd rtol])
          %sin  ~(sin rh:math [rnd rtol])
          %sqrt  ~(sqrt rh:math [rnd rtol])
        ==
      ==
    ==
  --
--
