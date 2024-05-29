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
  ::  innermost core of Saloon functionality
  ::  basically an implementation of /lib/math for Lagoon arrays
  +|  %uno
  ::
  ::  Comparison
  ::
  ::    +lth:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, less than
  ::    Examples
  ::      > (lth:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.0000.0000.0000.0000.0000.0000.3f80.0000.0000.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (lth:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.0] t=[i=~[.0] t=~[~[.0] ~[.1] ~[.0]]]]
  ::  Source
  ++  lth  lth:(lake rnd)
  ::    +lte:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, less than or equal to
  ::    Examples
  ::      > (lte:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.3f80.0000.0000.0000.3f80.0000.3f80.0000.0000.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (lte:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.0] t=~[~[.1] ~[.1] ~[.0]]]]
  ::  Source
  ++  lte  lte:(lake rnd)
  ::    +leq:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, less than or equal to
  ::  Alias for +lte.
  ::    Examples
  ::      > (leq:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.3f80.0000.0000.0000.3f80.0000.3f80.0000.0000.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (leq:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.0] t=~[~[.1] ~[.1] ~[.0]]]]
  ::  Source
  ++  leq  lte:(lake rnd)
  ::    +equ:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, equal to
  ::    Examples
  ::      > (equ:la:la (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.3f80.0000.0000.0000.3f80.0000.0000.0000.0000.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (equ:la:la (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.0] t=~[~[.1] ~[.0] ~[.0]]]]
  ::  Source
  ++  equ  equ:(lake rnd)
  ::    +gth:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, greater than
  ::    Examples
  ::      > (gth:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=data=0x1.0000.0000.3f80.0000.0000.0000.0000.0000.3f80.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (gth:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.0] t=[i=~[.1] t=~[~[.0] ~[.0] ~[.1]]]]
  ::  Source
  ++  gth  gth:(lake rnd)
  ::    +gte:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, greater than or equal to
  ::    Examples
  ::      > (gte:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=data=0x1.3f80.0000.3f80.0000.3f80.0000.0000.0000.3f80.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (gte:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.1] t=~[~[.1] ~[.0] ~[.1]]]]
  ::  Source
  ++  gte  gte:(lake rnd)
  ::    +geq:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, greater than or equal to
  ::  Alias for +gte.
  ::    Examples
  ::      > (geq:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=data=0x1.3f80.0000.3f80.0000.3f80.0000.0000.0000.3f80.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (geq:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.1] t=~[~[.1] ~[.0] ~[.1]]]]
  ::  Source
  ++  geq  gte:(lake rnd)
  ::    +neq:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, equal to
  ::    Examples
  ::      > (neq:la:la (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.0000.0000.3f80.0000.0000.0000.3f80.0000.3f80.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (neq:la:la (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.0] t=[i=~[.1] t=~[~[.0] ~[.1] ~[.1]]]]
  ::  Source
  ++  neq  neq:(lake rnd)
  ::    +is-close:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, close to
  ::    Examples
  ::      > (is-close:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.3f80.0000.0000.0000.3f80.0000.3f80.0000.0000.0000.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (is-close:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.0] t=~[~[.1] ~[.1] ~[.0]]]
  ::  Source
  ++  is-close
    |=  [a=ray:ls b=ray:ls]
    ^-  ray:ls
    (is-close:(lake rnd) a b [0x0 rtol])
  ::    +all-close:  [$ray $ray] -> ?
  ::
  ::  Returns the LOOBEAN comparison of two floating-point rays, all close to
  ::    Examples
  ::      > (all-close:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      %.n
  ::      > (all-close:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (ones:la:la [~[5 1] 5 %i754 ~]))
  ::      %.y
  ::  Source
  ++  all-close  all:(lake rnd)
  ::    +any-close:  [$ray $ray] -> ?
  ::
  ::  Returns the LOOBEAN comparison of two floating-point rays, any close to
  ::    Examples
  ::      > (any-close:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      %.y
  ::      > (any-close:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      %.y
  ::  Source
  ++  any-close  any:(lake rnd)
  ::
  ::  Algebraic
  ::
  ::    +add:  [$ray $ray] -> $ray
  ::
  ::  Returns the sum of two floating-point rays
  ::  Source
  ++  add  add:(lake rnd)
  ::    +sub:  [$ray $ray] -> $ray
  ::
  ::  Returns the difference of two floating-point rays
  ::  Source
  ++  sub  sub:(lake rnd)
  ::    +mul:  [$ray $ray] -> $ray
  ::
  ::  Returns the product of two floating-point rays
  ::  Source
  ++  mul  mul:(lake rnd)
  ::    +div:  [$ray $ray] -> $ray
  ::
  ::  Returns the quotient of two floating-point rays
  ::  Source
  ++  div  div:(lake rnd)
  ::    +fma:  [$ray $ray $ray] -> $ray
  ::
  ::  Returns the fused multiply-add of three floating-point rays
  ::  Examples
  ::    > (fma:sa:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (ones:la:la [~[5 1] 5 %i754 ~]) (ones:la:la [~[5 1] 5 %i754 ~]))
  ::    [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.4000.0000.4000.0000.4000.0000.4000.0000.4000.0000]
  ::  Source
  ++  fma  |=([a=ray:ls b=ray:ls c=ray:ls] (add:(lake rnd) (mul:(lake rnd) a b) c))
  ::    +neg:  $ray -> $ray
  ::
  ::  Returns the negation of each entry in a floating-point ray
  ::  Source
  ++  neg
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %neg))
  ::    +factorial:  $ray -> $ray
  ::
  ::  Returns the factorial of each entry in a floating-point ray
  ::  Source
  ++  factorial
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %factorial))
  ::    +abs: $ray -> $ray
  ::
  ::  Returns the absolute value of each entry in a floating-point ray
  ::  Source
  ++  abs  abs:(lake rnd)
  ::    +exp: $ray -> $ray
  ::
  ::  Returns the exponential of each entry in a floating-point ray
  ::  Source
  ++  exp
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %exp))
  ::    +sin: $ray -> $ray
  ::
  ::  Returns the sine of each entry in a floating-point ray
  ::  Source
  ++  sin
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %sin))
  ::    +cos: $ray -> $ray
  ::
  ::  Returns the cosine of each entry in a floating-point ray
  ::  Source
  ++  cos
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %cos))
  ::    +tan: $ray -> $ray
  ::
  ::  Returns the tangent of each entry in a floating-point ray
  ::  Source
  ++  tan
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %tan))
  ::    +pow-n: [$ray $ray] -> $ray
  ::
  ::  Returns the exponentiation of each entry in a floating-point ray by another ray
  ::  Source
  ++  pow-n
    |=  [a=ray:ls b=ray:ls]
    ^-  ray
    (bin-op:la a b (fun-scalar meta.a %pow-n))
  ::    +log: $ray -> $ray
  ::
  ::  Returns the natural logarithm of each entry in a floating-point ray
  ::  Source
  ++  log
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %log))
  ::    +log-10: $ray -> $ray
  ::
  ::  Returns the base-10 logarithm of each entry in a floating-point ray
  ::  Source
  ++  log-10
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %log-10))
  ::    +log-2: $ray -> $ray
  ::
  ::  Returns the base-2 logarithm of each entry in a floating-point ray
  ::  Source
  ++  log-2
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %log-2))
  ::    +pow: [$ray $ray] -> $ray
  ::
  ::  Returns the exponentiation of each entry in a floating-point ray by another ray
  ::  Source
  ++  pow
    |=  [a=ray:ls b=ray:ls]
    ^-  ray
    (bin-op:la a b (fun-scalar meta.a %pow))
  ::    +sqrt: $ray -> $ray
  ::
  ::  Returns the square root of each entry in a floating-point ray
  ::  Source
  ++  sqrt
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %sqrt))
  ::    +sqt: $ray -> $ray
  ::
  ::  Returns the square root of each entry in a floating-point ray
  ::  Alias for +sqrt.
  ::  Source
  ++  sqt  sqrt
  ::    +cbrt: $ray -> $ray
  ::
  ::  Returns the cube root of each entry in a floating-point ray
  ::  Source
  ++  cbrt
    |=  a=ray:ls
    ^-  ray
    (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %cbrt))
  ::    +cbt: $ray -> $ray
  ::
  ::  Returns the cube root of each entry in a floating-point ray
  ::  Alias for +cbrt.
  ::  Source
  ++  cbt  cbrt
  ::
  +$  unary-ops   $?  %neg
                      %factorial
                      %exp
                      %sin
                      %cos
                      %tan
                      %log
                      %log-10
                      %log-2
                      %sqrt
                      %cbrt
                  ==
  ::
  ++  trans-scalar
    |=  [=bloq =kind fun=unary-ops]
    ^-  $-(@ @)
    ?-    kind
        %int2  !!
        %uint
      ?-  fun
        %neg        !!
        %factorial  |=(x=@u ^-(@u =/(t 1 ?:(=(0 x) t ?:(=(1 x) t |-(?:(=(1 x) t $(x (^sub x 1), t (^mul t x)))))))))
        %exp        !!
        %sin        !!
        %cos        !!
        %tan        !!
        %log        !!
        %log-10     !!
        %log-2      !!
        %sqrt       !!
        %cbrt       !!
      ==  ::  fun
      ::
        %i754
      ?+    bloq  !!
          %7
        ?-  fun
          %neg        ~(neg rq:math [rnd rtol])
          %factorial  ~(factorial rq:math [rnd rtol])
          %exp        ~(exp rq:math [rnd rtol])
          %sin        ~(sin rq:math [rnd rtol])
          %cos        ~(cos rq:math [rnd rtol])
          %tan        ~(tan rq:math [rnd rtol])
          %log        ~(log rq:math [rnd rtol])
          %log-10     ~(log-10 rq:math [rnd rtol])
          %log-2      ~(log-2 rq:math [rnd rtol])
          %sqrt       ~(sqrt rq:math [rnd rtol])
          %cbrt       ~(cbrt rq:math [rnd rtol])
        ==  ::  fun
          %6
        ?-  fun
          %neg        ~(neg rd:math [rnd rtol])
          %factorial  ~(factorial rd:math [rnd rtol])
          %exp        ~(exp rd:math [rnd rtol])
          %sin        ~(sin rd:math [rnd rtol])
          %cos        ~(cos rd:math [rnd rtol])
          %tan        ~(tan rd:math [rnd rtol])
          %log        ~(log rd:math [rnd rtol])
          %log-10     ~(log-10 rd:math [rnd rtol])
          %log-2      ~(log-2 rd:math [rnd rtol])
          %sqrt       ~(sqrt rd:math [rnd rtol])
          %cbrt       ~(cbrt rd:math [rnd rtol])
        ==  ::  fun
          %5
        ?-  fun
          %neg        ~(neg rs:math [rnd rtol])
          %factorial  ~(factorial rs:math [rnd rtol])
          %exp        ~(exp rs:math [rnd rtol])
          %sin        ~(sin rs:math [rnd rtol])
          %cos        ~(cos rs:math [rnd rtol])
          %tan        ~(tan rs:math [rnd rtol])
          %log        ~(log rs:math [rnd rtol])
          %log-10     ~(log-10 rs:math [rnd rtol])
          %log-2      ~(log-2 rs:math [rnd rtol])
          %sqrt       ~(sqrt rs:math [rnd rtol])
          %cbrt       ~(cbrt rs:math [rnd rtol])
        ==  ::  fun
          %4
        ?-  fun
          %neg        ~(neg rh:math [rnd rtol])
          %factorial  ~(factorial rh:math [rnd rtol])
          %exp        ~(exp rh:math [rnd rtol])
          %sin        ~(sin rh:math [rnd rtol])
          %cos        ~(cos rh:math [rnd rtol])
          %tan        ~(tan rh:math [rnd rtol])
          %log        ~(log rh:math [rnd rtol])
          %log-10     ~(log-10 rh:math [rnd rtol])
          %log-2      ~(log-2 rh:math [rnd rtol])
          %sqrt       ~(sqrt rh:math [rnd rtol])
          %cbrt       ~(cbrt rh:math [rnd rtol])
        ==  ::  fun
      ==  ::  bloq
    ==  ::  kind
  ::
  +$  binary-ops  $?  %pow-n
                      %pow
                  ==
  ::
  ++  fun-scalar
    |=  [=meta fun=binary-ops]
    ^-  $-([@ @] @)
    ?-    kind.meta
        %int2  !!
        %uint
      ?-  fun
        %pow        (fun-scalar meta %pow-n)
        %pow-n      |=([x=@u n=@u] ^-(@u ?:(=(0 n) 1 =/(p x |-(?:((^lth n 2) p $(n (dec n), p (^mul p x))))))))
      ==  ::  fun
      ::
        %i754
      ?+    bloq.meta  !!
          %7
        ?-  fun
          %pow-n      ~(pow-n rq:math [rnd rtol])
          %pow        ~(pow rq:math [rnd rtol])
        ==  ::  fun
          %6
        ?-  fun
          %pow-n      ~(pow-n rd:math [rnd rtol])
          %pow        ~(pow rd:math [rnd rtol])
        ==  ::  fun
          %5
        ?-  fun
          %pow-n      ~(pow-n rs:math [rnd rtol])
          %pow        ~(pow rs:math [rnd rtol])
        ==  ::  fun
          %4
        ?-  fun
          %pow-n      ~(pow-n rh:math [rnd rtol])
          %pow        ~(pow rh:math [rnd rtol])
        ==  ::  fun
      ==  ::  bloq
    ==  ::  kind
  --
--
