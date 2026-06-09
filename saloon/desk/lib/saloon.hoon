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
    math,
    complex
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
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.0000.0000.3f80.0000.0000.0000.0000.0000.3f80.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (gth:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.0] t=[i=~[.1] t=~[~[.0] ~[.0] ~[.1]]]]
  ::  Source
  ++  gth  gth:(lake rnd)
  ::    +gte:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, greater than or equal to
  ::    Examples
  ::      > (gte:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.3f80.0000.3f80.0000.3f80.0000.0000.0000.3f80.0000]
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
  ::      [meta=[shape=~[5 1] bloq=5 kind=%i754 fxp=~] data=0x1.3f80.0000.3f80.0000.3f80.0000.0000.0000.3f80.0000]
  ::      > ;;((list (list @rs)) data:(de-ray:la:la (geq:sa:sa (ones:la:la [~[5 1] 5 %i754 ~]) (en-ray:la:la [~[5 1] 5 %i754 ~] ~[.1 .0 .1 .2 .0]))))
  ::      [i=~[.1] t=[i=~[.1] t=~[~[.1] ~[.0] ~[.1]]]]
  ::  Source
  ++  geq  gte:(lake rnd)
  ::    +neq:  [$ray $ray] -> $ray
  ::
  ::  Returns the BOOLEAN comparison of two floating-point rays, not equal to
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
    ?+    kind  !!
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
      ::
      ::  posits (/lib/unum): rpb/rph/rps/rpd are bloq 3/4/5/6 (NOT %i754's
      ::  4-7), and unum's square root is named +sqt.  Posit rounding is fixed
      ::  (round-to-nearest-even), so rnd/rtol don't apply.
        %unum
      ?+    bloq  !!
        %6  ?-(fun %neg neg:rpd:unum, %factorial factorial:rpd:unum, %exp exp:rpd:unum, %sin sin:rpd:unum, %cos cos:rpd:unum, %tan tan:rpd:unum, %log log:rpd:unum, %log-10 log-10:rpd:unum, %log-2 log-2:rpd:unum, %sqrt sqt:rpd:unum, %cbrt cbrt:rpd:unum)
        %5  ?-(fun %neg neg:rps:unum, %factorial factorial:rps:unum, %exp exp:rps:unum, %sin sin:rps:unum, %cos cos:rps:unum, %tan tan:rps:unum, %log log:rps:unum, %log-10 log-10:rps:unum, %log-2 log-2:rps:unum, %sqrt sqt:rps:unum, %cbrt cbrt:rps:unum)
        %4  ?-(fun %neg neg:rph:unum, %factorial factorial:rph:unum, %exp exp:rph:unum, %sin sin:rph:unum, %cos cos:rph:unum, %tan tan:rph:unum, %log log:rph:unum, %log-10 log-10:rph:unum, %log-2 log-2:rph:unum, %sqrt sqt:rph:unum, %cbrt cbrt:rph:unum)
        %3  ?-(fun %neg neg:rpb:unum, %factorial factorial:rpb:unum, %exp exp:rpb:unum, %sin sin:rpb:unum, %cos cos:rpb:unum, %tan tan:rpb:unum, %log log:rpb:unum, %log-10 log-10:rpb:unum, %log-2 log-2:rpb:unum, %sqrt sqt:rpb:unum, %cbrt cbrt:rpb:unum)
      ==
      ::
      ::  complex (/lib/complex): bloq 5/6/7/8 = ch/cs/cd/cq.  Real-only / ordered
      ::  functions (factorial, log-10, log-2, cbrt) are undefined on complex.
        %cplx
      ?+    bloq  !!
        %8  ?-(fun %neg ~(neg cq:complex rnd), %exp ~(cexp cq:complex rnd), %sin ~(csin cq:complex rnd), %cos ~(ccos cq:complex rnd), %tan ~(ctan cq:complex rnd), %log ~(clog cq:complex rnd), %sqrt ~(csqrt cq:complex rnd), %factorial |=(@ !!), %log-10 |=(@ !!), %log-2 |=(@ !!), %cbrt |=(@ !!))
        %7  ?-(fun %neg ~(neg cd:complex rnd), %exp ~(cexp cd:complex rnd), %sin ~(csin cd:complex rnd), %cos ~(ccos cd:complex rnd), %tan ~(ctan cd:complex rnd), %log ~(clog cd:complex rnd), %sqrt ~(csqrt cd:complex rnd), %factorial |=(@ !!), %log-10 |=(@ !!), %log-2 |=(@ !!), %cbrt |=(@ !!))
        %6  ?-(fun %neg ~(neg cs:complex rnd), %exp ~(cexp cs:complex rnd), %sin ~(csin cs:complex rnd), %cos ~(ccos cs:complex rnd), %tan ~(ctan cs:complex rnd), %log ~(clog cs:complex rnd), %sqrt ~(csqrt cs:complex rnd), %factorial |=(@ !!), %log-10 |=(@ !!), %log-2 |=(@ !!), %cbrt |=(@ !!))
        %5  ?-(fun %neg ~(neg ch:complex rnd), %exp ~(cexp ch:complex rnd), %sin ~(csin ch:complex rnd), %cos ~(ccos ch:complex rnd), %tan ~(ctan ch:complex rnd), %log ~(clog ch:complex rnd), %sqrt ~(csqrt ch:complex rnd), %factorial |=(@ !!), %log-10 |=(@ !!), %log-2 |=(@ !!), %cbrt |=(@ !!))
      ==
    ==  ::  kind
  ::
  +$  binary-ops  $?  %pow-n
                      %pow
                  ==
  ::
  ++  fun-scalar
    |=  [=meta fun=binary-ops]
    ^-  $-([@ @] @)
    ?+    kind.meta  !!
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
      ::  posits (/lib/unum): bloq 3/4/5/6.
        %unum
      ?+    bloq.meta  !!
        %6  ?-(fun %pow-n pow-n:rpd:unum, %pow pow:rpd:unum)
        %5  ?-(fun %pow-n pow-n:rps:unum, %pow pow:rps:unum)
        %4  ?-(fun %pow-n pow-n:rph:unum, %pow pow:rph:unum)
        %3  ?-(fun %pow-n pow-n:rpb:unum, %pow pow:rpb:unum)
      ==
      ::  complex: z^w via +cpow (both rays complex); integer +pow-n undefined.
        %cplx
      ?+    bloq.meta  !!
        %8  ?-(fun %pow ~(cpow cq:complex rnd), %pow-n |=([@ @] !!))
        %7  ?-(fun %pow ~(cpow cd:complex rnd), %pow-n |=([@ @] !!))
        %6  ?-(fun %pow ~(cpow cs:complex rnd), %pow-n |=([@ @] !!))
        %5  ?-(fun %pow ~(cpow ch:complex rnd), %pow-n |=([@ @] !!))
      ==
    ==  ::  kind
  ::
  +|  %linalg
  ::
  ::  Eigendecomposition of symmetric real matrices via the cyclic Jacobi
  ::  algorithm.  Phase A: %i754 single (bloq 5) and double (bloq 6) only.
  ::  Real eigenvalues, orthonormal eigenvectors, no complex arithmetic.
  ::
  ::    Scalar float helpers, dispatched on the component bloq (4=@rh, 5=@rs,
  ::    6=@rd, 7=@rq).  Each takes the bloq as its first argument and operates on
  ::    raw component atoms.
  ::
  ++  fadd  |=([b=@ x=@ y=@] ^-(@ ?:(=(4 b) (~(add rh:math [rnd rtol]) x y) ?:(=(5 b) (~(add rs:math [rnd rtol]) x y) ?:(=(6 b) (~(add rd:math [rnd rtol]) x y) (~(add rq:math [rnd rtol]) x y))))))
  ++  fsub  |=([b=@ x=@ y=@] ^-(@ ?:(=(4 b) (~(sub rh:math [rnd rtol]) x y) ?:(=(5 b) (~(sub rs:math [rnd rtol]) x y) ?:(=(6 b) (~(sub rd:math [rnd rtol]) x y) (~(sub rq:math [rnd rtol]) x y))))))
  ++  fmul  |=([b=@ x=@ y=@] ^-(@ ?:(=(4 b) (~(mul rh:math [rnd rtol]) x y) ?:(=(5 b) (~(mul rs:math [rnd rtol]) x y) ?:(=(6 b) (~(mul rd:math [rnd rtol]) x y) (~(mul rq:math [rnd rtol]) x y))))))
  ++  fdiv  |=([b=@ x=@ y=@] ^-(@ ?:(=(4 b) (~(div rh:math [rnd rtol]) x y) ?:(=(5 b) (~(div rs:math [rnd rtol]) x y) ?:(=(6 b) (~(div rd:math [rnd rtol]) x y) (~(div rq:math [rnd rtol]) x y))))))
  ++  fabs  |=([b=@ x=@] ^-(@ ?:(=(4 b) (~(abs rh:math [rnd rtol]) x) ?:(=(5 b) (~(abs rs:math [rnd rtol]) x) ?:(=(6 b) (~(abs rd:math [rnd rtol]) x) (~(abs rq:math [rnd rtol]) x))))))
  ++  fgte  |=([b=@ x=@ y=@] ^-(? ?:(=(4 b) (~(gte rh:math [rnd rtol]) x y) ?:(=(5 b) (~(gte rs:math [rnd rtol]) x y) ?:(=(6 b) (~(gte rd:math [rnd rtol]) x y) (~(gte rq:math [rnd rtol]) x y))))))
  ++  flte  |=([b=@ x=@ y=@] ^-(? ?:(=(4 b) (~(lte rh:math [rnd rtol]) x y) ?:(=(5 b) (~(lte rs:math [rnd rtol]) x y) ?:(=(6 b) (~(lte rd:math [rnd rtol]) x y) (~(lte rq:math [rnd rtol]) x y))))))
  ++  f0    |=(b=@ ^-(@ ?:(=(4 b) .~~0 ?:(=(5 b) .0 ?:(=(6 b) .~0 .~~~0)))))
  ++  f1    |=(b=@ ^-(@ ?:(=(4 b) .~~1 ?:(=(5 b) .1 ?:(=(6 b) .~1 .~~~1)))))
  ++  f2    |=(b=@ ^-(@ ?:(=(4 b) .~~2 ?:(=(5 b) .2 ?:(=(6 b) .~2 .~~~2)))))
  ::    width-fixed relative epsilon for sqrt convergence (decoupled from rtol).
  ++  feps  |=(b=@ ^-(@ ?:(=(4 b) .~~1e-2 ?:(=(5 b) .1e-6 ?:(=(6 b) .~1e-13 .~~~1e-30)))))
  ::    +fsqt: iteration-capped Newton sqrt.  The 50-step cap is a durability
  ::    backstop: a sub-ULP tolerance otherwise makes Newton oscillate forever
  ::    (which once OOM-crashed the ship).  Convergence uses feps, not rtol.
  ++  fsqt
    |=  [b=@ x=@]
    ^-  @
    ?:  (flte b x (f0 b))  (f0 b)
    =/  g  x
    =/  i  0
    |-  ^-  @
    ?:  =(50 i)  g
    =/  ng  (fmul b (fdiv b (f1 b) (f2 b)) (fadd b g (fdiv b x g)))
    ?:  (flte b (fabs b (fsub b g ng)) (fmul b (feps b) g))
      ng
    $(i +(i), g ng)
  ++  fneg  |=([b=@ x=@] ^-(@ (fsub b (f0 b) x)))
  ++  fsign  |=([b=@ x=@] ^-(@ ?:((fgte b x (f0 b)) (f1 b) (fneg b (f1 b)))))
  ::    Complex helpers, dispatched on the COMPLEX bloq (5=@ch, 6=@cs, 7=@cd,
  ::    8=@cq), over /lib/complex.  The real component bloq is always (dec cb).
  ++  cb-comp  |=(cb=@ ^-(@ (dec cb)))
  ++  cadd  |=([cb=@ p=@ q=@] ^-(@ ?:(=(5 cb) (~(add ch:complex rnd) p q) ?:(=(6 cb) (~(add cs:complex rnd) p q) ?:(=(7 cb) (~(add cd:complex rnd) p q) (~(add cq:complex rnd) p q))))))
  ++  csub  |=([cb=@ p=@ q=@] ^-(@ ?:(=(5 cb) (~(sub ch:complex rnd) p q) ?:(=(6 cb) (~(sub cs:complex rnd) p q) ?:(=(7 cb) (~(sub cd:complex rnd) p q) (~(sub cq:complex rnd) p q))))))
  ++  cmul  |=([cb=@ p=@ q=@] ^-(@ ?:(=(5 cb) (~(mul ch:complex rnd) p q) ?:(=(6 cb) (~(mul cs:complex rnd) p q) ?:(=(7 cb) (~(mul cd:complex rnd) p q) (~(mul cq:complex rnd) p q))))))
  ++  cdiv  |=([cb=@ p=@ q=@] ^-(@ ?:(=(5 cb) (~(div ch:complex rnd) p q) ?:(=(6 cb) (~(div cs:complex rnd) p q) ?:(=(7 cb) (~(div cd:complex rnd) p q) (~(div cq:complex rnd) p q))))))
  ++  cconj  |=([cb=@ p=@] ^-(@ ?:(=(5 cb) (~(conj ch:complex rnd) p) ?:(=(6 cb) (~(conj cs:complex rnd) p) ?:(=(7 cb) (~(conj cd:complex rnd) p) (~(conj cq:complex rnd) p))))))
  ++  cre  |=([cb=@ p=@] ^-(@ ?:(=(5 cb) (~(re ch:complex rnd) p) ?:(=(6 cb) (~(re cs:complex rnd) p) ?:(=(7 cb) (~(re cd:complex rnd) p) (~(re cq:complex rnd) p))))))
  ++  cpak  |=([cb=@ r=@ i=@] ^-(@ ?:(=(5 cb) (~(pak ch:complex rnd) r i) ?:(=(6 cb) (~(pak cs:complex rnd) r i) ?:(=(7 cb) (~(pak cd:complex rnd) r i) (~(pak cq:complex rnd) r i))))))
  ++  cabs-re  |=([cb=@ p=@] ^-(@ (cre cb ?:(=(5 cb) (~(abs ch:complex rnd) p) ?:(=(6 cb) (~(abs cs:complex rnd) p) ?:(=(7 cb) (~(abs cd:complex rnd) p) (~(abs cq:complex rnd) p)))))))
  ::    Indexed scalar access over the rounding-bound Lagoon door.
  ++  gi  |=([m=ray:ls ix=(list @)] ^-(@ (get-item:(lake rnd) m ix)))
  ++  si  |=([m=ray:ls ix=(list @) val=@] ^-(ray:ls (set-item:(lake rnd) m ix val)))
  ::    +stol/+near/+cnear: (skew-)symmetry tolerance.  An eig input is accepted
  ::    if it is symmetric/Hermitian WITHIN a relative tolerance, so a matrix
  ::    that is symmetric only up to rounding (e.g. a Gram matrix M^T*M from
  ::    Lagoon mmul, whose [i,j] and [j,i] accumulate in different orders) is not
  ::    rejected, while a genuinely asymmetric matrix (entries differing by O(1))
  ::    still crashes.  The check is magnitude-based, so a +-0.0 sign difference
  ::    in conjugate pairs (which arises under %d/%u rounding) also passes.
  ++  stol  |=(b=@ ^-(@ ?:(=(4 b) .~~1e-2 ?:(=(5 b) .1e-4 ?:(=(6 b) .~1e-9 .~~~1e-18)))))
  ++  near  |=([b=@ x=@ y=@] ^-(? (flte b (fabs b (fsub b x y)) (fadd b (stol b) (fmul b (stol b) (fadd b (fabs b x) (fabs b y)))))))
  ++  cnear
    |=  [cb=@ x=@ y=@]
    ^-  ?
    =/  rb  (cb-comp cb)
    (flte rb (cabs-re cb (csub cb x (cconj cb y))) (fadd rb (stol rb) (fmul rb (stol rb) (fadd rb (cabs-re cb x) (cabs-re cb y)))))
  ::    +symmetric: m equals its transpose within the symmetry tolerance (+near).
  ++  symmetric
    |=  m=ray:ls
    ^-  ?
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  k  0
    |-  ^-  ?
    ?:  =(k (^mul n n))  &
    =/  i  (^div k n)
    =/  j  (mod k n)
    ?:  (^lte i j)  $(k +(k))
    ?.  (near b (gi m ~[i j]) (gi m ~[j i]))  |
    $(k +(k))
  ::    +off-norm: Frobenius norm of the strictly off-diagonal part.
  ++  off-norm
    |=  m=ray:ls
    ^-  @
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  acc  (f0 b)
    =/  k  0
    |-  ^-  @
    ?:  =(k (^mul n n))  (fsqt b acc)
    =/  i  (^div k n)
    =/  j  (mod k n)
    ?:  =(i j)  $(k +(k))
    =/  mij  (gi m ~[i j])
    $(k +(k), acc (fadd b acc (fmul b mij mij)))
  ::    +frob: full Frobenius norm of m.
  ++  frob
    |=  m=ray:ls
    ^-  @
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  acc  (f0 b)
    =/  k  0
    |-  ^-  @
    ?:  =(k (^mul n n))  (fsqt b acc)
    =/  i  (^div k n)
    =/  j  (mod k n)
    =/  mij  (gi m ~[i j])
    $(k +(k), acc (fadd b acc (fmul b mij mij)))
  ::    +diag: extract the diagonal of a square matrix as a 1-D ray.
  ++  diag
    |=  m=ray:ls
    ^-  ray:ls
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  r  (zeros:(lake rnd) [~[n] b %i754 ~])
    =/  i  0
    |-  ^-  ray:ls
    ?:  =(i n)  r
    =.  r  (si r ~[i] (gi m ~[i i]))
    $(i +(i))
  ::    +rot-cols: post-multiply columns p,q of m by the Givens rotation J
  ::    (J_pp=J_qq=c, J_pq=s, J_qp=-s), i.e. m <- m*J.
  ++  rot-cols
    |=  [m=ray:ls p=@ q=@ c=@ s=@]
    ^-  ray:ls
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  i  0
    |-  ^-  ray:ls
    ?:  =(i n)  m
    =/  mip  (gi m ~[i p])
    =/  miq  (gi m ~[i q])
    =.  m  (si m ~[i p] (fsub b (fmul b c mip) (fmul b s miq)))
    =.  m  (si m ~[i q] (fadd b (fmul b s mip) (fmul b c miq)))
    $(i +(i))
  ::    +rot-rows: pre-multiply rows p,q of m by J-transpose, i.e. m <- J^T*m.
  ++  rot-rows
    |=  [m=ray:ls p=@ q=@ c=@ s=@]
    ^-  ray:ls
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  j  0
    |-  ^-  ray:ls
    ?:  =(j n)  m
    =/  mpj  (gi m ~[p j])
    =/  mqj  (gi m ~[q j])
    =.  m  (si m ~[p j] (fsub b (fmul b c mpj) (fmul b s mqj)))
    =.  m  (si m ~[q j] (fadd b (fmul b s mpj) (fmul b c mqj)))
    $(j +(j))
  ::    +sweep-once: one cyclic Jacobi sweep over every p<q pair, applying a
  ::    rotation that zeros a_pq to both the working matrix and the
  ::    accumulated eigenvector matrix.
  ++  sweep-once
    |=  [m=ray:ls v=ray:ls]
    ^-  [ray:ls ray:ls]
    =/  b  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    ?:  (^lte n 1)  [m v]
    =/  p  0
    =/  q  1
    |-  ^-  [ray:ls ray:ls]
    =/  apq  (gi m ~[p q])
    =/  mv=[ray:ls ray:ls]
      ?:  =(apq (f0 b))
        [m v]
      =/  app  (gi m ~[p p])
      =/  aqq  (gi m ~[q q])
      =/  theta  (fdiv b (fsub b aqq app) (fmul b (f2 b) apq))
      =/  t  (fdiv b (fsign b theta) (fadd b (fabs b theta) (fsqt b (fadd b (fmul b theta theta) (f1 b)))))
      =/  c  (fdiv b (f1 b) (fsqt b (fadd b (fmul b t t) (f1 b))))
      =/  s  (fmul b t c)
      :-  (rot-rows (rot-cols m p q c s) p q c s)
      (rot-cols v p q c s)
    =.  m  -.mv
    =.  v  +.mv
    ?:  =(+(q) n)
      ?:  =(+(p) (dec n))  [m v]
      $(p +(p), q (^add p 2))
    $(q +(q))
  ::  Hermitian (%cplx) Jacobi.  The complex rotation J that zeros a_pq has a
  ::  real diagonal c and complex off-diagonal b = s*(a_pq/|a_pq|), with
  ::  J[q,p] = -conj(b); updates are A <- J^H*A*J and V <- V*J.  Eigenvalues are
  ::  real (the diagonal at convergence), eigenvectors form a unitary matrix.
  ::
  ::    +hermitian: a_ij == conj(a_ji) within the symmetry tolerance (+cnear),
  ::    which implies a (near-)real diagonal.  Magnitude-based, so a real
  ::    diagonal stored with -0.0 imaginary parts is accepted.
  ++  hermitian
    |=  m=ray:ls
    ^-  ?
    =/  cb  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  k  0
    |-  ^-  ?
    ?:  =(k (^mul n n))  &
    =/  i  (^div k n)
    =/  j  (mod k n)
    ?:  (^lth i j)  $(k +(k))
    ?.  (cnear cb (gi m ~[i j]) (gi m ~[j i]))  |
    $(k +(k))
  ::    +off-norm-h: real Frobenius norm of the strictly off-diagonal part.
  ++  off-norm-h
    |=  m=ray:ls
    ^-  @
    =/  cb  bloq.meta.m
    =/  rb  (cb-comp cb)
    =/  n  (snag 0 shape.meta.m)
    =/  acc  (f0 rb)
    =/  k  0
    |-  ^-  @
    ?:  =(k (^mul n n))  (fsqt rb acc)
    =/  i  (^div k n)
    =/  j  (mod k n)
    ?:  =(i j)  $(k +(k))
    =/  mag  (cabs-re cb (gi m ~[i j]))
    $(k +(k), acc (fadd rb acc (fmul rb mag mag)))
  ::    +frob-h: full real Frobenius norm of a complex matrix.
  ++  frob-h
    |=  m=ray:ls
    ^-  @
    =/  cb  bloq.meta.m
    =/  rb  (cb-comp cb)
    =/  n  (snag 0 shape.meta.m)
    =/  acc  (f0 rb)
    =/  k  0
    |-  ^-  @
    ?:  =(k (^mul n n))  (fsqt rb acc)
    =/  i  (^div k n)
    =/  j  (mod k n)
    =/  mag  (cabs-re cb (gi m ~[i j]))
    $(k +(k), acc (fadd rb acc (fmul rb mag mag)))
  ::    +diag-real: real parts of the diagonal as a 1-D %i754 ray.
  ++  diag-real
    |=  m=ray:ls
    ^-  ray:ls
    =/  cb  bloq.meta.m
    =/  rb  (cb-comp cb)
    =/  n  (snag 0 shape.meta.m)
    =/  r  (zeros:(lake rnd) [~[n] rb %i754 ~])
    =/  i  0
    |-  ^-  ray:ls
    ?:  =(i n)  r
    =.  r  (si r ~[i] (cre cb (gi m ~[i i])))
    $(i +(i))
  ::    +rot-cols-h: m <- m*J  (complex Givens; J[q,p] = -conj(b)).
  ++  rot-cols-h
    |=  [m=ray:ls p=@ q=@ cc=@ b=@]
    ^-  ray:ls
    =/  cb  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  i  0
    |-  ^-  ray:ls
    ?:  =(i n)  m
    =/  mip  (gi m ~[i p])
    =/  miq  (gi m ~[i q])
    =.  m  (si m ~[i p] (csub cb (cmul cb cc mip) (cmul cb (cconj cb b) miq)))
    =.  m  (si m ~[i q] (cadd cb (cmul cb b mip) (cmul cb cc miq)))
    $(i +(i))
  ::    +rot-rows-h: m <- J^H*m.
  ++  rot-rows-h
    |=  [m=ray:ls p=@ q=@ cc=@ b=@]
    ^-  ray:ls
    =/  cb  bloq.meta.m
    =/  n  (snag 0 shape.meta.m)
    =/  j  0
    |-  ^-  ray:ls
    ?:  =(j n)  m
    =/  mpj  (gi m ~[p j])
    =/  mqj  (gi m ~[q j])
    =.  m  (si m ~[p j] (csub cb (cmul cb cc mpj) (cmul cb b mqj)))
    =.  m  (si m ~[q j] (cadd cb (cmul cb (cconj cb b) mpj) (cmul cb cc mqj)))
    $(j +(j))
  ::    +sweep-herm: one cyclic Hermitian-Jacobi sweep over every p<q pair.
  ++  sweep-herm
    |=  [m=ray:ls v=ray:ls]
    ^-  [ray:ls ray:ls]
    =/  cb  bloq.meta.m
    =/  rb  (cb-comp cb)
    =/  n  (snag 0 shape.meta.m)
    ?:  (^lte n 1)  [m v]
    =/  p  0
    =/  q  1
    |-  ^-  [ray:ls ray:ls]
    =/  apq  (gi m ~[p q])
    =/  mag  (cabs-re cb apq)
    =/  mv=[ray:ls ray:ls]
      ?:  =(mag (f0 rb))
        [m v]
      =/  app  (cre cb (gi m ~[p p]))
      =/  aqq  (cre cb (gi m ~[q q]))
      =/  tau  (fdiv rb (fsub rb aqq app) (fmul rb (f2 rb) mag))
      =/  t    (fdiv rb (fsign rb tau) (fadd rb (fabs rb tau) (fsqt rb (fadd rb (fmul rb tau tau) (f1 rb)))))
      =/  c    (fdiv rb (f1 rb) (fsqt rb (fadd rb (fmul rb t t) (f1 rb))))
      =/  s    (fmul rb t c)
      =/  phase  (cdiv cb apq (cpak cb mag (f0 rb)))
      =/  b      (cmul cb (cpak cb s (f0 rb)) phase)
      =/  cc     (cpak cb c (f0 rb))
      :-  (rot-rows-h (rot-cols-h m p q cc b) p q cc b)
      (rot-cols-h v p q cc b)
    =.  m  -.mv
    =.  v  +.mv
    ?:  =(+(q) n)
      ?:  =(+(p) (dec n))  [m v]
      $(p +(p), q (^add p 2))
    $(q +(q))
  ::    +eig-herm: Hermitian eig.  vals real (%i754), vecs unitary (%cplx).
  ++  eig-herm
    |=  a=ray:ls
    ^-  [vals=ray:ls vecs=ray:ls]
    =/  cb  bloq.meta.a
    ?>  ?|(=(5 cb) =(6 cb) =(7 cb) =(8 cb))
    =/  rb  (cb-comp cb)
    ::  rtol: same default + width guard as +eig, against the component (rb).
    =/  rtol  ?:(=(0x1 `@`rtol) `@r`(feps rb) rtol)
    ~|  'saloon eig (hermitian): rtol wider than the component; match +sake width to the component bloq'
    ?>  (^lte (met 3 `@`rtol) (bex (^sub rb 3)))
    ?>  =(2 (lent shape.meta.a))
    =/  n  (snag 0 shape.meta.a)
    ?>  =(n (snag 1 shape.meta.a))
    ?>  (hermitian a)
    =/  v  (eye:(lake rnd) [~[n n] cb %cplx ~])
    =/  m  a
    =/  thresh  (fmul rb `@`rtol (frob-h a))
    =/  sweep  0
    |-  ^-  [vals=ray:ls vecs=ray:ls]
    ?:  =(60 sweep)
      ~&  "saloon eig (hermitian): hit sweep cap (60) without converging to rtol"
      [(diag-real m) v]
    ?:  (flte rb (off-norm-h m) thresh)
      [(diag-real m) v]
    =/  mv  (sweep-herm m v)
    $(sweep +(sweep), m -.mv, v +.mv)
  ::    +eig:  $ray -> [vals=$ray vecs=$ray]
  ::
  ::  Returns the eigenvalues (1-D ray) and eigenvectors (columns of a square
  ::  ray) of a symmetric (%i754) or Hermitian (%cplx) matrix, via cyclic
  ::  Jacobi.  Dispatches on kind; the %cplx path returns real (%i754)
  ::  eigenvalues and %cplx (unitary) eigenvectors.
  ::
  ::  Asserts squareness and symmetry/Hermitian-ness WITHIN a relative tolerance
  ::  (+near/+cnear) -- a matrix symmetric only up to rounding (e.g. a Gram
  ::  matrix from +mmul) is accepted; a genuinely asymmetric one crashes.
  ::
  ::  CALLERS: prefer +sake to set rtol, matched in WIDTH to the component (@rs
  ::  for @cs etc.).  A tolerance WIDER than the component is now rejected (its
  ::  high bytes would be dropped, mis-scaling the threshold).  The bare `sa`
  ::  default rtol=0x1 (a denormal, not 1.0) is replaced by a width-appropriate
  ::  epsilon (+feps) so bare +eig still converges, though +sake is recommended.
  ::    Examples
  ::      > =sa  (sake %n .~1e-12)
  ::      > =a   (en-ray:la [[~[2 2] 6 %i754 ~] ~[~[.~2 .~1] ~[.~1 .~2]]])  ::  [[2 1] [1 2]]
  ::      > ;;((list @rd) data:(de-ray:la -:(eig:sa a)))
  ::      ~[.~1 .~2.9999999999999996]                          ::  eigenvalues {1,3}
  ::  Source
  ++  eig
    |=  a=ray:ls
    ^-  [vals=ray:ls vecs=ray:ls]
    ?:  =(%cplx kind.meta.a)  (eig-herm a)
    =/  b  bloq.meta.a
    ?>  ?|(=(4 b) =(5 b) =(6 b) =(7 b))
    ::  rtol: replace the unusable bare default (0x1, a denormal) with a width-
    ::  appropriate epsilon; reject a tolerance wider than the component, whose
    ::  high bytes would be dropped and silently mis-scale the threshold.
    =/  rtol  ?:(=(0x1 `@`rtol) `@r`(feps b) rtol)
    ~|  'saloon eig: rtol wider than the %i754 component; match +sake width to the array bloq'
    ?>  (^lte (met 3 `@`rtol) (bex (^sub b 3)))
    ?>  =(2 (lent shape.meta.a))
    =/  n  (snag 0 shape.meta.a)
    ?>  =(n (snag 1 shape.meta.a))
    ?>  (symmetric a)
    =/  v  (eye:(lake rnd) [~[n n] b %i754 ~])
    =/  m  a
    =/  thresh  (fmul b `@`rtol (frob a))
    =/  sweep  0
    |-  ^-  [vals=ray:ls vecs=ray:ls]
    ?:  =(60 sweep)
      ~&  "saloon eig: hit sweep cap (60) without converging to rtol"
      [(diag m) v]
    ?:  (flte b (off-norm m) thresh)
      [(diag m) v]
    =/  mv  (sweep-once m v)
    $(sweep +(sweep), m -.mv, v +.mv)
  ::    +eigvals:  $ray -> $ray
  ::
  ::  Returns just the eigenvalues (1-D ray) of a symmetric/Hermitian matrix;
  ::  real (%i754) even for a Hermitian (%cplx) input.  (Head of +eig, taken by
  ::  axis because the kind-dispatch ?: in +eig erodes the result's faces.)
  ::    Examples
  ::      > =sa  (sake %n .~1e-12)
  ::      > =a   (en-ray:la [[~[2 2] 6 %i754 ~] ~[~[.~2 .~1] ~[.~1 .~2]]])
  ::      > ;;((list @rd) data:(de-ray:la (eigvals:sa a)))
  ::      ~[.~1 .~2.9999999999999996]                          ::  eigenvalues {1,3}
  ::  Source
  ++  eigvals  |=(a=ray:ls ^-(ray:ls -:(eig a)))
  ::    +eigvecs:  $ray -> $ray
  ::
  ::  Returns just the eigenvectors as the columns of a square ray (orthonormal
  ::  for symmetric input, unitary for Hermitian).  (Tail of +eig, by axis.)
  ::    Examples
  ::      > =sa  (sake %n .~1e-12)
  ::      > =a   (en-ray:la [[~[2 2] 6 %i754 ~] ~[~[.~2 .~1] ~[.~1 .~2]]])
  ::      > shape.meta:(eigvecs:sa a)
  ::      ~[2 2]
  ::  Source
  ++  eigvecs  |=(a=ray:ls ^-(ray:ls +:(eig a)))
  --
--
