  ::  /tests/lib/unum
::::
::    Posits (2022 Posit Standard, es=2)
::
::  Test strategy: the heavy exhaustive cross-check of arithmetic against
::  SoftPosit (all 65,536 posit8 pairs for add/sub/mul/div) is done ONCE,
::  offline, in Python.  The on-ship suite stays lean: exhaustive
::  decode/encode round-trips, a light oracle-free property sweep, and
::  curated SoftPosit-verified spot values.
::
/+  *test,
    unum
^|
|%
::  +rt: exhaustively check bit(sea(p)) == p for every n-bit pattern.
::
++  test-round-trip-rpb  ^-  tang
  =|  i=@
  |-  ^-  tang
  ?:  =(256 i)  ~
  =/  rt  (bit:rpb:unum (sea:rpb:unum i))
  ?:  =(i rt)  $(i +(i))
  [(cat 3 'rpb round-trip failed at ' (scot %ux i)) $(i +(i))]
::
++  test-round-trip-rph  ^-  tang
  =|  i=@
  |-  ^-  tang
  ?:  =(65.536 i)  ~
  =/  rt  (bit:rph:unum (sea:rph:unum i))
  ?:  =(i rt)  $(i +(i))
  [(cat 3 'rph round-trip failed at ' (scot %ux i)) $(i +(i))]
::
::  posit8 spot values (es=2): hand-derived canonical bit patterns.
::
++  test-values-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x0)   !>(zero:rpb:unum)
    %+  expect-eq  !>(`@`0x80)  !>(nar:rpb:unum)
    %+  expect-eq  !>(`@`0x40)  !>(one:rpb:unum)             ::  1.0
    %+  expect-eq  !>(`@`0x7f)  !>(maxpos:rpb:unum)          ::  2^24
    %+  expect-eq  !>(`@`0x1)   !>(minpos:rpb:unum)          ::  2^-24
  ==
::
::  Mathematical constants, cross-checked against SoftPosit (pX2) and an
::  independent reference encoder at posit8/16/32.
::
++  test-consts-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4d)  !>(pi:rpb:unum)
    %+  expect-eq  !>(`@`0x55)  !>(tau:rpb:unum)
    %+  expect-eq  !>(`@`0x4b)  !>(e:rpb:unum)
    %+  expect-eq  !>(`@`0x45)  !>(phi:rpb:unum)
    %+  expect-eq  !>(`@`0x43)  !>(sqt2:rpb:unum)
    %+  expect-eq  !>(`@`0x3b)  !>(invsqt2:rpb:unum)
    %+  expect-eq  !>(`@`0x3b)  !>(log2:rpb:unum)
    %+  expect-eq  !>(`@`0x44)  !>(invlog2:rpb:unum)
    %+  expect-eq  !>(`@`0x49)  !>(log10:rpb:unum)
  ==
::
++  test-consts-rph  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4c91)  !>(pi:rph:unum)
    %+  expect-eq  !>(`@`0x5491)  !>(tau:rph:unum)
    %+  expect-eq  !>(`@`0x4ae0)  !>(e:rph:unum)
    %+  expect-eq  !>(`@`0x44f2)  !>(phi:rph:unum)
    %+  expect-eq  !>(`@`0x4350)  !>(sqt2:rph:unum)
    %+  expect-eq  !>(`@`0x3b50)  !>(invsqt2:rph:unum)
    %+  expect-eq  !>(`@`0x3b17)  !>(log2:rph:unum)
    %+  expect-eq  !>(`@`0x438b)  !>(invlog2:rph:unum)
    %+  expect-eq  !>(`@`0x4936)  !>(log10:rph:unum)
  ==
::
++  test-consts-rps  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4c90.fdaa)  !>(pi:rps:unum)
    %+  expect-eq  !>(`@`0x5490.fdaa)  !>(tau:rps:unum)
    %+  expect-eq  !>(`@`0x4adf.8546)  !>(e:rps:unum)
    %+  expect-eq  !>(`@`0x44f1.bbce)  !>(phi:rps:unum)
    %+  expect-eq  !>(`@`0x4350.4f33)  !>(sqt2:rps:unum)
    %+  expect-eq  !>(`@`0x3b50.4f33)  !>(invsqt2:rps:unum)
    %+  expect-eq  !>(`@`0x3b17.217f)  !>(log2:rps:unum)
    %+  expect-eq  !>(`@`0x438a.a3b3)  !>(invlog2:rps:unum)
    %+  expect-eq  !>(`@`0x4935.d8de)  !>(log10:rps:unum)
  ==
::
::  posit8 decode (sea) / encode (bit) spot checks against the g-layer form.
::
++  test-sea-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>([%z ~])             !>((sea:rpb:unum 0x0))
    %+  expect-eq  !>([%n ~])             !>((sea:rpb:unum 0x80))
    %+  expect-eq  !>([%p %.y -3 8])      !>((sea:rpb:unum 0x40))   ::  1.0 (8*2^-3)
    %+  expect-eq  !>([%p %.y -2 8])      !>((sea:rpb:unum 0x48))   ::  2.0
    %+  expect-eq  !>([%p %.y -3 10])     !>((sea:rpb:unum 0x42))   ::  1.25
  ==
::
++  test-bit-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x0)   !>((bit:rpb:unum [%z ~]))
    %+  expect-eq  !>(`@`0x80)  !>((bit:rpb:unum [%n ~]))
    %+  expect-eq  !>(`@`0x40)  !>((bit:rpb:unum [%p %.y --0 1]))   ::  1.0
    %+  expect-eq  !>(`@`0x48)  !>((bit:rpb:unum [%p %.y -2 8]))    ::  2.0
    %+  expect-eq  !>(`@`0xc0)  !>((bit:rpb:unum [%p %.n --0 1]))   ::  -1.0
  ==
::
::  comparisons and sign ops (signed two's-complement ordering, NaR lowest).
::
++  test-compare-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((lth:rpb:unum 0x38 0x40))   ::  0.5 < 1.0
    %+  expect-eq  !>(%.n)  !>((lth:rpb:unum 0x40 0x38))
    %+  expect-eq  !>(%.y)  !>((gth:rpb:unum 0x40 0xc0))   ::  1.0 > -1.0
    %+  expect-eq  !>(%.y)  !>((lth:rpb:unum 0x80 0xc0))   ::  NaR < -1.0
    %+  expect-eq  !>(%.y)  !>((equ:rpb:unum 0x80 0x80))   ::  NaR == NaR
    %+  expect-eq  !>(`@`0xc0)  !>((neg:rpb:unum 0x40))    ::  -(1.0)
    %+  expect-eq  !>(`@`0x40)  !>((abs:rpb:unum 0xc0))    ::  |-1.0|
    %+  expect-eq  !>(`@`0xc0)  !>((sgn:rpb:unum 0xa0))    ::  sign of a negative
  ==
::
::  Arithmetic spot checks (SoftPosit pX2, es=2).
::
++  test-arith-spot-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x48)  !>((add:rpb:unum 0x40 0x40))  ::  1+1=2
    %+  expect-eq  !>(`@`0x40)  !>((add:rpb:unum 0x38 0x38))  ::  .5+.5=1
    %+  expect-eq  !>(`@`0x4c)  !>((add:rpb:unum 0x40 0x48))  ::  1+2=3
    %+  expect-eq  !>(`@`0x48)  !>((sub:rpb:unum 0x4c 0x40))  ::  3-1=2
    %+  expect-eq  !>(`@`0x54)  !>((mul:rpb:unum 0x48 0x4c))  ::  2*3=6
    %+  expect-eq  !>(`@`0x38)  !>((div:rpb:unum 0x40 0x48))  ::  1/2=.5
    %+  expect-eq  !>(`@`0x80)  !>((div:rpb:unum 0x40 0x0))   ::  1/0=NaR
    %+  expect-eq  !>(`@`0x80)  !>((mul:rpb:unum 0x40 0x80))  ::  1*NaR=NaR
  ==
::
++  test-arith-spot-rph  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4c00)  !>((add:rph:unum 0x4000 0x4800))  ::  1+2=3
    %+  expect-eq  !>(`@`0x319a)  !>((mul:rph:unum 0x4c00 0x24cd))  ::  3*0.1
  ==
::
::  Light algebraic property sweep over posit8 (256 pseudo-random pairs):
::  add/mul commute, a+0=a, a*1=a, a*(-1)=neg(a).  No oracle needed; the
::  exhaustive arithmetic check lives in the offline Python harness.
::
++  test-arith-props-rpb  ^-  tang
  =/  mu  mul:rpb:unum
  =/  ad  add:rpb:unum
  =/  ng  neg:rpb:unum
  =/  on  one:rpb:unum
  =/  no  (ng on)
  =|  i=@
  |-  ^-  tang
  ?:  =(256 i)  ~
  =/  a  i
  =/  b  (mod (add (mul i 181) 67) 256)            :: pseudo-random partner
  =/  e1=tang  ?:(=((mu a b) (mu b a)) ~ ~[(cat 3 'mul comm at ' (scot %ux i))])
  =/  e2=tang  ?:(=((ad a b) (ad b a)) ~ ~[(cat 3 'add comm at ' (scot %ux i))])
  =/  e3=tang  ?:(=(a (ad a 0x0)) ~ ~[(cat 3 'a+0 at ' (scot %ux i))])
  =/  e4=tang  ?:(=(a (mu a on)) ~ ~[(cat 3 'a*1 at ' (scot %ux i))])
  =/  e5=tang  ?:(=((ng a) (mu a no)) ~ ~[(cat 3 'a*-1 at ' (scot %ux i))])
  :(weld e1 e2 e3 e4 e5 $(i +(i)))
::
::  sqrt, rounding, conversions, fma (SoftPosit pX2-verified, es=2).
::
++  test-sqt-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x48)  !>((sqt:rpb:unum 0x50))   ::  sqrt(4)=2
    %+  expect-eq  !>(`@`0x43)  !>((sqt:rpb:unum 0x48))   ::  sqrt(2)=sqt2
    %+  expect-eq  !>(`@`0x80)  !>((sqt:rpb:unum 0xc0))   ::  sqrt(-1)=NaR
    %+  expect-eq  !>(`@`0x0)   !>((sqt:rpb:unum 0x0))    ::  sqrt(0)=0
    %+  expect-eq  !>(`@`0x80)  !>((sqt:rpb:unum 0x80))   ::  sqrt(NaR)=NaR
  ==
::
++  test-round-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x40)  !>((rnd:rpb:unum 0x42))   ::  1.25 -> 1
    %+  expect-eq  !>(`@`0x48)  !>((rnd:rpb:unum 0x4a))   ::  2.5 -> 2 (even)
    %+  expect-eq  !>(`@`0x0)   !>((rnd:rpb:unum 0x38))   ::  0.5 -> 0 (even)
    %+  expect-eq  !>(`@`0x50)  !>((rnd:rpb:unum 0x4e))   ::  3.5 -> 4 (even)
    %+  expect-eq  !>(`@`0x40)  !>((flr:rpb:unum 0x42))   ::  floor 1.25 = 1
    %+  expect-eq  !>(`@`0x48)  !>((cel:rpb:unum 0x42))   ::  ceil 1.25 = 2
    %+  expect-eq  !>(`@`0xb8)  !>((flr:rpb:unum 0xbe))   ::  floor -1.25 = -2
    %+  expect-eq  !>(`@`0xc0)  !>((cel:rpb:unum 0xbe))   ::  ceil -1.25 = -1
  ==
::
++  test-convert-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4c)  !>((sun:rpb:unum 3))      ::  3 -> 3.0
    %+  expect-eq  !>(`@`0x0)   !>((sun:rpb:unum 0))
    %+  expect-eq  !>(`@`0x4c)  !>((san:rpb:unum --3))    ::  +3 -> 3.0
    %+  expect-eq  !>(`@`0xb4)  !>((san:rpb:unum -3))     ::  -3 -> -3.0
    %+  expect-eq  !>(`(unit @s)`[~ --1])  !>((toi:rpb:unum 0x42))  ::  1.25 -> 1
    %+  expect-eq  !>(`(unit @s)`[~ --4])  !>((toi:rpb:unum 0x4e))  ::  3.5 -> 4
    %+  expect-eq  !>(`(unit @s)`~)        !>((toi:rpb:unum 0x80))  ::  NaR -> ~
    %+  expect-eq  !>(`(unit @s)`[~ --0])  !>((toi:rpb:unum 0x0))
  ==
::
++  test-fma-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x56)  !>((fma:rpb:unum 0x48 0x4c 0x40))  ::  2*3+1=7
    %+  expect-eq  !>(`@`0x80)  !>((fma:rpb:unum 0x48 0x80 0x40))  ::  NaR
  ==
::
::  Quire (Phase 4): exact accumulation + fused dot product (SoftPosit qX2).
::
++  test-quire-rpb  ^-  tang
  ;:  weld
    ::  p -> q -> p round trips (exact)
    %+  expect-eq  !>(`@`0x42)  !>((q-to-p:rpb:unum (p-to-q:rpb:unum 0x42)))
    %+  expect-eq  !>(`@`0x38)  !>((q-to-p:rpb:unum (p-to-q:rpb:unum 0x38)))
    ::  q-add-p: 1 + 2 = 3
    %+  expect-eq  !>(`@`0x4c)
      !>((q-to-p:rpb:unum (q-add-p:rpb:unum (p-to-q:rpb:unum 0x40) 0x48)))
    ::  q-mul-add: 0 + 2*3 = 6
    %+  expect-eq  !>(`@`0x54)
      !>((q-to-p:rpb:unum (q-mul-add:rpb:unum q-zero:rpb:unum 0x48 0x4c)))
  ==
::
++  test-fdp-rpb  ^-  tang
  ;:  weld
    ::  [1,2] . [3,1] = 5
    %+  expect-eq  !>(`@`0x52)  !>((fdp:rpb:unum ~[0x40 0x48] ~[0x4c 0x40]))
    ::  [2,3] . [2,3] = 13
    %+  expect-eq  !>(`@`0x5d)  !>((fdp:rpb:unum ~[0x48 0x4c] ~[0x48 0x4c]))
    ::  exact accumulation: maxpos + 1 - maxpos = 1 (naive float sum gives 0)
    %+  expect-eq  !>(`@`0x40)
      !>((fdp:rpb:unum ~[0x7f 0x40 0x81] ~[0x40 0x40 0x40]))
  ==
--
