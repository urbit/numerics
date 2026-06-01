  ::  /tests/lib/unum-fns
::::
::    Posits (2022 Posit Standard, es=2) -- elementary functions, rounding,
::    integer/IEEE conversions, the quire, and transcendentals.  Split from a
::    single oversized file so each compiles quickly.
::
::  Note: the from-rh/rs/rd/rq gates are typed on the float aura, so their
::  arguments must be cast (`@rs`0x..) -- a bare 0x.. literal is @ux and will
::  not nest.  The to-* gates take a bare @, so their literals need no cast.
::
/+  *test,
    unum
^|
|%
++  test-sqt-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x48)  !>((sqt:rpb:unum 0x50))
    %+  expect-eq  !>(`@`0x43)  !>((sqt:rpb:unum 0x48))
    %+  expect-eq  !>(`@`0x80)  !>((sqt:rpb:unum 0xc0))
    %+  expect-eq  !>(`@`0x0)   !>((sqt:rpb:unum 0x0))
    %+  expect-eq  !>(`@`0x80)  !>((sqt:rpb:unum 0x80))
  ==
::
++  test-round-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x40)  !>((rnd:rpb:unum 0x42))
    %+  expect-eq  !>(`@`0x48)  !>((rnd:rpb:unum 0x4a))
    %+  expect-eq  !>(`@`0x0)   !>((rnd:rpb:unum 0x38))
    %+  expect-eq  !>(`@`0x50)  !>((rnd:rpb:unum 0x4e))
    %+  expect-eq  !>(`@`0x40)  !>((flr:rpb:unum 0x42))
    %+  expect-eq  !>(`@`0x48)  !>((cel:rpb:unum 0x42))
    %+  expect-eq  !>(`@`0xb8)  !>((flr:rpb:unum 0xbe))
    %+  expect-eq  !>(`@`0xc0)  !>((cel:rpb:unum 0xbe))
  ==
::
++  test-convert-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4c)  !>((sun:rpb:unum 3))
    %+  expect-eq  !>(`@`0x0)   !>((sun:rpb:unum 0))
    %+  expect-eq  !>(`@`0x4c)  !>((san:rpb:unum --3))
    %+  expect-eq  !>(`@`0xb4)  !>((san:rpb:unum -3))
    %+  expect-eq  !>(`(unit @s)`[~ --1])  !>((toi:rpb:unum 0x42))
    %+  expect-eq  !>(`(unit @s)`[~ --4])  !>((toi:rpb:unum 0x4e))
    %+  expect-eq  !>(`(unit @s)`~)        !>((toi:rpb:unum 0x80))
    %+  expect-eq  !>(`(unit @s)`[~ --0])  !>((toi:rpb:unum 0x0))
  ==
::
++  test-fma-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x56)  !>((fma:rpb:unum 0x48 0x4c 0x40))
    %+  expect-eq  !>(`@`0x80)  !>((fma:rpb:unum 0x48 0x80 0x40))
  ==
::
++  test-quire-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x42)  !>((q-to-p:rpb:unum (p-to-q:rpb:unum 0x42)))
    %+  expect-eq  !>(`@`0x38)  !>((q-to-p:rpb:unum (p-to-q:rpb:unum 0x38)))
    %+  expect-eq  !>(`@`0x4c)
      !>((q-to-p:rpb:unum (q-add-p:rpb:unum (p-to-q:rpb:unum 0x40) 0x48)))
    %+  expect-eq  !>(`@`0x54)
      !>((q-to-p:rpb:unum (q-mul-add:rpb:unum q-zero:rpb:unum 0x48 0x4c)))
  ==
::
++  test-fdp-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x52)  !>((fdp:rpb:unum ~[0x40 0x48] ~[0x4c 0x40]))
    %+  expect-eq  !>(`@`0x5d)  !>((fdp:rpb:unum ~[0x48 0x4c] ~[0x48 0x4c]))
    %+  expect-eq  !>(`@`0x40)
      !>((fdp:rpb:unum ~[0x7f 0x40 0x81] ~[0x40 0x40 0x40]))
  ==
::
::  posit32 <-> single (value-based; posit -2.0 = 0xb800.0000, NOT the
::  single -2.0 = 0xc000.0000 -- the two formats differ at the same width).
++  test-ieee-rps  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3f80.0000)  !>((to-rs:rps:unum 0x4000.0000))
    %+  expect-eq  !>(`@`0xc000.0000)  !>((to-rs:rps:unum 0xb800.0000))
    %+  expect-eq  !>(`@`0x3f00.0000)  !>((to-rs:rps:unum 0x3800.0000))
    %+  expect-eq  !>(`@`0x4000.0000)  !>((from-rs:rps:unum `@rs`0x3f80.0000))
    %+  expect-eq  !>(`@`0xb800.0000)  !>((from-rs:rps:unum `@rs`0xc000.0000))
    %+  expect-eq  !>(`@`0x3800.0000)  !>((from-rs:rps:unum `@rs`0x3f00.0000))
  ==
::
++  test-ieee-rph  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3c00)  !>((to-rh:rph:unum 0x4000))
    %+  expect-eq  !>(`@`0xc000)  !>((to-rh:rph:unum 0xb800))
    %+  expect-eq  !>(`@`0x4000)  !>((from-rh:rph:unum `@rh`0x3c00))
    %+  expect-eq  !>(`@`0xb800)  !>((from-rh:rph:unum `@rh`0xc000))
  ==
::
::  The matrix: ANY posit width <-> ANY float width (value-based).  posits
::  pack more accuracy per bit, so cross-width is the useful correspondence.
::
++  test-ieee-matrix  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3f80.0000)  !>((to-rs:rph:unum 0x4000))
    %+  expect-eq  !>(`@`0x3ff0.0000.0000.0000)  !>((to-rd:rps:unum 0x4000.0000))
    %+  expect-eq  !>(`@`0x3f80.0000)  !>((to-rs:rpb:unum 0x40))
    %+  expect-eq  !>(`@`0x4000)  !>((from-rs:rph:unum `@rs`0x3f80.0000))
    %+  expect-eq  !>(`@`0x4000.0000)  !>((from-rd:rps:unum `@rd`0x3ff0.0000.0000.0000))
  ==
::
::  Transcendentals (naive Taylor series, /lib/math style).  These posit8
::  values are the series outputs, which match the SoftPosit correctly-rounded
::  result (verified offline by a Python port of the same series).  pow-n is
::  exact integer power.
::
++  test-transcendental-rpb  ^-  tang
  =/  u  rpb:unum
  ;:  weld
    %+  expect-eq  !>(`@`0x40)  !>((exp:u 0x0))
    %+  expect-eq  !>(`@`0x4b)  !>((exp:u 0x40))
    %+  expect-eq  !>(`@`0x0)   !>((sin:u 0x0))
    %+  expect-eq  !>(`@`0x3d)  !>((sin:u 0x40))
    %+  expect-eq  !>(`@`0x40)  !>((cos:u 0x0))
    %+  expect-eq  !>(`@`0x39)  !>((cos:u 0x40))
    %+  expect-eq  !>(`@`0x44)  !>((tan:u 0x40))
    %+  expect-eq  !>(`@`0x0)   !>((log:u 0x40))
    %+  expect-eq  !>(`@`0x3b)  !>((log:u (sun:u 2)))
    %+  expect-eq  !>(`@`0x58)  !>((pow-n:u (sun:u 2) 3))
    %+  expect-eq  !>(`@`0x40)  !>((pow:u (sun:u 2) 0x0))
  ==
::
::  Domain guards: log/pow of out-of-domain inputs return NaR (not a divergent
::  series result), and NaR propagates through pow-n even when the exponent is 0.
::
++  test-domain-rpb  ^-  tang
  =/  u  rpb:unum
  ;:  weld
    %+  expect-eq  !>(`@`0x80)  !>((log:u 0x0))            ::  log(0)=NaR
    %+  expect-eq  !>(`@`0x80)  !>((log:u 0xc0))           ::  log(-1)=NaR
    %+  expect-eq  !>(`@`0x80)  !>((log:u 0x80))           ::  log(NaR)=NaR
    %+  expect-eq  !>(`@`0x80)  !>((log:u 0xb4))           ::  log(-3)=NaR
    %+  expect-eq  !>(`@`0x80)  !>((pow-n:u 0x80 0))       ::  NaR^0=NaR (not 1)
    %+  expect-eq  !>(`@`0x40)  !>((pow-n:u 0x40 0))       ::  1^0=1 still
    %+  expect-eq  !>(`@`0x80)  !>((pow:u 0x0 (sun:u 2)))  ::  pow(0,2)=NaR
    %+  expect-eq  !>(`@`0x80)  !>((pow:u 0xc0 (sun:u 2))) ::  pow(-1,2)=NaR
  ==
--
