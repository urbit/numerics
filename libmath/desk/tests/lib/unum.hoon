  ::  /tests/lib/unum
::::
::    Posits (2022 Posit Standard, es=2)
::
/+  *test,
    unum
^|
|%
::  +rt: exhaustively check that bit(sea(p)) == p for every n-bit pattern.
::
::  Posits are non-redundant, so decode-then-encode is the identity on every
::  bit pattern.  This catches the overwhelming majority of decode/encode and
::  rounding bugs for posit8 (256 patterns) and posit16 (65,536 patterns).
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
::  posit8 decode (sea) spot checks against the g-layer form.
::
++  test-sea-rpb  ^-  tang
  ;:  weld
    %+  expect-eq  !>([%z ~])             !>((sea:rpb:unum 0x0))
    %+  expect-eq  !>([%n ~])             !>((sea:rpb:unum 0x80))
    %+  expect-eq  !>([%p %.y --0 1])     !>((sea:rpb:unum 0x40))   ::  1.0
    %+  expect-eq  !>([%p %.y -2 8])      !>((sea:rpb:unum 0x48))   ::  2.0
    %+  expect-eq  !>([%p %.y -3 10])     !>((sea:rpb:unum 0x42))   ::  1.25
  ==
::
::  posit8 encode (bit) spot checks.
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
--
