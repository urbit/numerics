/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-unum -- posit (unum) arrays, /lib/unum
::
::  posit8 = bloq 3, kind %unum (rpb); posit16 = bloq 4 (rph).  Build rays
::  with `fill`, read back with get-item.  Expected bit patterns come from
::  the SoftPosit-cross-checked oracle in tools/posit_check.py (es=2).
::
::  Reference posit8 patterns: 1=0x40 2=0x48 3=0x4c 5=0x52 6=0x54 8=0x58
::  1.5=0x44 -2=0xb8 -5=0xae 1/4=0x30; zero=0x00 (one/zero = compare values).
::
^|
|%
::  apply a binary ray op to two posit8 scalars, read the result scalar
++  bin
  |=  [op=$-([ray ray] ray) b=@ c=@]  ^-  @
  =/  m1=meta  [~[1 1] 3 %unum ~]
  (get-item:la (op (fill:la m1 b) (fill:la m1 c)) ~[0 0])
::  same, for posit16
++  bin16
  |=  [op=$-([ray ray] ray) b=@ c=@]  ^-  @
  =/  m1=meta  [~[1 1] 4 %unum ~]
  (get-item:la (op (fill:la m1 b) (fill:la m1 c)) ~[0 0])
::
::  arithmetic, posit8 (correctly rounded, es=2)
++  test-unum-arith  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x52)  !>((bin add:la 0x48 0x4c))    ::  2+3=5
    %+  expect-eq  !>(`@`0xb8)  !>((bin sub:la 0x4c 0x52))    ::  3-5=-2
    %+  expect-eq  !>(`@`0x54)  !>((bin mul:la 0x48 0x4c))    ::  2*3=6
    %+  expect-eq  !>(`@`0x44)  !>((bin div:la 0x4c 0x48))    ::  3/2=1.5
    %+  expect-eq  !>(`@`0x30)  !>((bin div:la 0x40 0x50))    ::  1/4
  ==
::  mod (a - b*trunc(a/b)) via the ray op; pow via fun-scalar directly
::  (lagoon exposes no ray-level binary pow, only fun-scalar's %pow).
++  test-unum-mod-pow  ^-  tang
  =/  m1=meta  [~[1 1] 3 %unum ~]
  ;:  weld
    %+  expect-eq  !>(`@`0x40)  !>((bin mod:la 0x56 0x4c))            ::  7 mod 3 = 1
    %+  expect-eq  !>(`@`0x58)  !>(((fun-scalar:la m1 %pow) 0x48 0x4c))  ::  2^3 = 8
  ==
::  comparisons return posit 1 (one=0x40) / 0 (zero=0x00)
++  test-unum-compare  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x40)  !>((bin lth:la 0xae 0x40))    ::  -5 < 1 -> 1
    %+  expect-eq  !>(`@`0x0)  !>((bin gth:la 0xae 0x40))    ::  -5 > 1 -> 0
    %+  expect-eq  !>(`@`0x40)  !>((bin gth:la 0x40 0xae))    ::  1 > -5 -> 1
    %+  expect-eq  !>(`@`0x40)  !>((bin gte:la 0x48 0x48))    ::  2 >= 2 -> 1
    %+  expect-eq  !>(`@`0x40)  !>((bin equ:la 0x48 0x48))    ::  2 == 2 -> 1
    %+  expect-eq  !>(`@`0x0)  !>((bin equ:la 0x48 0x4c))    ::  2 == 3 -> 0
  ==
::  unary abs over a filled ray
++  test-unum-abs  ^-  tang
  =/  m1=meta  [~[1 1] 3 %unum ~]
  ;:  weld
    %+  expect-eq  !>(`@`0x52)  !>((get-item:la (abs:la (fill:la m1 0xae)) ~[0 0]))
    %+  expect-eq  !>(`@`0x52)  !>((get-item:la (abs:la (fill:la m1 0x52)) ~[0 0]))
  ==
::  ones builder uses the %unum constant one (0x40)
++  test-unum-ones  ^-  tang
  =/  m=meta  [~[2 2] 3 %unum ~]
  %+  expect-eq  !>(`@`0x40)  !>((get-item:la (ones:la m) ~[0 0]))
::  cumsum reduction over a 1x4 posit8 ray [2 3 -5 5] = 5
++  test-unum-cumsum  ^-  tang
  =/  r=ray  (fill:la [~[1 4] 3 %unum ~] 0x0)
  =.  r  (set-item:la r ~[0 0] 0x48)              ::  2
  =.  r  (set-item:la r ~[0 1] 0x4c)              ::  3
  =.  r  (set-item:la r ~[0 2] 0xae)              ::  -5
  =.  r  (set-item:la r ~[0 3] 0x52)              ::  5
  %+  expect-eq  !>(`@`0x52)  !>((get-item:la (cumsum:la r) ~[0 0]))
::  posit16 arithmetic (bloq 4, rph)
++  test-unum-posit16  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x5200)  !>((bin16 add:la 0x4800 0x4c00))   ::  2+3=5
    %+  expect-eq  !>(`@`0x5c00)  !>((bin16 mul:la 0x4c00 0x5000))   ::  3*4=12
  ==
--
