/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-fixp -- fixed-point (Q a.b) arrays, /lib/fixed
::
::  Element width is a power of two: a+b+1 = 2^bloq.  Here q3.4 (a=3 b=4,
::  N=8, bloq 3); the precision [a b] is carried in meta.tail.  Build rays
::  with `fill`, read back with get-item.
::
::  q3.4 scale = 2^4 = 16.  1.0=0x10 1.5=0x18 2.0=0x20 0.5=0x8 0.25=0x4
::  0.75=0x0c 3.0=0x30 4.0=0x40 3.5=0x38; -1.0=0xf0 -1.5=0xe8 -2.0=0xe0.
::  one (1.0) = 0x10.
::
^|
|%
::  apply a binary ray op to two q3.4 scalars, read the result scalar
++  bin
  |=  [op=$-([ray ray] ray) b=@ c=@]  ^-  @
  =/  m1=meta  [~[1 1] 3 %fixp [3 4]]
  (get-item:la (op (fill:la m1 b) (fill:la m1 c)) ~[0 0])
::
::  element-wise arithmetic (signed, rescaled for mul/div)
++  test-fixp-arith  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x18)  !>((bin add:la 0x10 0x8))   ::  1.0+0.5=1.5
    %+  expect-eq  !>(`@`0xf0)  !>((bin sub:la 0x10 0x20))   ::  1.0-2.0=-1.0
    %+  expect-eq  !>(`@`0x30)  !>((bin mul:la 0x18 0x20))   ::  1.5*2.0=3.0
    %+  expect-eq  !>(`@`0x4)  !>((bin mul:la 0x8 0x8))   ::  0.5*0.5=0.25
    %+  expect-eq  !>(`@`0x18)  !>((bin div:la 0x30 0x20))   ::  3.0/2.0=1.5
    %+  expect-eq  !>(`@`0x18)  !>((bin mod:la 0x38 0x20))   ::  3.5 mod 2.0=1.5
  ==
::  comparisons return fixed 1.0 (one=0x10) / 0.0 (0x0)
++  test-fixp-compare  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x10)  !>((bin gth:la 0x18 0x10))   ::  1.5 > 1.0 -> 1
    %+  expect-eq  !>(`@`0x0)   !>((bin gth:la 0x8 0x10))   ::  0.5 > 1.0 -> 0
    %+  expect-eq  !>(`@`0x10)  !>((bin lth:la 0xe8 0x10))   ::  -1.5 < 1.0 -> 1
    %+  expect-eq  !>(`@`0x10)  !>((bin gte:la 0x10 0x10))   ::  1.0 >= 1.0 -> 1
    %+  expect-eq  !>(`@`0x10)  !>((bin equ:la 0x10 0x10))   ::  1.0 == 1.0 -> 1
    %+  expect-eq  !>(`@`0x0)   !>((bin equ:la 0x10 0x18))   ::  1.0 == 1.5 -> 0
  ==
::  unary abs over a filled ray
++  test-fixp-abs  ^-  tang
  =/  m1=meta  [~[1 1] 3 %fixp [3 4]]
  ;:  weld
    %+  expect-eq  !>(`@`0x18)  !>((get-item:la (abs:la (fill:la m1 0xe8)) ~[0 0]))  ::  |-1.5|
    %+  expect-eq  !>(`@`0x18)  !>((get-item:la (abs:la (fill:la m1 0x18)) ~[0 0]))  ::  |1.5|
  ==
::  ones builder uses the %fixp constant 1.0 = 2^b = 0x10
++  test-fixp-ones  ^-  tang
  =/  m=meta  [~[2 2] 3 %fixp [3 4]]
  %+  expect-eq  !>(`@`0x10)  !>((get-item:la (ones:la m) ~[0 0]))
::  cumsum over [0.5 1.0 1.5 -1.0] = 2.0 (exact integer sum)
++  test-fixp-cumsum  ^-  tang
  =/  r=ray  (fill:la [~[1 4] 3 %fixp [3 4]] 0x0)
  =.  r  (set-item:la r ~[0 0] 0x8)              ::  0.5
  =.  r  (set-item:la r ~[0 1] 0x10)              ::  1.0
  =.  r  (set-item:la r ~[0 2] 0x18)              ::  1.5
  =.  r  (set-item:la r ~[0 3] 0xf0)              ::  -1.0
  %+  expect-eq  !>(`@`0x20)  !>((get-item:la (cumsum:la r) ~[0 0]))
::  exact dot product: 0.5*2.0 + 1.5*2.0 = 4.0 (0x40)
++  test-fixp-dot  ^-  tang
  =/  m=meta  [~[1 2] 3 %fixp [3 4]]
  =/  a=ray  (fill:la m 0x0)
  =.  a  (set-item:la a ~[0 0] 0x8)              ::  0.5
  =.  a  (set-item:la a ~[0 1] 0x18)              ::  1.5
  =/  b=ray  (fill:la m 0x0)
  =.  b  (set-item:la b ~[0 0] 0x20)              ::  2.0
  =.  b  (set-item:la b ~[0 1] 0x20)              ::  2.0
  %+  expect-eq  !>(`@`0x40)  !>((get-item:la (dot:la a b) ~[0 0]))
::  exact matrix multiply: [[1.0 0.5][0.0 2.0]] x [[2.0 0.0][0.0 0.5]]
::  = [[2.0 0.25][0.0 1.0]] -> 0x20 0x4 / 0x0 0x10
++  test-fixp-mmul  ^-  tang
  =/  m=meta  [~[2 2] 3 %fixp [3 4]]
  =/  a=ray  (fill:la m 0x0)
  =.  a  (set-item:la a ~[0 0] 0x10)              ::  1.0
  =.  a  (set-item:la a ~[0 1] 0x8)              ::  0.5
  =.  a  (set-item:la a ~[1 0] 0x0)              ::  0.0
  =.  a  (set-item:la a ~[1 1] 0x20)              ::  2.0
  =/  b=ray  (fill:la m 0x0)
  =.  b  (set-item:la b ~[0 0] 0x20)              ::  2.0
  =.  b  (set-item:la b ~[0 1] 0x0)              ::  0.0
  =.  b  (set-item:la b ~[1 0] 0x0)              ::  0.0
  =.  b  (set-item:la b ~[1 1] 0x8)              ::  0.5
  =/  c=ray  (mmul:la a b)
  ;:  weld
    %+  expect-eq  !>(`@`0x20)  !>((get-item:la c ~[0 0]))   ::  2.0
    %+  expect-eq  !>(`@`0x4)  !>((get-item:la c ~[0 1]))   ::  0.25
    %+  expect-eq  !>(`@`0x0)  !>((get-item:la c ~[1 0]))   ::  0.0
    %+  expect-eq  !>(`@`0x10)  !>((get-item:la c ~[1 1]))   ::  1.0
  ==
--
