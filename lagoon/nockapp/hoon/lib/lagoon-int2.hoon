/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-int2 -- two's-complement signed integer arrays
::
::  int8 = bloq 3, kind %int2.  Build rays with `fill`, read back with
::  get-item.  Verifies modular wrap, signed div/rem, abs, the 1=true/
::  0=false comparison convention, and a reduction (cumsum).
::
^|
|%
::  apply a binary ray op to two int8 scalars, read the result scalar
++  bin
  |=  [op=$-([ray ray] ray) b=@ c=@]  ^-  @
  =/  m1=meta  [~[1 1] 3 %int2 ~]
  (get-item:la (op (fill:la m1 b) (fill:la m1 c)) ~[0 0])
::
::  arithmetic, modular two's-complement (int8)
++  test-int2-arith  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x7)   !>((bin add:la 0x5 0x2))     ::  5+2=7
    %+  expect-eq  !>(`@`0x96)  !>((bin add:la 0x64 0x32))   ::  100+50 wraps -> -106
    %+  expect-eq  !>(`@`0xfb)  !>((bin sub:la 0x5 0xa))     ::  5-10=-5
    %+  expect-eq  !>(`@`0xf6)  !>((bin mul:la 0xfb 0x2))    ::  -5*2=-10
    %+  expect-eq  !>(`@`0xfd)  !>((bin div:la 0xf9 0x2))    ::  -7/2 trunc=-3
    %+  expect-eq  !>(`@`0xff)  !>((bin mod:la 0xf9 0x2))    ::  -7%2=-1
  ==
::  comparisons return value 1 (true) / 0 (false), matching %uint/%i754
++  test-int2-compare  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x1)  !>((bin lth:la 0xff 0x1))     ::  -1 < 1 -> 1
    %+  expect-eq  !>(`@`0x0)  !>((bin gth:la 0xff 0x1))     ::  -1 > 1 -> 0
    %+  expect-eq  !>(`@`0x1)  !>((bin gth:la 0x1 0xff))     ::  1 > -1 -> 1
    %+  expect-eq  !>(`@`0x1)  !>((bin gte:la 0x5 0x5))      ::  5 >= 5 -> 1
    %+  expect-eq  !>(`@`0x1)  !>((bin equ:la 0x5 0x5))      ::  5 == 5 -> 1
    %+  expect-eq  !>(`@`0x0)  !>((bin equ:la 0x5 0x6))      ::  5 == 6 -> 0
  ==
::  unary abs over a filled ray
++  test-int2-abs  ^-  tang
  =/  m1=meta  [~[1 1] 3 %int2 ~]
  ;:  weld
    %+  expect-eq  !>(`@`0x5)  !>((get-item:la (abs:la (fill:la m1 0xfb)) ~[0 0]))
    %+  expect-eq  !>(`@`0x5)  !>((get-item:la (abs:la (fill:la m1 0x5)) ~[0 0]))
  ==
::  ones builder uses the %int2 constant 1
++  test-int2-ones  ^-  tang
  =/  m=meta  [~[2 2] 3 %int2 ~]
  %+  expect-eq  !>(`@`0x1)  !>((get-item:la (ones:la m) ~[0 0]))
::  cumsum reduction over a 1x4 int2 ray [1 2 3 -1] = 5
++  test-int2-cumsum  ^-  tang
  =/  r=ray  (fill:la [~[1 4] 3 %int2 ~] 0x0)
  =.  r  (set-item:la r ~[0 0] 0x1)
  =.  r  (set-item:la r ~[0 1] 0x2)
  =.  r  (set-item:la r ~[0 2] 0x3)
  =.  r  (set-item:la r ~[0 3] 0xff)
  %+  expect-eq  !>(`@`0x5)  !>((get-item:la (cumsum:la r) ~[0 0]))
--
