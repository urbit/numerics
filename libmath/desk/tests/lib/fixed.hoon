  ::  /tests/lib/fixed
::::
::    Fixed-point arithmetic (signed, modular), built on /lib/twoc.
::    q8.8 = prec [8 8], N = 17.  Negative operands are the cases the old
::    unsigned mul/div/scale got wrong; the scale crash is also covered.
::
/+  *test,
    fixed
^|
=/  q88=prec:fixed  [8 8]
=/  n1-5  (neg:fixed 0x180 q88)   :: -1.5 in q8.8 = 0x1.fe80
=/  n3    (neg:fixed 0x300 q88)   :: -3.0 in q8.8 = 0x1.fd00
|%
++  test-add  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3c0)     !>((add:fixed 0x180 q88 0x240 q88))  ::  1.5+2.25=3.75
    %+  expect-eq  !>(`@`0x200)     !>((add:fixed 0x100 q88 0x100 q88))  ::  1+1=2
    %+  expect-eq  !>(`@`0x1.ff80)  !>((add:fixed 0x100 q88 n1-5 q88))   ::  1+(-1.5)=-0.5
  ==
::
++  test-sub  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x1.ff00)  !>((sub:fixed 0x100 q88 0x200 q88))  ::  1-2=-1
    %+  expect-eq  !>(`@`0x100)     !>((sub:fixed 0x200 q88 0x100 q88))  ::  2-1=1
  ==
::
++  test-neg  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x1.ff00)  !>((neg:fixed 0x100 q88))   ::  -(1.0)
    %+  expect-eq  !>(`@`0x1.fe80)  !>((neg:fixed 0x180 q88))   ::  -(1.5)
    %+  expect-eq  !>(`@`0x100)     !>((neg:fixed (neg:fixed 0x100 q88) q88))  ::  --1=1
  ==
::
::  mul widens precision to [17 16] (N=34).  Negative cases the old unsigned
::  mul returned garbage for.  -3.0 in q17.16 is 0x3.fffd.0000.
++  test-mul  ^-  tang
  =/  prodp=prec:fixed  [17 16]
  ;:  weld
    %+  expect-eq  !>([0x3.0000 prodp])       !>((mul:fixed 0x180 q88 0x200 q88))
    %+  expect-eq  !>([0x3.fffd.0000 prodp])  !>((mul:fixed n1-5 q88 0x200 q88))
    %+  expect-eq  !>([0x3.0000 prodp])       !>((mul:fixed n1-5 q88 (neg:fixed 0x200 q88) q88))
  ==
::
::  div keeps precision [8 8].  Negative cases the old unsigned div got wrong.
++  test-div  ^-  tang
  ;:  weld
    %+  expect-eq  !>([0x1.fe80 q88])  !>((div:fixed n3 q88 0x200 q88))   ::  -3/2=-1.5
    %+  expect-eq  !>([0x180 q88])     !>((div:fixed 0x300 q88 0x200 q88))  ::  3/2=1.5
    %+  expect-eq  !>([0x180 q88])     !>((div:fixed n3 q88 (neg:fixed 0x200 q88) q88))  ::  -3/-2=1.5
  ==
::
::  scale (the old version crashed on widen via a bad lsh; now correct).
::  -1.5 in q4.4 (N=9) is 0x1e8.
++  test-scale  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x100)  !>((scale:fixed 0x10 [4 4] [8 8]))   ::  1.0 q4.4 -> q8.8
    %+  expect-eq  !>(`@`0x10)   !>((scale:fixed 0x100 [8 8] [4 4]))  ::  1.0 q8.8 -> q4.4
    %+  expect-eq  !>(`@`0x1e8)  !>((scale:fixed n1-5 [8 8] [4 4]))   ::  -1.5 q8.8 -> q4.4
  ==
::
::  round-trip: (3.0/2.0) then *2.0 = 3.0 (exact here).
++  test-div-roundtrip  ^-  tang
  =/  q  (div:fixed 0x300 q88 0x200 q88)        ::  3.0/2.0 = 1.5 in q8.8
  =/  back  (mul:fixed -.q q88 0x200 q88)       ::  1.5*2.0 = 3.0 in q17.16
  %+  expect-eq  !>(`@`0x3.0000)  !>(-.back)
--
