/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-cplx -- complex (@cs) arrays, /lib/complex
::
::  @cs = bloq 6 (complex-single, two @rs); element packs real low, imag high.
::  Build rays with `fill`, read back with get-item.  Expected bit patterns
::  from NumPy complex64.  1+2i=0x4000.0000.3f80.0000, 3+4i=0x4080.0000.4040.0000,
::  one (1+0i)=0x3f80.0000, zero=0x0.
::
^|
|%
::  apply a binary ray op to two @cs scalars, read the result scalar
++  bin
  |=  [op=$-([ray ray] ray) b=@ c=@]  ^-  @
  =/  m1=meta  [~[1 1] 6 %cplx ~]
  (get-item:la (op (fill:la m1 b) (fill:la m1 c)) ~[0 0])
::
::  element-wise arithmetic (per-component IEEE single)
++  test-cplx-arith  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x40c0.0000.4080.0000)
      !>((bin add:la 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))   ::  (1+2i)+(3+4i)=4+6i
    %+  expect-eq  !>(`@`0xc000.0000.c000.0000)
      !>((bin sub:la 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))   ::  =-2-2i
    %+  expect-eq  !>(`@`0x4120.0000.c0a0.0000)
      !>((bin mul:la 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))   ::  =-5+10i
    %+  expect-eq  !>(`@`0x4000.0000)
      !>((bin div:la 0x4000.0000.4000.0000 0x3f80.0000.3f80.0000))   ::  (2+2i)/(1+1i)=2
  ==
::  abs (modulus, real-valued complex) and conj
++  test-cplx-unary  ^-  tang
  =/  m1=meta  [~[1 1] 6 %cplx ~]
  ;:  weld
    %+  expect-eq  !>(`@`0x40a0.0000)            !>((get-item:la (abs:la (fill:la m1 0x4080.0000.4040.0000)) ~[0 0]))   ::  |3+4i|=5
    %+  expect-eq  !>(`@`0xc000.0000.3f80.0000)  !>((get-item:la (conj:la (fill:la m1 0x4000.0000.3f80.0000)) ~[0 0]))  ::  conj(1+2i)=1-2i
  ==
::  comparisons return complex 1 (one=0x3f80.0000) / 0 (zero=0x0)
++  test-cplx-compare  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3f80.0000)  !>((bin equ:la 0x4000.0000.3f80.0000 0x4000.0000.3f80.0000))  ::  eq -> 1
    %+  expect-eq  !>(`@`0x0)          !>((bin equ:la 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))  ::  neq -> 0
    %+  expect-eq  !>(`@`0x3f80.0000)  !>((bin neq:la 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))  ::  neq -> 1
  ==
::  ones builder uses the %cplx constant 1+0i = 0x3f80.0000
++  test-cplx-ones  ^-  tang
  =/  m=meta  [~[2 2] 6 %cplx ~]
  %+  expect-eq  !>(`@`0x3f80.0000)  !>((get-item:la (ones:la m) ~[0 0]))
::  cumsum over [1+2i 3+4i -2-2i] = 2+4i
++  test-cplx-cumsum  ^-  tang
  =/  r=ray  (fill:la [~[1 3] 6 %cplx ~] 0x0)
  =.  r  (set-item:la r ~[0 0] 0x4000.0000.3f80.0000)   ::  1+2i
  =.  r  (set-item:la r ~[0 1] 0x4080.0000.4040.0000)   ::  3+4i
  =.  r  (set-item:la r ~[0 2] 0xc000.0000.c000.0000)   ::  -2-2i
  %+  expect-eq  !>(`@`0x4080.0000.4000.0000)  !>((get-item:la (cumsum:la r) ~[0 0]))   ::  2+4i
::  bilinear dot vs Hermitian dotc over [1+2i 3+4i].[5+6i 7+8i]
::  dot  = sum x*y       = -18+68i (0x4288.0000.c190.0000)
::  dotc = sum conj(x)*y =  70-8i  (0xc100.0000.428c.0000)
++  test-cplx-dot  ^-  tang
  =/  m=meta  [~[1 2] 6 %cplx ~]
  =/  a=ray  (fill:la m 0x0)
  =.  a  (set-item:la a ~[0 0] 0x4000.0000.3f80.0000)   ::  1+2i
  =.  a  (set-item:la a ~[0 1] 0x4080.0000.4040.0000)   ::  3+4i
  =/  b=ray  (fill:la m 0x0)
  =.  b  (set-item:la b ~[0 0] 0x40c0.0000.40a0.0000)   ::  5+6i
  =.  b  (set-item:la b ~[0 1] 0x4100.0000.40e0.0000)   ::  7+8i
  ;:  weld
    %+  expect-eq  !>(`@`0x4288.0000.c190.0000)  !>((get-item:la (dot:la a b) ~[0 0]))
    %+  expect-eq  !>(`@`0xc100.0000.428c.0000)  !>((get-item:la (dotc:la a b) ~[0 0]))
  ==
::  complex matrix multiply [[1+1i 2+0i][0 1+0i]] x [[1+0i 0][0 1+1i]]
::  = [[1+1i 2+2i][0 1+1i]]
++  test-cplx-mmul  ^-  tang
  =/  m=meta  [~[2 2] 6 %cplx ~]
  =/  a=ray  (fill:la m 0x0)
  =.  a  (set-item:la a ~[0 0] 0x3f80.0000.3f80.0000)   ::  1+1i
  =.  a  (set-item:la a ~[0 1] 0x4000.0000)             ::  2+0i
  =.  a  (set-item:la a ~[1 0] 0x0)                     ::  0
  =.  a  (set-item:la a ~[1 1] 0x3f80.0000)             ::  1+0i
  =/  b=ray  (fill:la m 0x0)
  =.  b  (set-item:la b ~[0 0] 0x3f80.0000)             ::  1+0i
  =.  b  (set-item:la b ~[0 1] 0x0)                     ::  0
  =.  b  (set-item:la b ~[1 0] 0x0)                     ::  0
  =.  b  (set-item:la b ~[1 1] 0x3f80.0000.3f80.0000)   ::  1+1i
  =/  c=ray  (mmul:la a b)
  ;:  weld
    %+  expect-eq  !>(`@`0x3f80.0000.3f80.0000)  !>((get-item:la c ~[0 0]))   ::  1+1i
    %+  expect-eq  !>(`@`0x4000.0000.4000.0000)  !>((get-item:la c ~[0 1]))   ::  2+2i
    %+  expect-eq  !>(`@`0x0)                     !>((get-item:la c ~[1 0]))   ::  0
    %+  expect-eq  !>(`@`0x3f80.0000.3f80.0000)  !>((get-item:la c ~[1 1]))   ::  1+1i
  ==
::  apply a binary ray op to two scalars of the given complex bloq
++  binw
  |=  [w=@ op=$-([ray ray] ray) b=@ c=@]  ^-  @
  =/  m1=meta  [~[1 1] w %cplx ~]
  (get-item:la (op (fill:la m1 b) (fill:la m1 c)) ~[0 0])
::  @ch (bloq 5): 1+2i=0x4000.3c00, 3+4i=0x4400.4200
++  test-cplx-ch  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4600.4400)  !>((binw 5 add:la 0x4000.3c00 0x4400.4200))   ::  =4+6i
    %+  expect-eq  !>(`@`0x4900.c500)  !>((binw 5 mul:la 0x4000.3c00 0x4400.4200))   ::  =-5+10i
    %+  expect-eq  !>(`@`0x4500)
      !>((get-item:la (abs:la (fill:la [~[1 1] 5 %cplx ~] 0x4400.4200)) ~[0 0]))     ::  |3+4i|=5
  ==
::  @cq (bloq 8): 1+2i / 3+4i in quad
++  test-cplx-cq  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4001.8000.0000.0000.0000.0000.0000.0000.4001.0000.0000.0000.0000.0000.0000.0000)
      !>  %-  binw
      :*  8  add:la
          0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000
          0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000           ::  =4+6i
      ==
    %+  expect-eq  !>(`@`0x4001.4000.0000.0000.0000.0000.0000.0000)
      !>((get-item:la (abs:la (fill:la [~[1 1] 8 %cplx ~] 0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000)) ~[0 0]))   ::  |3+4i|=5
  ==
--
