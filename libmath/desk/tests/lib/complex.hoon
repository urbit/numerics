/+  *test, complex
::::  /tests/lib/complex -- complex arithmetic over IEEE-754 components
::
::  @cs packs two @rs (real low, imag high); expected values from NumPy
::  complex64.  1+2i=0x4000.0000.3f80.0000, 3+4i=0x4080.0000.4040.0000.
::  @cd packs two @rd; values from NumPy complex128.  1+2i packs to
::  0x4000.0000.0000.0000.3ff0.0000.0000.0000, 3+4i to
::  0x4010.0000.0000.0000.4008.0000.0000.0000.
::
^|
|%
::  --- complex-single (@cs) ---
++  test-cs-components  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@rs`.1)  !>((~(re cs:complex %n) 0x4000.0000.3f80.0000))   ::  re(1+2i)=1
    %+  expect-eq  !>(`@rs`.2)  !>((~(im cs:complex %n) 0x4000.0000.3f80.0000))   ::  im(1+2i)=2
    %+  expect-eq  !>(`@`0x4000.0000.3f80.0000)  !>((~(pak cs:complex %n) .1 .2))  ::  1+2i
    %+  expect-eq  !>(`@`0x0)  !>(~(zero cs:complex %n))
  ==
++  test-cs-arith  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x40c0.0000.4080.0000)
      !>((~(add cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))   ::  (1+2i)+(3+4i)=4+6i
    %+  expect-eq  !>(`@`0xc000.0000.c000.0000)
      !>((~(sub cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))   ::  =-2-2i
    %+  expect-eq  !>(`@`0x4120.0000.c0a0.0000)
      !>((~(mul cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))   ::  =-5+10i
    %+  expect-eq  !>(`@`0x4000.0000)
      !>((~(div cs:complex %n) 0x4000.0000.4000.0000 0x3f80.0000.3f80.0000))   ::  (2+2i)/(1+1i)=2
  ==
++  test-cs-unary  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0xc000.0000.bf80.0000)  !>((~(neg cs:complex %n) 0x4000.0000.3f80.0000))   ::  -(1+2i)
    %+  expect-eq  !>(`@`0xc000.0000.3f80.0000)  !>((~(conj cs:complex %n) 0x4000.0000.3f80.0000))  ::  conj=1-2i
    %+  expect-eq  !>(`@`0x40a0.0000)            !>((~(abs cs:complex %n) 0x4080.0000.4040.0000))   ::  |3+4i|=5
    %+  expect-eq  !>(`@`0x3f80.0000)            !>(~(one cs:complex %n))                            ::  1+0i
  ==
++  test-cs-compare  ^-  tang
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(equ cs:complex %n) 0x4000.0000.3f80.0000 0x4000.0000.3f80.0000))
    %+  expect-eq  !>(%.n)  !>((~(equ cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))
    %+  expect-eq  !>(%.y)  !>((~(neq cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000))
  ==
::  --- complex-double (@cd) ---
++  test-cd-components  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@rd`.~1)  !>((~(re cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000))  ::  re(1+2i)=1
    %+  expect-eq  !>(`@rd`.~2)  !>((~(im cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000))  ::  im(1+2i)=2
    %+  expect-eq  !>(`@`0x4000.0000.0000.0000.3ff0.0000.0000.0000)  !>((~(pak cd:complex %n) .~1 .~2))
    %+  expect-eq  !>(`@`0x0)  !>(~(zero cd:complex %n))
  ==
++  test-cd-arith  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4018.0000.0000.0000.4010.0000.0000.0000)
      !>  %-  ~(add cd:complex %n)
          :-  0x4000.0000.0000.0000.3ff0.0000.0000.0000
          0x4010.0000.0000.0000.4008.0000.0000.0000                            ::  (1+2i)+(3+4i)=4+6i
    %+  expect-eq  !>(`@`0xc000.0000.0000.0000.c000.0000.0000.0000)
      !>  %-  ~(sub cd:complex %n)
          :-  0x4000.0000.0000.0000.3ff0.0000.0000.0000
          0x4010.0000.0000.0000.4008.0000.0000.0000                            ::  =-2-2i
    %+  expect-eq  !>(`@`0x4024.0000.0000.0000.c014.0000.0000.0000)
      !>  %-  ~(mul cd:complex %n)
          :-  0x4000.0000.0000.0000.3ff0.0000.0000.0000
          0x4010.0000.0000.0000.4008.0000.0000.0000                            ::  =-5+10i
    %+  expect-eq  !>(`@`0x4000.0000.0000.0000)
      !>  %-  ~(div cd:complex %n)
          :-  0x4000.0000.0000.0000.4000.0000.0000.0000
          0x3ff0.0000.0000.0000.3ff0.0000.0000.0000                            ::  (2+2i)/(1+1i)=2
  ==
++  test-cd-unary  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0xc000.0000.0000.0000.bff0.0000.0000.0000)
      !>((~(neg cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000))     ::  -(1+2i)
    %+  expect-eq  !>(`@`0xc000.0000.0000.0000.3ff0.0000.0000.0000)
      !>((~(conj cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000))    ::  conj=1-2i
    %+  expect-eq  !>(`@`0x4014.0000.0000.0000)
      !>((~(abs cd:complex %n) 0x4010.0000.0000.0000.4008.0000.0000.0000))     ::  |3+4i|=5
    %+  expect-eq  !>(`@`0x3ff0.0000.0000.0000)  !>(~(one cd:complex %n))      ::  1+0i
  ==
++  test-cd-compare  ^-  tang
  ;:  weld
    %+  expect-eq  !>(%.y)
      !>  %-  ~(equ cd:complex %n)
          :-  0x4000.0000.0000.0000.3ff0.0000.0000.0000
          0x4000.0000.0000.0000.3ff0.0000.0000.0000
    %+  expect-eq  !>(%.n)
      !>  %-  ~(equ cd:complex %n)
          :-  0x4000.0000.0000.0000.3ff0.0000.0000.0000
          0x4010.0000.0000.0000.4008.0000.0000.0000
    %+  expect-eq  !>(%.y)
      !>  %-  ~(neq cd:complex %n)
          :-  0x4000.0000.0000.0000.3ff0.0000.0000.0000
          0x4010.0000.0000.0000.4008.0000.0000.0000
  ==
--
