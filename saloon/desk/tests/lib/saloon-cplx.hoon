::  Saloon Tier-B transcendentals over %cplx (complex) rays.  Verifies the
::  trans-scalar/fun-scalar %cplx dispatch reaches /lib/complex per width.
::  Expected values are the same series outputs checked in /tests/lib/complex-fns
::  (validated vs numpy).  Inputs: z = 1+2i at @cs (bloq 6) and @ch (bloq 5).
::
/-  ls=lagoon
/+  *test, *saloon, *lagoon
|%
::  a 1-element %cplx ray at bloq .b holding packed complex pattern .p
++  ur  |=([b=@ p=@] ^-(ray:ls (en-ray:(lake %n) [[~[1] b %cplx ~] ~[p]])))
++  hd   |=(a=ray:ls ^-(@ -:(ravel:(lake %n) a)))
++  sad  (sake %n .1e-5)   :: %n rounding (bare +sa default rnd is not %n)
++  z6   `@`0x4000.0000.3f80.0000               ::  1+2i @cs
++  w26  `@`0x4000.0000                          ::  2+0i @cs
++  z5   `@`0x4000.3c00                          ::  1+2i @ch
++  test-cplx-exp-cs   (expect-eq !>(`@`0x401e.30c5.bf90.cb4e) !>((hd (exp:sad (ur 6 z6)))))
++  test-cplx-log-cs   (expect-eq !>(`@`0x3f8d.b70a.3f4e.0210) !>((hd (log:sad (ur 6 z6)))))
++  test-cplx-sqrt-cs  (expect-eq !>(`@`0x3f49.4138.3fa2.d18a) !>((hd (sqrt:sad (ur 6 z6)))))
++  test-cplx-sin-cs   (expect-eq !>(`@`0x3ffa.d435.404a.9c1e) !>((hd (sin:sad (ur 6 z6)))))
++  test-cplx-cos-cs   (expect-eq !>(`@`0xc043.524b.4002.1823) !>((hd (cos:sad (ur 6 z6)))))
++  test-cplx-tan-cs   (expect-eq !>(`@`0x3f81.e4c1.3d0a.7f5c) !>((hd (tan:sad (ur 6 z6)))))
++  test-cplx-pow-cs   (expect-eq !>(`@`0x4080.0004.c03f.fff5) !>((hd (pow:sad (ur 6 z6) (ur 6 w26)))))
++  test-cplx-exp-ch   (expect-eq !>(`@`0x40f1.bc86) !>((hd (exp:sad (ur 5 z5)))))
::  real-only / ordered functions are undefined on complex -> crash.
++  test-cplx-factorial-crashes  (expect-fail |.((factorial:sad (ur 6 z6))))
++  test-cplx-cbrt-crashes       (expect-fail |.((cbrt:sad (ur 6 z6))))
--
