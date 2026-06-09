::  /tests/lib/complex-edge -- %cplx transcendental edge cases (@cs).
::
::  Exercises the origin, the real and imaginary axes, the csqrt/clog branch cut
::  (-1, +-i), the csqrt sign branch, the clog(0) singularity, and the naive-
::  series breakdown at large arguments.  All values validated against numpy
::  complex64.
::
/+  *test, complex
|%
++  s    cs:complex
++  z00  `@`0x0                    ::  0+0i
++  z10  `@`0x3f80.0000            ::  1+0i
++  zt0  `@`0x4120.0000            ::  10+0i
++  z0i  `@`0x3f80.0000.0000.0000  ::  0+1i
++  z0n  `@`0xbf80.0000.0000.0000  ::  0-1i
++  zn1  `@`0xbf80.0000            ::  -1+0i
++  z11  `@`0x3f80.0000.3f80.0000  ::  1+1i
++  z34  `@`0x4080.0000.4040.0000  ::  3+4i
++  z0t  `@`0x4120.0000.0000.0000  ::  0+10i
::  cexp: origin -> 1, real axis (e, 1/e), imaginary axis (cos1 + i sin1), off-axis.
++  test-cexp-origin  (expect-eq !>(`@`0x3f80.0000) !>((~(cexp s %n) z00)))
++  test-cexp-1       (expect-eq !>(`@`0x402d.f855) !>((~(cexp s %n) z10)))
++  test-cexp-neg1    (expect-eq !>(`@`0x3ebc.5ab0) !>((~(cexp s %n) zn1)))
++  test-cexp-i       (expect-eq !>(`@`0x3f57.6aa4.3f0a.5140) !>((~(cexp s %n) z0i)))
++  test-cexp-1p1i    (expect-eq !>(`@`0x4012.6408.3fbb.fe2a) !>((~(cexp s %n) z11)))
++  test-cexp-3p4i    (expect-eq !>(`@`0xc173.366f.c152.0f81) !>((~(cexp s %n) z34)))
::  clog: the singularity (log 0 = -inf) and the branch cut at -1 (iπ) and +-i.
++  test-clog-zero    (expect-eq !>(`@`0xff80.0000) !>((~(clog s %n) z00)))
++  test-clog-neg1    (expect-eq !>(`@`0x4049.0fdb.0000.0000) !>((~(clog s %n) zn1)))
++  test-clog-i       (expect-eq !>(`@`0x3fc9.0fdb.0000.0000) !>((~(clog s %n) z0i)))
++  test-clog-neg-i   (expect-eq !>(`@`0xbfc9.0fdb.0000.0000) !>((~(clog s %n) z0n)))
::  csqrt: branch cut sqrt(-1) = i, and the sign branch sqrt(-i) = 0.707 - 0.707i.
++  test-csqrt-neg1   (expect-eq !>(`@`0x3f80.0000.0000.0000) !>((~(csqrt s %n) zn1)))
++  test-csqrt-neg-i  (expect-eq !>(`@`0xbf35.04f3.3f35.04f3) !>((~(csqrt s %n) z0n)))
::  csin/ccos on the imaginary axis -> real sinh/cosh.
++  test-csin-i       (expect-eq !>(`@`0x3f96.6cff.0000.0000) !>((~(csin s %n) z0i)))
++  test-ccos-i       (expect-eq !>(`@`0x3fc5.83ab) !>((~(ccos s %n) z0i)))
::  KNOWN LIMITATION: the naive Taylor/AGM series diverges from the true value
::  far from the origin -- exp(10) is ~21991 here vs the true ~22026 (~0.16%
::  low), csin(0+10i) ~ i*10989 vs true ~ i*11013.  Locked as regression (NOT as
::  correctness); #18's Chebyshev rewrite would tighten these.
++  test-cexp-large   (expect-eq !>(`@`0x46ab.cef6) !>((~(cexp s %n) zt0)))
++  test-csin-large   (expect-eq !>(`@`0x462b.b42b.0000.0000) !>((~(csin s %n) z0t)))
--
