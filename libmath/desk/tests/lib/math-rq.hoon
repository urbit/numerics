::::  /tests/lib/math-rq -- @rq (f128) transcendentals (PR #18).
::  Faithful native f128 (Cody-Waite/fdlibm at degree ~24).  Expected bits are
::  the SoftFloat-f128 output of libmath/tools/rq_check.c (the jet basis),
::  cross-checked against MPFR truth.
::
/+  *test, math
|%
++  eq  |=(x=@rq ^-(@ `@`(~(exp rq:math [%n .~~~1e-10]) x)))
::  ==== exp ====
++  test-exp-half  %+  expect-eq  !>(`@`0x3fff.a612.98e1.e069.bc97.2dfe.fab6.df34)
  !>((eq `@rq`0x3ffe.0000.0000.0000.0000.0000.0000.0000))
++  test-exp-1     %+  expect-eq  !>(`@`0x4000.5bf0.a8b1.4576.9535.5fb8.ac40.4e7a)
  !>((eq `@rq`0x3fff.0000.0000.0000.0000.0000.0000.0000))
++  test-exp-2     %+  expect-eq  !>(`@`0x4001.d8e6.4b8d.4dda.dcc3.3a3b.a206.b68b)
  !>((eq `@rq`0x4000.0000.0000.0000.0000.0000.0000.0000))
++  test-exp-n2    %+  expect-eq  !>(`@`0x3ffc.152a.aa3b.f81c.b9fd.b76e.ae12.d029)
  !>((eq `@rq`0xc000.0000.0000.0000.0000.0000.0000.0000))
++  test-exp-10    %+  expect-eq  !>(`@`0x400d.5829.dcf9.5055.f9f0.7ea8.c056.d135)
  !>((eq `@rq`0x4002.4000.0000.0000.0000.0000.0000.0000))
++  test-exp-0     %+  expect-eq  !>(`@`0x3fff.0000.0000.0000.0000.0000.0000.0000)
  !>((eq `@rq`0x0))
++  test-exp-inf   %+  expect-eq  !>(`@`0x7fff.0000.0000.0000.0000.0000.0000.0000)
  !>((eq `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000))
++  test-exp-ninf  %+  expect-eq  !>(`@`0x0)
  !>((eq `@rq`0xffff.0000.0000.0000.0000.0000.0000.0000))
++  test-exp-nan   %+  expect-eq  !>(`@`0x7fff.8000.0000.0000.0000.0000.0000.0000)
  !>((eq `@rq`0x7fff.8000.0000.0000.0000.0000.0000.0000))
--
