::::  /tests/lib/math-constants -- tau/pi/phi/sqt2/invsqt2 at all four
::::  precisions (@rs/@rd/@rh/@rq).
::
::  Regression for a real bug caught while building librand's /lib/complexrand
::  adapter: /lib/math's @rs +invsqt2 was written `.70710677` (missing the
::  leading `.0.` before the fractional digits), which Hoon parses as the
::  INTEGER 70,710,677.0, not 0.70710677 -- every other precision (rd/rh/rq)
::  had the correct `int.frac` form, so this only affected @rs. Nothing
::  exercised +invsqt2:rs:math until then. Expected bits below are ship-
::  verified (dojo `@ux` casts), not hand-derived.
::
/+  *test, math
|%
++  test-tau  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x40c9.0fdb)                                  !>(tau:rs:math)
    %+  expect-eq  !>(`@`0x4019.21fb.5444.2d18)                        !>(tau:rd:math)
    %+  expect-eq  !>(`@`0x4648)                                       !>(tau:rh:math)
    %+  expect-eq  !>(`@`0x4001.921f.b544.42d1.8469.898c.c517.01b8)    !>(tau:rq:math)
  ==
::
++  test-pi  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x4049.0fdb)                                  !>(pi:rs:math)
    %+  expect-eq  !>(`@`0x4009.21fb.5444.2d18)                        !>(pi:rd:math)
    %+  expect-eq  !>(`@`0x4248)                                       !>(pi:rh:math)
    %+  expect-eq  !>(`@`0x4000.921f.b544.42d1.8469.898c.c517.01b8)    !>(pi:rq:math)
  ==
::
++  test-phi  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3fcf.1bbd)                                  !>(phi:rs:math)
    %+  expect-eq  !>(`@`0x3ff9.e377.9b97.f4a8)                        !>(phi:rd:math)
    %+  expect-eq  !>(`@`0x3e79)                                       !>(phi:rh:math)
    %+  expect-eq  !>(`@`0x3fff.9e37.79b9.7f4a.7c15.f39c.c060.5cee)    !>(phi:rq:math)
  ==
::
++  test-sqt2  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3fb5.04f3)                                  !>(sqt2:rs:math)
    %+  expect-eq  !>(`@`0x3ff6.a09e.667f.3bcd)                        !>(sqt2:rd:math)
    %+  expect-eq  !>(`@`0x3da8)                                       !>(sqt2:rh:math)
    %+  expect-eq  !>(`@`0x3fff.6a09.e667.f3bc.c908.b2fb.1366.ea95)    !>(sqt2:rq:math)
  ==
::
::  +invsqt2:rs is the arm the bug above lived in -- checked against its
::  IEEE-754 value directly (0.7071067690849304 as a float32), not just
::  bit-equality with a value that could itself be wrong.
++  test-invsqt2  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0x3f35.04f3)                                  !>(invsqt2:rs:math)
    %+  expect-eq  !>(`@`0x3fe6.a09e.667f.3bcd)                        !>(invsqt2:rd:math)
    %+  expect-eq  !>(`@`0x39a8)                                       !>(invsqt2:rh:math)
    %+  expect-eq  !>(`@`0x3ffe.6a09.e667.f3bc.c908.b2fb.1366.ea95)    !>(invsqt2:rq:math)
  ==
--
