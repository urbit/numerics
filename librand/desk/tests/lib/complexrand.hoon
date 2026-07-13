  ::  /tests/lib/complexrand
::::
::    /lib/complexrand: value regressions for both width doors (+cd, the
::    @rd/@cd reference precision, and +cs, the @rs/@cs mirror) against
::    known-good ship output, plus sanity checks that don't depend on any
::    particular RNG draw: +on-circle lands exactly on the unit circle
::    (re^2+im^2 = 1) and +in-disk lands strictly inside the unit disk
::    (re^2+im^2 < 1).
::
::    +cnormal:cs is also a regression for a real bug this file's own
::    construction caught: /lib/math's @rs +invsqt2 was missing a leading
::    `.0.` (parsed as the integer 70710677 instead of 0.70710677) -- see
::    libmath's tests/lib/math-constants.hoon for the fix itself.
::
/+  *test,
    rand,
    i754rand,
    complex,
    math,
    complexrand
|%
++  test-cuniform  ^-  tang
  ;:  weld
    %+  expect-eq
      !>(`@`0x3fe8.9e6a.a1b9.65f4.3f95.072f.63b9.b5e0)
      !>(out:(cuniform:cd:complexrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x3f39.65f4.3dee.6d78)
      !>(out:(cuniform:cs:complexrand (from-atom:seed:rand %sm64 0)))
  ==
::
++  test-normal-parts  ^-  tang
  ;:  weld
    %+  expect-eq
      !>(`@`0x3ff2.e308.2bd4.88f8.bfee.55f7.4af6.d8b9)
      !>(out:(normal-parts:cd:complexrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x4009.430e.bf17.e80f)
      !>(out:(normal-parts:cs:complexrand (from-atom:seed:rand %sm64 0)))
  ==
::
++  test-cnormal  ^-  tang
  ;:  weld
    %+  expect-eq
      !>(`@`0x3fea.b5c4.88e8.bf41.bfe5.735e.01a0.2b7b)
      !>(out:(cnormal:cd:complexrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x3fc2.1e20.bed6.d405)
      !>(out:(cnormal:cs:complexrand (from-atom:seed:rand %sm64 0)))
  ==
::
++  test-on-circle  ^-  tang
  ;:  weld
    %+  expect-eq
      !>(`@`0x3fc0.7838.f0db.1ef6.3fef.bbe7.a757.e0e1)
      !>(out:(on-circle:cd:complexrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x3f2b.0087.3f3e.82b8)
      !>(out:(on-circle:cs:complexrand (from-atom:seed:rand %sm64 0)))
  ==
::
++  test-in-disk  ^-  tang
  ;:  weld
    %+  expect-eq
      !>(`@`0xbfd1.1d5e.36cd.f850.bfe7.45ce.ffed.7562)
      !>(out:(in-disk:cd:complexrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x3ee5.97d0.bf44.64a2)
      !>(out:(in-disk:cs:complexrand (from-atom:seed:rand %sm64 0)))
  ==
::
::  +on-circle lands exactly on the unit circle regardless of seed: no
::  rejection loop, no accumulated error beyond the underlying cos/sin's
::  own correctly-rounded precision.  Checked across several distinct
::  seeds, not just the one regression value above.
++  on-circle-mag2
  |=  s=@
  ^-  @rd
  =/  m  ~(. rd:math [%n .~1e-13 .~0])
  =/  z  out:(on-circle:cd:complexrand (from-atom:seed:rand %sm64 s))
  =/  re  (~(re cd:complex %n) z)
  =/  im  (~(im cd:complex %n) z)
  (add:m (mul:m re re) (mul:m im im))
++  test-on-circle-unit-magnitude  ^-  tang
  =/  m  ~(. rd:math [%n .~1e-13 .~0])
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((gth:m (on-circle-mag2 0) .~0.999999999))
    %+  expect-eq  !>(%.y)  !>((lth:m (on-circle-mag2 0) .~1.000000001))
    %+  expect-eq  !>(%.y)  !>((gth:m (on-circle-mag2 1) .~0.999999999))
    %+  expect-eq  !>(%.y)  !>((lth:m (on-circle-mag2 1) .~1.000000001))
    %+  expect-eq  !>(%.y)  !>((gth:m (on-circle-mag2 2) .~0.999999999))
    %+  expect-eq  !>(%.y)  !>((lth:m (on-circle-mag2 2) .~1.000000001))
    %+  expect-eq  !>(%.y)  !>((gth:m (on-circle-mag2 3) .~0.999999999))
    %+  expect-eq  !>(%.y)  !>((lth:m (on-circle-mag2 3) .~1.000000001))
    %+  expect-eq  !>(%.y)  !>((gth:m (on-circle-mag2 4) .~0.999999999))
    %+  expect-eq  !>(%.y)  !>((lth:m (on-circle-mag2 4) .~1.000000001))
  ==
::
::  +in-disk always lands strictly inside the unit disk by construction
::  (the rejection loop's own accept condition), across several seeds.
++  in-disk-mag2
  |=  s=@
  ^-  @rd
  =/  m  ~(. rd:math [%n .~1e-13 .~0])
  =/  z  out:(in-disk:cd:complexrand (from-atom:seed:rand %sm64 s))
  =/  re  (~(re cd:complex %n) z)
  =/  im  (~(im cd:complex %n) z)
  (add:m (mul:m re re) (mul:m im im))
++  test-in-disk-inside-unit-disk  ^-  tang
  =/  m  ~(. rd:math [%n .~1e-13 .~0])
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((lth:m (in-disk-mag2 0) .~1))
    %+  expect-eq  !>(%.y)  !>((lth:m (in-disk-mag2 1) .~1))
    %+  expect-eq  !>(%.y)  !>((lth:m (in-disk-mag2 2) .~1))
    %+  expect-eq  !>(%.y)  !>((lth:m (in-disk-mag2 3) .~1))
    %+  expect-eq  !>(%.y)  !>((lth:m (in-disk-mag2 4) .~1))
  ==
--
