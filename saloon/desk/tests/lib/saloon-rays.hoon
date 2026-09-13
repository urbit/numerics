::  Tests for Saloon ray-level transcendental operations (+sa core).
::
::  The existing saloon.hoon test file covers *scalar* math via +rs:saloon
::  (which resolves to the scalar math door in /lib/math).  This file covers
::  the *ray-level* operations -- ++exp, ++sin, ++cos, ++tan, ++log,
::  ++log-2, ++log-10, ++sqrt, ++cbrt, ++pow-n, ++pow -- verifying that
::  ++trans-scalar dispatches correctly through ++el-wise-op for %i754 rays.
::
::  Expected values match the scalar tests in saloon.hoon so the two suites
::  cross-validate: if a scalar test and its ray counterpart disagree, the
::  dispatch is broken.
::
/-  ls=lagoon
/+  *test, *saloon, *lagoon
|%
::
::  ++close-rs: are two @rs (bloq=5) rays elementwise within ~0.1%?
::  Uses is-close with atol=rtol=0x3a83126f ≈ 0.001 (IEEE @rs).
::
++  close-rs
  |=  [x=ray:ls y=ray:ls]  ^-  ?
  (all:(lake %n) (is-close:(lake %n) x y [0x3a83.126f 0x3a83.126f]))
::
::  ++rs1: single-element @rs (bloq=5) ray.
::  ++rs2: two-element @rs ray.
::
++  rs1
  |=  [a=@]  ^-  ray:ls
  (en-ray:(lake %n) [[~[1] 5 %i754 ~] ~[a]])
++  rs2
  |=  [a=@ b=@]  ^-  ray:ls
  (en-ray:(lake %n) [[~[2] 5 %i754 ~] ~[a b]])
::
::  ++close-rh/rd/rq + ++rh1/rd1/rq1: the same pattern as ++close-rs/++rs1
::  above, at the other three %i754 bloqs (4=@rh, 6=@rd, 7=@rq).  Added
::  alongside the saloon.hoon +sa scalar-dispatch rtol fix (rand-spec.md
::  milestone 8 prerequisite): the existing tests here only ever exercised
::  bloq 5 (@rs), so the fix's bloq 4/6/7 branches -- mechanically
::  identical, but otherwise untested at runtime -- get one spot check
::  each below (+test-exp-rh/rd/rq).
::
++  close-rh
  |=  [x=ray:ls y=ray:ls]  ^-  ?
  (all:(lake %n) (is-close:(lake %n) x y [.~~1e-2 .~~1e-2]))
++  close-rd
  |=  [x=ray:ls y=ray:ls]  ^-  ?
  (all:(lake %n) (is-close:(lake %n) x y [.~1e-9 .~1e-9]))
++  close-rq
  |=  [x=ray:ls y=ray:ls]  ^-  ?
  (all:(lake %n) (is-close:(lake %n) x y [.~~~1e-15 .~~~1e-15]))
::
++  rh1
  |=  [a=@]  ^-  ray:ls
  (en-ray:(lake %n) [[~[1] 4 %i754 ~] ~[a]])
++  rd1
  |=  [a=@]  ^-  ray:ls
  (en-ray:(lake %n) [[~[1] 6 %i754 ~] ~[a]])
++  rq1
  |=  [a=@]  ^-  ray:ls
  (en-ray:(lake %n) [[~[1] 7 %i754 ~] ~[a]])
::
::  Unary ops: input and expected-output rays share the same values as the
::  scalar tests in saloon.hoon so the two suites cross-validate.
::
++  test-exp-rs
  ^-  tang
  =/  a  (rs2 .4 .5)
  =/  w  (rs2 .54.598 .148.413)
  (expect !>((close-rs (exp:sa a) w)))
::
++  test-sin-rs
  ^-  tang
  =/  a  (rs2 .4 .0.55)
  =/  w  (rs2 .-0.756802 .0.522687)
  (expect !>((close-rs (sin:sa a) w)))
::
++  test-cos-rs
  ^-  tang
  =/  a  (rs2 .4 .0.55)
  =/  w  (rs2 .-0.65364 .0.852525)
  (expect !>((close-rs (cos:sa a) w)))
::
++  test-tan-rs
  ^-  tang
  =/  a  (rs2 .4 .0.55)
  =/  w  (rs2 .1.15782 .0.613105)
  (expect !>((close-rs (tan:sa a) w)))
::
++  test-log-rs
  ^-  tang
  =/  a  (rs2 .0.1 .60)
  =/  w  (rs2 .-2.30259 .4.094345)
  (expect !>((close-rs (log:sa a) w)))
::
++  test-log-2-rs
  ^-  tang
  =/  a  (rs2 .30 .0.55)
  =/  w  (rs2 .4.9069 .-0.8625)
  (expect !>((close-rs (log-2:sa a) w)))
::
++  test-log-10-rs
  ^-  tang
  =/  a  (rs2 .30 .0.66)
  =/  w  (rs2 .1.477121 .-0.180456)
  (expect !>((close-rs (log-10:sa a) w)))
::
++  test-sqrt-rs
  ^-  tang
  =/  a  (rs2 .4 .2)
  =/  w  (rs2 .2 .1.41421)
  (expect !>((close-rs (sqrt:sa a) w)))
::
++  test-cbrt-rs
  ^-  tang
  =/  a  (rs2 .27 .4.4)
  =/  w  (rs2 .3 .1.63864)
  (expect !>((close-rs (cbrt:sa a) w)))
::
::  Binary ops: pow-n and pow take matching-shape rays.
::
++  test-pow-n-rs
  ^-  tang
  =/  a  (rs2 .5.1 .-3)
  =/  b  (rs2 .3 .3)
  =/  w  (rs2 .132.651 .-27)
  (expect !>((close-rs (pow-n:sa a b) w)))
::
++  test-pow-rs
  ^-  tang
  =/  a  (rs2 .5 .3.3)
  =/  b  (rs2 .3.3 .5)
  =/  w  (rs2 .202.582 .391.35393)
  (expect !>((close-rs (pow:sa a b) w)))
::
::  Shape-preservation: a [3 2] matrix stays [3 2] after element-wise op.
::
++  test-shape-preserved-rs
  ^-  tang
  =/  a  (en-ray:(lake %n) [[~[3 2] 5 %i754 ~] ~[.1 .2 .3 .4 .5 .6]])
  =/  res  (exp:sa a)
  (expect !>(=(shape.meta.res ~[3 2])))
::
::  +test-exp-rh/rd/rq: see the ++close-rh/rd/rq header comment above.
::
++  test-exp-rh
  ^-  tang
  =/  a  (rh1 .~~4)
  =/  w  (rh1 .~~54.6)
  (expect !>((close-rh (exp:sa a) w)))
::
++  test-exp-rd
  ^-  tang
  =/  a  (rd1 .~4)
  =/  w  (rd1 .~54.598150033144236)
  (expect !>((close-rd (exp:sa a) w)))
::
++  test-exp-rq
  ^-  tang
  =/  a  (rq1 .~~~4)
  =/  w  (rq1 .~~~54.59815003314423907811026120286088)
  (expect !>((close-rq (exp:sa a) w)))
--
