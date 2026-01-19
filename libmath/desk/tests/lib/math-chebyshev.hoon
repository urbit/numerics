/+  *test
/+  *math-chebyshev
::
::::  Tests for Chebyshev/fixed-polynomial math library
::
^|
|_  $:  atol=_.1e-3          :: absolute tolerance
        rtol=_.1e-5          :: relative tolerance
    ==
::
::  Helper functions
::
++  expect-near-rs
  |=  [expected=@rs actual=@rs]  ^-  tang
  =/  diff  (abs:rs:^math-chebyshev (sub:^rs expected actual))
  ?:  (lth:^rs diff atol)
    ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<expected>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual" "{<actual>}"]]
      [%palm [": " ~ ~ ~] [leaf+"diff" "{<diff>}"]]
  ==
::
++  expect-near-rd
  |=  [expected=@rd actual=@rd]  ^-  tang
  =/  diff  `@rd`(dis (sub:^rd expected actual) 0x7fff.ffff.ffff.ffff)
  ?:  (lth:^rd diff .~1e-10)
    ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<expected>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual" "{<actual>}"]]
  ==
::
::  ============================================================
::  SINGLE PRECISION (@rs) TESTS
::  ============================================================
::
++  test-rs-sin  ^-  tang
  ;:  weld
    %+  expect-near-rs  .0  (sin:rs .0)
    %+  expect-near-rs  .0.84147096  (sin:rs .1)
    %+  expect-near-rs  .0.9092974  (sin:rs .2)
    %+  expect-near-rs  .-0.7568025  (sin:rs .4)
    %+  expect-near-rs  .0  (sin:rs .3.1415927)  ::  sin(pi) ~ 0
  ==
::
++  test-rs-cos  ^-  tang
  ;:  weld
    %+  expect-near-rs  .1  (cos:rs .0)
    %+  expect-near-rs  .0.5403023  (cos:rs .1)
    %+  expect-near-rs  .-0.4161468  (cos:rs .2)
    %+  expect-near-rs  .-0.6536436  (cos:rs .4)
    %+  expect-near-rs  .-1  (cos:rs .3.1415927)  ::  cos(pi) ~ -1
  ==
::
++  test-rs-tan  ^-  tang
  ;:  weld
    %+  expect-near-rs  .0  (tan:rs .0)
    %+  expect-near-rs  .1.5574077  (tan:rs .1)
  ==
::
++  test-rs-exp  ^-  tang
  ;:  weld
    %+  expect-near-rs  .1  (exp:rs .0)
    %+  expect-near-rs  .2.7182817  (exp:rs .1)
    %+  expect-near-rs  .7.389056  (exp:rs .2)
    %+  expect-near-rs  .0.36787945  (exp:rs .-1)
  ==
::
++  test-rs-log  ^-  tang
  ;:  weld
    %+  expect-near-rs  .0  (log:rs .1)
    %+  expect-near-rs  .0.6931472  (log:rs .2)
    %+  expect-near-rs  .2.302585  (log:rs .10)
    %+  expect-near-rs  .-0.6931472  (log:rs .0.5)
  ==
::
++  test-rs-sqrt  ^-  tang
  ;:  weld
    %+  expect-near-rs  .0  (sqt:rs .0)
    %+  expect-near-rs  .1  (sqt:rs .1)
    %+  expect-near-rs  .1.4142135  (sqt:rs .2)
    %+  expect-near-rs  .2  (sqt:rs .4)
    %+  expect-near-rs  .3  (sqt:rs .9)
    %+  expect-near-rs  .10  (sqt:rs .100)
  ==
::
++  test-rs-atan  ^-  tang
  ;:  weld
    %+  expect-near-rs  .0  (atan:rs .0)
    %+  expect-near-rs  .0.7853982  (atan:rs .1)        ::  atan(1) = pi/4
    %+  expect-near-rs  .1.1071488  (atan:rs .2)
    %+  expect-near-rs  .-0.7853982  (atan:rs .-1)
  ==
::
++  test-rs-asin  ^-  tang
  ;:  weld
    %+  expect-near-rs  .0  (asin:rs .0)
    %+  expect-near-rs  .0.5235988  (asin:rs .0.5)      ::  pi/6
    %+  expect-near-rs  .1.5707964  (asin:rs .1)        ::  pi/2
    %+  expect-near-rs  .-1.5707964  (asin:rs .-1)
  ==
::
++  test-rs-acos  ^-  tang
  ;:  weld
    %+  expect-near-rs  .1.5707964  (acos:rs .0)       ::  pi/2
    %+  expect-near-rs  .0  (acos:rs .1)
    %+  expect-near-rs  .3.1415927  (acos:rs .-1)      ::  pi
  ==
::
++  test-rs-pow-n  ^-  tang
  ;:  weld
    %+  expect-near-rs  .1  (pow-n:rs .2 .0)
    %+  expect-near-rs  .2  (pow-n:rs .2 .1)
    %+  expect-near-rs  .4  (pow-n:rs .2 .2)
    %+  expect-near-rs  .8  (pow-n:rs .2 .3)
    %+  expect-near-rs  .27  (pow-n:rs .3 .3)
  ==
::
++  test-rs-pow  ^-  tang
  ;:  weld
    %+  expect-near-rs  .1  (pow:rs .2 .0)
    %+  expect-near-rs  .8  (pow:rs .2 .3)
    %+  expect-near-rs  .2.828427  (pow:rs .2 .1.5)    ::  2^1.5 = 2*sqrt(2)
    %+  expect-near-rs  .0.5  (pow:rs .2 .-1)
  ==
::
++  test-rs-cheb  ^-  tang
  ::  Chebyshev polynomials: T_0=1, T_1=x, T_2=2x^2-1, T_3=4x^3-3x
  ;:  weld
    %+  expect-near-rs  .1  ((cheb:rs 0) .0.5)
    %+  expect-near-rs  .0.5  ((cheb:rs 1) .0.5)
    %+  expect-near-rs  .-0.5  ((cheb:rs 2) .0.5)      ::  2*(0.5)^2 - 1 = -0.5
    %+  expect-near-rs  .-1  ((cheb:rs 3) .0.5)        ::  4*(0.5)^3 - 3*0.5 = -1
  ==
::
::  ============================================================
::  DOUBLE PRECISION (@rd) TESTS
::  ============================================================
::
++  test-rd-sin  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~0  (sin:rd .~0)
    %+  expect-near-rd  .~0.8414709848078965  (sin:rd .~1)
    %+  expect-near-rd  .~0.9092974268256817  (sin:rd .~2)
  ==
::
++  test-rd-cos  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~1  (cos:rd .~0)
    %+  expect-near-rd  .~0.5403023058681398  (cos:rd .~1)
  ==
::
++  test-rd-exp  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~1  (exp:rd .~0)
    %+  expect-near-rd  .~2.718281828459045  (exp:rd .~1)
  ==
::
++  test-rd-log  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~0  (log:rd .~1)
    %+  expect-near-rd  .~0.6931471805599453  (log:rd .~2)
    %+  expect-near-rd  .~2.302585092994046  (log:rd .~10)
  ==
::
++  test-rd-sqrt  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~1  (sqt:rd .~1)
    %+  expect-near-rd  .~1.4142135623730951  (sqt:rd .~2)
    %+  expect-near-rd  .~2  (sqt:rd .~4)
  ==
::
++  test-rd-atan  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~0  (atan:rd .~0)
    %+  expect-near-rd  .~0.7853981633974483  (atan:rd .~1)  ::  pi/4
  ==
::
++  test-rd-pow-n  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~8  (pow-n:rd .~2 .~3)
    %+  expect-near-rd  .~27  (pow-n:rd .~3 .~3)
  ==
::
++  test-rd-cheb  ^-  tang
  ;:  weld
    %+  expect-near-rd  .~1  ((cheb:rd 0) .~0.5)
    %+  expect-near-rd  .~0.5  ((cheb:rd 1) .~0.5)
    %+  expect-near-rd  .~-0.5  ((cheb:rd 2) .~0.5)
    %+  expect-near-rd  .~-1  ((cheb:rd 3) .~0.5)
  ==
::
::  ============================================================
::  IDENTITY TESTS - verify mathematical relationships
::  ============================================================
::
++  test-rs-sin-cos-identity  ^-  tang
  ::  sin^2(x) + cos^2(x) = 1
  =/  x  .1.234
  =/  s  (sin:rs x)
  =/  c  (cos:rs x)
  =/  sum  (add:^rs (mul:^rs s s) (mul:^rs c c))
  (expect-near-rs .1 sum)
::
++  test-rs-exp-log-identity  ^-  tang
  ::  exp(log(x)) = x for x > 0
  =/  x  .2.5
  =/  result  (exp:rs (log:rs x))
  (expect-near-rs x result)
::
++  test-rs-log-exp-identity  ^-  tang
  ::  log(exp(x)) = x
  =/  x  .1.5
  =/  result  (log:rs (exp:rs x))
  (expect-near-rs x result)
::
++  test-rs-sqrt-square-identity  ^-  tang
  ::  sqrt(x)^2 = x for x >= 0
  =/  x  .7.5
  =/  s  (sqt:rs x)
  =/  result  (mul:^rs s s)
  (expect-near-rs x result)
::
++  test-rs-atan-tan-identity  ^-  tang
  ::  atan(tan(x)) = x for |x| < pi/2
  =/  x  .0.5
  =/  result  (atan:rs (tan:rs x))
  (expect-near-rs x result)
--
