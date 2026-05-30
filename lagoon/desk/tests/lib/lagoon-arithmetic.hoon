/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-arithmetic -- elementwise arithmetic
::
::  Table-driven: each op is asserted once across all precisions and a few
::  representative shapes (scalar, rectangular, square) via `fill`-built
::  rays, rather than hand-written per (op, bloq, shape, kind).
::
^|
|_  $:  atol=_.1e-3          :: absolute tolerance for precision of operations
        rtol=_.1e-5          :: relative tolerance for precision of operations
    ==
++  is-equal
  |=  [a=ray b=ray]  ^-  tang
  ?:  =(a b)  ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<`ray`a>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual  " "{<`ray`b>}"]]
  ==
::
::  Representative shapes: a 1x1 scalar, a 2x3 rectangle, a 3x3 square.
++  shapes  ^-  (list (list @))  ~[~[1 1] ~[2 3] ~[3 3]]
::  Per-precision sample values, written once here instead of inline in
::  every arm.  [bloq zero one two three five].
+$  rrow  [bloq=@ zero=@ one=@ two=@ three=@ five=@]
++  reals
  ^-  (list rrow)
  :~  [4 .~~0 .~~1 .~~2 .~~3 .~~5]
      [5 .0 .1 .2 .3 .5]
      [6 .~0 .~1 .~2 .~3 .~5]
      [7 .~~~0 .~~~1 .~~~2 .~~~3 .~~~5]
  ==
++  uint-bloqs  ^-  (list @)  ~[3 4 5 6]
::  Run an %i754 assertion over reals x shapes.
++  each-real
  |=  fun=$-([rrow meta] tang)
  ^-  tang
  %-  zing
  %+  turn  reals
  |=  r=rrow
  ^-  tang
  %-  zing
  %+  turn  shapes
  |=  shp=(list @)
  (fun r [shp bloq.r %i754 ~])
::  Run a %uint assertion over bloqs x shapes.
++  each-uint
  |=  fun=$-(meta tang)
  ^-  tang
  %-  zing
  %+  turn  uint-bloqs
  |=  b=@
  ^-  tang
  %-  zing
  %+  turn  shapes
  |=  shp=(list @)
  (fun [shp b %uint ~])
::
::  ones + ones = twos
++  test-add  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m two.r) (add:la (fill:la m one.r) (fill:la m one.r)))
++  test-add-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 2) (add:la (fill:la m 1) (fill:la m 1)))
::  ones - ones = zeros
++  test-sub  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (sub:la (fill:la m one.r) (fill:la m one.r)))
++  test-sub-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 0) (sub:la (fill:la m 1) (fill:la m 1)))
::  fives - threes = twos  (asymmetric: distinguishes x-y from y-x and from
::  returning an unmodified operand -- the bug fixed in the sub jet)
++  test-sub-asym  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m two.r) (sub:la (fill:la m five.r) (fill:la m three.r)))
::  ones * ones = ones
++  test-mul  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (mul:la (fill:la m one.r) (fill:la m one.r)))
++  test-mul-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 1) (mul:la (fill:la m 1) (fill:la m 1)))
::  ones / ones = ones
++  test-div  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (div:la (fill:la m one.r) (fill:la m one.r)))
++  test-div-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 1) (div:la (fill:la m 1) (fill:la m 1)))
::  ones mod ones = zeros
++  test-mod  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (mod:la (fill:la m one.r) (fill:la m one.r)))
++  test-mod-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 0) (mod:la (fill:la m 1) (fill:la m 1)))
::
::  Regression: mod must round the quotient with the active mode (default
::  %n, round-nearest), not truncate.  5 mod 3 = 5 - 3*round(5/3) = -1.
++  test-mod-nearest-6r  ^-  tang
  =/  m  `meta`[~[1 1] 6 %i754 ~]
  (is-equal (fill:la m .~-1) (mod:la (fill:la m .~5) (fill:la m .~3)))
::  Regression: quad mod-scalar reciprocal constant (was 0.0, returned x).
::  7 mod 3 = 1, 5 mod 3 = 2.
++  test-mods-7r  ^-  tang
  =/  input-mods-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~7.0 .~~~5.0]]])
  =/  canon-mods-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~2.0]]])
  %+  is-equal
    canon-mods-1x2-7r
  (mod-scalar:la input-mods-1x2-7r .~~~3.0)
--
