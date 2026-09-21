/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-builders -- array builders
::
::  Table-driven across precisions and representative shapes.  zeros/ones
::  are cross-checked against `fill`; eye/range/linspace against small
::  explicitly-built canons.
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
++  shapes  ^-  (list (list @))  ~[~[1 1] ~[2 3] ~[3 3]]
+$  rrow  [bloq=@ zero=@ one=@ two=@ three=@ five=@]
++  reals
  ^-  (list rrow)
  :~  [4 .~~0 .~~1 .~~2 .~~3 .~~5]
      [5 .0 .1 .2 .3 .5]
      [6 .~0 .~1 .~2 .~3 .~5]
      [7 .~~~0 .~~~1 .~~~2 .~~~3 .~~~5]
  ==
++  uint-bloqs  ^-  (list @)  ~[3 4 5 6]
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
::  iterate just the precisions (for builders whose shape is fixed)
++  each-prec
  |=  fun=$-(rrow tang)
  ^-  tang
  %-  zing
  (turn reals fun)
::
::  zeros == a fill with the zero value; ones == a fill with the one value.
++  test-zeros  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (zeros:la m))
++  test-zeros-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 0) (zeros:la m))
++  test-ones  ^-  tang
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (ones:la m))
++  test-ones-uint  ^-  tang
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 1) (ones:la m))
::  eye 2x2 == [[1 0] [0 1]] and eye 3x3 == identity, per precision.
++  test-eye-2x2  ^-  tang
  %-  each-prec
  |=  r=rrow
  =/  m  `meta`[~[2 2] bloq.r %i754 ~]
  (is-equal (en-ray:la [m ~[~[one.r zero.r] ~[zero.r one.r]]]) (eye:la m))
++  test-eye-3x3  ^-  tang
  %-  each-prec
  |=  r=rrow
  =/  m  `meta`[~[3 3] bloq.r %i754 ~]
  %+  is-equal
    (en-ray:la [m ~[~[one.r zero.r zero.r] ~[zero.r one.r zero.r] ~[zero.r zero.r one.r]]])
  (eye:la m)
::  range [0, 3) step 1 == [0 1 2], per precision.
++  test-range  ^-  tang
  %-  each-prec
  |=  r=rrow
  =/  m  `meta`[~[3] bloq.r %i754 ~]
  (is-equal (en-ray:la [m ~[zero.r one.r two.r]]) (range:la m [zero.r three.r] one.r))
::  linspace [0, 2] in 3 steps == [0 1 2], per precision.
++  test-linspace  ^-  tang
  %-  each-prec
  |=  r=rrow
  =/  m  `meta`[~[3] bloq.r %i754 ~]
  (is-equal (en-ray:la [m ~[zero.r one.r two.r]]) (linspace:la m [zero.r two.r] 3))
--
