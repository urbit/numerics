/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-linalg -- linear algebra
::
::  Named ++test-<op>-<shape>-<kind>.  A `canon` is the reference result,
::  an `assay` is the result of the operation under test.
::
^|
|_  $:  atol=_.1e-3          :: absolute tolerance for precision of operations
        rtol=_.1e-5          :: relative tolerance for precision of operations
    ==
::  Auxiliary tools
++  is-equal
  |=  [a=ray b=ray]  ^-  tang
  ?:  =(a b)  ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<`ray`a>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual  " "{<`ray`b>}"]]
  ==
::
++  is-close
  |=  [a=ray b=ray]  ^-  tang
  ?:  (all:la (is-close:la a b [atol rtol]))  ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<a>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual  " "{<b>}"]]
  ==
::

++  test-mmul-3x3-3u  ^-  tang
  =/  meta-3x3-3  [~[3 3] 3 %uint ~]
  =/  assay-3x3-3-i  (eye:la meta-3x3-3)
  =/  assay-3x3-3  (en-ray:la [meta-3x3-3 ~[~[1 2 3] ~[4 5 6] ~[7 8 9]]])
  =/  canon-3x3-3  (en-ray:la [meta-3x3-3 ~[~[30 36 42] ~[66 81 96] ~[102 126 150]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x3-3)
      !>((mmul:la assay-3x3-3 assay-3x3-3))
    %+  expect-eq
      !>(assay-3x3-3)
      !>((mmul:la assay-3x3-3 assay-3x3-3-i))
    %+  expect-eq
      !>((en-ray:la [[~[3 1] 3 %uint ~] ~[~[6] ~[15] ~[24]]]))
      !>((mmul:la assay-3x3-3 (ones:la [~[3 1] 3 %uint ~])))
  ==
::

++  test-mmul-3x3-4u  ^-  tang
  =/  meta-3x3-4  [~[3 3] 4 %uint ~]
  =/  assay-3x3-4-i  (eye:la meta-3x3-4)
  =/  assay-3x3-4  (en-ray:la [meta-3x3-4 ~[~[1 2 3] ~[4 5 6] ~[7 8 9]]])
  =/  canon-3x3-4  (en-ray:la [meta-3x3-4 ~[~[30 36 42] ~[66 81 96] ~[102 126 150]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x3-4)
      !>((mmul:la assay-3x3-4 assay-3x3-4))
    %+  expect-eq
      !>(assay-3x3-4)
      !>((mmul:la assay-3x3-4 assay-3x3-4-i))
    %+  expect-eq
      !>((en-ray:la [[~[3 1] 4 %uint ~] ~[~[6] ~[15] ~[24]]]))
      !>((mmul:la assay-3x3-4 (ones:la [~[3 1] 4 %uint ~])))
  ==
::

++  test-mmul-3x3-5u  ^-  tang
  =/  meta-3x3-5  [~[3 3] 5 %uint ~]
  =/  assay-3x3-5-i  (eye:la meta-3x3-5)
  =/  assay-3x3-5  (en-ray:la [meta-3x3-5 ~[~[1 2 3] ~[4 5 6] ~[7 8 9]]])
  =/  canon-3x3-5  (en-ray:la [meta-3x3-5 ~[~[30 36 42] ~[66 81 96] ~[102 126 150]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x3-5)
      !>((mmul:la assay-3x3-5 assay-3x3-5))
    %+  expect-eq
      !>(assay-3x3-5)
      !>((mmul:la assay-3x3-5 assay-3x3-5-i))
    %+  expect-eq
      !>((en-ray:la [[~[3 1] 5 %uint ~] ~[~[6] ~[15] ~[24]]]))
      !>((mmul:la assay-3x3-5 (ones:la [~[3 1] 5 %uint ~])))
  ==
::

++  test-mmul-3x4x5-4r  ^-  tang
  =/  meta-3x4-4  [~[3 4] 4 %i754 ~]
  =/  meta-4x5-4  [~[4 5] 4 %i754 ~]
  =/  meta-3x5-4  [~[3 5] 4 %i754 ~]
  =/  assay-3x4-4  (en-ray:la [meta-3x4-4 ~[~[.~~1 .~~2 .~~3 .~~4] ~[.~~5 .~~6 .~~7 .~~8] ~[.~~9 .~~10 .~~11 .~~12]]])
  =/  assay-4x5-4  (en-ray:la [meta-4x5-4 ~[~[.~~1 .~~2 .~~3 .~~4 .~~5] ~[.~~4 .~~5 .~~6 .~~7 .~~8] ~[.~~7 .~~8 .~~9 .~~10 .~~11] ~[.~~10 .~~11 .~~12 .~~13 .~~14]]])
  =/  canon-3x5-4  (en-ray:la [meta-3x5-4 ~[~[.~~70 .~~80 .~~90 .~~100 .~~110] ~[.~~158 .~~184 .~~210 .~~236 .~~262] ~[.~~246 .~~288 .~~330 .~~372 .~~414]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x5-4)
      !>((mmul:la assay-3x4-4 assay-4x5-4))
  ==
::

++  test-mmul-3x4x5-5r  ^-  tang
  =/  meta-3x4-5  [~[3 4] 5 %i754 ~]
  =/  meta-4x5-5  [~[4 5] 5 %i754 ~]
  =/  meta-3x5-5  [~[3 5] 5 %i754 ~]
  =/  assay-3x4-5  (en-ray:la [meta-3x4-5 ~[~[.1 .2 .3 .4] ~[.5 .6 .7 .8] ~[.9 .10 .11 .12]]])
  =/  assay-4x5-5  (en-ray:la [meta-4x5-5 ~[~[.1 .2 .3 .4 .5] ~[.4 .5 .6 .7 .8] ~[.7 .8 .9 .10 .11] ~[.10 .11 .12 .13 .14]]])
  =/  canon-3x5-5  (en-ray:la [meta-3x5-5 ~[~[.70 .80 .90 .100 .110] ~[.158 .184 .210 .236 .262] ~[.246 .288 .330 .372 .414]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x5-5)
      !>((mmul:la assay-3x4-5 assay-4x5-5))
  ==
::

++  test-mmul-3x4x5-6r  ^-  tang
  =/  meta-3x4-6  [~[3 4] 6 %i754 ~]
  =/  meta-4x5-6  [~[4 5] 6 %i754 ~]
  =/  meta-3x5-6  [~[3 5] 6 %i754 ~]
  =/  assay-3x4-6  (en-ray:la [meta-3x4-6 ~[~[.~1 .~2 .~3 .~4] ~[.~5 .~6 .~7 .~8] ~[.~9 .~10 .~11 .~12]]])
  =/  assay-4x5-6  (en-ray:la [meta-4x5-6 ~[~[.~1 .~2 .~3 .~4 .~5] ~[.~4 .~5 .~6 .~7 .~8] ~[.~7 .~8 .~9 .~10 .~11] ~[.~10 .~11 .~12 .~13 .~14]]])
  =/  canon-3x5-6  (en-ray:la [meta-3x5-6 ~[~[.~70 .~80 .~90 .~100 .~110] ~[.~158 .~184 .~210 .~236 .~262] ~[.~246 .~288 .~330 .~372 .~414]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x5-6)
      !>((mmul:la assay-3x4-6 assay-4x5-6))
  ==
::

++  test-mmul-3x4x5-7r  ^-  tang
  =/  meta-3x4-7  [~[3 4] 7 %i754 ~]
  =/  meta-4x5-7  [~[4 5] 7 %i754 ~]
  =/  meta-3x5-7  [~[3 5] 7 %i754 ~]
  =/  assay-3x4-7  (en-ray:la [meta-3x4-7 ~[~[.~~~1 .~~~2 .~~~3 .~~~4] ~[.~~~5 .~~~6 .~~~7 .~~~8] ~[.~~~9 .~~~10 .~~~11 .~~~12]]])
  =/  assay-4x5-7  (en-ray:la [meta-4x5-7 ~[~[.~~~1 .~~~2 .~~~3 .~~~4 .~~~5] ~[.~~~4 .~~~5 .~~~6 .~~~7 .~~~8] ~[.~~~7 .~~~8 .~~~9 .~~~10 .~~~11] ~[.~~~10 .~~~11 .~~~12 .~~~13 .~~~14]]])
  =/  canon-3x5-7  (en-ray:la [meta-3x5-7 ~[~[.~~~70 .~~~80 .~~~90 .~~~100 .~~~110] ~[.~~~158 .~~~184 .~~~210 .~~~236 .~~~262] ~[.~~~246 .~~~288 .~~~330 .~~~372 .~~~414]]])
  ;:  weld
    %+  expect-eq
      !>(canon-3x5-7)
      !>((mmul:la assay-3x4-7 assay-4x5-7))
  ==
::

++  test-dot-1-4r  ^-  tang
  =/  meta-1x1-4  [~[1 1] 4 %i754 ~]
  =/  assay-1x1-4  (en-ray:la [meta-1x1-4 ~[~[.~~10]]])
  =/  canon-1x1-4  (en-ray:la [meta-1x1-4 ~[~[.~~100]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-4)
      !>((dot:la assay-1x1-4 assay-1x1-4))
  ==
::

++  test-dot-4-4r  ^-  tang
  =/  meta-1x1-4  [~[1 1] 4 %i754 ~]
  =/  meta-1x4-4  [~[1 4] 4 %i754 ~]
  =/  assay-1x4-a-4  (en-ray:la [meta-1x4-4 ~[~[.~~1 .~~2 .~~3 .~~4]]])
  =/  assay-1x4-b-4  (en-ray:la [meta-1x4-4 ~[~[.~~0.5 .~~0.25 .~~0.125 .~~0.0625]]])
  =/  canon-1x1-4  (en-ray:la [meta-1x1-4 ~[~[.~~1.625]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-4)
      !>((dot:la assay-1x4-a-4 assay-1x4-b-4))
  ==
::

++  test-dot-1-5r  ^-  tang
  =/  meta-1x1-5  [~[1 1] 5 %i754 ~]
  =/  assay-1x1-5  (en-ray:la [meta-1x1-5 ~[~[.10]]])
  =/  canon-1x1-5  (en-ray:la [meta-1x1-5 ~[~[.100]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-5)
      !>((dot:la assay-1x1-5 assay-1x1-5))
  ==
::

++  test-dot-4-5r  ^-  tang
  =/  meta-1x1-5  [~[1 1] 5 %i754 ~]
  =/  meta-1x4-5  [~[1 4] 5 %i754 ~]
  =/  assay-1x4-a-5  (en-ray:la [meta-1x4-5 ~[~[.1 .2 .3 .4 .5]]])
  =/  assay-1x4-b-5  (en-ray:la [meta-1x4-5 ~[~[.0.5 .0.25 .0.125 .0.0625 .0.03125]]])
  =/  canon-1x1-5  (en-ray:la [meta-1x1-5 ~[~[.1.625]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-5)
      !>((dot:la assay-1x4-a-5 assay-1x4-b-5))
  ==
::

++  test-dot-1-6r  ^-  tang
  =/  meta-1x1-6  [~[1 1] 6 %i754 ~]
  =/  assay-1x1-6  (en-ray:la [meta-1x1-6 ~[~[.~10]]])
  =/  canon-1x1-6  (en-ray:la [meta-1x1-6 ~[~[.~100]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-6)
      !>((dot:la assay-1x1-6 assay-1x1-6))
  ==
::

++  test-dot-4-6r  ^-  tang
  =/  meta-1x1-6  [~[1 1] 6 %i754 ~]
  =/  meta-1x4-6  [~[1 4] 6 %i754 ~]
  =/  assay-1x4-a-6  (en-ray:la [meta-1x4-6 ~[~[.~1 .~2 .~3 .~4 .~5 .~6]]])
  =/  assay-1x4-b-6  (en-ray:la [meta-1x4-6 ~[~[.~0.5 .~0.25 .~0.125 .~0.0625 .~0.03125 .~0.015625]]])
  =/  canon-1x1-6  (en-ray:la [meta-1x1-6 ~[~[.~1.625]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-6)
      !>((dot:la assay-1x4-a-6 assay-1x4-b-6))
  ==
::

++  test-dot-1-7r  ^-  tang
  =/  meta-1x1-7  [~[1 1] 7 %i754 ~]
  =/  assay-1x1-7  (en-ray:la [meta-1x1-7 ~[~[.~~~10]]])
  =/  canon-1x1-7  (en-ray:la [meta-1x1-7 ~[~[.~~~100]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-7)
      !>((dot:la assay-1x1-7 assay-1x1-7))
  ==
::

++  test-dot-4-7r  ^-  tang
  =/  meta-1x1-7  [~[1 1] 7 %i754 ~]
  =/  meta-1x4-7  [~[1 4] 7 %i754 ~]
  =/  assay-1x4-a-7  (en-ray:la [meta-1x4-7 ~[~[.~~~1 .~~~2 .~~~3 .~~~4 .~~~5 .~~~6 .~~~7]]])
  =/  assay-1x4-b-7  (en-ray:la [meta-1x4-7 ~[~[.~~~0.5 .~~~0.25 .~~~0.125 .~~~0.0625 .~~~0.03125 .~~~0.015625 .~~~0.0078125]]])
  =/  canon-1x1-7  (en-ray:la [meta-1x1-7 ~[~[.~~~1.625]]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x1-7)
      !>((dot:la assay-1x4-a-7 assay-1x4-b-7))
  ==
--
