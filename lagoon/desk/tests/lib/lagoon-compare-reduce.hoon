/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-compare-reduce -- comparisons and reductions
::
::  Table-driven.  Comparisons run over precisions x shapes and cover all
::  three branches (>, <, =) per op rather than only the symmetric ones-vs-
::  ones case the old suite used.  Reductions run over precisions on a
::  distinct-valued vector and read the reduced scalar via its ravel head
::  (independent of min/max's output shape).
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
::  elementwise: %i754 over reals x shapes
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
::  elementwise: %uint over bloqs x shapes
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
::  reductions: one distinct-valued 1-D vector [zero one two] per precision
++  each-vec
  |=  fun=$-([rrow ray] tang)
  ^-  tang
  %-  zing
  %+  turn  reals
  |=  r=rrow
  (fun r (en-ray:la [`meta`[~[3] bloq.r %i754 ~] ~[zero.r one.r two.r]]))
::  head scalar of a reduced ray, rank-independent
++  hed  |=(a=ray ^-(@ -:(ravel:la a)))
::
::  Comparisons (numeric convention true=1, false=0).  Both branches of
::  each op are covered -- one true-producing and one false-producing
::  arm -- across all precisions and shapes.  (The old suite only tested
::  the symmetric ones-vs-ones case, exercising one branch per op.)
++  test-gth-t  ^-  tang  :: 2 > 1
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (gth:la (fill:la m two.r) (fill:la m one.r)))
++  test-gth-f  ^-  tang  :: 1 > 1
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (gth:la (fill:la m one.r) (fill:la m one.r)))
++  test-gte-t  ^-  tang  :: 1 >= 1
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (gte:la (fill:la m one.r) (fill:la m one.r)))
++  test-gte-f  ^-  tang  :: 1 >= 2
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (gte:la (fill:la m one.r) (fill:la m two.r)))
++  test-lth-t  ^-  tang  :: 1 < 2
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (lth:la (fill:la m one.r) (fill:la m two.r)))
++  test-lth-f  ^-  tang  :: 1 < 1
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (lth:la (fill:la m one.r) (fill:la m one.r)))
++  test-lte-t  ^-  tang  :: 1 <= 1
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m one.r) (lte:la (fill:la m one.r) (fill:la m one.r)))
++  test-lte-f  ^-  tang  :: 2 <= 1
  %-  each-real
  |=  [r=rrow m=meta]
  (is-equal (fill:la m zero.r) (lte:la (fill:la m two.r) (fill:la m one.r)))
++  test-gth-uint-t  ^-  tang  :: 2 > 1
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 1) (gth:la (fill:la m 2) (fill:la m 1)))
++  test-gth-uint-f  ^-  tang  :: 1 > 1
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 0) (gth:la (fill:la m 1) (fill:la m 1)))
++  test-lte-uint-t  ^-  tang  :: 1 <= 1
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 1) (lte:la (fill:la m 1) (fill:la m 1)))
++  test-lte-uint-f  ^-  tang  :: 2 <= 1
  %-  each-uint
  |=  m=meta
  (is-equal (fill:la m 0) (lte:la (fill:la m 2) (fill:la m 1)))
::
::  Reductions over [zero one two]: min=zero, max=two, cumsum=three (0+1+2),
::  argmin=0, argmax=2.  argmin/argmax assert the forward (ravel) index.
++  test-min  ^-  tang
  %-  each-vec
  |=  [r=rrow a=ray]
  (expect-eq !>(`@`zero.r) !>((hed (min:la a))))
++  test-max  ^-  tang
  %-  each-vec
  |=  [r=rrow a=ray]
  (expect-eq !>(`@`two.r) !>((hed (max:la a))))
++  test-cumsum  ^-  tang
  %-  each-vec
  |=  [r=rrow a=ray]
  (expect-eq !>(`@`three.r) !>((hed (cumsum:la a))))
++  test-argmin  ^-  tang
  %-  each-vec
  |=  [r=rrow a=ray]
  (expect-eq !>(`@ud`0) !>(`@ud`(argmin:la a)))
++  test-argmax  ^-  tang
  %-  each-vec
  |=  [r=rrow a=ray]
  (expect-eq !>(`@ud`2) !>(`@ud`(argmax:la a)))
::
::  any/all over boolean rays (numeric true=1.0, false=0.0).
++  test-any-all  ^-  tang
  =/  all-true   (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  has-false  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .0.0 .1.0]]])
  =/  all-false  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0]]])
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((all:la all-true))
    %+  expect-eq  !>(%.y)  !>((any:la all-true))
    %+  expect-eq  !>(%.n)  !>((all:la has-false))
    %+  expect-eq  !>(%.y)  !>((any:la has-false))
    %+  expect-eq  !>(%.n)  !>((all:la all-false))
    %+  expect-eq  !>(%.n)  !>((any:la all-false))
  ==
++  test-any-all-1d  ^-  tang
  =/  all-true   (en-ray:la [meta=[shape=~[3] bloq=5 kind=%i754 tail=~] baum=~[.1.0 .1.0 .1.0]])
  =/  has-false  (en-ray:la [meta=[shape=~[3] bloq=5 kind=%i754 tail=~] baum=~[.1.0 .0.0 .1.0]])
  =/  all-false  (en-ray:la [meta=[shape=~[3] bloq=5 kind=%i754 tail=~] baum=~[.0.0 .0.0 .0.0]])
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((all:la all-true))
    %+  expect-eq  !>(%.y)  !>((any:la all-true))
    %+  expect-eq  !>(%.n)  !>((all:la has-false))
    %+  expect-eq  !>(%.y)  !>((any:la has-false))
    %+  expect-eq  !>(%.n)  !>((all:la all-false))
    %+  expect-eq  !>(%.n)  !>((any:la all-false))
  ==
::  %cplx has no total order: the ordering reductions (max/min, and any/all
::  which read them) crash via +fun-scalar %gth -> +ord, rather than silently
::  mis-ordering.  Test complex rays with +is-close/+equ instead.
++  test-cplx-reduce-crashes  ^-  tang
  =/  c  (fill:la [shape=~[1 2] bloq=6 kind=%cplx tail=~] 0x3f80.0000)
  ;:  weld
    (expect-fail |.((max:la c)))
    (expect-fail |.((min:la c)))
    (expect-fail |.((any:la c)))
    (expect-fail |.((all:la c)))
  ==
--
