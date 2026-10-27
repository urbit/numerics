/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-unary -- unary / elementwise transforms
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

::  el-wise-op must preserve element order.  It previously flopped the
::  raveled list ("compensate for LSB"), which reversed every Hoon-only
::  transcendental (exp/sin/cos/tan/abs) for n>1 elements.  Identity on
::  [1 2 3 4] must return [1 2 3 4], not [4 3 2 1].
++  test-el-wise-op-order-5r  ^-  tang
  =/  meta-1x4-5  [~[4] 5 %i754 ~]
  =/  input-1x4-5  (en-ray:la [meta-1x4-5 ~[.1.0 .2.0 .3.0 .4.0]])
  ;:  weld
    %+  expect-eq
      !>(input-1x4-5)
      !>((el-wise-op:la input-1x4-5 |=(e=@ e)))
  ==
::

::  abs flows through el-wise-op; distinct magnitudes and mixed signs make
::  any order reversal visible in the public single-arg path.
++  test-abs-1x4-5r  ^-  tang
  =/  meta-1x4-5  [~[4] 5 %i754 ~]
  =/  input-1x4-5  (en-ray:la [meta-1x4-5 ~[.-1.0 .2.0 .-3.0 .4.0]])
  =/  canon-1x4-5  (en-ray:la [meta-1x4-5 ~[.1.0 .2.0 .3.0 .4.0]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x4-5)
      !>((abs:la input-1x4-5))
  ==
--
