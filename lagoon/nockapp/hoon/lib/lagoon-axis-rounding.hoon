/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-axis-rounding -- order and rounding of the axis arms
::
::  The arms in /tests/lib/lagoon-axis and /tests/lib/lagoon-sort-dist use
::  exactly-representable integers, so they pass whatever the accumulation order
::  and rounding mode are.  These arms pin the two decisions those cannot see:
::
::    1. reductions fold LEFT TO RIGHT along the axis, seeded with the slice's
::       first element (a right fold gives different bits here), and
::    2. +cdist-sq sums squared differences DIRECTLY rather than via the Gram
::       identity |x|^2 + |y|^2 - 2*x.y (the two differ by an ulp here).
::
::  Both are jet contracts: a jet that folds the other way, or that reaches for
::  the Gram identity because it is faster, fails these.
::
::  Inputs and expectations are raw bit patterns.  Oracle: exact rational
::  arithmetic (fractions.Fraction) with an IEEE-754 encoder cross-checked
::  against NumPy on 4120 round-to-nearest probes; see
::  lagoon/tools/gen-axis-oracle.py.
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
::  build a ray from raw bit patterns (no conversion)
++  bits
  |=  [shp=(list @) b=@ vals=(list @)]
  ^-  ray
  (reshape:la (en-ray:la [`meta`[~[(lent vals)] b %i754 ~] vals]) shp)
::
::  ORDER.  v = [1, 2^-24, 2^-24] in binary32.  1 + 2^-24 is exactly the
::  midpoint between 1 and 1+2^-23, so folding LEFT TO RIGHT ties down to 1.0
::  twice under %n (and truncates to 1.0 under %z), while folding RIGHT TO LEFT
::  first forms the exact 2^-23 and yields 1+2^-23 = 0x3f80.0001 in every mode.
::  So 0x3f80.0000 here is positive evidence of the left fold.
++  v32  ^-  ray  (bits ~[3] 5 ~[0x3f80.0000 0x3380.0000 0x3380.0000])
::
++  test-sum-dim-order-nearest  ^-  tang
  ::  left fold, %n: 1.0 (a right fold would give 0x3f80.0001)
  %+  is-equal
    `ray`[[~[1] 5 %i754 ~] 0x1.3f80.0000]
  (sum-dim:(lake %n) v32 0)
++  test-sum-dim-order-zero  ^-  tang
  ::  left fold, %z: both additions truncate back to 1.0
  %+  is-equal
    `ray`[[~[1] 5 %i754 ~] 0x1.3f80.0000]
  (sum-dim:(lake %z) v32 0)
++  test-sum-dim-order-down  ^-  tang
  ::  left fold, %d: toward -inf, also 1.0
  %+  is-equal
    `ray`[[~[1] 5 %i754 ~] 0x1.3f80.0000]
  (sum-dim:(lake %d) v32 0)
++  test-sum-dim-order-up  ^-  tang
  ::  left fold, %u: 1 -> 1+2^-23 -> 1+2^-22 = 0x3f80.0002
  %+  is-equal
    `ray`[[~[1] 5 %i754 ~] 0x1.3f80.0002]
  (sum-dim:(lake %u) v32 0)
++  test-sum-dim-mode-is-honored  ^-  tang
  ::  the door's mode reaches the reduction: %u and %z disagree
  ?.  =((sum-dim:(lake %u) v32 0) (sum-dim:(lake %z) v32 0))  ~
  ~[leaf+"sum-dim ignored the rounding mode"]
::
::  DIRECT vs GRAM.  A = [1.375, 1.875] and B = [16/3, 5/3] rounded to the
::  width.  The direct sum of squared differences and the Gram identity differ
::  by one ulp in both binary32 and binary64.
::
++  test-cdist-sq-direct-not-gram-32  ^-  tang
  ::  direct 0x417b.638f; the Gram identity would give 0x417b.6390
  %+  is-equal
    `ray`[[~[1 1] 5 %i754 ~] 0x1.417b.638f]
  %+  cdist-sq:(lake %n)
    (bits ~[1 2] 5 ~[0x3fb0.0000 0x3ff0.0000])
  (bits ~[1 2] 5 ~[0x40aa.aaab 0x3fd5.5555])
++  test-cdist-sq-direct-not-gram-64  ^-  tang
  ::  direct 0x402f.6c71.c71c.71c6; the Gram identity would give ...71c8
  %+  is-equal
    `ray`[[~[1 1] 6 %i754 ~] 0x1.402f.6c71.c71c.71c6]
  %+  cdist-sq:(lake %n)
    (bits ~[1 2] 6 ~[0x3ff6.0000.0000.0000 0x3ffe.0000.0000.0000])
  (bits ~[1 2] 6 ~[0x4015.5555.5555.5555 0x3ffa.aaaa.aaaa.aaab])
--
