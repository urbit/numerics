/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-axis -- axis-wise reductions, moments, and norms
::
::  Table-driven over all four %i754 widths (bloq 4/5/6/7) and, where the arm is
::  kind-generic, over %uint widths too.
::
::  BIT-EXACTNESS.  Every input and every expected output here is a small
::  integer, exactly representable in binary16 through binary128, and every
::  intermediate (sums, one division for a mean, one square root for a std) is
::  exact as well.  So these expectations hold in all four widths AND in all
::  four rounding modes, and they are compared with `=` on the whole ray, not
::  within a tolerance.  Order- and rounding-sensitive behavior is pinned
::  separately in /tests/lib/lagoon-axis-rounding.
::
::  Oracle: NumPy (numpy 1.26.4), via lagoon/tools/gen-axis-oracle.py.
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
::  Builders.  +ints makes an %i754 ray of width .b from integers (via +change,
::  so no per-width literal table is needed); +uints makes a %uint ray.
++  ints
  |=  [shp=(list @) b=@ vals=(list @)]
  ^-  ray
  %+  reshape:la
    (change:la (en-ray:la [`meta`[~[(lent vals)] 6 %uint ~] vals]) %i754 b)
  shp
++  uints
  |=  [shp=(list @) b=@ vals=(list @)]
  ^-  ray
  (reshape:la (en-ray:la [`meta`[~[(lent vals)] b %uint ~] vals]) shp)
::  indices always come back as %uint bloq 6
++  idxs  |=([shp=(list @) vals=(list @)] ^-(ray (uints shp 6 vals)))
::
++  i754-bloqs  ^-  (list @)  ~[4 5 6 7]
++  uint-bloqs  ^-  (list @)  ~[3 4 5 6]
++  each-real  |=(fun=$-(@ tang) ^-(tang (zing (turn i754-bloqs fun))))
++  each-uint  |=(fun=$-(@ tang) ^-(tang (zing (turn uint-bloqs fun))))
::
::  Fixtures.
::    +mat-a:  2x3  [[1 2 3] [3 4 5]]  -- column means 2,3,4 and row means 2,4
::             are integers, so +mean-dim stays exact in every width.
++  mat-a  |=(b=@ ^-(ray (ints ~[2 3] b ~[1 2 3 3 4 5])))
++  mat-au  |=(b=@ ^-(ray (uints ~[2 3] b ~[1 2 3 3 4 5])))
::    +mat-b:  2x2  [[1 3] [2 6]]  -- row variances 1,4 (ddof 0) and 2,8
::             (ddof 1), row stds 1,2.
++  mat-b  |=(b=@ ^-(ray (ints ~[2 2] b ~[1 3 2 6])))
::    +vec-c:  [1 3]  -- whole-array mean 2, var 1 (ddof 0) / 2 (ddof 1), std 1.
++  vec-c  |=(b=@ ^-(ray (ints ~[2] b ~[1 3])))
::    +vec-d:  [3 4]  -- l1 7, l2 5, linf 4.
++  vec-d  |=(b=@ ^-(ray (ints ~[2] b ~[3 4])))
::    +mat-e:  2x2  [[3 4] [6 8]]  -- row l2 norms 5,10; row l1 7,14.
++  mat-e  |=(b=@ ^-(ray (ints ~[2 2] b ~[3 4 6 8])))
::    +mat-f:  2x2  [[1 2] [2 4]]  -- Frobenius norm exactly 5.
++  mat-f  |=(b=@ ^-(ray (ints ~[2 2] b ~[1 2 2 4])))
::
::  +sum-dim / +prod-dim
::
++  test-sum-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[4 6 8]) (sum-dim:la (mat-a b) 0))
++  test-sum-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[6 12]) (sum-dim:la (mat-a b) 1))
++  test-sum-dim-rank-1  ^-  tang
  ::  reducing a rank-1 ray gives shape ~[1], not a rank-0 ray
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[4]) (sum-dim:la (vec-c b) 0))
++  test-sum-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (uints ~[3] b ~[4 6 8]) (sum-dim:la (mat-au b) 0))
++  test-prod-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[3 8 15]) (prod-dim:la (mat-a b) 0))
++  test-prod-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[6 60]) (prod-dim:la (mat-a b) 1))
++  test-prod-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (uints ~[2] b ~[6 60]) (prod-dim:la (mat-au b) 1))
::
::  +max-dim / +min-dim
::
++  test-max-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[3 4 5]) (max-dim:la (mat-a b) 0))
++  test-max-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[3 5]) (max-dim:la (mat-a b) 1))
++  test-min-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[1 2 3]) (min-dim:la (mat-a b) 0))
++  test-min-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[1 3]) (min-dim:la (mat-a b) 1))
++  test-max-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (uints ~[3] b ~[3 4 5]) (max-dim:la (mat-au b) 0))
++  test-min-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (uints ~[2] b ~[1 3]) (min-dim:la (mat-au b) 1))
::
::  +argmax-dim / +argmin-dim (results are %uint bloq 6 regardless of input)
::
++  test-argmax-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[3] ~[1 1 1]) (argmax-dim:la (mat-a b) 0))
++  test-argmax-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2] ~[2 2]) (argmax-dim:la (mat-a b) 1))
++  test-argmin-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[3] ~[0 0 0]) (argmin-dim:la (mat-a b) 0))
++  test-argmin-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2] ~[0 0]) (argmin-dim:la (mat-a b) 1))
++  test-argmax-dim-first-wins  ^-  tang
  ::  a tie must report the FIRST maximum, like +argmax
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[1] ~[0]) (argmax-dim:la (ints ~[3] b ~[5 5 1]) 0))
++  test-argmin-dim-first-wins  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[1] ~[0]) (argmin-dim:la (ints ~[3] b ~[1 1 5]) 0))
++  test-argmax-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (idxs ~[3] ~[1 1 1]) (argmax-dim:la (mat-au b) 0))
::
::  rank 3.  Reducing the MIDDLE dimension of a ~[2 2 2] ray strides over the
::  data (inner = 2), which the rank-2 cases above never exercise.
::  a[i][j][k] = 4i + 2j + k.
::
++  cub-i  |=(b=@ ^-(ray (ints ~[2 2 2] b ~[0 1 2 3 4 5 6 7])))
++  test-sum-dim-rank-3-mid  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[2 4 10 12]) (sum-dim:la (cub-i b) 1))
++  test-sum-dim-rank-3-first  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[4 6 8 10]) (sum-dim:la (cub-i b) 0))
++  test-sum-dim-rank-3-last  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[1 5 9 13]) (sum-dim:la (cub-i b) 2))
++  test-max-dim-rank-3-mid  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[2 3 6 7]) (max-dim:la (cub-i b) 1))
++  test-min-dim-rank-3-mid  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[0 1 4 5]) (min-dim:la (cub-i b) 1))
++  test-argmax-dim-rank-3-mid  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2 2] ~[1 1 1 1]) (argmax-dim:la (cub-i b) 1))
++  test-mean-dim-rank-3-mid  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[1 2 5 6]) (mean-dim:la (cub-i b) 1))
++  test-broadcast-rank-3-mid  ^-  tang
  ::  ~[2 1 2] -> ~[2 3 2]: the middle dimension repeats
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3 2] b ~[0 1 0 1 0 1 4 5 4 5 4 5])
  (broadcast-to:la (ints ~[2 1 2] b ~[0 1 4 5]) ~[2 3 2])
::
::  +mean-dim / +mean
::
++  test-mean-dim-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[2 3 4]) (mean-dim:la (mat-a b) 0))
++  test-mean-dim-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[2 4]) (mean-dim:la (mat-a b) 1))
++  test-mean-whole  ^-  tang
  ::  (1+2+3+3+4+5)/6 = 3
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1 1] b ~[3]) (mean:la (mat-a b)))
++  test-mean-rank-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[2]) (mean:la (vec-c b)))
::
::  +var-dim / +std-dim / +var / +std
::
++  test-var-dim-ddof-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[1 4]) (var-dim:la (mat-b b) 1 0))
++  test-var-dim-ddof-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[2 8]) (var-dim:la (mat-b b) 1 1))
++  test-std-dim-ddof-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[1 2]) (std-dim:la (mat-b b) 1 0))
++  test-var-whole-ddof-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[1]) (var:la (vec-c b) 0))
++  test-var-whole-ddof-1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[2]) (var:la (vec-c b) 1))
++  test-std-whole-ddof-0  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[1]) (std:la (vec-c b) 0))
::
::  +norm / +norm-dim
::
++  test-norm-l1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[7]) (norm:la (vec-d b) %l1))
++  test-norm-l2  ^-  tang
  ::  sqrt(9+16) = 5, exact
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[5]) (norm:la (vec-d b) %l2))
++  test-norm-linf  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1] b ~[4]) (norm:la (vec-d b) %linf))
++  test-norm-fro  ^-  tang
  ::  sqrt(1+4+4+16) = 5, exact
  %-  each-real
  |=  b=@
  (is-equal (ints ~[1 1] b ~[5]) (norm:la (mat-f b) %fro))
++  test-norm-l1-uint  ^-  tang
  ::  %l1 needs only +abs and +add, so it is kind-generic
  %-  each-uint
  |=  b=@
  (is-equal (uints ~[1] b ~[7]) (norm:la (uints ~[2] b ~[3 4]) %l1))
++  test-norm-dim-l2  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[5 10]) (norm-dim:la (mat-e b) 1 %l2))
++  test-norm-dim-l1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[7 14]) (norm-dim:la (mat-e b) 1 %l1))
++  test-norm-dim-linf  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[4 8]) (norm-dim:la (mat-e b) 1 %linf))
++  test-norm-dim-0-l1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2] b ~[9 12]) (norm-dim:la (mat-e b) 0 %l1))
++  test-norm-negative-abs  ^-  tang
  ::  norms take absolute values: |-3|,|4| gives l1 7 and l2 5
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[1] b ~[5])
  (norm:la (sub:la (ints ~[2] b ~[0 4]) (ints ~[2] b ~[3 0])) %l2)
::
::  +broadcast-to
::
++  test-broadcast-row  ^-  tang
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3] b ~[1 2 3 1 2 3])
  (broadcast-to:la (ints ~[1 3] b ~[1 2 3]) ~[2 3])
++  test-broadcast-col  ^-  tang
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3] b ~[1 1 1 2 2 2])
  (broadcast-to:la (ints ~[2 1] b ~[1 2]) ~[2 3])
++  test-broadcast-rank-up  ^-  tang
  ::  a rank-1 ray right-aligns against the target shape
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3] b ~[1 2 3 1 2 3])
  (broadcast-to:la (ints ~[3] b ~[1 2 3]) ~[2 3])
++  test-broadcast-scalar  ^-  tang
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 2] b ~[7 7 7 7])
  (broadcast-to:la (ints ~[1 1] b ~[7]) ~[2 2])
++  test-broadcast-identity  ^-  tang
  ::  broadcasting to the same shape is the identity
  %-  each-real
  |=  b=@
  (is-equal (mat-a b) (broadcast-to:la (mat-a b) ~[2 3]))
++  test-broadcast-uint  ^-  tang
  %-  each-uint
  |=  b=@
  %+  is-equal  (uints ~[2 3] b ~[1 1 1 2 2 2])
  (broadcast-to:la (uints ~[2 1] b ~[1 2]) ~[2 3])
++  test-broadcast-then-add  ^-  tang
  ::  the point of broadcasting: elementwise ops against a row vector
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3] b ~[2 4 6 4 6 8])
  (add:la (mat-a b) (broadcast-to:la (ints ~[1 3] b ~[1 2 3]) ~[2 3]))
--
