/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-sort-dist -- sorting, selection, and distances
::
::  Table-driven over all four %i754 widths (bloq 4/5/6/7) and, for the
::  comparison-only arms, over %uint widths too.
::
::  BIT-EXACTNESS.  Sorting and selection move values without arithmetic, so
::  they are exact by construction.  The distance cases use small integers whose
::  squares and sums are exact in binary16 through binary128, so the
::  expectations hold in every width and every rounding mode.  The cases that
::  turn on rounding live in /tests/lib/lagoon-axis-rounding.
::
::  Oracle: NumPy (numpy 1.26.4), via lagoon/tools/gen-axis-oracle.py; the
::  stable-tie expectations match np.argsort(kind='stable').
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
++  idxs  |=([shp=(list @) vals=(list @)] ^-(ray (uints shp 6 vals)))
::
++  i754-bloqs  ^-  (list @)  ~[4 5 6 7]
++  uint-bloqs  ^-  (list @)  ~[3 4 5 6]
++  each-real  |=(fun=$-(@ tang) ^-(tang (zing (turn i754-bloqs fun))))
++  each-uint  |=(fun=$-(@ tang) ^-(tang (zing (turn uint-bloqs fun))))
::
::  Fixtures.
++  vec-s   |=(b=@ ^-(ray (ints ~[3] b ~[3 1 2])))
++  vec-su  |=(b=@ ^-(ray (uints ~[3] b ~[3 1 2])))
::    +vec-t:  [2 1 2] -- the duplicate pins tie-breaking in both directions
++  vec-t   |=(b=@ ^-(ray (ints ~[3] b ~[2 1 2])))
++  mat-g   |=(b=@ ^-(ray (ints ~[2 2] b ~[3 1 1 2])))
++  mat-h   |=(b=@ ^-(ray (ints ~[2 3] b ~[1 5 3 9 2 7])))
::    rank-3, a[i][j][k] = 4i+2j+k: exercises inner > 1 (dim 1 of ~[2 2 2])
++  cub-i   |=(b=@ ^-(ray (ints ~[2 2 2] b ~[0 1 2 3 4 5 6 7])))
::    +pts-a / +pts-b:  2x2 and 3x2 point sets with integer squared distances
++  pts-a   |=(b=@ ^-(ray (ints ~[2 2] b ~[0 0 1 0])))
++  pts-b   |=(b=@ ^-(ray (ints ~[3 2] b ~[0 0 0 1 1 1])))
::
::  +sort-dim
::
++  test-sort-dim-asc  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[1 2 3]) (sort-dim:la (vec-s b) 0 %asc))
++  test-sort-dim-des  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (ints ~[3] b ~[3 2 1]) (sort-dim:la (vec-s b) 0 %des))
++  test-sort-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (uints ~[3] b ~[1 2 3]) (sort-dim:la (vec-su b) 0 %asc))
++  test-sort-dim-cols  ^-  tang
  ::  dim 0 sorts each COLUMN: [[3 1] [1 2]] -> [[1 1] [3 2]]
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[1 1 3 2]) (sort-dim:la (mat-g b) 0 %asc))
++  test-sort-dim-rows  ^-  tang
  ::  dim 1 sorts each ROW: [[3 1] [1 2]] -> [[1 3] [1 2]]
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[1 3 1 2]) (sort-dim:la (mat-g b) 1 %asc))
++  test-sort-dim-rank-3  ^-  tang
  ::  dim 1 of a ~[2 2 2] ray, so the slices are strided (inner = 2)
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 2 2] b ~[2 3 0 1 6 7 4 5])
  (sort-dim:la (cub-i b) 1 %des)
++  test-sort-dim-shape-kept  ^-  tang
  ::  sorting does not change the shape or the kind
  %-  each-real
  |=  b=@
  =/  a=ray  (mat-h b)
  =/  s=ray  (sort-dim:la a 1 %asc)
  ?:  =(meta.s meta.a)  ~
  ~[leaf+"sort-dim changed the meta"]
::
::  +argsort-dim (indices are %uint bloq 6)
::
++  test-argsort-dim-asc  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[3] ~[1 2 0]) (argsort-dim:la (vec-s b) 0 %asc))
++  test-argsort-dim-des  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[3] ~[0 2 1]) (argsort-dim:la (vec-s b) 0 %des))
++  test-argsort-dim-tie-asc  ^-  tang
  ::  [2 1 2] ascending: the tied 2s keep their order, so index 0 before 2
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[3] ~[1 0 2]) (argsort-dim:la (vec-t b) 0 %asc))
++  test-argsort-dim-tie-des  ^-  tang
  ::  [2 1 2] descending: the tie STILL resolves to the lower index first
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[3] ~[0 2 1]) (argsort-dim:la (vec-t b) 0 %des))
++  test-argsort-dim-cols  ^-  tang
  ::  dim 0 of [[3 1] [1 2]]: column 0 is [3 1] -> [1 0], column 1 is [1 2] -> [0 1]
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2 2] ~[1 0 0 1]) (argsort-dim:la (mat-g b) 0 %asc))
++  test-argsort-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (idxs ~[3] ~[1 2 0]) (argsort-dim:la (vec-su b) 0 %asc))
::
::  +argtop-dim (replaces .dim with .k rather than dropping it)
::
++  test-argtop-dim-k2  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2 2] ~[1 2 0 2]) (argtop-dim:la (mat-h b) 1 2))
++  test-argtop-dim-k1  ^-  tang
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2 1] ~[1 0]) (argtop-dim:la (mat-h b) 1 1))
++  test-argtop-dim-dim-0  ^-  tang
  ::  down the columns of [[1 5 3] [9 2 7]]: 9 beats 1, 5 beats 2, 7 beats 3
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[1 3] ~[1 0 1]) (argtop-dim:la (mat-h b) 0 1))
++  test-argtop-dim-full  ^-  tang
  ::  k = the axis length is just a descending argsort
  %-  each-real
  |=  b=@
  (is-equal (argsort-dim:la (mat-h b) 1 %des) (argtop-dim:la (mat-h b) 1 3))
++  test-argtop-dim-tie  ^-  tang
  ::  tied maxima: the lower index comes first
  %-  each-real
  |=  b=@
  (is-equal (idxs ~[2] ~[0 1]) (argtop-dim:la (ints ~[3] b ~[5 5 1]) 0 2))
++  test-argtop-dim-uint  ^-  tang
  %-  each-uint
  |=  b=@
  (is-equal (idxs ~[2 2] ~[1 2 0 2]) (argtop-dim:la (uints ~[2 3] b ~[1 5 3 9 2 7]) 1 2))
::
::  +take-dim
::
++  test-take-dim-top-k  ^-  tang
  ::  the two largest of each row of [[1 5 3] [9 2 7]], largest first
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 2] b ~[5 3 9 7])
  (take-dim:la (mat-h b) (argtop-dim:la (mat-h b) 1 2) 1)
++  test-take-dim-argsort-round-trip  ^-  tang
  ::  take along the argsort permutation == sort
  %-  each-real
  |=  b=@
  %+  is-equal  (sort-dim:la (mat-h b) 1 %asc)
  (take-dim:la (mat-h b) (argsort-dim:la (mat-h b) 1 %asc) 1)
++  test-take-dim-cols  ^-  tang
  ::  pick row 1 of column 0 and row 0 of columns 1,2 out of [[1 5 3] [9 2 7]]
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[1 3] b ~[9 5 3])
  (take-dim:la (mat-h b) (idxs ~[1 3] ~[1 0 0]) 0)
++  test-take-dim-repeat  ^-  tang
  ::  indices may repeat and need not be a permutation
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3] b ~[5 5 1 2 2 9])
  (take-dim:la (mat-h b) (idxs ~[2 3] ~[1 1 0 1 1 0]) 1)
++  test-take-dim-rank-3  ^-  tang
  ::  strided slices: reverse dim 1 of the ~[2 2 2] ray by index
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 2 2] b ~[2 3 0 1 6 7 4 5])
  (take-dim:la (cub-i b) (idxs ~[2 2 2] ~[1 1 0 0 1 1 0 0]) 1)
::
::  +cdist-sq / +pdist-sq
::
++  test-cdist-sq  ^-  tang
  ::  rows of [[0 0] [1 0]] against [[0 0] [0 1] [1 1]]
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[2 3] b ~[0 1 2 1 2 1])
  (cdist-sq:la (pts-a b) (pts-b b))
++  test-cdist-sq-single  ^-  tang
  ::  3-4-5: the squared distance is 25
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[1 1] b ~[25])
  (cdist-sq:la (ints ~[1 2] b ~[0 0]) (ints ~[1 2] b ~[3 4]))
++  test-pdist-sq  ^-  tang
  ::  self-distances: zero diagonal, symmetric
  %-  each-real
  |=  b=@
  (is-equal (ints ~[2 2] b ~[0 1 1 0]) (pdist-sq:la (pts-a b)))
++  test-pdist-sq-diagonal-zero  ^-  tang
  ::  Read the diagonal with +get-item rather than +diag.  On a ship whose vere
  ::  still carries the PRE-urbit/vere#1057 lagoon jets, +diag crashes: that old
  ::  C jet crashed unconditionally (see numerics #75 and urbit/urbit#7388).  It
  ::  is the jet, not this desk's Hoon -- an inline copy of +diag's own body,
  ::  outside the jetted core, runs fine on the same ship -- and it reproduces on
  ::  main, so it has nothing to do with the arm under test here.
  %-  each-real
  |=  b=@
  =/  p  (pdist-sq:la (pts-b b))
  ?:  ?&  =(0 (get-item:la p ~[0 0]))
          =(0 (get-item:la p ~[1 1]))
          =(0 (get-item:la p ~[2 2]))
      ==
    ~
  ~[leaf+"pdist-sq diagonal is not +0.0"]
++  test-cdist-sq-wide  ^-  tang
  ::  d = 3, so the inner sum has three terms: |(1,2,3)-(0,0,0)|^2 = 14
  %-  each-real
  |=  b=@
  %+  is-equal  (ints ~[1 1] b ~[14])
  (cdist-sq:la (ints ~[1 3] b ~[1 2 3]) (ints ~[1 3] b ~[0 0 0]))
--
