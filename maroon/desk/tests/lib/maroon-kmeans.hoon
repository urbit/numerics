/-  ls=lagoon
/+  *test, *lagoon, maroon, rand
::::  /tests/lib/maroon-kmeans -- assignment, Lloyd's algorithm, k-means++
::
::  Lloyd's algorithm is deterministic given its starting centroids, so most of
::  this compares whole rays exactly: the fixtures are small integers whose
::  cluster means are exact halves, verified in exact rational arithmetic by
::  maroon/tools/ml_check.py (which also checks the partition against
::  sklearn.cluster.KMeans).
::
::  The seeding is random, so it is tested for the properties that must hold
::  regardless of the draw: same seed gives the same centroids, the seeds are
::  distinct samples, and the fit recovers a well-separated partition.  Labels
::  are compared up to permutation, since cluster numbering is arbitrary.
::
|%
++  mm  (make:maroon %n .~1e-12)
++  lk  (lake %n)
++  mat
  |=  [r=@ c=@ v=(list (list @))]
  ^-  ray:ls
  (en-ray:lk [[~[r c] 6 %i754 ~] v])
++  vec
  |=  v=(list @)
  ^-  ray:ls
  (en-ray:lk [[~[(lent v)] 6 %i754 ~] v])
::  labels come back as %uint bloq 6
++  labs
  |=  v=(list @)
  ^-  ray:ls
  (en-ray:lk [[~[(lent v)] 6 %uint ~] v])
++  seed  (from-atom:seed:rand %sm64 0)
::  fixtures
::    +line4:  1-D, two obvious clusters {0,1} and {10,11}
++  line4  (mat 4 1 ~[~[.~0] ~[.~1] ~[.~10] ~[.~11]])
++  init2  (mat 2 1 ~[~[.~0] ~[.~10]])
::    +rect4:  2-D, two clusters two units apart in y, ten apart in x
++  rect4  (mat 4 2 ~[~[.~0 .~0] ~[.~0 .~2] ~[.~10 .~0] ~[.~10 .~2]])
++  init2d  (mat 2 2 ~[~[.~0 .~0] ~[.~10 .~0]])
::
::  +assign
::
++  test-assign  ^-  tang
  (expect-eq !>((labs ~[0 0 1 1])) !>((assign:mm line4 init2)))
++  test-assign-tie  ^-  tang
  ::  a sample equidistant from two centroids takes the lower index
  (expect-eq !>((labs ~[0])) !>((assign:mm (mat 1 1 ~[~[.~5]]) init2)))
++  test-assign-2d  ^-  tang
  (expect-eq !>((labs ~[0 0 1 1])) !>((assign:mm rect4 init2d)))
::
::  +update
::
++  test-update  ^-  tang
  ::  cluster means: (0+1)/2 and (10+11)/2, both exact
  %+  expect-eq  !>((mat 2 1 ~[~[.~0.5] ~[.~10.5]]))
  !>((update:mm line4 (labs ~[0 0 1 1]) init2 2))
++  test-update-empty-cluster  ^-  tang
  ::  an empty cluster keeps its previous centroid
  %+  expect-eq  !>((mat 2 1 ~[~[.~5.5] ~[.~99]]))
  !>((update:mm line4 (labs ~[0 0 0 0]) (mat 2 1 ~[~[.~0] ~[.~99]]) 2))
::
::  +kmeans
::
++  test-kmeans-labels  ^-  tang
  (expect-eq !>((labs ~[0 0 1 1])) !>(labels:(kmeans:mm line4 init2 20)))
++  test-kmeans-centroids  ^-  tang
  %+  expect-eq  !>((mat 2 1 ~[~[.~0.5] ~[.~10.5]]))
  !>(centroids:(kmeans:mm line4 init2 20))
++  test-kmeans-iterations  ^-  tang
  ::  this configuration reaches its fixed point in one iteration
  (expect-eq !>(1) !>(iter:(kmeans:mm line4 init2 20)))
++  test-kmeans-2d  ^-  tang
  =/  res  (kmeans:mm rect4 init2d 20)
  %+  weld  (expect-eq !>((labs ~[0 0 1 1])) !>(labels.res))
  (expect-eq !>((mat 2 2 ~[~[.~0 .~1] ~[.~10 .~1]])) !>(centroids.res))
++  test-kmeans-maxit-zero  ^-  tang
  ::  no iterations allowed: the initial assignment, unmoved centroids
  =/  res  (kmeans:mm line4 init2 0)
  %+  weld  (expect-eq !>((labs ~[0 0 1 1])) !>(labels.res))
  %+  weld  (expect-eq !>(init2) !>(centroids.res))
  (expect-eq !>(0) !>(iter.res))
++  test-kmeans-idempotent  ^-  tang
  ::  restarting from the converged centroids changes nothing
  =/  res  (kmeans:mm line4 init2 20)
  (expect-eq !>(centroids.res) !>(centroids:(kmeans:mm line4 centroids.res 20)))
++  test-kmeans-crashes-feature-mismatch  ^-  tang
  (expect-fail |.((kmeans:mm line4 init2d 5)))
::
::  +inertia
::
++  test-inertia  ^-  tang
  ::  every point sits 0.5 from its centroid, so 4 * 0.25
  =/  res  (kmeans:mm line4 init2 20)
  (expect-eq !>(`@rd`.~1) !>(`@rd`(inertia:mm line4 labels.res centroids.res)))
++  test-inertia-zero-at-centroids  ^-  tang
  ::  samples that ARE the centroids have zero inertia
  %+  expect-eq  !>(`@rd`.~0)
  !>(`@rd`(inertia:mm init2 (labs ~[0 1]) init2))
::
::  +kmeans-pp and +kmeans-fit
::
++  test-kmeans-pp-deterministic  ^-  tang
  ::  the same seed gives the same centroids
  %+  expect-eq  !>(centroids:(kmeans-pp:mm line4 2 seed))
  !>(centroids:(kmeans-pp:mm line4 2 seed))
++  test-kmeans-pp-distinct  ^-  tang
  ::  k-means++ never picks the same sample twice
  =/  c  centroids:(kmeans-pp:mm line4 2 seed)
  (expect !>(!=((get-item:lk c ~[0 0]) (get-item:lk c ~[1 0]))))
++  test-kmeans-pp-picks-samples  ^-  tang
  ::  every seed is one of the samples
  =/  c  centroids:(kmeans-pp:mm line4 2 seed)
  =/  seen  (ravel:lk line4)
  =/  a  (get-item:lk c ~[0 0])
  =/  b  (get-item:lk c ~[1 0])
  (expect !>(?&(!=(~ (find ~[a] seen)) !=(~ (find ~[b] seen)))))
++  test-kmeans-pp-shape  ^-  tang
  =/  c  centroids:(kmeans-pp:mm rect4 2 seed)
  (expect !>(=(~[2 2] shape.meta.c)))
++  test-kmeans-fit-partition  ^-  tang
  ::  a well-separated dataset is recovered whatever the seeding, up to the
  ::  arbitrary numbering of the clusters
  =/  res  (kmeans-fit:mm line4 2 20 seed)
  =/  l  labels.res
  =/  a  (get-item:lk l ~[0])
  =/  b  (get-item:lk l ~[1])
  =/  c  (get-item:lk l ~[2])
  =/  d  (get-item:lk l ~[3])
  (expect !>(?&(=(a b) =(c d) !=(a c))))
++  test-kmeans-fit-centroids  ^-  tang
  ::  and it finds the same two centroids, in some order
  =/  res  (kmeans-fit:mm line4 2 20 seed)
  =/  c  centroids.res
  =/  lo  (get-item:lk c ~[0 0])
  =/  hi  (get-item:lk c ~[1 0])
  =/  half  `@ux`.~0.5
  =/  ten-half  `@ux`.~10.5
  (expect !>(?|(?&(=(lo half) =(hi ten-half)) ?&(=(lo ten-half) =(hi half)))))
++  test-kmeans-pp-crashes-k-too-big  ^-  tang
  (expect-fail |.((kmeans-pp:mm line4 5 seed)))
--
