/-  ls=lagoon
/+  *test, *saloon, *lagoon
::::  /tests/lib/saloon-linalg -- Cholesky, conjugate gradient, one-sided SVD
::
::  Two kinds of case.
::
::  EXACT: whole rays compared with +expect-eq.  Every input and every
::  intermediate (sums, the one division per row, the one square root per pivot)
::  is exactly representable in binary16 through binary128, so these hold at any
::  width and in any rounding mode.  Verified operation by operation in exact
::  rational arithmetic by saloon/tools/linalg_check.py.
::
::  APPROXIMATE: compared with +close against NumPy references (numpy.linalg.svd
::  and a float64 CG mirror), also from that script.
::
::  +tpose is local rather than +transpose:la because the lagoon transpose jet
::  crashes on runtimes older than urbit/vere#1057.
::
|%
++  sad  (sake %n .~1e-12)
++  lk   (lake %n)
++  close
  |=  [x=ray:ls y=ray:ls]
  ^-  ?
  (all:lk (is-close:lk x y [.~1e-9 .~1e-9]))
::  builders
++  mat
  |=  [r=@ c=@ v=(list (list @))]
  ^-  ray:ls
  (en-ray:lk [[~[r c] 6 %i754 ~] v])
++  vec
  |=  v=(list @)
  ^-  ray:ls
  (en-ray:lk [[~[(lent v)] 6 %i754 ~] v])
++  mat5
  |=  [r=@ c=@ v=(list (list @))]
  ^-  ray:ls
  (en-ray:lk [[~[r c] 5 %i754 ~] v])
++  vec5
  |=  v=(list @)
  ^-  ray:ls
  (en-ray:lk [[~[(lent v)] 5 %i754 ~] v])
::  transpose without the broken jet
++  tpose
  |=  a=ray:ls
  ^-  ray:ls
  =/  r  (snag 0 shape.meta.a)
  =/  c  (snag 1 shape.meta.a)
  =/  o  (zeros:lk [~[c r] bloq.meta.a %i754 ~])
  =/  i  0
  |-  ^-  ray:ls
  ?:  =(i r)  o
  =/  row
    =/  j  0
    |-  ^-  ray:ls
    ?:  =(j c)  o
    $(j +(j), o (set-item:lk o ~[j i] (get-item:lk a ~[i j])))
  $(i +(i), o row)
::  n x n diagonal matrix from a rank-1 ray
++  mk-diag
  |=  [d=ray:ls n=@]
  ^-  ray:ls
  =/  z  (zeros:lk [~[n n] 6 %i754 ~])
  =/  i  0
  |-  ^-  ray:ls
  ?:  =(i n)  z
  =.  z  (set-item:lk z ~[i i] (get-item:lk d ~[i]))
  $(i +(i))
::
::  fixtures
++  spd2   (mat 2 2 ~[~[.~4 .~2] ~[.~2 .~5]])        ::  L = [[2 0] [1 2]]
++  spd2b  (mat 2 2 ~[~[.~9 .~3] ~[.~3 .~5]])        ::  L = [[3 0] [1 2]]
++  spd3   (mat 3 3 ~[~[.~4 .~2 .~2] ~[.~2 .~5 .~3] ~[.~2 .~3 .~6]])
++  lo3    (mat 3 3 ~[~[.~2 .~0 .~0] ~[.~1 .~2 .~0] ~[.~1 .~1 .~2]])
++  eye2   (mat 2 2 ~[~[.~1 .~0] ~[.~0 .~1]])
++  dg24   (mat 2 2 ~[~[.~2 .~0] ~[.~0 .~4]])
++  gen22  (mat 2 2 ~[~[.~1 .~2] ~[.~3 .~4]])
++  gen43  (mat 4 3 ~[~[.~1 .~2 .~3] ~[.~4 .~5 .~6] ~[.~7 .~8 .~10] ~[.~2 .~0 .~1]])
::
::  +chol
::
++  test-chol-2x2  ^-  tang
  (expect-eq !>((mat 2 2 ~[~[.~2 .~0] ~[.~1 .~2]])) !>((chol:sad spd2)))
++  test-chol-2x2-other  ^-  tang
  (expect-eq !>((mat 2 2 ~[~[.~3 .~0] ~[.~1 .~2]])) !>((chol:sad spd2b)))
++  test-chol-3x3  ^-  tang
  (expect-eq !>(lo3) !>((chol:sad spd3)))
++  test-chol-bloq-5  ^-  tang
  ::  the same factorization is exact in binary32 too
  %+  expect-eq  !>((mat5 2 2 ~[~[.2 .0] ~[.1 .2]]))
  !>((chol:(sake %n .1e-6) (mat5 2 2 ~[~[.4 .2] ~[.2 .5]])))
++  test-chol-reconstructs  ^-  tang
  ::  L*L^T = A, exactly
  =/  l  (chol:sad spd3)
  (expect-eq !>(spd3) !>((mmul:lk l (tpose l))))
++  test-chol-unit-not-pd  ^-  tang
  (expect !>(=(~ (chol-unit:sad (mat 2 2 ~[~[.~1 .~2] ~[.~2 .~1]])))))
++  test-chol-unit-zero-pivot  ^-  tang
  ::  a pivot of exactly 0 is rejected: [[4 2] [2 1]] is singular
  (expect !>(=(~ (chol-unit:sad (mat 2 2 ~[~[.~4 .~2] ~[.~2 .~1]])))))
++  test-chol-unit-zeros  ^-  tang
  (expect !>(=(~ (chol-unit:sad (zeros:lk [~[2 2] 6 %i754 ~])))))
++  test-chol-crashes-not-pd  ^-  tang
  (expect-fail |.((chol:sad (mat 2 2 ~[~[.~1 .~2] ~[.~2 .~1]]))))
++  test-chol-crashes-asymmetric  ^-  tang
  (expect-fail |.((chol:sad (mat 2 2 ~[~[.~1 .~2] ~[.~3 .~1]]))))
::
::  triangular solves and +chol-solve
::
++  test-trsv-lo  ^-  tang
  ::  [[2 0] [1 2]] y = [10 9]  ->  y = [5 2]
  (expect-eq !>((vec ~[.~5 .~2])) !>((trsv-lo:sad (mat 2 2 ~[~[.~2 .~0] ~[.~1 .~2]]) (vec ~[.~10 .~9]))))
++  test-trsv-up  ^-  tang
  ::  [[2 0] [1 2]]^T x = [5 2]  ->  x = [2 1]
  (expect-eq !>((vec ~[.~2 .~1])) !>((trsv-up:sad (mat 2 2 ~[~[.~2 .~0] ~[.~1 .~2]]) (vec ~[.~5 .~2]))))
++  test-chol-solve  ^-  tang
  (expect-eq !>((vec ~[.~2 .~1])) !>((chol-solve:sad spd2 (vec ~[.~10 .~9]))))
++  test-chol-solve-3x3  ^-  tang
  ::  A x = A*[1 1 1] with A = spd3, so x = [1 1 1]
  %+  expect-eq  !>((vec ~[.~1 .~1 .~1]))
  !>((chol-solve:sad spd3 (matvec:sad spd3 (vec ~[.~1 .~1 .~1]))))
::
::  vector helpers
::
++  test-matvec  ^-  tang
  (expect-eq !>((vec ~[.~10 .~9])) !>((matvec:sad spd2 (vec ~[.~2 .~1]))))
++  test-matvec-rectangular  ^-  tang
  ::  a 4x3 times a 3-vector is a 4-vector
  %+  expect-eq  !>((vec ~[.~14 .~32 .~53 .~5]))
  !>((matvec:sad gen43 (vec ~[.~1 .~2 .~3])))
++  test-dotv  ^-  tang
  (expect-eq !>(`@rd`.~25) !>(`@rd`(dotv:sad (vec ~[.~3 .~4]) (vec ~[.~3 .~4]))))
++  test-nrm2  ^-  tang
  (expect-eq !>(`@rd`.~5) !>(`@rd`(nrm2:sad (vec ~[.~3 .~4]))))
++  test-axpyv  ^-  tang
  ::  y + 2*x with x = [1 2], y = [3 4]
  (expect-eq !>((vec ~[.~5 .~8])) !>((axpyv:sad .~2 (vec ~[.~1 .~2]) (vec ~[.~3 .~4]))))
::
::  +cg / +pcg
::
++  test-cg-identity  ^-  tang
  ::  A = I converges in one iteration, exactly
  =/  res  (cg:sad eye2 (vec ~[.~3 .~4]) 20)
  %+  weld  (expect-eq !>((vec ~[.~3 .~4])) !>(x.res))
  (expect-eq !>(1) !>(iter.res))
++  test-cg-diagonal  ^-  tang
  ::  diag(2,4) x = e1 -> x = [0.5 0], one iteration, exactly
  =/  res  (cg:sad dg24 (vec ~[.~1 .~0]) 20)
  %+  weld  (expect-eq !>((vec ~[.~0.5 .~0])) !>(x.res))
  (expect-eq !>(1) !>(iter.res))
++  test-cg-maxit-zero  ^-  tang
  ::  no iterations allowed: x stays 0 and the residual is |v|
  =/  res  (cg:sad dg24 (vec ~[.~2 .~4]) 0)
  %+  weld  (expect-eq !>((zeros:lk [~[2] 6 %i754 ~])) !>(x.res))
  (expect-eq !>(0) !>(iter.res))
++  test-cg-zero-rhs  ^-  tang
  ::  v = 0 is already solved
  =/  res  (cg:sad spd2 (zeros:lk [~[2] 6 %i754 ~]) 20)
  %+  weld  (expect-eq !>((zeros:lk [~[2] 6 %i754 ~])) !>(x.res))
  (expect-eq !>(0) !>(iter.res))
++  test-cg-converges  ^-  tang
  ::  [[4 2] [2 5]] x = [10 9] -> [2 1] in at most n=2 iterations
  =/  res  (cg:sad spd2 (vec ~[.~10 .~9]) 20)
  %+  weld  (expect !>((close (vec ~[.~2 .~1]) x.res)))
  (expect !>((lte iter.res 2)))
++  test-cg-matches-chol-solve  ^-  tang
  ::  the two solvers agree on the same system
  =/  v  (matvec:sad spd3 (vec ~[.~1 .~-2 .~3]))
  (expect !>((close (chol-solve:sad spd3 v) x:(cg:sad spd3 v 50))))
++  test-cg-crashes-asymmetric  ^-  tang
  (expect-fail |.((cg:sad (mat 2 2 ~[~[.~1 .~2] ~[.~3 .~1]]) (vec ~[.~1 .~1]) 5)))
++  test-pcg-identity  ^-  tang
  =/  res  (pcg:sad eye2 (vec ~[.~3 .~4]) 20)
  %+  weld  (expect-eq !>((vec ~[.~3 .~4])) !>(x.res))
  (expect-eq !>(1) !>(iter.res))
++  test-pcg-diagonal  ^-  tang
  ::  the Jacobi preconditioner solves a diagonal system in one step, exactly
  =/  res  (pcg:sad dg24 (vec ~[.~1 .~0]) 20)
  %+  weld  (expect-eq !>((vec ~[.~0.5 .~0])) !>(x.res))
  (expect-eq !>(1) !>(iter.res))
++  test-pcg-converges  ^-  tang
  =/  res  (pcg:sad spd2 (vec ~[.~10 .~9]) 20)
  (expect !>((close (vec ~[.~2 .~1]) x.res)))
++  test-pcg-badly-scaled  ^-  tang
  ::  diag(1e6, 1) is where preconditioning earns its keep; both must converge
  =/  a  (mat 2 2 ~[~[.~1e6 .~0] ~[.~0 .~1]])
  =/  v  (vec ~[.~1e6 .~1])
  (expect !>((close (vec ~[.~1 .~1]) x:(pcg:sad a v 50))))
::
::  +svd -- exact cases (columns already orthogonal, so no rotation is needed)
::
++  test-svd-diagonal  ^-  tang
  =/  res  (svd:sad dg24)
  ::  s = [4 2] descending, so the columns swap
  %+  weld  (expect-eq !>((vec ~[.~4 .~2])) !>(s.res))
  %+  weld  (expect-eq !>((mat 2 2 ~[~[.~0 .~1] ~[.~1 .~0]])) !>(u.res))
  (expect-eq !>((mat 2 2 ~[~[.~0 .~1] ~[.~1 .~0]])) !>(v.res))
++  test-svd-sorted-identity  ^-  tang
  ::  already descending: U = V = I
  =/  res  (svd:sad (mat 2 2 ~[~[.~2 .~0] ~[.~0 .~1]]))
  %+  weld  (expect-eq !>((vec ~[.~2 .~1])) !>(s.res))
  %+  weld  (expect-eq !>(eye2) !>(u.res))
  (expect-eq !>(eye2) !>(v.res))
++  test-svd-equal-values  ^-  tang
  ::  [[3 4] [4 -3]] has orthogonal columns of equal norm 5
  (expect-eq !>((vec ~[.~5 .~5])) !>(s:(svd:sad (mat 2 2 ~[~[.~3 .~4] ~[.~4 .~-3]]))))
++  test-svd-rectangular  ^-  tang
  ::  3x2 thin SVD: U is 3x2, V is 2x2
  =/  res  (svd:sad (mat 3 2 ~[~[.~3 .~0] ~[.~0 .~4] ~[.~0 .~0]]))
  %+  weld  (expect-eq !>((vec ~[.~4 .~3])) !>(s.res))
  %+  weld  (expect-eq !>((mat 3 2 ~[~[.~0 .~1] ~[.~1 .~0] ~[.~0 .~0]])) !>(u.res))
  (expect-eq !>((mat 2 2 ~[~[.~0 .~1] ~[.~1 .~0]])) !>(v.res))
++  test-svd-zero-column  ^-  tang
  ::  a zero singular value leaves that column of U zero rather than dividing
  =/  res  (svd:sad (mat 2 2 ~[~[.~3 .~0] ~[.~0 .~0]]))
  %+  weld  (expect-eq !>((vec ~[.~3 .~0])) !>(s.res))
  (expect-eq !>((mat 2 2 ~[~[.~1 .~0] ~[.~0 .~0]])) !>(u.res))
++  test-svd-vals  ^-  tang
  (expect-eq !>((vec ~[.~4 .~2])) !>((svd-vals:sad dg24)))
++  test-svd-crashes-wide  ^-  tang
  (expect-fail |.((svd:sad (mat 2 3 ~[~[.~1 .~2 .~3] ~[.~4 .~5 .~6]]))))
::
::  +svd -- approximate cases against numpy.linalg.svd
::
++  test-svd-values-2x2  ^-  tang
  =/  want  (vec ~[.~5.464985704219043 .~0.3659661906262575])
  (expect !>((close want (svd-vals:sad gen22))))
++  test-svd-values-4x3  ^-  tang
  =/  want  (vec ~[.~17.488318893441512 .~1.6812882949946584 .~0.5761700707347798])
  (expect !>((close want (svd-vals:sad gen43))))
++  test-svd-reconstructs-2x2  ^-  tang
  =/  res  (svd:sad gen22)
  (expect !>((close gen22 (mmul:lk (mmul:lk u.res (mk-diag s.res 2)) (tpose v.res)))))
++  test-svd-reconstructs-4x3  ^-  tang
  =/  res  (svd:sad gen43)
  (expect !>((close gen43 (mmul:lk (mmul:lk u.res (mk-diag s.res 3)) (tpose v.res)))))
++  test-svd-v-orthogonal  ^-  tang
  =/  res  (svd:sad gen43)
  (expect !>((close (eye:lk [~[3 3] 6 %i754 ~]) (mmul:lk (tpose v.res) v.res))))
++  test-svd-u-orthonormal  ^-  tang
  ::  thin U has orthonormal columns: U^T U = I
  =/  res  (svd:sad gen43)
  (expect !>((close (eye:lk [~[3 3] 6 %i754 ~]) (mmul:lk (tpose u.res) u.res))))
++  test-svd-descending  ^-  tang
  ::  singular values come out sorted
  =/  s  (svd-vals:sad gen43)
  =/  a  `@rd`(get-item:lk s ~[0])
  =/  b  `@rd`(get-item:lk s ~[1])
  =/  c  `@rd`(get-item:lk s ~[2])
  =/  ge  ~(gte rd:math [%n .~1e-12 `@rd`0])
  (expect !>(?&((ge a b) (ge b c))))
--
