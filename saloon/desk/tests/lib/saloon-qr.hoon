/-  ls=lagoon
/+  *test, *saloon, *lagoon
::::  /tests/lib/saloon-qr -- Householder QR, least squares, and their helpers
::
::  EXACT cases compare whole rays: inputs whose column norms are powers of two
::  (or zero), so every Householder vector, scale and reflection is exactly
::  representable.  Verified op-by-op in exact rational arithmetic by
::  saloon/tools/linalg_check.py.
::
::  APPROXIMATE cases compare within a tolerance against numpy.linalg.qr and
::  numpy.linalg.lstsq.  Signs are compared too, not just magnitudes: +house
::  uses LAPACK's convention (alpha = -sign(x0)*|x|), and the oracle confirms
::  the factors match numpy's including signs.
::
::  Q^T*Q is checked with +gram rather than a transpose, which also exercises
::  +gram and avoids the lagoon transpose jet (broken before urbit/vere#1057).
::
|%
++  sad  (sake %n .~1e-12)
++  lk   (lake %n)
++  close
  |=  [x=ray:ls y=ray:ls]
  ^-  ?
  (all:lk (is-close:lk x y [.~1e-9 .~1e-9]))
++  mat
  |=  [r=@ c=@ v=(list (list @))]
  ^-  ray:ls
  (en-ray:lk [[~[r c] 6 %i754 ~] v])
++  vec
  |=  v=(list @)
  ^-  ray:ls
  (en-ray:lk [[~[(lent v)] 6 %i754 ~] v])
::  fixtures
++  diag2  (mat 2 2 ~[~[.~2 .~0] ~[.~0 .~2]])
++  tall   (mat 3 2 ~[~[.~2 .~0] ~[.~0 .~2] ~[.~0 .~0]])
++  gen    (mat 3 2 ~[~[.~1 .~2] ~[.~3 .~4] ~[.~5 .~7]])
++  m32    (mat 3 2 ~[~[.~1 .~2] ~[.~3 .~4] ~[.~5 .~6]])
::
::  +matvec-t / +gram
::
++  test-matvec-t  ^-  tang
  ::  M^T*[1 1 1] is the column sums
  (expect-eq !>((vec ~[.~9 .~12])) !>((matvec-t:sad m32 (vec ~[.~1 .~1 .~1]))))
++  test-gram  ^-  tang
  (expect-eq !>((mat 2 2 ~[~[.~35 .~44] ~[.~44 .~56]])) !>((gram:sad m32)))
++  test-gram-symmetric  ^-  tang
  ::  exactly symmetric, since (i,j) and (j,i) multiply the same pairs
  =/  g  (gram:sad gen)
  (expect !>(=((get-item:lk g ~[0 1]) (get-item:lk g ~[1 0]))))
::
::  +qr -- exact
::
++  test-qr-diagonal  ^-  tang
  ::  LAPACK signs: Q = -I, R = -2I
  =/  f  (qr:sad diag2)
  %+  weld  (expect-eq !>((mat 2 2 ~[~[.~-1 .~0] ~[.~0 .~-1]])) !>(q.f))
  (expect-eq !>((mat 2 2 ~[~[.~-2 .~0] ~[.~0 .~-2]])) !>(r.f))
++  test-qr-tall  ^-  tang
  ::  thin QR of a 3x2: Q is 3x2, R is 2x2
  =/  f  (qr:sad tall)
  %+  weld  (expect-eq !>((mat 3 2 ~[~[.~-1 .~0] ~[.~0 .~-1] ~[.~0 .~0]])) !>(q.f))
  (expect-eq !>((mat 2 2 ~[~[.~-2 .~0] ~[.~0 .~-2]])) !>(r.f))
++  test-qr-single-column  ^-  tang
  ::  [3 4]^T has norm 5, so R is exactly -5 (Q holds 0.6/0.8, which are not)
  (expect-eq !>((mat 1 1 ~[~[.~-5]])) !>(r:(qr:sad (mat 2 1 ~[~[.~3] ~[.~4]]))))
++  test-qr-zero-column  ^-  tang
  ::  a zero first column needs no reflection and is skipped, not divided by
  =/  f  (qr:sad (mat 2 2 ~[~[.~0 .~1] ~[.~0 .~1]]))
  %+  weld  (expect-eq !>((mat 2 2 ~[~[.~1 .~0] ~[.~0 .~-1]])) !>(q.f))
  (expect-eq !>((mat 2 2 ~[~[.~0 .~1] ~[.~0 .~-1]])) !>(r.f))
++  test-qr-subdiagonal-exactly-zero  ^-  tang
  ::  R is upper triangular with EXACT zeros below the diagonal, not residue
  (expect !>(=(0 (get-item:lk r:(qr:sad gen) ~[1 0]))))
::
::  +qr -- approximate, against numpy (signs included)
::
++  test-qr-r-matches-numpy  ^-  tang
  =/  want  (mat 2 2 ~[~[.~-5.916079783099615 .~-8.282511696339462] ~[.~0 .~0.6324555320336759]])
  (expect !>((close want r:(qr:sad gen))))
++  test-qr-q-matches-numpy  ^-  tang
  =/  want
    %^  mat  3  2
    :~  ~[.~-0.16903085094570325 .~0.9486832980505133]
        ~[.~-0.50709255283711 .~-0.31622776601683833]
        ~[.~-0.8451542547285166 .~2.3208668709370606e-16]
    ==
  (expect !>((close want q:(qr:sad gen))))
++  test-qr-reconstructs  ^-  tang
  =/  f  (qr:sad gen)
  (expect !>((close gen (mmul:lk q.f r.f))))
++  test-qr-q-orthonormal  ^-  tang
  (expect !>((close (eye:lk [~[2 2] 6 %i754 ~]) (gram:sad q:(qr:sad gen)))))
++  test-qr-crashes-wide  ^-  tang
  (expect-fail |.((qr:sad (mat 2 3 ~[~[.~1 .~2 .~3] ~[.~4 .~5 .~6]]))))
::
::  +trsv-r / +lstsq
::
++  test-trsv-r  ^-  tang
  ::  [[2 1] [0 4]] x = [5 8]  ->  x = [1.5 2]
  (expect-eq !>((vec ~[.~1.5 .~2])) !>((trsv-r:sad (mat 2 2 ~[~[.~2 .~1] ~[.~0 .~4]]) (vec ~[.~5 .~8]))))
++  test-lstsq-tall  ^-  tang
  ::  the third row of v is unreachable, so the best fit ignores it
  (expect-eq !>((vec ~[.~2 .~3])) !>((lstsq:sad tall (vec ~[.~4 .~6 .~1]))))
++  test-lstsq-square  ^-  tang
  ::  on a square full-rank system least squares is just the solve
  (expect-eq !>((vec ~[.~2 .~3])) !>((lstsq:sad diag2 (vec ~[.~4 .~6]))))
++  test-lstsq-matches-numpy  ^-  tang
  ::  numpy.linalg.lstsq: [-3/14, 1/2] to rounding
  =/  want  (vec ~[.~-0.2142857142857133 .~0.4999999999999992])
  (expect !>((close want (lstsq:sad gen (vec ~[.~1 .~2 .~2])))))
++  test-lstsq-matches-chol-normal-equations  ^-  tang
  ::  QR and the normal equations agree on a well-conditioned problem
  =/  v  (vec ~[.~1 .~2 .~2])
  =/  ne  (chol-solve:sad (gram:sad gen) (matvec-t:sad gen v))
  (expect !>((close ne (lstsq:sad gen v))))
++  test-lstsq-crashes-length  ^-  tang
  (expect-fail |.((lstsq:sad tall (vec ~[.~1 .~2]))))
--
