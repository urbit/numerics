/-  ls=lagoon
/+  *test, *lagoon, maroon, math
::::  /tests/lib/maroon-pca -- centering, covariance, and PCA
::
::  EXACT cases compare whole rays: the inputs are small integers whose means,
::  centred values and Gram entries are all exactly representable, verified
::  op-by-op in exact rational arithmetic by maroon/tools/ml_check.py.
::
::  APPROXIMATE cases compare within a tolerance against scikit-learn.  The
::  explained VARIANCE is always approximate even when the data is exact,
::  because it comes back through Saloon's Newton-iteration square root: for
::  data whose true variance is 1 the ship returns 0.9999999999999998.
::
::  Component signs are arbitrary in any PCA, so sign-sensitive checks compare
::  absolute values.
::
|%
++  mm  (make:maroon %n .~1e-12)
++  lk  (lake %n)
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
::    +x-axis:  three points along the x axis; the first component is e1
++  x-axis  (mat 3 2 ~[~[.~1 .~0] ~[.~2 .~0] ~[.~3 .~0]])
::    +box:  a 4-point rectangle centred at (2,1)
++  box  (mat 4 2 ~[~[.~0 .~0] ~[.~4 .~0] ~[.~0 .~2] ~[.~4 .~2]])
::    +diag-line:  four points on the 45-degree line through (5,4)
++  diag-line  (mat 4 2 ~[~[.~2 .~1] ~[.~4 .~3] ~[.~6 .~5] ~[.~8 .~7]])
::
::  +col-mean / +center
::
++  test-col-mean  ^-  tang
  (expect-eq !>((vec ~[.~2 .~0])) !>((col-mean:mm x-axis)))
++  test-col-mean-box  ^-  tang
  (expect-eq !>((vec ~[.~2 .~1])) !>((col-mean:mm box)))
++  test-center  ^-  tang
  %+  expect-eq  !>((mat 3 2 ~[~[.~-1 .~0] ~[.~0 .~0] ~[.~1 .~0]]))
  !>((center:mm x-axis))
++  test-center-box  ^-  tang
  %+  expect-eq
    !>((mat 4 2 ~[~[.~-2 .~-1] ~[.~2 .~-1] ~[.~-2 .~1] ~[.~2 .~1]]))
  !>((center:mm box))
++  test-center-sums-to-zero  ^-  tang
  ::  a centred dataset has zero column means
  (expect !>((close (vec ~[.~0 .~0]) (col-mean:mm (center:mm diag-line)))))
::
::  +cov
::
++  test-cov  ^-  tang
  ::  the box's centred columns are orthogonal: Xc^T*Xc/n = diag(4,1), exactly
  (expect-eq !>((mat 2 2 ~[~[.~4 .~0] ~[.~0 .~1]])) !>((cov:mm box 0)))
++  test-cov-symmetric  ^-  tang
  ::  a covariance matrix is symmetric, which is what makes it +eig-able
  =/  c  (cov:mm diag-line 1)
  (expect !>(=((get-item:lk c ~[0 1]) (get-item:lk c ~[1 0]))))
++  test-cov-crashes-ddof  ^-  tang
  (expect-fail |.((cov:mm (mat 1 2 ~[~[.~1 .~2]]) 1)))
::
::  +pca
::
++  test-pca-component  ^-  tang
  ::  data on the x axis: the first component is e1, exactly (the centred
  ::  columns are already orthogonal, so Jacobi performs no rotation)
  (expect-eq !>((mat 2 1 ~[~[.~1] ~[.~0]])) !>(comp:(pca:mm x-axis 1)))
++  test-pca-mean  ^-  tang
  (expect-eq !>((vec ~[.~2 .~0])) !>(mean:(pca:mm x-axis 1)))
++  test-pca-variance  ^-  tang
  ::  variance 1, to within the sqrt-then-square round trip
  (expect !>((close (vec ~[.~1]) vals:(pca:mm x-axis 1))))
++  test-pca-transform  ^-  tang
  ::  scores of the training data on e1: the centred x coordinates
  =/  p  (pca:mm x-axis 1)
  %+  expect-eq  !>((mat 3 1 ~[~[.~-1] ~[.~0] ~[.~1]]))
  !>((pca-transform:mm x-axis comp.p mean.p))
++  test-pca-two-components  ^-  tang
  ::  asking for both components of the box gives the two axes, largest first
  =/  p  (pca:mm box 2)
  =/  want  (mat 2 2 ~[~[.~1 .~0] ~[.~0 .~1]])
  (expect !>((close want (abs:lk comp.p))))
++  test-pca-rotated  ^-  tang
  ::  points on the 45-degree line: the component is (1,1)/sqrt(2), sign-free
  =/  p  (pca:mm diag-line 1)
  =/  want  (mat 2 1 ~[~[.~0.7071067811865475] ~[.~0.7071067811865475]])
  (expect !>((close want (abs:lk comp.p))))
++  test-pca-rotated-variance  ^-  tang
  ::  sklearn: explained_variance_ = 13.333333333333336
  (expect !>((close (vec ~[.~13.333333333333336]) vals:(pca:mm diag-line 1))))
++  test-pca-descending  ^-  tang
  ::  variances come out largest first
  =/  v  vals:(pca:mm box 2)
  =/  ge  ~(gte rd:math [%n .~1e-12 `@rd`0])
  (expect !>((ge `@rd`(get-item:lk v ~[0]) `@rd`(get-item:lk v ~[1]))))
++  test-pca-crashes-k-zero  ^-  tang
  (expect-fail |.((pca:mm x-axis 0)))
++  test-pca-crashes-k-too-big  ^-  tang
  (expect-fail |.((pca:mm x-axis 3)))
++  test-pca-crashes-one-sample  ^-  tang
  (expect-fail |.((pca:mm (mat 1 2 ~[~[.~1 .~2]]) 1)))
--
