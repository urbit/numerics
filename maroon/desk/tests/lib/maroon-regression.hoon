/-  ls=lagoon
/+  *test, *lagoon, maroon, math
::::  /tests/lib/maroon-regression -- linear and ridge regression, metrics
::
::  EXACT cases use x = [0 0 2 2], y = 3x + 1: the centred feature [-1 -1 1 1]
::  has norm exactly 2, so the Householder step, the back-substitution, and
::  the ridge Cholesky (with alpha = 12, making the system [[16]]) all stay
::  exact.  Verified in exact rational arithmetic by maroon/tools/ml_check.py,
::  which also checks every value against scikit-learn.
::
::  APPROXIMATE cases fit two features against sklearn's LinearRegression and
::  Ridge.  Intercepts are raw %i754 scalars, compared as @rd.
::
|%
++  mm  (make:maroon %n .~1e-12)
++  lk  (lake %n)
++  close
  |=  [x=ray:ls y=ray:ls]
  ^-  ?
  (all:lk (is-close:lk x y [.~1e-9 .~1e-9]))
++  near
  |=  [x=@ y=@]
  ^-  ?
  (close (en-ray:lk [[~[1] 6 %i754 ~] ~[x]]) (en-ray:lk [[~[1] 6 %i754 ~] ~[y]]))
++  mat
  |=  [r=@ c=@ v=(list (list @))]
  ^-  ray:ls
  (en-ray:lk [[~[r c] 6 %i754 ~] v])
++  vec
  |=  v=(list @)
  ^-  ray:ls
  (en-ray:lk [[~[(lent v)] 6 %i754 ~] v])
::  fixtures
++  rx  (mat 4 1 ~[~[.~0] ~[.~0] ~[.~2] ~[.~2]])
++  ry  (vec ~[.~1 .~1 .~7 .~7])
++  mx  (mat 5 2 ~[~[.~1 .~2] ~[.~2 .~1] ~[.~3 .~4] ~[.~4 .~3] ~[.~5 .~6]])
++  my  (vec ~[.~3.1 .~2.9 .~7.2 .~6.8 .~11.1])
::
::  +linreg
::
++  test-linreg-coef  ^-  tang
  (expect-eq !>((vec ~[.~3])) !>(coef:(linreg:mm rx ry)))
++  test-linreg-intercept  ^-  tang
  (expect-eq !>(`@rd`.~1) !>(`@rd`intercept:(linreg:mm rx ry)))
++  test-linreg-perfect-fit  ^-  tang
  ::  y is exactly linear in x, so the residual vanishes
  =/  f  (linreg:mm rx ry)
  =/  p  (predict:mm rx coef.f intercept.f)
  %+  weld  (expect-eq !>(ry) !>(p))
  %+  weld  (expect-eq !>(`@rd`.~0) !>(`@rd`(mse:mm ry p)))
  (expect-eq !>(`@rd`.~1) !>(`@rd`(r2:mm ry p)))
++  test-linreg-two-features  ^-  tang
  ::  sklearn LinearRegression: coef [0.84833.. 1.14166..], intercept 0.021666..
  =/  f  (linreg:mm mx my)
  %+  weld  (expect !>((close (vec ~[.~0.8483333333333334 .~1.1416666666666668]) coef.f)))
  (expect !>((near .~0.0216666666666665 intercept.f)))
++  test-linreg-two-features-r2  ^-  tang
  =/  f  (linreg:mm mx my)
  (expect !>((near .~0.9997674486206797 (r2:mm my (predict:mm mx coef.f intercept.f)))))
++  test-linreg-crashes-mismatch  ^-  tang
  (expect-fail |.((linreg:mm rx (vec ~[.~1 .~2]))))
::
::  +ridge
::
++  test-ridge-coef  ^-  tang
  ::  (4 + 12) coef = 12, so coef = 0.75
  (expect-eq !>((vec ~[.~0.75])) !>(coef:(ridge:mm rx ry .~12)))
++  test-ridge-intercept  ^-  tang
  (expect-eq !>(`@rd`.~3.25) !>(`@rd`intercept:(ridge:mm rx ry .~12)))
++  test-ridge-shrinks  ^-  tang
  ::  the penalty pulls the slope toward zero: 0.75 < 3
  =/  lt  ~(lth rd:math [%n .~1e-12 `@rd`0])
  (expect !>((lt `@rd`(get-item:lk coef:(ridge:mm rx ry .~12) ~[0]) .~3)))
++  test-ridge-alpha-zero-is-ols  ^-  tang
  ::  alpha = 0 reduces to ordinary least squares
  (expect-eq !>(coef:(linreg:mm rx ry)) !>(coef:(ridge:mm rx ry .~0)))
++  test-ridge-mse  ^-  tang
  =/  f  (ridge:mm rx ry .~12)
  (expect-eq !>(`@rd`.~5.0625) !>(`@rd`(mse:mm ry (predict:mm rx coef.f intercept.f))))
++  test-ridge-r2  ^-  tang
  ::  1 - 20.25/36 = 0.4375, exactly
  =/  f  (ridge:mm rx ry .~12)
  (expect-eq !>(`@rd`.~0.4375) !>(`@rd`(r2:mm ry (predict:mm rx coef.f intercept.f))))
++  test-ridge-two-features  ^-  tang
  ::  sklearn Ridge(alpha=1.0)
  =/  f  (ridge:mm mx my .~1)
  %+  weld  (expect !>((close (vec ~[.~0.8214092140921403 .~1.0864498644986456]) coef.f)))
  (expect !>((near .~0.279132791327914 intercept.f)))
++  test-ridge-crashes-negative-alpha  ^-  tang
  (expect-fail |.((ridge:mm rx ry .~-1)))
::
::  +predict / +mse / +r2
::
++  test-predict  ^-  tang
  (expect-eq !>(ry) !>((predict:mm rx (vec ~[.~3]) .~1)))
++  test-r2-mean-predictor  ^-  tang
  ::  predicting the mean everywhere scores exactly 0
  (expect-eq !>(`@rd`.~0) !>(`@rd`(r2:mm ry (vec ~[.~4 .~4 .~4 .~4]))))
++  test-r2-crashes-constant  ^-  tang
  (expect-fail |.((r2:mm (vec ~[.~2 .~2 .~2]) (vec ~[.~1 .~2 .~3]))))
--
