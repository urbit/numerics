/-  ls=lagoon
/+  *lagoon
|^
~
::  conjugate gradient descent method
::    a  is a 2D system of linear equations
::    b  is the corresponding solutions
::    x0 is a starting guess at the unknowns
::
++  conjgrad
  |=  [a=ray:ls b=ray:ls x0=ray:ls]
  ^-  ray:ls
  =/  r=ray:ls  (sub b (mmul a x0))
  =/  p=ray:ls  r
  =/  rsold=ray:ls  (dot r r)
  ::
  =/  i=@ud  0
  =/  x=ray:ls  x0
  |-  ^-  ray
  ?:  =(i (&1.shape.meta.b))  x
  =/  ap=ray  (mmul a p)
  =/  alpha=ray  (div rsold (dot p ap))
  =/  x  (add x (mul (fill meta.p alpha) p))
  =/  r  (sub r (mul (fill meta.ap alpha) ap))
  =/  rsnew=ray  (dot r r)
  ?:  (lth rsnew (scalarize meta.rsnew %rs .1e-5))  x
  %=  $
    i      +(i)
    p      (add r (mul (fill meta.p (div rsnew rsold)) p))
    rsold  rsnew
  ==
--

