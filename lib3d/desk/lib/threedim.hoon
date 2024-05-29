  ::  /lib/threedim
::::  Version ~2024.5.29 by ~lagrev-nocfep
::
::
/-  ls=lagoon,
    ts=threedim
/+  *lagoon,
    *saloon
|%
::  Render 3D points in 2D space.
++  render-3d-points
  |=  [a=ray:ls c=ray:ls theta=ray:ls e=ray:ls]
  ^-  ray:ls
  =/  asin  (sin:sa theta)
  =/  acos  (cos:sa theta)
  =/  xyz   (sub:la a c)
  =/  d  (zeros:la meta.a)
  =.  d  (set-col:la d ~[0] (sub:la (mul:la (get-col:la acos ~[1]) (add:la (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[0])))) (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[2]))))
  =.  d  (set-col:la d ~[1] (add:la (mul:la (get-col:la asin ~[0]) (add:la (mul:la (get-col:la acos ~[1]) (get-col:la xyz ~[2])) (mul:la (get-col:la asin ~[1]) (add:la (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[0])))))) (mul:la (get-col:la acos ~[0]) (add:la (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[0]))))))
  =.  d  (set-col:la d ~[2] (sub:la (mul:la (get-col:la acos ~[0]) (add:la (mul:la (get-col:la acos ~[1]) (get-col:la xyz ~[2])) (mul:la (get-col:la asin ~[1]) (add:la (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[0])))))) (mul:la (get-col:la asin ~[0]) (add:la (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[0]))))))
  =/  b  (zeros:la [~[(snag 0 shape.meta.a) 2] 3 %uint ~])
  =.  b  (set-col:la b ~[0] (add:la (mul:la (get-col:la e ~[2]) (div:la (get-col:la d ~[0]) (get-col:la d ~[2]))) (get-col:la e ~[0])))
  (set-col:la b ~[1] (add:la (mul:la (get-col:la e ~[2]) (div:la (get-col:la d ~[1]) (get-col:la d ~[2]))) (get-col:la e ~[1])))
::  Produce a grid of points for the canvas with rendered points indicated.
++  gen-3d
  |=  [ipts=(list (list @rs)) =view:ts]
  ^-  ray:ls
  =/  m  (lent ipts)
  =/  n  (lent (snag 0 ipts))
  ?>  =(3 n)
  =/  dx  (div:rs (sub:rs xmax.view xmin.view) (sun:rs 100))
  =/  dy  (div:rs (sub:rs ymax.view ymin.view) (sun:rs 100))
  =/  rpts   (en-ray:la [~[m n] 3 %uint ~] ipts)
  =/  c      (en-ray:la [~[m n] 3 %i754 ~] (reap m ~[.114 .-83 .100]))
  =/  theta  (en-ray:la [~[m n] 3 %i754 ~] (reap m ~[.0.0 .0.0 .0.0]))
  =/  e      (en-ray:la [~[m n] 3 %i754 ~] (reap m ~[.10.0 .10.0 .10.0]))
  =/  pts    (render-3d-points rpts c theta e)
  ~&  pts
  =/  grid   (zeros:la [~[m n] 3 %uint ~])
  =|  idx=@
  |-
  ?:  =(m idx)  grid
  =/  pt=(list @rs)  (snag idx ipts)
  =/  sx=@rs  (snag 0 pt)
  =/  sy=@rs  (snag 1 pt)
  ::  Locate the pi-th bin from the floating-point coordinates of the point.
  =/  px=@ud  (bin sx xmin.view dx m)
  =/  py=@ud  (bin sy ymin.view dy m)
  $(idx +(idx), grid (set-item:la grid ~[px py] .1))
++  bin
  ::  Bin the value.
  |=  [p=@rs min=@rs ds=@rs np=@ud]
  =/  idx=@ud  0
  |-  ^-  @ud
    ::  The naive algorithm is to count up in range until we find the right one.
    =/  current-bin  (add:rs min (mul:rs ds (sun:rs idx)))
    ?:  (lth:rs p current-bin)  idx
  $(idx +(idx))
--
