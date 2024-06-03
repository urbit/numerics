  ::  /lib/threedim
::::  Version ~2024.5.29 by ~lagrev-nocfep
::
::
/-  ls=lagoon,
    ts=threedim
/+  *lagoon,
    math,
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
  =.  d  (set-col:la d ~[0] (sub:la (mul:la (get-col:la acos ~[1]) (add:la (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[0])))) (mul:la (get-col:la asin ~[1]) (get-col:la xyz ~[2]))))
  =.  d  (set-col:la d ~[1] (add:la (mul:la (get-col:la asin ~[0]) (add:la (mul:la (get-col:la acos ~[1]) (get-col:la xyz ~[2])) (mul:la (get-col:la asin ~[1]) (add:la (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[0])))))) (mul:la (get-col:la acos ~[0]) (add:la (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[0]))))))
  =.  d  (set-col:la d ~[2] (sub:la (mul:la (get-col:la acos ~[0]) (add:la (mul:la (get-col:la acos ~[1]) (get-col:la xyz ~[2])) (mul:la (get-col:la asin ~[1]) (add:la (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[0])))))) (mul:la (get-col:la asin ~[0]) (add:la (mul:la (get-col:la acos ~[2]) (get-col:la xyz ~[1])) (mul:la (get-col:la asin ~[2]) (get-col:la xyz ~[0]))))))
  =/  b  (zeros:la [~[(snag 0 shape.meta.a) 2] 5 %i754 ~])
  =.  b  (set-col:la b ~[0] (add:la (mul:la (get-col:la e ~[2]) (div:la (get-col:la d ~[0]) (get-col:la d ~[2]))) (get-col:la e ~[0])))
  =.  b  (set-col:la b ~[1] (add:la (mul:la (get-col:la e ~[2]) (div:la (get-col:la d ~[1]) (get-col:la d ~[2]))) (get-col:la e ~[1])))
  b
::  Produce a grid of points for the canvas with rendered points indicated.
++  gen-3d
  |=  [ipts=(list (list @rs)) =view:ts =lens:ts]
  ^-  ray:ls
  =/  m  (lent ipts)
  =/  n  (lent (snag 0 ipts))
  ?>  =(3 n)
  ::  Project 3D Points to 2D: Transform the 3D points to 2D coordinates on the viewport.
  ::  Calculate Depth (Distance): Compute the distance of each point from the camera to determine the shading.
  ::  Rasterize Points: Place the projected points on the 2D grid and apply shading based on depth.
  =/  dx  (div:rs (sub:rs xmax.view xmin.view) (sun:rs 100))
  =/  dy  (div:rs (sub:rs ymax.view ymin.view) (sun:rs 100))
  =/  rpts   (change:la (en-ray:la [~[m n] 3 %uint ~] ipts) %i754 5)
  =/  c      (en-ray:la [~[m n] 5 %i754 ~] (reap m ~[.10 .-8 .10]))
  =/  theta  (en-ray:la [~[m n] 5 %i754 ~] (reap m ~[.0.0 .0.0 .0.0]))
  =/  e      (en-ray:la [~[m n] 5 %i754 ~] (reap m ~[.3.0 .3.0 .3.0]))
  =/  pts    (render-3d-points rpts c theta e)
  =/  grid   (zeros:la [~[nx.view ny.view] 5 %i754 ~])
  =|  idx=@
  |-
  ?:  =(m idx)  grid
  =/  pt=(list @rs)  (snag idx ipts)
  =/  sx=@rs  (snag 0 pt)
  =/  sy=@rs  (snag 1 pt)
  ::  Locate the p-th bin from the floating-point coordinates of the point.
  =/  px=@ud  (min (dec nx.view) (bin sx xmin.view dx m))
  =/  py=@ud  (min (dec ny.view) (bin sy ymin.view dy m))
  ~&  >>>  [sx sy px py]
  $(idx +(idx), grid (set-item:la grid ~[px py] .255))
++  bin
  ::  Bin the value.
  |=  [p=@rs min=@rs ds=@rs np=@ud]
  ^-  @ud
  =|  idx=@ud
  |-
  ::  The naive algorithm is to count up in range until we find the right one.
  =/  current-bin  (add:rs min (mul:rs ds (sun:rs idx)))
  ?:  (lth:rs p current-bin)  idx
  $(idx +(idx))
::  +triplicate
::
::  Convert single points into triples of the point.
::  For this case, we are working with @rs values only.
++  triplicate
  |=  =ray:ls
  ^-  ray:ls
  =,  meta.ray
  =/  pts  (ravel:la ray)
  =/  mn  (lent pts)
  =|  pxl=(list @rs)
  |-
  ?~  pts  (en-ray:la [~[(mul 3 mn)] +.meta.ray] pxl)
  $(pts t.pts, pxl [i.pts i.pts i.pts pxl])
::
++  look-at-matrix
  |=  =lens:ts
  ^-  ray:ls
  =/  c-posn  (en-ray:la [~[3 1] 5 %i754 ~] ~[~[x.lens] ~[y.lens] ~[z.lens]])
  =/  t-posn  (en-ray:la [~[3 1] 5 %i754 ~] ~[~[tx.lens] ~[ty.lens] ~[tz.lens]])
  =/  u  (en-ray:la [~[3 1] 5 %i754 ~] ~[~[.0.0] ~[.1.0] ~[.0.0]])
  ::  Create the look-at view matrix.
  =/  f  (sub:la c-posn t-posn)
  =.  f  (div-scalar:la f (get-item:la (sqrt:sa (dot:la f f)) ~[0 0]))
  =/  r1  `@rs`(sub:rs (mul:rs `@rs`(get-item:la u ~[1 0]) `@rs`(get-item:la f ~[2 0])) (mul:rs `@rs`(get-item:la u ~[2 0]) `@rs`(get-item:la f ~[1 0])))
  =/  r2  `@rs`(sub:rs (mul:rs `@rs`(get-item:la u ~[2 0]) `@rs`(get-item:la f ~[0 0])) (mul:rs `@rs`(get-item:la u ~[0 0]) `@rs`(get-item:la f ~[2 0])))
  =/  r3  `@rs`(sub:rs (mul:rs `@rs`(get-item:la u ~[0 0]) `@rs`(get-item:la f ~[1 0])) (mul:rs `@rs`(get-item:la u ~[1 0]) `@rs`(get-item:la f ~[0 0])))
  =/  r  (en-ray:la [~[3 1] 5 %i754 ~] ~[~[r1] ~[r2] ~[r3]])
  =.  r  (div-scalar:la r (get-item:la (sqrt:sa (dot:la r r)) ~[0 0]))
  =/  u1  (sub:rs (mul:rs `@rs`(get-item:la f ~[1 0]) `@rs`(get-item:la r ~[2 0])) (mul:rs `@rs`(get-item:la f ~[2 0]) `@rs`(get-item:la r ~[1 0])))
  =/  u2  (sub:rs (mul:rs `@rs`(get-item:la f ~[2 0]) `@rs`(get-item:la r ~[0 0])) (mul:rs `@rs`(get-item:la f ~[0 0]) `@rs`(get-item:la r ~[2 0])))
  =/  u3  (sub:rs (mul:rs `@rs`(get-item:la f ~[0 0]) `@rs`(get-item:la r ~[1 0])) (mul:rs `@rs`(get-item:la f ~[1 0]) `@rs`(get-item:la r ~[0 0])))
  =/  u  (en-ray:la [~[3 1] 5 %i754 ~] ~[~[u1] ~[u2] ~[u3]])
  =/  m  (eye:la [~[4 4] 5 %i754 ~])
  =.  m  (set-item:la m ~[0 0] (get-item:la r ~[0 0]))
  =.  m  (set-item:la m ~[0 1] (get-item:la r ~[1 0]))
  =.  m  (set-item:la m ~[0 2] (get-item:la r ~[2 0]))
  =.  m  (set-item:la m ~[1 0] (get-item:la u ~[0 0]))
  =.  m  (set-item:la m ~[1 1] (get-item:la u ~[1 0]))
  =.  m  (set-item:la m ~[1 2] (get-item:la u ~[2 0]))
  =.  m  (set-item:la m ~[2 0] (get-item:la f ~[0 0]))
  =.  m  (set-item:la m ~[2 1] (get-item:la f ~[1 0]))
  =.  m  (set-item:la m ~[2 2] (get-item:la f ~[2 0]))
  ::  status:  need to fix +submatrix
  =/  mr  (mmul:la (sub:la (zeros:la [~[2 2] 5 %i754 ~]) (submatrix:la ~[`[[`0 `0]] `[[`2 `2]]] m)) c-posn)
  =.  m  (set-item:la m ~[3 0] (get-item:la mr ~[0 0]))
  =.  m  (set-item:la m ~[3 1] (get-item:la mr ~[1 0]))
  =.  m  (set-item:la m ~[3 2] (get-item:la mr ~[2 0]))
  =.  m  (set-item:la m ~[3 3] (get-item:la mr ~[3 0]))
  m
::
++  project-points
  |=  [pts=ray:ls =view:ts =lens:ts]
  ^-  ray:ls
  =/  fov  .60  ::  field of view, degrees
  =/  asp  (div:rs (sun:rs nx.view) (sun:rs ny.view))  ::  aspect ratio
  =/  near  .0.1  ::  near clipping plane
  =/  far  .1000  ::  far clipping plane
  ::  Compute the projection matrix.
  =/  f  (div:rs .1 (tan:rs:math (mul:rs (deg2rad fov) .0.5)))
  =/  m  (zeros:la [~[4 4] 5 %i754 ~])
  =.  m  (set-item:la m ~[0 0] (div:rs f asp))
  =.  m  (set-item:la m ~[1 1] f)
  =.  m  (set-item:la m ~[2 2] (div:rs (add:rs near far) (sub:rs near far)))
  =.  m  (set-item:la m ~[2 3] (mul:rs :(mul:rs .2 near far) (div:rs near far)))
  =.  m  (set-item:la m ~[3 2] -1)
  ::  Compute the view matrix.
  =/  l  (look-at-matrix lens)
  ::  Project the points to camera space.
  =/  ph  (ones:la [~[(snag 0 shape.meta.pts) 4] 5 %i754 ~])
  =.  ph  (set-col:la ph ~[0] (get-col:la pts ~[0]))
  =.  ph  (set-col:la ph ~[1] (get-col:la pts ~[1]))
  =.  ph  (set-col:la ph ~[2] (get-col:la pts ~[2]))
  =/  pc  (mmul:la (mmul:la ph (transpose:la l)) m)
  ::  Perspective divide.  (seems to be unnecessary)
  :: =.  pc  (div:la pc (get-col:la pc ~[3]))
  ::  Transform the points to screen space.
  =/  pix  (add-scalar:la (mul-scalar:la (submatrix:la ~[`[[~ ~]] `[[~ `3]]] pc) .5) .5)
  =.  pix  (set-col:la pix ~[0] (mul-scalar:la (get-col:la pix ~[0]) (sun:rs nx.view)))
  =.  pix  (set-col:la pix ~[1] (mul-scalar:la (get-col:la pix ~[1]) (sun:rs ny.view)))
  pix
::
++  render-point-cloud
  |=  [pts=ray:ls =view:ts =lens:ts]
  ^-  ray:ls
  ::  Initialize with black background.
  =/  img  (zeros:la [~[nx.view ny.view 3] 5 %i754 ~])
  ::  Project points to 2D.
  =/  p2d  (project-points pts view lens)
  =/  m  (snag 0 shape.meta.p2d)
  ::  Calculate depth (distance from camera).
  =/  dd  (sub-scalar:la pts .2)
  =/  d  (sqrt:sa (dot:la dd dd))
  =/  dn  %+  div:la
            (sub-scalar:la d (get-item:la (min:la d) ~[0]))
          (sub:la (max:la d) (min:la d))
  ::  Rasterize points.
  =|  idx=@
  |-
  ?:  =(m idx)  img
  =/  px  (get-item:la p2d ~[idx 0])
  =/  py  (get-item:la p2d ~[idx 1])
  ?:  &(&((lte .0 px) (lth px nx.view)) &((lte .0 py) (lth py ny.view)))
    $(idx +(idx))
  =/  shd  (get-item:la dn ~[idx])
  %=  $
    idx  +(idx)
    img  (set-item:la (set-item:la (set-item:la img ~[py px 0] shd) ~[py px 1] shd) ~[py px 2] shd)
  ==
::
++  deg2rad
  |=  deg=@rs
  ^-  @rs
  =/  pi  pi:rs:math
  =/  rad  (mul:rs (mul:rs pi .5.555555e-3) deg)
  rad
--
