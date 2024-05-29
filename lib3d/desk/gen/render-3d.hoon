  ::  /gen/render-3d
::::  Version ~2024.5.29 by ~lagrev-nocfep
::
::  Load a point cloud and render to a bitmap.
/-  ls=lagoon
/+  bmp,
    *lagoon,
    threedim
:-  %say
|=  $:  $:  now=@da             ::  timestamp
            eny=@uvJ            ::  entropy
            bec=beak            ::  clay beak
        ==
        $:                      ::  required arguments
            $:  data-file=@ta   ::  data file
            ==
            ~
        ==
        $:                      ::  optional arguments
            $:  angle=ray:ls    ::  camera view angle
            ==
            ~
        ==
    ==
~&  >  "Rendering /dat/{(trip data-file)}.hoon at {<now>}."
:-  %noun
::  Acquire 3D point cloud.
=/  fil  .^(@t %cx /(scot %p p.bec)/[q.bec]/(scot %da p.r.bec)/dat/[data-file]/hoon)
~&  >  "Loaded /dat/{(trip data-file)}.hoon."
=/  dat  (ride %noun fil)
=/  pts  ;;((list (list @rs)) +.q.dat)
~&  >  "Loaded {<(lent pts)>} points."
::  Define the Camera and Viewport: Set up the camera position and the parameters of the viewport (the 2D grid).
=/  img  (gen-3d:threedim pts [.-1 .1 .-1 .1 100 100])
~&  >>  meta.img
~&  >  "Rendered points."
=/  raw  (change:la (mul-scalar:la img .255) %uint 3)
~&  >>  meta.raw
~&  >  "Converted to bitmap."
(write-task:bmp /dat/[data-file]/bmp raw)
