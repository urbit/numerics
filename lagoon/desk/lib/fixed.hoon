  ::  /lib/fixed
::::  Version ~2024.1.15 by ~lagrev-nocfep and ~mopfel-winrux
::
|%
::
::  \frac{1}{2^b}\left(-2^{N-1}x_{N-1}+\sum_{n=0}^{N-2} 2^n x_n \right)
::
::  where
::    N = a + b + 1, 2^bloq size
::    a = integer bits
::    b = fractional bits
::
::  Yates2023 http://www.digitalsignallabs.com/downloads/fp.pdf
::
+$  prec  [a=@ b=@]   ::  fixed-point precision, a+b+1=bloq
::
++  add
  |=  [x=@ xprec=prec y=@ yprec=prec]
  ?>  =(xprec yprec)
  =/  sum  (^add x y)
  =/  sump=prec  [+(a.xprec) b.xprec]
  (lo sum sump :(^add a.xprec b.xprec 1))
++  sub
  |=  [x=@ xprec=prec y=@ yprec=prec]
  ?>  =(xprec yprec)
  =/  n  :(^add a.xprec b.xprec 1)
  (add x xprec (twoc y n) yprec)
++  mul
  |=  [x=@ xprec=prec y=@ yprec=prec]
  =/  prod  (^mul x y)
  =/  prodp=prec  [:(^add a.xprec a.yprec 1) (^add b.xprec b.yprec)]
  [prod prodp]
++  div  :: TODO This is not correct yet.
  |=  [x=@ xprec=prec y=@ yprec=prec]
  =/  quot  (^div x y)
  =/  quotp=prec  [:(^add a.xprec b.yprec 1) (^add a.yprec b.xprec)]
  [quot quotp]
++  hi
  |=  [x=@ =prec n=@]
  ?>  (lte n :(^add a.prec b.prec 1))
  =/  a  a.prec
  =/  b  b.prec
  ?:  =(0b0 (dis x (lsh [0 (^add a.prec b.prec)] 0b1)))
    :: positive number or zero
    (rsh [0 (^sub n +(a))] x)
  :: negative number
  =/  mask  (lsh [0 n] (rep [0 1] (reap (^sub :(^add a b 1) n) 0b1)))
  (con (rsh [0 (^sub n +(a))] x) mask)
++  lo
  |=  [x=@ =prec n=@]
  ?>  (lte n :(^add a.prec b.prec 1))
  =/  b  b.prec
  =/  a  (^sub n (^add b 1))
  =/  mask  (rep [0 1] (reap n 0b1))
  (dis x mask)
::  Produce twos-complement negation of x.
++  twoc
  |=  [x=@ n=@]
  %+  dis
    (rep [0 1] (reap n 0b1))
  (^add (mix x (rep [0 1] (reap n 0b1))) 0b1)
++  scale
  |=  [x=@ xprec=prec yprec=prec]
  =/  nx  :(^add a.xprec b.xprec 1)
  =/  lox  (lo x xprec nx)
  =/  hix  (hi x xprec nx)
  =/  loy  ?:  (gth b.xprec b.yprec)
    (rsh [0 (^sub b.xprec b.yprec)] x)
  (lsh [0 (^sub b.yprec b.xprec)] x)
  =/  hiy  ?:  (gth a.xprec a.yprec)
    (rsh [0 (^sub a.xprec a.yprec)] x)
  (lsh [0 (^sub a.yprec a.xprec)] x)
  (con (lsh hiy (^sub a.xprec b.xprec)) loy)
--
