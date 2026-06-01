/+  twoc
::::  /lib/fixed -- fixed-point arithmetic
::::  Version ~2024.1.15 by ~lagrev-nocfep and ~mopfel-winrux
::::  ~2026.6.1: rebuilt on /lib/twoc for correct signed mul/div/scale.
::
::  A fixed-point number is a two's-complement integer of N = a + b + 1 bits
::  interpreted as value / 2^b:
::
::    value = 1/2^b * ( -2^(N-1) x_(N-1) + sum_{n=0}^{N-2} 2^n x_n )
::
::  where a = integer bits, b = fractional bits, and the leading bit is the
::  sign.  (Yates2023, http://www.digitalsignallabs.com/downloads/fp.pdf)
::
::  Arithmetic is signed and MODULAR (wraps mod 2^N on overflow), delegating
::  to the width-keyed +twid door of /lib/twoc so the sign handling is shared
::  with the rest of numerics rather than reimplemented here.
::
|%
+$  prec  [a=@ b=@]   ::  fixed-point precision, a+b+1=bloq
::  +wid: total bit width N = a + b + 1 of a precision.
++  wid  |=(=prec ^-(@ :(^add a.prec b.prec 1)))
::  +ng: the /lib/twoc width-keyed door at a given precision's N.
++  ng   |=(=prec ~(. twid:twoc (wid prec)))
::    +add: [@ prec @ prec] -> @
::
::  Add two fixed-point numbers of equal precision (signed, wrapping).
::    Examples
::      > (add 0x180 [8 8] 0x240 [8 8])   :: 1.5 + 2.25
::      0x3c0                             :: 3.75
::  Source
++  add
  |=  [x=@ xprec=prec y=@ yprec=prec]
  ?>  =(xprec yprec)
  (add:(ng xprec) x y)
::    +sub: [@ prec @ prec] -> @
::
::  Subtract two fixed-point numbers of equal precision (signed, wrapping).
::    Examples
::      > (sub 0x100 [8 8] 0x200 [8 8])   :: 1.0 - 2.0
::      0x1.ff00                          :: -1.0 in N=17 bits
::  Source
++  sub
  |=  [x=@ xprec=prec y=@ yprec=prec]
  ?>  =(xprec yprec)
  (sub:(ng xprec) x y)
::    +mul: [@ prec @ prec] -> [@ prec]
::
::  Multiply two fixed-point numbers (signed).  The product precision widens:
::  integer bits a.x + a.y + 1, fractional bits b.x + b.y.
::    Examples
::      > (mul 0x180 [8 8] 0x200 [8 8])   :: 1.5 * 2.0
::      [0x3.0000 17 16]                  :: 3.0 in q17.16
::      > (mul (neg 0x180 [8 8]) [8 8] 0x200 [8 8])   :: -1.5 * 2.0
::      [0x3.fffd.0000 17 16]             :: -3.0
::  Source
++  mul
  |=  [x=@ xprec=prec y=@ yprec=prec]
  ^-  [@ prec]
  =/  prodp=prec  [:(^add a.xprec a.yprec 1) (^add b.xprec b.yprec)]
  ::  re-encode each operand at the product width, then signed multiply.
  =/  px  (s-to-twoc:(ng prodp) (to-s x xprec))
  =/  py  (s-to-twoc:(ng prodp) (to-s y yprec))
  [(mul:(ng prodp) px py) prodp]
::    +div: [@ prec @ prec] -> [@ prec]
::
::  Divide two fixed-point numbers of equal precision (signed), keeping the
::  input precision.  Division by zero crashes.
::    Examples
::      > (div (neg 0x300 [8 8]) [8 8] 0x200 [8 8])   :: -3.0 / 2.0
::      [0x1.fe80 8 8]                                :: -1.5
::  Source
++  div
  |=  [x=@ xprec=prec y=@ yprec=prec]
  ^-  [@ prec]
  ?>  ~|(%argument-precision-must-match =(xprec yprec))
  ::  scale the dividend up by 2^b before the signed divide so the quotient
  ::  lands back in xprec's fixed-point scale: (x * 2^b) / y.
  =/  sx   (to-s x xprec)
  =/  num  (new:si (syn:si sx) (lsh [0 b.xprec] (abs:si sx)))
  =/  qa   (fra:si num (to-s y yprec))
  [(s-to-twoc:(ng xprec) qa) xprec]
::    +to-s: [@ prec] -> @s
::
::  Decode a fixed-point bit pattern to a hoon signed integer (the stored
::  two's-complement value, fractional bits included).
::  Source
++  to-s
  |=  [x=@ =prec]
  ^-  @s
  (twoc-to-s:(ng prec) x)
::    +neg: [@ prec] -> @
::
::  Two's-complement negation of a fixed-point number at its own width.
::  (Renamed from the old +twoc, which shadowed the /lib/twoc import face
::  and duplicated a name; behavior unchanged: -x.)
::    Examples
::      > (neg 0x100 [8 8])   :: -(1.0)
::      0x1.ff00
::  Source
++  neg
  |=  [x=@ =prec]
  ^-  @
  (neg:(ng prec) x)
::    +hi: [@ prec @] -> @   (high n bits)
::  Source
++  hi
  |=  [x=@ =prec n=@]
  ?>  (lte n (wid prec))
  (rsh [0 (^sub (wid prec) n)] x)
::    +lo: [@ prec @] -> @   (low n bits)
::  Source
++  lo
  |=  [x=@ =prec n=@]
  ?>  (lte n (wid prec))
  (dis x (dec (bex n)))
::    +scale: [@ prec prec] -> @
::
::  Re-quantize x from xprec to yprec (signed): shift the fractional part to
::  the new b, then re-encode at the new width.  Saturation/rounding are not
::  applied; high bits beyond yprec's range are dropped (wrap).
::    Examples
::      > (scale 0x10 [4 4] [8 8])   :: 1.0 in q4.4 -> q8.8
::      0x100
::      > (scale 0x100 [8 8] [4 4])  :: 1.0 in q8.8 -> q4.4
::      0x10
::  Source
++  scale
  |=  [x=@ xprec=prec yprec=prec]
  ^-  @
  ::  decode to the signed value, rescale the fractional precision, re-encode.
  =/  sx  (to-s x xprec)
  =/  sg  (syn:si sx)
  =/  mg  (abs:si sx)
  =/  mg2  ?:  (gth b.yprec b.xprec)
             (lsh [0 (^sub b.yprec b.xprec)] mg)
           (rsh [0 (^sub b.xprec b.yprec)] mg)
  (s-to-twoc:(ng yprec) (new:si sg mg2))
--
