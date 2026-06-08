::::  /lib/complex -- complex-number arithmetic over IEEE-754 component floats
::::  ~2026.6.5
::
::  A complex value is a packed atom with the REAL component in the low
::  2^(bloq-1) bits and the IMAGINARY component in the high bits.  This matches
::  SoftBLAS `complexN_t { real; imag }` over little-endian bytes, so an array
::  of these packed elements IS a BLAS interleaved complex array (re0, im0,
::  re1, im1, ...) and a future jet can cast the data straight to complexN_t*.
::
::  Aura family @c (complex), parallel to @r (real); each @c? is a packed pair
::  of the corresponding @r? component float:
::
::    @ch  complex-half    = 2 x @rh (16b)  = 32-bit element
::    @cs  complex-single  = 2 x @rs (32b)  = 64-bit element   (BLAS  c )
::    @cd  complex-double  = 2 x @rd (64b)  = 128-bit element  (BLAS  z )
::    @cq  complex-quad    = 2 x @rq (128b) = 256-bit element
::
::  Each width is a door keyed on the rounding mode; arms take and return
::  PACKED complex atoms.  @cs and @cd ship first (BLAS c/z); @ch/@cq are
::  additive copies later.  Complex has no total order, so there is
::  intentionally no gth/gte/lth/lte -- only equ/neq.
::
|%
+$  rounding-mode  ?(%n %u %d %z)
::    +cs:  complex-single (@cs), two @rs (32-bit) components.
::
::  Door keyed on the rounding mode.  In the examples below 1+2i packs to
::  0x4000.0000.3f80.0000 and 3+4i to 0x4080.0000.4040.0000 (real in the low
::  32 bits).
::
++  cs
  |_  rnd=rounding-mode
  ::
  ::  Components
  ::
  ::    +re:  @cs -> @rs
  ::
  ::  Returns the real component (low 32 bits) of a packed @cs value.
  ::    Examples
  ::      > (~(re cs:complex %n) 0x4000.0000.3f80.0000)  ::  re(1+2i)
  ::      .1
  ::  Source
  ++  re    |=(p=@ `@rs`(end [0 32] p))
  ::    +im:  @cs -> @rs
  ::
  ::  Returns the imaginary component (high 32 bits) of a packed @cs value.
  ::    Examples
  ::      > (~(im cs:complex %n) 0x4000.0000.3f80.0000)  ::  im(1+2i)
  ::      .2
  ::  Source
  ++  im    |=(p=@ `@rs`(rsh [0 32] p))
  ::    +pak:  [@rs @rs] -> @cs
  ::
  ::  Packs real and imaginary @rs components into one @cs atom (real low).
  ::    Examples
  ::      > (~(pak cs:complex %n) .1 .2)
  ::      0x4000.0000.3f80.0000
  ::  Source
  ++  pak   |=([r=@rs i=@rs] ^-(@ (con `@`r (lsh [0 32] `@`i))))
  ::  component helpers (internal)
  ++  fneg  |=(x=@rs (~(sub rs rnd) .0 x))
  ++  fabs  |=(x=@rs ?:((~(lth rs rnd) x .0) (~(sub rs rnd) .0 x) x))
  ::
  ::  Constants
  ::
  ::    +zero:  @cs
  ::
  ::  The additive identity 0+0i.
  ::    Examples
  ::      > ~(zero cs:complex %n)
  ::      0x0
  ::  Source
  ++  zero  ^-(@ 0)
  ::    +one:  @cs
  ::
  ::  The multiplicative identity 1+0i.
  ::    Examples
  ::      > ~(one cs:complex %n)
  ::      0x3f80.0000
  ::  Source
  ++  one   (pak .1 .0)
  ::
  ::  Arithmetic
  ::
  ::    +add:  [@cs @cs] -> @cs
  ::
  ::  Returns the sum of two complex values (component-wise).
  ::    Examples
  ::      > (~(add cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000)
  ::      0x40c0.0000.4080.0000
  ::  Source
  ++  add   |=([p=@ q=@] (pak (~(add rs rnd) (re p) (re q)) (~(add rs rnd) (im p) (im q))))
  ::    +sub:  [@cs @cs] -> @cs
  ::
  ::  Returns the difference of two complex values (component-wise).
  ::    Examples
  ::      > (~(sub cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000)
  ::      0xc000.0000.c000.0000
  ::  Source
  ++  sub   |=([p=@ q=@] (pak (~(sub rs rnd) (re p) (re q)) (~(sub rs rnd) (im p) (im q))))
  ::    +neg:  @cs -> @cs
  ::
  ::  Returns the additive inverse -z.
  ::    Examples
  ::      > (~(neg cs:complex %n) 0x4000.0000.3f80.0000)  ::  -(1+2i)
  ::      0xc000.0000.bf80.0000
  ::  Source
  ++  neg   |=(p=@ (pak (fneg (re p)) (fneg (im p))))
  ::    +conj:  @cs -> @cs
  ::
  ::  Returns the complex conjugate (negate the imaginary component).
  ::    Examples
  ::      > (~(conj cs:complex %n) 0x4000.0000.3f80.0000)  ::  conj(1+2i)=1-2i
  ::      0xc000.0000.3f80.0000
  ::  Source
  ++  conj  |=(p=@ (pak (re p) (fneg (im p))))
  ::    +mul:  [@cs @cs] -> @cs
  ::
  ::  Returns the complex product (ar+ai i)(br+bi i), rounded per component.
  ::    Examples
  ::      > (~(mul cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000)
  ::      0x4120.0000.c0a0.0000                                ::  (1+2i)(3+4i)=-5+10i
  ::  Source
  ++  mul
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    %+  pak
      (~(sub rs rnd) (~(mul rs rnd) ar br) (~(mul rs rnd) ai bi))
    (~(add rs rnd) (~(mul rs rnd) ar bi) (~(mul rs rnd) ai br))
  ::    +div:  [@cs @cs] -> @cs
  ::
  ::  Returns the complex quotient via Smith's algorithm (scale by the larger
  ::  denominator component for numerical stability).
  ::    Examples
  ::      > (~(div cs:complex %n) 0x4000.0000.4000.0000 0x3f80.0000.3f80.0000)
  ::      0x4000.0000                                          ::  (2+2i)/(1+1i)=2
  ::  Source
  ++  div
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    ?:  (~(gte rs rnd) (fabs br) (fabs bi))
      =/  r   (~(div rs rnd) bi br)
      =/  dn  (~(add rs rnd) br (~(mul rs rnd) bi r))
      %+  pak
        (~(div rs rnd) (~(add rs rnd) ar (~(mul rs rnd) ai r)) dn)
      (~(div rs rnd) (~(sub rs rnd) ai (~(mul rs rnd) ar r)) dn)
    =/  r   (~(div rs rnd) br bi)
    =/  dn  (~(add rs rnd) (~(mul rs rnd) br r) bi)
    %+  pak
      (~(div rs rnd) (~(add rs rnd) (~(mul rs rnd) ar r) ai) dn)
    (~(div rs rnd) (~(sub rs rnd) (~(mul rs rnd) ai r) ar) dn)
  ::    +abs:  @cs -> @cs
  ::
  ::  Returns the modulus |z| as a real-valued complex [|z|, 0], computed via
  ::  hypot (scale by the larger component to avoid overflow).
  ::    Examples
  ::      > (~(abs cs:complex %n) 0x4080.0000.4040.0000)       ::  |3+4i|=5
  ::      0x40a0.0000
  ::  Source
  ++  abs
    |=  p=@
    =/  xr  (fabs (re p))  =/  xi  (fabs (im p))
    ?:  &(=(`@rs`.0 xr) =(`@rs`.0 xi))  (pak .0 .0)
    ?:  (~(gte rs rnd) xr xi)
      =/  t  (~(div rs rnd) xi xr)
      (pak (~(mul rs rnd) xr (~(sqt rs rnd) (~(add rs rnd) .1 (~(mul rs rnd) t t)))) .0)
    =/  t  (~(div rs rnd) xr xi)
    (pak (~(mul rs rnd) xi (~(sqt rs rnd) (~(add rs rnd) .1 (~(mul rs rnd) t t)))) .0)
  ::
  ::  Comparison (no total order on complex -- only equality)
  ::
  ::    +equ:  [@cs @cs] -> ?
  ::
  ::  Loobean: bit-equality of both components (matches the %i754 convention).
  ::    Examples
  ::      > (~(equ cs:complex %n) 0x4000.0000.3f80.0000 0x4000.0000.3f80.0000)
  ::      %.y
  ::      > (~(equ cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000)
  ::      %.n
  ::  Source
  ++  equ   |=([p=@ q=@] &(=((re p) (re q)) =((im p) (im q))))
  ::    +neq:  [@cs @cs] -> ?
  ::
  ::  Loobean negation of +equ.
  ::    Examples
  ::      > (~(neq cs:complex %n) 0x4000.0000.3f80.0000 0x4080.0000.4040.0000)
  ::      %.y
  ::  Source
  ++  neq   |=([p=@ q=@] =(| (equ p q)))
  --
::    +cd:  complex-double (@cd), two @rd (64-bit) components.
::
::  The @rd/@cd analogue of +cs: identical API and semantics over double-
::  precision components.  Same arms (re/im/pak/zero/one/add/sub/neg/conj/mul/
::  div/abs/equ/neq); see +cs for per-arm documentation.  Here 1+2i packs to
::  0x4000.0000.0000.0000.3ff0.0000.0000.0000.
::    Examples
::      > (~(mul cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4010.0000.0000.0000.4008.0000.0000.0000)
::      0x4024.0000.0000.0000.c014.0000.0000.0000               ::  (1+2i)(3+4i)=-5+10i
::      > ~(one cd:complex %n)
::      0x3ff0.0000.0000.0000                                   ::  1+0i
::  Source
++  cd
  |_  rnd=rounding-mode
  ++  re    |=(p=@ `@rd`(end [0 64] p))
  ++  im    |=(p=@ `@rd`(rsh [0 64] p))
  ++  pak   |=([r=@rd i=@rd] ^-(@ (con `@`r (lsh [0 64] `@`i))))
  ++  fneg  |=(x=@rd (~(sub rd rnd) .~0 x))
  ++  fabs  |=(x=@rd ?:((~(lth rd rnd) x .~0) (~(sub rd rnd) .~0 x) x))
  ++  zero  ^-(@ 0)
  ++  one   (pak .~1 .~0)
  ++  add   |=([p=@ q=@] (pak (~(add rd rnd) (re p) (re q)) (~(add rd rnd) (im p) (im q))))
  ++  sub   |=([p=@ q=@] (pak (~(sub rd rnd) (re p) (re q)) (~(sub rd rnd) (im p) (im q))))
  ++  neg   |=(p=@ (pak (fneg (re p)) (fneg (im p))))
  ++  conj  |=(p=@ (pak (re p) (fneg (im p))))
  ++  mul
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    %+  pak
      (~(sub rd rnd) (~(mul rd rnd) ar br) (~(mul rd rnd) ai bi))
    (~(add rd rnd) (~(mul rd rnd) ar bi) (~(mul rd rnd) ai br))
  ++  div
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    ?:  (~(gte rd rnd) (fabs br) (fabs bi))
      =/  r   (~(div rd rnd) bi br)
      =/  dn  (~(add rd rnd) br (~(mul rd rnd) bi r))
      %+  pak
        (~(div rd rnd) (~(add rd rnd) ar (~(mul rd rnd) ai r)) dn)
      (~(div rd rnd) (~(sub rd rnd) ai (~(mul rd rnd) ar r)) dn)
    =/  r   (~(div rd rnd) br bi)
    =/  dn  (~(add rd rnd) (~(mul rd rnd) br r) bi)
    %+  pak
      (~(div rd rnd) (~(add rd rnd) (~(mul rd rnd) ar r) ai) dn)
    (~(div rd rnd) (~(sub rd rnd) (~(mul rd rnd) ai r) ar) dn)
  ++  abs
    |=  p=@
    =/  xr  (fabs (re p))  =/  xi  (fabs (im p))
    ?:  &(=(`@rd`.~0 xr) =(`@rd`.~0 xi))  (pak .~0 .~0)
    ?:  (~(gte rd rnd) xr xi)
      =/  t  (~(div rd rnd) xi xr)
      (pak (~(mul rd rnd) xr (~(sqt rd rnd) (~(add rd rnd) .~1 (~(mul rd rnd) t t)))) .~0)
    =/  t  (~(div rd rnd) xr xi)
    (pak (~(mul rd rnd) xi (~(sqt rd rnd) (~(add rd rnd) .~1 (~(mul rd rnd) t t)))) .~0)
  ++  equ   |=([p=@ q=@] &(=((re p) (re q)) =((im p) (im q))))
  ++  neq   |=([p=@ q=@] =(| (equ p q)))
  --
--
