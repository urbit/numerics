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
::  TODO (PENCILLED IN -- complex transcendentals, follow-on to the %unum
::  transcendental work).  This lib currently provides arithmetic + abs/conj/
::  re/im/pak only.  To let Saloon's Tier B (exp/sin/cos/log/sqrt/...) accept
::  %cplx rays the way it now accepts %unum, each width door needs the complex
::  elementary functions, then a %cplx branch in saloon +trans-scalar/+fun-scalar
::  (mirroring the %unum branch).  Closed forms over the real component ops:
::    cexp(a+bi)  = e^a * (cos b + i sin b)
::    clog(a+bi)  = log|z| + i*atan2(b, a)            (needs a real atan2)
::    csqrt(z)    = sqrt((|z|+a)/2) + i*sgn(b)*sqrt((|z|-a)/2)
::    csin/ccos   via the real sin/cos/sinh/cosh, cpow via cexp(w*clog z)
::  Depends on the real transcendentals (which #18's Chebyshev rewrite will
::  harden); accuracy will inherit whatever the real layer provides.
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
::  The @rd/@cd analogue of +cs over double-precision components.  Here 1+2i
::  packs to 0x4000.0000.0000.0000.3ff0.0000.0000.0000 and 3+4i to
::  0x4010.0000.0000.0000.4008.0000.0000.0000 (real in the low 64 bits).
::
++  cd
  |_  rnd=rounding-mode
  ::
  ::  Components
  ::
  ::    +re:  @cd -> @rd
  ::
  ::  Returns the real component (low 64 bits) of a packed @cd value.
  ::    Examples
  ::      > (~(re cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000)  ::  re(1+2i)
  ::      .~1
  ::  Source
  ++  re    |=(p=@ `@rd`(end [0 64] p))
  ::    +im:  @cd -> @rd
  ::
  ::  Returns the imaginary component (high 64 bits) of a packed @cd value.
  ::    Examples
  ::      > (~(im cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000)  ::  im(1+2i)
  ::      .~2
  ::  Source
  ++  im    |=(p=@ `@rd`(rsh [0 64] p))
  ::    +pak:  [@rd @rd] -> @cd
  ::
  ::  Packs real and imaginary @rd components into one @cd atom (real low).
  ::    Examples
  ::      > (~(pak cd:complex %n) .~1 .~2)
  ::      0x4000.0000.0000.0000.3ff0.0000.0000.0000
  ::  Source
  ++  pak   |=([r=@rd i=@rd] ^-(@ (con `@`r (lsh [0 64] `@`i))))
  ::  component helpers (internal)
  ++  fneg  |=(x=@rd (~(sub rd rnd) .~0 x))
  ++  fabs  |=(x=@rd ?:((~(lth rd rnd) x .~0) (~(sub rd rnd) .~0 x) x))
  ::
  ::  Constants
  ::
  ::    +zero:  @cd
  ::
  ::  The additive identity 0+0i.
  ::    Examples
  ::      > ~(zero cd:complex %n)
  ::      0x0
  ::  Source
  ++  zero  ^-(@ 0)
  ::    +one:  @cd
  ::
  ::  The multiplicative identity 1+0i.
  ::    Examples
  ::      > ~(one cd:complex %n)
  ::      0x3ff0.0000.0000.0000
  ::  Source
  ++  one   (pak .~1 .~0)
  ::
  ::  Arithmetic
  ::
  ::    +add:  [@cd @cd] -> @cd
  ::
  ::  Returns the sum of two complex values (component-wise).
  ::    Examples
  ::      > (~(add cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4010.0000.0000.0000.4008.0000.0000.0000)
  ::      0x4018.0000.0000.0000.4010.0000.0000.0000
  ::  Source
  ++  add   |=([p=@ q=@] (pak (~(add rd rnd) (re p) (re q)) (~(add rd rnd) (im p) (im q))))
  ::    +sub:  [@cd @cd] -> @cd
  ::
  ::  Returns the difference of two complex values (component-wise).
  ::    Examples
  ::      > (~(sub cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4010.0000.0000.0000.4008.0000.0000.0000)
  ::      0xc000.0000.0000.0000.c000.0000.0000.0000
  ::  Source
  ++  sub   |=([p=@ q=@] (pak (~(sub rd rnd) (re p) (re q)) (~(sub rd rnd) (im p) (im q))))
  ::    +neg:  @cd -> @cd
  ::
  ::  Returns the additive inverse -z.
  ::    Examples
  ::      > (~(neg cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000)  ::  -(1+2i)
  ::      0xc000.0000.0000.0000.bff0.0000.0000.0000
  ::  Source
  ++  neg   |=(p=@ (pak (fneg (re p)) (fneg (im p))))
  ::    +conj:  @cd -> @cd
  ::
  ::  Returns the complex conjugate (negate the imaginary component).
  ::    Examples
  ::      > (~(conj cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000)  ::  conj(1+2i)=1-2i
  ::      0xc000.0000.0000.0000.3ff0.0000.0000.0000
  ::  Source
  ++  conj  |=(p=@ (pak (re p) (fneg (im p))))
  ::    +mul:  [@cd @cd] -> @cd
  ::
  ::  Returns the complex product (ar+ai i)(br+bi i), rounded per component.
  ::    Examples
  ::      > (~(mul cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4010.0000.0000.0000.4008.0000.0000.0000)
  ::      0x4024.0000.0000.0000.c014.0000.0000.0000               ::  (1+2i)(3+4i)=-5+10i
  ::  Source
  ++  mul
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    %+  pak
      (~(sub rd rnd) (~(mul rd rnd) ar br) (~(mul rd rnd) ai bi))
    (~(add rd rnd) (~(mul rd rnd) ar bi) (~(mul rd rnd) ai br))
  ::    +div:  [@cd @cd] -> @cd
  ::
  ::  Returns the complex quotient via Smith's algorithm (scale by the larger
  ::  denominator component for numerical stability).
  ::    Examples
  ::      > (~(div cd:complex %n) 0x4000.0000.0000.0000.4000.0000.0000.0000 0x3ff0.0000.0000.0000.3ff0.0000.0000.0000)
  ::      0x4000.0000.0000.0000                                   ::  (2+2i)/(1+1i)=2
  ::  Source
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
  ::    +abs:  @cd -> @cd
  ::
  ::  Returns the modulus |z| as a real-valued complex [|z|, 0], computed via
  ::  hypot (scale by the larger component to avoid overflow).
  ::    Examples
  ::      > (~(abs cd:complex %n) 0x4010.0000.0000.0000.4008.0000.0000.0000)  ::  |3+4i|=5
  ::      0x4014.0000.0000.0000
  ::  Source
  ++  abs
    |=  p=@
    =/  xr  (fabs (re p))  =/  xi  (fabs (im p))
    ?:  &(=(`@rd`.~0 xr) =(`@rd`.~0 xi))  (pak .~0 .~0)
    ?:  (~(gte rd rnd) xr xi)
      =/  t  (~(div rd rnd) xi xr)
      (pak (~(mul rd rnd) xr (~(sqt rd rnd) (~(add rd rnd) .~1 (~(mul rd rnd) t t)))) .~0)
    =/  t  (~(div rd rnd) xr xi)
    (pak (~(mul rd rnd) xi (~(sqt rd rnd) (~(add rd rnd) .~1 (~(mul rd rnd) t t)))) .~0)
  ::
  ::  Comparison (no total order on complex -- only equality)
  ::
  ::    +equ:  [@cd @cd] -> ?
  ::
  ::  Loobean: bit-equality of both components (matches the %i754 convention).
  ::    Examples
  ::      > (~(equ cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4000.0000.0000.0000.3ff0.0000.0000.0000)
  ::      %.y
  ::      > (~(equ cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4010.0000.0000.0000.4008.0000.0000.0000)
  ::      %.n
  ::  Source
  ++  equ   |=([p=@ q=@] &(=((re p) (re q)) =((im p) (im q))))
  ::    +neq:  [@cd @cd] -> ?
  ::
  ::  Loobean negation of +equ.
  ::    Examples
  ::      > (~(neq cd:complex %n) 0x4000.0000.0000.0000.3ff0.0000.0000.0000 0x4010.0000.0000.0000.4008.0000.0000.0000)
  ::      %.y
  ::  Source
  ++  neq   |=([p=@ q=@] =(| (equ p q)))
  --
::    +ch:  complex-half (@ch), two @rh (16-bit) components.
::
::  The @rh/@ch analogue of +cs over half-precision components (32-bit element).
::  Here 1+2i packs to 0x4000.3c00 and 3+4i to 0x4400.4200.
::
++  ch
  |_  rnd=rounding-mode
  ::
  ::  Components
  ::
  ::    +re:  @ch -> @rh
  ::
  ::  Returns the real component (low 16 bits) of a packed @ch value.
  ::    Examples
  ::      > (~(re ch:complex %n) 0x4000.3c00)  ::  re(1+2i)
  ::      .~~1
  ::  Source
  ++  re    |=(p=@ `@rh`(end [0 16] p))
  ::    +im:  @ch -> @rh
  ::
  ::  Returns the imaginary component (high 16 bits) of a packed @ch value.
  ::    Examples
  ::      > (~(im ch:complex %n) 0x4000.3c00)  ::  im(1+2i)
  ::      .~~2
  ::  Source
  ++  im    |=(p=@ `@rh`(rsh [0 16] p))
  ::    +pak:  [@rh @rh] -> @ch
  ::
  ::  Packs real and imaginary @rh components into one @ch atom (real low).
  ::    Examples
  ::      > (~(pak ch:complex %n) .~~1 .~~2)
  ::      0x4000.3c00
  ::  Source
  ++  pak   |=([r=@rh i=@rh] ^-(@ (con `@`r (lsh [0 16] `@`i))))
  ::  component helpers (internal)
  ++  fneg  |=(x=@rh (~(sub rh rnd) .~~0 x))
  ++  fabs  |=(x=@rh ?:((~(lth rh rnd) x .~~0) (~(sub rh rnd) .~~0 x) x))
  ::
  ::  Constants
  ::
  ::    +zero:  @ch
  ::
  ::  The additive identity 0+0i.
  ::    Examples
  ::      > ~(zero ch:complex %n)
  ::      0x0
  ::  Source
  ++  zero  ^-(@ 0)
  ::    +one:  @ch
  ::
  ::  The multiplicative identity 1+0i.
  ::    Examples
  ::      > ~(one ch:complex %n)
  ::      0x3c00
  ::  Source
  ++  one   (pak .~~1 .~~0)
  ::
  ::  Arithmetic
  ::
  ::    +add:  [@ch @ch] -> @ch
  ::
  ::  Returns the sum of two complex values (component-wise).
  ::    Examples
  ::      > (~(add ch:complex %n) 0x4000.3c00 0x4400.4200)  ::  (1+2i)+(3+4i)=4+6i
  ::      0x4600.4400
  ::  Source
  ++  add   |=([p=@ q=@] (pak (~(add rh rnd) (re p) (re q)) (~(add rh rnd) (im p) (im q))))
  ::    +sub:  [@ch @ch] -> @ch
  ::
  ::  Returns the difference of two complex values (component-wise).
  ::    Examples
  ::      > (~(sub ch:complex %n) 0x4000.3c00 0x4400.4200)  ::  =-2-2i
  ::      0xc000.c000
  ::  Source
  ++  sub   |=([p=@ q=@] (pak (~(sub rh rnd) (re p) (re q)) (~(sub rh rnd) (im p) (im q))))
  ::    +neg:  @ch -> @ch
  ::
  ::  Returns the additive inverse -z.
  ::    Examples
  ::      > (~(neg ch:complex %n) 0x4000.3c00)  ::  -(1+2i)
  ::      0xc000.bc00
  ::  Source
  ++  neg   |=(p=@ (pak (fneg (re p)) (fneg (im p))))
  ::    +conj:  @ch -> @ch
  ::
  ::  Returns the complex conjugate (negate the imaginary component).
  ::    Examples
  ::      > (~(conj ch:complex %n) 0x4000.3c00)  ::  conj(1+2i)=1-2i
  ::      0xc000.3c00
  ::  Source
  ++  conj  |=(p=@ (pak (re p) (fneg (im p))))
  ::    +mul:  [@ch @ch] -> @ch
  ::
  ::  Returns the complex product (ar+ai i)(br+bi i), rounded per component.
  ::    Examples
  ::      > (~(mul ch:complex %n) 0x4000.3c00 0x4400.4200)  ::  (1+2i)(3+4i)=-5+10i
  ::      0x4900.c500
  ::  Source
  ++  mul
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    %+  pak
      (~(sub rh rnd) (~(mul rh rnd) ar br) (~(mul rh rnd) ai bi))
    (~(add rh rnd) (~(mul rh rnd) ar bi) (~(mul rh rnd) ai br))
  ::    +div:  [@ch @ch] -> @ch
  ::
  ::  Returns the complex quotient via Smith's algorithm (scale by the larger
  ::  denominator component for numerical stability).
  ::    Examples
  ::      > (~(div ch:complex %n) 0x4000.4000 0x3c00.3c00)  ::  (2+2i)/(1+1i)=2
  ::      0x4000
  ::  Source
  ++  div
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    ?:  (~(gte rh rnd) (fabs br) (fabs bi))
      =/  r   (~(div rh rnd) bi br)
      =/  dn  (~(add rh rnd) br (~(mul rh rnd) bi r))
      %+  pak
        (~(div rh rnd) (~(add rh rnd) ar (~(mul rh rnd) ai r)) dn)
      (~(div rh rnd) (~(sub rh rnd) ai (~(mul rh rnd) ar r)) dn)
    =/  r   (~(div rh rnd) br bi)
    =/  dn  (~(add rh rnd) (~(mul rh rnd) br r) bi)
    %+  pak
      (~(div rh rnd) (~(add rh rnd) (~(mul rh rnd) ar r) ai) dn)
    (~(div rh rnd) (~(sub rh rnd) (~(mul rh rnd) ai r) ar) dn)
  ::    +abs:  @ch -> @ch
  ::
  ::  Returns the modulus |z| as a real-valued complex [|z|, 0], via hypot.
  ::    Examples
  ::      > (~(abs ch:complex %n) 0x4400.4200)  ::  |3+4i|=5
  ::      0x4500
  ::  Source
  ++  abs
    |=  p=@
    =/  xr  (fabs (re p))  =/  xi  (fabs (im p))
    ?:  &(=(`@rh`.~~0 xr) =(`@rh`.~~0 xi))  (pak .~~0 .~~0)
    ?:  (~(gte rh rnd) xr xi)
      =/  t  (~(div rh rnd) xi xr)
      (pak (~(mul rh rnd) xr (~(sqt rh rnd) (~(add rh rnd) .~~1 (~(mul rh rnd) t t)))) .~~0)
    =/  t  (~(div rh rnd) xr xi)
    (pak (~(mul rh rnd) xi (~(sqt rh rnd) (~(add rh rnd) .~~1 (~(mul rh rnd) t t)))) .~~0)
  ::
  ::  Comparison (no total order on complex -- only equality)
  ::
  ::    +equ:  [@ch @ch] -> ?
  ::
  ::  Loobean: bit-equality of both components.
  ::    Examples
  ::      > (~(equ ch:complex %n) 0x4000.3c00 0x4000.3c00)
  ::      %.y
  ::      > (~(equ ch:complex %n) 0x4000.3c00 0x4400.4200)
  ::      %.n
  ::  Source
  ++  equ   |=([p=@ q=@] &(=((re p) (re q)) =((im p) (im q))))
  ::    +neq:  [@ch @ch] -> ?
  ::
  ::  Loobean negation of +equ.
  ::    Examples
  ::      > (~(neq ch:complex %n) 0x4000.3c00 0x4400.4200)
  ::      %.y
  ::  Source
  ++  neq   |=([p=@ q=@] =(| (equ p q)))
  --
::    +cq:  complex-quad (@cq), two @rq (128-bit) components.
::
::  The @rq/@cq analogue of +cs over quad-precision components (256-bit element).
::  Here 1+2i packs to
::  0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000
::  and 3+4i to
::  0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000.
::
++  cq
  |_  rnd=rounding-mode
  ::
  ::  Components
  ::
  ::    +re:  @cq -> @rq
  ::
  ::  Returns the real component (low 128 bits) of a packed @cq value.
  ::    Examples
  ::      > (~(re cq:complex %n) 0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000)
  ::      .~~~1
  ::  Source
  ++  re    |=(p=@ `@rq`(end [0 128] p))
  ::    +im:  @cq -> @rq
  ::
  ::  Returns the imaginary component (high 128 bits) of a packed @cq value.
  ::    Examples
  ::      > (~(im cq:complex %n) 0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000)
  ::      .~~~2
  ::  Source
  ++  im    |=(p=@ `@rq`(rsh [0 128] p))
  ::    +pak:  [@rq @rq] -> @cq
  ::
  ::  Packs real and imaginary @rq components into one @cq atom (real low).
  ::    Examples
  ::      > (~(pak cq:complex %n) .~~~1 .~~~2)
  ::      0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  pak   |=([r=@rq i=@rq] ^-(@ (con `@`r (lsh [0 128] `@`i))))
  ::  component helpers (internal)
  ++  fneg  |=(x=@rq (~(sub rq rnd) .~~~0 x))
  ++  fabs  |=(x=@rq ?:((~(lth rq rnd) x .~~~0) (~(sub rq rnd) .~~~0 x) x))
  ::
  ::  Constants
  ::
  ::    +zero:  @cq
  ::
  ::  The additive identity 0+0i.
  ::    Examples
  ::      > ~(zero cq:complex %n)
  ::      0x0
  ::  Source
  ++  zero  ^-(@ 0)
  ::    +one:  @cq
  ::
  ::  The multiplicative identity 1+0i.
  ::    Examples
  ::      > ~(one cq:complex %n)
  ::      0x3fff.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  one   (pak .~~~1 .~~~0)
  ::
  ::  Arithmetic
  ::
  ::    +add:  [@cq @cq] -> @cq
  ::
  ::  Returns the sum of two complex values (component-wise).  (1+2i)+(3+4i)=4+6i.
  ::    Examples
  ::      > (~(add cq:complex %n) `@`0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000 0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000)
  ::      0x4001.8000.0000.0000.0000.0000.0000.0000.4001.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  add   |=([p=@ q=@] (pak (~(add rq rnd) (re p) (re q)) (~(add rq rnd) (im p) (im q))))
  ::    +sub:  [@cq @cq] -> @cq
  ::
  ::  Returns the difference of two complex values (component-wise).  =-2-2i.
  ::    Examples
  ::      > (~(sub cq:complex %n) `@`0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000 0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000)
  ::      0xc000.0000.0000.0000.0000.0000.0000.0000.c000.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  sub   |=([p=@ q=@] (pak (~(sub rq rnd) (re p) (re q)) (~(sub rq rnd) (im p) (im q))))
  ::    +neg:  @cq -> @cq
  ::
  ::  Returns the additive inverse -z.  -(1+2i)=-1-2i.
  ::    Examples
  ::      > (~(neg cq:complex %n) 0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000)
  ::      0xc000.0000.0000.0000.0000.0000.0000.0000.bfff.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  neg   |=(p=@ (pak (fneg (re p)) (fneg (im p))))
  ::    +conj:  @cq -> @cq
  ::
  ::  Returns the complex conjugate (negate the imaginary component).  conj(1+2i)=1-2i.
  ::    Examples
  ::      > (~(conj cq:complex %n) 0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000)
  ::      0xc000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  conj  |=(p=@ (pak (re p) (fneg (im p))))
  ::    +mul:  [@cq @cq] -> @cq
  ::
  ::  Returns the complex product, rounded per component.  (1+2i)(3+4i)=-5+10i.
  ::    Examples
  ::      > (~(mul cq:complex %n) `@`0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000 0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000)
  ::      0x4002.4000.0000.0000.0000.0000.0000.0000.c001.4000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  mul
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    %+  pak
      (~(sub rq rnd) (~(mul rq rnd) ar br) (~(mul rq rnd) ai bi))
    (~(add rq rnd) (~(mul rq rnd) ar bi) (~(mul rq rnd) ai br))
  ::    +div:  [@cq @cq] -> @cq
  ::
  ::  Returns the complex quotient via Smith's algorithm.  (2+2i)/(1+1i)=2.
  ::    Examples
  ::      > (~(div cq:complex %n) `@`0x4000.0000.0000.0000.0000.0000.0000.0000.4000.0000.0000.0000.0000.0000.0000.0000 0x3fff.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000)
  ::      0x4000.0000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  div
    |=  [p=@ q=@]
    =/  ar  (re p)  =/  ai  (im p)
    =/  br  (re q)  =/  bi  (im q)
    ?:  (~(gte rq rnd) (fabs br) (fabs bi))
      =/  r   (~(div rq rnd) bi br)
      =/  dn  (~(add rq rnd) br (~(mul rq rnd) bi r))
      %+  pak
        (~(div rq rnd) (~(add rq rnd) ar (~(mul rq rnd) ai r)) dn)
      (~(div rq rnd) (~(sub rq rnd) ai (~(mul rq rnd) ar r)) dn)
    =/  r   (~(div rq rnd) br bi)
    =/  dn  (~(add rq rnd) (~(mul rq rnd) br r) bi)
    %+  pak
      (~(div rq rnd) (~(add rq rnd) (~(mul rq rnd) ar r) ai) dn)
    (~(div rq rnd) (~(sub rq rnd) (~(mul rq rnd) ai r) ar) dn)
  ::    +abs:  @cq -> @cq
  ::
  ::  Returns the modulus |z| as a real-valued complex [|z|, 0], via hypot.  |3+4i|=5.
  ::    Examples
  ::      > (~(abs cq:complex %n) 0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000)
  ::      0x4001.4000.0000.0000.0000.0000.0000.0000
  ::  Source
  ++  abs
    |=  p=@
    =/  xr  (fabs (re p))  =/  xi  (fabs (im p))
    ?:  &(=(`@rq`.~~~0 xr) =(`@rq`.~~~0 xi))  (pak .~~~0 .~~~0)
    ?:  (~(gte rq rnd) xr xi)
      =/  t  (~(div rq rnd) xi xr)
      (pak (~(mul rq rnd) xr (~(sqt rq rnd) (~(add rq rnd) .~~~1 (~(mul rq rnd) t t)))) .~~~0)
    =/  t  (~(div rq rnd) xr xi)
    (pak (~(mul rq rnd) xi (~(sqt rq rnd) (~(add rq rnd) .~~~1 (~(mul rq rnd) t t)))) .~~~0)
  ::
  ::  Comparison (no total order on complex -- only equality)
  ::
  ::    +equ:  [@cq @cq] -> ?
  ::
  ::  Loobean: bit-equality of both components.
  ::    Examples
  ::      > =z  0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000  ::  1+2i
  ::      > (~(equ cq:complex %n) z z)
  ::      %.y
  ::  Source
  ++  equ   |=([p=@ q=@] &(=((re p) (re q)) =((im p) (im q))))
  ::    +neq:  [@cq @cq] -> ?
  ::
  ::  Loobean negation of +equ.
  ::    Examples
  ::      > (~(neq cq:complex %n) `@`0x4000.0000.0000.0000.0000.0000.0000.0000.3fff.0000.0000.0000.0000.0000.0000.0000 0x4001.0000.0000.0000.0000.0000.0000.0000.4000.8000.0000.0000.0000.0000.0000.0000)
  ::      %.y
  ::  Source
  ++  neq   |=([p=@ q=@] =(| (equ p q)))
  --
--
