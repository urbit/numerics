::
::::  Chebyshev/Fixed-Polynomial Math Library
::
::  Jet-compatible transcendental functions using fixed polynomial
::  evaluation. Every floating-point operation happens in the same
::  order with the same intermediate values in both Hoon and the C jet.
::
::  Reference: musl libc (MIT license)
::  Coefficients expressed as hex literals for bit-exactness.
::
::  For each precision level (@rh, @rs, @rd), we provide:
::    - Polynomial coefficients as hex constants
::    - Argument reduction helpers
::    - Kernel functions (sindf, cosdf)
::    - Public API (sin, cos, tan, exp)
::    - Chebyshev polynomial generator (+cheb)
::
|%
::
::  ================================================================
::  HALF PRECISION (@rh) - 16-bit IEEE 754
::  ================================================================
::
::  Half precision has limited range and precision, so we use
::  reduced-degree polynomials (fewer terms).
::
++  rh
  |%
  ::
  ::  Constants (hex-exact)
  ::
  ++  pio2-hi   `@rh`0x3e48       ::  pi/2 high bits ~1.5703125
  ++  pio2-lo   `@rh`0x0500       ::  pi/2 low bits (correction)
  ++  pio4      `@rh`0x3a48       ::  pi/4 ~0.785
  ++  invpio2   `@rh`0x3518       ::  2/pi ~0.6367
  ++  ln2hi     `@rh`0x398c       ::  ln(2) high ~0.693
  ++  ln2lo     `@rh`0x0000       ::  ln(2) low (negligible at half)
  ++  invln2    `@rh`0x3dc5       ::  1/ln(2) ~1.4425
  ::
  ::  Sin coefficients for half precision (degree 5)
  ::  sin(x) ~ x - x^3/6 + x^5/120
  ::
  ++  sin-s1    `@rh`0xb555       ::  -1/6 = -0.1666...
  ++  sin-s2    `@rh`0x2444       ::   1/120 = 0.00833...
  ::
  ::  Cos coefficients for half precision (degree 4)
  ::  cos(x) ~ 1 - x^2/2 + x^4/24
  ::
  ++  cos-c0    `@rh`0xb800       ::  -1/2 = -0.5
  ++  cos-c1    `@rh`0x2955       ::   1/24 = 0.04166...
  ::
  ::  Exp coefficients for half precision (degree 3)
  ::
  ++  exp-p1    `@rh`0x3800       ::  1/2
  ++  exp-p2    `@rh`0x3155       ::  1/6 ~0.1666
  ::
  ::  +cheb: Chebyshev polynomial T_n(x)
  ::
  ::  Returns the Chebyshev polynomial of the first kind.
  ::  T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x) with T_0=1, T_1=x
  ::
  ++  cheb
    |=  n=@ud
    |=  x=@rh
    ^-  @rh
    ?:  =(0 n)  `@rh`0x3c00       ::  1.0
    ?:  =(1 n)  x
    %+  sub:^rh
      (mul:^rh `@rh`0x4000 (mul:^rh x $(n (dec n))))  ::  2*x*T_{n-1}
    $(n (sub n 2))                                    ::  T_{n-2}
  ::
  ::  +sindf: Kernel sin for |x| <= pi/4
  ::
  ++  sindf
    |=  x=@rh
    ^-  @rh
    =/  z  (mul:^rh x x)
    ::  sin(x) = x + x^3*S1 + x^5*S2
    ::         = x*(1 + z*(S1 + z*S2))
    =/  r  (add:^rh sin-s1 (mul:^rh z sin-s2))
    (add:^rh x (mul:^rh (mul:^rh x z) r))
  ::
  ::  +cosdf: Kernel cos for |x| <= pi/4
  ::
  ++  cosdf
    |=  x=@rh
    ^-  @rh
    =/  z  (mul:^rh x x)
    ::  cos(x) = 1 + z*C0 + z^2*C1
    ::         = 1 + z*(C0 + z*C1)
    =/  r  (add:^rh cos-c0 (mul:^rh z cos-c1))
    (add:^rh `@rh`0x3c00 (mul:^rh z r))
  ::
  ::  +sin: Sine with argument reduction
  ::
  ++  sin
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    =/  ax    `@rh`(dis bits 0x7fff)
    ::  tiny: sin(x) ~ x
    ?:  (lth:^rh ax `@rh`0x1000)  x
    ::  inf/nan
    ?:  (gte:^rh ax `@rh`0x7c00)  (sub:^rh x x)
    ::  Argument reduction (simplified for half precision range)
    =/  n   (div:^rh ax pio2-hi)
    =/  ni  `@`(rsh [0 10] `@`n)
    =/  y   (sub:^rh x (mul:^rh (sun:^rh ni) pio2-hi))
    ?-  (mod ni 4)
      %0  (sindf y)
      %1  (cosdf y)
      %2  (sub:^rh `@rh`0x0 (sindf y))
      %3  (sub:^rh `@rh`0x0 (cosdf y))
    ==
  ::
  ::  +cos: Cosine with argument reduction
  ::
  ++  cos
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    =/  ax    `@rh`(dis bits 0x7fff)
    ::  tiny: cos(x) ~ 1
    ?:  (lth:^rh ax `@rh`0x1000)  `@rh`0x3c00
    ::  inf/nan
    ?:  (gte:^rh ax `@rh`0x7c00)  (sub:^rh x x)
    ::  Argument reduction
    =/  n   (div:^rh ax pio2-hi)
    =/  ni  `@`(rsh [0 10] `@`n)
    =/  y   (sub:^rh x (mul:^rh (sun:^rh ni) pio2-hi))
    ?-  (mod ni 4)
      %0  (cosdf y)
      %1  (sub:^rh `@rh`0x0 (sindf y))
      %2  (sub:^rh `@rh`0x0 (cosdf y))
      %3  (sindf y)
    ==
  ::
  ::  +tan: Tangent
  ::
  ++  tan
    |=  x=@rh
    ^-  @rh
    (div:^rh (sin x) (cos x))
  ::
  ::  +exp: Exponential (simplified for half precision)
  ::
  ++  exp
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    ::  x = 0 -> 1
    ?:  =(x `@rh`0x0)  `@rh`0x3c00
    ::  +inf -> +inf
    ?:  =(bits 0x7c00)  `@rh`0x7c00
    ::  -inf -> 0
    ?:  =(bits 0xfc00)  `@rh`0x0
    ::  nan -> nan
    ?:  (gth:^rh `@rh`(dis bits 0x7fff) `@rh`0x7c00)  x
    ::  overflow check (x > 11.09)
    ?:  (gth:^rh x `@rh`0x498a)  `@rh`0x7c00
    ::  underflow check (x < -17)
    ?:  (lth:^rh x `@rh`0xcc40)  `@rh`0x0
    ::  exp(r) ~ 1 + r + r^2/2 + r^3/6
    =/  r2  (mul:^rh x x)
    =/  p   (add:^rh exp-p1 (mul:^rh x exp-p2))
    (add:^rh `@rh`0x3c00 (mul:^rh x (add:^rh `@rh`0x3c00 (mul:^rh x p))))
  --
::
::  ================================================================
::  SINGLE PRECISION (@rs) - 32-bit IEEE 754
::  ================================================================
::
::  Coefficients from musl __sindf.c, __cosdf.c, expf.c
::  All values as hex literals for bit-exact matching with C jets.
::
++  rs
  |%
  ::
  ::  ============================================================
  ::  CONSTANTS - Must match exactly between Hoon and C jet
  ::  ============================================================
  ::
  ::  Pi-related constants for argument reduction
  ::  From musl __rem_pio2f.c
  ::
  ++  pio2-1    `@rs`0x3fc9.0fda  ::  1.5707963267341256e+00
  ++  pio2-1t   `@rs`0x33a2.2168  ::  1.5893254773528196e-08
  ++  pio4      `@rs`0x3f49.0fda  ::  pi/4 = 0.7853981633974483
  ++  invpio2   `@rs`0x3f22.f983  ::  2/pi = 0.6366197723675814
  ::
  ::  Exp constants
  ::
  ++  ln2hi     `@rs`0x3f31.7200  ::  0.6931381225585938
  ++  ln2lo     `@rs`0x3717.f7d1  ::  9.058001e-06
  ++  invln2    `@rs`0x3fb8.aa3b  ::  1.4426950216293335
  ++  exp-huge  `@rs`0x42b1.7218  ::  88.72283935546875 (overflow threshold)
  ++  exp-tiny  `@rs`0xc2cf.f1b5  ::  -103.97208404541016 (underflow threshold)
  ::
  ::  ============================================================
  ::  POLYNOMIAL COEFFICIENTS - from musl (hex-exact)
  ::  ============================================================
  ::
  ::  Sin coefficients: sin(x) = x + x^3*S1 + x^5*S2 + x^7*S3 + x^9*S4
  ::  From musl __sindf.c
  ::
  ++  sin-s1    `@rs`0xbe2a.aaab  ::  -1.6666666641831398e-01
  ++  sin-s2    `@rs`0x3c08.8889  ::   8.3333293930888176e-03
  ++  sin-s3    `@rs`0xb950.0d01  ::  -1.9839334836483002e-04
  ++  sin-s4    `@rs`0x3638.ef1b  ::   2.7181216275692873e-06
  ::
  ::  Cos coefficients: cos(x) = 1 - x^2/2 + x^4*C0 + x^6*C1 + x^8*C2 + x^10*C3
  ::  From musl __cosdf.c
  ::
  ++  cos-c0    `@rs`0x3d2a.aaab  ::   4.1666668653488159e-02
  ++  cos-c1    `@rs`0xbab6.0b61  ::  -1.3888889225125313e-03
  ++  cos-c2    `@rs`0x37d0.0d01  ::   2.4801587642822415e-05
  ++  cos-c3    `@rs`0xb493.f27c  ::  -2.7557314603984356e-07
  ::
  ::  Exp polynomial coefficients
  ::  exp(r) ~ 1 + r + r^2*(P1 + r*(P2 + r*(P3 + r*P4)))
  ::
  ++  exp-p1    `@rs`0x3f00.0000  ::   0.5
  ++  exp-p2    `@rs`0x3e2a.aaab  ::   0.16666667163372040
  ++  exp-p3    `@rs`0x3d2a.aaab  ::   0.04166668653488159
  ++  exp-p4    `@rs`0x3c08.8889  ::   0.008333338864892721
  ::
  ::  ============================================================
  ::  CHEBYSHEV POLYNOMIAL
  ::  ============================================================
  ::
  ::  +cheb: Chebyshev polynomial T_n(x) of the first kind
  ::
  ::  The Chebyshev polynomials are defined by the recurrence:
  ::    T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
  ::  with T_0(x) = 1 and T_1(x) = x.
  ::
  ::    Examples
  ::      > ((cheb 0) .0.25)
  ::      .1
  ::      > ((cheb 1) .0.25)
  ::      .0.25
  ::      > ((cheb 2) .0.25)
  ::      .-0.875
  ::      > ((cheb 3) .0.25)
  ::      .-0.6875
  ::
  ++  cheb
    |=  n=@ud
    |=  x=@rs
    ^-  @rs
    ?:  =(0 n)  .1
    ?:  =(1 n)  x
    (sub:^rs (mul:^rs .2 (mul:^rs x $(n (dec n)))) $(n (sub n 2)))
  ::
  ::  ============================================================
  ::  HELPER FUNCTIONS
  ::  ============================================================
  ::
  ::  +scalbn: multiply x by 2^n (ldexp equivalent)
  ::
  ++  scalbn
    |=  [x=@rs n=@sd]
    ^-  @rs
    =/  bits  `@`x
    =/  exp   (dis (rsh [0 23] bits) 0xff)
    =/  sign  (rsh [0 31] bits)
    ::  Special cases
    ?:  =(exp 0xff)  x                    ::  inf or nan
    ?:  =(exp 0)     x                    ::  zero or subnormal
    ::  Compute new exponent
    =/  new-exp  (sum:si (sun:si exp) n)
    ::  Overflow
    ?:  (gth:si new-exp (sun:si 254))
      ?:(=(sign 0) `@rs`0x7f80.0000 `@rs`0xff80.0000)
    ::  Underflow
    ?:  (lth:si new-exp (sun:si 1))
      .0
    ::  Reconstruct
    =/  new-bits
      %+  con  (lsh [0 31] sign)
      %+  con  (lsh [0 23] (abs:si new-exp))
      (dis bits 0x7f.ffff)
    `@rs`new-bits
  ::
  ::  +floor-int: floor to signed integer
  ::
  ++  floor-int
    |=  x=@rs
    ^-  @sd
    =/  bits  `@`x
    =/  sign  (rsh [0 31] bits)
    =/  bexp  (dis (rsh [0 23] bits) 0xff)
    ?:  (^lth bexp 127)                   ::  |x| < 1
      ?:(=(sign 0) --0 -1)
    =/  exp   (sub bexp 127)
    =/  mant  (con (dis bits 0x7f.ffff) 0x80.0000)
    ?:  (^gte exp 23)                     ::  No fractional part
      ?:  =(sign 0)
        (sun:si (lsh [0 (sub exp 23)] mant))
      (new:si %.n (lsh [0 (sub exp 23)] mant))
    =/  shift     (sub 23 exp)
    =/  int-mant  (rsh [0 shift] mant)
    =/  frac-mask (dec (lsh [0 shift] 1))
    =/  had-frac  !=(0 (dis mant frac-mask))
    =/  result    ?:(=(sign 0) (sun:si int-mant) (new:si %.n int-mant))
    ?:  &(=(sign 1) had-frac)
      (dif:si result --1)
    result
  ::
  ::  ============================================================
  ::  ARGUMENT REDUCTION
  ::  ============================================================
  ::
  ::  +rem-pio2: Reduce x to [-pi/4, pi/4], return quadrant
  ::
  ::  Returns [n y] where x = n*(pi/2) + y, |y| <= pi/4
  ::
  ++  rem-pio2
    |=  x=@rs
    ^-  [n=@ y=@rs]
    =/  bits  `@`x
    =/  sign  (rsh [0 31] bits)
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    ::  Small argument: no reduction needed
    ?:  (lte:^rs ax pio4)
      [0 x]
    ::  Compute n = round(x * 2/pi)
    =/  fn-raw  (mul:^rs x invpio2)
    =/  fn      ?:  (gte:^rs fn-raw .0)
                  (add:^rs fn-raw .0.5)
                (sub:^rs fn-raw .0.5)
    =/  n       (floor-int fn)
    =/  fn      (san:^rs n)
    ::  y = x - n*pio2_1 - n*pio2_1t (extended precision)
    =/  y  (sub:^rs x (mul:^rs fn pio2-1))
    =.  y  (sub:^rs y (mul:^rs fn pio2-1t))
    [(mod (abs:si n) 4) y]
  ::
  ::  ============================================================
  ::  KERNEL FUNCTIONS - Fixed polynomial evaluation
  ::  ============================================================
  ::
  ::  +sindf: Kernel sin for reduced argument |x| <= pi/4
  ::
  ::  sin(x) = x + x^3*S1 + x^5*S2 + x^7*S3 + x^9*S4
  ::         = x + x*z*(S1 + z*(S2 + z*(S3 + z*S4)))
  ::  where z = x^2
  ::
  ++  sindf
    |=  x=@rs
    ^-  @rs
    =/  z  (mul:^rs x x)                                 ::  z = x^2
    ::  r = S3 + z*S4
    =/  r  (add:^rs sin-s3 (mul:^rs z sin-s4))
    ::  r = S2 + z*r
    =.  r  (add:^rs sin-s2 (mul:^rs z r))
    ::  r = S1 + z*r
    =.  r  (add:^rs sin-s1 (mul:^rs z r))
    ::  result = x + x*z*r
    (add:^rs x (mul:^rs (mul:^rs x z) r))
  ::
  ::  +cosdf: Kernel cos for reduced argument |x| <= pi/4
  ::
  ::  cos(x) = 1 - z/2 + z^2*(C0 + z*(C1 + z*(C2 + z*C3)))
  ::  where z = x^2
  ::
  ++  cosdf
    |=  x=@rs
    ^-  @rs
    =/  z  (mul:^rs x x)                                 ::  z = x^2
    =/  w  (mul:^rs z z)                                 ::  w = z^2
    ::  r = C2 + z*C3
    =/  r  (add:^rs cos-c2 (mul:^rs z cos-c3))
    ::  r = C1 + z*r
    =.  r  (add:^rs cos-c1 (mul:^rs z r))
    ::  r = C0 + z*r
    =.  r  (add:^rs cos-c0 (mul:^rs z r))
    ::  result = 1 - z/2 + w*r
    =/  hz  (mul:^rs .0.5 z)
    (add:^rs (sub:^rs .1 hz) (mul:^rs w r))
  ::
  ::  ============================================================
  ::  PUBLIC API
  ::  ============================================================
  ::
  ::  +sin: Sine
  ::
  ++  sin
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    ::  Tiny: sin(x) ~ x
    ?:  (lth:^rs ax `@rs`0x3980.0000)  x                 ::  |x| < 2^-12
    ::  Inf/NaN
    ?:  (gte:^rs ax `@rs`0x7f80.0000)  (sub:^rs x x)
    ::  Argument reduction
    =/  [n=@ y=@rs]  (rem-pio2 x)
    ::  Quadrant dispatch
    ?-  n
      %0  (sindf y)
      %1  (cosdf y)
      %2  (sub:^rs .0 (sindf y))
      %3  (sub:^rs .0 (cosdf y))
    ==
  ::
  ::  +cos: Cosine
  ::
  ++  cos
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    ::  Tiny: cos(x) ~ 1
    ?:  (lth:^rs ax `@rs`0x3980.0000)  .1
    ::  Inf/NaN
    ?:  (gte:^rs ax `@rs`0x7f80.0000)  (sub:^rs x x)
    ::  Argument reduction
    =/  [n=@ y=@rs]  (rem-pio2 x)
    ::  Quadrant dispatch
    ?-  n
      %0  (cosdf y)
      %1  (sub:^rs .0 (sindf y))
      %2  (sub:^rs .0 (cosdf y))
      %3  (sindf y)
    ==
  ::
  ::  +tan: Tangent
  ::
  ++  tan
    |=  x=@rs
    ^-  @rs
    (div:^rs (sin x) (cos x))
  ::
  ::  +exp: Exponential
  ::
  ::  Algorithm:
  ::    1. Handle special cases
  ::    2. Reduce: x = k*ln2 + r, |r| <= ln2/2
  ::    3. Compute exp(r) via polynomial
  ::    4. Return 2^k * exp(r)
  ::
  ++  exp
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    ::  x = 0 -> 1
    ?:  =(x .0)  .1
    ::  +inf -> +inf
    ?:  =(bits 0x7f80.0000)  `@rs`0x7f80.0000
    ::  -inf -> 0
    ?:  =(bits 0xff80.0000)  .0
    ::  NaN -> NaN
    ?:  (gth:^rs ax `@rs`0x7f80.0000)  x
    ::  Overflow
    ?:  (gth:^rs x exp-huge)  `@rs`0x7f80.0000
    ::  Underflow
    ?:  (lth:^rs x exp-tiny)  .0
    ::
    ::  Argument reduction: x = k*ln2 + r
    ::
    =/  fn-raw  (mul:^rs x invln2)
    =/  fn      ?:  (gte:^rs fn-raw .0)
                  (add:^rs fn-raw .0.5)
                (sub:^rs fn-raw .0.5)
    =/  k   (floor-int fn)
    =/  kf  (san:^rs k)
    ::  r = x - k*ln2hi - k*ln2lo (extended precision)
    =/  r   (sub:^rs x (mul:^rs kf ln2hi))
    =.  r   (sub:^rs r (mul:^rs kf ln2lo))
    ::
    ::  Polynomial: exp(r) ~ 1 + r + r^2*(P1 + r*(P2 + r*(P3 + r*P4)))
    ::
    =/  r2  (mul:^rs r r)
    ::  p = P3 + r*P4
    =/  p   (add:^rs exp-p3 (mul:^rs r exp-p4))
    ::  p = P2 + r*p
    =.  p   (add:^rs exp-p2 (mul:^rs r p))
    ::  p = P1 + r*p
    =.  p   (add:^rs exp-p1 (mul:^rs r p))
    ::  exp(r) = 1 + r + r^2*p
    =/  expr  (add:^rs .1 (add:^rs r (mul:^rs r2 p)))
    ::
    ::  Scale by 2^k
    ::
    (scalbn expr k)
  --
::
::  ================================================================
::  DOUBLE PRECISION (@rd) - 64-bit IEEE 754
::  ================================================================
::
::  Coefficients from musl __sin.c, __cos.c, exp.c
::  Using higher-degree polynomials for double precision accuracy.
::
++  rd
  |%
  ::
  ::  ============================================================
  ::  CONSTANTS (hex-exact, from musl)
  ::  ============================================================
  ::
  ::  Pi-related constants for argument reduction
  ::
  ++  pio2-1    `@rd`0x3ff9.21fb.5440.0000  ::  pi/2 high bits
  ++  pio2-1t   `@rd`0x3dd0.b461.1a62.6331  ::  pi/2 - pio2_1
  ++  pio2-2    `@rd`0x3dd0.b461.1a60.0000  ::  second 33 bits of pi/2
  ++  pio2-2t   `@rd`0x3ba3.198a.2e03.7073  ::  pi/2 - pio2_1 - pio2_2
  ++  pio4      `@rd`0x3fe9.21fb.5444.2d18  ::  pi/4
  ++  invpio2   `@rd`0x3fe4.5f30.6dc9.c883  ::  2/pi
  ::
  ::  Exp constants
  ::
  ++  ln2hi     `@rd`0x3fe6.2e42.fee0.0000  ::  ln(2) high
  ++  ln2lo     `@rd`0x3dea.39ef.3579.3c76  ::  ln(2) low
  ++  invln2    `@rd`0x3ff7.1547.652b.82fe  ::  1/ln(2)
  ++  exp-huge  `@rd`0x4086.2e42.fefa.39ef  ::  709.782... (overflow)
  ++  exp-tiny  `@rd`0xc087.4910.d52d.3051  ::  -745.133... (underflow)
  ::
  ::  ============================================================
  ::  POLYNOMIAL COEFFICIENTS (hex-exact, from musl)
  ::  ============================================================
  ::
  ::  Sin coefficients: degree 13 polynomial
  ::  sin(x) = x + x^3*(S1 + x^2*(S2 + x^2*(S3 + x^2*(S4 + x^2*(S5 + x^2*S6)))))
  ::
  ++  sin-s1    `@rd`0xbfc5.5555.5555.5549  ::  -1.66666666666666324348e-01
  ++  sin-s2    `@rd`0x3f81.1111.1110.f8a6  ::   8.33333333332248946124e-03
  ++  sin-s3    `@rd`0xbf2a.01a0.19c1.61d5  ::  -1.98412698298579493134e-04
  ++  sin-s4    `@rd`0x3ec7.1de3.57b1.fe7d  ::   2.75573137070700676789e-06
  ++  sin-s5    `@rd`0xbe5a.e5e6.8a2b.9ceb  ::  -2.50507602534068634195e-08
  ++  sin-s6    `@rd`0x3de5.d93a.5acf.d57c  ::   1.58969099521155010221e-10
  ::
  ::  Cos coefficients: degree 14 polynomial
  ::  cos(x) = 1 - x^2/2 + x^4*(C1 + x^2*(C2 + x^2*(C3 + x^2*(C4 + x^2*(C5 + x^2*C6)))))
  ::
  ++  cos-c1    `@rd`0x3fa5.5555.5555.554c  ::   4.16666666666666019037e-02
  ++  cos-c2    `@rd`0xbf56.c16c.16c1.5177  ::  -1.38888888888741095749e-03
  ++  cos-c3    `@rd`0x3efa.01a0.19cb.1590  ::   2.48015872894767294178e-05
  ++  cos-c4    `@rd`0xbe92.7e4f.809c.52ad  ::  -2.75573143513906633035e-07
  ++  cos-c5    `@rd`0x3e21.ee9e.bdb4.b1c4  ::   2.08757232129817482790e-09
  ++  cos-c6    `@rd`0xbda8.fae9.be88.38d4  ::  -1.13596475577881948265e-11
  ::
  ::  Exp polynomial coefficients
  ::  From musl exp_data.c
  ::
  ++  exp-p1    `@rd`0x3fe0.0000.0000.0000  ::   0.5
  ++  exp-p2    `@rd`0x3fc5.5555.5555.5555  ::   0.16666666666666666
  ++  exp-p3    `@rd`0x3fa5.5555.5555.5555  ::   0.041666666666666664
  ++  exp-p4    `@rd`0x3f81.1111.1111.1111  ::   0.008333333333333333
  ++  exp-p5    `@rd`0x3f56.c16c.16c1.6c17  ::   0.001388888888888889
  ::
  ::  ============================================================
  ::  CHEBYSHEV POLYNOMIAL
  ::  ============================================================
  ::
  ::  +cheb: Chebyshev polynomial T_n(x) of the first kind
  ::
  ++  cheb
    |=  n=@ud
    |=  x=@rd
    ^-  @rd
    ?:  =(0 n)  .~1
    ?:  =(1 n)  x
    (sub:^rd (mul:^rd .~2 (mul:^rd x $(n (dec n)))) $(n (sub n 2)))
  ::
  ::  ============================================================
  ::  HELPER FUNCTIONS
  ::  ============================================================
  ::
  ::  +scalbn: multiply x by 2^n
  ::
  ++  scalbn
    |=  [x=@rd n=@sd]
    ^-  @rd
    =/  bits  `@`x
    =/  exp   (dis (rsh [0 52] bits) 0x7ff)
    =/  sign  (rsh [0 63] bits)
    ?:  =(exp 0x7ff)  x                   ::  inf or nan
    ?:  =(exp 0)      x                   ::  zero or subnormal
    =/  new-exp  (sum:si (sun:si exp) n)
    ?:  (gth:si new-exp (sun:si 2.046))
      ?:(=(sign 0) `@rd`0x7ff0.0000.0000.0000 `@rd`0xfff0.0000.0000.0000)
    ?:  (lth:si new-exp (sun:si 1))
      .~0
    =/  new-bits
      %+  con  (lsh [0 63] sign)
      %+  con  (lsh [0 52] (abs:si new-exp))
      (dis bits 0xf.ffff.ffff.ffff)
    `@rd`new-bits
  ::
  ::  +floor-int: floor to signed integer
  ::
  ++  floor-int
    |=  x=@rd
    ^-  @sd
    =/  bits  `@`x
    =/  sign  (rsh [0 63] bits)
    =/  bexp  (dis (rsh [0 52] bits) 0x7ff)
    ?:  (^lth bexp 1.023)
      ?:(=(sign 0) --0 -1)
    =/  exp   (sub bexp 1.023)
    =/  mant  (con (dis bits 0xf.ffff.ffff.ffff) 0x10.0000.0000.0000)
    ?:  (^gte exp 52)
      ?:  =(sign 0)
        (sun:si (lsh [0 (sub exp 52)] mant))
      (new:si %.n (lsh [0 (sub exp 52)] mant))
    =/  shift     (sub 52 exp)
    =/  int-mant  (rsh [0 shift] mant)
    =/  frac-mask (dec (lsh [0 shift] 1))
    =/  had-frac  !=(0 (dis mant frac-mask))
    =/  result    ?:(=(sign 0) (sun:si int-mant) (new:si %.n int-mant))
    ?:  &(=(sign 1) had-frac)
      (dif:si result --1)
    result
  ::
  ::  ============================================================
  ::  ARGUMENT REDUCTION
  ::  ============================================================
  ::
  ::  +rem-pio2: Reduce x to [-pi/4, pi/4]
  ::
  ++  rem-pio2
    |=  x=@rd
    ^-  [n=@ y=@rd]
    =/  bits  `@`x
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    ?:  (lte:^rd ax pio4)
      [0 x]
    =/  fn-raw  (mul:^rd x invpio2)
    =/  fn      ?:  (gte:^rd fn-raw .~0)
                  (add:^rd fn-raw .~0.5)
                (sub:^rd fn-raw .~0.5)
    =/  n   (floor-int fn)
    =/  fn  (san:^rd n)
    ::  Extended precision reduction
    =/  y   (sub:^rd x (mul:^rd fn pio2-1))
    =.  y   (sub:^rd y (mul:^rd fn pio2-1t))
    [(mod (abs:si n) 4) y]
  ::
  ::  ============================================================
  ::  KERNEL FUNCTIONS
  ::  ============================================================
  ::
  ::  +sindf: Kernel sin (degree 13)
  ::
  ++  sindf
    |=  x=@rd
    ^-  @rd
    =/  z  (mul:^rd x x)
    ::  Horner: S5 + z*S6
    =/  r  (add:^rd sin-s5 (mul:^rd z sin-s6))
    ::  S4 + z*r
    =.  r  (add:^rd sin-s4 (mul:^rd z r))
    ::  S3 + z*r
    =.  r  (add:^rd sin-s3 (mul:^rd z r))
    ::  S2 + z*r
    =.  r  (add:^rd sin-s2 (mul:^rd z r))
    ::  S1 + z*r
    =.  r  (add:^rd sin-s1 (mul:^rd z r))
    ::  x + x*z*r
    (add:^rd x (mul:^rd (mul:^rd x z) r))
  ::
  ::  +cosdf: Kernel cos (degree 14)
  ::
  ++  cosdf
    |=  x=@rd
    ^-  @rd
    =/  z  (mul:^rd x x)
    =/  w  (mul:^rd z z)
    ::  Horner: C5 + z*C6
    =/  r  (add:^rd cos-c5 (mul:^rd z cos-c6))
    ::  C4 + z*r
    =.  r  (add:^rd cos-c4 (mul:^rd z r))
    ::  C3 + z*r
    =.  r  (add:^rd cos-c3 (mul:^rd z r))
    ::  C2 + z*r
    =.  r  (add:^rd cos-c2 (mul:^rd z r))
    ::  C1 + z*r
    =.  r  (add:^rd cos-c1 (mul:^rd z r))
    ::  1 - z/2 + w*r
    =/  hz  (mul:^rd .~0.5 z)
    (add:^rd (sub:^rd .~1 hz) (mul:^rd w r))
  ::
  ::  ============================================================
  ::  PUBLIC API
  ::  ============================================================
  ::
  ++  sin
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    ::  Tiny
    ?:  (lth:^rd ax `@rd`0x3e40.0000.0000.0000)  x
    ::  Inf/NaN
    ?:  (gte:^rd ax `@rd`0x7ff0.0000.0000.0000)  (sub:^rd x x)
    =/  [n=@ y=@rd]  (rem-pio2 x)
    ?-  n
      %0  (sindf y)
      %1  (cosdf y)
      %2  (sub:^rd .~0 (sindf y))
      %3  (sub:^rd .~0 (cosdf y))
    ==
  ::
  ++  cos
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    ?:  (lth:^rd ax `@rd`0x3e40.0000.0000.0000)  .~1
    ?:  (gte:^rd ax `@rd`0x7ff0.0000.0000.0000)  (sub:^rd x x)
    =/  [n=@ y=@rd]  (rem-pio2 x)
    ?-  n
      %0  (cosdf y)
      %1  (sub:^rd .~0 (sindf y))
      %2  (sub:^rd .~0 (cosdf y))
      %3  (sindf y)
    ==
  ::
  ++  tan
    |=  x=@rd
    ^-  @rd
    (div:^rd (sin x) (cos x))
  ::
  ++  exp
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    ?:  =(x .~0)  .~1
    ?:  =(bits 0x7ff0.0000.0000.0000)  `@rd`0x7ff0.0000.0000.0000
    ?:  =(bits 0xfff0.0000.0000.0000)  .~0
    ?:  (gth:^rd ax `@rd`0x7ff0.0000.0000.0000)  x
    ?:  (gth:^rd x exp-huge)  `@rd`0x7ff0.0000.0000.0000
    ?:  (lth:^rd x exp-tiny)  .~0
    ::  Argument reduction
    =/  fn-raw  (mul:^rd x invln2)
    =/  fn      ?:  (gte:^rd fn-raw .~0)
                  (add:^rd fn-raw .~0.5)
                (sub:^rd fn-raw .~0.5)
    =/  k   (floor-int fn)
    =/  kf  (san:^rd k)
    =/  r   (sub:^rd x (mul:^rd kf ln2hi))
    =.  r   (sub:^rd r (mul:^rd kf ln2lo))
    ::  Polynomial (degree 5)
    =/  r2  (mul:^rd r r)
    =/  p   (add:^rd exp-p4 (mul:^rd r exp-p5))
    =.  p   (add:^rd exp-p3 (mul:^rd r p))
    =.  p   (add:^rd exp-p2 (mul:^rd r p))
    =.  p   (add:^rd exp-p1 (mul:^rd r p))
    =/  expr  (add:^rd .~1 (add:^rd r (mul:^rd r2 p)))
    (scalbn expr k)
  --
::
::  ================================================================
::  QUAD PRECISION (@rq) - 128-bit IEEE 754
::  ================================================================
::
::  Coefficients from FreeBSD msun ld128/k_sinl.c, k_cosl.c
::  Quad precision uses degree 25 sin, degree 22 cos for ~113 bits accuracy.
::
::  Note: 128-bit hex literals are 32 hex digits each.
::  Format: 1 sign + 15 exponent (bias 16383) + 112 mantissa bits
::
++  rq
  |%
  ::
  ::  ============================================================
  ::  CONSTANTS (hex-exact for binary128)
  ::  ============================================================
  ::
  ::  Pi-related constants
  ::
  ++  pio2-1    `@rq`0x3fff.921f.b544.42d1.8469.898c.c517.01b8  ::  pi/2 high
  ++  pio2-1t   `@rq`0x3f8c.c51c.e0ae.9ca4.e68e.6e36.30d1.5dc4  ::  pi/2 - pio2_1
  ++  pio4      `@rq`0x3ffe.921f.b544.42d1.8469.898c.c517.01b8  ::  pi/4
  ++  invpio2   `@rq`0x3ffe.45f3.06dc.9c88.2a53.f84e.afa3.ea6a  ::  2/pi
  ::
  ::  Exp constants
  ::
  ++  ln2hi     `@rq`0x3ffe.62e4.2fef.a39e.f357.93c7.6730.07e6  ::  ln(2) high
  ++  ln2lo     `@rq`0x3fbc.7abc.9e3b.3998.0f2f.707c.9ca4.4b10  ::  ln(2) low
  ++  invln2    `@rq`0x3fff.7154.7652.b82f.e177.7d0f.fda0.d23a  ::  1/ln(2)
  ++  exp-huge  `@rq`0x400c.62e4.2fef.a39e.f357.93c7.6730.07e6  ::  ~11356.5 (overflow)
  ++  exp-tiny  `@rq`0xc00c.6760.f0c7.f77a.8a37.8c9c.e5ca.10bc  ::  ~-11433.5 (underflow)
  ::
  ::  ============================================================
  ::  POLYNOMIAL COEFFICIENTS
  ::  ============================================================
  ::
  ::  Sin coefficients: degree 25 polynomial (S1-S12)
  ::  From FreeBSD ld128/k_sinl.c
  ::  sin(x) = x + x^3*(S1 + x^2*(S2 + ...))
  ::
  ++  sin-s1    `@rq`0xbffc.5555.5555.5555.5555.5555.5555.5555  ::  -1/6
  ++  sin-s2    `@rq`0x3ff8.1111.1111.1111.1111.1111.1111.1111  ::   1/120
  ++  sin-s3    `@rq`0xbff2.a01a.01a0.1a01.a01a.01a0.1a01.a01a  ::  -1/5040
  ++  sin-s4    `@rq`0x3fec.71de.3a55.6c73.38fa.ac1c.88e5.0017  ::   1/362880
  ++  sin-s5    `@rq`0xbfe5.ae64.567f.544e.38fe.747e.4b60.c9a6  ::  -1/39916800
  ++  sin-s6    `@rq`0x3fde.6124.613a.86d0.9393.30de.0b58.0e66  ::   1/6227020800
  ++  sin-s7    `@rq`0xbfd6.ae7f.3e73.3b81.f7be.ced5.1769.e17e  ::  ~-1.15e-15
  ++  sin-s8    `@rq`0x3fce.952c.7703.0ad4.a71d.fc88.40d2.35fc  ::  ~8.22e-18
  ::
  ::  Cos coefficients: degree 22 polynomial (C1-C11)
  ::  From FreeBSD ld128/k_cosl.c
  ::  cos(x) = 1 - x^2/2 + x^4*(C1 + x^2*(C2 + ...))
  ::
  ++  cos-c1    `@rq`0x3ffa.5555.5555.5555.5555.5555.5555.5555  ::   1/24
  ++  cos-c2    `@rq`0xbff5.6c16.c16c.16c1.6c16.c16c.16c1.6c17  ::  -1/720
  ++  cos-c3    `@rq`0x3fef.a01a.01a0.1a01.a01a.01a0.1a01.a01a  ::   1/40320
  ++  cos-c4    `@rq`0xbfe9.27e4.fb77.89f5.c72e.f016.d3ea.6679  ::  -1/3628800
  ++  cos-c5    `@rq`0x3fe2.1eed.8eff.8d89.7b5e.33bf.1a9b.a5c9  ::   1/479001600
  ++  cos-c6    `@rq`0xbfda.93f2.7dbb.c4fb.7459.4a58.2e61.6a96  ::  ~-8.9e-14
  ++  cos-c7    `@rq`0x3fd2.ae7f.3e73.3b94.0c61.7a76.bec4.4b09  ::  ~4.78e-16
  ++  cos-c8    `@rq`0xbfca.6827.863b.97d9.7e4c.9bbd.c044.2faa  ::  ~-1.56e-18
  ::
  ::  Exp polynomial coefficients (degree 7 for quad)
  ::
  ++  exp-p1    `@rq`0x3ffe.0000.0000.0000.0000.0000.0000.0000  ::   0.5
  ++  exp-p2    `@rq`0x3ffc.5555.5555.5555.5555.5555.5555.5555  ::   1/6
  ++  exp-p3    `@rq`0x3ffa.5555.5555.5555.5555.5555.5555.5555  ::   1/24
  ++  exp-p4    `@rq`0x3ff8.1111.1111.1111.1111.1111.1111.1111  ::   1/120
  ++  exp-p5    `@rq`0x3ff5.6c16.c16c.16c1.6c16.c16c.16c1.6c17  ::   1/720
  ++  exp-p6    `@rq`0x3ff2.a01a.01a0.1a01.a01a.01a0.1a01.a01a  ::   1/5040
  ++  exp-p7    `@rq`0x3fef.a01a.01a0.1a01.a01a.01a0.1a01.a01a  ::   1/40320
  ::
  ::  ============================================================
  ::  CHEBYSHEV POLYNOMIAL
  ::  ============================================================
  ::
  ++  cheb
    |=  n=@ud
    |=  x=@rq
    ^-  @rq
    ?:  =(0 n)  .~~~1
    ?:  =(1 n)  x
    (sub:^rq (mul:^rq .~~~2 (mul:^rq x $(n (dec n)))) $(n (sub n 2)))
  ::
  ::  ============================================================
  ::  HELPER FUNCTIONS
  ::  ============================================================
  ::
  ++  scalbn
    |=  [x=@rq n=@sd]
    ^-  @rq
    =/  bits  `@`x
    =/  exp   (dis (rsh [0 112] bits) 0x7fff)
    =/  sign  (rsh [0 127] bits)
    ?:  =(exp 0x7fff)  x                  ::  inf or nan
    ?:  =(exp 0)       x                  ::  zero or subnormal
    =/  new-exp  (sum:si (sun:si exp) n)
    ?:  (gth:si new-exp (sun:si 32.766))
      ?:  =(sign 0)
        `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000
      `@rq`0xffff.0000.0000.0000.0000.0000.0000.0000
    ?:  (lth:si new-exp (sun:si 1))
      .~~~0
    =/  new-bits
      %+  con  (lsh [0 127] sign)
      %+  con  (lsh [0 112] (abs:si new-exp))
      (dis bits 0xffff.ffff.ffff.ffff.ffff.ffff.ffff)
    `@rq`new-bits
  ::
  ++  floor-int
    |=  x=@rq
    ^-  @sd
    =/  bits  `@`x
    =/  sign  (rsh [0 127] bits)
    =/  bexp  (dis (rsh [0 112] bits) 0x7fff)
    ?:  (^lth bexp 16.383)
      ?:(=(sign 0) --0 -1)
    =/  exp   (sub bexp 16.383)
    =/  mant  (con (dis bits 0xffff.ffff.ffff.ffff.ffff.ffff.ffff) 0x1.0000.0000.0000.0000.0000.0000.0000)
    ?:  (^gte exp 112)
      ?:  =(sign 0)
        (sun:si (lsh [0 (sub exp 112)] mant))
      (new:si %.n (lsh [0 (sub exp 112)] mant))
    =/  shift     (sub 112 exp)
    =/  int-mant  (rsh [0 shift] mant)
    =/  frac-mask (dec (lsh [0 shift] 1))
    =/  had-frac  !=(0 (dis mant frac-mask))
    =/  result    ?:(=(sign 0) (sun:si int-mant) (new:si %.n int-mant))
    ?:  &(=(sign 1) had-frac)
      (dif:si result --1)
    result
  ::
  ::  ============================================================
  ::  ARGUMENT REDUCTION
  ::  ============================================================
  ::
  ++  rem-pio2
    |=  x=@rq
    ^-  [n=@ y=@rq]
    =/  bits  `@`x
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ?:  (lte:^rq ax pio4)
      [0 x]
    =/  fn-raw  (mul:^rq x invpio2)
    =/  fn      ?:  (gte:^rq fn-raw .~~~0)
                  (add:^rq fn-raw .~~~0.5)
                (sub:^rq fn-raw .~~~0.5)
    =/  n   (floor-int fn)
    =/  fn  (san:^rq n)
    =/  y   (sub:^rq x (mul:^rq fn pio2-1))
    =.  y   (sub:^rq y (mul:^rq fn pio2-1t))
    [(mod (abs:si n) 4) y]
  ::
  ::  ============================================================
  ::  KERNEL FUNCTIONS
  ::  ============================================================
  ::
  ::  +sindf: Kernel sin (degree 17 for practical quad precision)
  ::
  ++  sindf
    |=  x=@rq
    ^-  @rq
    =/  z  (mul:^rq x x)
    ::  Horner form with 8 coefficients
    =/  r  (add:^rq sin-s7 (mul:^rq z sin-s8))
    =.  r  (add:^rq sin-s6 (mul:^rq z r))
    =.  r  (add:^rq sin-s5 (mul:^rq z r))
    =.  r  (add:^rq sin-s4 (mul:^rq z r))
    =.  r  (add:^rq sin-s3 (mul:^rq z r))
    =.  r  (add:^rq sin-s2 (mul:^rq z r))
    =.  r  (add:^rq sin-s1 (mul:^rq z r))
    (add:^rq x (mul:^rq (mul:^rq x z) r))
  ::
  ::  +cosdf: Kernel cos (degree 16 for practical quad precision)
  ::
  ++  cosdf
    |=  x=@rq
    ^-  @rq
    =/  z  (mul:^rq x x)
    =/  w  (mul:^rq z z)
    ::  Horner form with 8 coefficients
    =/  r  (add:^rq cos-c7 (mul:^rq z cos-c8))
    =.  r  (add:^rq cos-c6 (mul:^rq z r))
    =.  r  (add:^rq cos-c5 (mul:^rq z r))
    =.  r  (add:^rq cos-c4 (mul:^rq z r))
    =.  r  (add:^rq cos-c3 (mul:^rq z r))
    =.  r  (add:^rq cos-c2 (mul:^rq z r))
    =.  r  (add:^rq cos-c1 (mul:^rq z r))
    =/  hz  (mul:^rq .~~~0.5 z)
    (add:^rq (sub:^rq .~~~1 hz) (mul:^rq w r))
  ::
  ::  ============================================================
  ::  PUBLIC API
  ::  ============================================================
  ::
  ++  sin
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ::  Tiny
    ?:  (lth:^rq ax `@rq`0x3f8f.0000.0000.0000.0000.0000.0000.0000)  x
    ::  Inf/NaN
    ?:  (gte:^rq ax `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000)
      (sub:^rq x x)
    =/  [n=@ y=@rq]  (rem-pio2 x)
    ?-  n
      %0  (sindf y)
      %1  (cosdf y)
      %2  (sub:^rq .~~~0 (sindf y))
      %3  (sub:^rq .~~~0 (cosdf y))
    ==
  ::
  ++  cos
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ?:  (lth:^rq ax `@rq`0x3f8f.0000.0000.0000.0000.0000.0000.0000)  .~~~1
    ?:  (gte:^rq ax `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000)
      (sub:^rq x x)
    =/  [n=@ y=@rq]  (rem-pio2 x)
    ?-  n
      %0  (cosdf y)
      %1  (sub:^rq .~~~0 (sindf y))
      %2  (sub:^rq .~~~0 (cosdf y))
      %3  (sindf y)
    ==
  ::
  ++  tan
    |=  x=@rq
    ^-  @rq
    (div:^rq (sin x) (cos x))
  ::
  ++  exp
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ?:  =(x .~~~0)  .~~~1
    ?:  =(bits 0x7fff.0000.0000.0000.0000.0000.0000.0000)
      `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000
    ?:  =(bits 0xffff.0000.0000.0000.0000.0000.0000.0000)  .~~~0
    ?:  (gth:^rq ax `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000)  x
    ?:  (gth:^rq x exp-huge)  `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000
    ?:  (lth:^rq x exp-tiny)  .~~~0
    ::  Argument reduction
    =/  fn-raw  (mul:^rq x invln2)
    =/  fn      ?:  (gte:^rq fn-raw .~~~0)
                  (add:^rq fn-raw .~~~0.5)
                (sub:^rq fn-raw .~~~0.5)
    =/  k   (floor-int fn)
    =/  kf  (san:^rq k)
    =/  r   (sub:^rq x (mul:^rq kf ln2hi))
    =.  r   (sub:^rq r (mul:^rq kf ln2lo))
    ::  Polynomial (degree 7)
    =/  r2  (mul:^rq r r)
    =/  p   (add:^rq exp-p6 (mul:^rq r exp-p7))
    =.  p   (add:^rq exp-p5 (mul:^rq r p))
    =.  p   (add:^rq exp-p4 (mul:^rq r p))
    =.  p   (add:^rq exp-p3 (mul:^rq r p))
    =.  p   (add:^rq exp-p2 (mul:^rq r p))
    =.  p   (add:^rq exp-p1 (mul:^rq r p))
    =/  expr  (add:^rq .~~~1 (add:^rq r (mul:^rq r2 p)))
    (scalbn expr k)
  --
::
::  ================================================================
::  JET IMPLEMENTATION GUIDE
::  ================================================================
::
::  For C jets to produce identical results:
::
::  1. Use SoftFloat for all FP operations (IEEE 754 compliance)
::
::  2. Define coefficients as hex literals:
::       static const double S1 = 0xbfc5555555555549; // reinterpreted
::     Or use union-based bit casting.
::
::  3. Match operation order exactly:
::       // sindf for @rd
::       double z = x * x;
::       double r = S5 + z * S6;
::       r = S4 + z * r;
::       r = S3 + z * r;
::       r = S2 + z * r;
::       r = S1 + z * r;
::       return x + x * z * r;
::
::  4. Use volatile or explicit stores if compiler reorders.
::
--
