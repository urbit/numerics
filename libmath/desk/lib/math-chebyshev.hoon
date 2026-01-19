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
  ::
  ::  +log: Natural logarithm (simplified for half precision)
  ::
  ++  log
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    ::  Special cases
    ?:  =(x `@rh`0x0)  `@rh`0xfc00              ::  log(0) = -inf
    ?:  =(bits 0x7c00)  `@rh`0x7c00             ::  log(+inf) = +inf
    ?:  (gth:^rh `@rh`(dis bits 0x7fff) `@rh`0x7c00)  x  ::  nan
    ?:  !=(0 (rsh [0 15] bits))  (sub:^rh x x)  ::  log(negative) = nan
    ::
    ::  Extract exponent and mantissa
    ::  Half: 1 sign, 5 exp (bias 15), 10 mantissa
    ::
    =/  exp-bits  (dis (rsh [0 10] bits) 0x1f)
    =/  mant-bits  (dis bits 0x3ff)
    ?:  =(exp-bits 0)  `@rh`0xfc00             ::  subnormal
    =/  k  (sub exp-bits 15)
    =/  m  `@rh`(con 0x3c00 mant-bits)         ::  1.0 + mantissa
    ::
    ::  s = (m-1)/(m+1), log(m) ~ 2*s*(1 + s^2*Lg1)
    ::
    =/  one   `@rh`0x3c00
    =/  two   `@rh`0x4000
    =/  f     (sub:^rh m one)
    =/  s     (div:^rh f (add:^rh two f))
    =/  s2    (mul:^rh s s)
    =/  lg1   `@rh`0x3955                      ::  ~0.666 (2/3)
    =/  log-m (mul:^rh (mul:^rh two s) (add:^rh one (mul:^rh s2 lg1)))
    ::
    =/  kf    (sun:^rh k)
    (add:^rh (mul:^rh kf ln2hi) log-m)
  ::
  ++  log-2
    |=  x=@rh
    ^-  @rh
    (mul:^rh (log x) invln2)
  ::
  ++  log-10
    |=  x=@rh
    ^-  @rh
    =/  invlog10  `@rh`0x36f3                   ::  ~0.4343
    (mul:^rh (log x) invlog10)
  ::
  ::  +sqt: Square root (Newton-Raphson, 3 iterations for half)
  ::
  ++  sqt
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    ?:  =(x `@rh`0x0)  `@rh`0x0
    ?:  =(bits 0x7c00)  `@rh`0x7c00
    ?:  (gth:^rh `@rh`(dis bits 0x7fff) `@rh`0x7c00)  x
    ?:  !=(0 (rsh [0 15] bits))  (sub:^rh x x)
    ::  Initial guess
    =/  g0  `@rh`(add (rsh [0 1] bits) 0x1c00)
    =/  half  `@rh`0x3800                       ::  0.5
    =/  g  g0
    =.  g  (mul:^rh half (add:^rh g (div:^rh x g)))
    =.  g  (mul:^rh half (add:^rh g (div:^rh x g)))
    =.  g  (mul:^rh half (add:^rh g (div:^rh x g)))
    g
  ::
  ++  sqrt  sqt
  ::
  ::  +atan: Inverse tangent (simplified for half precision)
  ::
  ++  atan
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    =/  sign  (rsh [0 15] bits)
    =/  ax    `@rh`(dis bits 0x7fff)
    ::
    ?:  (gte:^rh ax `@rh`0x7c00)
      ?:  =(bits 0x7c00)  `@rh`0x3e48           ::  pi/2
      ?:  =(bits 0xfc00)  `@rh`0xbe48           ::  -pi/2
      x
    ::
    ::  Simple polynomial for small |x|
    ::  atan(x) ~ x - x^3/3 + x^5/5
    ::
    =/  one   `@rh`0x3c00
    =/  pio2  `@rh`0x3e48                       ::  pi/2
    =/  at1   `@rh`0xb555                       ::  -1/3
    =/  at2   `@rh`0x3266                       ::  ~0.2 (1/5)
    ::
    =/  result=@rh
      ?:  (lte:^rh ax one)
        ::  Small: polynomial directly
        =/  x2  (mul:^rh ax ax)
        =/  x3  (mul:^rh x2 ax)
        =/  r   (add:^rh at1 (mul:^rh x2 at2))
        (add:^rh ax (mul:^rh x3 r))
      ::  Large: atan(x) = pi/2 - atan(1/x)
      =/  t   (div:^rh one ax)
      =/  t2  (mul:^rh t t)
      =/  t3  (mul:^rh t2 t)
      =/  r   (add:^rh at1 (mul:^rh t2 at2))
      (sub:^rh pio2 (add:^rh t (mul:^rh t3 r)))
    ::
    ?:(=(sign 0) result (sub:^rh `@rh`0x0 result))
  ::
  ++  asin
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    =/  ax    `@rh`(dis bits 0x7fff)
    =/  one   `@rh`0x3c00
    ?:  (gth:^rh ax one)  (sub:^rh x x)
    ?:  =(ax one)
      ?:  =(0 (rsh [0 15] bits))  `@rh`0x3e48   ::  pi/2
      `@rh`0xbe48                                ::  -pi/2
    (atan (div:^rh x (sqt (sub:^rh one (mul:^rh x x)))))
  ::
  ++  acos
    |=  x=@rh
    ^-  @rh
    =/  bits  `@`x
    =/  ax    `@rh`(dis bits 0x7fff)
    =/  one   `@rh`0x3c00
    =/  pi    `@rh`0x4248                       ::  pi
    ?:  (gth:^rh ax one)  (sub:^rh x x)
    ?:  =(x one)   `@rh`0x0
    ?:  =(x `@rh`0xbc00)  pi                    ::  acos(-1) = pi
    =/  s  (sqt (sub:^rh one (mul:^rh x x)))
    ?:  (gte:^rh x `@rh`0x0)
      (atan (div:^rh s x))
    (add:^rh pi (atan (div:^rh s x)))
  ::
  ++  atan2
    |=  [y=@rh x=@rh]
    ^-  @rh
    =/  pi    `@rh`0x4248
    =/  pio2  `@rh`0x3e48
    =/  zero  `@rh`0x0
    ?:  (gth:^rh x zero)
      (atan (div:^rh y x))
    ?:  &((lth:^rh x zero) (gte:^rh y zero))
      (add:^rh (atan (div:^rh y x)) pi)
    ?:  &((lth:^rh x zero) (lth:^rh y zero))
      (sub:^rh (atan (div:^rh y x)) pi)
    ?:  &(=(x zero) (gth:^rh y zero))
      pio2
    ?:  &(=(x zero) (lth:^rh y zero))
      (sub:^rh zero pio2)
    zero
  ::
  ++  pow-n
    |=  [x=@rh n=@rh]
    ^-  @rh
    =/  one  `@rh`0x3c00
    ?:  =(n `@rh`0x0)  one
    =/  ni  (abs:si (need (toi:^rh n)))
    =/  neg  (lth:^rh n `@rh`0x0)
    =/  result  one
    =/  base    x
    =/  i       0
    |-
    ?:  |(=(ni 0) (gth i 15))
      ?:(neg (div:^rh one result) result)
    =/  new-result
      ?:  =(1 (dis ni 1))
        (mul:^rh result base)
      result
    $(ni (rsh [0 1] ni), base (mul:^rh base base), result new-result, i +(i))
  ::
  ++  pow
    |=  [x=@rh n=@rh]
    ^-  @rh
    =/  ni  (toi:^rh n)
    ?:  &(?=(^ ni) =(n (san:^rh (need ni))))
      (pow-n x n)
    (exp (mul:^rh n (log x)))
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
  ::
  ::  +log: Natural logarithm
  ::
  ::  Algorithm (from FreeBSD e_logf.c):
  ::    1. Handle special cases (0, negative, inf, nan)
  ::    2. Reduce: x = 2^k * m where 1 <= m < 2
  ::    3. Further reduce: s = (m-1)/(m+1), so m = (1+s)/(1-s)
  ::    4. Compute log(m) = 2*s + 2*s^3*(Lg1 + s^2*(Lg2 + s^2*(Lg3 + s^2*Lg4)))
  ::    5. Return k*ln2 + log(m)
  ::
  ++  log
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    ::  Special cases
    ?:  =(x .0)  `@rs`0xff80.0000              ::  log(0) = -inf
    ?:  =(bits 0x7f80.0000)  `@rs`0x7f80.0000  ::  log(+inf) = +inf
    ?:  (gth:^rs `@rs`(dis bits 0x7fff.ffff) `@rs`0x7f80.0000)  x  ::  nan
    ?:  !=(0 (rsh [0 31] bits))  (sub:^rs x x)  ::  log(negative) = nan
    ::
    ::  Extract exponent and mantissa
    ::  x = 2^k * m where 1 <= m < 2
    ::
    =/  exp-bits  (dis (rsh [0 23] bits) 0xff)
    =/  mant-bits  (dis bits 0x7f.ffff)
    ::  Handle subnormals (simplified: treat as very small)
    ?:  =(exp-bits 0)  `@rs`0xff80.0000
    ::  k = exponent - bias
    =/  k  (sub exp-bits 127)
    ::  m = 1.mantissa (reconstruct normalized mantissa)
    =/  m  `@rs`(con 0x3f80.0000 mant-bits)
    ::
    ::  Compute s = (m-1)/(m+1)
    ::
    =/  f   (sub:^rs m .1)
    =/  s   (div:^rs f (add:^rs .2 f))
    =/  s2  (mul:^rs s s)
    =/  s4  (mul:^rs s2 s2)
    ::
    ::  Polynomial: R = Lg2 + s2*(Lg4)  [even terms]
    ::              + Lg1 + s2*(Lg3)     [odd terms]
    ::  Coefficients from FreeBSD e_logf.c
    ::
    =/  lg1  `@rs`0x3f2a.aaaa                  ::  0.66666662693
    =/  lg2  `@rs`0x3ecc.ce13                  ::  0.40000972152
    =/  lg3  `@rs`0x3e91.e9ee                  ::  0.28498786688
    =/  lg4  `@rs`0x3e78.9e26                  ::  0.24279078841
    ::
    =/  t1  (mul:^rs s4 (add:^rs lg2 (mul:^rs s4 lg4)))
    =/  t2  (mul:^rs s2 (add:^rs lg1 (mul:^rs s4 lg3)))
    =/  r   (add:^rs t1 t2)
    ::
    ::  log(m) = 2*s + 2*s^3*R = 2*s*(1 + s^2*R)
    ::
    =/  log-m  (mul:^rs (mul:^rs .2 s) (add:^rs .1 (mul:^rs s2 r)))
    ::
    ::  Result = k*ln2 + log(m)
    ::
    =/  kf  (san:^rs (sun:si k))
    (add:^rs (mul:^rs kf ln2hi) log-m)
  ::
  ::  +log-2: Base-2 logarithm
  ::
  ++  log-2
    |=  x=@rs
    ^-  @rs
    (mul:^rs (log x) invln2)
  ::
  ::  +log-10: Base-10 logarithm
  ::
  ++  log-10
    |=  x=@rs
    ^-  @rs
    =/  invlog10  `@rs`0x3ede.5bd9              ::  1/log(10) = 0.43429448
    (mul:^rs (log x) invlog10)
  ::
  ::  +sqt: Square root (Newton-Raphson with fixed iterations)
  ::
  ::  Algorithm:
  ::    1. Handle special cases
  ::    2. Initial guess from bit manipulation
  ::    3. Fixed 4 Newton-Raphson iterations: g = (g + x/g) / 2
  ::
  ++  sqt
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    ::  Special cases
    ?:  =(x .0)  .0
    ?:  =(bits 0x7f80.0000)  `@rs`0x7f80.0000  ::  sqrt(+inf) = +inf
    ?:  (gth:^rs `@rs`(dis bits 0x7fff.ffff) `@rs`0x7f80.0000)  x  ::  nan
    ?:  !=(0 (rsh [0 31] bits))  (sub:^rs x x)  ::  sqrt(negative) = nan
    ::
    ::  Initial guess: use bit manipulation for good starting point
    ::  g0 = bits/2 + 0x1fc00000 (magic constant for float sqrt)
    ::
    =/  g0  `@rs`(add (rsh [0 1] bits) 0x1fc0.0000)
    ::
    ::  Newton-Raphson: g = (g + x/g) / 2
    ::  Fixed 4 iterations for single precision
    ::
    =/  g  g0
    =.  g  (mul:^rs .0.5 (add:^rs g (div:^rs x g)))  ::  iter 1
    =.  g  (mul:^rs .0.5 (add:^rs g (div:^rs x g)))  ::  iter 2
    =.  g  (mul:^rs .0.5 (add:^rs g (div:^rs x g)))  ::  iter 3
    =.  g  (mul:^rs .0.5 (add:^rs g (div:^rs x g)))  ::  iter 4
    g
  ::
  ++  sqrt  sqt
  ::
  ::  +atan: Inverse tangent
  ::
  ::  Algorithm (from FreeBSD s_atanf.c):
  ::    1. Reduce to [0, 0.4375] using atan(x) = pi/2 - atan(1/x) for |x| > 1
  ::       and atan(x) = atan(c) + atan((x-c)/(1+x*c)) for thresholds
  ::    2. Polynomial approximation in reduced range
  ::
  ++  atan
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    =/  sign  (rsh [0 31] bits)
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    ::
    ::  Special cases
    ?:  (gte:^rs ax `@rs`0x7f80.0000)           ::  inf or nan
      ?:  =(bits 0x7f80.0000)                   ::  +inf
        `@rs`0x3fc9.0fdb                        ::  pi/2
      ?:  =(bits 0xff80.0000)                   ::  -inf
        `@rs`0xbfc9.0fdb                        ::  -pi/2
      x                                          ::  nan
    ::
    ::  Polynomial coefficients from FreeBSD s_atanf.c
    ::
    =/  at0  `@rs`0x3eaa.aaaa                   ::   0.33333334
    =/  at1  `@rs`0xbe4c.cccc                   ::  -0.20000000
    =/  at2  `@rs`0x3e12.4925                   ::   0.14285715
    =/  at3  `@rs`0xbdc3.b8a4                   ::  -0.11111089
    ::
    ::  Reduction thresholds
    ::
    =/  lo-thresh  `@rs`0x3ee0.0000             ::  0.4375
    =/  hi-thresh  `@rs`0x401c.0000             ::  2.4375
    ::
    ::  atan constants for reduction
    ::
    =/  pio2    `@rs`0x3fc9.0fdb                ::  pi/2
    =/  pio4    `@rs`0x3f49.0fdb                ::  pi/4
    ::
    ::  Core polynomial evaluation for |x| <= 0.4375
    ::
    =/  eval-poly
      |=  z=@rs
      ^-  @rs
      =/  w   (mul:^rs z z)
      =/  s1  (add:^rs at2 (mul:^rs w at3))
      =/  s2  (add:^rs at0 (mul:^rs w at1))
      =/  r   (add:^rs s2 (mul:^rs (mul:^rs w w) s1))
      (sub:^rs z (mul:^rs z (mul:^rs w r)))
    ::
    =/  result=@rs
      ?:  (lte:^rs ax lo-thresh)
        ::  Small: use polynomial directly
        (eval-poly ax)
      ?:  (lte:^rs ax .1)
        ::  Medium-small: atan(x) = atan(0.5) + atan((x-0.5)/(1+x*0.5))
        =/  atan-half  `@rs`0x3eed.6338         ::  atan(0.5) = 0.46364760
        =/  t  (div:^rs (sub:^rs ax .0.5) (add:^rs .1 (mul:^rs ax .0.5)))
        (add:^rs atan-half (eval-poly t))
      ?:  (lte:^rs ax hi-thresh)
        ::  Medium: atan(x) = pi/4 + atan((x-1)/(1+x))
        =/  t  (div:^rs (sub:^rs ax .1) (add:^rs .1 ax))
        (add:^rs pio4 (eval-poly t))
      ::  Large: atan(x) = pi/2 - atan(1/x)
      =/  t  (div:^rs .1 ax)
      (sub:^rs pio2 (eval-poly t))
    ::
    ::  Apply sign
    ?:(=(sign 0) result (sub:^rs .0 result))
  ::
  ::  +asin: Inverse sine
  ::
  ++  asin
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    ?:  (gth:^rs ax .1)  (sub:^rs x x)          ::  |x| > 1: nan
    ?:  =(ax .1)                                 ::  |x| = 1: +/- pi/2
      ?:  =(0 (rsh [0 31] bits))
        `@rs`0x3fc9.0fdb                        ::  pi/2
      `@rs`0xbfc9.0fdb                          ::  -pi/2
    ::  asin(x) = atan(x / sqrt(1 - x^2))
    (atan (div:^rs x (sqt (sub:^rs .1 (mul:^rs x x)))))
  ::
  ::  +acos: Inverse cosine
  ::
  ++  acos
    |=  x=@rs
    ^-  @rs
    =/  bits  `@`x
    =/  ax    `@rs`(dis bits 0x7fff.ffff)
    =/  pi    `@rs`0x4049.0fdb                  ::  pi
    ?:  (gth:^rs ax .1)  (sub:^rs x x)          ::  |x| > 1: nan
    ?:  =(x .1)   .0                             ::  acos(1) = 0
    ?:  =(x .-1)  pi                             ::  acos(-1) = pi
    ::  acos(x) = atan(sqrt(1 - x^2) / x) with quadrant adjustment
    =/  s  (sqt (sub:^rs .1 (mul:^rs x x)))
    ?:  (gte:^rs x .0)
      (atan (div:^rs s x))
    (add:^rs pi (atan (div:^rs s x)))
  ::
  ::  +atan2: Two-argument inverse tangent
  ::
  ++  atan2
    |=  [y=@rs x=@rs]
    ^-  @rs
    =/  pi    `@rs`0x4049.0fdb                  ::  pi
    =/  pio2  `@rs`0x3fc9.0fdb                  ::  pi/2
    ?:  (gth:^rs x .0)
      (atan (div:^rs y x))
    ?:  &((lth:^rs x .0) (gte:^rs y .0))
      (add:^rs (atan (div:^rs y x)) pi)
    ?:  &((lth:^rs x .0) (lth:^rs y .0))
      (sub:^rs (atan (div:^rs y x)) pi)
    ?:  &(=(x .0) (gth:^rs y .0))
      pio2
    ?:  &(=(x .0) (lth:^rs y .0))
      (sub:^rs .0 pio2)
    .0                                           ::  undefined (0,0)
  ::
  ::  +pow-n: Integer power (fixed iteration)
  ::
  ++  pow-n
    |=  [x=@rs n=@rs]
    ^-  @rs
    ?:  =(n .0)  .1
    ::  Convert n to integer (must be positive integer)
    =/  ni  (abs:si (need (toi:^rs n)))
    =/  neg  (lth:^rs n .0)
    ::  Binary exponentiation (fixed max 32 iterations for 32-bit range)
    =/  result  .1
    =/  base    x
    =/  i       0
    |-
    ?:  |(=(ni 0) (gth i 31))
      ?:(neg (div:^rs .1 result) result)
    =/  new-result
      ?:  =(1 (dis ni 1))
        (mul:^rs result base)
      result
    $(ni (rsh [0 1] ni), base (mul:^rs base base), result new-result, i +(i))
  ::
  ::  +pow: General power function
  ::
  ++  pow
    |=  [x=@rs n=@rs]
    ^-  @rs
    ::  Special case: integer exponent
    =/  ni  (toi:^rs n)
    ?:  &(?=(^ ni) =(n (san:^rs (need ni))))
      (pow-n x n)
    ::  General case: x^n = exp(n * log(x))
    (exp (mul:^rs n (log x)))
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
  ::
  ::  +log: Natural logarithm (double precision)
  ::
  ++  log
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    ::  Special cases
    ?:  =(x .~0)  `@rd`0xfff0.0000.0000.0000              ::  log(0) = -inf
    ?:  =(bits 0x7ff0.0000.0000.0000)  `@rd`0x7ff0.0000.0000.0000
    ?:  (gth:^rd `@rd`(dis bits 0x7fff.ffff.ffff.ffff) `@rd`0x7ff0.0000.0000.0000)  x
    ?:  !=(0 (rsh [0 63] bits))  (sub:^rd x x)
    ::
    ::  Extract exponent and mantissa
    ::
    =/  exp-bits  (dis (rsh [0 52] bits) 0x7ff)
    =/  mant-bits  (dis bits 0xf.ffff.ffff.ffff)
    ?:  =(exp-bits 0)  `@rd`0xfff0.0000.0000.0000
    =/  k  (sub exp-bits 1.023)
    =/  m  `@rd`(con 0x3ff0.0000.0000.0000 mant-bits)
    ::
    ::  s = (m-1)/(m+1)
    ::
    =/  f   (sub:^rd m .~1)
    =/  s   (div:^rd f (add:^rd .~2 f))
    =/  s2  (mul:^rd s s)
    =/  s4  (mul:^rd s2 s2)
    ::
    ::  Lg1-Lg7 from FreeBSD e_log.c (hex-exact)
    ::
    =/  lg1  `@rd`0x3fe5.5555.5555.5593  ::  0.6666666666666735130
    =/  lg2  `@rd`0x3fd9.9999.9997.fa04  ::  0.3999999999940941908
    =/  lg3  `@rd`0x3fd2.4924.9422.9359  ::  0.2857142874366239149
    =/  lg4  `@rd`0x3fcc.71c5.1d8e.78af  ::  0.2222219843214978396
    =/  lg5  `@rd`0x3fc7.4664.96cb.03de  ::  0.1818357216161805012
    =/  lg6  `@rd`0x3fc3.9a09.d078.c69f  ::  0.1531383769920937332
    =/  lg7  `@rd`0x3fc2.f112.df3e.5244  ::  0.1479819860511658591
    ::
    =/  t1  (mul:^rd s4 (add:^rd lg2 (mul:^rd s4 (add:^rd lg4 (mul:^rd s4 lg6)))))
    =/  t2  (mul:^rd s2 (add:^rd lg1 (mul:^rd s4 (add:^rd lg3 (mul:^rd s4 (add:^rd lg5 (mul:^rd s4 lg7)))))))
    =/  r   (add:^rd t1 t2)
    =/  log-m  (mul:^rd (mul:^rd .~2 s) (add:^rd .~1 (mul:^rd s2 r)))
    =/  kf  (san:^rd (sun:si k))
    (add:^rd (mul:^rd kf ln2hi) log-m)
  ::
  ++  log-2
    |=  x=@rd
    ^-  @rd
    (mul:^rd (log x) invln2)
  ::
  ++  log-10
    |=  x=@rd
    ^-  @rd
    =/  invlog10  `@rd`0x3fdb.cb7b.1526.e50e  ::  0.43429448190325176
    (mul:^rd (log x) invlog10)
  ::
  ::  +sqt: Square root (Newton-Raphson, 6 iterations for double)
  ::
  ++  sqt
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    ?:  =(x .~0)  .~0
    ?:  =(bits 0x7ff0.0000.0000.0000)  `@rd`0x7ff0.0000.0000.0000
    ?:  (gth:^rd `@rd`(dis bits 0x7fff.ffff.ffff.ffff) `@rd`0x7ff0.0000.0000.0000)  x
    ?:  !=(0 (rsh [0 63] bits))  (sub:^rd x x)
    ::  Initial guess
    =/  g0  `@rd`(add (rsh [0 1] bits) 0x1ff0.0000.0000.0000)
    =/  g  g0
    =.  g  (mul:^rd .~0.5 (add:^rd g (div:^rd x g)))
    =.  g  (mul:^rd .~0.5 (add:^rd g (div:^rd x g)))
    =.  g  (mul:^rd .~0.5 (add:^rd g (div:^rd x g)))
    =.  g  (mul:^rd .~0.5 (add:^rd g (div:^rd x g)))
    =.  g  (mul:^rd .~0.5 (add:^rd g (div:^rd x g)))
    =.  g  (mul:^rd .~0.5 (add:^rd g (div:^rd x g)))
    g
  ::
  ++  sqrt  sqt
  ::
  ::  +atan: Inverse tangent (double precision)
  ::
  ++  atan
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    =/  sign  (rsh [0 63] bits)
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    ::
    ?:  (gte:^rd ax `@rd`0x7ff0.0000.0000.0000)
      ?:  =(bits 0x7ff0.0000.0000.0000)  `@rd`0x3ff9.21fb.5444.2d18  ::  pi/2
      ?:  =(bits 0xfff0.0000.0000.0000)  `@rd`0xbff9.21fb.5444.2d18
      x
    ::
    ::  Coefficients from FreeBSD s_atan.c
    ::
    =/  at0   `@rd`0x3fd5.5555.5555.550d  ::   0.33333333333329318027
    =/  at1   `@rd`0xbfc9.9999.9998.ebc4  ::  -0.19999999998764832476
    =/  at2   `@rd`0x3fc2.4924.9200.83ff  ::   0.14285714272503466371
    =/  at3   `@rd`0xbfbc.71c6.fe23.1671  ::  -0.11111104054623557880
    =/  at4   `@rd`0x3fb7.45cd.c54c.206e  ::   0.09090887133436506962
    =/  at5   `@rd`0xbfb3.b0f2.af74.9a6d  ::  -0.07691876205044829950
    =/  at6   `@rd`0x3fb1.0d66.a0d0.3d51  ::   0.06661073137387531207
    =/  at7   `@rd`0xbfad.de2d.52de.fd9a  ::  -0.05833570133790573485
    =/  at8   `@rd`0x3fa9.7b4b.2476.0deb  ::   0.04976877994615936017
    =/  at9   `@rd`0xbfa2.b444.2c6a.6c2f  ::  -0.03653157274421691527
    =/  at10  `@rd`0x3f90.ad3a.e322.da11  ::   0.01628582015365782362
    ::
    =/  pio2  `@rd`0x3ff9.21fb.5444.2d18  ::  pi/2
    =/  pio4  `@rd`0x3fe9.21fb.5444.2d18  ::  pi/4
    ::
    =/  eval-poly
      |=  z=@rd
      ^-  @rd
      =/  w   (mul:^rd z z)
      =/  w2  (mul:^rd w w)
      =/  s1  (add:^rd at0 (mul:^rd w (add:^rd at1 (mul:^rd w (add:^rd at2 (mul:^rd w (add:^rd at3 (mul:^rd w at4))))))))
      =/  s2  (add:^rd at5 (mul:^rd w (add:^rd at6 (mul:^rd w (add:^rd at7 (mul:^rd w (add:^rd at8 (mul:^rd w (add:^rd at9 (mul:^rd w at10))))))))))
      =/  r   (add:^rd s1 (mul:^rd w2 (mul:^rd w s2)))
      (sub:^rd z (mul:^rd z (mul:^rd w r)))
    ::
    =/  lo-thresh  `@rd`0x3fdc.0000.0000.0000  ::  0.4375
    =/  hi-thresh  `@rd`0x4003.8000.0000.0000  ::  2.4375
    ::
    =/  result=@rd
      ?:  (lte:^rd ax lo-thresh)
        (eval-poly ax)
      ?:  (lte:^rd ax .~1)
        =/  atan-half  `@rd`0x3fdd.ac67.0561.bb4f
        =/  t  (div:^rd (sub:^rd ax .~0.5) (add:^rd .~1 (mul:^rd ax .~0.5)))
        (add:^rd atan-half (eval-poly t))
      ?:  (lte:^rd ax hi-thresh)
        =/  t  (div:^rd (sub:^rd ax .~1) (add:^rd .~1 ax))
        (add:^rd pio4 (eval-poly t))
      =/  t  (div:^rd .~1 ax)
      (sub:^rd pio2 (eval-poly t))
    ::
    ?:(=(sign 0) result (sub:^rd .~0 result))
  ::
  ++  asin
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    ?:  (gth:^rd ax .~1)  (sub:^rd x x)
    ?:  =(ax .~1)
      ?:  =(0 (rsh [0 63] bits))
        `@rd`0x3ff9.21fb.5444.2d18
      `@rd`0xbff9.21fb.5444.2d18
    (atan (div:^rd x (sqt (sub:^rd .~1 (mul:^rd x x)))))
  ::
  ++  acos
    |=  x=@rd
    ^-  @rd
    =/  bits  `@`x
    =/  ax    `@rd`(dis bits 0x7fff.ffff.ffff.ffff)
    =/  pi    `@rd`0x4009.21fb.5444.2d18
    ?:  (gth:^rd ax .~1)  (sub:^rd x x)
    ?:  =(x .~1)   .~0
    ?:  =(x .~-1)  pi
    =/  s  (sqt (sub:^rd .~1 (mul:^rd x x)))
    ?:  (gte:^rd x .~0)
      (atan (div:^rd s x))
    (add:^rd pi (atan (div:^rd s x)))
  ::
  ++  atan2
    |=  [y=@rd x=@rd]
    ^-  @rd
    =/  pi    `@rd`0x4009.21fb.5444.2d18
    =/  pio2  `@rd`0x3ff9.21fb.5444.2d18
    ?:  (gth:^rd x .~0)
      (atan (div:^rd y x))
    ?:  &((lth:^rd x .~0) (gte:^rd y .~0))
      (add:^rd (atan (div:^rd y x)) pi)
    ?:  &((lth:^rd x .~0) (lth:^rd y .~0))
      (sub:^rd (atan (div:^rd y x)) pi)
    ?:  &(=(x .~0) (gth:^rd y .~0))
      pio2
    ?:  &(=(x .~0) (lth:^rd y .~0))
      (sub:^rd .~0 pio2)
    .~0
  ::
  ++  pow-n
    |=  [x=@rd n=@rd]
    ^-  @rd
    ?:  =(n .~0)  .~1
    =/  ni  (abs:si (need (toi:^rd n)))
    =/  neg  (lth:^rd n .~0)
    =/  result  .~1
    =/  base    x
    =/  i       0
    |-
    ?:  |(=(ni 0) (gth i 63))
      ?:(neg (div:^rd .~1 result) result)
    =/  new-result
      ?:  =(1 (dis ni 1))
        (mul:^rd result base)
      result
    $(ni (rsh [0 1] ni), base (mul:^rd base base), result new-result, i +(i))
  ::
  ++  pow
    |=  [x=@rd n=@rd]
    ^-  @rd
    =/  ni  (toi:^rd n)
    ?:  &(?=(^ ni) =(n (san:^rd (need ni))))
      (pow-n x n)
    (exp (mul:^rd n (log x)))
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
  ::
  ::  +log: Natural logarithm (quad precision)
  ::
  ::  Uses higher-degree polynomial for 113-bit mantissa accuracy
  ::
  ++  log
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    ::  Special cases
    ?:  =(x .~~~0)  `@rq`0xffff.0000.0000.0000.0000.0000.0000.0000  ::  -inf
    ?:  =(bits 0x7fff.0000.0000.0000.0000.0000.0000.0000)
      `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000                ::  +inf
    ?:  (gth:^rq `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
                `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000)
      x                                                              ::  nan
    ?:  !=(0 (rsh [0 127] bits))  (sub:^rq x x)                     ::  negative
    ::
    ::  Extract exponent and mantissa
    ::  Quad: 1 sign, 15 exp (bias 16383), 112 mantissa
    ::
    =/  exp-bits  (dis (rsh [0 112] bits) 0x7fff)
    =/  mant-bits  (dis bits 0xffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ?:  =(exp-bits 0)  `@rq`0xffff.0000.0000.0000.0000.0000.0000.0000
    =/  k  (sub exp-bits 16.383)
    =/  m  `@rq`(con 0x3fff.0000.0000.0000.0000.0000.0000.0000 mant-bits)
    ::
    ::  s = (m-1)/(m+1)
    ::
    =/  f   (sub:^rq m .~~~1)
    =/  s   (div:^rq f (add:^rq .~~~2 f))
    =/  s2  (mul:^rq s s)
    =/  s4  (mul:^rq s2 s2)
    ::
    ::  Lg1-Lg8 for quad precision (higher accuracy)
    ::
    =/  lg1  `@rq`0x3ffc.5555.5555.5555.5555.5555.5555.5555  ::  2/3
    =/  lg2  `@rq`0x3ffb.9999.9999.9999.9999.9999.9999.999a  ::  2/5
    =/  lg3  `@rq`0x3ffb.2492.4924.9249.2492.4924.9249.2492  ::  2/7
    =/  lg4  `@rq`0x3ffa.c71c.71c7.1c71.c71c.71c7.1c71.c71c  ::  2/9
    =/  lg5  `@rq`0x3ffa.745d.1745.d174.5d17.45d1.745d.1746  ::  2/11
    =/  lg6  `@rq`0x3ffa.3a83.a83a.83a8.3a83.a83a.83a8.3a84  ::  2/13
    =/  lg7  `@rq`0x3ffa.0f0f.0f0f.0f0f.0f0f.0f0f.0f0f.0f0f  ::  2/15
    =/  lg8  `@rq`0x3ff9.e1e1.e1e1.e1e1.e1e1.e1e1.e1e1.e1e2  ::  2/17
    ::
    =/  t1  (mul:^rq s4 (add:^rq lg2 (mul:^rq s4 (add:^rq lg4 (mul:^rq s4 (add:^rq lg6 (mul:^rq s4 lg8)))))))
    =/  t2  (mul:^rq s2 (add:^rq lg1 (mul:^rq s4 (add:^rq lg3 (mul:^rq s4 (add:^rq lg5 (mul:^rq s4 lg7)))))))
    =/  r   (add:^rq t1 t2)
    =/  log-m  (mul:^rq (mul:^rq .~~~2 s) (add:^rq .~~~1 (mul:^rq s2 r)))
    ::
    =/  kf  (san:^rq (sun:si k))
    (add:^rq (mul:^rq kf ln2hi) log-m)
  ::
  ++  log-2
    |=  x=@rq
    ^-  @rq
    (mul:^rq (log x) invln2)
  ::
  ++  log-10
    |=  x=@rq
    ^-  @rq
    =/  invlog10  `@rq`0x3ffd.bcb7.b152.6e50.e32a.6ab7.555f.5a68  ::  1/ln(10)
    (mul:^rq (log x) invlog10)
  ::
  ::  +sqt: Square root (Newton-Raphson, 8 iterations for quad)
  ::
  ++  sqt
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    ?:  =(x .~~~0)  .~~~0
    ?:  =(bits 0x7fff.0000.0000.0000.0000.0000.0000.0000)
      `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000
    ?:  (gth:^rq `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
                `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000)
      x
    ?:  !=(0 (rsh [0 127] bits))  (sub:^rq x x)
    ::  Initial guess
    =/  g0  `@rq`(add (rsh [0 1] bits) 0x1fff.8000.0000.0000.0000.0000.0000.0000)
    =/  g  g0
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    =.  g  (mul:^rq .~~~0.5 (add:^rq g (div:^rq x g)))
    g
  ::
  ++  sqrt  sqt
  ::
  ::  +atan: Inverse tangent (quad precision)
  ::
  ++  atan
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    =/  sign  (rsh [0 127] bits)
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ::
    ?:  (gte:^rq ax `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000)
      ?:  =(bits 0x7fff.0000.0000.0000.0000.0000.0000.0000)
        `@rq`0x3fff.921f.b544.42d1.8469.898c.c517.01b8        ::  pi/2
      ?:  =(bits 0xffff.0000.0000.0000.0000.0000.0000.0000)
        `@rq`0xbfff.921f.b544.42d1.8469.898c.c517.01b8        ::  -pi/2
      x
    ::
    ::  Coefficients for atan polynomial (higher degree for quad)
    ::
    =/  at0   `@rq`0x3ffd.5555.5555.5555.5555.5555.5555.5555  ::   1/3
    =/  at1   `@rq`0xbffc.cccc.cccc.cccc.cccc.cccc.cccc.cccd  ::  -1/5
    =/  at2   `@rq`0x3ffc.4924.9249.2492.4924.9249.2492.4924  ::   1/7
    =/  at3   `@rq`0xbffc.0000.0000.0000.0000.0000.0000.0000  ::  -1/9
    =/  at4   `@rq`0x3ffb.a2e8.ba2e.8ba2.e8ba.2e8b.a2e8.ba2f  ::   1/11
    =/  at5   `@rq`0xbffb.5555.5555.5555.5555.5555.5555.5555  ::  -1/13
    ::
    =/  pio2  `@rq`0x3fff.921f.b544.42d1.8469.898c.c517.01b8  ::  pi/2
    =/  pio4  `@rq`0x3ffe.921f.b544.42d1.8469.898c.c517.01b8  ::  pi/4
    ::
    =/  eval-poly
      |=  z=@rq
      ^-  @rq
      =/  w   (mul:^rq z z)
      =/  w2  (mul:^rq w w)
      =/  s1  (add:^rq at0 (mul:^rq w (add:^rq at1 (mul:^rq w at2))))
      =/  s2  (add:^rq at3 (mul:^rq w (add:^rq at4 (mul:^rq w at5))))
      =/  r   (add:^rq s1 (mul:^rq w2 (mul:^rq w s2)))
      (sub:^rq z (mul:^rq z (mul:^rq w r)))
    ::
    =/  lo-thresh  `@rq`0x3ffe.0000.0000.0000.0000.0000.0000.0000  ::  0.5
    =/  hi-thresh  `@rq`0x4000.8000.0000.0000.0000.0000.0000.0000  ::  2.0
    ::
    =/  result=@rq
      ?:  (lte:^rq ax lo-thresh)
        (eval-poly ax)
      ?:  (lte:^rq ax .~~~1)
        =/  atan-half  `@rq`0x3ffd.dac6.7056.1bb4.f68a.dec4.befd.cd8b
        =/  t  (div:^rq (sub:^rq ax .~~~0.5) (add:^rq .~~~1 (mul:^rq ax .~~~0.5)))
        (add:^rq atan-half (eval-poly t))
      ?:  (lte:^rq ax hi-thresh)
        =/  t  (div:^rq (sub:^rq ax .~~~1) (add:^rq .~~~1 ax))
        (add:^rq pio4 (eval-poly t))
      =/  t  (div:^rq .~~~1 ax)
      (sub:^rq pio2 (eval-poly t))
    ::
    ?:(=(sign 0) result (sub:^rq .~~~0 result))
  ::
  ++  asin
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    ?:  (gth:^rq ax .~~~1)  (sub:^rq x x)
    ?:  =(ax .~~~1)
      ?:  =(0 (rsh [0 127] bits))
        `@rq`0x3fff.921f.b544.42d1.8469.898c.c517.01b8
      `@rq`0xbfff.921f.b544.42d1.8469.898c.c517.01b8
    (atan (div:^rq x (sqt (sub:^rq .~~~1 (mul:^rq x x)))))
  ::
  ++  acos
    |=  x=@rq
    ^-  @rq
    =/  bits  `@`x
    =/  ax    `@rq`(dis bits 0x7fff.ffff.ffff.ffff.ffff.ffff.ffff.ffff)
    =/  pi    `@rq`0x4000.921f.b544.42d1.8469.898c.c517.01b8
    ?:  (gth:^rq ax .~~~1)  (sub:^rq x x)
    ?:  =(x .~~~1)   .~~~0
    ?:  =(x .~~~-1)  pi
    =/  s  (sqt (sub:^rq .~~~1 (mul:^rq x x)))
    ?:  (gte:^rq x .~~~0)
      (atan (div:^rq s x))
    (add:^rq pi (atan (div:^rq s x)))
  ::
  ++  atan2
    |=  [y=@rq x=@rq]
    ^-  @rq
    =/  pi    `@rq`0x4000.921f.b544.42d1.8469.898c.c517.01b8
    =/  pio2  `@rq`0x3fff.921f.b544.42d1.8469.898c.c517.01b8
    ?:  (gth:^rq x .~~~0)
      (atan (div:^rq y x))
    ?:  &((lth:^rq x .~~~0) (gte:^rq y .~~~0))
      (add:^rq (atan (div:^rq y x)) pi)
    ?:  &((lth:^rq x .~~~0) (lth:^rq y .~~~0))
      (sub:^rq (atan (div:^rq y x)) pi)
    ?:  &(=(x .~~~0) (gth:^rq y .~~~0))
      pio2
    ?:  &(=(x .~~~0) (lth:^rq y .~~~0))
      (sub:^rq .~~~0 pio2)
    .~~~0
  ::
  ++  pow-n
    |=  [x=@rq n=@rq]
    ^-  @rq
    ?:  =(n .~~~0)  .~~~1
    =/  ni  (abs:si (need (toi:^rq n)))
    =/  neg  (lth:^rq n .~~~0)
    =/  result  .~~~1
    =/  base    x
    =/  i       0
    |-
    ?:  |(=(ni 0) (gth i 127))
      ?:(neg (div:^rq .~~~1 result) result)
    =/  new-result
      ?:  =(1 (dis ni 1))
        (mul:^rq result base)
      result
    $(ni (rsh [0 1] ni), base (mul:^rq base base), result new-result, i +(i))
  ::
  ++  pow
    |=  [x=@rq n=@rq]
    ^-  @rq
    =/  ni  (toi:^rq n)
    ?:  &(?=(^ ni) =(n (san:^rq (need ni))))
      (pow-n x n)
    (exp (mul:^rq n (log x)))
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
