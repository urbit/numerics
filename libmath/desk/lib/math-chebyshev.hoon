::
::::  Chebyshev/Fixed-Polynomial Math Library Prototype
::
::  This file demonstrates the approach for jet-compatible transcendental
::  functions. The key principle is that every floating-point operation
::  must happen in the same order with the same intermediate values in
::  both Hoon and the C jet.
::
::  Reference: musl libc (MIT license)
::  - https://git.musl-libc.org/cgit/musl/tree/src/math/
::
::  For jetting, the C code should use the exact same algorithm with
::  the exact same coefficients. SoftFloat should be used for IEEE 754
::  compliance in intermediate operations.
::
|%
::
::  Single-precision (@rs) implementations
::
++  rs
  |%
  ::
  ::  ============================================================
  ::  CONSTANTS - These must match exactly between Hoon and C jet
  ::  ============================================================
  ::
  ::  Pi-related constants for argument reduction
  ::
  ++  pio2-1    .1.5707963         ::  first 25 bits of pi/2 (0x3fc90fda)
  ++  pio2-1t   .1.5893255e-8      ::  pi/2 - pio2-1 (0x33a22168)
  ++  pio4      .0.78539816        ::  pi/4 (0x3f490fda)
  ++  invpio2   .0.6366198         ::  2/pi (0x3f22f984)
  ::
  ::  Exp constants
  ::
  ++  ln2hi     .0.6931381         ::  high bits of ln(2) (0x3f317200)
  ++  ln2lo     .9.058001e-6       ::  low bits of ln(2) (0x3717f7d1)
  ++  invln2    .1.442695          ::  1/ln(2) (0x3fb8aa3b)
  ++  exp-overflow   .88.72284     ::  value above which exp overflows
  ++  exp-underflow  .-103.97208   ::  value below which exp underflows to 0
  ::
  ::  ============================================================
  ::  POLYNOMIAL COEFFICIENTS - from musl __sindf.c, __cosdf.c
  ::  ============================================================
  ::
  ::  Sin coefficients: approximates sin(x)/x - 1 on [-pi/4, pi/4]
  ::  sin(x) = x + x^3*S1 + x^5*S2 + x^7*S3 + x^9*S4
  ::
  ++  sin-s1  .-0.16666667        ::  -0x1.555556p-3 = -1/6
  ++  sin-s2  .8.3333291e-3       ::   0x1.11110ap-7 =  1/120
  ++  sin-s3  .-1.9839358e-4      ::  -0x1.a013a2p-13 = -1/5040
  ++  sin-s4  .2.718311e-6        ::   0x1.6d2d70p-19 = 1/362880
  ::
  ::  Cos coefficients: approximates cos(x) - 1 + x^2/2 on [-pi/4, pi/4]
  ::  cos(x) = 1 - x^2/2 + x^4*C0 + x^6*C1 + x^8*C2 + x^10*C3
  ::
  ++  cos-c0  .0.041666668        ::   0x1.555556p-5 =  1/24
  ++  cos-c1  .-1.388889e-3       ::  -0x1.6c16c2p-10 = -1/720
  ++  cos-c2  .2.4801588e-5       ::   0x1.a00eb0p-16 = 1/40320
  ++  cos-c3  .-2.7557314e-7      ::  -0x1.27e4a8p-22 = -1/3628800
  ::
  ::  Exp coefficients: approximates (exp(r)-1-r)/r^2 on [-ln2/2, ln2/2]
  ::  exp(r) = 1 + r + r^2/2 + r^3*P1 + r^4*P2 + r^5*P3
  ::
  ++  exp-p1  .0.16666667         ::  1/6
  ++  exp-p2  .0.041666668        ::  1/24  (actually uses different form)
  ::
  ::  ============================================================
  ::  HELPER FUNCTIONS
  ::  ============================================================
  ::
  ::  +scalbn: multiply by 2^n (ldexp equivalent)
  ::  This is used for exp reconstruction
  ::
  ++  scalbn
    |=  [x=@rs n=@sd]
    ^-  @rs
    ::  Extract sign, exponent, mantissa from x
    =/  bits=@  `@`x
    =/  exp=@   (dis (rsh [0 23] bits) 0xff)
    =/  sign=@  (rsh [0 31] bits)
    ::  Handle special cases
    ?:  =(exp 0xff)  x                    ::  inf or nan
    ?:  =(exp 0)     x                    ::  zero or subnormal (simplified)
    ::  Compute new exponent
    =/  n-val=@sd  n
    =/  new-exp=@sd  (sum:si (sun:si exp) n-val)
    ::  Check for overflow/underflow
    ?:  (gth:si new-exp (sun:si 254))
      ?:(=(sign 0) `@rs`0x7f80.0000 `@rs`0xff80.0000)  ::  +/- inf
    ?:  (lth:si new-exp (sun:si 1))
      .0                                  ::  underflow to zero (simplified)
    ::  Reconstruct float
    =/  new-bits=@
      %+  con  (lsh [0 31] sign)
      %+  con  (lsh [0 23] (abs:si new-exp))
      (dis bits 0x7f.ffff)
    `@rs`new-bits
  ::
  ::  +floor-int: floor to integer, return as @sd
  ::
  ++  floor-int
    |=  x=@rs
    ^-  @sd
    =/  bits=@   `@`x
    =/  sign=@   (rsh [0 31] bits)
    =/  exp=@    (sub (dis (rsh [0 23] bits) 0xff) 127)
    =/  mant=@   (con (dis bits 0x7f.ffff) 0x80.0000)
    ::  If exponent < 0, result is 0 or -1
    ?:  (^lth exp 0)
      ?:(=(sign 0) --0 -1)
    ::  If exponent >= 23, no fractional part
    ?:  (^gte exp 23)
      ?:  =(sign 0)
        (sun:si (lsh [0 (sub exp 23)] mant))
      (new:si %.n (lsh [0 (sub exp 23)] mant))
    ::  Otherwise, mask off fractional bits
    =/  int-mant=@  (rsh [0 (sub 23 exp)] mant)
    =/  result=@sd
      ?:(=(sign 0) (sun:si int-mant) (new:si %.n int-mant))
    ::  Floor: if negative and had fractional part, subtract 1
    =/  had-frac=?  !=(0 (dis mant (dec (lsh [0 (sub 23 exp)] 1))))
    ?:  &(=(sign 1) had-frac)
      (dif:si result --1)
    result
  ::
  ::  ============================================================
  ::  ARGUMENT REDUCTION FOR SIN/COS
  ::  ============================================================
  ::
  ::  +rem-pio2: Reduce x to y in [-pi/4, pi/4], return quadrant
  ::
  ::  Returns [n y] where:
  ::    n = quadrant (0-3)
  ::    y = reduced argument
  ::    x = n*(pi/2) + y
  ::
  ::  This is a simplified version for |x| < 2^7 * pi/2
  ::  A full implementation would handle larger arguments
  ::
  ++  rem-pio2
    |=  x=@rs
    ^-  [n=@ y=@rs]
    ::  Compute n = round(x * 2/pi)
    =/  fn=@rs  (add:^rs (mul:^rs x invpio2) .0.5)
    =/  n=@sd   (floor-int fn)
    =/  fn=@rs  (san:^rs n)
    ::  Compute y = x - n*pi/2 using extended precision
    ::  y = x - fn*pio2_1 - fn*pio2_1t
    =/  y=@rs   (sub:^rs x (mul:^rs fn pio2-1))
    =.  y       (sub:^rs y (mul:^rs fn pio2-1t))
    ::  Return quadrant mod 4 and reduced argument
    [(mod (abs:si n) 4) y]
  ::
  ::  ============================================================
  ::  KERNEL FUNCTIONS - Fixed polynomial evaluation
  ::  ============================================================
  ::
  ::  +sindf: Kernel sin for reduced argument |x| <= pi/4
  ::
  ::  Uses degree-9 polynomial:
  ::    sin(x) = x + x^3*(S1 + x^2*(S2 + x^2*(S3 + x^2*S4)))
  ::
  ::  Written in Horner form for efficiency and determinism:
  ::    sin(x) = x * (1 + z*(S1 + z*(S2 + z*(S3 + z*S4))))
  ::  where z = x^2
  ::
  ++  sindf
    |=  x=@rs
    ^-  @rs
    =/  z=@rs   (mul:^rs x x)                              ::  z = x^2
    =/  w=@rs   (mul:^rs z z)                              ::  w = x^4
    ::  r = S3 + z*S4
    =/  r=@rs   (add:^rs sin-s3 (mul:^rs z sin-s4))
    ::  r = S2 + z*r = S2 + z*S3 + z^2*S4
    =.  r       (add:^rs sin-s2 (mul:^rs z r))
    ::  Compute: x + x^3*S1 + x^5*(S2 + z*(S3 + z*S4))
    ::         = x + x*z*S1 + x*z^2*r
    ::         = x + x*z*(S1 + z*r)
    =/  s=@rs   (mul:^rs z x)                              ::  s = x^3
    (add:^rs x (mul:^rs s (add:^rs sin-s1 (mul:^rs z r))))
  ::
  ::  +cosdf: Kernel cos for reduced argument |x| <= pi/4
  ::
  ::  Uses degree-10 polynomial:
  ::    cos(x) = 1 - x^2/2 + x^4*(C0 + x^2*(C1 + x^2*(C2 + x^2*C3)))
  ::
  ::  Written for numerical stability:
  ::    cos(x) = 1 - hz + (z*z*(C0 + z*(C1 + z*(C2 + z*C3))))
  ::  where z = x^2, hz = z/2
  ::
  ++  cosdf
    |=  x=@rs
    ^-  @rs
    =/  z=@rs   (mul:^rs x x)                              ::  z = x^2
    =/  w=@rs   (mul:^rs z z)                              ::  w = x^4
    ::  r = C2 + z*C3
    =/  r=@rs   (add:^rs cos-c2 (mul:^rs z cos-c3))
    ::  r = C1 + z*r
    =.  r       (add:^rs cos-c1 (mul:^rs z r))
    ::  r = C0 + z*r
    =.  r       (add:^rs cos-c0 (mul:^rs z r))
    ::  cos = 1 - z/2 + z^2*r
    =/  hz=@rs  (mul:^rs .0.5 z)
    (add:^rs (sub:^rs .1 hz) (mul:^rs w r))
  ::
  ::  ============================================================
  ::  PUBLIC API - sin and exp with fixed-polynomial evaluation
  ::  ============================================================
  ::
  ::  +sin: Sine using argument reduction + polynomial
  ::
  ::  Algorithm:
  ::    1. Handle special cases (tiny, inf, nan)
  ::    2. Reduce argument to [-pi/4, pi/4]
  ::    3. Based on quadrant, compute sin or cos with sign adjustment
  ::
  ++  sin
    |=  x=@rs
    ^-  @rs
    =/  bits=@  `@`x
    =/  sign=@  (rsh [0 31] bits)
    =/  ax=@rs  `@rs`(dis bits 0x7fff.ffff)               ::  |x|
    ::
    ::  Special case: tiny argument |x| < 2^-12
    ::  sin(x) ~ x for small x
    ?:  (lth:^rs ax .2.4414063e-4)                        ::  0x39800000
      x
    ::
    ::  Special case: inf or nan
    ?:  (gte:^rs ax `@rs`0x7f80.0000)
      (sub:^rs x x)                                        ::  return nan
    ::
    ::  Argument reduction
    =/  [n=@ y=@rs]  (rem-pio2 x)
    ::
    ::  Based on quadrant n mod 4:
    ::    n=0: sin(y)
    ::    n=1: cos(y)
    ::    n=2: -sin(y)
    ::    n=3: -cos(y)
    ::
    ?-  n
      %0  (sindf y)
      %1  (cosdf y)
      %2  (sub:^rs .0 (sindf y))
      %3  (sub:^rs .0 (cosdf y))
    ==
  ::
  ::  +cos: Cosine using argument reduction + polynomial
  ::
  ++  cos
    |=  x=@rs
    ^-  @rs
    =/  bits=@  `@`x
    =/  ax=@rs  `@rs`(dis bits 0x7fff.ffff)               ::  |x|
    ::
    ::  Special case: tiny argument |x| < 2^-12
    ::  cos(x) ~ 1 for small x
    ?:  (lth:^rs ax .2.4414063e-4)
      .1
    ::
    ::  Special case: inf or nan
    ?:  (gte:^rs ax `@rs`0x7f80.0000)
      (sub:^rs x x)                                        ::  return nan
    ::
    ::  Argument reduction
    =/  [n=@ y=@rs]  (rem-pio2 x)
    ::
    ::  Based on quadrant n mod 4:
    ::    n=0: cos(y)
    ::    n=1: -sin(y)
    ::    n=2: -cos(y)
    ::    n=3: sin(y)
    ::
    ?-  n
      %0  (cosdf y)
      %1  (sub:^rs .0 (sindf y))
      %2  (sub:^rs .0 (cosdf y))
      %3  (sindf y)
    ==
  ::
  ::  +exp: Exponential using argument reduction + polynomial
  ::
  ::  Algorithm (simplified polynomial version, not table-based):
  ::    1. Handle special cases (0, inf, nan, overflow, underflow)
  ::    2. Reduce: x = k*ln2 + r where |r| <= ln2/2
  ::    3. Compute exp(r) using polynomial
  ::    4. Return 2^k * exp(r)
  ::
  ::  For full musl compatibility, a table-based version should be used.
  ::  This simplified version demonstrates the fixed-operation approach.
  ::
  ++  exp
    |=  x=@rs
    ^-  @rs
    =/  bits=@  `@`x
    =/  sign=@  (rsh [0 31] bits)
    =/  ax=@rs  `@rs`(dis bits 0x7fff.ffff)
    ::
    ::  Special case: x = 0
    ?:  =(x .0)  .1
    ::
    ::  Special case: +inf
    ?:  =(bits 0x7f80.0000)  `@rs`0x7f80.0000
    ::
    ::  Special case: -inf
    ?:  =(bits 0xff80.0000)  .0
    ::
    ::  Special case: nan
    ?:  (gth:^rs ax `@rs`0x7f80.0000)
      x                                                    ::  return nan
    ::
    ::  Check overflow: exp(x) for x > 88.72 overflows
    ?:  (gth:^rs x exp-overflow)
      `@rs`0x7f80.0000                                     ::  +inf
    ::
    ::  Check underflow: exp(x) for x < -103.97 underflows
    ?:  (lth:^rs x exp-underflow)
      .0
    ::
    ::  Argument reduction: x = k*ln2 + r
    ::  k = round(x / ln2), r = x - k*ln2
    ::
    ::  Use extended precision for ln2 to minimize error:
    ::  r = (x - k*ln2hi) - k*ln2lo
    ::
    =/  fn=@rs  (add:^rs (mul:^rs x invln2) .0.5)
    =/  k=@sd   (floor-int fn)
    =/  kf=@rs  (san:^rs k)
    ::
    ::  r = x - k*ln2hi - k*ln2lo
    =/  r=@rs   (sub:^rs x (mul:^rs kf ln2hi))
    =.  r       (sub:^rs r (mul:^rs kf ln2lo))
    ::
    ::  Compute exp(r) using polynomial approximation
    ::  exp(r) = 1 + r + r^2/2 + r^3/6 + r^4/24 + r^5/120
    ::
    ::  Horner form: exp(r) = 1 + r*(1 + r*(1/2 + r*(1/6 + r*(1/24 + r/120))))
    ::
    =/  r2=@rs  (mul:^rs r r)                              ::  r^2
    ::
    ::  Compute polynomial: 1/2 + r*(1/6 + r*(1/24 + r/120))
    ::  P = 0.5 + r*(0.16666667 + r*(0.041666668 + r*0.0083333336))
    ::
    =/  p=@rs   .0.0083333336                              ::  1/120
    =.  p       (add:^rs .0.041666668 (mul:^rs r p))       ::  1/24 + r/120
    =.  p       (add:^rs .0.16666667 (mul:^rs r p))        ::  1/6 + ...
    =.  p       (add:^rs .0.5 (mul:^rs r p))               ::  1/2 + ...
    ::
    ::  exp(r) = 1 + r + r^2*P = 1 + r*(1 + r*P)
    =/  expr=@rs  (add:^rs .1 (mul:^rs r (add:^rs .1 (mul:^rs r p))))
    ::
    ::  Result = 2^k * exp(r)
    ::  Use scalbn to multiply by 2^k
    (scalbn expr k)
  ::
  ::  +tan: Tangent as sin/cos
  ::
  ++  tan
    |=  x=@rs
    ^-  @rs
    (div:^rs (sin x) (cos x))
  --
::
::  ============================================================
::  JET IMPLEMENTATION NOTES
::  ============================================================
::
::  For the C jet to produce identical results:
::
::  1. Use SoftFloat for all floating-point operations
::     - Ensures IEEE 754 compliance regardless of hardware
::     - Handles rounding modes correctly
::
::  2. Use the EXACT same polynomial coefficients
::     - Define as static const float/double with hex literals
::     - Example: static const float S1 = -0x1.555556p-3f;
::
::  3. Perform operations in the EXACT same order
::     - Follow the Horner evaluation exactly as written
::     - Don't let the C compiler reorder operations
::     - Use volatile or explicit temporaries if needed
::
::  4. Handle special cases identically
::     - Same checks for inf, nan, overflow, underflow
::     - Same return values for edge cases
::
::  Example C jet for sindf:
::
::  ```c
::  static const float S1 = -0.16666667f;  // -1/6
::  static const float S2 =  0.008333329f; //  1/120
::  static const float S3 = -0.00019839358f;
::  static const float S4 =  0.0000027183114f;
::
::  float __sindf(float x) {
::      float z = x * x;
::      float r = S3 + z * S4;
::      r = S2 + z * r;
::      float s = z * x;
::      return x + s * (S1 + z * r);
::  }
::  ```
::
--
