/+  *math   :: notably, this shadows floating-point doors
  ::
::::  Mathematical library
::
::  UNIVERSAL NUMBERS (UNUMS)
::
::  Type III universal numbers (unums) are divided into three classes:
::  1. Posits, or precise real numbers (well, as real as floats anyway).
::  2. Quires, or compensated auxiliary sums.
::  3. Valids, or bounded intervals.
::
::  We represent these using Hoon auras at three bitwidths:
::  - @rpb, @rph, @rps  :: posits
::  - @rqb, @rqh, @rqs  :: quires
::  - @rvb, @rvh, @rvs  :: valids
::
::  The bitwidths are 8-bit Bytes, 16-bit Halfs, and 32-bit Singles.
::  While we square on quires and valids, the current implementation is
::  only for posits.
::
::  A posit has four fields:
::  - Sign bit (1), 0 for positive, 1 for negative.
::  - Regime bits (1--(ps-1)), unary run-length encoded.
::    - Regime scale is 2^2^es, where es is the max exponent size.
::  - Exponent bits (0--(ps-2)), fixed if available but can be truncated.
::  - Fraction bits (the rest), remaining bits to total bitwidth.
::
::  - Gustafson & Yonemoto (2017), "Beating Floating Point at its Own Game:
::    Posit Arithmetic", Supercomputing Frontiers and Innovations. 4 (2).
::    Publishing Center of South Ural State University, Chelyabinsk, Russia.
::    http://www.johngustafson.net/pdfs/BeatingFloatingPoint.pdf
|%
::  Representation of a Type III Unum Posit
::  (uq, uv)
+$  up
  $%  $:  %p    :: real-valued posit
          b=@u  :: bitwidth (bloq), ∈ 3 (byte), 4 (half), 5 (single)
          s=?   :: sign, 0 (+) or 1 (-)
          r=@s  :: regime
          e=@s  :: exponent
          f=@u  :: fraction
      ==
      [%n ~]    :: Not a Real (NaR), unum NaN
  ==
::  TODO possibly superfluous
::  8-bit posit tuple
+$  pb
  $%  s=?         :: sign bit
      r=@         :: regime bits, bitwidth first (1--7 bits)
      e=@         :: exponent bits, bitwidth first (0 for posit8)
      f=@         :: fraction bits, remaining to total bitwidth (the rest)
  ==
::  16-bit posit tuple
+$  ph
  $%  s=?         :: sign bit
      r=@         :: regime bits, bitwidth first (1--15 bits)
      e=@         :: exponent bits, bitwidth first (1 for posit16)
      f=@         :: fraction bits, remaining to total bitwidth (the rest)
  ==
::  32-bit posit tuple
+$  ps
  $%  s=?         :: sign bit
      r=@         :: regime bits, bitwidth first (1--31 bits)
      e=@         :: exponent bits, bitwidth first (2 for posit32)
      f=@         :: fraction bits, remaining to total bitwidth (the rest)
  ==
++  rpb
  |%
  ::
  ::    +huge:  @rpb
  ::
  ::  Returns the value of the largest representable number.
  ::    Examples
  ::      > `@ub`huge
  ::      0b111.1111
  ::  Source
  ++  huge  `@rpb`0x7f         ::  2**24
  ::    +tiny:  @rpb
  ::
  ::  Returns the value of the smallest representable normal number.
  ::    Examples
  ::      > `@ub`tiny
  ::      0b1
  ::  Source
  ++  tiny  `@rpb`0x1          ::  2**-24
  ::    +nar:  @rpb
  ::
  ::  Returns the value of Not a Real (NaR).
  ::    Examples
  ::      > `@ub`nar
  ::      0b1000.0000
  ::  Source
  ++  nar  `@rpb`0x80          ::  NaR, Not a Real (NaN)
  ::
  ::  Operations
  ::
  ::    +sea:  @rpb -> $up
  ::
  ::  Returns the +$up representation of a posit atom.
  ::    Examples
  ::      :: posit 1.0
  ::      > (sea 0b100.0000)
  ::      [%p s=%.y r=1 e=0 f=0b0]
  ::      :: posit 0.5
  ::      > (sea 0b10.0000)
  ::      [%p s=%.y r=1 e= f=0b0]
  ::      :: posit largest possible 8-bit value, 2**24
  ::      > (sea 0b111.1111)
  ::      [%p s=%.y r=0b111.1111 e=0b0 f=0b0]
  ::      :: posit largest possible negative 8-bit value, -2**24+1
  ::  Source
  ++  sea
    :: =seaa !:
    |=  =@rpb
    ^-  up
    |^
    ::  Sign bit at MSB.
    =/  s=?  ;;(? (rsh [0 7] rpb))
    ::  Regime bits, unary run-length encoding.
    =+  [r0 k]=(get-regime rpb)
    =/  r=@s  ?:(r0 (dif:si --0 (sun:si k)) (sun:si (dec k)))
    ::  Exponent bits, zero in posit8.
    =/  e=@  0
    ::  Fraction bits, remaining bits to total bitwidth.
    =/  f=@  (dis (dec (bex +(+(k)))) rpb)
    [%p 8 s r e f]
    ::  Retrieve unary run-length encoded regime bits.
    ::  k in Gustafson's notation.
    ++  get-regime
      |=  p=@
      ^-  [r0=? k=@]
      ::  get RLE bit polarity
      =/  r0=?  ;;(? (rsh [0 6] (dis 0x40 p)))
      ::  get RLE bit count
      =|  b=_5
      =|  k=@
      :-  r0
      |-
      =/  r  ;;(? (rsh [0 b] (dis (bex b) p)))
      ?:  !=(r0 r)  +(k)
      ?:  =(0 b)  +(k)
      $(k +(k), b (dec b))
    --
  ::
  ::    +to-rs:  @pb -> @rs
  ::
  ::  Returns the closest @rs equivalent to a posit.
  ::  $((1-3s)+f)\times 2^{(1-2s)\times(4r+e+s)}$
  ::
  ::  Source
  :: ++  to-rs
  ::   |=  =@rpb
  ::   ^-  @rs
  ::   =/  =pb  (sea rpb)
  ::   =/  s=@  `@`s.pb
  ::   =/  r=@  (get-regime r.pb)
  ::   =/  e=@  (sun:^rs (to-sd:twoc e.pb))
  ::   =/  f=@  (sun:^rs (to-sd:twoc f.pb))
  ::   %+  mul:rs
  ::     %+  add:rs
  ::       (sub:rs .1 (mul:rs .3 s))
  ::     f
  ::   %-  pow-n:rs
  ::   %+  mul:rs
  ::     :(sub:rs .1 s s)
  ::   :(add:rs (mul:rs .4 r) e s)
  :: ::
  :: ++  to-rpb  !!
  :: ::
  :: ++  add
  ::   |=  [a=@rpb b=@rpb]
  ::   ^-  @rpb
  :: - `++add`, $+$ addition
  :: - `++sub`, $-$ subtraction
  :: - `++mul`, $\times$ multiplication
  :: - `++div`, $/$ division
  :: - `++mod`, modulo (remainder after division)
  :: - `++fma`, $\text{fma}$ fused multiply-add
  :: - `++sgn`, $\text{sgn}$ signum (also `++sig`)
  :: - `++neg`, $-$ unary negation
  :: - `++factorial`, $!$ factorial
  :: - `++abs`, $\text{abs}$
  :: - `++exp`, $\exp$
  :: - `++sin`, $\sin$
  :: - `++cos`, $\cos$
  :: - `++tan`, $\tan$
  :: - `++pow-n`, $\text{pow}$ to integer power
  :: - `++log`, $\log$ (natural logarithm)
  :: - `++log-10`, $\log_{10}$ (log base-10)
  :: - `++log-2`, $\log_{2}$ (log base-2)
  :: - `++pow`, $\text{pow}$
  :: - `++sqrt`, $\sqrt{}$ (also `++sqt`)
  :: - `++cbrt`, $\sqrt[3]{}$ (also `++cbt`)
  :: - `++arg` (alias for `++abs` in real values)

  :: Logical functions:

  :: - `++lth`, $<$
  :: - `++lte` $\leq$ (also `++leq`)
  :: - `++gth`, $>$
  :: - `++gte`, $\geq$ (also `++geq`)
  :: - `++equ`, $=$
  :: - `++neq`, $\neq$
  :: - `is-close`
  :: - `all-close`
  :: - `is-int`

  :: Constants (to machine accuracy):

  :: - `++tau`, $\tau = 2 \pi$
  :: - `++pi`, $\pi$
  :: - `++e`, $e$ (base of natural logarithm)
  :: - `++phi`, $\phi$, Euler's constant
  :: - `++sqt2`, $\sqrt{2}$
  :: - `++invsqt2`, $\frac{1}{\sqrt{2}}$
  :: - `++log2`, $\log 2$
  :: - `++invlog2`, $\frac{1}{\log 2}$
  :: - `++log10`, $\log 10$
  :: - `++huge`, largest valid number in `bloq` width
  :: - `++tiny`, smallest valid number in `bloq` size

  --
::  Type III Unum Posit, 16-bit width ("half")
::  Type III Unum Posit, 32-bit width ("single")
--
