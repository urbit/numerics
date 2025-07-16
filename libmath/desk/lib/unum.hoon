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
::
::  While posits can be written in a general purpose form, we are interested
::  in standard posit8, posit16, and posit32 representations.  For these, the
::  following conventions apply:
::  - posit8:  8 bits total (posit<8,0>)
::    - Sign bit s (1)
::    - Regime bits k (1--7), unary run-length encoded (one different to end)
::    - Exponent bits e (0), fixed
::    - Fraction bits f (the rest), remaining to total bitwidth
::  - posit16: 16 bits total (posit<16,1>)
::    - Sign bit s (1)
::    - Regime bits k (1--15), unary run-length encoded (one different to end)
::    - Exponent bits e (1), fixed (but may be occluded by a full regime)
::    - Fraction bits f (the rest), remaining to total bitwidth
::  - posit32: 32 bits total (posit<32,2>)
::    - Sign bit s (1)
::    - Regime bits k (1--31), unary run-length encoded (one different to end)
::    - Exponent bits e (2), fixed (but may be occluded by a full regime)
::    - Fraction bits f (the rest), remaining to total bitwidth
::
::  In general, a posit<bw,es> has a bitwidth of bw and a useed = 2^2^es (or
::  regime multiplier) of $2^{2^\text{es}}$.  These two equations are
::  equivalent:
::
::  $$(-1)^s \times \text{useed}^r \times 2^e \times (1+f/{2^{fs}})$$
::  $$(1-3s+f) \times 2^{(-1)^s (r \times 2^\text{es} + e + s)}
::
::  (Note that Gustafson equivocates between r and k in notation.)
::
::  - Cook (2018), "Eight-bit floating point",
::    https://www.johndcook.com/blog/2018/04/15/eight-bit-floating-point/
::  - Cook (2018), "All posit<8,0> representable numbers",
::    https://www.johndcook.com/eightbit.html
::  - Gustafson & Yonemoto (2017), "Beating Floating Point at its Own Game:
::    Posit Arithmetic", Supercomputing Frontiers and Innovations. 4 (2).
::    Publishing Center of South Ural State University, Chelyabinsk, Russia.
::    http://www.johngustafson.net/pdfs/BeatingFloatingPoint.pdf
::  - Gustafson (2019), "Posit Arithmetic Lecture 1", CoNGa 2019,
::    https://posithub.org/conga/2019/docs/13/1430-John-Introductory.pdf
|%
::  Representation of a Type III Unum Posit
::
::  A posit has four fields:
::  - Sign bit (1), 0 for positive, 1 for negative.
::  - Regime bits (1--(ps-1)), unary run-length encoded.
::    - Regime scale is 2^2^es, where es is the max exponent size.
::  - Exponent bits (0--(ps-2)), fixed if available but can be truncated.
::  - Fraction bits (the rest), remaining bits to total bitwidth.
::
+$  up
  $%  $:  %p    :: real-valued posit
          b=@u  :: bitwidth (bloq), ∈ 3 (byte), 4 (half), 5 (single)
          s=?   :: sign, 0 (+) or 1 (-)
          x=@s  :: total exponent (to be broken into regime and exponent)
          f=@u  :: fraction
      ==
      [%n b=@u ~]   :: Not a Real (NaR), unum NaN
      [%z b=@u ~]   :: Zero, unum 0
  ==
::
::  Representation of a Type III Unum Quire
::
::  > The other major component of the definition of posit arithmetic is the
::  quire, the fixed-size set of bits used for scratch operations that are
::  mathematically exact in the g-layer. A quire can be thought of as a
::  dedicated register that permits dot products, sums, and other operations to
::  be performed with rounding error deferred to the very end of the
::  calculation. All computers use a hidden scratch area for temporary
::  calculations. In posit arithmetic, the quire data type is accessible to the
::  programmer, which is what makes possible for posits to follow the rules of
::  algebra much more closely than floats do.  (Posits4.nb)
+$  uq
  $:  p=@ux
  ==
::
++  rpb
  |%
  ::  mathematics constants to posit8 precision
  ::  TODO
  ++  es  0
  ::    +useed:  @
  ::
  ::  Returns the value of useed, $2^{2^{es}}$, as used in this posit.
  ::    Examples
  ::      > `@ub`useed
  ::      0b100.0000
  ::  Source
  ++  useed  2
  ::    +zero:  @rpb
  ::
  ::  Returns the value of zero.
  ::    Examples
  ::      > `@ub`zero
  ::      0b0
  ::  Source
  ++  zero  `@rpb`0x0          ::  0
  ::    +one:  @rpb
  ::
  ::  Returns the value of one.
  ::    Examples
  ::      > `@ub`one
  ::      0b0
  ::  Source
  ++  one  `@rpb`0x40          ::  1
  ::    +neg-one:  @rpb
  ::
  ::  Returns the value of negative one.
  ::    Examples
  ::      > `@ub`neg-one
  ::      0b1100.0000
  ::  Source
  ++  neg-one  `@rpb`0xc0      ::  -1
  ::    +nar:  @rpb
  ::
  ::  Returns the value of Not a Real (NaR).
  ::    Examples
  ::      > `@ub`nar
  ::      0b1000.0000
  ::  Source
  ++  nar  `@rpb`0x80          ::  NaR, Not a Real (NaN)
  ::    +pi:  @rpb
  ::
  ::  Returns the value of pi.
  ::    Examples
  ::      > `@ub`pi
  ::      0b110.1001
  ::  Source
  ++  pi  `@rpb`0x69           ::  π
  ::    +tau:  @rpb
  ::
  ::  Returns the value of tau, $2\pi$.
  ::    Examples
  ::      > `@ub`tau
  ::      0b111.0101
  ::  Source
  ++  tau  `@rpb`0x75          ::  τ
  ::    +e:  @rpb
  ::
  ::  Returns the value of e.
  ::    Examples
  ::      > `@ub`e
  ::      0b110.0110
  ::  Source
  ++  e  `@rpb`0x66            ::  e
  ::    +phi:  @rpb
  ::
  ::  Returns the value of phi, $\frac{1+\sqrt{5}}{2}$.
  ::    Examples
  ::      > `@ub`phi
  ::      0b101.0100
  ::  Source
  ++  phi  `@rpb`0x54          ::  φ
  ::
  ::    +sqt2:  @rpb
  ::
  ::  Returns the value of sqrt(2).
  ::    Examples
  ::      > `@ub`sqt2
  ::      0b100.1101
  ::  Source
  ++  sqt2  `@rpb`0x4d         ::  √2
  ::
  ::  TODO other constants
  ::
  ::    +huge:  @rpb
  ::
  ::  Returns the value of the largest representable number.
  ::    Examples
  ::      > `@ub`huge
  ::      0b111.1111
  ::  Source
  ++  huge  `@rpb`0x7f         ::  2**6
  ::
  ::    +neghuge:  @rpb
  ::
  ::  Returns the value of the largest representable number.
  ::    Examples
  ::      > `@ub`neg-huge
  ::      0b1000.0001
  ::  Source
  ++  neg-huge  `@rpb`0x81     ::  -2**6
  ::    +tiny:  @rpb
  ::
  ::  Returns the value of the smallest representable normal number.
  ::    Examples
  ::      > `@ub`tiny
  ::      0b1
  ::  Source
  ++  tiny  `@rpb`0x1          ::  2**-6
  ::    +neg-tiny:  @rpb
  ::
  ::  Returns the value of the smallest representable negative normal number.
  ::    Examples
  ::      > `@ub`neg-tiny
  ::      0b1111.1111
  ::  Source
  ++  neg-tiny  `@rpb`0xff     ::  -2**-6
  ::    +invsqt2:  @rs
  ::
  ::  Returns the value 1/sqrt(2) (OEIS A010503).
  ::    Examples
  ::      > invsqt2
  ::      0b10.1101
  ::  Source
  ++  invsqt2  `@rpb`0x2d    ::  1/√2
  ::    +log2:  @rs
  ::
  ::  Returns the value log(2) (OEIS A002162).
  ::    Examples
  ::      > log2
  ::      0b10.1100
  ::  Source
  ++  log2  `@rpb`0x2c    ::  log(2)
  ::    +invlog2:  @rs
  ::
  ::  Returns the value 1/log(2).
  ::    Examples
  ::      > invlog2
  ::      0b100.1110
  ::  Source
  ++  invlog2  `@rpb`0x4e  ::  1/log(2)
  ::    +log10:  @rs
  ::
  ::  Returns the value log(10) (OEIS A002392).
  ::    Examples
  ::      > log10
  ::      0b110.0010
  ::  Source
  ++  log10  `@rpb`0x62    ::  log(10)
  ::
  ::  Operations
  ::
  ::    +from:  @rpb -> $up
  ::
  ::  Returns the +$up representation of a posit atom.
  ::    Examples
  ::      :: posit 1.0
  ::      > (from 0b100.0000)
  ::      [%p b=3 s=%.y r=--0 e=--0 f=0]
  ::      :: posit 0.5
  ::      > (from 0b10.0000)
  ::      [%p b=3 s=%.y r=-1 e=--0 f=0]
  ::      :: posit largest possible 8-bit value, 2**24
  ::      > (from 0b111.1111)
  ::      [%p b=3 s=%.y r=--5 e=--0 f=127]
  ::      :: posit largest possible negative 8-bit value, -2**24+1
  ::  Source
  ++  from
    |=  =@rpb
    ^-  up
    |^
    ?:  =(0x0 rpb)   [%z 3 ~]
    ?:  =(0x80 rpb)  [%n 3 ~]
    ::  Sign bit at MSB.
    =/  s=?  ;;(? (rsh [0 7] rpb))
    ::  Regime bits, unary run-length encoding.
    =+  [r0 k]=(get-regime rpb)
    =/  r=@s  ?:(r0 (dif:si --0 (sun:si k)) (sun:si (dec k)))
    ::  Exponent bits, zero in posit8.
    =/  e=@  0
    =/  x=@s  (re-to-exp r e)
    ::  Fraction bits, remaining bits to total bitwidth.
    =/  f=@  (dis (dec (bex (sub 6 k))) rpb)
    [%p 3 s x f]
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
      ?:  =(0 b)    +(+(k))  :: no empty unary terminator
      $(k +(k), b (dec b))
    ::  x = r*2^es + e
    ++  re-to-exp
      |=  [r=@s e=@]
      ^-  @s
      (sum:si (pro:si r (new:si & (bex es))) e)
    --
  ::
  ::    +into:  $up -> @rpb
  ::
  ::  Returns the closest @rpb equivalent to a posit tuple.
  ::  $((1-3s)+f)\times 2^{(1-2s)\times(2^0r+e+s)}$
  ::
  ::  Source
  ++  into
    |=  =up
    ^-  @rpb
    |^
    =|  rpb=@rpbD
    ?:  ?=(%z -.up)  `@rpb`0x0
    ?:  ?=(%n -.up)  `@rpb`0x80
    ?>  ?=(%p -.up)
    ::  s sign bit
    =.  rpb  (lsh [0 7] s.up)
    ::
    =+  [r e]=(exp-to-er x.up)
    ::  r regime bits
    =+  [sg ab]=(old:si r)
    =/  r0  ?:(sg %| %&)
    =/  k   ?:(sg +(ab) ab)
    ~&  rke+[r k e]
    =.  rpb
      %+  con  rpb
      =/  rs  (fil 0 k r0)
      (lsh [0 (sub 7 k)] rs)
    ?:  (gte k 7)  rpb  :: regime bits are full, no fraction
    ::  regime terminator
    =.  rpb
      %+  con  rpb
      (lsh [0 (sub 7 +(k))] !r0)
    ::  no exponent in posit8
    ?>  =(0 e)
    ::  f fraction bits
    ?:  (gte +(k) 7)  rpb  :: regime bits are full, no fraction
    %+  con  rpb
    =/  sp  (sub 7 +(k))
    ?:  (gth sp (met 0 f.up))
      f.up
    =/  sh  (sub sp (met 0 f.up))
    (rsh [0 sh] f.up)
    ::  Convert base exponent into regime and exponent.
    ::  2^(r*2^es + e) = 2^x
    ::  x = r*2^es+e
    ++  exp-to-er
      |=  x=@s
      ^-  [r=@s e=@]
      ?:  =(0 es)  [x 0]
      =/  lg2useed  ^~((new:si %& (rsh [0 1] useed)))
      :-  (fra:si x lg2useed)
      =+  [sg ab]=(old:si (rem:si x lg2useed))
      ?>  sg
      ab
    --
  :: ++  to-rs
  ::   |=  =up
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
  ::  See the 2022 Posit Standard, section 5, for a full list.
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
::  
++  rqb
  |%
  ++  qsize  !!
  ++  qextra  !!
  --
::  Type III Unum Posit, 16-bit width ("half")
::  Type III Unum Posit, 32-bit width ("single")
--
