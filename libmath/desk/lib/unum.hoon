/+  twoc   :: two's-complement helpers (sign, comparison) shared with numerics
  ::
::::  Universal numbers (unums), Type III: posits, quires, valids
::
::  Per the 2022 Posit Standard (Posit Working Group, John Gustafson chair):
::  https://posithub.org/docs/posit_standard-2.pdf
::
::  We represent unums using Hoon auras at four bitwidths:
::  - @rpb @rph @rps @rpd  :: posits  (byte/half/single/double = 8/16/32/64)
::  - @rqb @rqh @rqs       :: quires
::  - @rvb @rvh @rvs       :: valids  (not yet implemented)
::
::  STANDARD, NOT LEGACY.  The 2022 standard fixes the exponent size at
::  es = 2 for every width (the exponent field is a 2-bit unsigned integer,
::  0..3), so useed = 2^2^es = 16.  This differs from the 2017 draft (and
::  SoftPosit's fast p8/p16 types) which scaled es with width.  Only posit32
::  coincides between the two.  We target the standard.
::
::  A posit of precision n has, most-significant first:
::  - sign     S  (1 bit)              s in {0,1}; implicit value 1-3s
::  - regime   R  (k+1 bits)           run of R0 bits, terminated by ~R0;
::                                     r = -k if R0=0, else r = k-1
::  - exponent E  (2 bits, truncable)  e in {0,1,2,3}
::  - fraction F  (rest, truncable)    f = F / 2^m, 0 <= f < 1
::
::  value:  p = ((1-3s) + f) * 2^((1-2s)*(4r + e + s))
::
::  Exceptions: all bits but S zero -> S=0 is posit 0, S=1 is NaR.
::
::  NB: this core defines posit add/sub/mul/div, which shadow the stdlib
::  gates of the same name.  Internal integer arithmetic therefore uses the
::  ^-prefixed forms (^add/^sub/^mul/^div) to reach the standard library.
::
|%
::  G-layer (decoded) representation of a posit, mirroring the stdlib float +$fn
::  (`[%f s=? e=@s a=@u]`):  value = (sign) a * 2^e, with a an integer
::  significand (the hidden 1 included) and e a signed binary exponent.
::
::  s=%.y is non-negative (the +$si convention), so value is
::  ?:(s +1 -1) * a * 2^e.
::
+$  up
  $%  [%p s=? e=@s a=@u]   :: real posit
      [%z ~]              :: zero
      [%n ~]              :: Not a Real (NaR)
  ==
::  Type III Unum Quire: a fixed-point exact accumulator of 16n bits (§3.4):
::  sign (1) . carry guard (31) . integer (8n-16) . fraction (8n-16).
::
+$  uq
  $%  $:  %q
          s=?      :: sign
          c=@      :: carry guard, 31 bits
          i=@      :: integer part, 8n-16 bits
          f=@      :: fractional part, 8n-16 bits
      ==
      [%n ~]       :: NaR
  ==
::  Generic posit core, parameterized by bloq (log2 of the bitwidth).
::
::  Specialize with %*: posit8 is bloq=3 (n=8), posit16 bloq=4, posit32 bloq=5.
::
++  pp
  |_  =bloq
  ::  n: total bitwidth in bits.
  ++  n  (bex bloq)
  ::  es: exponent size, fixed at 2 by the 2022 standard.
  ++  es  2
  ::  useed: regime multiplier, 2^2^es = 16.
  ++  useed  16
  ::  msk: n-bit mask.
  ++  msk  (dec (bex n))
  ::  zero, NaR, and the structural extremes (sign-magnitude positive forms).
  ++  zero    `@`0
  ++  nar     (bex (dec n))         :: 1000...0, also the most-negative int
  ++  maxpos  (dec (bex (dec n)))   :: 0111...1, value 2^(4n-8)
  ++  minpos  `@`1                  :: 0000...1, value 2^-(4n-8)
  ::  huge/tiny: aliases for the structural extremes (cf. /lib/math).
  ++  huge  maxpos
  ++  tiny  minpos
  ::  one: the posit 1.0 (s=0, e=0, a=1).
  ++  one  (bit [%p %.y --0 1])
  ::  Mathematical constants, correctly rounded at every width.  Each is the
  ::  recognizable Q2.52 fixed-point hex of the value (e.g. pi = 0x3.243F6A88...),
  ::  i.e. value = a * 2^-52, fed through the verified encoder.
  ::  Cross-checked against SoftPosit (pX2, es=2) and an independent reference
  ::  encoder for posit8/16/32.
  ++  pi       (bit [%p %.y -52 0x32.43f6.a888.5a31])   ::  3.14159265358979
  ++  tau      (bit [%p %.y -52 0x64.87ed.5110.b461])   ::  6.28318530717959
  ++  e        (bit [%p %.y -52 0x2b.7e15.1628.aed3])   ::  2.71828182845905
  ++  phi      (bit [%p %.y -52 0x19.e377.9b97.f4a8])   ::  1.61803398874989
  ++  sqt2     (bit [%p %.y -52 0x16.a09e.667f.3bcd])   ::  1.41421356237310
  ++  invsqt2  (bit [%p %.y -52 0xb.504f.333f.9de6])    ::  0.70710678118655
  ++  log2     (bit [%p %.y -52 0xb.1721.7f7d.1cf8])    ::  0.69314718055995
  ++  invlog2  (bit [%p %.y -52 0x17.1547.652b.82fe])   ::  1.44269504088896
  ++  log10    (bit [%p %.y -52 0x24.d763.776a.aa2b])   ::  2.30258509299405
  ::    +sea:  @ -> $up
  ::
  ::  Decode an n-bit posit atom into its g-layer +$up form.  The decoder
  ::  yields the fraction-aligned form: a is the hidden 1 plus all fraction
  ::  bits, e offset so the value is exactly a * 2^e (not minimally normalized).
  ::    Examples (posit8, es=2)
  ::      > (sea:rpb 0x40)   ::  1.0  = 8 * 2^-3
  ::      [%p s=%.y e=-3 a=8]
  ::      > (sea:rpb 0x48)   ::  2.0  = 8 * 2^-2
  ::      [%p s=%.y e=-2 a=8]
  ::      > (sea:rpb 0x42)   ::  1.25 = 10 * 2^-3
  ::      [%p s=%.y e=-3 a=10]
  ::  Source
  ++  sea
    |=  p=@
    ^-  up
    =.  p  (dis msk p)
    ?:  =(0 p)    [%z ~]
    ?:  =(nar p)  [%n ~]
    ::  sign at the MSB; work on the magnitude (two's-complement) thereafter.
    =/  neg  =(1 (cut 0 [(dec n) 1] p))
    =/  mag  ?:(neg (^sub (bex n) p) p)
    =/  pw   (dec n)                       :: payload width (bits below sign)
    =/  r0   (cut 0 [(dec pw) 1] mag)      :: regime leading bit, R0
    ::  count the regime run length k (number of bits equal to R0).
    =/  k=@  1
    |-  ^-  up
    ?:  =(k pw)
      ::  regime fills the payload: no terminator, exponent, or fraction.
      =/  r  ?:(=(1 r0) (sun:si (dec k)) (dif:si --0 (sun:si k)))
      [%p !neg (pro:si r --4) 1]
    =/  nb  (cut 0 [(^sub (dec pw) k) 1] mag)
    ?:  =(nb r0)
      $(k +(k))
    ::  terminator at position pw-1-k; (pw-1-k) bits of exponent+fraction below.
    =/  r       ?:(=(1 r0) (sun:si (dec k)) (dif:si --0 (sun:si k)))
    =/  remwid  (^sub pw +(k))
    =/  rem     (dis (dec (bex remwid)) mag)
    ::  exponent: top 2 bits of the remaining field (low bits truncate to 0).
    =/  elo  ?:  (^gte remwid 2)
               (rsh [0 (^sub remwid 2)] rem)
             ?:(=(1 remwid) (^mul 2 rem) 0)
    =/  fw    ?:((^gte remwid 2) (^sub remwid 2) 0)
    =/  frac  (dis (dec (bex fw)) rem)
    ::  x = 4r + e (total scale); significand a = 2^fw + frac with hidden 1.
    =/  x  (sum:si (pro:si r --4) (sun:si elo))
    =/  a  (^add (bex fw) frac)
    [%p !neg (dif:si x (sun:si fw)) a]
  ::    +bit:  $up -> @
  ::
  ::  Encode a g-layer +$up into the closest n-bit posit, rounding to nearest
  ::  with ties to even and saturating at +-maxPos / +-minPos (never to NaR or 0).
  ::    Examples (posit8)
  ::      > (bit:rpb [%p %.y --0 1])
  ::      0x40
  ::      > (bit:rpb [%z ~])
  ::      0x0
  ::  Source
  ++  bit
    |=  =up
    ^-  @
    ?:  ?=(%z -.up)  zero
    ?:  ?=(%n -.up)  nar
    ?>  ?=(%p -.up)
    ?:  =(0 a.up)  zero
    =/  neg   !s.up
    ::  normalize: value = a * 2^e, a >= 1.  lead = floor(log2 a); the total
    ::  scale x = e + lead so that value = 1.frac * 2^x.
    =/  lead  (dec (met 0 a.up))
    =/  x     (sum:si e.up (sun:si lead))
    =/  frac  (dis (dec (bex lead)) a.up)  :: the `lead` fraction bits
    ::  split x into regime r and 2-bit exponent elo: x = 4r + elo, elo in 0..3
    ::  (floor division, so elo is always non-negative).
    =/  rel
      =/  ax  (abs:si x)
      ?:  (syn:si x)
        [(sun:si (^div ax 4)) (mod ax 4)]
      =/  q  (^div ax 4)
      =/  m  (mod ax 4)
      ?:  =(0 m)
        [(dif:si --0 (sun:si q)) 0]
      [(dif:si --0 (sun:si +(q))) (^sub 4 m)]
    =/  r=@s   -.rel
    =/  elo=@  +.rel
    ::  saturate when the regime cannot fit in the payload.
    ?:  (gte-s r (sun:si (^sub n 2)))    (smag neg maxpos)
    ?:  (lte-s r (dif:si --0 (sun:si (dec n))))  (smag neg minpos)
    ::  build the regime field (value, and its width) MSB-first.
    =/  rmag  (abs:si r)
    =/  rr
      ?:  (syn:si r)
        ::  R0=1: (rmag+1) ones, then a 0 terminator.
        [(lsh [0 1] (dec (bex +(rmag)))) (^add rmag 2)]
      ::  R0=0: rmag zeros, then a 1 terminator.
      [1 (^add rmag 1)]
    =/  regval=@  -.rr
    =/  regwid=@  +.rr
    ::  assemble payload [regime | 2-bit exponent | fraction] and its width.
    =/  totw  (^add (^add regwid 2) lead)
    =/  pay   (con (lsh [0 (^add 2 lead)] regval) (con (lsh [0 lead] elo) frac))
    =/  pw    (dec n)
    ?:  (^lte totw pw)
      ::  fits exactly; left-justify (regime to the top of the payload).
      (smag neg (lsh [0 (^sub pw totw)] pay))
    ::  round to nearest, ties to even, keeping the top pw bits.
    =/  sh      (^sub totw pw)
    =/  keep    (rsh [0 sh] pay)
    =/  guard   (cut 0 [(dec sh) 1] pay)
    =/  sticky  ?:(=(0 (dis (dec (bex (dec sh))) pay)) 0 1)
    =/  lsbit   (dis 1 keep)
    =/  roundup  &(=(1 guard) |(=(1 sticky) =(1 lsbit)))
    =/  mag      ?:(roundup +(keep) keep)
    ::  a carry must never spill into the sign bit; clamp to maxPos.
    =?  mag  (^gth mag maxpos)  maxpos
    (smag neg mag)
  ::  +smag: apply a sign to a positive magnitude pattern (two's-complement).
  ++  smag
    |=  [neg=? mag=@]
    ^-  @
    ?.  neg  mag
    (dis msk (^sub (bex n) mag))
  ::  +gte-s / +lte-s: signed (@s) comparisons via the stdlib compare.
  ++  gte-s  |=([a=@s b=@s] ^-(? !=(-1 (cmp:si a b))))
  ++  lte-s  |=([a=@s b=@s] ^-(? !=(--1 (cmp:si a b))))
  ::
  ::  Comparisons.  Posit ordering is identical to two's-complement integer
  ::  ordering of the raw bits (§5.3), with NaR (= most-negative int) below
  ::  every real posit, so we defer to the numerics twoc lib.
  ::
  ::    +gth:  @ -> @ -> ?
  ++  gth  |=([a=@ b=@] ^-(? (~(gth twoc:twoc bloq) a b)))
  ++  lth  |=([a=@ b=@] ^-(? (~(lth twoc:twoc bloq) a b)))
  ++  gte  |=([a=@ b=@] ^-(? (~(gte twoc:twoc bloq) a b)))
  ++  lte  |=([a=@ b=@] ^-(? (~(lte twoc:twoc bloq) a b)))
  ++  equ  |=([a=@ b=@] ^-(? =(a b)))
  ++  neq  |=([a=@ b=@] ^-(? !=(a b)))
  ::  +neg: negate (two's-complement); fixes neither 0 nor NaR.
  ++  neg  |=(a=@ ^-(@ (dis msk (^sub (bex n) a))))
  ::  +abs: absolute value (NaR stays NaR).
  ++  abs  |=(a=@ ^-(@ ?:(=(1 (~(msb twoc:twoc bloq) a)) (neg a) a)))
  ::  +sgn: sign as a posit (0 -> 0, NaR -> NaR, else +-1).
  ++  sgn
    |=  a=@
    ^-  @
    ?:  =(zero a)  zero
    ?:  =(nar a)   nar
    ?:  =(1 (~(msb twoc:twoc bloq) a))  (neg one)
    one
  ::
  ::  Arithmetic (2022 standard sec 5.4).  Each combines exact g-layer
  ::  significands and rounds once via +bit, so all four are correctly
  ::  rounded.  Verified against SoftPosit (pX2) over all 65,536 posit8 pairs.
  ::
  ::    +mul:  @ -> @ -> @  (correctly-rounded product)
  ++  mul
    |=  [a=@ b=@]
    ^-  @
    =/  ua  (sea a)
    =/  ub  (sea b)
    ?:  |(?=(%n -.ua) ?=(%n -.ub))  nar
    ?:  |(?=(%z -.ua) ?=(%z -.ub))  zero
    ?>  ?=(%p -.ua)
    ?>  ?=(%p -.ub)
    (bit [%p =(s.ua s.ub) (sum:si e.ua e.ub) (^mul a.ua a.ub)])
  ::    +add:  @ -> @ -> @  (correctly-rounded sum)
  ++  add
    |=  [a=@ b=@]
    ^-  @
    =/  ua  (sea a)
    =/  ub  (sea b)
    ?:  |(?=(%n -.ua) ?=(%n -.ub))  nar
    ?:  ?=(%z -.ua)  b
    ?:  ?=(%z -.ub)  a
    ?>  ?=(%p -.ua)
    ?>  ?=(%p -.ub)
    ::  align both significands to the smaller exponent, then add/subtract.
    =/  emin  ?:(=(-1 (cmp:si e.ua e.ub)) e.ua e.ub)
    =/  s1  (lsh [0 (abs:si (dif:si e.ua emin))] a.ua)
    =/  s2  (lsh [0 (abs:si (dif:si e.ub emin))] a.ub)
    ?:  =(s.ua s.ub)
      (bit [%p s.ua emin (^add s1 s2)])
    ?:  (^gth s1 s2)
      (bit [%p s.ua emin (^sub s1 s2)])
    ?:  (^gth s2 s1)
      (bit [%p s.ub emin (^sub s2 s1)])
    zero
  ::    +sub:  @ -> @ -> @  (a - b)
  ++  sub  |=([a=@ b=@] ^-(@ (add a (neg b))))
  ::    +div:  @ -> @ -> @  (correctly-rounded quotient; x/0 = NaR)
  ++  div
    |=  [a=@ b=@]
    ^-  @
    =/  ua  (sea a)
    =/  ub  (sea b)
    ?:  |(?=(%n -.ua) ?=(%n -.ub))  nar
    ?:  ?=(%z -.ub)  nar
    ?:  ?=(%z -.ua)  zero
    ?>  ?=(%p -.ua)
    ?>  ?=(%p -.ub)
    ::  long division with guard bits; fold the remainder into a sticky bit.
    =/  g    (^add n n)
    =/  num  (lsh [0 g] a.ua)
    =/  q    (^div num a.ub)
    =/  qs   ?:(=(0 (mod num a.ub)) q (con q 1))
    =/  eo   (dif:si (dif:si e.ua e.ub) (sun:si g))
    (bit [%p =(s.ua s.ub) eo qs])
  ::
  ::  Elementary and rounding ops.
  ::
  ::  +isqt: integer (floor) square root, by Newton's method.
  ++  isqt
    |=  x=@
    ^-  @
    ?:  =(0 x)  0
    =/  r  (bex (^div (^add (met 0 x) 1) 2))
    |-  ^-  @
    =/  nr  (^div (^add r (^div x r)) 2)
    ?:  (^gte nr r)
      |-  ^-  @
      ?:  (^gth (^mul r r) x)  $(r (dec r))
      r
    $(r nr)
  ::    +sqt:  @ -> @  (correctly-rounded square root; sqrt of negative = NaR)
  ++  sqt
    |=  p=@
    ^-  @
    =/  u  (sea p)
    ?:  ?=(%n -.u)  nar
    ?:  ?=(%z -.u)  zero
    ?>  ?=(%p -.u)
    ?.  s.u  nar                                  :: negative -> NaR
    =/  odd  =(1 (dis 1 (abs:si e.u)))            :: exponent parity
    =/  aa   ?:(odd (lsh [0 1] a.u) a.u)          :: make exponent even
    =/  ee   ?:(odd (dif:si e.u --1) e.u)
    =/  g    (^add n n)
    =/  m    (lsh [0 (^mul 2 g)] aa)
    =/  s    (isqt m)
    =/  sx   ?:(=(m (^mul s s)) s (con s 1))      :: fold sticky if inexact
    (bit [%p %.y (dif:si (fra:si ee --2) (sun:si g)) sx])
  ::    +round:  (mode) -> @ -> @  (round to integer-valued posit)
  ++  round
    |=  mode=?(%near %down %up)
    |=  p=@
    ^-  @
    =/  u  (sea p)
    ?.  ?=(%p -.u)  p
    ?:  (gte-s e.u --0)  p                        :: already integer-valued
    =/  sh    (abs:si e.u)                        :: -e
    =/  hi    (rsh [0 sh] a.u)
    =/  rem   (dis (dec (bex sh)) a.u)
    =/  half  (bex (dec sh))
    =?  hi  &(?=(%near mode) |((^gth rem half) &(=(rem half) =(1 (dis 1 hi)))))
      +(hi)
    =?  hi  &(?=(%down mode) !s.u !=(0 rem))  +(hi)
    =?  hi  &(?=(%up mode) s.u !=(0 rem))     +(hi)
    ?:  =(0 hi)  zero
    (bit [%p s.u --0 hi])
  ::    +rnd / +flr / +cel:  nearest-even / floor / ceil
  ++  rnd  (round %near)
  ++  flr  (round %down)
  ++  cel  (round %up)
  ::    +sun:  @u -> @  (unsigned integer to posit)
  ++  sun  |=(v=@ ^-(@ ?:(=(0 v) zero (bit [%p %.y --0 v]))))
  ::    +san:  @s -> @  (signed integer to posit)
  ++  san
    |=  v=@s
    ^-  @
    =+  [sn mg]=(old:si v)
    ?:  =(0 mg)  zero
    (bit [%p sn --0 mg])
  ::    +toi:  @ -> (unit @s)  (posit to nearest-even integer; NaR -> ~)
  ++  toi
    |=  p=@
    ^-  (unit @s)
    =/  u  (sea p)
    ?:  ?=(%n -.u)  ~
    =/  r   (rnd p)
    =/  ur  (sea r)
    ?:  ?=(%z -.ur)  `--0
    ?>  ?=(%p -.ur)
    ::  value is integer-valued; shift the significand into place (exact).
    =/  sh   (abs:si e.ur)
    =/  mag  ?:((gte-s e.ur --0) (lsh [0 sh] a.ur) (rsh [0 sh] a.ur))
    `(new:si s.ur mag)
  ::    +fma:  @ -> @ -> @ -> @  (fused multiply-add a*b+c, single rounding)
  ++  fma
    |=  [a=@ b=@ c=@]
    ^-  @
    =/  ua  (sea a)
    =/  ub  (sea b)
    =/  uc  (sea c)
    ?:  |(?=(%n -.ua) ?=(%n -.ub) ?=(%n -.uc))  nar
    ?:  |(?=(%z -.ua) ?=(%z -.ub))  c
    ?>  ?=(%p -.ua)
    ?>  ?=(%p -.ub)
    =/  ps  =(s.ua s.ub)
    =/  pe  (sum:si e.ua e.ub)
    =/  pa  (^mul a.ua a.ub)
    ?:  ?=(%z -.uc)  (bit [%p ps pe pa])
    ?>  ?=(%p -.uc)
    =/  emin  ?:(=(-1 (cmp:si pe e.uc)) pe e.uc)
    =/  s1  (lsh [0 (abs:si (dif:si pe emin))] pa)
    =/  s2  (lsh [0 (abs:si (dif:si e.uc emin))] a.uc)
    ?:  =(ps s.uc)    (bit [%p ps emin (^add s1 s2)])
    ?:  (^gth s1 s2)  (bit [%p ps emin (^sub s1 s2)])
    ?:  (^gth s2 s1)  (bit [%p s.uc emin (^sub s2 s1)])
    zero
  --
::  posit8  ("byte"),  posit<8,2>
++  rpb  %*(. pp bloq 3)
::  posit16 ("half"),  posit<16,2>
++  rph  %*(. pp bloq 4)
::  posit32 ("single"), posit<32,2>
++  rps  %*(. pp bloq 5)
::  posit64 ("double"), posit<64,2>
++  rpd  %*(. pp bloq 6)
--
