|%
::
++  twoc
  |_  =bloq
  ::
  ++  len  (bex bloq)
  ++  msb
    |=  a=@
    ?:  (^lth (xeb a) len)
      0
    1
  ++  ones  (dec len)
  ::
  ::
  :: https://gist.github.com/mfuerstenau/ba870a29e16536fdbaba#file-zigzag-encoding-readme-L53
  ::  (i >>> 1) ^ (~(i & 1) + 1)
  ::
  ::  Test cases
  :: > (~(s-to-twoc twoc:twoc 2) --2)
  ::   2
  :: > (~(s-to-twoc twoc:twoc 2) --1)
  ::   1
  :: > (~(s-to-twoc twoc:twoc 2) --8)
  ::  /lib/twoc/hoon:<[19 5].[22 45]>
  :: > (~(s-to-twoc twoc:twoc 3) -14)
  ::  242
  ++  s-to-twoc
    |=  a=@s
    ^-  @
    ?>  (^lte `@`a (dec ~(out fe bloq)))
    %+  mix
      (rsh 0 a)
    (~(sum fe bloq) (not 0 len (dis a 1)) 1)
  ::
  :: 1001 is -7 in int4
  :: 1111 - 1001 = 0110
  :: 0110 + 0001 = 0111 (7)
  ++  twoc-to-s
    |=  a=@
    ^-  @s
    ?:  =(1 (msb a))
      (new:si | +((^sub (dec (bex len)) a)))
    (new:si & a)
  ::
  ::
  ::  Arithmetic is MODULAR two's-complement: results wrap mod 2^len rather
  ::  than crashing on overflow.  This DIVERGES from a checked-overflow
  ::  signed-integer type, but matches hardware two's-complement and the
  ::  wrapping behavior of Lagoon's %uint (fe) scalars.  Because the low len
  ::  bits of a sum/product are independent of sign, the same bit-level add
  ::  and multiply serve signed and unsigned alike.
  ::
  ++  add  |=([a=@ b=@] ^-(@ (mod (^add a b) len-mod)))
  ::  +neg: two's-complement negation, ~a + 1 (flip the len bits, add one),
  ::  the standard XOR-then-increment.  neg(min) = min (wraps).
  ++  neg  |=(a=@ ^-(@ (mod +((not 0 len (mod a len-mod))) len-mod)))
  ++  sub  |=([a=@ b=@] ^-(@ (add a (neg b))))
  ++  mul  |=([a=@ b=@] ^-(@ (mod (^mul (mod a len-mod) (mod b len-mod)) len-mod)))
  ::  +abs: |a| as a two's-complement value (abs(min) wraps back to min).
  ++  abs  |=(a=@ ^-(@ ?:(=(1 (msb a)) (neg a) (mod a len-mod))))
  ::  +len-mod: 2^len, the modulus (len is the bit width = 2^bloq).
  ++  len-mod  (bex len)
  ::  +overflow: would a + b overflow a CHECKED signed add?  Retained for
  ::  callers that want to detect, rather than wrap; +add no longer uses it.
  ++  overflow
    |=  [a=@ b=@ c=@]
    ?|  &(=(0 (msb c)) =(1 (msb a)) =(1 (msb b)))
        &(=(1 (msb c)) =(0 (msb a)) =(0 (msb b)))
    ==
  ::  +div: signed division, truncating toward zero (C / Hoon `div` style).
  ::  Division by zero crashes (as `div` does).  div(min, -1) wraps to min.
  ++  div
    |=  [a=@ b=@]
    ^-  @
    =/  sa  (msb a)
    =/  sb  (msb b)
    =/  q   (^div (abs a) (abs b))
    ?:  =(sa sb)  (mod q len-mod)            :: same sign -> non-negative
    (neg q)                                  :: opposite sign -> negative
  ::  +rem: signed remainder, sign follows the DIVIDEND (C `%` / truncated).
  ::  a == (add (mul (div a b) b) (rem a b)).
  ++  rem
    |=  [a=@ b=@]
    ^-  @
    =/  r  (mod (abs a) (abs b))
    ?:  =(0 (msb a))  (mod r len-mod)        :: dividend non-negative
    (neg r)                                  :: dividend negative
  ::  +pow: a to a NON-NEGATIVE integer power n (n a raw @, not two's-comp),
  ::  by modular squaring.  pow(a, 0) = 1.
  ++  pow
    |=  [a=@ n=@]
    ^-  @
    =/  base  (mod a len-mod)
    =/  acc   (mod 1 len-mod)
    |-  ^-  @
    ?:  =(0 n)  acc
    =?  acc  =(1 (dis n 1))  (mul acc base)
    $(n (rsh 0 n), base (mul base base))
  ::
  ::
  ++  extend
    |=  a=@
    ^-  @
    ?:  =((msb a) 0)
      0
    (dec (bex len))
  ::
  ::
  ::  > (~(gth twoc:twoc 3) `@ub`0b1010.1111 `@ub`0b1010.1110)
  ::    %.y
  ::
  ::  > (~(gth twoc:twoc 3) 256 126)
  ::    %.n
  ::
  ++  gth
    |=  [a=@ b=@]
    ::
    ::  check for different signs
    ?:  =(1 (mix (msb a) (msb b)))
      ::
      ::  if different, choose the one that is positive
      =(0 (msb a))
    ::
    ::  if signs same, use the default gth
    (^gth a b)
  ::
  ++  lth  |=([a=@ b=@] (gth b a))
  ++  lte  |=([a=@ b=@] !(gth a b))
  ++  gte  |=([a=@ b=@] !(gth b a))
  --
--
