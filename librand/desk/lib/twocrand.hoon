/+  rand, twoc
::::  /lib/twocrand -- two's-complement integer output adapter (rand-spec.md
::::  section 12.2)
::
::  Routes /lib/rand's engines through /lib/twoc's width-keyed two's-
::  complement representation.  No new engine work: both arms are thin
::  compositions of /lib/rand's +bits/+below and /lib/twoc's +twid door.
::  This is /lib/rand's LEANEST adapter -- no floats anywhere, so no
::  dependency on /lib/i754rand or /lib/math, unlike /lib/complexrand.
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::    +twoc-full:  [r=rng w=@] -> [@ rng]
::
::  .w raw bits, interpreted directly as the width-.w two's-complement
::  pattern.  This is the IDENTITY on the raw bits -- a two's-complement
::  value IS just a bit pattern, so a uniform draw over raw bits already
::  IS a uniform draw over two's-complement values at that width.  No
::  encoding step needed (rand-spec.md section 12.2 marks this "Done").
::    Examples
::      > (twoc-full:twocrand (from-atom:seed:rand %sm64 0) 8)
::      [out=175 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
::  Source
++  twoc-full
  |=  [r=rng:rand w=@]
  ^-  [out=@ r=rng:rand]
  (bits:uni:rand r w)
::    +twoc-between:  [r=rng w=@ a=@ b=@] -> [@ rng]
::
::  Inclusive range [a,b], where .a and .b are width-.w two's-complement
::  PATTERNS (raw bits, not @s) and the range is inclusive in TWOC order
::  (negatives below non-negatives, per /lib/twoc's own +gth/+lth).
::  Crashes unless a<=b in that order.
::
::  Bias to unsigned: the span b-a+1 is computed via /lib/twoc's own
::  modular +sub, so it wraps correctly regardless of a/b's individual
::  signs; drawn uniformly via Lemire's +below (rand-spec.md section
::  5.1); biased back via /lib/twoc's modular +add.  /lib/twoc is the
::  primitive /lib/fixedrand's +fixed-between delegates to, per
::  rand-spec.md section 12.1.
::
::  JET WARNING (rand-spec.md section 12.2, verbatim): the span of a
::  full-width-64 range (a=minint, b=maxint) needs 65 bits, not 64 --
::  Hoon's arbitrary-precision @ doesn't care (this Hoon implementation
::  is correct at every width with no special case), but a future C jet
::  computing the span in a machine uint64_t WILL overflow for that one
::  case.  The jet must special-case span=2^w, or compute the span in a
::  wider-than-w integer type.
::    Examples
::      > (twoc-between:twocrand (from-atom:seed:rand %sm64 0) 8 0xfb 0x5)
::      [out=4 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
::  Source
++  twoc-between
  |=  [r=rng:rand w=@ a=@ b=@]
  ^-  [out=@ r=rng:rand]
  =/  tw  ~(. twid:twoc w)
  ~|  %rand-bad-range
  ?>  !(gth:tw a b)
  =/  span  +((sub:tw b a))
  =^  off  r  (below:uni:rand r span)
  [(add:tw a off) r]
--
