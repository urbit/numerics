/+  rand, unum
::::  /lib/unumrand -- posit output adapter (rand-spec.md section 12.4)
::
::  Two inequivalent uniform semantics, named so they cannot be confused --
::  posits are TAPERED, so consecutive bit patterns are NOT evenly spaced in
::  value:
::
::    +posit-lattice -- uniform over raw BIT PATTERNS excluding NaR.  The
::    induced VALUE distribution is only approximately log-uniform (dense
::    near +-1, sparse toward the extremes).  A fuzzing/property-testing
::    primitive, not a statistics one.  Ships at all five width doors
::    (rpb/rph/rps/rpd/rpq) -- no bit-count subtlety, exact at any width.
::
::    +posit-unit -- uniform VALUE on [0,1), exact.  Ships at posit8/16/32
::    only (see librand/NEXT-STEPS.md's milestone-7-pass-2 entry for the
::    full k=4n derivation and the scope decision).  Draw k raw bits u,
::    encode the dyadic u*2^-k -- the g-layer value [%p %.y -k u] (%z when
::    u=0) -- via /lib/unum's existing +bit (round-to-nearest-even,
::    saturating; confirmed no rounding-mode parameter needed).  For this
::    to be EXACTLY the round-to-nearest image of continuous uniform (not
::    approximately), k must exceed the finest rounding-cell width anywhere
::    in [0,1) -- which sits near zero, at exponent -(4(n-2)+1), since
::    minpos = 2^-4(n-2).  k=4n (32/64/128 for posit8/16/32) gives a
::    CONSTANT 7-bit safety margin at any width: 4n - (4(n-2)+1) = 7 always.
::    Verified against an exact mpmath-free rational oracle (every quantity
::    here is already dyadic, so plain Fraction arithmetic is exact) in
::    librand/tools/posit_unit_check.py -- see NEXT-STEPS.md for how to run
::    the chi-square check against ship-drawn draws.
::
::  Distributions ("sample at @rd, convert") need NO new /lib/unum plumbing,
::  unlike fixedrand's +from-rd: /lib/unum already ships +from-rh/rs/rd/rq
::  at every width door.  Not shipped as dedicated wrapper arms here, per
::  the same decision as fixedrand/complexrand -- documented composition +
::  one test proving it (tests/lib/unumrand.hoon).
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::  +pl: shared +posit-lattice implementation, keyed on bloq -- mirrors
::  /lib/unum's own +pp idiom (one generic door, thin per-width forwarders
::  below), so the algorithm exists exactly once regardless of width.
++  pl
  |=  [r=rng:rand =bloq]
  ^-  [out=@ r=rng:rand]
  =/  w  ~(. pp:unum bloq)
  |-  ^-  [out=@ r=rng:rand]
  =^  b  r  (bits:uni:rand r n:w)
  ?:  =(b nar:w)
    $
  [b r]
::  +pu: shared +posit-unit implementation, keyed on bloq and the drawn
::  bit-count k.  Only called from rpb/rph/rps below (the widths in scope).
++  pu
  |=  [r=rng:rand =bloq k=@]
  ^-  [out=@ r=rng:rand]
  =/  w  ~(. pp:unum bloq)
  =^  u  r  (bits:uni:rand r k)
  :_  r
  %-  bit:w
  ?:  =(u 0)
    [%z ~]
  [%p %.y (dif:si --0 (sun:si k)) u]
::    +rpb:  posit8 (n=8) adapters.
::  Source
++  rpb
  |%
  ::    +posit-lattice:  rng -> [@ rng]
  ::    Examples
  ::      > (posit-lattice:rpb:unumrand (from-atom:seed:rand %sm64 0))
  ::  Source
  ++  posit-lattice  |=(r=rng:rand (pl r 3))
  ::    +posit-unit:  rng -> [@ rng]
  ::    Examples
  ::      > (posit-unit:rpb:unumrand (from-atom:seed:rand %sm64 0))
  ::  Source
  ++  posit-unit     |=(r=rng:rand (pu r 3 32))
  --
::    +rph:  posit16 (n=16) adapters.
::  Source
++  rph
  |%
  ++  posit-lattice  |=(r=rng:rand (pl r 4))
  ++  posit-unit     |=(r=rng:rand (pu r 4 64))
  --
::    +rps:  posit32 (n=32) adapters.
::  Source
++  rps
  |%
  ++  posit-lattice  |=(r=rng:rand (pl r 5))
  ++  posit-unit     |=(r=rng:rand (pu r 5 128))
  --
::    +rpd:  posit64 (n=64) adapter -- +posit-lattice only (see header).
::  Source
++  rpd
  |%
  ++  posit-lattice  |=(r=rng:rand (pl r 6))
  --
::    +rpq:  posit128 (n=128) adapter -- +posit-lattice only (see header).
::  Source
++  rpq
  |%
  ++  posit-lattice  |=(r=rng:rand (pl r 7))
  --
--
