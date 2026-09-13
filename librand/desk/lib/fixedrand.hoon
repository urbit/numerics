/+  rand, twocrand, fx=fixed
::::  /lib/fixedrand -- fixed-point output adapter (rand-spec.md section
::::  12.1)
::
::  Routes /lib/rand's engines through /lib/fixed's two's-complement
::  fixed-point representation.  No new engine work, and no floats: a
::  fixed-point lattice is UNIFORM, so uniform generation is exact by
::  construction (stronger than the float case in rand-spec.md section
::  5.2 -- no rounding argument needed here at all).
::
::  Imports /lib/fixed aliased to .fx, not .fixed: this file's own
::  +fixed arm (named to match rand-spec.md's public API) would
::  otherwise shadow the library face of the same name within its own
::  body and every later sibling arm's -- the exact footgun /lib/fixed
::  itself hit and fixed by renaming +twoc to +neg (see its own +neg
::  doc comment).  Here the arm name is the spec-mandated one, so the
::  import is what gets renamed instead.
::
::  Distributions (rand-spec.md 12.1: "sample at @rd, quantize round-to-
::  nearest") are NOT shipped as dedicated per-distribution wrapper arms
::  here -- that would be a dozen-plus nearly-identical one-liners doing
::  nothing beyond composing two already-existing arms.  The composition
::  is exactly: draw a value from /lib/i754rand's ++dist, then quantize
::  via /lib/fixed's +from-rd (added alongside this file, mirroring
::  +from-rs, since it didn't exist yet).  See tests/lib/fixedrand.hoon
::  for that composition proven end to end, and librand/NEXT-STEPS.md
::  for the reasoning.
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::    +fixed:  [r=rng p=prec:fx] -> [@ rng]
::
::  Full-range uniform: draw N=(wid p) raw bits, return as the two's-
::  complement pattern at that width.  Like +twoc-full, this is the
::  IDENTITY on the raw bits -- a fixed-point value's on-the-wire
::  representation at a given precision IS just a two's-complement bit
::  pattern, so a uniform draw over raw bits already is a uniform draw
::  over fixed-point values at that precision.
::    Examples
::      > (fixed:fixedrand (from-atom:seed:rand %sm64 0) [8 8])
::      [out=118.191 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
::  Source
++  fixed
  |=  [r=rng:rand p=prec:fx]
  ^-  [out=@ r=rng:rand]
  (bits:uni:rand r (wid:fx p))
::    +fixed-unit:  [r=rng p=prec:fx] -> [@ rng]
::
::  Uniform in [0,1): draw .b.p bits as the fraction field, with the
::  integer field (and sign) forced to zero.  Since the upper bits of a
::  Hoon atom are implicitly zero, drawing exactly .b.p bits already
::  produces the correct full-width pattern -- no shifting or masking
::  needed on top.
::    Examples
::      > (fixed-unit:fixedrand (from-atom:seed:rand %sm64 0) [8 8])
::      [out=175 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
::  Source
++  fixed-unit
  |=  [r=rng:rand p=prec:fx]
  ^-  [out=@ r=rng:rand]
  (bits:uni:rand r b.p)
::    +fixed-between:  [r=rng p=prec:fx a=@ b=@] -> [@ rng]
::
::  Inclusive range [a,b], where .a and .b are width-N fixed-point
::  PATTERNS (N = wid p) in the SAME sense +twoc-between's .a/.b are raw
::  two's-complement patterns, not decoded values -- convert numeric
::  bounds to patterns yourself first (e.g. via /lib/fixed's +from-s/
::  +from-rd) if that's what you have.  A thin delegation to
::  +twoc-between at width N: fixed-point's bit-level representation IS
::  two's-complement, so once the width matches there's nothing left to
::  "reinterpret" -- the same bits mean a fixed-point value instead of a
::  plain integer purely by the caller's own convention.  Crashes unless
::  a<=b in twoc order (see +twoc-between).
::    Examples
::      > (fixed-between:fixedrand (from-atom:seed:rand %sm64 0) [8 8] 0x1.ff00 0x300)
::      ::  range [-1.0, 3.0] in q8.8
::      [out=649 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
::  Source
++  fixed-between
  |=  [r=rng:rand p=prec:fx a=@ b=@]
  ^-  [out=@ r=rng:rand]
  (twoc-between:twocrand r (wid:fx p) a b)
--
