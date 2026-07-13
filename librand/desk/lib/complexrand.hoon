/+  rand, i754rand, complex, math
::::  /lib/complexrand -- complex-number output adapter (rand-spec.md
::::  section 12.3)
::
::  One arm set per component-width door, matching /lib/complex's own
::  ship order: `++cd` (double, @cd) first, `++cs` (single, @cs) second.
::  Every arm here draws floats directly (there is no raw-bit passthrough
::  the way twocrand/fixedrand have), so -- unlike those two -- this
::  adapter DOES depend on /lib/i754rand (for +rd:uni/+rd-oo:uni and
::  ++dist's +normal), plus /lib/math (for +cos/+sin/+tau/+invsqt2) and
::  /lib/complex (for +pak/+re/+im, to assemble/inspect the packed atom).
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::    +cd:  complex-double (@cd) adapters.
::
::  Uses /lib/math's rd door and /lib/i754rand's @rd generators
::  throughout.
::  Source
++  cd
  |%
  ++  m  ~(. rd:math [%n .~1e-13 .~0])
  ::    +cuniform:  rng -> [@cd rng]
  ::
  ::  Uniform over the unit square [0,1) x [0,1): two independent [0,1)
  ::  component draws, packed.
  ::    Examples
  ::      > `@ux`out:(cuniform:cd:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3fe8.9e6a.a1b9.65f4.3f95.072f.63b9.b5e0
  ::  Source
  ++  cuniform
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  u  r  (rd:uni:i754rand r)
    =^  v  r  (rd:uni:i754rand r)
    [(~(pak cd:complex %n) u v) r]
  ::    +normal-parts:  rng -> [@cd rng]
  ::
  ::  Independent N(0,1) real and imaginary components (two draws from
  ::  /lib/i754rand's +normal, one per component).
  ::    Examples
  ::      > `@ux`out:(normal-parts:cd:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3ff2.e308.2bd4.88f8.bfee.55f7.4af6.d8b9
  ::  Source
  ++  normal-parts
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  u  r  (normal:rd:dist:i754rand r)
    =^  v  r  (normal:rd:dist:i754rand r)
    [(~(pak cd:complex %n) u v) r]
  ::    +cnormal:  rng -> [@cd rng]
  ::
  ::  Standard circularly-symmetric complex normal CN(0,1): components
  ::  N(0, 1/2), i.e. +normal-parts scaled by 1/sqrt(2).  Deliberately
  ::  distinct from +normal-parts (N(0,1) components) -- "complex
  ::  Gaussian" means different things to signal-processing users (who
  ::  usually mean CN(0,1), total variance 1) and statistics users (who
  ::  usually mean iid N(0,1) parts, total variance 2).  Both are named
  ::  unambiguously so neither reading silently gets the wrong one.
  ::    Examples
  ::      > `@ux`out:(cnormal:cd:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3fea.b5c4.88e8.bf41.bfe5.735e.01a0.2b7b
  ::  Source
  ++  cnormal
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  z  r  (normal-parts r)
    =/  u  (mul:m invsqt2:m (~(re cd:complex %n) z))
    =/  v  (mul:m invsqt2:m (~(im cd:complex %n) z))
    [(~(pak cd:complex %n) u v) r]
  ::    +on-circle:  rng -> [@cd rng]
  ::
  ::  Uniform on the unit circle: theta = tau*u for u drawn uniform
  ::  [0,1), components (cos theta, sin theta) via /lib/math's
  ::  bit-specified @rd kernels (jet parity holds, unlike this adapter's
  ::  own naive Taylor loops it borrows nothing from).
  ::    Examples
  ::      > `@ux`out:(on-circle:cd:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3fc0.7838.f0db.1ef6.3fef.bbe7.a757.e0e1  ::  re^2+im^2 = 1
  ::  Source
  ++  on-circle
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  u  r  (rd:uni:i754rand r)
    =/  theta  (mul:m tau:m u)
    [(~(pak cd:complex %n) (cos:m theta) (sin:m theta)) r]
  ::    +in-disk:  rng -> [@cd rng]
  ::
  ::  Uniform on the unit disk via rejection from the enclosing square:
  ::  draw x,y uniform on (-1,1) (the same +rd:uni-mapped-2x-1 trick
  ::  +normal:rd:dist:i754rand uses for its own rejection setup), reject
  ::  if x^2+y^2 >= 1, else accept [x y] as the point.  Variable-
  ::  consumption (rejection loop): acceptance probability pi/4 (~78.5%),
  ::  astronomically bounded in practice -- same rand-spec.md section 5
  ::  policy as every other rejection-based arm in this codebase.
  ::    Examples
  ::      > `@ux`out:(in-disk:cd:complexrand (from-atom:seed:rand %sm64 0))
  ::      0xbfd1.1d5e.36cd.f850.bfe7.45ce.ffed.7562  ::  re^2+im^2 < 1
  ::  Source
  ++  in-disk
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    |-  ^-  [out=@ r=rng:rand]
    =^  u1  r  (rd:uni:i754rand r)
    =^  u2  r  (rd:uni:i754rand r)
    =/  x  (sub:m (mul:m .~2 u1) .~1)
    =/  y  (sub:m (mul:m .~2 u2) .~1)
    ?:  (gte:m (add:m (mul:m x x) (mul:m y y)) .~1)
      $
    [(~(pak cd:complex %n) x y) r]
  --
::    +cs:  complex-single (@cs) adapters.
::
::  The @rs/@cs mirror of +cd, same five arms, same algorithms -- a
::  mechanical re-instantiation over /lib/math's rs door and /lib/
::  i754rand's @rs generators, per rand-spec.md's "each arm exists at
::  @rd (reference) and @rs" convention (already used throughout
::  /lib/i754rand's own ++dist).
::  Source
++  cs
  |%
  ++  m  ~(. rs:math [%n .1e-5 .0])
  ::    +cuniform:  rng -> [@cs rng]  (see +cuniform:cd for algorithm)
  ::    Examples
  ::      > `@ux`out:(cuniform:cs:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3f39.65f4.3dee.6d78
  ::  Source
  ++  cuniform
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  u  r  (rs:uni:i754rand r)
    =^  v  r  (rs:uni:i754rand r)
    [(~(pak cs:complex %n) u v) r]
  ::    +normal-parts:  rng -> [@cs rng]  (see +normal-parts:cd for algorithm)
  ::    Examples
  ::      > `@ux`out:(normal-parts:cs:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x4009.430e.bf17.e80f
  ::  Source
  ++  normal-parts
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  u  r  (normal:rs:dist:i754rand r)
    =^  v  r  (normal:rs:dist:i754rand r)
    [(~(pak cs:complex %n) u v) r]
  ::    +cnormal:  rng -> [@cs rng]  (see +cnormal:cd for algorithm)
  ::    Examples
  ::      > `@ux`out:(cnormal:cs:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3fc2.1e20.bed6.d405
  ::  Source
  ++  cnormal
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  z  r  (normal-parts r)
    =/  u  (mul:m invsqt2:m (~(re cs:complex %n) z))
    =/  v  (mul:m invsqt2:m (~(im cs:complex %n) z))
    [(~(pak cs:complex %n) u v) r]
  ::    +on-circle:  rng -> [@cs rng]  (see +on-circle:cd for algorithm)
  ::    Examples
  ::      > `@ux`out:(on-circle:cs:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3f2b.0087.3f3e.82b8  ::  re^2+im^2 = 1
  ::  Source
  ++  on-circle
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    =^  u  r  (rs:uni:i754rand r)
    =/  theta  (mul:m tau:m u)
    [(~(pak cs:complex %n) (cos:m theta) (sin:m theta)) r]
  ::    +in-disk:  rng -> [@cs rng]  (see +in-disk:cd for algorithm)
  ::    Examples
  ::      > `@ux`out:(in-disk:cs:complexrand (from-atom:seed:rand %sm64 0))
  ::      0x3ee5.97d0.bf44.64a2  ::  re^2+im^2 < 1
  ::  Source
  ++  in-disk
    |=  r=rng:rand
    ^-  [out=@ r=rng:rand]
    |-  ^-  [out=@ r=rng:rand]
    =^  u1  r  (rs:uni:i754rand r)
    =^  u2  r  (rs:uni:i754rand r)
    =/  x  (sub:m (mul:m .2 u1) .1)
    =/  y  (sub:m (mul:m .2 u2) .1)
    ?:  (gte:m (add:m (mul:m x x) (mul:m y y)) .1)
      $
    [(~(pak cs:complex %n) x y) r]
  --
--
