/+  rand, math
::::  /lib/i754rand -- IEEE-754 float generation + distributions
::
::  The porcelain layer on top of /lib/rand's plumbing: uniform floats
::  (++uni: +rs/+rd/+rh/+rq/+rs-oo/+rd-oo), Vose's alias method (++alias,
::  moved here from /lib/rand's ++sample since +draw needs an actual
::  uniform @rd float draw), and the continuous/discrete distributions
::  built on both (++dist).  Named for the numeric domain it serves (IEEE
::  754), mirroring how /lib/twoc/fixed/complex/unum are kept separate
::  from /lib/math -- this is /lib/rand's own "math.hoon".  Non-float
::  output adapters (twocrand, fixedrand, complexrand, unumrand) are this
::  file's siblings, not its dependents: none of them need floats for
::  their own core uniform-generation arms, only (optionally) for a
::  "sample at @rd, then quantize/convert" distribution wrapper -- see
::  NEXT-STEPS.md.
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::::                    ++uni                          ::  (5.2) float generation
::
::  NAMING FOOTGUN (same class as /lib/rand's ++seed +mix): this core
::  defines arms named +rs/+rd/+rh/+rq, which shadow the Hoon standard
::  library's own +rs/+rd/+rh/+rq IEEE-754 doors.  /lib/math has the exact
::  same collision (its own +rs/+rd/+rh/+rq doors wrap the stdlib ones) and
::  handles it the same way: any reference to the real stdlib door from
::  inside ++uni uses ^rs/^rd/^rh/^rq, e.g. `(~(mul ^rs %n) a b)`.  This
::  core also shadows /lib/rand's OWN ++uni (the integer-only one) by
::  name -- they're siblings in DIFFERENT files, not nested, so there is
::  no actual collision, but the float generators below reach the integer
::  primitives they're built on via `bits:uni:rand` (qualified through the
::  imported `rand` face), never bare `bits`.
::
++  uni
  |%
  ::    +rs:  rng -> [@rs rng]
  ::
  ::  Uniform in [0,1): draw 24 bits, multiply by the exact constant
  ::  2^-24.  +sun of a 24-bit unsigned integer and multiplying by an
  ::  exact power of two are BOTH exact IEEE-754 operations (no rounding,
  ::  regardless of mode), so every representable output is a multiple of
  ::  2^-24 -- exactly uniform over that lattice.  0 is possible; 1 is not.
  ::    Examples
  ::      > (rs:uni:i754rand (from-atom:seed:rand %sm64 0))
  ::      [out=.0.1164197 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
  ::  Source
  ++  rs
    |=  r=rng:rand
    ^-  [out=@rs r=rng:rand]
    =^  b  r  (bits:uni:rand r 24)
    [(~(mul ^rs %n) (~(sun ^rs %n) b) `@rs`0x3380.0000) r]
  ::    +rd:  rng -> [@rd rng]
  ::
  ::  Same as +rs at double precision: 53 bits, times the exact constant
  ::  2^-53.
  ::    Source
  ++  rd
    |=  r=rng:rand
    ^-  [out=@rd r=rng:rand]
    =^  b  r  (bits:uni:rand r 53)
    [(~(mul ^rd %n) (~(sun ^rd %n) b) `@rd`0x3ca0.0000.0000.0000) r]
  ::    +rh:  rng -> [@rh rng]
  ::
  ::  Same as +rs at half precision: 11 bits, times the exact constant
  ::  2^-11.
  ::    Source
  ++  rh
    |=  r=rng:rand
    ^-  [out=@rh r=rng:rand]
    =^  b  r  (bits:uni:rand r 11)
    [(~(mul ^rh %n) (~(sun ^rh %n) b) `@rh`0x1000) r]
  ::    +rq:  rng -> [@rq rng]
  ::
  ::  Same as +rs at quad precision: 113 bits, times the exact constant
  ::  2^-113.
  ::    Source
  ++  rq
    |=  r=rng:rand
    ^-  [out=@rq r=rng:rand]
    =^  b  r  (bits:uni:rand r 113)
    [(~(mul ^rq %n) (~(sun ^rq %n) b) `@rq`0x3f8e.0000.0000.0000.0000.0000.0000.0000) r]
  ::    +rs-oo:  rng -> [@rs rng]
  ::
  ::  Uniform in the OPEN interval (0,1): draw the same 24 bits as +rs, but
  ::  construct (2*bits+1) * 2^-25 instead of bits * 2^-24.  With bits
  ::  ranging over [0,2^24), (2*bits+1) ranges over the odd integers in
  ::  [1,2^25), so the result spans (0,1) on a lattice of spacing 2^-24,
  ::  symmetric about 0.5, never touching either endpoint.  Needed by
  ::  Box-Muller/log-based transforms (milestone 5) so log(0) never fires.
  ::  (2*bits+1) is exact @ arithmetic; the multiply by 2^-25 is exact for
  ::  the same reason it is in +rs.
  ::    Source
  ++  rs-oo
    |=  r=rng:rand
    ^-  [out=@rs r=rng:rand]
    =^  b  r  (bits:uni:rand r 24)
    [(~(mul ^rs %n) (~(sun ^rs %n) +((mul 2 b))) `@rs`0x3300.0000) r]
  ::    +rd-oo:  rng -> [@rd rng]
  ::
  ::  Same construction as +rs-oo at double precision: 53 bits, (2*bits+1)
  ::  * 2^-54.
  ::    Source
  ++  rd-oo
    |=  r=rng:rand
    ^-  [out=@rd r=rng:rand]
    =^  b  r  (bits:uni:rand r 53)
    [(~(mul ^rd %n) (~(sun ^rd %n) +((mul 2 b))) `@rd`0x3c90.0000.0000.0000) r]
  --
::
::::                    ++alias                        ::  (6.1) alias method
::
::  Vose's alias method, moved here from /lib/rand's ++sample: +draw needs
::  an actual uniform @rd float draw to decide accept-vs-redirect, so it
::  is not float-free the way the rest of ++sample is.
::
::    $alias-table:  [n=@ prob=(list @rd) alias=(list @ud)]
::
::  Vose's alias method (1991) table: O(n) to build, O(1) to draw from.
::  .prob and .alias are indexed 0..n-1, parallel to the input weight
::  list's own order.
+$  alias-table  [n=@ prob=(list @rd) alias=(list @ud)]
::
++  alias
  |%
  ::  +mth: the shared @rd math door for +build's arithmetic, forced %n
  ::  per section 6.1's "require @rd %n arithmetic only."
  ++  mth  ~(. rd:math [%n .~1e-13 .~0])
  ::    +build:  (list @rd) -> alias-table
  ::
  ::  Vose's O(n) construction.  Pure (no rng).  Weights need not sum to
  ::  1 -- they're normalized here (scaled so their mean is 1, which is
  ::  what Vose's small/large partition compares against).  Weights must
  ::  be non-negative with at least one positive; crashes otherwise.
  ::
  ::  The small/large worklist (the classic place for off-by-one and
  ::  float-compare bugs, per section 6.1): +build-loop pops one
  ::  under-weight index (l, scaled prob < 1) and one over-weight index
  ::  (g, scaled prob >= 1) each iteration. l's own probability becomes
  ::  its final +prob entry, and l's +alias entry is set to g -- so
  ::  drawing l falls through to g whenever the per-draw uniform check
  ::  fails (see +draw).  g gives up exactly the mass l was short by
  ::  (1 - wt.l), and whatever's left of g's weight goes back on
  ::  whichever worklist it now belongs to.  When either worklist runs
  ::  dry, everything left on the other is a floating-point leftover
  ::  (mathematically it should be exactly weight 1, off by rounding) --
  ::  +build-cleanup gives each remaining index probability 1 outright,
  ::  with its own index as a harmless unused alias.
  ::    Examples
  ::      > (build:alias:i754rand ~[.~1 .~1 .~2])
  ::      [n=3 prob=~[.~0.75 .~1 .~0.5] alias=~[2 1 2]]
  ::  Source
  ++  build
    |=  weights=(list @rd)
    ^-  alias-table
    =/  n  (lent weights)
    ~|  %rand-empty-list
    ?>  (gth n 0)
    =/  total  (roll weights add:mth)
    ~|  %rand-bad-prob
    ?>  (gth:mth total .~0)
    =/  scale   (div:mth (sun:mth n) total)
    =/  scaled  (turn weights |=(w=@rd (mul:mth w scale)))
    =/  idxd  ^-  (list [idx=@ud wt=@rd])
      (turn (gulf 0 (dec n)) |=(k=@ud [idx=k wt=(snag k scaled)]))
    =/  parted  (skid idxd |=([idx=@ud wt=@rd] (lth:mth wt .~1)))
    =/  fin  (build-loop -.parted +.parted *(map @ud @rd) *(map @ud @ud))
    :+  n
      (turn (gulf 0 (dec n)) |=(k=@ud (~(got by probm.fin) k)))
    (turn (gulf 0 (dec n)) |=(k=@ud (~(got by aliasm.fin) k)))
  ++  build-loop
    |=  $:  small=(list [idx=@ud wt=@rd])
            large=(list [idx=@ud wt=@rd])
            probm=(map @ud @rd)
            aliasm=(map @ud @ud)
        ==
    ^-  [probm=(map @ud @rd) aliasm=(map @ud @ud)]
    ?:  |(?=(~ small) ?=(~ large))
      (build-cleanup (weld small large) probm aliasm)
    =/  l  i.small
    =/  g  i.large
    =.  probm   (~(put by probm) idx.l wt.l)
    =.  aliasm  (~(put by aliasm) idx.l idx.g)
    =/  pg  (sub:mth (add:mth wt.g wt.l) .~1)
    ::  pg<1: g is now under-weight -> joins small.  pg>=1: g stays
    ::  over-weight -> joins large.  (Reversing these was a real bug
    ::  caught during on-ship verification: prob/alias came out wrong
    ::  for a 3-element table, traced back to g landing in the opposite
    ::  worklist from where its own weight said it belonged.)
    ?:  (lth:mth pg .~1)
      (build-loop [[idx=idx.g wt=pg] t.small] t.large probm aliasm)
    (build-loop t.small [[idx=idx.g wt=pg] t.large] probm aliasm)
  ++  build-cleanup
    |=  $:  rest=(list [idx=@ud wt=@rd])
            probm=(map @ud @rd)
            aliasm=(map @ud @ud)
        ==
    ^-  [probm=(map @ud @rd) aliasm=(map @ud @ud)]
    ?~  rest
      [probm aliasm]
    %=  $
      rest    t.rest
      probm   (~(put by probm) idx.i.rest .~1)
      aliasm  (~(put by aliasm) idx.i.rest idx.i.rest)
    ==
  ::    +draw:  [t=alias-table r=rng] -> [out=@ud rng]
  ::
  ::  O(1): pick an index uniformly, then a coin flip (biased by
  ::  .prob.t at that index) decides whether to keep it or redirect to
  ::  its alias.
  ::    Source
  ++  draw
    |=  [t=alias-table r=rng:rand]
    ^-  [out=@ud r=rng:rand]
    =^  i  r  (below:uni:rand r n.t)
    =^  u  r  (rd:uni r)
    :-  ?:((lth:mth u (snag i prob.t)) i (snag i alias.t))
    r
  --
::
::::                    ++dist                          ::  (6) distributions
::
::  Nested ++rd / ++rs sub-cores, one per precision (see each for detail).
::  Internally forced to round-to-nearest (%n) via +m in each, matching
::  /lib/math's own doors -- callers never see or choose a rounding mode
::  here.  chi2/student-t aren't named in rand-spec.md's milestone punch
::  list but are one-line compositions of gamma/normal with no new
::  machinery, so they land here too rather than waiting on an unscheduled
::  slot.
::
::  Crashes (`?>` with `~|` tags) on out-of-domain parameters, per section
::  9's uniform error policy: domain violations are programmer error, not
::  data, so no NaN-returning "soft" errors.
::
++  dist
  |%
  ::    ++rd:  the reference distributions, at @rd (double precision)
  ::
  ::  Everything rand-spec.md section 6 specifies, at the reference
  ::  precision.  See ++rs below for the routed-through-the-same-
  ::  algorithm single-precision mirror (rand-spec.md: "each arm exists
  ::  at @rd (reference) and @rs").
  ++  rd
    |%
    ::  +m: the shared @rd math door, forced %n, rtol=1e-13 (the same "sane
    ::  default" convergence tolerance Saloon's own +feps uses for @rd).
    ::  Internal helper, mirroring /lib/fixed's +ng pattern.
    ++  m  ~(. rd:math [%n .~1e-13 .~0])
    ::  +mz: same door, forced %z (truncate toward zero) -- used only by
    ::  +geometric's ceiling, since +toi respects the door's rounding mode
    ::  and %n would round instead of truncate.
    ++  mz  ~(. rd:math [%z .~1e-13 .~0])
    ::    +normal:  rng -> [@rd rng]
    ::
    ::  Standard normal N(0,1) via the Marsaglia polar method (the Box-Muller
    ::  polar variant): draw u,v uniform on (-1,1) (via +rd:uni mapped
    ::  2x-1), reject if s=u^2+v^2 is >=1 or =0, else return u*sqrt(-2 ln(s)/s).
    ::  Returns ONE deviate and discards the pair-mate v*sqrt(-2 ln(s)/s) --
    ::  caching it would make the rng state opaque (the next +normal call
    ::  would need to "remember" a pending mate outside the plain +$rng
    ::  noun), so v's factor is thrown away at v1.  A documented waste, not a
    ::  bug; Ziggurat is a NEXT-STEPS optimization.  Variable-consumption
    ::  (rejection loop): expected iterations ~1.27 (rejection probability
    ::  1 - pi/4), astronomically bounded in practice.
    ::    Examples
    ::      > (normal:rd:dist:i754rand (from-atom:seed:rand %sm64 0))
    ::      [out=.~-0.9479938949723624 r=[%sm64 s=0x78dd.e6e5.fd29.f054]]
    ::  Source
    ++  normal
      |=  r=rng:rand
      ^-  [out=@rd r=rng:rand]
      |-  ^-  [out=@rd r=rng:rand]
      =^  u1  r  (rd:uni r)
      =^  u2  r  (rd:uni r)
      =/  u  (sub:m (mul:m .~2 u1) .~1)
      =/  v  (sub:m (mul:m .~2 u2) .~1)
      =/  s  (add:m (mul:m u u) (mul:m v v))
      ?:  |((gte:m s .~1) (equ:m s .~0))
        $
      =/  factor  (sqt:m (div:m (mul:m .~-2 (log:m s)) s))
      [(mul:m u factor) r]
    ::    +normal-mv:  [r=rng mu=@rd sigma=@rd] -> [@rd rng]
    ::
    ::  N(mu, sigma^2): mu + sigma*z where z ~ N(0,1).  Crashes if sigma < 0.
    ::    Source
    ++  normal-mv
      |=  [r=rng:rand mu=@rd sigma=@rd]
      ^-  [out=@rd r=rng:rand]
      ~|  %rand-bad-sigma
      ?>  !(lth:m sigma .~0)
      =^  z  r  (normal r)
      [(add:m mu (mul:m sigma z)) r]
    ::    +expon:  [r=rng lambda=@rd] -> [@rd rng]
    ::
    ::  Exponential(lambda) via inversion: -ln(u)/lambda, u drawn from the
    ::  OPEN (0,1) (+rd-oo, not +rd) specifically so log(0) never fires --
    ::  this is exactly the case rand-spec.md section 5.2 built +rd-oo for.
    ::  Crashes if lambda <= 0.
    ::    Source
    ++  expon
      |=  [r=rng:rand lambda=@rd]
      ^-  [out=@rd r=rng:rand]
      ~|  %rand-bad-rate
      ?>  (gth:m lambda .~0)
      =^  u  r  (rd-oo:uni r)
      [(div:m (neg:m (log:m u)) lambda) r]
    ::    +gamma:  [r=rng alpha=@rd] -> [@rd rng]
    ::
    ::  Gamma(alpha, scale=1) via Marsaglia-Tsang (2000).  alpha>=1 direct
    ::  (+gamma-ge1); alpha<1 via the standard boost gamma(alpha) =
    ::  gamma(alpha+1) * u^(1/alpha), u drawn from the open (0,1) so the
    ::  u=0 lattice point (probability 2^-53, not truly 0 as it would be for
    ::  a continuous uniform) never manufactures a spurious exact-zero
    ::  sample.  Crashes if alpha <= 0.
    ::    Source
    ++  gamma
      |=  [r=rng:rand alpha=@rd]
      ^-  [out=@rd r=rng:rand]
      ~|  %rand-bad-shape
      ?>  (gth:m alpha .~0)
      ?:  (gte:m alpha .~1)
        (gamma-ge1 r alpha)
      =^  g  r  (gamma-ge1 r (add:m alpha .~1))
      =^  u  r  (rd-oo:uni r)
      [(mul:m g (pow:m u (div:m .~1 alpha))) r]
    ::  +gamma-ge1: Marsaglia-Tsang squeeze for alpha>=1.  d=alpha-1/3,
    ::  c=1/sqrt(9d); draw x~N(0,1), v=(1+cx)^3 (reject if v<=0), draw
    ::  u~(0,1) open (so log(u) never fires on 0), accept d*v if
    ::  ln(u) < x^2/2 + d - d*v + d*ln(v), else reject and redraw both x,u.
    ::  Variable-consumption (rejection loop), astronomically bounded.
    ++  gamma-ge1
      |=  [r=rng:rand alpha=@rd]
      ^-  [out=@rd r=rng:rand]
      =/  d  (sub:m alpha (div:m .~1 .~3))
      =/  c  (div:m .~1 (sqt:m (mul:m .~9 d)))
      |-  ^-  [out=@rd r=rng:rand]
      =^  x  r  (normal r)
      =/  t  (add:m .~1 (mul:m c x))
      =/  v  (mul:m t (mul:m t t))
      ?:  !(gth:m v .~0)
        $
      =^  u  r  (rd-oo:uni r)
      =/  rhs
        %+  add:m
          (add:m (mul:m .~0.5 (mul:m x x)) d)
        (sub:m (mul:m d (log:m v)) (mul:m d v))
      ?:  (lth:m (log:m u) rhs)
        [(mul:m d v) r]
      $
    ::    +beta:  [r=rng a=@rd b=@rd] -> [@rd rng]
    ::
    ::  Beta(a,b) via two independent gammas: x/(x+y), x~gamma(a), y~gamma(b).
    ::  Crashes if a <= 0 or b <= 0 (via +gamma's own precondition).
    ::    Source
    ++  beta
      |=  [r=rng:rand a=@rd b=@rd]
      ^-  [out=@rd r=rng:rand]
      =^  x  r  (gamma r a)
      =^  y  r  (gamma r b)
      [(div:m x (add:m x y)) r]
    ::    +chi2:  [r=rng k=@] -> [@rd rng]
    ::
    ::  Chi-squared with k degrees of freedom: gamma(k/2, scale=2) ==
    ::  2*gamma(k/2, scale=1) (+gamma is scale=1, so the factor of 2 is
    ::  applied directly -- a standard gamma scaling property).  Crashes if
    ::  k = 0.
    ::    Source
    ++  chi2
      |=  [r=rng:rand k=@]
      ^-  [out=@rd r=rng:rand]
      ~|  %rand-bad-df
      ?>  !=(k 0)
      =^  x  r  (gamma r (div:m (sun:m k) .~2))
      [(mul:m .~2 x) r]
    ::    +student-t:  [r=rng k=@] -> [@rd rng]
    ::
    ::  Student's t with k degrees of freedom: z / sqrt(chi2(k)/k), z~N(0,1).
    ::  Crashes if k = 0.
    ::    Source
    ++  student-t
      |=  [r=rng:rand k=@]
      ^-  [out=@rd r=rng:rand]
      ~|  %rand-bad-df
      ?>  !=(k 0)
      =^  z  r  (normal r)
      =^  c  r  (chi2 r k)
      [(div:m z (sqt:m (div:m c (sun:m k)))) r]
    ::    +bernoulli:  [r=rng p=@rd] -> [? rng]
    ::
    ::  %.y with probability p, else %.n: draw u~[0,1), return u<p.  Crashes
    ::  unless 0 <= p <= 1.
    ::    Source
    ++  bernoulli
      |=  [r=rng:rand p=@rd]
      ^-  [out=? r=rng:rand]
      ~|  %rand-bad-prob
      ?>  &((gte:m p .~0) (lte:m p .~1))
      =^  u  r  (rd:uni r)
      [(lth:m u p) r]
    ::    +geometric:  [r=rng p=@rd] -> [@ud rng]
    ::
    ::  Number of Bernoulli(p) trials up to and including the first success:
    ::  ceil(ln(u)/ln(1-p)), u drawn from the open (0,1) (avoids ln(0)).
    ::  p=1 is a special-cased edge that returns 1 directly (ln(1-p) would
    ::  divide by ln(0) = -inf otherwise).  ceil is computed via +toi under
    ::  a SEPARATE %z-rounding (truncate-toward-zero) door instance, +mz --
    ::  +toi respects the door's own rounding mode, and this core's shared
    ::  +m door is forced %n (round-to-nearest), which would round instead
    ::  of truncate.  Crashes unless 0 < p <= 1.
    ::    Source
    ++  geometric
      |=  [r=rng:rand p=@rd]
      ^-  [out=@ud r=rng:rand]
      ~|  %rand-bad-prob
      ?>  &((gth:m p .~0) (lte:m p .~1))
      ?:  (equ:m p .~1)
        [1 r]
      =^  u  r  (rd-oo:uni r)
      =/  raw  (div:m (log:m u) (log:m (sub:m .~1 p)))
      =/  fl   (abs:si (need (toi:mz raw)))
      [?:(=(raw (sun:m fl)) fl +(fl)) r]
    ::    +categorical:  [r=rng t=alias-table] -> [@ud rng]
    ::
    ::  One draw from a pre-built alias table (rand-spec.md section 6.1's
    ::  ++alias, a sibling core in this file).  Deliberately takes an
    ::  ALREADY-BUILT table,
    ::  not raw weights: alias's whole point is O(n) build, O(1) per draw,
    ::  amortized over many draws from the same distribution -- rebuilding
    ::  the table on every single draw would defeat that.  Build once via
    ::  `(build:alias weights)`, draw many times via this arm (or
    ::  `draw:alias` directly, which this is a thin rename of).
    ::  Lives ONLY here, not mirrored in ++rs: alias-table's .prob is fixed
    ::  at @rd (rand-spec.md section 6.1), so this arm is precision-
    ::  invariant -- an "@rs categorical" would be byte-for-byte identical
    ::  code, not a real variant.
    ::    Source
    ++  categorical
      |=  [r=rng:rand t=alias-table]
      ^-  [out=@ud r=rng:rand]
      (draw:alias t r)
    ::    +poisson:  [r=rng lambda=@rd] -> [@ud rng]
    ::
    ::  Poisson(lambda): Knuth's product method for lambda<10 (simple,
    ::  O(lambda) expected multiplications -- fine at this scale, unusable
    ::  above it); Hörmann's PTRS (1993) transformed rejection for
    ::  lambda>=10 (O(1) expected, needed because Knuth's method's cost
    ::  scales linearly with lambda and would silently become the wrong
    ::  choice above the threshold rather than crashing -- so this arm
    ::  switches automatically rather than leaving the choice to the
    ::  caller).  Crashes if lambda <= 0.
    ::    Source
    ++  poisson
      |=  [r=rng:rand lambda=@rd]
      ^-  [out=@ud r=rng:rand]
      ~|  %rand-bad-rate
      ?>  (gth:m lambda .~0)
      ?:  (lth:m lambda .~10)
        (poisson-knuth r lambda)
      (poisson-ptrs r lambda)
    ::  +poisson-knuth: k=0, p=1; loop: k+=1, p*=uniform(0,1); accept
    ::  (return k-1, i.e. the count BEFORE the last increment) once
    ::  p <= exp(-lambda).  Variable-consumption (rejection-free here, but
    ::  a variable number of multiplications -- expected lambda of them).
    ++  poisson-knuth
      |=  [r=rng:rand lambda=@rd]
      ^-  [out=@ud r=rng:rand]
      =/  bigl  (exp:m (neg:m lambda))
      =/  k  0
      =/  p  .~1
      |-  ^-  [out=@ud r=rng:rand]
      =^  u  r  (rd:uni r)
      =/  p2  (mul:m p u)
      ?:  (lte:m p2 bigl)
        [k r]
      $(k +(k), p p2)
    ::  +poisson-ptrs: Hörmann 1993's transformed rejection with squeeze,
    ::  verified against NumPy's random_poisson_ptrs (src/distributions/
    ::  distributions.c) -- an independent re-implementation of the exact
    ::  same reference this spec cites, not reconstructed from memory.
    ::  Draws a candidate k from a scaled/shifted uniform (the "transformed"
    ::  part), accepts immediately if a cheap squeeze test passes (avoiding
    ::  the exact log-probability computation on the common path), and
    ::  falls back to the exact accept/reject test (via +loggam, log(k!))
    ::  otherwise.  Rejection loop, astronomically bounded in practice.
    ++  poisson-ptrs
      |=  [r=rng:rand lambda=@rd]
      ^-  [out=@ud r=rng:rand]
      =/  slam      (sqt:m lambda)
      =/  loglam    (log:m lambda)
      =/  b         (add:m .~0.931 (mul:m .~2.53 slam))
      =/  a         (add:m .~-0.059 (mul:m .~0.02483 b))
      =/  invalpha  (add:m .~1.1239 (div:m .~1.1328 (sub:m b .~3.4)))
      =/  vr        (sub:m .~0.9277 (div:m .~3.6224 (sub:m b .~2)))
      |-  ^-  [out=@ud r=rng:rand]
      =^  u1  r  (rd:uni r)
      =^  v   r  (rd:uni r)
      =/  u   (sub:m u1 .~0.5)
      =/  us  (sub:m .~0.5 (abs:m u))
      =/  kk  (dfloor (add:m (add:m (mul:m (add:m (div:m (mul:m .~2 a) us) b) u) lambda) .~0.43))
      ?:  &((gte:m us .~0.07) (lte:m v vr))
        [(abs:si kk) r]
      ?:  |(=(-1 (cmp:si kk --0)) &((lth:m us .~0.013) (gth:m v us)))
        $
      =/  lhs  (sub:m (add:m (log:m v) (log:m invalpha)) (log:m (add:m (div:m a (mul:m us us)) b)))
      =/  rhs  (sub:m (add:m (mul:m .~-1 lambda) (mul:m (san:m kk) loglam)) (loggam (add:m (san:m kk) .~1)))
      ?:  (lte:m lhs rhs)
        [(abs:si kk) r]
      $
    ::  +dfloor: @rd -> @s, floor (round toward -infinity, NOT toward
    ::  zero) -- +poisson-ptrs's candidate k can be negative, where
    ::  truncate-toward-zero (what +toi under +mz gives) would round the
    ::  wrong way.  Positive/exact values: truncation IS the floor.
    ::  Negative non-integers: truncation rounds up, so subtract 1.
    ++  dfloor
      |=  x=@rd
      ^-  @s
      =/  t  (need (toi:mz x))
      ?:  |((gte:m x .~0) (equ:m x (san:m t)))
        t
      (dif:si t --1)
    ::  +loggam: @rd -> @rd, natural log of the gamma function (so log(k!)
    ::  = loggam(k+1)).  Stirling's series with a small-argument recurrence
    ::  for x<7 (shifts x up by whole integers until >=7, where the series
    ::  is accurate, then subtracts the shifted-away log-factors back out)
    ::  -- transcribed from NumPy's random_loggam (src/distributions/
    ::  distributions.c), itself from Zhang & Jin's SPECFUN algorithm, not
    ::  derived independently: log-gamma accuracy is exactly the kind of
    ::  detail worth matching a vetted reference bit-for-bit rather than
    ::  reconstructing from the general shape of Stirling's series.
    ++  loggam
      |=  x=@rd
      ^-  @rd
      ?:  |((equ:m x .~1) (equ:m x .~2))
        .~0
      =/  a  ^-  (list @rd)
        :~  .~0.08333333333333333  .~-0.002777777777777778  .~0.0007936507936507937
            .~-0.0005952380952380952  .~0.0008417508417508417  .~-0.001917526917526918
            .~0.006410256410256411  .~-0.02955065359477124  .~0.1796443723688307
            .~-1.39243221690590
        ==
      =/  lg2pi  .~1.8378770664093453
      =/  n  ?:((lth:m x .~7) (abs:si (need (toi:mz (sub:m .~7 x)))) 0)
      =/  x0  (add:m x (sun:m n))
      =/  x2  (div:m .~1 (mul:m x0 x0))
      =/  gl0
        =/  acc  (snag 9 a)
        =/  i  0
        |-  ^-  @rd
        ?:  =(i 9)
          acc
        $(acc (add:m (mul:m acc x2) (snag (sub 8 i) a)), i +(i))
      =/  gl
        %+  sub:m
          %+  add:m
            (add:m (div:m gl0 x0) (mul:m .~0.5 lg2pi))
          (mul:m (sub:m x0 .~0.5) (log:m x0))
        x0
      ?.  (lth:m x .~7)
        gl
      =/  fin
        =/  gl2  gl
        =/  x02  x0
        =/  k  1
        |-  ^-  @rd
        ?:  (gth k n)
          gl2
        $(gl2 (sub:m gl2 (log:m (sub:m x02 .~1))), x02 (sub:m x02 .~1), k +(k))
      fin
    ::    +binomial:  [r=rng n=@ p=@rd] -> [@ud rng]
    ::
    ::  Binomial(n,p) via inversion by recurrence: P(X=0)=(1-p)^n, then
    ::  P(X=k+1) = P(X=k) * (n-k)/(k+1) * p/(1-p), accumulating the CDF
    ::  until it exceeds a single uniform draw.  Only valid for
    ::  n*min(p,1-p) < 30 (Kachitvichyanukul & Schmeiser 1988) -- above
    ::  that, inversion needs too many terms and BTPE is the right
    ::  algorithm, deferred to NEXT-STEPS; this arm crashes rather than
    ::  running slow (many terms) or, worse, silently truncating.
    ::    Source
    ++  binomial
      |=  [r=rng:rand n=@ p=@rd]
      ^-  [out=@ud r=rng:rand]
      ~|  %rand-bad-prob
      ?>  &((gte:m p .~0) (lte:m p .~1))
      =/  q     (sub:m .~1 p)
      =/  qmin  ?:((lth:m p .~0.5) p q)
      ~|  %rand-binomial-n-too-large
      ?>  (lth:m (mul:m (sun:m n) qmin) .~30)
      =^  u  r  (rd:uni r)
      =/  f    (pow:m q (sun:m n))
      =/  cdf  f
      =/  k    0
      |-  ^-  [out=@ud r=rng:rand]
      ?:  (lth:m u cdf)
        [k r]
      =/  f2  %+  mul:m  f
        (mul:m (div:m (sun:m (sub n k)) (sun:m +(k))) (div:m p q))
      $(f f2, cdf (add:m cdf f2), k +(k))
    ::    +dirichlet:  [r=rng alphas=(list @rd)] -> [(list @rd) rng]
    ::
    ::  Dirichlet(alphas): draw one gamma(alpha_i) per entry, normalize by
    ::  their sum.  Crashes on an empty .alphas (and, via +gamma's own
    ::  precondition, on any non-positive alpha).
    ::    Source
    ++  dirichlet
      |=  [r=rng:rand alphas=(list @rd)]
      ^-  [out=(list @rd) r=rng:rand]
      ~|  %rand-empty-list
      ?>  !=(~ alphas)
      =/  gd
        |-  ^-  [(list @rd) rng:rand]
        ?~  alphas
          [~ r]
        =^  x     r  (gamma r i.alphas)
        =^  rest  r  $(alphas t.alphas)
        [[x rest] r]
      =/  total  (roll -.gd add:m)
      [(turn -.gd |=(g=@rd (div:m g total))) +.gd]
  --
  ::    ++rs:  the same distributions, routed through @rs (single
  ::  precision)
  ::
  ::  Identical algorithms to ++rd, just built on ++uni's @rs generators
  ::  and /lib/math's @rs door instead of @rd's -- not a separate
  ::  derivation, a mechanical re-instantiation at the other precision
  ::  (mirroring how /lib/math itself has independent per-precision
  ::  doors rather than one generic implementation).  +categorical has
  ::  no @rs mirror: alias-table's .prob is fixed at @rd, so it's
  ::  precision-invariant (see its doccord in ++rd).
  ++  rs
    |%
    ::  +m: the shared @rs math door, forced %n, rtol=1e-6 (the same "sane
    ::  default" convergence tolerance Saloon's own +feps uses for @rs).
    ::  Internal helper, mirroring /lib/fixed's +ng pattern.
    ++  m  ~(. rs:math [%n .1e-6 .0])
    ::  +mz: same door, forced %z (truncate toward zero) -- used only by
    ::  +geometric's ceiling, since +toi respects the door's rounding mode
    ::  and %n would round instead of truncate.
    ++  mz  ~(. rs:math [%z .1e-6 .0])
    ::    +normal:  rng -> [@rs rng]
    ::
    ::  Standard normal N(0,1) via the Marsaglia polar method (the Box-Muller
    ::  polar variant): draw u,v uniform on (-1,1) (via +rs:uni mapped
    ::  2x-1), reject if s=u^2+v^2 is >=1 or =0, else return u*sqrt(-2 ln(s)/s).
    ::  Returns ONE deviate and discards the pair-mate v*sqrt(-2 ln(s)/s) --
    ::  caching it would make the rng state opaque (the next +normal call
    ::  would need to "remember" a pending mate outside the plain +$rng
    ::  noun), so v's factor is thrown away at v1.  A documented waste, not a
    ::  bug; Ziggurat is a NEXT-STEPS optimization.  Variable-consumption
    ::  (rejection loop): expected iterations ~1.27 (rejection probability
    ::  1 - pi/4), astronomically bounded in practice.
    ::    Examples
    ::      > (normal:rs:dist:i754rand (from-atom:seed:rand %sm64 0))
    ::      [out=.-0.5933847 r=[%sm64 s=0x3c6e.f372.fe94.f82a]]
    ::  Source
    ++  normal
      |=  r=rng:rand
      ^-  [out=@rs r=rng:rand]
      |-  ^-  [out=@rs r=rng:rand]
      =^  u1  r  (rs:uni r)
      =^  u2  r  (rs:uni r)
      =/  u  (sub:m (mul:m .2 u1) .1)
      =/  v  (sub:m (mul:m .2 u2) .1)
      =/  s  (add:m (mul:m u u) (mul:m v v))
      ?:  |((gte:m s .1) (equ:m s .0))
        $
      =/  factor  (sqt:m (div:m (mul:m .-2 (log:m s)) s))
      [(mul:m u factor) r]
    ::    +normal-mv:  [r=rng mu=@rs sigma=@rs] -> [@rs rng]
    ::
    ::  N(mu, sigma^2): mu + sigma*z where z ~ N(0,1).  Crashes if sigma < 0.
    ::    Source
    ++  normal-mv
      |=  [r=rng:rand mu=@rs sigma=@rs]
      ^-  [out=@rs r=rng:rand]
      ~|  %rand-bad-sigma
      ?>  !(lth:m sigma .0)
      =^  z  r  (normal r)
      [(add:m mu (mul:m sigma z)) r]
    ::    +expon:  [r=rng lambda=@rs] -> [@rs rng]
    ::
    ::  Exponential(lambda) via inversion: -ln(u)/lambda, u drawn from the
    ::  OPEN (0,1) (+rd-oo, not +rd) specifically so log(0) never fires --
    ::  this is exactly the case rand-spec.md section 5.2 built +rd-oo for.
    ::  Crashes if lambda <= 0.
    ::    Source
    ++  expon
      |=  [r=rng:rand lambda=@rs]
      ^-  [out=@rs r=rng:rand]
      ~|  %rand-bad-rate
      ?>  (gth:m lambda .0)
      =^  u  r  (rs-oo:uni r)
      [(div:m (neg:m (log:m u)) lambda) r]
    ::    +gamma:  [r=rng alpha=@rs] -> [@rs rng]
    ::
    ::  Gamma(alpha, scale=1) via Marsaglia-Tsang (2000).  alpha>=1 direct
    ::  (+gamma-ge1); alpha<1 via the standard boost gamma(alpha) =
    ::  gamma(alpha+1) * u^(1/alpha), u drawn from the open (0,1) so the
    ::  u=0 lattice point (probability 2^-53, not truly 0 as it would be for
    ::  a continuous uniform) never manufactures a spurious exact-zero
    ::  sample.  Crashes if alpha <= 0.
    ::    Source
    ++  gamma
      |=  [r=rng:rand alpha=@rs]
      ^-  [out=@rs r=rng:rand]
      ~|  %rand-bad-shape
      ?>  (gth:m alpha .0)
      ?:  (gte:m alpha .1)
        (gamma-ge1 r alpha)
      =^  g  r  (gamma-ge1 r (add:m alpha .1))
      =^  u  r  (rs-oo:uni r)
      [(mul:m g (pow:m u (div:m .1 alpha))) r]
    ::  +gamma-ge1: Marsaglia-Tsang squeeze for alpha>=1.  d=alpha-1/3,
    ::  c=1/sqrt(9d); draw x~N(0,1), v=(1+cx)^3 (reject if v<=0), draw
    ::  u~(0,1) open (so log(u) never fires on 0), accept d*v if
    ::  ln(u) < x^2/2 + d - d*v + d*ln(v), else reject and redraw both x,u.
    ::  Variable-consumption (rejection loop), astronomically bounded.
    ++  gamma-ge1
      |=  [r=rng:rand alpha=@rs]
      ^-  [out=@rs r=rng:rand]
      =/  d  (sub:m alpha (div:m .1 .3))
      =/  c  (div:m .1 (sqt:m (mul:m .9 d)))
      |-  ^-  [out=@rs r=rng:rand]
      =^  x  r  (normal r)
      =/  t  (add:m .1 (mul:m c x))
      =/  v  (mul:m t (mul:m t t))
      ?:  !(gth:m v .0)
        $
      =^  u  r  (rs-oo:uni r)
      =/  rhs
        %+  add:m
          (add:m (mul:m .0.5 (mul:m x x)) d)
        (sub:m (mul:m d (log:m v)) (mul:m d v))
      ?:  (lth:m (log:m u) rhs)
        [(mul:m d v) r]
      $
    ::    +beta:  [r=rng a=@rs b=@rs] -> [@rs rng]
    ::
    ::  Beta(a,b) via two independent gammas: x/(x+y), x~gamma(a), y~gamma(b).
    ::  Crashes if a <= 0 or b <= 0 (via +gamma's own precondition).
    ::    Source
    ++  beta
      |=  [r=rng:rand a=@rs b=@rs]
      ^-  [out=@rs r=rng:rand]
      =^  x  r  (gamma r a)
      =^  y  r  (gamma r b)
      [(div:m x (add:m x y)) r]
    ::    +chi2:  [r=rng k=@] -> [@rs rng]
    ::
    ::  Chi-squared with k degrees of freedom: gamma(k/2, scale=2) ==
    ::  2*gamma(k/2, scale=1) (+gamma is scale=1, so the factor of 2 is
    ::  applied directly -- a standard gamma scaling property).  Crashes if
    ::  k = 0.
    ::    Source
    ++  chi2
      |=  [r=rng:rand k=@]
      ^-  [out=@rs r=rng:rand]
      ~|  %rand-bad-df
      ?>  !=(k 0)
      =^  x  r  (gamma r (div:m (sun:m k) .2))
      [(mul:m .2 x) r]
    ::    +student-t:  [r=rng k=@] -> [@rs rng]
    ::
    ::  Student's t with k degrees of freedom: z / sqrt(chi2(k)/k), z~N(0,1).
    ::  Crashes if k = 0.
    ::    Source
    ++  student-t
      |=  [r=rng:rand k=@]
      ^-  [out=@rs r=rng:rand]
      ~|  %rand-bad-df
      ?>  !=(k 0)
      =^  z  r  (normal r)
      =^  c  r  (chi2 r k)
      [(div:m z (sqt:m (div:m c (sun:m k)))) r]
    ::    +bernoulli:  [r=rng p=@rs] -> [? rng]
    ::
    ::  %.y with probability p, else %.n: draw u~[0,1), return u<p.  Crashes
    ::  unless 0 <= p <= 1.
    ::    Source
    ++  bernoulli
      |=  [r=rng:rand p=@rs]
      ^-  [out=? r=rng:rand]
      ~|  %rand-bad-prob
      ?>  &((gte:m p .0) (lte:m p .1))
      =^  u  r  (rs:uni r)
      [(lth:m u p) r]
    ::    +geometric:  [r=rng p=@rs] -> [@ud rng]
    ::
    ::  Number of Bernoulli(p) trials up to and including the first success:
    ::  ceil(ln(u)/ln(1-p)), u drawn from the open (0,1) (avoids ln(0)).
    ::  p=1 is a special-cased edge that returns 1 directly (ln(1-p) would
    ::  divide by ln(0) = -inf otherwise).  ceil is computed via +toi under
    ::  a SEPARATE %z-rounding (truncate-toward-zero) door instance, +mz --
    ::  +toi respects the door's own rounding mode, and this core's shared
    ::  +m door is forced %n (round-to-nearest), which would round instead
    ::  of truncate.  Crashes unless 0 < p <= 1.
    ::    Source
    ++  geometric
      |=  [r=rng:rand p=@rs]
      ^-  [out=@ud r=rng:rand]
      ~|  %rand-bad-prob
      ?>  &((gth:m p .0) (lte:m p .1))
      ?:  (equ:m p .1)
        [1 r]
      =^  u  r  (rs-oo:uni r)
      =/  raw  (div:m (log:m u) (log:m (sub:m .1 p)))
      =/  fl   (abs:si (need (toi:mz raw)))
      [?:(=(raw (sun:m fl)) fl +(fl)) r]
    ::    +poisson:  [r=rng lambda=@rs] -> [@ud rng]
    ::
    ::  Poisson(lambda): Knuth's product method for lambda<10 (simple,
    ::  O(lambda) expected multiplications -- fine at this scale, unusable
    ::  above it); Hörmann's PTRS (1993) transformed rejection for
    ::  lambda>=10 (O(1) expected, needed because Knuth's method's cost
    ::  scales linearly with lambda and would silently become the wrong
    ::  choice above the threshold rather than crashing -- so this arm
    ::  switches automatically rather than leaving the choice to the
    ::  caller).  Crashes if lambda <= 0.
    ::    Source
    ++  poisson
      |=  [r=rng:rand lambda=@rs]
      ^-  [out=@ud r=rng:rand]
      ~|  %rand-bad-rate
      ?>  (gth:m lambda .0)
      ?:  (lth:m lambda .10)
        (poisson-knuth r lambda)
      (poisson-ptrs r lambda)
    ::  +poisson-knuth: k=0, p=1; loop: k+=1, p*=uniform(0,1); accept
    ::  (return k-1, i.e. the count BEFORE the last increment) once
    ::  p <= exp(-lambda).  Variable-consumption (rejection-free here, but
    ::  a variable number of multiplications -- expected lambda of them).
    ++  poisson-knuth
      |=  [r=rng:rand lambda=@rs]
      ^-  [out=@ud r=rng:rand]
      =/  bigl  (exp:m (neg:m lambda))
      =/  k  0
      =/  p  .1
      |-  ^-  [out=@ud r=rng:rand]
      =^  u  r  (rs:uni r)
      =/  p2  (mul:m p u)
      ?:  (lte:m p2 bigl)
        [k r]
      $(k +(k), p p2)
    ::  +poisson-ptrs: Hörmann 1993's transformed rejection with squeeze,
    ::  verified against NumPy's random_poisson_ptrs (src/distributions/
    ::  distributions.c) -- an independent re-implementation of the exact
    ::  same reference this spec cites, not reconstructed from memory.
    ::  Draws a candidate k from a scaled/shifted uniform (the "transformed"
    ::  part), accepts immediately if a cheap squeeze test passes (avoiding
    ::  the exact log-probability computation on the common path), and
    ::  falls back to the exact accept/reject test (via +loggam, log(k!))
    ::  otherwise.  Rejection loop, astronomically bounded in practice.
    ++  poisson-ptrs
      |=  [r=rng:rand lambda=@rs]
      ^-  [out=@ud r=rng:rand]
      =/  slam      (sqt:m lambda)
      =/  loglam    (log:m lambda)
      =/  b         (add:m .0.931 (mul:m .2.53 slam))
      =/  a         (add:m .-0.059 (mul:m .0.02483 b))
      =/  invalpha  (add:m .1.1239 (div:m .1.1328 (sub:m b .3.4)))
      =/  vr        (sub:m .0.9277 (div:m .3.6224 (sub:m b .2)))
      |-  ^-  [out=@ud r=rng:rand]
      =^  u1  r  (rs:uni r)
      =^  v   r  (rs:uni r)
      =/  u   (sub:m u1 .0.5)
      =/  us  (sub:m .0.5 (abs:m u))
      =/  kk  (dfloor (add:m (add:m (mul:m (add:m (div:m (mul:m .2 a) us) b) u) lambda) .0.43))
      ?:  &((gte:m us .0.07) (lte:m v vr))
        [(abs:si kk) r]
      ?:  |(=(-1 (cmp:si kk --0)) &((lth:m us .0.013) (gth:m v us)))
        $
      =/  lhs  (sub:m (add:m (log:m v) (log:m invalpha)) (log:m (add:m (div:m a (mul:m us us)) b)))
      =/  rhs  (sub:m (add:m (mul:m .-1 lambda) (mul:m (san:m kk) loglam)) (loggam (add:m (san:m kk) .1)))
      ?:  (lte:m lhs rhs)
        [(abs:si kk) r]
      $
    ::  +dfloor: @rs -> @s, floor (round toward -infinity, NOT toward
    ::  zero) -- +poisson-ptrs's candidate k can be negative, where
    ::  truncate-toward-zero (what +toi under +mz gives) would round the
    ::  wrong way.  Positive/exact values: truncation IS the floor.
    ::  Negative non-integers: truncation rounds up, so subtract 1.
    ++  dfloor
      |=  x=@rs
      ^-  @s
      =/  t  (need (toi:mz x))
      ?:  |((gte:m x .0) (equ:m x (san:m t)))
        t
      (dif:si t --1)
    ::  +loggam: @rs -> @rs, natural log of the gamma function (so log(k!)
    ::  = loggam(k+1)).  Stirling's series with a small-argument recurrence
    ::  for x<7 (shifts x up by whole integers until >=7, where the series
    ::  is accurate, then subtracts the shifted-away log-factors back out)
    ::  -- transcribed from NumPy's random_loggam (src/distributions/
    ::  distributions.c), itself from Zhang & Jin's SPECFUN algorithm, not
    ::  derived independently: log-gamma accuracy is exactly the kind of
    ::  detail worth matching a vetted reference bit-for-bit rather than
    ::  reconstructing from the general shape of Stirling's series.
    ++  loggam
      |=  x=@rs
      ^-  @rs
      ?:  |((equ:m x .1) (equ:m x .2))
        .0
      =/  a  ^-  (list @rs)
        :~  .0.08333333333333333  .-0.002777777777777778  .0.0007936507936507937
            .-0.0005952380952380952  .0.0008417508417508417  .-0.001917526917526918
            .0.006410256410256411  .-0.02955065359477124  .0.1796443723688307
            .-1.39243221690590
        ==
      =/  lg2pi  .1.8378770664093453
      =/  n  ?:((lth:m x .7) (abs:si (need (toi:mz (sub:m .7 x)))) 0)
      =/  x0  (add:m x (sun:m n))
      =/  x2  (div:m .1 (mul:m x0 x0))
      =/  gl0
        =/  acc  (snag 9 a)
        =/  i  0
        |-  ^-  @rs
        ?:  =(i 9)
          acc
        $(acc (add:m (mul:m acc x2) (snag (sub 8 i) a)), i +(i))
      =/  gl
        %+  sub:m
          %+  add:m
            (add:m (div:m gl0 x0) (mul:m .0.5 lg2pi))
          (mul:m (sub:m x0 .0.5) (log:m x0))
        x0
      ?.  (lth:m x .7)
        gl
      =/  fin
        =/  gl2  gl
        =/  x02  x0
        =/  k  1
        |-  ^-  @rs
        ?:  (gth k n)
          gl2
        $(gl2 (sub:m gl2 (log:m (sub:m x02 .1))), x02 (sub:m x02 .1), k +(k))
      fin
    ::    +binomial:  [r=rng n=@ p=@rs] -> [@ud rng]
    ::
    ::  Binomial(n,p) via inversion by recurrence: P(X=0)=(1-p)^n, then
    ::  P(X=k+1) = P(X=k) * (n-k)/(k+1) * p/(1-p), accumulating the CDF
    ::  until it exceeds a single uniform draw.  Only valid for
    ::  n*min(p,1-p) < 30 (Kachitvichyanukul & Schmeiser 1988) -- above
    ::  that, inversion needs too many terms and BTPE is the right
    ::  algorithm, deferred to NEXT-STEPS; this arm crashes rather than
    ::  running slow (many terms) or, worse, silently truncating.
    ::    Source
    ++  binomial
      |=  [r=rng:rand n=@ p=@rs]
      ^-  [out=@ud r=rng:rand]
      ~|  %rand-bad-prob
      ?>  &((gte:m p .0) (lte:m p .1))
      =/  q     (sub:m .1 p)
      =/  qmin  ?:((lth:m p .0.5) p q)
      ~|  %rand-binomial-n-too-large
      ?>  (lth:m (mul:m (sun:m n) qmin) .30)
      =^  u  r  (rs:uni r)
      =/  f    (pow:m q (sun:m n))
      =/  cdf  f
      =/  k    0
      |-  ^-  [out=@ud r=rng:rand]
      ?:  (lth:m u cdf)
        [k r]
      =/  f2  %+  mul:m  f
        (mul:m (div:m (sun:m (sub n k)) (sun:m +(k))) (div:m p q))
      $(f f2, cdf (add:m cdf f2), k +(k))
    ::    +dirichlet:  [r=rng alphas=(list @rs)] -> [(list @rs) rng]
    ::
    ::  Dirichlet(alphas): draw one gamma(alpha_i) per entry, normalize by
    ::  their sum.  Crashes on an empty .alphas (and, via +gamma's own
    ::  precondition, on any non-positive alpha).
    ::    Source
    ++  dirichlet
      |=  [r=rng:rand alphas=(list @rs)]
      ^-  [out=(list @rs) r=rng:rand]
      ~|  %rand-empty-list
      ?>  !=(~ alphas)
      =/  gd
        |-  ^-  [(list @rs) rng:rand]
        ?~  alphas
          [~ r]
        =^  x     r  (gamma r i.alphas)
        =^  rest  r  $(alphas t.alphas)
        [[x rest] r]
      =/  total  (roll -.gd add:m)
      [(turn -.gd |=(g=@rs (div:m g total))) +.gd]
  --
--
--
