/+  math
::::  /lib/rand -- deterministic pseudo-random number generation
::
::  NOT CRYPTOGRAPHIC.  This library produces deterministic pseudo-randomness
::  for numerical work (Monte Carlo, ML init, shuffles, simulation).  It must
::  NOT be used for key material, nonces, or anything adversarial.  Entropy
::  acquisition is Arvo's job (`eny`); cryptographic randomness is Zuse's job.
::
::  Status: milestones 1-5 of rand-spec.md -- ++philox, ++split-mix, ++seed,
::  +step, +fork, ++gen, ++uni, ++pcg, ++dist.  ++sample lands later.
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::    $phil:  Philox key/counter state
::
::  ctr is the 128-bit counter as a single @, key is the 64-bit key as a
::  single @.  Both are plain atoms; width discipline is enforced by masking,
::  never by aura tricks.
+$  phil  [key=@ ctr=@]
::    $sm64:  SplitMix64 state (64 bits)
+$  sm64  @
::    $pcg64:  PCG64 state (128-bit state, 128-bit odd increment)
::
::  (++pcg itself lands in milestone 4.)
+$  pcg64  [state=@ inc=@]
::    $rng:  a generic stream -- a tagged union so distribution code (later
::  milestones) is engine-agnostic.
+$  rng
  $%  [%phil p=phil]
      [%sm64 s=sm64]
      [%pcg p=pcg64]
  ==
::
::::                    ++philox                       ::  (3.1) Philox4x32-10
::
::  Reference: Salmon, Manuson, Jung, Shaw, "Parallel Random Numbers: As Easy
::  as 1, 2, 3" (SC'11); the Random123 library is the reference
::  implementation.  The primary generator: output is a pure function of
::  (key, counter), so a Lagoon jet can fill a ray in parallel and still
::  match this sequential Hoon loop exactly (rand-spec.md section 8).
::
++  philox
  |%
  ::    +block:  [key=@ ctr=@] -> @
  ::
  ::  The pure Philox4x32-10 function.  .ctr packs four 32-bit
  ::  little-endian lanes c0/c1/c2/c3 (c0 = bits [0,32), ... c3 = bits
  ::  [96,128)); .key packs two 32-bit lanes k0/k1 (k0 = bits [0,32), k1 =
  ::  bits [32,64)).  Ten rounds of the Philox S-box (mulhilo32 is exact:
  ::  a plain multiply, then hi/lo = high/low 32 bits of the product):
  ::
  ::    hi0,lo0 = mulhilo32(M0, c0) ;  hi1,lo1 = mulhilo32(M1, c2)
  ::    c0' = hi1 ^ c1 ^ k0 ;  c1' = lo1
  ::    c2' = hi0 ^ c3 ^ k1 ;  c3' = lo0
  ::    k0 += W0 ;  k1 += W1   (mod 2^32)
  ::
  ::  Output is the four final lanes repacked as one 128-bit atom in the
  ::  same little-endian lane order as .ctr (c0' lowest).  KAT vectors
  ::  (Random123's kat_vectors file -- all-zero, all-0xffffffff, and the
  ::  pi-digits vector -- cross-checked against an independent Python
  ::  re-implementation of this exact algorithm) are in tests/lib/rand.hoon.
  ::    Examples
  ::      > (block:philox 0 0)
  ::      0x9b00.dbd8.bc57.ac4c.e169.c58d.6627.e8d5
  ::  Source
  ++  block
    ~/  %block
    |=  [key=@ ctr=@]
    ^-  @
    =/  c0  (cut 0 [0 32] ctr)
    =/  c1  (cut 0 [32 32] ctr)
    =/  c2  (cut 0 [64 32] ctr)
    =/  c3  (cut 0 [96 32] ctr)
    =/  k0  (cut 0 [0 32] key)
    =/  k1  (cut 0 [32 32] key)
    =/  n   0
    |-  ^-  @
    ?:  =(n 10)
      :(add c0 (lsh [0 32] c1) (lsh [0 64] c2) (lsh [0 96] c3))
    =/  m0   (mul 0xd251.1f53 c0)
    =/  hi0  (rsh [0 32] m0)
    =/  lo0  (end [0 32] m0)
    =/  m1   (mul 0xcd9e.8d57 c2)
    =/  hi1  (rsh [0 32] m1)
    =/  lo1  (end [0 32] m1)
    %=  $
      c0  (mix hi1 (mix c1 k0))
      c1  lo1
      c2  (mix hi0 (mix c3 k1))
      c3  lo0
      k0  (end [0 32] (add k0 0x9e37.79b9))
      k1  (end [0 32] (add k1 0xbb67.ae85))
      n   +(n)
    ==
  ::    +next:  phil -> [out=@ p=phil]
  ::
  ::  One Philox4x32-10 draw: out = (block key ctr), new state advances the
  ::  counter by 1 (mod 2^128), key unchanged.  See +step (below) for how
  ::  callers reduce the 128-bit .out to a 64-bit draw.
  ::    Examples
  ::      > (next:philox [key=0 ctr=0])
  ::      [out=0x9b00.dbd8.bc57.ac4c.e169.c58d.6627.e8d5 p=[key=0 ctr=1]]
  ::  Source
  ++  next
    ~/  %next
    |=  p=phil
    ^-  [out=@ p=phil]
    [(block key.p ctr.p) [key.p (mod +(ctr.p) (bex 128))]]
  --
::
::::                    ++split-mix                    ::  (3.2) SplitMix64
::
::  Reference: Steele, Lea, Flood (OOPSLA'14); Vigna's public-domain C.  The
::  primary generator (++philox, milestone 2) is counter-based; SplitMix64
::  here serves two roles: a small sequential generator in its own right
::  (%sm64 in +$rng), and the seed-mixing primitive every engine's +seed
::  path and (later) +fork depend on -- see +finalize below.
::
++  split-mix
  |%
  ::    +finalize:  @ -> @
  ::
  ::  The z-mixing tail of SplitMix64's step function, exposed as its own
  ::  arm because +seed's +mix (below) reuses it on arbitrary atoms that are
  ::  NOT a sequential-generator state (no golden-ratio increment applies to
  ::  those inputs -- only +next's state-advance gets that step).  Maps
  ::  0 -> 0 (a fixed point of the finalizer alone, without the increment);
  ::  this is expected and is why +next always adds the increment before
  ::  finalizing, and why +mix's degenerate mix(0,0) = 0 is documented
  ::  rather than "fixed" (see +mix).
  ::    Examples
  ::      > (finalize:split-mix 0x9e37.79b9.7f4a.7c15)
  ::      0xe220.a839.7b1d.cdaf
  ::  Source
  ++  finalize
    |=  z=@
    ^-  @
    =.  z  (end [0 64] (mul (mix z (rsh [0 30] z)) 0xbf58.476d.1ce4.e5b9))
    =.  z  (end [0 64] (mul (mix z (rsh [0 27] z)) 0x94d0.49bb.1331.11eb))
    (mix z (rsh [0 31] z))
  ::    +next:  sm64 -> [out=@ s=sm64]
  ::
  ::  One SplitMix64 draw: advance state by the golden-ratio increment
  ::  (mod 2^64), then finalize.  Matches Vigna's reference C bit-for-bit;
  ::  see the KAT vectors in tests/lib/rand.hoon.
  ::    Examples
  ::      > (next:split-mix 0)
  ::      [out=0xe220.a839.7b1d.cdaf s=0x9e37.79b9.7f4a.7c15]
  ::  Source
  ++  next
    ~/  %next
    |=  s=sm64
    ^-  [out=@ s=sm64]
    =/  s2  (end [0 64] (add s 0x9e37.79b9.7f4a.7c15))
    [(finalize s2) s2]
  ::    +split:  sm64 -> [a=@ b=@ s=sm64]
  ::
  ::  Two decorrelated child states from one parent state, via two +next
  ::  draws used directly as child seeds.  Gamma-stepping per the original
  ::  paper is overkill for this library's needs and adds a second tunable
  ::  parameter for no benefit here; two draws is adequate and simpler.
  ::    Examples
  ::      > (split:split-mix 0)
  ::      [a=0xe220.a839.7b1d.cdaf b=0x6e78.9e6a.a1b9.65f4 s=0x3c6e.f372.fe94.f82a]
  ::  Source
  ++  split
    |=  s=sm64
    ^-  [a=@ b=@ s=sm64]
    =^  a  s  (next s)
    =^  b  s  (next s)
    [a b s]
  --
::
::::                    ++seed                          ::  (4) seeding
::
::  No implicit global seeding.  There is no "default rng"; callers own
::  their state.  NOT CRYPTOGRAPHIC -- see the file-level warning above;
::  +from-eny carries the same warning explicitly since it is the arm most
::  likely to be reached for by someone wanting "real" randomness.
::
++  seed
  |%
  ::    +from-atom:  [eng=?(%phil %sm64 %pcg) a=@] -> rng
  ::
  ::  Seed any engine deterministically from a single atom, treating .a as
  ::  raw SplitMix64 state and drawing as many outputs as that engine's
  ::  state needs.  %sm64 needs no pre-draw: SplitMix64's own state can be
  ::  seeded with any value (Vigna's reference explicitly allows this), and
  ::  +next's golden-ratio increment does the mixing on first use -- drawing
  ::  an extra output here would just be redundant hashing.  %phil takes the
  ::  first output as .key with .ctr reset to 0.  %pcg takes two outputs as
  ::  .state and .inc, forcing .inc odd (PCG64 requires an odd increment).
  ::    Examples
  ::      > (from-atom:seed %sm64 0)
  ::      [%sm64 s=0]
  ::      > (from-atom:seed %phil 0)
  ::      [%phil p=[key=0xe220.a839.7b1d.cdaf ctr=0]]
  ::      > (from-atom:seed %pcg 0)
  ::      [%pcg p=[state=0xe220.a839.7b1d.cdaf inc=0x6e78.9e6a.a1b9.65f5]]
  ::  Source
  ++  from-atom
    |=  [eng=?(%phil %sm64 %pcg) a=@]
    ^-  rng
    ?-  eng
      %sm64  [%sm64 s=(end [0 64] a)]
      %phil
        =^  key  a  (next:split-mix (end [0 64] a))
        [%phil p=[key=key ctr=0]]
      %pcg
        =^  st   a  (next:split-mix (end [0 64] a))
        =^  inc  a  (next:split-mix a)
        [%pcg p=[state=st inc=(con inc 1)]]
    ==
  ::    +from-eny:  [eng=?(%phil %sm64 %pcg) eny=@uvJ] -> rng
  ::
  ::  Same pipeline as +from-atom, over Arvo entropy.  NOT CRYPTOGRAPHIC --
  ::  see the file-level warning.  `eny` is 512 bits; +fold-wide below
  ::  compresses it to the 64 bits +from-atom expects.
  ::    Source
  ++  from-eny
    |=  [eng=?(%phil %sm64 %pcg) eny=@uvJ]
    ^-  rng
    (from-atom eng (fold-wide eny))
  ::    +fold-wide:  @ -> @
  ::
  ::  Compress an atom of ANY width to 64 bits: split into 64-bit words
  ::  (low word first, via +met bloq 6 for the word count) and fold them
  ::  through +mix in order (acc starts at 0; acc = (mix acc word) per
  ::  word).  Not itself an rng seed -- feeds +from-atom (for @uvJ eny,
  ::  always 8 words) and +fork (for salts wider than 64 bits, section
  ::  2.1c) alike, so there is exactly one fold algorithm in the library.
  ::    Examples
  ::      > (fold-wide:seed 0)
  ::      0
  ::  Source
  ++  fold-wide
    |=  a=@
    ^-  @
    =/  n    (met 6 a)
    =/  i    0
    =/  acc  0
    |-  ^-  @
    ?:  =(i n)
      acc
    $(acc (mix acc (cut 6 [i 1] a)), i +(i))
  ::    +mix:  [a=@ b=@] -> @
  ::
  ::  Combine two 64-bit seed words (e.g. a ship + a per-agent salt) into
  ::  one.  Fully specified as a two-word compression, NOT a single
  ::  finalizer call on the raw 128-bit concatenation of .a and .b -- the
  ::  finalizer is a function of one 64-bit word, so a wider input must be
  ::  folded, not fed in directly:
  ::
  ::    w0 = a mod 2^64  ;  w1 = b mod 2^64
  ::    return finalize(finalize(w0) XOR w1)
  ::
  ::  Order-sensitive by construction: (mix a b) differs from (mix b a) in
  ::  general, since .a is absorbed through one extra +finalize application
  ::  relative to .b.  This is the primitive +fork (milestone 2) uses to
  ::  make (fork (fork r a) b) differ from (fork (fork r b) a): each nested
  ::  +fork call adds another +finalize layer at a different depth.
  ::  Deterministic, no rng threading (pure @ -> @).
  ::
  ::  Degenerate case, documented rather than "fixed": mix(0,0) = 0, because
  ::  +finalize maps 0 -> 0 when no golden-ratio increment has been mixed in
  ::  first (see +finalize).  Real seeds and salts are essentially never
  ::  both exactly zero, so this is not a practical concern.
  ::
  ::  NAMING FOOTGUN: this arm's name collides with the Hoon standard
  ::  library's `++mix` (bitwise XOR).  Inside `++seed`, an unqualified
  ::  `(mix a b)` call resolves to THIS arm, not stdlib XOR -- which is why
  ::  the body below reaches for `^mix` to get bitwise XOR.  Any other code
  ::  added to `++seed` that wants stdlib XOR must do the same; callers
  ::  outside `++seed` are unaffected.
  ::    Examples
  ::      > (mix:seed 1 2)
  ::      0xef30.b01c.2974.aeeb
  ::      > (mix:seed 2 1)
  ::      0x3ec2.d42f.3a45.cc6e
  ::  Source
  ++  mix
    |=  [a=@ b=@]
    ^-  @
    =/  w0  (end [0 64] a)
    =/  w1  (end [0 64] b)
    (finalize:split-mix (^mix (finalize:split-mix w0) w1))
  --
::
::::                    +step                          ::  (2) generic draw
::
::  One 64-bit draw, engine-dispatched.  Lives at the top level (not nested
::  under any one engine core) because it dispatches across the +$rng
::  tagged union -- this is what lets ++uni/++dist/++sample (later
::  milestones) be written once each, not once per engine.
::
++  step
  |=  r=rng
  ^-  [out=@ r=rng]
  ?-  -.r
    %sm64
    =^  out  s.r  (next:split-mix s.r)
    [out r]
  ::
    %phil
    ::  +next:philox returns the full 128-bit block; +step keeps bits
    ::  [0,64) -- lanes c0'/c1', the low two lanes of the little-endian
    ::  repack in +block:philox.  The top 64 bits (c2'/c3') are discarded;
    ::  wasting them is acceptable at v1 (rand-spec.md section 2), and a
    ::  buffered variant that reuses both halves is a NEXT-STEPS item.
    =^  blk  p.r  (next:philox p.r)
    [(end [0 64] blk) r]
  ::
    %pcg
    =^  out  p.r  (next:pcg p.r)
    [out r]
  ==
::
::::                    +fork                          ::  (2.1c) key derivation
::
::  Deterministic, independent child stream from a parent stream and a
::  salt -- the JAX-style key-derivation payoff of counter-based-first
::  (rand-spec.md section 2.1c).  Callers fork by path instead of
::  threading one sequential state through a computation tree: draws are
::  independent of sibling evaluation order, and a draw inserted in one
::  branch doesn't perturb any other branch.
::    Examples
::      > =/  a  (fork (from-atom:seed %phil 0) 1)
::      > =/  b  (fork (from-atom:seed %phil 0) 2)
::      > =(a b)
::      %.n
::  Source
++  fork
  |=  [r=rng salt=@]
  ^-  rng
  ::  Salts up to 64 bits feed +mix directly (the common case: small
  ::  integer salts, e.g. per-event or per-layer indices).  Wider salts
  ::  (e.g. a caller-combined `(cat 6 a b)`, or an @uvJ) are folded to 64
  ::  bits first via +fold-wide, so no salt bits are silently dropped --
  ::  +mix itself only ever consumes one 64-bit word per side.
  =/  sw  ?:((lte (met 6 salt) 1) salt (fold-wide:seed salt))
  ?-  -.r
    %phil  r(key.p (mix:seed key.p.r sw), ctr.p 0)
  ::
    ::  "child increment derived the same way [as %phil's key], forced
    ::  odd; state re-mixed" (section 2.1c) -- read here as: both fields
    ::  pushed through +mix with the same salt-derived word.
    %pcg
    r(state.p (mix:seed state.p.r sw), inc.p (con (mix:seed inc.p.r sw) 1))
  ::
    ::  "the SplitMix split construction" (section 2.1c) is read here as
    ::  reusing the same +mix primitive +split's decorrelation relies on,
    ::  salt-directed rather than sequential -- NOT a literal call to
    ::  +split (which takes no salt and can't be path-directed).  This is
    ::  this implementation's resolution of that spec ambiguity.
    %sm64  r(s (mix:seed s.r sw))
  ==
::
::::                    ++gen                          ::  (2.1b) door facade
::
::  Ergonomic sugar over the functional core, never a storage format.
::  HARD RULE (section 2.1b): persisting this door in agent state pins a
::  battery across library upgrades (the classic +og misuse) -- persist
::  the +$rng noun, reconstruct the door locally.  The facade adds no
::  capability; it is a strict wrapper over +step/+fork, so jets register
::  on those functional arms only.
::
++  gen
  |_  r=rng
  ::    +draw:  gen -> [@ _..draw]
  ::
  ::  Door wrapper over +step.  `..draw` (not `+>`/`+>.$`) is the reliable
  ::  way to reference "this door, before any local rebinding" from an
  ::  arm's body regardless of whether the arm itself has a `|=` sample
  ::  (`+>`/`+>.$` axis arithmetic differs depending on that, which is a
  ::  footgun in itself -- `..<arm-name>` sidesteps it).
  ::    Examples
  ::      > =/  g  ~(. gen (from-atom:seed:rand %sm64 0))
  ::      > =^  x  g  draw:g
  ::      > x
  ::      16.294.208.416.658.607.535
  ::  Source
  ++  draw
    ^-  [@ _..draw]
    =^  out  r  (step r)
    [out ..draw(r r)]
  ::    +fork:  [gen @] -> _..fork
  ::
  ::  Door wrapper over the functional +fork.  NAMING FOOTGUN (same class
  ::  as ++seed's +mix, above): this arm's name shadows the top-level
  ::  +fork gate, so calling it unqualified from inside this arm's own
  ::  body would recurse into itself with the wrong arity -- the body
  ::  below reaches for `^fork` to reach the top-level gate.
  ::  Source
  ++  fork
    |=  salt=@
    ^-  _..fork
    ..fork(r (^fork r salt))
  --
::
::::                    ++uni                          ::  (5) uniform deviates
::
::  All arms `[r=rng ...] -> [out new-rng]` unless noted, per section 5.
::
::  NAMING FOOTGUN (same class as ++seed's +mix and ++gen's +fork): this
::  core defines arms named +rs/+rd/+rh/+rq, which shadow the Hoon standard
::  library's own +rs/+rd/+rh/+rq IEEE-754 doors.  /lib/math has the exact
::  same collision (its own +rs/+rd/+rh/+rq doors wrap the stdlib ones) and
::  handles it the same way: any reference to the real stdlib door from
::  inside ++uni uses ^rs/^rd/^rh/^rq, e.g. `(~(mul ^rs %n) a b)`.
::
++  uni
  |%
  ::    +bits:  [r=rng n=@] -> [@ rng]
  ::
  ::  .n uniform bits, assembled from ceil(n/64) +step draws, little-endian
  ::  (the first draw is bits [0,64) of the result, the second is bits
  ::  [64,128), and so on); the LAST draw is masked down to n's remainder
  ::  mod 64 when that remainder is nonzero (i.e. when n isn't itself a
  ::  multiple of 64).
  ::    Examples
  ::      > (bits:uni:rand (from-atom:seed:rand %sm64 0) 4)
  ::      [out=15 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
  ::  Source
  ++  bits
    |=  [r=rng n=@]
    ^-  [out=@ r=rng]
    =/  full   (div n 64)
    =/  rem    (mod n 64)
    =/  words  ?:(=(rem 0) full +(full))
    =/  i      0
    =/  acc    0
    |-  ^-  [out=@ r=rng]
    ?:  =(i words)
      [acc r]
    =^  out  r  (step r)
    =/  bc  ?:(&(=(i (dec words)) !=(rem 0)) rem 64)
    %=  $
      acc  (add acc (lsh [0 (mul i 64)] (end [0 bc] out)))
      i    +(i)
    ==
  ::    +below:  [r=rng n=@] -> [@ rng]
  ::
  ::  Uniform in [0,n), UNBIASED, via Lemire's multiply-shift with
  ::  rejection (Lemire 2019, "Fast Random Integer Generation in an
  ::  Interval").  .t depends only on .n, so it's computed once outside the
  ::  rejection loop, not per draw.  The primitive everything else
  ::  (shuffle, choice, +between) uses -- never ship modulo-bias
  ::  `(mod x n)` instead, here or anywhere else in this library.
  ::  Crashes (`?>`) on n=0.
  ::    Examples
  ::      > (below:uni:rand (from-atom:seed:rand %sm64 0) 10)
  ::      [out=8 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
  ::  Source
  ++  below
    ~/  %below
    |=  [r=rng n=@]
    ^-  [out=@ r=rng]
    ~|  %rand-zero-modulus
    ?>  !=(n 0)
    =/  t  (mod (sub (bex 64) n) n)
    |-  ^-  [out=@ r=rng]
    =^  x  r  (step r)
    =/  m  (mul x n)
    =/  l  (end [0 64] m)
    ?:  &((lth l n) (lth l t))
      $
    [(rsh [0 64] m) r]
  ::    +between:  [r=rng a=@s b=@s] -> [@s rng]
  ::
  ::  Inclusive signed range [a,b], via +below on the span (b - a + 1),
  ::  offset back by .a.  Crashes (`?>`) if a > b.
  ::    Examples
  ::      > (between:uni:rand (from-atom:seed:rand %sm64 0) -5 --5)
  ::      [out=--4 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
  ::  Source
  ++  between
    |=  [r=rng a=@s b=@s]
    ^-  [out=@s r=rng]
    ~|  %rand-bad-range
    ?>  !=(--1 (cmp:si a b))
    =/  span  +((abs:si (dif:si b a)))
    =^  off  r  (below r span)
    [(sum:si a (sun:si off)) r]
  ::    +rs:  rng -> [@rs rng]
  ::
  ::  Uniform in [0,1): draw 24 bits, multiply by the exact constant
  ::  2^-24.  +sun of a 24-bit unsigned integer and multiplying by an
  ::  exact power of two are BOTH exact IEEE-754 operations (no rounding,
  ::  regardless of mode), so every representable output is a multiple of
  ::  2^-24 -- exactly uniform over that lattice.  0 is possible; 1 is not.
  ::    Examples
  ::      > (rs:uni:rand (from-atom:seed:rand %sm64 0))
  ::      [out=.0.1164197 r=[%sm64 s=0x9e37.79b9.7f4a.7c15]]
  ::  Source
  ++  rs
    |=  r=rng
    ^-  [out=@rs r=rng]
    =^  b  r  (bits r 24)
    [(~(mul ^rs %n) (~(sun ^rs %n) b) `@rs`0x3380.0000) r]
  ::    +rd:  rng -> [@rd rng]
  ::
  ::  Same as +rs at double precision: 53 bits, times the exact constant
  ::  2^-53.
  ::    Source
  ++  rd
    |=  r=rng
    ^-  [out=@rd r=rng]
    =^  b  r  (bits r 53)
    [(~(mul ^rd %n) (~(sun ^rd %n) b) `@rd`0x3ca0.0000.0000.0000) r]
  ::    +rh:  rng -> [@rh rng]
  ::
  ::  Same as +rs at half precision: 11 bits, times the exact constant
  ::  2^-11.
  ::    Source
  ++  rh
    |=  r=rng
    ^-  [out=@rh r=rng]
    =^  b  r  (bits r 11)
    [(~(mul ^rh %n) (~(sun ^rh %n) b) `@rh`0x1000) r]
  ::    +rq:  rng -> [@rq rng]
  ::
  ::  Same as +rs at quad precision: 113 bits, times the exact constant
  ::  2^-113.
  ::    Source
  ++  rq
    |=  r=rng
    ^-  [out=@rq r=rng]
    =^  b  r  (bits r 113)
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
    |=  r=rng
    ^-  [out=@rs r=rng]
    =^  b  r  (bits r 24)
    [(~(mul ^rs %n) (~(sun ^rs %n) +((mul 2 b))) `@rs`0x3300.0000) r]
  ::    +rd-oo:  rng -> [@rd rng]
  ::
  ::  Same construction as +rs-oo at double precision: 53 bits, (2*bits+1)
  ::  * 2^-54.
  ::    Source
  ++  rd-oo
    |=  r=rng
    ^-  [out=@rd r=rng]
    =^  b  r  (bits r 53)
    [(~(mul ^rd %n) (~(sun ^rd %n) +((mul 2 b))) `@rd`0x3c90.0000.0000.0000) r]
  --
::
::::                    ++pcg                          ::  (3.3) PCG64 XSL-RR
::
::  Reference: O'Neill 2014, pcg-random.org; the pcg-c library
::  (`include/pcg_variants.h`) is the reference implementation.
::
::  **Spec correction, verified against pcg-c's `pcg_setseq_128_xsl_rr_64_
::  random_r` and independently against NumPy's vendored `pcg64.orig.h`
::  (both consistent): the real reference ADVANCES the state first, THEN
::  computes the output permutation from the NEW state** -- not "output
::  the current state, then advance" as an earlier draft of rand-spec.md
::  claimed. That earlier claim was backwards; see rand-spec.md section
::  3.3 for the correction. This is the implementation of the corrected
::  order: +next below steps first, then outputs.
::
++  pcg
  |%
  ::    +step-lcg:  pcg64 -> pcg64
  ::
  ::  The bare LCG advance, state' = state * MULT + inc (mod 2^128), no
  ::  output.  MULT is PCG's fixed 128-bit default multiplier; .inc is the
  ::  caller's odd stream increment (enforced odd at seed/fork time, not
  ::  here).
  ::  Source
  ++  step-lcg
    |=  p=pcg64
    ^-  pcg64
    %=  p
      state  (end [0 128] (add (mul state.p 0x2360.ed05.1fc6.5da4.4385.df64.9fcc.f645) inc.p))
    ==
  ::    +rotr64:  [value=@ rot=@] -> @
  ::
  ::  64-bit right-rotate: (value >> rot) | (value << (64-rot)), matching
  ::  pcg_rotr_64.  .rot is always < 64 in practice (it's the top 6 bits of
  ::  a 128-bit state, so already in [0,63]); the formula is correct at
  ::  rot=0 too (the left-shift-by-64 term vanishes under the 64-bit mask,
  ::  leaving value unchanged) so no special case is needed.
  ::  Source
  ++  rotr64
    |=  [value=@ rot=@]
    ^-  @
    (mix (rsh [0 rot] value) (end [0 64] (lsh [0 (sub 64 rot)] value)))
  ::    +output:  @ -> @
  ::
  ::  The xsl-rr output permutation of a 128-bit state to a 64-bit output:
  ::  rot = state >> 122; xored = (state >> 64) ^ (state & mask64); output
  ::  = rotr64(xored, rot).  Matches pcg_output_xsl_rr_128_64.
  ::  Source
  ++  output
    |=  state=@
    ^-  @
    =/  hi   (rsh [0 64] state)
    =/  lo   (end [0 64] state)
    =/  rot  (rsh [0 122] state)
    (rotr64 (mix hi lo) rot)
  ::    +next:  pcg64 -> [out=@ p=pcg64]
  ::
  ::  One PCG64 draw.  Per the spec correction above: advance first
  ::  (+step-lcg), THEN compute the output from the resulting NEW state --
  ::  matching `pcg_setseq_128_xsl_rr_64_random_r`'s actual order, not the
  ::  order an earlier draft of this library's spec claimed.
  ::    Examples
  ::      > (next:pcg [state=0xde2b.ce05.be01.3be3.d3f6.c45a.41e5.4320 inc=0x6d])
  ::      [out=0x86b1.da1d.7206.2b68 p=[state=0x10af.065f.4ea9.6e85.7bb2.a788.6ecb.d80d inc=0x6d]]
  ::  Source
  ++  next
    ~/  %next
    |=  p=pcg64
    ^-  [out=@ p=pcg64]
    =/  p2  (step-lcg p)
    [(output state.p2) p2]
  ::    +advance:  [p=pcg64 delta=@] -> pcg64
  ::
  ::  Advance .p by .delta LCG steps in O(log delta) via the standard PCG
  ::  skip-ahead trick: the affine map step(s) = s*mult + inc composes
  ::  under repeated squaring, so .delta steps can be applied as one
  ::  affine map (acc-mult, acc-plus) built by binary exponentiation,
  ::  rather than .delta sequential steps.  +jump is the .delta = 2^64
  ::  case this library actually needs (for stream partitioning); +advance
  ::  is exposed separately so the algorithm can be tested at a small,
  ::  tractable .delta instead of only at the astronomical one +jump uses:
  ::  `(advance p 3)` must equal three sequential `next:pcg` steps (see
  ::  tests/lib/rand.hoon).
  ::  Source
  ++  advance
    |=  [p=pcg64 delta=@]
    ^-  pcg64
    =/  cur-mult  0x2360.ed05.1fc6.5da4.4385.df64.9fcc.f645
    =/  cur-plus  inc.p
    =/  acc-mult  1
    =/  acc-plus  0
    |-  ^-  pcg64
    ?:  =(delta 0)
      %=  p
        state  (end [0 128] (add (mul acc-mult state.p) acc-plus))
      ==
    ?:  =(1 (dis delta 1))
      %=  $
        acc-mult  (end [0 128] (mul acc-mult cur-mult))
        acc-plus  (end [0 128] (add (mul acc-plus cur-mult) cur-plus))
        cur-plus  (end [0 128] (mul (add cur-mult 1) cur-plus))
        cur-mult  (end [0 128] (mul cur-mult cur-mult))
        delta     (rsh [0 1] delta)
      ==
    %=  $
      cur-plus  (end [0 128] (mul (add cur-mult 1) cur-plus))
      cur-mult  (end [0 128] (mul cur-mult cur-mult))
      delta     (rsh [0 1] delta)
    ==
  ::    +jump:  pcg64 -> pcg64
  ::
  ::  Advance by exactly 2^64 steps, for stream partitioning (give
  ::  different logical streams non-overlapping windows of the same
  ::  underlying sequence).  A thin +advance call at the one .delta this
  ::  library needs.
  ::  Source
  ++  jump
    |=  p=pcg64
    ^-  pcg64
    (advance p (bex 64))
  --
::
::::                    ++dist                          ::  (6) distributions
::
::  All arms at @rd only (the @rs routing and the poisson/binomial/
::  categorical/dirichlet arms land in a later milestone, per rand-spec.md's
::  milestone order).  Internally forced to round-to-nearest (%n) via +m
::  below, matching /lib/math's own doors -- callers never see or choose a
::  rounding mode here.  chi2/student-t aren't named in the milestone list
::  but are one-line compositions of gamma/normal with no new machinery, so
::  they land here too rather than waiting on an unscheduled slot.
::
::  Crashes (`?>` with `~|` tags) on out-of-domain parameters, per section
::  9's uniform error policy: domain violations are programmer error, not
::  data, so no NaN-returning "soft" errors.
::
++  dist
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
  ::      > (normal:dist:rand (from-atom:seed:rand %sm64 0))
  ::      [out=.~-0.9479938949723624 r=[%sm64 s=0x78dd.e6e5.fd29.f054]]
  ::  Source
  ++  normal
    |=  r=rng
    ^-  [out=@rd r=rng]
    |-  ^-  [out=@rd r=rng]
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
    |=  [r=rng mu=@rd sigma=@rd]
    ^-  [out=@rd r=rng]
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
    |=  [r=rng lambda=@rd]
    ^-  [out=@rd r=rng]
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
    |=  [r=rng alpha=@rd]
    ^-  [out=@rd r=rng]
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
    |=  [r=rng alpha=@rd]
    ^-  [out=@rd r=rng]
    =/  d  (sub:m alpha (div:m .~1 .~3))
    =/  c  (div:m .~1 (sqt:m (mul:m .~9 d)))
    |-  ^-  [out=@rd r=rng]
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
    |=  [r=rng a=@rd b=@rd]
    ^-  [out=@rd r=rng]
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
    |=  [r=rng k=@]
    ^-  [out=@rd r=rng]
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
    |=  [r=rng k=@]
    ^-  [out=@rd r=rng]
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
    |=  [r=rng p=@rd]
    ^-  [out=? r=rng]
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
    |=  [r=rng p=@rd]
    ^-  [out=@ud r=rng]
    ~|  %rand-bad-prob
    ?>  &((gth:m p .~0) (lte:m p .~1))
    ?:  (equ:m p .~1)
      [1 r]
    =^  u  r  (rd-oo:uni r)
    =/  raw  (div:m (log:m u) (log:m (sub:m .~1 p)))
    =/  fl   (abs:si (need (toi:mz raw)))
    [?:(=(raw (sun:m fl)) fl +(fl)) r]
  --
--
