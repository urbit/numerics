::::  /lib/rand -- deterministic pseudo-random number generation
::
::  NOT CRYPTOGRAPHIC.  This library produces deterministic pseudo-randomness
::  for numerical work (Monte Carlo, ML init, shuffles, simulation).  It must
::  NOT be used for key material, nonces, or anything adversarial.  Entropy
::  acquisition is Arvo's job (`eny`); cryptographic randomness is Zuse's job.
::
::  Status: milestones 1-2 of rand-spec.md -- ++philox, ++split-mix, ++seed,
::  +step, +fork, ++gen.  ++pcg, ++uni, ++dist, and ++sample land later.
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
    ::  ++pcg (milestone 4) isn't implemented yet.  Unlike +fork's %pcg
    ::  branch below (which only remixes the state/inc tuple via +mix and
    ::  needs no engine-specific logic), a real draw needs the xsl-rr
    ::  output permutation -- so this branch crashes rather than silently
    ::  returning something wrong.
    ~|  %rand-pcg-step-not-yet-implemented
    !!
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
--
