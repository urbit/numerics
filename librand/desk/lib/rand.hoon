/-  rand
/+  math
=+  rand
::::  /lib/rand -- deterministic pseudo-random number generation (plumbing)
::
::  NOT CRYPTOGRAPHIC.  This library produces deterministic pseudo-randomness
::  for numerical work (Monte Carlo, ML init, shuffles, simulation).  It must
::  NOT be used for key material, nonces, or anything adversarial.  Entropy
::  acquisition is Arvo's job (`eny`); cryptographic randomness is Zuse's job.
::
::  This is the FOUNDATION every other rand library imports: the three
::  counter/sequential engines (++philox, ++split-mix, ++pcg), seeding
::  (++seed), the generic engine-dispatched draw/fork (+step, +fork), the
::  door facade (++gen), integer-only uniform deviates (++uni: +bits,
::  +below, +between -- no floats here), and generic sampling (++sample:
::  shuffle/permutation/choice/reservoir -- +alias moved out, see below).
::  +$rng and its component engine-state types live in /sur/rand, not
::  here, so any file that only needs the type (not this library's logic)
::  can get it via a lightweight `/-` import instead of pulling in this
::  whole file.
::
::  IEEE-754 float generation and the distributions built on it (++dist)
::  live in a separate "porcelain" library, /lib/i754rand -- matching how
::  /lib/twoc, /lib/fixed, /lib/complex, and /lib/unum are already kept
::  separate from /lib/math in this codebase.  ++sample's own +alias
::  (Vose's alias method) moved there too, alongside ++dist's +categorical
::  which uses it: +draw needs an actual uniform @rd float draw to decide
::  accept-vs-redirect, so it is NOT float-free the way the rest of
::  ++sample is, and belongs with the porcelain, not the plumbing.
::  Non-float output adapters (twocrand, fixedrand, complexrand,
::  unumrand) are i754rand's siblings, each importing this file for the
::  same engine/uni/sample primitives.
::
::  Status: milestones 1-6 of rand-spec.md, split into this file plus
::  i754rand per the above.  See rand-spec.md and each library's own
::  NEXT-STEPS.md for the full milestone order.
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
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
::
::::                    ++uni                          ::  (5) uniform deviates (integer)
::
::  +bits/+below/+between only -- the integer/bit-pattern primitives every
::  adapter (twocrand, fixedrand, unumrand) and i754rand's own float
::  generators are built on.  IEEE-754 float generation (+rs/+rd/+rh/+rq/
::  +rs-oo/+rd-oo) lives in i754rand's OWN ++uni instead of here, kept
::  separate so non-float adapters never need to import /lib/math.
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
  --
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
::::                    ++sample                        ::  (7) sampling
::
::  Shuffles, permutations, choice, reservoir sampling.  +below:uni is the
::  shared unbiased primitive throughout, per section 5.1's policy: never
::  modulo-bias.  Vose's alias method (section 6.1) moved to i754rand: its
::  +draw needs an actual uniform @rd float draw (to decide accept vs.
::  redirect), so it isn't float-free the way the rest of this core is.
::
++  sample
  |%
  ::    +shuffle:  [r=rng l=(list)] -> [(list) rng]
  ::
  ::  Fisher-Yates, iterating from the end down: for i from n-1 to 1, swap
  ::  position i with a uniform position in [0,i] (+below:uni, inclusive
  ::  via i+1).  Implemented over a (map @ud _elem) rather than repeated
  ::  +snag/+oust on the list itself, which would be O(n^2); the map
  ::  round-trip is O(n log n).  Wet gate so the element type is
  ::  preserved for the caller (a dry `(list)` would erase it to `*`).
  ::    Examples
  ::      > (shuffle:sample (from-atom:seed:rand %sm64 0) ~[1 2 3 4 5])
  ::      [~[3 4 1 2 5] [%sm64 s=0x78dd.e6e5.fd29.f054]]
  ::  Source
  ++  shuffle
    |*  [r=rng l=(list)]
    ^+  [l r]
    =/  n  (lent l)
    ?:  (lth n 2)
      [l r]
    =/  elem  ?>(?=(^ l) i.l)
    =/  m  (~(gas by *(map @ud _elem)) (turn (gulf 0 (dec n)) |=(k=@ud [k (snag k l)])))
    =/  rr  r
    =/  ix  (dec n)
    |-  ^+  [l r]
    ?:  =(ix 0)
      [(turn (gulf 0 (dec n)) |=(k=@ud (~(got by m) k))) rr]
    =^  j  rr  (below:uni rr +(ix))
    =/  vi  (~(got by m) ix)
    =/  vj  (~(got by m) j)
    %=  $
      m   (~(put by (~(put by m) ix vj)) j vi)
      ix  (dec ix)
    ==
  ::    +permutation:  [r=rng n=@] -> [(list @ud) rng]
  ::
  ::  A uniform random permutation of 0..n-1: +shuffle of (gulf 0 (dec n)).
  ::    Source
  ++  permutation
    |=  [r=rng n=@]
    ^-  [(list @ud) rng]
    ?:  =(n 0)
      [~ r]
    (shuffle r (gulf 0 (dec n)))
  ::    +choice:  [r=rng l=(list)] -> [* rng]
  ::
  ::  One uniformly-chosen element.  Wet gate (rand-spec.md section 7
  ::  writes this arm's signature dry, `[* rng]` -- but a bare `*` return
  ::  loses the element's type at every call site for no benefit, and a
  ::  dry `(list)` argument runs into a real Hoon type-inference wall at
  ::  `+snag` -- mull-grow/nest-fail trying to prove a `(list *)` is
  ::  non-null after the `?~` guard.  `|*` sidesteps both: each call site
  ::  gets its own precise element type, same reasoning as +shuffle.
  ::  Crashes on an empty list.
  ::    Source
  ++  choice
    |*  [r=rng l=(list)]
    ~|  %rand-empty-list
    ?>  !=(~ l)
    =^  i  r  (below:uni r (lent l))
    [(snag i l) r]
  ::    +choices:  [r=rng n=@ l=(list)] -> [(list) rng]
  ::
  ::  .n elements chosen uniformly WITH replacement (independent +choice
  ::  draws).  Crashes on an empty .l if n > 0.
  ::    Source
  ++  choices
    |*  [r=rng n=@ l=(list)]
    ^+  [l r]
    ?:  =(n 0)
      [~ r]
    ~|  %rand-empty-list
    =/  elem  ?>(?=(^ l) i.l)
    =/  i    0
    =/  acc  *(list _elem)
    =/  rr   r
    |-  ^+  [l r]
    ?:  =(i n)
      [(flop acc) rr]
    =^  x  rr  (choice rr l)
    %=  $
      i    +(i)
      acc  [x acc]
    ==
  ::    +sample-n:  [r=rng k=@ l=(list)] -> [(list) rng]
  ::
  ::  .k elements WITHOUT replacement, via partial Fisher-Yates (the first
  ::  -- here, for implementation convenience, the LAST -- k positions of
  ::  a full shuffle): iterate i from n-1 down to n-k, swapping position i
  ::  with a uniform position in [0,i], same as +shuffle but stopping
  ::  early: this is the spec'd "first n of a permutation," just taken
  ::  from whichever end +shuffle itself iterates from, and either end is
  ::  equally uniform. NOT repeated-rejection sampling. Crashes if k > the
  ::  list's length. Wet gate, same reasoning as +shuffle.
  ::    Source
  ++  sample-n
    |*  [r=rng k=@ l=(list)]
    ^+  [l r]
    =/  n  (lent l)
    ~|  %rand-bad-count
    ?>  (lte k n)
    ?:  =(k 0)
      [~ r]
    =/  elem  ?>(?=(^ l) i.l)
    =/  m  (~(gas by *(map @ud _elem)) (turn (gulf 0 (dec n)) |=(j=@ud [j (snag j l)])))
    =/  rr  r
    =/  ix  (dec n)
    =/  stop  (sub n k)
    |-  ^+  [l r]
    ?:  =(ix stop)
      [(turn (gulf stop (dec n)) |=(j=@ud (~(got by m) j))) rr]
    =^  j  rr  (below:uni rr +(ix))
    =/  vi  (~(got by m) ix)
    =/  vj  (~(got by m) j)
    %=  $
      m   (~(put by (~(put by m) ix vj)) j vi)
      ix  (dec ix)
    ==
  ::    +reservoir:  [r=rng n=@ l=(list)] -> [(list) rng]
  ::
  ::  Algorithm R (Vitter 1985): a uniform sample of .n items from .l,
  ::  processed one at a time (the algorithm this library's own +below is
  ::  built on doesn't need true streaming, but the same algorithm serves
  ::  agents that DO stream events one at a time).  The first .n items
  ::  seed the reservoir; each later item at index i replaces a uniformly
  ::  chosen reservoir slot with probability n/(i+1) (drawn as "is the
  ::  uniform draw in [0,i] less than n").  Crashes if n > the list's
  ::  length.  Wet gate, same reasoning as +shuffle.
  ::    Source
  ++  reservoir
    |*  [r=rng n=@ l=(list)]
    =/  len  (lent l)
    ~|  %rand-bad-count
    ?>  (lte n len)
    ^+  [(scag n l) r]
    ?:  =(n 0)
      [~ r]
    =/  elem  ?>(?=(^ l) i.l)
    =/  m  (~(gas by *(map @ud _elem)) (turn (gulf 0 (dec n)) |=(k=@ud [k (snag k l)])))
    =/  rr  r
    =/  ix  n
    |-  ^+  [(scag n l) r]
    ?:  =(ix len)
      [(turn (gulf 0 (dec n)) |=(k=@ud (~(got by m) k))) rr]
    =^  j  rr  (below:uni rr +(ix))
    ?:  (lth j n)
      %=  $
        m   (~(put by m) j (snag ix l))
        ix  +(ix)
      ==
    $(ix +(ix))
  --
--
