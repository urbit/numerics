::::  /lib/rand -- deterministic pseudo-random number generation
::
::  NOT CRYPTOGRAPHIC.  This library produces deterministic pseudo-randomness
::  for numerical work (Monte Carlo, ML init, shuffles, simulation).  It must
::  NOT be used for key material, nonces, or anything adversarial.  Entropy
::  acquisition is Arvo's job (`eny`); cryptographic randomness is Zuse's job.
::
::  Status: milestone 1 of rand-spec.md -- ++split-mix and ++seed only.
::  ++philox, ++pcg, ++uni, ++dist, and ++sample land in later milestones.
::
~%  %non  ..part  ~  :: jet registration; nest non in hex (cf /lib/math, /lib/twoc)
|%
::    $phil:  Philox key/counter state
::
::  ctr is the 128-bit counter as a single @, key is the 64-bit key as a
::  single @.  Both are plain atoms; width discipline is enforced by masking,
::  never by aura tricks.  (++philox itself lands in milestone 2.)
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
  ::  see the file-level warning.  `eny` is 512 bits; +fold-eny below
  ::  compresses it to the 64 bits +from-atom expects.
  ::    Source
  ++  from-eny
    |=  [eng=?(%phil %sm64 %pcg) eny=@uvJ]
    ^-  rng
    (from-atom eng (fold-eny eny))
  ::    +fold-eny:  @uvJ -> @
  ::
  ::  Compress a wide atom (512-bit .eny) to 64 bits: split into eight
  ::  64-bit words, low word first, and fold them through +mix in order
  ::  (acc starts at 0; acc = (mix acc word) per word).  Not itself an rng
  ::  seed -- feeds +from-atom.  The same chunk-and-fold rule applies to any
  ::  seed/salt wider than 64 bits (e.g. a wide +fork salt).
  ::  Source
  ++  fold-eny
    |=  eny=@uvJ
    ^-  @
    =/  i    0
    =/  acc  0
    |-  ^-  @
    ?:  =(i 8)
      acc
    $(acc (mix acc (cut 6 [i 1] eny)), i +(i))
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
--
