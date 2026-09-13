  ::  /sur/rand
::::  Types for /lib/rand and its adapters (i754rand, twocrand, fixedrand,
::::  complexrand, unumrand): the engine states and the tagged +$rng union.
::
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
+$  pcg64  [state=@ inc=@]
::    $rng:  a generic stream -- a tagged union so distribution code is
::  engine-agnostic.
+$  rng
  $%  [%phil p=phil]
      [%sm64 s=sm64]
      [%pcg p=pcg64]
  ==
--
