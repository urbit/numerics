  ::  /tests/lib/fixedrand
::::
::    /lib/fixedrand: +fixed (a trivial re-export of +bits:uni:rand at the
::    precision's own width, checked for exact equality with it), +fixed-unit
::    and +fixed-between (value regressions against known-good ship output),
::    and one test proving the "sample at @rd via /lib/i754rand, quantize via
::    /lib/fixed's +from-rd" composition end to end (rand-spec.md section
::    12.1's distribution story for fixed-point, per the "document + test,
::    don't ship dedicated wrapper arms" decision -- see librand/NEXT-STEPS.md).
::
/+  *test,
    rand,
    twocrand,
    fixedrand,
    i754rand,
    fixed
^|
=/  q88=prec:fixed  [8 8]
|%
::  +fixed is exactly +bits:uni:rand at width (wid q88) -- same seed, same
::  result and resulting rng, by construction (not a coincidence to regress
::  against, but worth asserting so a future refactor that breaks the
::  passthrough is caught, same rationale as twocrand's own +twoc-full test).
++  test-fixed-is-bits  ^-  tang
  %+  expect-eq
    !>((bits:uni:rand (from-atom:seed:rand %sm64 0) (wid:fixed q88)))
    !>((fixed:fixedrand (from-atom:seed:rand %sm64 0) q88))
::  seed %sm64 0, q8.8: known-good ship output.
++  test-fixed-unit  ^-  tang
  %+  expect-eq
    !>(`@`175)
    !>(out:(fixed-unit:fixedrand (from-atom:seed:rand %sm64 0) q88))
::  range [-1.0, 3.0] in q8.8 (0x1.ff00, 0x300): known-good ship output,
::  delegates to +twoc-between at width 17 (rand-spec.md section 12.2).
++  test-fixed-between  ^-  tang
  %+  expect-eq
    !>(`@`649)
    !>(out:(fixed-between:fixedrand (from-atom:seed:rand %sm64 0) q88 0x1.ff00 0x300))
::  a > b in twoc order crashes, same as the +twoc-between it delegates to.
++  test-fixed-between-bad-range  ^-  tang
  %-  expect-fail
  |.((fixed-between:fixedrand (from-atom:seed:rand %sm64 0) q88 0x300 0x1.ff00))
::  The distribution composition rand-spec.md 12.1 describes but this library
::  doesn't wrap: draw a standard normal deviate at @rd via /lib/i754rand,
::  then quantize it to q8.8 via /lib/fixed's +from-rd.  Both steps are
::  independently tested elsewhere (i754rand's own +test-normal, fixed's own
::  +test-rd-bridge); this proves the COMPOSITION against a known-good
::  ship-verified result, not either step in isolation.
++  test-fixed-dist-composition  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  z  r  (normal:rd:dist:i754rand r)
  %+  expect-eq
    !>(`@`130.829)
    !>((from-rd:fixed z q88))
--
