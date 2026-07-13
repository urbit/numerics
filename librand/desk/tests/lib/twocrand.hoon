  ::  /tests/lib/twocrand
::::
::    /lib/twocrand: +twoc-full (a trivial re-export of +bits:uni:rand,
::    checked for exact equality with it -- no independent value needed
::    since it IS that arm) and +twoc-between (value regression against a
::    by-hand Lemire trace, plus a crash test for a>b in twoc order).
::
/+  *test,
    rand,
    twocrand
|%
::  +twoc-full is exactly +bits:uni:rand -- same seed, same width, same
::  result and same resulting rng, by construction (not a coincidence to
::  regress against, but worth asserting so a future refactor that
::  breaks the passthrough is caught).
++  test-twoc-full-is-bits  ^-  tang
  %+  expect-eq
    !>((bits:uni:rand (from-atom:seed:rand %sm64 0) 8))
    !>((twoc-full:twocrand (from-atom:seed:rand %sm64 0) 8))
::  w=8, a=0xfb (-5), b=0x5 (5): span=11, below(11)=9 (the same known
::  Lemire trace as ++uni's own +test-between), biased back by +9 from
::  -5 lands on 4.
++  test-twoc-between  ^-  tang
  %+  expect-eq
    !>(`@`4)
    !>(out:(twoc-between:twocrand (from-atom:seed:rand %sm64 0) 8 0xfb 0x5))
::  a > b in twoc order (0x5 = 5, 0xfb = -5) crashes.
++  test-twoc-between-bad-range  ^-  tang
  %-  expect-fail
  |.((twoc-between:twocrand (from-atom:seed:rand %sm64 0) 8 0x5 0xfb))
--
