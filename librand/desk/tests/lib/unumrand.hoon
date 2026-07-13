  ::  /tests/lib/unumrand
::::
::    /lib/unumrand: +posit-lattice regressions at all five width doors (one
::    proving the identity-on-non-NaR-bits passthrough, one proving the
::    NaR-rejection redraw actually redraws and threads rng state correctly
::    through both draws), +posit-unit regressions at posit8/16/32 (two
::    seeds each), and one test proving the "sample at @rd via i754rand,
::    convert via /lib/unum's existing +from-rd" distribution composition
::    end to end.
::
::    The RIGOROUS chi-square verification against librand/tools/
::    posit_unit_check.py's exact oracle is a separate, offline check (see
::    NEXT-STEPS.md) -- not reproduced here, since embedding a 256- or
::    65536-row expected-probability table as Hoon literals isn't practical
::    for a self-contained regression test.
::
/+  *test,
    rand,
    i754rand,
    unum,
    unumrand
|%
::  +posit-lattice is the identity on non-NaR raw bits -- same seed, same
::  width, same result as a plain +bits:uni:rand draw, whenever that draw
::  doesn't happen to be NaR (as it isn't at seed 0, any width here).
++  test-posit-lattice-is-bits  ^-  tang
  ;:  weld
    %+  expect-eq
      !>((bits:uni:rand (from-atom:seed:rand %sm64 0) 8))
      !>((posit-lattice:rpb:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>((bits:uni:rand (from-atom:seed:rand %sm64 0) 16))
      !>((posit-lattice:rph:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>((bits:uni:rand (from-atom:seed:rand %sm64 0) 32))
      !>((posit-lattice:rps:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>((bits:uni:rand (from-atom:seed:rand %sm64 0) 64))
      !>((posit-lattice:rpd:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>((bits:uni:rand (from-atom:seed:rand %sm64 0) 128))
      !>((posit-lattice:rpq:unumrand (from-atom:seed:rand %sm64 0)))
  ==
::  Seed %sm64 149's first 8-bit draw is EXACTLY NaR (0x80) -- ship-found by
::  brute-force search, not hand-picked.  +posit-lattice:rpb must reject it,
::  redraw (0x9b, the second 8-bit draw from the post-first-draw state), and
::  return an rng state matching two draws consumed, not one.
++  test-posit-lattice-rejects-nar  ^-  tang
  %+  expect-eq
    !>([out=0x9b r=[%sm64 s=0x3c6e.f372.fe94.f8bf]])
    !>((posit-lattice:rpb:unumrand (from-atom:seed:rand %sm64 149)))
::  +posit-unit value regressions, two seeds per width (posit8/16/32 only --
::  posit64/128 are out of scope, see NEXT-STEPS.md).  Each value is exactly
::  encode(False, -k, u, n) for u = the corresponding k-bit +bits:uni:rand
::  draw, cross-checked against librand/tools/posit_unit_check.py's encode()
::  before being taken as ground truth here (not hand-derived).
++  test-posit-unit  ^-  tang
  ;:  weld
    %+  expect-eq
      !>(`@`0x37)
      !>(out:(posit-unit:rpb:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x39)
      !>(out:(posit-unit:rpb:unumrand (from-atom:seed:rand %sm64 1)))
    %+  expect-eq
      !>(`@`0x3e22)
      !>(out:(posit-unit:rph:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x3911)
      !>(out:(posit-unit:rph:unumrand (from-atom:seed:rand %sm64 1)))
    %+  expect-eq
      !>(`@`0x35cf.13cd)
      !>(out:(posit-unit:rps:unumrand (from-atom:seed:rand %sm64 0)))
    %+  expect-eq
      !>(`@`0x3bee.b8da)
      !>(out:(posit-unit:rps:unumrand (from-atom:seed:rand %sm64 1)))
  ==
::  The distribution composition rand-spec.md 12.4 describes but this
::  library doesn't wrap: draw a standard normal deviate at @rd via
::  /lib/i754rand, then convert to posit32 via /lib/unum's EXISTING
::  +from-rd (no new /lib/unum plumbing needed, unlike fixedrand).  Both
::  steps are independently tested elsewhere (i754rand's own +test-normal,
::  /lib/unum's own from-rd round-trip tests); this proves the COMPOSITION
::  against a known-good ship-verified result.
++  test-unumrand-dist-composition  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  z  r  (normal:rd:dist:i754rand r)
  %+  expect-eq
    !>(`@`0xc0d5.045b)
    !>((from-rd:rps:unum z))
--
