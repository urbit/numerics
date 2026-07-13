::  Tests for Saloon's +rand-ray layer (rand-spec.md section 8): +fill-
::  uniform/+fill-normal/+fill-expon/+fill-below in the +sa core.
::
::  Every value below is cross-checked against a DIRECT call to the
::  underlying /lib/rand or /lib/i754rand primitive at the expected
::  counter/state, not hand-derived -- proving the per-element counter
::  math (rand-spec.md section 8's %phil "ctr0+i" / "ctr0+i*2^32 window"
::  design), not just that the ray-filling loop runs without crashing.
::
/-  ls=lagoon
/+  *test, *saloon, *lagoon, rand, i754rand
|%
::  +fill-uniform, %phil: element i must equal a plain rs:uni:i754rand
::  draw at counter i (rand-spec.md: "element i uses counter ctr0+i"),
::  and the post-state counter must be exactly ctr0+n (here 0+3=3).
++  test-fill-uniform-phil  ^-  tang
  =/  res  (fill-uniform:sa [~[3] 5 %i754 ~] (from-atom:seed:rand %phil 0))
  ?>  ?=(%phil -.r.res)
  =/  key0  key.p.r.res
  ;:  weld
    %+  expect-eq  !>(`@`3)                       !>(ctr.p.r.res)
    %+  expect-eq
      !>(out:(rs:uni:i754rand [%phil p=[key0 0]]))
      !>((end [0 32] data.ray.res))
    %+  expect-eq
      !>(out:(rs:uni:i754rand [%phil p=[key0 1]]))
      !>((cut 0 [32 32] data.ray.res))
    %+  expect-eq
      !>(out:(rs:uni:i754rand [%phil p=[key0 2]]))
      !>((cut 0 [64 32] data.ray.res))
  ==
::  +fill-uniform, %sm64 (sequential): elements thread the rng state
::  exactly the way two direct, sequential rs:uni:i754rand calls would.
++  test-fill-uniform-sm64  ^-  tang
  =/  res  (fill-uniform:sa [~[2] 5 %i754 ~] (from-atom:seed:rand %sm64 0))
  =/  a0  (from-atom:seed:rand %sm64 0)
  =^  e0  a0  (rs:uni:i754rand a0)
  =^  e1  a0  (rs:uni:i754rand a0)
  ;:  weld
    %+  expect-eq  !>(e0)  !>((end [0 32] data.ray.res))
    %+  expect-eq  !>(e1)  !>((cut 0 [32 32] data.ray.res))
    %+  expect-eq  !>(a0)  !>(r.res)
  ==
::  +fill-normal, %phil: element i must equal a normal:rd:dist:i754rand
::  draw at counter i*2^32 (rand-spec.md's counter-WINDOW design for
::  rejection-based transforms), and the post-state must be exactly
::  ctr0 + n*2^32 (here 0 + 2*2^32), regardless of how many sub-draws
::  each element's rejection loop actually used.
++  test-fill-normal-phil  ^-  tang
  =/  res  (fill-normal:sa [~[2] 6 %i754 ~] (from-atom:seed:rand %phil 0))
  ?>  ?=(%phil -.r.res)
  =/  key0  key.p.r.res
  =/  win  (bex 32)
  ;:  weld
    %+  expect-eq  !>(`@`(mul 2 win))             !>(ctr.p.r.res)
    %+  expect-eq
      !>(out:(normal:rd:dist:i754rand [%phil p=[key0 0]]))
      !>((end [0 64] data.ray.res))
    %+  expect-eq
      !>(out:(normal:rd:dist:i754rand [%phil p=[key0 win]]))
      !>((cut 0 [64 64] data.ray.res))
  ==
::  +fill-expon, %phil: same window scheme as +fill-normal, with the
::  distribution's own lambda parameter threaded through.
++  test-fill-expon-phil  ^-  tang
  =/  res  (fill-expon:sa [~[2] 6 %i754 ~] (from-atom:seed:rand %phil 0) .~2)
  ?>  ?=(%phil -.r.res)
  =/  key0  key.p.r.res
  =/  win  (bex 32)
  ;:  weld
    %+  expect-eq  !>(`@`(mul 2 win))             !>(ctr.p.r.res)
    %+  expect-eq
      !>(out:(expon:rd:dist:i754rand [%phil p=[key0 0]] .~2))
      !>((end [0 64] data.ray.res))
  ==
::  +fill-below, %sm64: elements are plain sequential Lemire draws.
++  test-fill-below-sm64  ^-  tang
  =/  res  (fill-below:sa [~[1] 5 %uint ~] 100 (from-atom:seed:rand %sm64 0))
  %+  expect-eq
    !>(out:(below:uni:rand (from-atom:seed:rand %sm64 0) 100))
    !>((end [0 32] data.ray.res))
::  Domain violations crash, tagged (rand-spec.md section 9 policy).
++  test-fill-uniform-bad-kind-crashes  ^-  tang
  %-  expect-fail
  |.((fill-uniform:sa [~[2] 5 %uint ~] (from-atom:seed:rand %sm64 0)))
++  test-fill-uniform-bad-bloq-crashes  ^-  tang
  %-  expect-fail
  |.((fill-uniform:sa [~[2] 7 %i754 ~] (from-atom:seed:rand %sm64 0)))
++  test-fill-below-bad-kind-crashes  ^-  tang
  %-  expect-fail
  |.((fill-below:sa [~[2] 5 %i754 ~] 100 (from-atom:seed:rand %sm64 0)))
--
