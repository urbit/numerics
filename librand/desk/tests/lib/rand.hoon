  ::  /tests/lib/rand
::::
::    Milestone 1 (rand-spec.md): ++split-mix + ++seed.  The SplitMix64 KAT
::    vectors are checked against Vigna's reference algorithm (verified via
::    an independent Python re-implementation, not transcribed by hand);
::    +from-atom, +fold-eny, and +mix have no external reference (they are
::    this library's own seeding design), so those are checked by direct
::    computation of the spec'd formula and by determinism/order-sensitivity
::    properties instead of KAT.
::
/+  *test,
    rand
|%
::  SplitMix64 KAT: first 5 outputs from seed 0, against Vigna's reference C.
++  test-splitmix-kat-seed-0  ^-  tang
  =/  s0=@  0
  =^  o0  s0  (next:split-mix:rand s0)
  =^  o1  s0  (next:split-mix:rand s0)
  =^  o2  s0  (next:split-mix:rand s0)
  =^  o3  s0  (next:split-mix:rand s0)
  =^  o4  s0  (next:split-mix:rand s0)
  ;:  weld
    %+  expect-eq  !>(`@`0xe220.a839.7b1d.cdaf)  !>(o0)
    %+  expect-eq  !>(`@`0x6e78.9e6a.a1b9.65f4)  !>(o1)
    %+  expect-eq  !>(`@`0x6c4.5d18.8009.454f)   !>(o2)
    %+  expect-eq  !>(`@`0xf88b.b8a8.724c.81ec)  !>(o3)
    %+  expect-eq  !>(`@`0x1b39.896a.51a8.749b)  !>(o4)
  ==
::  SplitMix64 KAT: first 5 outputs from seed 0xdeadbeef.
++  test-splitmix-kat-seed-deadbeef  ^-  tang
  =/  s0=@  0xdead.beef
  =^  o0  s0  (next:split-mix:rand s0)
  =^  o1  s0  (next:split-mix:rand s0)
  =^  o2  s0  (next:split-mix:rand s0)
  =^  o3  s0  (next:split-mix:rand s0)
  =^  o4  s0  (next:split-mix:rand s0)
  ;:  weld
    %+  expect-eq  !>(`@`0x4adf.b90f.68c9.eb9b)  !>(o0)
    %+  expect-eq  !>(`@`0xde58.6a31.41a1.0922)  !>(o1)
    %+  expect-eq  !>(`@`0x21f.bc2f.8e1c.fc1d)   !>(o2)
    %+  expect-eq  !>(`@`0x7466.ce73.7be1.6790)  !>(o3)
    %+  expect-eq  !>(`@`0x3bfa.8764.f685.bd1c)  !>(o4)
  ==
::  +split is two +next draws used directly as child seeds.
++  test-split  ^-  tang
  =/  r  (split:split-mix:rand 0)
  ;:  weld
    %+  expect-eq  !>(`@`0xe220.a839.7b1d.cdaf)  !>(a.r)
    %+  expect-eq  !>(`@`0x6e78.9e6a.a1b9.65f4)  !>(b.r)
    %+  expect-eq  !>(`@`0x3c6e.f372.fe94.f82a)  !>(s.r)
  ==
::  +from-atom %sm64: state = seed directly, no pre-draw.
++  test-from-atom-sm64  ^-  tang
  %+  expect-eq
    !>(`rng:rand`[%sm64 s=0])
    !>((from-atom:seed:rand %sm64 0))
::  +from-atom %phil: key = first SplitMix output, ctr reset to 0.
++  test-from-atom-phil  ^-  tang
  %+  expect-eq
    !>(`rng:rand`[%phil p=[key=0xe220.a839.7b1d.cdaf ctr=0]])
    !>((from-atom:seed:rand %phil 0))
::  +from-atom %pcg: state/inc from two outputs, inc forced odd.
++  test-from-atom-pcg  ^-  tang
  %+  expect-eq
    !>(`rng:rand`[%pcg p=[state=0xe220.a839.7b1d.cdaf inc=0x6e78.9e6a.a1b9.65f5]])
    !>((from-atom:seed:rand %pcg 0))
::  same seed -> identical rng noun, for all three engines (determinism).
++  test-from-atom-determinism  ^-  tang
  ;:  weld
    %+  expect-eq  !>((from-atom:seed:rand %sm64 42))  !>((from-atom:seed:rand %sm64 42))
    %+  expect-eq  !>((from-atom:seed:rand %phil 42))  !>((from-atom:seed:rand %phil 42))
    %+  expect-eq  !>((from-atom:seed:rand %pcg 42))   !>((from-atom:seed:rand %pcg 42))
  ==
::  +mix: known values from the spec'd two-word compression, plus the
::  documented mix(0,0) = 0 degenerate point (see +mix's doc comment).
++  test-mix-values  ^-  tang
  ;:  weld
    %+  expect-eq  !>(`@`0xef30.b01c.2974.aeeb)  !>((mix:seed:rand 1 2))
    %+  expect-eq  !>(`@`0x3ec2.d42f.3a45.cc6e)  !>((mix:seed:rand 2 1))
    %+  expect-eq  !>(`@`0x0)                    !>((mix:seed:rand 0 0))
  ==
::  order-sensitivity: (mix a b) != (mix b a) in general -- the whole point
::  of the two-word compression (rand-spec.md section 4), load-bearing for
::  +fork's path-sensitivity property once +fork lands in milestone 2.
++  test-mix-order-sensitive  ^-  tang
  %+  expect-eq  !>(%.n)  !>(=((mix:seed:rand 42 7) (mix:seed:rand 7 42)))
::  +from-eny: determinism, and a +fold-eny KAT against the same reference
::  fold computed independently in Python over 8 known 64-bit words.
++  test-from-eny  ^-  tang
  =/  e=@uvJ
    `@uvJ`0x8888.8888.8888.8888.7777.7777.7777.7777.6666.6666.6666.6666.5555.5555.5555.5555.4444.4444.4444.4444.3333.3333.3333.3333.2222.2222.2222.2222.1111.1111.1111.1111
  ;:  weld
    %+  expect-eq  !>((from-eny:seed:rand %sm64 e))  !>((from-eny:seed:rand %sm64 e))
    %+  expect-eq  !>(`@`0xdeb.8a1b.ec35.d57d)  !>((fold-eny:seed:rand e))
  ==
--
