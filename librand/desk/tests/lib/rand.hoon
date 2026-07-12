  ::  /tests/lib/rand
::::
::    Milestones 1-2 (rand-spec.md): ++split-mix, ++philox, ++seed, +fork.
::    The SplitMix64 and Philox4x32-10 KAT vectors are checked against their
::    published reference values (Vigna's reference C for SplitMix64; the
::    Random123 kat_vectors file -- all-zero, all-0xffffffff, and the
::    pi-digits vector -- for Philox), each cross-checked against an
::    independent Python re-implementation before being transcribed here.
::    +from-atom, +fold-wide, +mix, and +fork have no external reference
::    (they are this library's own seeding/forking design), so those are
::    checked by direct computation of the spec'd formula and by
::    determinism/order-sensitivity/path-sensitivity properties instead.
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
::  +from-eny: determinism, and a +fold-wide KAT against the same reference
::  fold computed independently in Python over 8 known 64-bit words.
++  test-from-eny  ^-  tang
  =/  e=@uvJ
    `@uvJ`0x8888.8888.8888.8888.7777.7777.7777.7777.6666.6666.6666.6666.5555.5555.5555.5555.4444.4444.4444.4444.3333.3333.3333.3333.2222.2222.2222.2222.1111.1111.1111.1111
  ;:  weld
    %+  expect-eq  !>((from-eny:seed:rand %sm64 e))  !>((from-eny:seed:rand %sm64 e))
    %+  expect-eq  !>(`@`0xdeb.8a1b.ec35.d57d)  !>((fold-wide:seed:rand e))
  ==
::  Philox4x32-10 KAT: Random123 kat_vectors, all-zero case.
++  test-philox-kat-zero  ^-  tang
  %+  expect-eq
    !>(`@`0x9b00.dbd8.bc57.ac4c.e169.c58d.6627.e8d5)
    !>((block:philox:rand 0 0))
::  Philox4x32-10 KAT: Random123 kat_vectors, all-0xffffffff case.
++  test-philox-kat-ones  ^-  tang
  %+  expect-eq
    !>(`@`0x6d54.51fd.a20b.c7c6.41c8.3b0e.408f.276d)
    !>((block:philox:rand 0xffff.ffff.ffff.ffff 0xffff.ffff.ffff.ffff.ffff.ffff.ffff.ffff))
::  Philox4x32-10 KAT: Random123 kat_vectors, pi-digits case.
++  test-philox-kat-pi  ^-  tang
  %+  expect-eq
    !>(`@`0x2412.6ea1.5001.e420.94fd.cceb.d16c.fe09)
    !>((block:philox:rand 0x299f.31d0.a409.3822 0x370.7344.1319.8a2e.85a3.08d3.243f.6a88))
::  +next:philox: out = (block key ctr), ctr advances by 1, key unchanged.
++  test-philox-next  ^-  tang
  %+  expect-eq
    !>(`[out=@ p=phil:rand]`[0x9b00.dbd8.bc57.ac4c.e169.c58d.6627.e8d5 [key=0 ctr=1]])
    !>((next:philox:rand [key=0 ctr=0]))
::  +step dispatches to the right engine and keeps only bits [0,64) of a
::  Philox block.
++  test-step-phil  ^-  tang
  =/  r  (from-atom:seed:rand %phil 0)
  =^  out  r  (step:rand r)
  ;:  weld
    %+  expect-eq  !>(`@`0x24cb.d2fb.a9e3.9636)  !>(out)
    %+  expect-eq  !>(`rng:rand`[%phil p=[key=0xe220.a839.7b1d.cdaf ctr=1]])  !>(r)
  ==
++  test-step-sm64  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  out  r  (step:rand r)
  ;:  weld
    %+  expect-eq  !>(`@`0xe220.a839.7b1d.cdaf)  !>(out)
    %+  expect-eq  !>(`rng:rand`[%sm64 s=0x9e37.79b9.7f4a.7c15])  !>(r)
  ==
::  +fork: deterministic (same parent + salt -> same child).
++  test-fork-determinism  ^-  tang
  =/  r  (from-atom:seed:rand %phil 0)
  %+  expect-eq
    !>((fork:rand r 1))
    !>((fork:rand r 1))
::  +fork: different salts -> different children, for all three engines.
++  test-fork-distinct-salts  ^-  tang
  ;:  weld
    %+  expect-eq  !>(%.n)
      !>  =((fork:rand (from-atom:seed:rand %phil 0) 1) (fork:rand (from-atom:seed:rand %phil 0) 2))
    %+  expect-eq  !>(%.n)
      !>  =((fork:rand (from-atom:seed:rand %pcg 0) 1) (fork:rand (from-atom:seed:rand %pcg 0) 2))
    %+  expect-eq  !>(%.n)
      !>  =((fork:rand (from-atom:seed:rand %sm64 0) 1) (fork:rand (from-atom:seed:rand %sm64 0) 2))
  ==
::  +fork: %phil's ctr always resets to 0 in the child, regardless of the
::  parent's ctr.
++  test-fork-phil-ctr-reset  ^-  tang
  =/  r  [%phil p=[key=0x2a ctr=99]]
  =/  child  (fork:rand r 7)
  ?>  ?=(%phil -.child)
  %+  expect-eq  !>(0)  !>(ctr.p.child)
::  +fork: %pcg's child increment is always forced odd.
++  test-fork-pcg-inc-odd  ^-  tang
  =/  r  (from-atom:seed:rand %pcg 0)
  =/  child  (fork:rand r 5)
  ?>  ?=(%pcg -.child)
  %+  expect-eq  !>(1)  !>((dis 1 inc.p.child))
::  Nesting rule (rand-spec.md section 2.1c): the mix is genuinely
::  path-sensitive.  (fork (fork r a) b) must differ from
::  (fork (fork r b) a) and from (fork r (cat 6 a b)).
++  test-fork-path-sensitive  ^-  tang
  =/  r    (from-atom:seed:rand %phil 0)
  =/  a    11
  =/  b    22
  =/  fab  (fork:rand (fork:rand r a) b)
  =/  fba  (fork:rand (fork:rand r b) a)
  =/  fcat  (fork:rand r (cat 6 a b))
  ;:  weld
    %+  expect-eq  !>(%.n)  !>(=(fab fba))
    %+  expect-eq  !>(%.n)  !>(=(fab fcat))
    %+  expect-eq  !>(%.n)  !>(=(fba fcat))
  ==
::  ++gen door facade: +draw mirrors the functional +step, +fork mirrors
::  the functional +fork.
++  test-gen-draw  ^-  tang
  =/  g  ~(. gen:rand (from-atom:seed:rand %sm64 0))
  =^  x  g  draw:g
  %+  expect-eq  !>(`@`0xe220.a839.7b1d.cdaf)  !>(x)
++  test-gen-fork  ^-  tang
  =/  g  ~(. gen:rand (from-atom:seed:rand %phil 0))
  %+  expect-eq
    !>((fork:rand (from-atom:seed:rand %phil 0) 3))
    !>(r:(fork:g 3))
--
