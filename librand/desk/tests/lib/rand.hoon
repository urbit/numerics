  ::  /tests/lib/rand
::::
::    Milestones 1-5 (rand-spec.md): ++split-mix, ++philox, ++seed, +fork,
::    ++uni, ++pcg, ++dist.  The SplitMix64 and Philox4x32-10 KAT vectors
::    are checked against their published reference values (Vigna's
::    reference C for SplitMix64; the Random123 kat_vectors file --
::    all-zero, all-0xffffffff, and the pi-digits vector -- for Philox),
::    each cross-checked against an independent Python re-implementation
::    before being transcribed here.  ++uni's float arms (+rs/+rd/+rh/+rq/
::    +rs-oo/+rd-oo) are exact bit constructions (+sun then a power-of-two
::    multiply), so their expected values are likewise checked against an
::    independent Python IEEE-754 encoder, not the Hoon float printer.
::    ++pcg's KAT (seed=42, seq=54, pcg-c's own demo convention) is
::    checked against an independent Python re-implementation of pcg-c's
::    reference C (include/pcg_variants.h) -- this also exercises the
::    CORRECTED advance-then-output order (rand-spec.md section 3.3; an
::    earlier spec draft had it backwards).  ++dist's value spot-checks
::    involve /lib/math transcendentals (log/sqrt), so their expected
::    values are read off this same ship's actual output (sanity-checked
::    by hand against the closed-form algorithm first) rather than
::    independently rederived in Python, which isn't guaranteed to match
::    this codebase's correctly-rounded kernels to the last bit; ++dist's
::    moment tests are regression checks against the THEORETICAL targets
::    (mean/variance) at a fixed seed, with tolerances wide enough to
::    absorb real sampling noise at n=50000 but tight enough to catch an
::    actual algorithm bug.  +from-atom, +fold-wide, +mix, +fork, +bits,
::    and +below have no external reference (they are this library's own
::    design, or Lemire's algorithm applied to this library's own +step),
::    so those are checked by direct computation of the spec'd formula and
::    by determinism/order-sensitivity/path-sensitivity/unbiasedness
::    properties instead.
::
/+  *test,
    rand,
    math
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
::  +bits: n=4 -> low 4 bits of the first +step draw (0xe220a8397b1dcdaf).
++  test-bits  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  out  r  (bits:uni:rand r 4)
  ;:  weld
    %+  expect-eq  !>(`@`15)  !>(out)
    %+  expect-eq  !>(`rng:rand`[%sm64 s=0x9e37.79b9.7f4a.7c15])  !>(r)
  ==
::  +bits: n=130 (>64, spans 3 draws: two full 64-bit words, one 2-bit
::  remainder) round-trips against a manual reconstruction from the same
::  three +step draws taken independently.
++  test-bits-wide  ^-  tang
  =/  r0      (from-atom:seed:rand %sm64 0)
  =/  result  (bits:uni:rand r0 130)
  =^  w0  r0  (step:rand r0)
  =^  w1  r0  (step:rand r0)
  =^  w2  r0  (step:rand r0)
  =/  expected  :(add w0 (lsh [0 64] w1) (lsh [0 128] (end [0 2] w2)))
  ;:  weld
    %+  expect-eq  !>(expected)  !>(out.result)
    %+  expect-eq  !>(r0)  !>(r.result)
  ==
::  +below: Lemire, n=10, against an independently-computed reference
::  (same algorithm, computed in Python before transcription here).
++  test-below  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  out  r  (below:uni:rand r 10)
  ;:  weld
    %+  expect-eq  !>(`@`8)  !>(out)
    %+  expect-eq  !>(`rng:rand`[%sm64 s=0x9e37.79b9.7f4a.7c15])  !>(r)
  ==
::  +below: unbiasedness smoke test, n=6, 1200 draws, chi-square below a
::  fixed threshold (section 10 item 3: a smoke test, not TestU01).  Uses
::  ++fork to get 1200 independent draws without threading one rng
::  sequentially (each fork is its own single-draw +below call).  Chi-
::  square is computed as an exact integer comparison (sum-sq < k*expected)
::  rather than a division, so this needs no float/math import.
++  test-below-unbiased  ^-  tang
  =/  base      (from-atom:seed:rand %phil 0)
  =/  n         6
  =/  reps      1.200
  =/  expected  (div reps n)
  =/  counts
    %+  roll  (gulf 0 (dec reps))
    |=  [i=@ counts=(map @ @)]
    =/  child  (fork:rand base i)
    =^  draw   child  (below:uni:rand child n)
    (~(put by counts) draw +((fall (~(get by counts) draw) 0)))
  =/  sum-sq
    %+  roll  (gulf 0 (dec n))
    |=  [k=@ acc=@]
    =/  obs   (fall (~(get by counts) k) 0)
    =/  diff  (abs:si (dif:si (sun:si obs) (sun:si expected)))
    (add acc (mul diff diff))
  ::  chi-square = sum-sq / expected; threshold chosen generously (df=5,
  ::  p~0.001 critical value is ~20.5) since this is a smoke test, not a
  ::  rigorous statistical test -- compared as sum-sq < 30*expected to
  ::  avoid a division.
  %+  expect-eq  !>(%.y)  !>((lth sum-sq (mul 30 expected)))
::  +between: span=11 (a=-5,b=--5), same underlying draw as +test-below.
++  test-between  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  out  r  (between:uni:rand r -5 --5)
  %+  expect-eq  !>(`@s`--4)  !>(out)
::  +between: crashes if a > b.
++  test-between-bad-range  ^-  tang
  %-  expect-fail
  |.((between:uni:rand (from-atom:seed:rand %sm64 0) --5 -5))
::  +below: crashes on n=0.
++  test-below-zero  ^-  tang
  %-  expect-fail
  |.((below:uni:rand (from-atom:seed:rand %sm64 0) 0))
::  +rs/+rd/+rh/+rq: exact bit constructions, checked against an
::  independent Python IEEE-754 encoder (not the Hoon float printer).
++  test-float-rs  ^-  tang
  %+  expect-eq
    !>(`@rs`0x3dee.6d78)
    !>(out:(rs:uni:rand (from-atom:seed:rand %sm64 0)))
++  test-float-rd  ^-  tang
  %+  expect-eq
    !>(`@rd`0x3f95.072f.63b9.b5e0)
    !>(out:(rd:uni:rand (from-atom:seed:rand %sm64 0)))
++  test-float-rh  ^-  tang
  %+  expect-eq
    !>(`@rh`0x39af)
    !>(out:(rh:uni:rand (from-atom:seed:rand %sm64 0)))
++  test-float-rq  ^-  tang
  %+  expect-eq
    !>(`@rq`0x3ffd.3cd5.4372.cbe9.c441.5072.f63b.9b5e)
    !>(out:(rq:uni:rand (from-atom:seed:rand %sm64 0)))
++  test-float-rs-oo  ^-  tang
  %+  expect-eq
    !>(`@rs`0x3dee.6d7c)
    !>(out:(rs-oo:uni:rand (from-atom:seed:rand %sm64 0)))
++  test-float-rd-oo  ^-  tang
  %+  expect-eq
    !>(`@rd`0x3f95.072f.63b9.b5f0)
    !>(out:(rd-oo:uni:rand (from-atom:seed:rand %sm64 0)))
::  PCG64 KAT: pcg-c's own seeding (pcg_setseq_128_srandom_r, seed=42,
::  seq=54 -- the library's classic demo parameters), first 5 outputs,
::  against an independent Python re-implementation of pcg-c's reference
::  C (include/pcg_variants.h).  This exercises the CORRECTED operation
::  order (advance then output -- see rand-spec.md section 3.3), not the
::  order an earlier draft of the spec claimed.  Seeded directly with the
::  literal state/inc here rather than via +from-atom, since this KAT
::  checks +next:pcg's own algorithm against pcg-c's native seeding, not
::  this library's (unrelated, SplitMix64-derived) +from-atom scheme.
++  test-pcg-kat  ^-  tang
  =/  p  [state=0xde2b.ce05.be01.3be3.d3f6.c45a.41e5.4320 inc=0x6d]
  =^  o0  p  (next:pcg:rand p)
  =^  o1  p  (next:pcg:rand p)
  =^  o2  p  (next:pcg:rand p)
  =^  o3  p  (next:pcg:rand p)
  =^  o4  p  (next:pcg:rand p)
  ;:  weld
    %+  expect-eq  !>(`@`0x86b1.da1d.7206.2b68)  !>(o0)
    %+  expect-eq  !>(`@`0x1304.aa46.c985.3d39)  !>(o1)
    %+  expect-eq  !>(`@`0xa367.0e9e.0dd5.0358)  !>(o2)
    %+  expect-eq  !>(`@`0xf909.0e52.9a7d.ae00)  !>(o3)
    %+  expect-eq  !>(`@`0xc85b.9fd8.3799.6f2c)  !>(o4)
  ==
::  +advance: at a small, tractable delta (3), must equal 3 sequential
::  +next:pcg steps -- the property that justifies the O(log delta)
::  skip-ahead algorithm +jump relies on at the (untestable-by-brute-
::  force) delta = 2^64 scale.
++  test-pcg-advance  ^-  tang
  =/  p0  [state=5 inc=0x6d]
  =/  p1  (advance:pcg:rand p0 3)
  =^  s1  p0  (next:pcg:rand p0)
  =^  s2  p0  (next:pcg:rand p0)
  =^  s3  p0  (next:pcg:rand p0)
  %+  expect-eq  !>(state.p0)  !>(state.p1)
::  +jump: deterministic, and actually changes the state (catches a
::  no-op mistake).
++  test-pcg-jump  ^-  tang
  =/  p  [state=5 inc=0x6d]
  ;:  weld
    %+  expect-eq  !>((jump:pcg:rand p))  !>((jump:pcg:rand p))
    %+  expect-eq  !>(%.n)  !>(=(p (jump:pcg:rand p)))
  ==
::  +step dispatches %pcg correctly (next:pcg wired in, not the earlier
::  crash stub).
++  test-step-pcg  ^-  tang
  =/  r  (from-atom:seed:rand %pcg 0)
  =^  out  r  (step:rand r)
  ;:  weld
    %+  expect-eq  !>(`@`0x517a.36a9.6d93.79b8)  !>(out)
    %+  expect-eq
      !>(`rng:rand`[%pcg p=[state=0x15ab.4f3e.beb3.372a.3aed.9a13.0cdc.0020 inc=0x6e78.9e6a.a1b9.65f5]])
      !>(r)
  ==
::  ++dist: value spot-checks against actual on-ship computation (these
::  involve /lib/math transcendentals -- log/sqrt -- so the expected
::  values below are read off this same ship's output, not independently
::  rederived in Python, which isn't guaranteed to match this codebase's
::  correctly-rounded kernels to the last bit).  Each is sanity-checked
::  against the closed-form algorithm by hand (e.g. -ln(u)/lambda for
::  +expon) before being trusted as a regression constant.
++  test-normal  ^-  tang
  %+  expect-eq
    !>(`@rd`.~-0.9479938949723624)
    !>(out:(normal:dist:rand (from-atom:seed:rand %sm64 0)))
::  determinism: same seed -> same output (spot check; +normal's rejection
::  loop is the one arm in ++dist where a threading bug would most likely
::  show up as nondeterminism).
++  test-normal-determinism  ^-  tang
  %+  expect-eq
    !>((normal:dist:rand (from-atom:seed:rand %phil 7)))
    !>((normal:dist:rand (from-atom:seed:rand %phil 7)))
++  test-normal-mv  ^-  tang
  %+  expect-eq
    !>(`@rd`.~8.104012210055275)
    !>(out:(normal-mv:dist:rand (from-atom:seed:rand %sm64 0) .~10 .~2))
++  test-normal-mv-bad-sigma  ^-  tang
  %-  expect-fail
  |.((normal-mv:dist:rand (from-atom:seed:rand %sm64 0) .~0 .~-1))
::  -ln(u)/lambda, u = +rd-oo of the first SplitMix64 output (0.0205...):
::  -ln(0.0205352...)/2 ~ 1.9428 -- matches to displayed precision.
++  test-expon  ^-  tang
  %+  expect-eq
    !>(`@rd`.~1.942806871605541)
    !>(out:(expon:dist:rand (from-atom:seed:rand %sm64 0) .~2))
++  test-expon-bad-rate  ^-  tang
  %-  expect-fail
  |.((expon:dist:rand (from-atom:seed:rand %sm64 0) .~0))
::  alpha=2 (>=1, direct Marsaglia-Tsang path).
++  test-gamma-ge1  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.7179344171650754)
    !>(out:(gamma:dist:rand (from-atom:seed:rand %sm64 0) .~2))
::  alpha=0.5 (<1, boost path: gamma(1.5) * u^(1/0.5)).
++  test-gamma-lt1  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.05447897345403265)
    !>(out:(gamma:dist:rand (from-atom:seed:rand %sm64 0) .~0.5))
++  test-gamma-bad-shape  ^-  tang
  %-  expect-fail
  |.((gamma:dist:rand (from-atom:seed:rand %sm64 0) .~0))
++  test-beta  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.23622515961139348)
    !>(out:(beta:dist:rand (from-atom:seed:rand %sm64 0) .~2 .~3))
++  test-chi2  ^-  tang
  %+  expect-eq
    !>(`@rd`.~0.8261342997997718)
    !>(out:(chi2:dist:rand (from-atom:seed:rand %sm64 0) 3))
++  test-chi2-bad-df  ^-  tang
  %-  expect-fail
  |.((chi2:dist:rand (from-atom:seed:rand %sm64 0) 0))
++  test-student-t  ^-  tang
  %+  expect-eq
    !>(`@rd`.~-0.7137613421937223)
    !>(out:(student-t:dist:rand (from-atom:seed:rand %sm64 0) 5))
++  test-student-t-bad-df  ^-  tang
  %-  expect-fail
  |.((student-t:dist:rand (from-atom:seed:rand %sm64 0) 0))
::  u ~ 0.0205 < 0.5 -> %.y.
++  test-bernoulli  ^-  tang
  %+  expect-eq
    !>(`?`%.y)
    !>(out:(bernoulli:dist:rand (from-atom:seed:rand %sm64 0) .~0.5))
++  test-bernoulli-bad-prob  ^-  tang
  %-  expect-fail
  |.((bernoulli:dist:rand (from-atom:seed:rand %sm64 0) .~1.5))
++  test-geometric  ^-  tang
  %+  expect-eq
    !>(`@ud`11)
    !>(out:(geometric:dist:rand (from-atom:seed:rand %sm64 0) .~0.3))
::  p=1 edge case: returns 1 without consuming a draw (rng unchanged).
++  test-geometric-p1  ^-  tang
  =/  r  (from-atom:seed:rand %sm64 0)
  =^  out  r  (geometric:dist:rand r .~1)
  ;:  weld
    %+  expect-eq  !>(`@ud`1)  !>(out)
    %+  expect-eq  !>((from-atom:seed:rand %sm64 0))  !>(r)
  ==
++  test-geometric-bad-prob  ^-  tang
  %-  expect-fail
  |.((geometric:dist:rand (from-atom:seed:rand %sm64 0) .~0))
::  Moment tests (rand-spec.md section 10 item 4): sample mean/variance at
::  50k draws, FIXED seed -- deterministic, so this is a regression test
::  against the theoretical constants, not a flaky statistical one.
::  Tolerances are generous relative to the actual sampling error at
::  n=50000 (stderr(mean) ~ sigma/sqrt(n), stderr(var) ~ sigma^2*sqrt(2/n))
::  so a real algorithm bug (wrong constant, sign error, wrong formula)
::  fails loudly while ordinary run-to-run noise from a fixed seed does not
::  (there being only one fixed seed, "run-to-run" here really means
::  "robust to this exact computation shifting by a few ULP if the
::  underlying math kernels ever change").
++  test-moments-normal  ^-  tang
  =/  mv  (moments (from-atom:seed:rand %phil 0) normal:dist:rand 50.000)
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.02]) -.mv .~0))
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) +.mv .~1))
  ==
++  test-moments-expon  ^-  tang
  =/  mv  (moments (from-atom:seed:rand %phil 0) |=(r=rng:rand (expon:dist:rand r .~1)) 50.000)
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) -.mv .~1))
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) +.mv .~1))
  ==
++  test-moments-gamma  ^-  tang
  =/  mv  (moments (from-atom:seed:rand %phil 0) |=(r=rng:rand (gamma:dist:rand r .~2)) 50.000)
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) -.mv .~2))
    %+  expect-eq  !>(%.y)  !>((~(is-close rd:math [%n .~0 .~0.05]) +.mv .~2))
  ==
::  +moments: draw .n samples via .f from .r0, return [mean=@rd var=@rd]
::  (population variance).  Shared by the three moment tests above.
++  moments
  |=  [r0=rng:rand f=$-(rng:rand [@rd rng:rand]) n=@]
  ^-  [mean=@rd var=@rd]
  =/  mth  ~(. rd:math [%n .~1e-13 .~0])
  =/  r    r0
  =/  i    0
  =/  sum  .~0
  =/  ssq  .~0
  |-  ^-  [@rd @rd]
  ?:  =(i n)
    =/  ct    (sun:mth n)
    =/  mean  (div:mth sum ct)
    [mean (sub:mth (div:mth ssq ct) (mul:mth mean mean))]
  =^  x  r  (f r)
  %=  $
    i    +(i)
    sum  (add:mth sum x)
    ssq  (add:mth ssq (mul:mth x x))
  ==
--
