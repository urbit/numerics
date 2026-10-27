::  /tests/lib/lagoon-pathological: jet/Hoon parity on edge-case inputs
::
::  Each arm runs one call through the jetted +la and through
::  /lib/lagoon-ref, a copy of /lib/lagoon with its ~% and ~/ hints
::  stripped so it never reaches a jet.  +par passes when both
::  produce the same noun or both crash.  A jet that crashes where the
::  Hoon succeeds, or returns where the Hoon crashes, fails here even
::  when the expected value is hard to write by hand.
::
::  The asymmetric +sub and +div arms also pin literal oracles: both
::  jets once swapped their operands (urbit/vere#1057, #1139), which
::  symmetric inputs such as 1/1 cannot detect.
::
::  Not covered, because each ends in a %fail that +mule cannot catch:
::  a +range that cannot terminate (the jet bails %fail by design), and
::  a zero-product shape with an indirect dimension such as
::  ~[0 (bex 64)] (the stdlib jets the Hoon reaches bail %fail on the
::  indirect size, so there is no deterministic Hoon result to match).
::  A M*P overflow in +mmul needs ~12GB of input and is not testable.
::
/-  *lagoon
/+  *test
/+  *lagoon
/+  ref=lagoon-ref
^|
|_  ~
::  +par: the jet and the reference agree, crash for crash.
++  par
  |*  [nam=@tas jet=(trap *) hon=(trap *)]
  ^-  tang
  =/  j  (mule jet)
  =/  h  (mule hon)
  ?:  &(?=(%| -.j) ?=(%| -.h))  ~
  ?:  &(?=(%& -.j) ?=(%& -.h))
    ?:  =(p.j p.h)  ~
    :~  leaf+"{<nam>} value mismatch"
        leaf+"  jet:  {<p.j>}"
        leaf+"  hoon: {<p.h>}"
    ==
  :~  leaf+"{<nam>} crash mismatch"
      leaf+"jet:  {?:(?=(%& -.j) "value" "crash")}"
      leaf+"hoon: {?:(?=(%& -.h) "value" "crash")}"
  ==
::  +bin: every binary elementwise op on (a, b).
++  bin
  |=  [a=ray b=ray]
  ^-  tang
  ;:  weld
    (par %add |.((add:la a b)) |.((add:la:ref a b)))
    (par %sub |.((sub:la a b)) |.((sub:la:ref a b)))
    (par %mul |.((mul:la a b)) |.((mul:la:ref a b)))
    (par %div |.((div:la a b)) |.((div:la:ref a b)))
    (par %mod |.((mod:la a b)) |.((mod:la:ref a b)))
    (par %gth |.((gth:la a b)) |.((gth:la:ref a b)))
    (par %gte |.((gte:la a b)) |.((gte:la:ref a b)))
    (par %lth |.((lth:la a b)) |.((lth:la:ref a b)))
    (par %lte |.((lte:la a b)) |.((lte:la:ref a b)))
    (par %dot |.((dot:la a b)) |.((dot:la:ref a b)))
  ==
::  +una: every unary op, and every ray-scalar op with scalar .n.
++  una
  |=  [a=ray n=@]
  ^-  tang
  ;:  weld
    (par %min |.((min:la a)) |.((min:la:ref a)))
    (par %max |.((max:la a)) |.((max:la:ref a)))
    (par %argmin |.((argmin:la a)) |.((argmin:la:ref a)))
    (par %argmax |.((argmax:la a)) |.((argmax:la:ref a)))
    (par %cumsum |.((cumsum:la a)) |.((cumsum:la:ref a)))
    (par %abs |.((abs:la a)) |.((abs:la:ref a)))
    (par %ravel |.((ravel:la a)) |.((ravel:la:ref a)))
    (par %transpose |.((transpose:la a)) |.((transpose:la:ref a)))
    (par %add-scalar |.((add-scalar:la a n)) |.((add-scalar:la:ref a n)))
    (par %sub-scalar |.((sub-scalar:la a n)) |.((sub-scalar:la:ref a n)))
    (par %mul-scalar |.((mul-scalar:la a n)) |.((mul-scalar:la:ref a n)))
    (par %div-scalar |.((div-scalar:la a n)) |.((div-scalar:la:ref a n)))
    (par %mod-scalar |.((mod-scalar:la a n)) |.((mod-scalar:la:ref a n)))
  ==
++  mmul-par
  |=  [a=ray b=ray]
  ^-  tang
  (par %mmul |.((mmul:la a b)) |.((mmul:la:ref a b)))
::  +flt: 1..n as %i754 scalars at .bloq.
++  flt
  |=  [n=@ =bloq]
  ^-  (list @)
  =/  ns  (gulf 1 n)
  ?+  bloq  !!
    %4  (turn ns sun:rh)
    %5  (turn ns sun:rs)
    %6  (turn ns sun:rd)
    %7  (turn ns sun:rq)
  ==
::  +grid: the m x n %i754 ray of 1..mn, row-major, at .bloq.
++  grid
  |=  [m=@ n=@ =bloq]
  ^-  ray
  =/  xs  (flt (^mul m n) bloq)
  %+  en-ray:la  [~[m n] bloq %i754 ~]
  %+  turn  (gulf 0 (dec m))
  |=(i=@ (swag [(^mul i n) n] xs))
::  +vec: the rank-1 %i754 ray of 1..n at .bloq.
++  vec
  |=  [n=@ =bloq]
  ^-  ray
  (en-ray:la [~[n] bloq %i754 ~] (flt n bloq))
::
::  Asymmetric operands, every %i754 bloq.  [6 1] - [3 4] = [3 -3];
::  [6 1] / [3 4] = [2 0.25].  The bloq-7 +div jet computed y/x.
++  test-sub-asym-all-bloqs  ^-  tang
  ;:  weld
    %+  expect-eq
      !>((en-ray:la [~[2] 4 %i754 ~] ~[.~~3 .~~-3]))
    !>((sub:la (en-ray:la [~[2] 4 %i754 ~] ~[.~~6 .~~1]) (en-ray:la [~[2] 4 %i754 ~] ~[.~~3 .~~4])))
    %+  expect-eq
      !>((en-ray:la [~[2] 5 %i754 ~] ~[.3 .-3]))
    !>((sub:la (en-ray:la [~[2] 5 %i754 ~] ~[.6 .1]) (en-ray:la [~[2] 5 %i754 ~] ~[.3 .4])))
    %+  expect-eq
      !>((en-ray:la [~[2] 6 %i754 ~] ~[.~3 .~-3]))
    !>((sub:la (en-ray:la [~[2] 6 %i754 ~] ~[.~6 .~1]) (en-ray:la [~[2] 6 %i754 ~] ~[.~3 .~4])))
    %+  expect-eq
      !>((en-ray:la [~[2] 7 %i754 ~] ~[.~~~3 .~~~-3]))
    !>((sub:la (en-ray:la [~[2] 7 %i754 ~] ~[.~~~6 .~~~1]) (en-ray:la [~[2] 7 %i754 ~] ~[.~~~3 .~~~4])))
  ==
::
++  test-div-asym-all-bloqs  ^-  tang
  ;:  weld
    %+  expect-eq
      !>((en-ray:la [~[2] 4 %i754 ~] ~[.~~2 .~~0.25]))
    !>((div:la (en-ray:la [~[2] 4 %i754 ~] ~[.~~6 .~~1]) (en-ray:la [~[2] 4 %i754 ~] ~[.~~3 .~~4])))
    %+  expect-eq
      !>((en-ray:la [~[2] 5 %i754 ~] ~[.2 .0.25]))
    !>((div:la (en-ray:la [~[2] 5 %i754 ~] ~[.6 .1]) (en-ray:la [~[2] 5 %i754 ~] ~[.3 .4])))
    %+  expect-eq
      !>((en-ray:la [~[2] 6 %i754 ~] ~[.~2 .~0.25]))
    !>((div:la (en-ray:la [~[2] 6 %i754 ~] ~[.~6 .~1]) (en-ray:la [~[2] 6 %i754 ~] ~[.~3 .~4])))
    %+  expect-eq
      !>((en-ray:la [~[2] 7 %i754 ~] ~[.~~~2 .~~~0.25]))
    !>((div:la (en-ray:la [~[2] 7 %i754 ~] ~[.~~~6 .~~~1]) (en-ray:la [~[2] 7 %i754 ~] ~[.~~~3 .~~~4])))
  ==
::
::  The same asymmetric pair through every binary op, every bloq.
++  test-asym-parity-all-bloqs  ^-  tang
  ;:  weld
    (bin (en-ray:la [~[2] 4 %i754 ~] ~[.~~7 .~~1]) (en-ray:la [~[2] 4 %i754 ~] ~[.~~2 .~~4]))
    (bin (en-ray:la [~[2] 5 %i754 ~] ~[.7 .1]) (en-ray:la [~[2] 5 %i754 ~] ~[.2 .4]))
    (bin (en-ray:la [~[2] 6 %i754 ~] ~[.~7 .~1]) (en-ray:la [~[2] 6 %i754 ~] ~[.~2 .~4]))
    (bin (en-ray:la [~[2] 7 %i754 ~] ~[.~~~7 .~~~1]) (en-ray:la [~[2] 7 %i754 ~] ~[.~~~2 .~~~4]))
  ==
::
::  Empty rays: a zero dimension leaves only the pin.  +bin-op,
::  +min/+max and +diag crash on them; the rest must agree too.
++  test-empty-binary  ^-  tang
  =/  z  (zeros:la [~[0 3] 5 %i754 ~])
  ;:  weld
    (bin z z)
    (bin (zeros:la [~[3 0] 6 %i754 ~]) (zeros:la [~[3 0] 6 %i754 ~]))
    (bin (zeros:la [~[0] 7 %i754 ~]) (zeros:la [~[0] 7 %i754 ~]))
  ==
::
++  test-empty-unary  ^-  tang
  ;:  weld
    (una (zeros:la [~[0 3] 5 %i754 ~]) .2)
    (una (zeros:la [~[3 0] 6 %i754 ~]) .~2)
    (una (zeros:la [~[0] 4 %i754 ~]) .~~2)
  ==
::
++  test-empty-diag-trace  ^-  tang
  =/  z  (zeros:la [~[0 0] 5 %i754 ~])
  ;:  weld
    (par %diag |.((diag:la z)) |.((diag:la:ref z)))
    (par %trace |.((trace:la z)) |.((trace:la:ref z)))
  ==
::
::  +check counts blocks, so an empty ray's single block may hold
::  junk below the pin: data 0x1234 at bloq 5 is one block, (dec 1)=0.
++  test-empty-junk-data  ^-  tang
  =/  j=ray  [[~[0 3] 5 %i754 ~] 0x1234]
  ;:  weld
    (bin j j)
    (una j .2)
  ==
::
::  The same on a non-empty ray: the pin block is 0x1234, not 0x1.
::  +check still passes (one extra block), so every op must agree on
::  what it does with the junk.
++  test-junk-pin  ^-  tang
  =/  g  (grid 2 2 5)
  =/  j=ray  g(data (con (lsh [5 4] 0x1234) (end [5 4] data.g)))
  =/  v  (vec 3 6)
  =/  k=ray  v(data (con (lsh [6 3] 0xdead) (end [6 3] data.v)))
  ;:  weld
    (bin j j)
    (una j .2)
    (mmul-par j j)
    (par %diag |.((diag:la j)) |.((diag:la:ref j)))
    (par %trace |.((trace:la j)) |.((trace:la:ref j)))
    (bin k k)
    (una k .~2)
  ==
::
::  %int2 rays take their own kernels: a junk pin, an empty ray, and
::  asymmetric operands.
++  test-int2-edges  ^-  tang
  =/  g  (en-ray:la [~[2 2] 5 %int2 ~] ~[~[7 1] ~[2 4]])
  =/  h  (en-ray:la [~[2 2] 5 %int2 ~] ~[~[2 4] ~[7 1]])
  =/  j=ray  g(data (con (lsh [5 4] 0x1234) (end [5 4] data.g)))
  =/  z  (zeros:la [~[0 3] 5 %int2 ~])
  ;:  weld
    (bin g h)
    (bin j j)
    (bin z z)
    (una j 2)
    (una z 2)
  ==
::
::  An empty %unum ray: +dot takes the quire path, not +bin-op, so the
::  jets must punt on empty rays rather than crash for every kind.
++  test-unum-empty-dot  ^-  tang
  =/  z  (zeros:la [~[0] 5 %unum ~])
  (par %dot |.((dot:la z z)) |.((dot:la:ref z z)))
::
::  +mmul with a zero dimension in each position.
++  test-mmul-zero-dims  ^-  tang
  ;:  weld
    (mmul-par (zeros:la [~[0 3] 5 %i754 ~]) (grid 3 2 5))
    (mmul-par (grid 2 3 5) (zeros:la [~[3 0] 5 %i754 ~]))
    (mmul-par (zeros:la [~[2 0] 5 %i754 ~]) (zeros:la [~[0 3] 5 %i754 ~]))
  ==
::
::  +mmul shape errors: inner mismatch, wrong rank, mixed bloqs.
++  test-mmul-bad-shapes  ^-  tang
  ;:  weld
    (mmul-par (grid 2 3 5) (grid 2 3 5))
    (mmul-par (vec 3 5) (grid 3 2 5))
    (mmul-par (grid 2 3 5) (vec 3 5))
    (mmul-par (reshape:la (vec 8 5) ~[2 2 2]) (grid 2 2 5))
    (mmul-par (grid 2 2 5) (grid 2 2 6))
  ==
::
::  Non-square +mmul (2x3 * 3x4) at every bloq: M != P.
++  test-mmul-nonsquare-all-bloqs  ^-  tang
  ;:  weld
    (mmul-par (grid 2 3 4) (grid 3 4 4))
    (mmul-par (grid 2 3 5) (grid 3 4 5))
    (mmul-par (grid 2 3 6) (grid 3 4 6))
    (mmul-par (grid 2 3 7) (grid 3 4 7))
  ==
::
::  +dot needs equal shapes, not just equal element counts.
++  test-dot-shapes  ^-  tang
  ;:  weld
    (par %dot |.((dot:la (grid 2 3 5) (grid 3 2 5))) |.((dot:la:ref (grid 2 3 5) (grid 3 2 5))))
    (par %dot |.((dot:la (grid 2 3 5) (grid 2 3 5))) |.((dot:la:ref (grid 2 3 5) (grid 2 3 5))))
    (par %dot |.((dot:la (grid 2 3 5) (grid 2 3 6))) |.((dot:la:ref (grid 2 3 5) (grid 2 3 6))))
  ==
::
::  Shape/data inconsistency: one element short and one too many.
++  test-inconsistent-data  ^-  tang
  =/  s=ray  [[~[2 2] 5 %i754 ~] (lsh [5 3] 1)]
  =/  l=ray  [[~[2 2] 5 %i754 ~] (lsh [5 5] 1)]
  ;:  weld
    (bin s s)
    (bin l l)
    (una s .2)
    (una l .2)
    (mmul-par s s)
  ==
::
::  A cell where a dimension belongs.  The vase is the only way past
::  the type system.  Jet side only: the Hoon +check reaches the
::  stdlib +mul jet with a cell, and u3wa_mul bails %fail, which
::  +mule cannot catch.  The lagoon jets must bail %exit instead.
++  test-cell-dimension  ^-  tang
  =/  c=ray  !<(ray [-:!>(*ray) [[[[1 2] 3 ~] 5 %i754 ~] 0x1]])
  ;:  weld
    (expect-fail |.((add:la c c)))
    (expect-fail |.((dot:la c c)))
    (expect-fail |.((mmul:la c c)))
    (expect-fail |.((cumsum:la c)))
    (expect-fail |.((min:la c)))
    (expect-fail |.((abs:la c)))
    (expect-fail |.((transpose:la c)))
    (expect-fail |.((add-scalar:la c .2)))
  ==
--
