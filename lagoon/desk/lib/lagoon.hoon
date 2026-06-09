/-  lagoon
/+  twoc, unum, complex, fixed
=+  lagoon
::                                                    ::
::::                    ++la                          ::  (2v) vector/matrix ops
~%  %non  ..part  ~  :: nest non in hex for now
|%
::  +lake: set +la params
::
::    rnd: rounding mode
::
::  Limited to the four modes the hoon @rs/@rd/@rq/@rh precision doors
::  accept: %n (nearest, ties to even), %u (toward +inf), %d (toward
::  -inf), %z (toward zero).  The underlying ++fl engine — and the C
::  jets' _set_rounding — also implement %a (nearest, ties away from
::  zero), but the precision doors do not expose it, so %a is
::  intentionally excluded here; the jets handle it defensively only.
+$  rounding-mode  ?(%n %u %d %z)
++  lake
  |=  [inrnd=rounding-mode]
  %*(. la rnd inrnd)
::
++  la
  ^|
  =+  [rnd=*rounding-mode]
  ~/  %lagoon
  |%
  ::
  ::  Utilities
  ::
  ++  print  |=(a=ray ~>(%slog.1^(to-tank (ravel a) shape.meta.a kind.meta.a) ~))
  ::
  ++  slog   |=(a=ray (^slog (to-tank (ravel a) shape.meta.a kind.meta.a) ~))
  ::
  ::
  ++  to-tank
    |=  [dat=(list @) shape=(list @) =kind]
    ^-  tank
    ::  1D vector case
    ?:  =(1 (lent shape))
      :+  %rose  [" " "[" "]"]
      %+  turn  dat
      |=  i=@
      ^-  tank
      (sell [%atom kind ~] i)
    ::  general tensor case
    =/  width  (^div (lent dat) -.shape)
    :+  %rose  [" " "[" "]\0a"]
    ^-  (list tank)
    %+  turn  (gulf 0 (dec -.shape))
    |=  i=@
    (to-tank (swag [(^mul width i) width] dat) +.shape kind)
  ::
  ++  get-term
    |=  =meta
    ?-    kind.meta
        %uint
      %ud
      ::
        %int2
      %sd
      ::
        %i754
      ?+    bloq.meta  ~|(bloq.meta !!)
        %7  %rq
        %6  %rd
        %5  %rs
        %4  %rh
      ==
      ::
        %unum
      ::  no decimal pretty-printer yet; the @rpb/@rph/@rps/@rpd auras are
      ::  not registered with the printer, so values render as raw hex.
      ?+    bloq.meta  ~|(bloq.meta !!)
        %6  %rpd
        %5  %rps
        %4  %rph
        %3  %rpb
      ==
      ::
        %cplx
      ::  packed complex; render as the @c aura (no decimal printer yet).
      ?+    bloq.meta  ~|(bloq.meta !!)
        %8  %cq
        %7  %cd
        %6  %cs
        %5  %ch
      ==
        %fixp
      ::  no fixed-point printer; render the raw stored bits as hex.
      %ux
    ==
  ::
  ++  squeeze  |*(a=* `(list _a)`~[a])
  ++  submatrix
    ::
    ::  +submatrix: grab a submatrix using numpy slice notation, except indices are inclusive.
    ::  If you aren't slicing the last few dimensions, you can omit them.
    ::
    ~/  %submatrix
    |=  [sli=(list slice) a=ray]
    ::
    ::  example: sli=~[`[`1 `3] `[`1 ~] ~] is equivalent to a[1:3,1:,:]
    ::
    ^-  ray
    ::
    ::  pad slice with sigs if necessary
    =?  sli  (^lth (lent sli) (lent shape.meta.a))
      (weld sli (reap (^sub (lent shape.meta.a) (lent sli)) *(unit [(unit @) (unit @)])))
    ::
    ::  calculate indices to grab using cartesian product
    =/  out-indices=(list (list (list @)))
    %+  turn
      (gulf 0 (dec (lent shape.meta.a)))
    |=  i=@
    =/  s  (snag i sli)
    =/  dim  (snag i shape.meta.a)
    ?~  s
      (turn (gulf 0 (dec dim)) squeeze)
    =/  s2=[(unit @) (unit @)]  (need s)
    =/  c=^  [(fall -.s2 ~) (fall +.s2 ~)]
    ^-  (list (list @))
    ?+    c  !!
        [j=@ k=~]  (turn (gulf j.c (dec dim)) squeeze)
        [j=@ k=@]  (turn (gulf j.c k.c) squeeze)
        [j=~ k=@]  (turn (gulf 0 k.c) squeeze)
        [~ ~]  (turn (gulf 0 (dec dim)) squeeze)
    ==
    ::
    ::  calculate the shape of the result
    =/  out-shape=(list @)
    %+  turn
      out-indices
    |=(inds=(list (list @)) (lent inds))
    :: 
    ::  grab submatrix entries from cartesian product
    =/  new-dat=@ux
    %+  rep  bloq.meta.a
    %+  turn
      (gather out-indices)
    |=  dex=(list @)
    (get-item a dex)
    ::
    ::  construct new ray
    %-  spac
    =,  meta.a
    :-  [out-shape bloq kind ~]
    new-dat
  ::
  ++  product  ::  cartesian product
    |*  [a=(list) b=(list)]
    ?~  a
      b
    %-  zing
    %+  turn  a
    |=  ai=_-.a
    (turn b |=(bi=_-.b (welp ai bi)))
  ::
  ++  gather
    |=  [a=(list (list (list @)))]
    ^-  (list (list @))
    =/  i  0
    =|  c=(list (list @)) 
    |-
    ?:  =(i (lent a))
      c
    $(i +(i), c `(list (list @))`(product c (snag i a)))
  ::
  ++  get-item  ::  extract item at index .dex
    |=  [=ray dex=(list @)]
    ^-  @ux
    =/  len  (^sub (roll shape.meta.ray ^mul) 1)
    %^    cut
        bloq.meta.ray
      :: [(^sub len (get-bloq-offset meta.ray dex)) 1]
      [(get-bloq-offset meta.ray dex) 1]
    data.ray
  ::
  ++  set-item  ::  set item at index .dex to .val
    |=  [=ray dex=(list @) val=@]
    ^+  ray
    =/  len  (^sub (roll shape.meta.ray ^mul) 1)
    :-  meta.ray
    %^    sew
        bloq.meta.ray
      :: [(^sub len (get-bloq-offset meta.ray dex)) 1 val]
      [(get-bloq-offset meta.ray dex) 1 val]
    data.ray
  ::
  ++  get-row
    |=  [a=ray dex=(list @)]
    ^-  ray
    =,  meta.a
    ?:  =(1 (lent shape))
      (spac [~[1] bloq kind ~] (get-item a dex))
    ?>  =(+((lent dex)) (lent shape))
    =/  res
      %-  zeros
      :*  ~[1 (snag 0 (flop shape))]
          bloq
          kind
          ~
      ==
    =/  idx  0
    |-  ^-  ray
    ?:  =((snag 0 (flop shape.meta.res)) idx)  res
    =/  val  (get-item a (weld dex ~[idx]))
    $(idx +(idx), res (set-item res ~[0 idx] val))
  ::
  ++  set-row
    |=  [a=ray dex=(list @) row=ray]
    ^-  ray
    ?:  &(=(1 (lent dex)) =(1 (lent shape.meta.row)))  (set-item a dex (get-item row ~[0]))
    ?>  =(+((lent dex)) (lent shape.meta.a))
    ?>  =((lent shape.meta.row) 2)
    ?>  =(1 (snag 0 shape.meta.row))
    ?>  =((snag 1 shape.meta.row) (snag 0 (flop shape.meta.a)))
    =/  idx  0
    |-  ^-  ray
    ?:  =((snag 1 shape.meta.row) idx)  a
    %=  $
      idx  +(idx)
      a    (set-item a (weld dex ~[idx]) (get-item row ~[0 idx]))
    ==
  :: only works for 2D
  ++  get-col
    |=  [a=ray dex=(list @)]
    ^-  ray
    (get-row (transpose a) dex)
  ::
  ++  set-col
    |=  [a=ray dex=(list @) col=ray]
    ^-  ray
    (transpose (set-row (transpose a) dex col))
  ::
  ++  get-bloq-offset  ::  get bloq offset of n-dimensional index
    |=  [=meta dex=(list @)]
    ^-  @
    (get-item-number shape.meta dex)
  ::
  ++  get-item-number  ::  convert n-dimensional index to scalar index
    |=  [sap=(list @) dex=(list @)]
    ^-  @
    =/  cof  1
    =/  ret  0
    =.  sap  (flop sap)
    =.  dex  (flop dex)
    |-  ^+  ret
    ?~  sap  ret :: out of indices, return
    ?~  dex  !!  :: no indices past size
    ?>  (^lth i.dex i.sap)  :: index must be less than size
    %=  $
      sap  t.sap
      dex  t.dex
      cof  (^mul cof i.sap)
      ret  (^add ret (^mul i.dex cof))
    ==
  ::  Return the stride in each dimension:  row, col, layer, &c.
  ::  The stride is reported in units of bits.
  ::
  ++  strides
    |=  =meta
    ^-  (list @)
    =/  idx  0
    =|  res=(list @)
    |-  ^-  (list @)
    ?~  shape.meta  (flop res)
    =/  stride  (roll (scag idx `(list @)`shape.meta) ^mul)
    %=  $
      idx         +(idx)
      shape.meta  t.shape.meta
      res         [(^mul (pow 2 bloq.meta) stride) res]
    ==
  ::
  ++  get-dim  :: convert scalar index to n-dimensional index
    |=  [shape=(list @) ind=@]
    =/  shape  (flop shape)
    =/  i=@  0
    =|  res=(list @)
    ?>  (^lth ind (roll shape ^mul))
    |-  ^-  (list @)
    ?:  (^gte i (lent shape))  (flop res)
    %=    $
      res  `(list @)`(snoc res (^mod ind (snag i shape)))
      ind  (^div ind (snag i shape))
      i    (^add i 1)
    ==
  ::
  ++  get-item-index
    |=  [shape=(list @) num=@]
    ^-  @
    =/  len  (roll shape ^mul)
    =-  (roll - ^add)
    ^-  (list @)
    %+  turn  shape
    |=  wid=@
    (^mod (^div len wid) num)
  ::
  ++  ravel
    :: ~/  %ravel
    |=  a=ray
    ^-  (list @)
    (snip (rip bloq.meta.a data.a))
  ::
  ++  en-ray    :: baum to ray
    |=  =baum
    ^-  ray
    =/  a=ray  (zeros meta.baum)
    =/  i  0
    =/  n  (roll shape.meta.a ^mul)
    |-
    ?:  =(i n)  a
    %=  $
      i  +(i)
      data.a
        %+  con
          data.a
        %+  lsh
          [bloq.meta.a i]
        (get-item-baum baum (get-dim shape.meta.a i))
    ==
  ::
  ++  de-ray    :: ray to baum
    |=  =ray
    ^-  baum
    |^
    :-  meta.ray
    ^-  ndray
    ::
    =,  meta.ray
    ?:  =(1 (lent shape))
      ::  Snip off tail which is the pinned 0x1 MSB
      (snip (rip bloq data.ray))
    ::
    ?:  =(2 (lent shape))
      =/  dims  (flop shape)
      =|  fin=(list ndray)
      =|  els=ndray
      |-
      ?:  =(0x1 data.ray)  (welp ;;((list ndray) fin) ~[;;((list ndray) els)])
      %=  $
        els   (snip (rip bloq (cut bloq [0 +((snag 0 dims))] data:(spac `^ray`[[~[(snag 0 dims) 1] bloq kind tail] `@ux`data.ray]))))
        fin   ?~  els  fin  :: skip on first row
              ?~  fin  `(list (list ndray))`~[;;((list ndray) els)]
              (welp ;;((list (list ndray)) fin) ~[;;((list ndray) els)])
        data.ray  (rsh [bloq (snag 0 dims)] data.ray)
      ==
    !!
    ::  cut off end
    ++  rip
      |=  [a=bite b=@]
      ^-  (list @)
      =/  cnt  ?@(a 1 +.a)
      ?~  (^rip a b)  (reap cnt 0)
      (^rip a b)
    --
  ::
  ::  +check: validate a ray's data length against its shape.  The data atom
  ::  carries a 0x1 most-significant "pin" above the elements (so leading-zero
  ::  elements aren't lost), hence (met bloq data) counts prod(shape)+1 blocks;
  ::  this asserts prod(shape) == (met ...) - 1.  Every public arm calls it.
  ++  check
    |=  =ray
    ^-  ?
    .=  (roll shape.meta.ray ^mul)
    (dec (met bloq.meta.ray data.ray))
  ::
  ++  get-item-baum
    |=  [=baum dex=(list @)]
    ^-  @
    =/  a=ndray  data.baum
    |-
    ?~  a  !!
    ?@  (snag -.dex ;;((list ndray) a))
      ;;(@ (snag -.dex ((list ndray) a)))
    %=  $
      dex  +.dex
      a    (snag -.dex ;;((list ndray) a))
    ==
  ::
  ++  fill
    |=  [=meta x=@]
    ^-  ray
    =/  len  (roll shape.meta ^mul)
    :-  meta
    (con +:(zeros meta) (fil bloq.meta len x))
  ::
  ++  spac
    |=  =ray
    ^-  ^ray
    :-  meta.ray
    (con data:(zeros meta.ray) data.ray)
  ::
  ++  unspac
    |=  =ray
    ^-  ^ray
    :-  meta.ray
    (cut bloq.meta.ray [0 (roll shape.meta.ray ^mul)] data.ray)
  ::
  ++  scalar-to-ray
    |=  [=meta data=@]
    ^-  ray
    =.  shape.meta  (reap (lent shape.meta) 1)
    %-  spac
    [meta data]
  ::
  ++  change
    |=  [=ray =kind =bloq]
    ^-  ^ray
    ?+    kind.meta.ray  !!
        %uint
      ?+    kind  !!
          :: %uint -> %uint
          %uint
        %-  en-ray
        :-  [shape.meta.ray bloq %uint tail.meta.ray]
        data:(de-ray ray)
          :: %uint -> %i754
          %i754
        %-  en-ray
        :-  [shape.meta.ray bloq %i754 tail.meta.ray]
        %+  turn  (ravel ray)
        ?+  bloq  !!
          %7  ~(sun rq rnd)
          %6  ~(sun rd rnd)
          %5  ~(sun rs rnd)
          %4  ~(sun rh rnd)
        ==
      ==
      ::
        %i754
      ?+    kind  !!
          :: %i754 -> %uint.  KNOWN LIMITATION: a negative %i754 value has no
          :: %uint representation; +toi yields a negative @s and the (div ... 2)
          :: wraps it rather than clamping/erroring.  Callers must pre-clamp.
          %uint
        %-  en-ray
        :-  [shape.meta.ray bloq %uint tail.meta.ray]
        %+  turn  (ravel ray)
        ?+  bloq.meta.ray  !!
          %7  |=(a=@rq ^-(@u (^div (need (~(toi rq rnd) a)) 2)))
          %6  |=(a=@rd ^-(@u (^div (need (~(toi rd rnd) a)) 2)))
          %5  |=(a=@rs ^-(@u (^div (need (~(toi rs rnd) a)) 2)))
          %4  |=(a=@rh ^-(@u (^div (need (~(toi rh rnd) a)) 2)))
        ==
          :: %i754 -> %i754
          %i754
        ?>  &((^gte bloq %4) (^lte bloq %7))
        %-  en-ray
        :-  [shape.meta.ray bloq %i754 tail.meta.ray]
        data:(de-ray ray)
      ==
      ::  %unum conversion (to/from posits) is not yet wired; see the unum
      ::  to-r*/from-r* matrix in /lib/unum for the eventual mapping.
        %unum
      ~|('lagoon: change/convert not yet implemented for %unum' !!)
      ::  out of complex is lossy and ambiguous (real part? modulus?); refuse
      ::  and make the caller pick (abs for modulus).  Into-complex TODO.
        %cplx
      ~|('lagoon: change out of %cplx is lossy; use abs or an explicit re/im' !!)
      ::  %fixp conversion not yet wired; use /lib/fixed's to-rs/from-rs.
        %fixp
      ~|('lagoon: change/convert not yet implemented for %fixp' !!)
    ==
  ::
  ::  Builders
  ::
  ::
  ++  eye      ::  produces identity matrix of shape nxn.
    |=  =meta
    ^-  ray
    ~_  leaf+"lagoon-fail"
    ?>  =(2 (lent shape.meta))
    ?>  =((snag 0 shape.meta) (snag 1 shape.meta))
    =/  n  (snag 0 shape.meta)
    =<  +
    %^    spin
        (gulf 0 (dec n))
      ^-  ray  (zeros [~[n n] bloq.meta kind.meta ~])
    |=  [i=@ r=ray]
    :: [i (set-item r ~[i i] 1)]
    :-  i
    %^  set-item
        r
      ~[i i]
    ^-  @
    ?-    kind.meta
        %uint  `@`1
      ::
        %int2  `@`1
      ::
        %i754
      ?+  bloq.meta  ~|(bloq.meta !!)
        %7  .~~~1
        %6  .~1
        %5  .1
        %4  .~~1
      ==
      ::
        %unum
      ?+  bloq.meta  ~|(bloq.meta !!)
        %6  one:rpd:unum
        %5  one:rps:unum
        %4  one:rph:unum
        %3  one:rpb:unum
      ==
      ::
        %cplx
      ?+  bloq.meta  ~|(bloq.meta !!)
        %8  ~(one cq:complex rnd)
        %7  ~(one cd:complex rnd)
        %6  ~(one cs:complex rnd)
        %5  ~(one ch:complex rnd)
      ==
      ::  fixed-point 1.0 = 2^b, with b the fractional bits from meta.tail.
        %fixp  (bex +:;;([a=@ b=@] tail.meta))
    ==
  ::    Zeroes
  ++  zeros
    |=  =meta  ^-  ray
    ~_  leaf+"lagoon-fail"
    :-  meta
    (lsh [bloq.meta (roll shape.meta ^mul)] 1)
  ::    Ones
  ++  ones
    |=  =meta  ^-  ray
    ~_  leaf+"lagoon-fail"
    =/  one
      ?-    kind.meta
          %uint  `@`1
        ::
          %int2  `@`1
        ::
          %i754
        ?+  bloq.meta  !!
          %7  .~~~1
          %6  .~1
          %5  .1
          %4  .~~1
        ==
        ::
          %unum
        ?+  bloq.meta  !!
          %6  one:rpd:unum
          %5  one:rps:unum
          %4  one:rph:unum
          %3  one:rpb:unum
        ==
        ::
          %cplx
        ?+  bloq.meta  !!
          %8  ~(one cq:complex rnd)
          %7  ~(one cd:complex rnd)
          %6  ~(one cs:complex rnd)
          %5  ~(one ch:complex rnd)
        ==
        ::  fixed-point 1.0 = 2^b (b = fractional bits from meta.tail).
          %fixp  (bex +:;;([a=@ b=@] tail.meta))
      ==
    (fill meta one)
  ::  Produce a 1-dimensional index array.
  ::  Only produces %uint.
  ::  Note that this runs from 0 to n-1.  The point of ++iota is to be an index,
  ::  so it needs to pattern-match the context rather than slavishly follow APL.
  ::
  ++  iota
    |=  =meta
    ^-  ray
    ?>  =((lent shape.meta) 1)
    =/  n  (snag 0 shape.meta)
    =.  kind.meta  %uint
    (en-ray meta (gulf 0 (dec n)))
  ::  Produce a magic square in nD.
  ++  magic
    |=  =meta
    ^-  ray
    =/  n  (roll shape.meta ^mul)
    %+  reshape
      (en-ray [~[n] bloq.meta %uint ~] (gulf 0 (dec n)))
    shape.meta
  ::  Produce a 1-dimensional range along one dimension
  ::  as [a, b) with interval d.
  ::  Only produces %i754.
  ::
  ++  range
    ~/  %range
    |=  [=meta [a=@ b=@] d=@]
    ^-  ray
    =.  kind.meta  %i754
    %-  spac
    %-  en-ray
    ::
    ?+    bloq.meta  !!
        %7
      =/  ba  (~(sub rq rnd) b a)
      =/  bad  `(list @)`~[a]
      |-  ^-  baum
      ?:  (?:((~(lth rq rnd) ba .~~~0) ~(lte rq rnd) ~(gte rq rnd)) (~(add rq rnd) (snag 0 bad) d) b)
        =.  shape.meta  ~[(lent bad)]
        [meta (flop bad)]
      $(bad [(~(add rq rnd) (snag 0 bad) d) bad])
      ::
        %6
      =/  ba  (~(sub rd rnd) b a)
      =/  bad  `(list @)`~[a]
      |-  ^-  baum
      ?:  (?:((~(lth rd rnd) ba .~0) ~(lte rd rnd) ~(gte rd rnd)) (~(add rd rnd) (snag 0 bad) d) b)
        =.  shape.meta  ~[(lent bad)]
        [meta (flop bad)]
      $(bad [(~(add rd rnd) (snag 0 bad) d) bad])
      ::
        %5
      =/  ba  (~(sub rs rnd) b a)
      =/  bad  `(list @)`~[a]
      |-  ^-  baum
      ?:  (?:((~(lth rs rnd) ba .0) ~(lte rs rnd) ~(gte rs rnd)) (~(add rs rnd) (snag 0 bad) d) b)
        =.  shape.meta  ~[(lent bad)]
        [meta (flop bad)]
      $(bad [(~(add rs rnd) (snag 0 bad) d) bad])
      ::
        %4
      =/  ba  (~(sub rh rnd) b a)
      =/  bad  `(list @)`~[a]
      |-  ^-  baum
      ?:  (?:((~(lth rh rnd) ba .~~0) ~(lte rh rnd) ~(gte rh rnd)) (~(add rh rnd) (snag 0 bad) d) b)
        [meta (flop bad)]
      $(bad [(~(add rh rnd) (snag 0 bad) d) bad])
    ==
  ::  Produce a 1-dimensional range along one dimension
  ::  as [a b] with number of steps n.
  ::  Only produces %i754.
  ::
  ++  linspace
    ~/  %linspace
    |=  [=meta [a=@ b=@] n=@ud]
    ^-  ray
    =.  shape.meta  ~[n]
    =.  kind.meta  %i754
    ?<  =(n 0)
    ?:  =(n 1)  (en-ray meta ~[a])
    %-  en-ray
    :-  meta
    ::
    ?+    bloq.meta  !!
        %7
      =/  ba  (~(sub rq rnd) b a)
      =/  d  (~(div rq rnd) ba (~(sun rq rnd) (dec n)))
      =|  bad=(list @)
      =|  i=@ud
      |-  ^-  ndray
      ?:  (^lte n +(i))  (flop [b bad])
      $(i +(i), bad [(~(add rq rnd) a (~(mul rq rnd) (~(sun rq rnd) i) d)) bad])
      ::
        %6
      =/  ba  (~(sub rd rnd) b a)
      =/  d  (~(div rd rnd) ba (~(sun rd rnd) (dec n)))
      =|  bad=(list @)
      =|  i=@ud
      |-  ^-  ndray
      ?:  (^lte n +(i))  (flop [b bad])
      $(i +(i), bad [(~(add rd rnd) a (~(mul rd rnd) (~(sun rd rnd) i) d)) bad])
      ::
        %5
      =/  ba  (~(sub rs rnd) b a)
      =/  d  (~(div rs rnd) ba (~(sun rs rnd) (dec n)))
      =|  bad=(list @)
      =|  i=@ud
      |-  ^-  ndray
      ?:  (^lte n +(i))  (flop [b bad])
      $(i +(i), bad [(~(add rs rnd) a (~(mul rs rnd) (~(sun rs rnd) i) d)) bad])
      ::
        %4
      =/  ba  (~(sub rh rnd) b a)
      =/  d  (~(div rh rnd) ba (~(sun rh rnd) (dec n)))
      =|  bad=(list @)
      =|  i=@ud
      |-  ^-  ndray
      ?:  (^lte n +(i))  (flop [b bad])
      $(i +(i), bad [(~(add rh rnd) a (~(mul rh rnd) (~(sun rh rnd) i) d)) bad])
    ==
  ::  Coerce 1D array along specified dimension with given overall
  ::  dimensionality.
  ::
  ++  urge
    |=  [=ray [i=@ud n=@ud]]
    ^-  ^ray
    ?>  =(1 (lent shape.meta.ray))
    =.  shape.meta.ray  `(list @)`(zing (reap n ~[1]))
    =.  shape.meta.ray  (snap shape.meta.ray i n)
    |-
    ray
  ::  Produce an n-dimensional array containing a single value.
  ::
  ++  scale
    |=  [=meta data=@]
    ^-  ray
    =.  shape.meta  `(list @)`(zing (reap (lent shape.meta) ~[1]))
    ?-    kind.meta
        %uint
      (spac [meta `@ux`data])
      ::
        %int2
      (spac [meta `@ux`data])
      ::
        %unum
      ::  data is already a posit bit pattern; pack it raw, like %int2.
      (spac [meta `@ux`data])
      ::
        %cplx
      ::  data is already a packed complex bit pattern; pack it raw.
      (spac [meta `@ux`data])
      ::
        %fixp
      ::  data is already a fixed-point bit pattern; pack it raw.
      (spac [meta `@ux`data])
      ::
        %i754
      ::  convert date to fl to @r XXX TODO REVISIT whether we want to specify input type
      =/  fin
        ?+    bloq.meta  !!
          %7  (bit:ma:rq (sea:ma:rq data))
          %6  (bit:ma:rd (sea:ma:rd data))
          %5  (bit:ma:rs (sea:ma:rs data))
          %4  (bit:ma:rh (sea:ma:rh data))
        ==
      (spac [meta `@ux`fin])
    ==
  ::
  ::  Operators
  ::
  ::  +max/+min (and +argmax/+argmin, which call them) reduce with the kind's
  ::  %gth/%lth ordering, so they require a TOTALLY ORDERED kind.  On %cplx the
  ::  ordering op deliberately crashes with 'lagoon: %cplx has no total order;
  ::  use abs/equ' (the complex numbers have no total order).  +any/+all inherit
  ::  this, since they read the reduced +max/+min scalar.
  ++  max
    ~/  %max
    |=  a=ray
    ?>  (check a)
    =/  fun
      |:  [b=1 c=-:(ravel a)] 
      ?.  =(((fun-scalar meta.a %gth) b c) 0)
        b  c 
    (scalar-to-ray meta.a (reel (ravel a) fun))
  ::
  ++  argmax :: Only returns first match
    ~/  %argmax
    |=  a=ray
    ^-  @ud
    ?>  (check a)
    +:(find ~[(get-item (max a) ~[0 0])] (ravel a))
  ::
  ++  min
    ~/  %min
    |=  a=ray
    ?>  (check a)
    =/  fun
      |:  [b=1 c=-:(ravel a)] 
      ?.  =(((fun-scalar meta.a %lth) b c) 0)
        b  c 
    (scalar-to-ray meta.a (reel (ravel a) fun))
  ::
  ++  argmin :: Only returns first match
    ~/  %argmin
    |=  a=ray
    ^-  @ud
    ?>  (check a)
    +:(find ~[(get-item (min a) ~[0 0])] (ravel a))
  ::
  ::  NB: despite the name, this is a full reduction to a scalar (the total
  ::  sum), not a running/prefix cumulative sum.
  ++  cumsum
    ~/  %cumsum
    |=  a=ray
    ^-  ray
    ?>  (check a)
    ::  %unum: exact quire sum (single rounding), see +unum-sum.
    ?:  ?=(%unum kind.meta.a)
      (scalar-to-ray meta.a (unum-sum bloq.meta.a (ravel a)))
    %+  scalar-to-ray  meta.a
    (reel (ravel a) |=([b=@ c=@] ((fun-scalar meta.a %add) b c)))
  ::
  ++  prod
    |=  a=ray
    ^-  ray
    ?>  (check a)
    %+  scalar-to-ray  meta.a
    (reel (ravel a) |=([b=_1 c=_1] ((fun-scalar meta.a %mul) b c)))
  ::
  ++  reshape
    |=  [a=ray shape=(list @)]
    ^-  ray
    ?>  (check a)
    =/  in-cnt  (reel shape.meta.a ^mul)
    =/  out-cnt  (reel shape ^mul)
    ?>  =(in-cnt out-cnt)
    =.  shape.meta.a  shape
    a
  ::  stack along dimension (0 row, 1 col, 2 lay, etc.)
  ++  stack
    ~/  %stack
    |=  [a=ray b=ray dim=@ud]
    ^-  ray
    ?>  (check a)
    ?>  (check b)
    ::  check same dims overall
    ?>  =((lent shape.meta.a) (lent shape.meta.b))
    ::  check same dims other than target dim
    ?>  =/  idx  0
      |-  ^-  ?
      ?:  =(idx (lent shape.meta.a))  %.y
      ?:  =(dim idx)  $(idx +(idx))
      ?.  =((snag idx shape.meta.a) (snag idx shape.meta.b))
        %.n
      $(idx +(idx))
    ::  TODO revisit this assumption/requirement
    ?>  (^lte dim (lent shape.meta.a))
    =|  c=ray
    ?:  =(0 dim)
      =.  meta.c  meta.a
      =.  shape.meta.c
        :-  (^add (snag dim shape.meta.a) (snag dim shape.meta.b))
            +.shape.meta.a
      =.  data.c  (con (lsh [bloq.meta.a (roll shape.meta.a ^mul)] data.a) data:(unspac b))
      c
    ?:  =(1 dim)
      =.  meta.c  meta.a
      =.  shape.meta.c
        (snap shape.meta.c dim (^add (snag dim shape.meta.a) (snag dim shape.meta.b)))
      =/  c  (zeros meta.c)
      =/  idx  0
      |-  ^-  ray
      ?:  =((snag 0 (flop shape.meta.a)) idx)  c
      =/  off  (weld (snip (snip shape.meta.a)) ~[idx])
      =/  row-a  (get-row a off)
      =/  row-b  (get-row b off)
      =/  data-c  (con (lsh [bloq.meta.a (snag 1 shape.meta.row-a)] data.row-a) data:(unspac row-b))
      =/  meta-c=meta  meta.row-a
      =.  shape.meta-c  (snap shape.meta-c dim (^add (snag dim shape.meta.row-a) (snag dim shape.meta.row-b)))
      =/  row-c=ray  (spac [meta-c data-c])
      %=  $
        idx  +(idx)
        c    (set-row c off row-c)
      ==
    !!
  ::
  ++  hstack
    |=  [a=ray b=ray]
    ^-  ray
    (stack a b 1)
  ::
  ++  vstack
    |=  [a=ray b=ray]
    ^-  ray
    (stack a b 0)
  ::
  ::
  ++  transpose
    ~/  %transpose
    |=  a=ray  ^-  ray
    ?>  (check a)
    =,  meta.a
    ?>  =(2 (lent shape))
    =/  i  0
    =/  j  0
    =/  shape=(list @)  ~[(snag 1 shape) (snag 0 shape)]
    =/  prod=ray  (zeros [shape bloq kind ~])
    |-
      ?:  =(i (snag 0 shape.meta.a))
        prod
      %=  $
        i  +(i)
        prod
      |-
        ?:  =(j (snag 1 shape.meta.a))
          prod
        %=  $
          j  +(j)
          prod  (set-item prod ~[j i] (get-item a ~[i j]))
        ==
    ==
  ::  Returns diagonal of an array.
  ::
  ++  diag
    ~/  %diag
    |=  a=ray
    ^-  ray
    ?>  (check a)
    =,  meta.a
    ?>  =(2 (lent shape))
    ?>  =(-.shape +<.shape)
    ^-  ray
    %-  en-ray
    ^-  baum
    :-  `meta`[~[-.shape 1] bloq kind tail]
    ^-  ndray
    %+  turn
      `(list @)`(flop (gulf 0 (dec -.shape)))
    |=(i=@ (get-item a ~[i i]))
  ::
  ++  trace
    ~/  %trace
    |=  a=ray
    ^-  ray
    (cumsum (diag a))
  ::
  ++  dot
    ~/  %dot
    |=  [a=ray b=ray]
    ^-  ray
    ?>  =(shape.meta.a shape.meta.b)
    ::  %unum: fuse the multiply-accumulate in the quire (single rounding)
    ::  rather than rounding each product and each partial sum.
    ?:  ?=(%unum kind.meta.a)
      (scalar-to-ray meta.a (unum-fdp bloq.meta.a (ravel a) (ravel b)))
    ::  %fixp: accumulate the integer products exactly, rescale once.
    ?:  ?=(%fixp kind.meta.a)
      (scalar-to-ray meta.a (fixp-fdp ;;([a=@ b=@] tail.meta.a) (ravel a) (ravel b)))
    (cumsum (mul a b))
  ::    +dotc:  [$ray $ray] -> $ray
  ::
  ::  Returns the Hermitian (conjugated) dot product Sum conj(a_i)*b_i.  This is
  ::  the inner product whose norm <a,a> = Sum |a_i|^2 is real and nonnegative
  ::  (what Saloon's eig/QR want).  On real kinds conj is identity, so dotc
  ::  coincides with +dot.
  ::    Examples
  ::      > =a  (en-ray:la [[~[1 2] 6 %cplx ~] ~[~[0x4000.0000.3f80.0000 0x4080.0000.4040.0000]]])  ::  [1+2i 3+4i]
  ::      > =b  (en-ray:la [[~[1 2] 6 %cplx ~] ~[~[0x40c0.0000.40a0.0000 0x4100.0000.40e0.0000]]])  ::  [5+6i 7+8i]
  ::      > (get-item:la (dotc:la a b) ~[0 0])
  ::      0xc100.0000.428c.0000                                ::  Sum conj(x)*y = 70-8i
  ::      > (get-item:la (dot:la a b) ~[0 0])
  ::      0x4288.0000.c190.0000                                ::  Sum x*y       = -18+68i
  ::  Source
  ++  dotc
    ~/  %dotc
    |=  [a=ray b=ray]
    ^-  ray
    ?>  =(shape.meta.a shape.meta.b)
    (dot (conj a) b)
  ::
  ++  mmul
    ~/  %mmul
    |=  [a=ray b=ray]
    ?>  (check a)
    ?>  (check b)
    =/  i  0
    =/  j  0
    =/  k  0
    =/  shape=(list @)  ~[(snag 0 shape.meta.a) (snag 1 shape.meta.b)]
    =/  prod=ray  =,(meta.a (zeros [^shape bloq kind ~]))
    ::  
    ::  multiplication conditions
    ?>
    ?&  =(2 (lent shape.meta.b))
        =(2 (lent shape.meta.a))
        =((snag 1 shape.meta.a) (snag 0 shape.meta.b))
    ==
    ::  %unum: accumulate each output cell exactly in the quire, see +mmul-unum.
    ?:  ?=(%unum kind.meta.a)  (mmul-unum a b)
    ::  %fixp: accumulate each output cell's integer products exactly.
    ?:  ?=(%fixp kind.meta.a)  (mmul-fixp a b)
    |-
      ?:   =(i (snag 0 shape.meta.prod))
        prod
      %=    $
        i  +(i)
        prod
      |-
        ?:  =(j (snag 1 shape.meta.prod))
          prod
        =/  cume  0
        %=    $
            j  +(j)
            prod
          |-  
          ?:   =(k (snag 1 shape.meta.a))
            (set-item prod `(list @)`~[i j] cume)
          %=    $
              k  +(k)
              cume
            %+  (fun-scalar meta.a %add)
              cume
            %+  (fun-scalar meta.a %mul)
              (get-item a ~[i k])
            (get-item b ~[k j])
          ==
        ==
      ==
  ::  +mmul-unum: posit matrix multiply with quire-exact accumulation.  Each
  ::  output cell [i j] is the fused dot product of row i of a and column j of
  ::  b, rounded once.  Assumes the shape checks in +mmul already passed.
  ++  mmul-unum
    ~/  %mmul-unum
    |=  [a=ray b=ray]
    ^-  ray
    =/  m   (snag 0 shape.meta.a)
    =/  nn  (snag 1 shape.meta.b)
    =/  kk  (snag 1 shape.meta.a)
    =/  bl  bloq.meta.a
    =/  prod=ray  =,(meta.a (zeros [~[m nn] bloq kind ~]))
    =/  i  0
    |-  ^-  ray
    ?:  =(i m)  prod
    =/  j  0
    |-  ^-  ray
    ?:  =(j nn)  ^$(i +(i), prod prod)
    =/  av=(list @)  (turn (gulf 0 (dec kk)) |=(p=@ `@`(get-item a ~[i p])))
    =/  bv=(list @)  (turn (gulf 0 (dec kk)) |=(p=@ `@`(get-item b ~[p j])))
    $(j +(j), prod (set-item prod ~[i j] (unum-fdp bl av bv)))
  ::  +mmul-fixp: fixed-point matrix multiply, each cell an exact integer
  ::  dot product rescaled once (see +fixp-fdp), preserving the array
  ::  precision (meta.tail).  Assumes +mmul's shape checks already passed.
  ++  mmul-fixp
    ~/  %mmul-fixp
    |=  [a=ray b=ray]
    ^-  ray
    =/  prc  ;;([a=@ b=@] tail.meta.a)
    =/  m   (snag 0 shape.meta.a)
    =/  nn  (snag 1 shape.meta.b)
    =/  kk  (snag 1 shape.meta.a)
    =/  prod=ray  =,(meta.a (zeros [~[m nn] bloq kind tail]))
    =/  i  0
    |-  ^-  ray
    ?:  =(i m)  prod
    =/  j  0
    |-  ^-  ray
    ?:  =(j nn)  ^$(i +(i), prod prod)
    =/  av=(list @)  (turn (gulf 0 (dec kk)) |=(p=@ `@`(get-item a ~[i p])))
    =/  bv=(list @)  (turn (gulf 0 (dec kk)) |=(p=@ `@`(get-item b ~[p j])))
    $(j +(j), prod (set-item prod ~[i j] (fixp-fdp prc av bv)))
::
  ++  abs
    ~/  %abs
    |=  a=ray
    ^-  ray
    (el-wise-op a (trans-scalar bloq.meta.a kind.meta.a %abs))
::
  ::    +conj:  $ray -> $ray
  ::
  ::  Returns the elementwise complex conjugate of a ray.  Identity on real
  ::  kinds, which is what makes +dotc reduce to +dot on reals.
  ::    Examples
  ::      > (get-item:la (conj:la (fill:la [~[1 1] 6 %cplx ~] 0x4000.0000.3f80.0000)) ~[0 0])
  ::      0xc000.0000.3f80.0000                                ::  conj(1+2i)=1-2i
  ::  Source
  ++  conj
    ~/  %conj
    |=  a=ray
    ^-  ray
    (el-wise-op a (trans-scalar bloq.meta.a kind.meta.a %conj))
::
  ++  add-scalar
    ~/  %add-scal
    |=  [a=ray n=@]
    ^-  ray
    =/  b=ray  (fill meta.a n)
    (add a b)
  ::
  ++  sub-scalar
    ~/  %sub-scal
    |=  [a=ray n=@]
    ^-  ray
    =/  b=ray  (fill meta.a n)
    (sub a b)
  ::
  ++  mul-scalar
    ~/  %mul-scal
    |=  [a=ray n=@]
    ^-  ray
    =/  b=ray  (fill meta.a n)
    (mul a b)
  ::
  ++  div-scalar
    ~/  %div-scal
    |=  [a=ray n=@]
    ^-  ray
    =/  b=ray  (fill meta.a n)
    (div a b)
  ::
  ++  mod-scalar
    ~/  %mod-scal
    |=  [a=ray n=@]
    ^-  ray
    =/  b=ray  (fill meta.a n)
    (mod a b)
  ::
  ++  add
    ~/  %add-rays
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %add))
  ::
  ++  sub
    ~/  %sub-rays
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %sub))
  ::
  ++  mul
    ~/  %mul-rays
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %mul))
  ::
  ++  div
    ~/  %div-rays
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %div))
  ::
  ++  mod
    ~/  %mod-rays
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %mod))
  ::
  ++  pow-n
    |=  [a=ray n=@ud]
    ^-  ray
    ?>  (check a)
    ?~  =(0 n)  (ones meta.a)
    =/  b=ray   a
    |-  ^-  ray
    ?~  =(1 n)  b
    $(b (mul a b), n (dec n))
  ::
  ++  gth
    ~/  %gth
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %gth))
  ::
  ++  gte
    ~/  %gte
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %gte))
  ::
  ++  lth
    ~/  %lth
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %lth))
  ::
  ++  lte
    ~/  %lte
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %lte))
  ::
  ++  equ
    :: ~/  %equ
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %equ))
  ::
  ++  neq
    :: ~/  %equ
    |=  [a=ray b=ray]
    ^-  ray
    (bin-op a b (fun-scalar meta.a %neq))
  ::
  ++  mpow-n
    |=  [a=ray n=@ud]
    ^-  ray
    ?~  =(0 n)  (ones meta.a)
    =/  b=ray   a
    |-  ^-  ray
    ?~  =(1 n)  b
    $(b (mmul a b), n (dec n))
  ::  proximity check [abs rel]
  ::
  ++  is-close
    |=  [a=ray b=ray tol=[@ @]]
    ^-  ray
    ?>  =(shape.meta.a shape.meta.b)
    =/  atol  (fill meta.a data:(scale meta.a -.tol))
    =/  rtol  (fill meta.a data:(scale meta.a +.tol))
    (lte (abs (sub a b)) (add atol (mul rtol (abs b))))
  ::
  ::  Booleans use the numeric convention true=1, false=0 (see fun-scalar).
  ::  An element is truthy iff its encoded value is nonzero (also +0.0=0x0).
  ::  any: some element truthy  <=>  the max element is nonzero.
  ::  all: every element truthy  <=>  the min element is nonzero.
  ::  Read the reduced scalar as the head of its ravel, so this is
  ::  independent of min/max's output shape and works for any rank.
  ::  Because they reduce via +max/+min, +any/+all need a totally ordered kind
  ::  and crash on %cplx (see +max) -- test complex rays with +is-close/+equ.
  ++  any
    |=  [a=ray]
    ^-  ?(%.y %.n)
    ?!(=(-:(ravel (max a)) 0))
  ::
  ++  all
    |=  [a=ray]
    ^-  ?(%.y %.n)
    ?!(=(-:(ravel (min a)) 0))
  ::
  ::  Quire-exact posit reductions (/lib/unum fdp).  Sums of posit products
  ::  accumulate in the 16n-bit quire and round exactly once, so dot/mmul/
  ::  cumsum over %unum rays are correctly-rounded rather than rounding at
  ::  every intermediate step.  bloq selects width: 3=rpb 4=rph 5=rps 6=rpd.
  ::
  ++  unum-one
    |=  =bloq
    ^-  @
    ?+  bloq  !!
      %3  one:rpb:unum
      %4  one:rph:unum
      %5  one:rps:unum
      %6  one:rpd:unum
    ==
  ::  +unum-fdp: fused dot product of two equal-length posit vectors.
  ++  unum-fdp
    |=  [=bloq av=(list @) bv=(list @)]
    ^-  @
    ?+  bloq  !!
      %3  (fdp:rpb:unum av bv)
      %4  (fdp:rph:unum av bv)
      %5  (fdp:rps:unum av bv)
      %6  (fdp:rpd:unum av bv)
    ==
  ::  +unum-sum: exact sum of a posit vector (fdp against a vector of ones).
  ++  unum-sum
    |=  [=bloq v=(list @)]
    ^-  @
    (unum-fdp bloq v (reap (lent v) (unum-one bloq)))
  ::  +fixp-fdp: exact fixed-point dot product at precision [a b].  Decode each
  ::  operand to its stored signed integer, accumulate the integer products
  ::  exactly (no per-product truncation), then rescale by 2^b once and
  ::  re-encode at the element width N = a+b+1.
  ++  fixp-fdp
    |=  [prc=[a=@ b=@] av=(list @) bv=(list @)]
    ^-  @
    =/  sm=@s
      =|  sm=@s
      |-  ^-  @s
      ?~  av  sm
      ?~  bv  sm
      %=  $
        av  t.av
        bv  t.bv
        sm  (sum:si sm (pro:si (to-s:fixed i.av prc) (to-s:fixed i.bv prc)))
      ==
    (s-to-twoc:(ng:fixed prc) (fra:si sm (sun:si (bex b.prc))))
  ::  +unum-mod: posit remainder a - b*trunc(a/b), the quotient truncated
  ::  TOWARD ZERO (C fmod / remainder-with-sign-of-dividend).  Returns nar on
  ::  a NaR operand or division by zero (b = posit 0) rather than crashing.
  ++  unum-mod
    |=  [=bloq a=@ b=@]
    ^-  @
    ?+  bloq  !!
        %3
      =/  q  (div:rpb:unum a b)
      ?:  =(nar:rpb:unum q)  nar:rpb:unum
      =/  t  ?:((gte:rpb:unum q zero:rpb:unum) (flr:rpb:unum q) (cel:rpb:unum q))
      (sub:rpb:unum a (mul:rpb:unum b t))
        %4
      =/  q  (div:rph:unum a b)
      ?:  =(nar:rph:unum q)  nar:rph:unum
      =/  t  ?:((gte:rph:unum q zero:rph:unum) (flr:rph:unum q) (cel:rph:unum q))
      (sub:rph:unum a (mul:rph:unum b t))
        %5
      =/  q  (div:rps:unum a b)
      ?:  =(nar:rps:unum q)  nar:rps:unum
      =/  t  ?:((gte:rps:unum q zero:rps:unum) (flr:rps:unum q) (cel:rps:unum q))
      (sub:rps:unum a (mul:rps:unum b t))
        %6
      =/  q  (div:rpd:unum a b)
      ?:  =(nar:rpd:unum q)  nar:rpd:unum
      =/  t  ?:((gte:rpd:unum q zero:rpd:unum) (flr:rpd:unum q) (cel:rpd:unum q))
      (sub:rpd:unum a (mul:rpd:unum b t))
    ==
  ::  +unum-pow: a raised to an integer power.  The exponent posit b is taken
  ::  to its nearest integer (via toi); a NEGATIVE exponent yields the
  ::  reciprocal a^-n = 1/(a^n).  Returns nar if b is NaR.
  ++  unum-pow
    |=  [=bloq a=@ b=@]
    ^-  @
    ?+  bloq  !!
        %3
      =/  ui  (toi:rpb:unum b)
      ?~  ui  nar:rpb:unum
      =/  p  (pow-n:rpb:unum a `@u`(abs:si u.ui))
      ?:((lth:rpb:unum b zero:rpb:unum) (div:rpb:unum one:rpb:unum p) p)
        %4
      =/  ui  (toi:rph:unum b)
      ?~  ui  nar:rph:unum
      =/  p  (pow-n:rph:unum a `@u`(abs:si u.ui))
      ?:((lth:rph:unum b zero:rph:unum) (div:rph:unum one:rph:unum p) p)
        %5
      =/  ui  (toi:rps:unum b)
      ?~  ui  nar:rps:unum
      =/  p  (pow-n:rps:unum a `@u`(abs:si u.ui))
      ?:((lth:rps:unum b zero:rps:unum) (div:rps:unum one:rps:unum p) p)
        %6
      =/  ui  (toi:rpd:unum b)
      ?~  ui  nar:rpd:unum
      =/  p  (pow-n:rpd:unum a `@u`(abs:si u.ui))
      ?:((lth:rpd:unum b zero:rpd:unum) (div:rpd:unum one:rpd:unum p) p)
    ==
  ::
  +$  ops   $?  %add
                %sub
                %mul
                %div
                %mod
                %pow
                %gth
                %gte
                %lth
                %lte
                %equ
                %neq
                %abs
                %conj
            ==
  ::
  ++  fun-scalar
    |=  [=meta fun=ops]
    ^-  $-([@ @] @)
    =,  meta
    ?-    `^kind`kind
        %uint
      ?+  fun  !!
        %add  ~(sum fe bloq)
        %sub  ~(dif fe bloq)
        %mul  |=([b=_1 c=_1] (~(sit fe bloq) (^mul b c)))
        %div  |=([b=_1 c=_1] (~(sit fe bloq) (^div b c)))
        %mod  |=([b=@ c=@] (~(sit fe bloq) (^mod b c)))
        %pow  |=([b=@ c=@] (~(sit fe bloq) (^pow b c)))
        %gth  |=([b=@ c=@] !(^gth b c))
        %gte  |=([b=@ c=@] !(^gte b c))
        %lth  |=([b=@ c=@] !(^lth b c))
        %lte  |=([b=@ c=@] !(^lte b c))
        %equ  |=([b=@ c=@] ?:(.=(b c) 1 0))
        %neq  |=([b=@ c=@] ?:(.=(b c) 0 1))
      ==
      ::
        %int2
      ::  modular two's-complement (wraps on overflow), via /lib/twoc.
      ::  comparisons return 1=true / 0=false as a value, matching %uint.
      ?+  fun  !!
        %add  ~(add twoc:twoc bloq)
        %sub  ~(sub twoc:twoc bloq)
        %mul  ~(mul twoc:twoc bloq)
        %div  ~(div twoc:twoc bloq)
        %mod  ~(rem twoc:twoc bloq)
        %pow  ~(pow twoc:twoc bloq)
        %gth  |=([b=@ c=@] ?:((~(gth twoc:twoc bloq) b c) 1 0))
        %gte  |=([b=@ c=@] ?:((~(gte twoc:twoc bloq) b c) 1 0))
        %lth  |=([b=@ c=@] ?:((~(lth twoc:twoc bloq) b c) 1 0))
        %lte  |=([b=@ c=@] ?:((~(lte twoc:twoc bloq) b c) 1 0))
        %equ  |=([b=@ c=@] ?:(.=(b c) 1 0))
        %neq  |=([b=@ c=@] ?:(.=(b c) 0 1))
      ==
      ::
        %i754
      ::  %mod is C fmod: a - b*trunc(a/b), truncating the quotient TOWARD ZERO
      ::  (so the result takes the sign of the dividend), matching %unum.  A zero
      ::  divisor or a non-finite operand makes the quotient non-representable as
      ::  an integer (toi -> ~) and returns NaN rather than crashing.
      ::  NOTE (Vere follow-on): the f32/f64 %mod jet in vere's lagoon.c still
      ::  does round-nearest (IEEE remainder) and returns the dividend on /0, so
      ::  on a jetted ship f32/f64 disagree with this Hoon (and with the
      ::  un-jetted @rq path).  The jet must be updated to C fmod + NaN-on-/0.
      ?+    `^bloq`bloq  !!
          %7
        ?+  fun  !!
          %add  ~(add rq rnd)
          %sub  ~(sub rq rnd)
          %mul  ~(mul rq rnd)
          %div  ~(div rq rnd)
          %mod  |=([a=@rq b=@rq] ^-(@rq ?:((~(equ rq rnd) b .~~~0) `@rq`0x7fff.8000.0000.0000.0000.0000.0000.0000 =/(u (~(toi rq %z) (~(div rq rnd) a b)) ?~(u `@rq`0x7fff.8000.0000.0000.0000.0000.0000.0000 (~(sub rq rnd) a (~(mul rq rnd) b (~(san rq rnd) u.u))))))))
          %gth  |=([a=@rq b=@rq] ?:((~(gth rq rnd) a b) .~~~1 .~~~0))
          %gte  |=([a=@rq b=@rq] ?:((~(gte rq rnd) a b) .~~~1 .~~~0))
          %lth  |=([a=@rq b=@rq] ?:((~(lth rq rnd) a b) .~~~1 .~~~0))
          %lte  |=([a=@rq b=@rq] ?:((~(lte rq rnd) a b) .~~~1 .~~~0))
          %equ  |=([a=@rq b=@rq] ?:(.=(a b) .~~~1 .~~~0))
          %neq  |=([a=@rq b=@rq] ?:(.=(a b) .~~~0 .~~~1))
        ==
          %6
        ?+  fun  !!
          %add  ~(add rd rnd)
          %sub  ~(sub rd rnd)
          %mul  ~(mul rd rnd)
          %div  ~(div rd rnd)
          %mod  |=([a=@rd b=@rd] ^-(@rd ?:((~(equ rd rnd) b .~0) `@rd`0x7ff8.0000.0000.0000 =/(u (~(toi rd %z) (~(div rd rnd) a b)) ?~(u `@rd`0x7ff8.0000.0000.0000 (~(sub rd rnd) a (~(mul rd rnd) b (~(san rd rnd) u.u))))))))
          %gth  |=([a=@rd b=@rd] ?:((~(gth rd rnd) a b) .~1 .~0))
          %gte  |=([a=@rd b=@rd] ?:((~(gte rd rnd) a b) .~1 .~0))
          %lth  |=([a=@rd b=@rd] ?:((~(lth rd rnd) a b) .~1 .~0))
          %lte  |=([a=@rd b=@rd] ?:((~(lte rd rnd) a b) .~1 .~0))
          %equ  |=([a=@rd b=@rd] ?:(.=(a b) .~1 .~0))
          %neq  |=([a=@rd b=@rd] ?:(.=(a b) .~0 .~1))
        ==
          %5
        ?+  fun  !!
          %add  ~(add rs rnd)
          %sub  ~(sub rs rnd)
          %mul  ~(mul rs rnd)
          %div  ~(div rs rnd)
          %mod  |=([a=@rs b=@rs] ^-(@rs ?:((~(equ rs rnd) b .0) `@rs`0x7fc0.0000 =/(u (~(toi rs %z) (~(div rs rnd) a b)) ?~(u `@rs`0x7fc0.0000 (~(sub rs rnd) a (~(mul rs rnd) b (~(san rs rnd) u.u))))))))
          %gth  |=([a=@rs b=@rs] ?:((~(gth rs rnd) a b) .1 .0))
          %gte  |=([a=@rs b=@rs] ?:((~(gte rs rnd) a b) .1 .0))
          %lth  |=([a=@rs b=@rs] ?:((~(lth rs rnd) a b) .1 .0))
          %lte  |=([a=@rs b=@rs] ?:((~(lte rs rnd) a b) .1 .0))
          %equ  |=([a=@rs b=@rs] ?:(.=(a b) .1 .0))
          %neq  |=([a=@rs b=@rs] ?:(.=(a b) .0 .1))
        ==
          %4
        ?+  fun  !!
          %add  ~(add rh rnd)
          %sub  ~(sub rh rnd)
          %mul  ~(mul rh rnd)
          %div  ~(div rh rnd)
          %mod  |=([a=@rh b=@rh] ^-(@rh ?:((~(equ rh rnd) b .~~0) `@rh`0x7e00 =/(u (~(toi rh %z) (~(div rh rnd) a b)) ?~(u `@rh`0x7e00 (~(sub rh rnd) a (~(mul rh rnd) b (~(san rh rnd) u.u))))))))
          %gth  |=([a=@rh b=@rh] ?:((~(gth rh rnd) a b) .~~1 .~~0))
          %gte  |=([a=@rh b=@rh] ?:((~(gte rh rnd) a b) .~~1 .~~0))
          %lth  |=([a=@rh b=@rh] ?:((~(lth rh rnd) a b) .~~1 .~~0))
          %lte  |=([a=@rh b=@rh] ?:((~(lte rh rnd) a b) .~~1 .~~0))
          %equ  |=([a=@rh b=@rh] ?:(.=(a b) .~~1 .~~0))
          %neq  |=([a=@rh b=@rh] ?:(.=(a b) .~~0 .~~1))
        ==
      ==  :: bloq real
      ::
        %unum
      ::  posits (/lib/unum), bloq selects width: 3=rpb 4=rph 5=rps 6=rpd.
      ::  arithmetic is correctly-rounded per the 2022 standard; comparisons
      ::  return posit 1 / 0 (i.e. `one` / `zero`) to match the value
      ::  convention of the other kinds.  %mod is the truncated remainder
      ::  a - b*trunc(a/b) (see +unum-mod); %pow raises to an integer exponent,
      ::  negative -> reciprocal (see +unum-pow).  Both return nar on a NaR
      ::  operand or division by zero rather than crashing.
      ?+    `^bloq`bloq  !!
          %3
        ?+  fun  !!
          %add  add:rpb:unum
          %sub  sub:rpb:unum
          %mul  mul:rpb:unum
          %div  div:rpb:unum
          %mod  |=([a=@ b=@] (unum-mod %3 a b))
          %pow  |=([a=@ b=@] (unum-pow %3 a b))
          %gth  |=([a=@ b=@] ?:((gth:rpb:unum a b) one:rpb:unum zero:rpb:unum))
          %gte  |=([a=@ b=@] ?:((gte:rpb:unum a b) one:rpb:unum zero:rpb:unum))
          %lth  |=([a=@ b=@] ?:((lth:rpb:unum a b) one:rpb:unum zero:rpb:unum))
          %lte  |=([a=@ b=@] ?:((lte:rpb:unum a b) one:rpb:unum zero:rpb:unum))
          %equ  |=([a=@ b=@] ?:((equ:rpb:unum a b) one:rpb:unum zero:rpb:unum))
          %neq  |=([a=@ b=@] ?:((neq:rpb:unum a b) one:rpb:unum zero:rpb:unum))
        ==
          %4
        ?+  fun  !!
          %add  add:rph:unum
          %sub  sub:rph:unum
          %mul  mul:rph:unum
          %div  div:rph:unum
          %mod  |=([a=@ b=@] (unum-mod %4 a b))
          %pow  |=([a=@ b=@] (unum-pow %4 a b))
          %gth  |=([a=@ b=@] ?:((gth:rph:unum a b) one:rph:unum zero:rph:unum))
          %gte  |=([a=@ b=@] ?:((gte:rph:unum a b) one:rph:unum zero:rph:unum))
          %lth  |=([a=@ b=@] ?:((lth:rph:unum a b) one:rph:unum zero:rph:unum))
          %lte  |=([a=@ b=@] ?:((lte:rph:unum a b) one:rph:unum zero:rph:unum))
          %equ  |=([a=@ b=@] ?:((equ:rph:unum a b) one:rph:unum zero:rph:unum))
          %neq  |=([a=@ b=@] ?:((neq:rph:unum a b) one:rph:unum zero:rph:unum))
        ==
          %5
        ?+  fun  !!
          %add  add:rps:unum
          %sub  sub:rps:unum
          %mul  mul:rps:unum
          %div  div:rps:unum
          %mod  |=([a=@ b=@] (unum-mod %5 a b))
          %pow  |=([a=@ b=@] (unum-pow %5 a b))
          %gth  |=([a=@ b=@] ?:((gth:rps:unum a b) one:rps:unum zero:rps:unum))
          %gte  |=([a=@ b=@] ?:((gte:rps:unum a b) one:rps:unum zero:rps:unum))
          %lth  |=([a=@ b=@] ?:((lth:rps:unum a b) one:rps:unum zero:rps:unum))
          %lte  |=([a=@ b=@] ?:((lte:rps:unum a b) one:rps:unum zero:rps:unum))
          %equ  |=([a=@ b=@] ?:((equ:rps:unum a b) one:rps:unum zero:rps:unum))
          %neq  |=([a=@ b=@] ?:((neq:rps:unum a b) one:rps:unum zero:rps:unum))
        ==
          %6
        ?+  fun  !!
          %add  add:rpd:unum
          %sub  sub:rpd:unum
          %mul  mul:rpd:unum
          %div  div:rpd:unum
          %mod  |=([a=@ b=@] (unum-mod %6 a b))
          %pow  |=([a=@ b=@] (unum-pow %6 a b))
          %gth  |=([a=@ b=@] ?:((gth:rpd:unum a b) one:rpd:unum zero:rpd:unum))
          %gte  |=([a=@ b=@] ?:((gte:rpd:unum a b) one:rpd:unum zero:rpd:unum))
          %lth  |=([a=@ b=@] ?:((lth:rpd:unum a b) one:rpd:unum zero:rpd:unum))
          %lte  |=([a=@ b=@] ?:((lte:rpd:unum a b) one:rpd:unum zero:rpd:unum))
          %equ  |=([a=@ b=@] ?:((equ:rpd:unum a b) one:rpd:unum zero:rpd:unum))
          %neq  |=([a=@ b=@] ?:((neq:rpd:unum a b) one:rpd:unum zero:rpd:unum))
        ==
      ==  :: bloq unum
      ::
        %cplx
      ::  complex (/lib/complex), bloq selects width: 5=ch (2x@rh) 6=cs (2x@rs)
      ::  7=cd (2x@rd) 8=cq (2x@rq).  add/sub/mul/div are per-component float ops;
      ::  comparisons are equ/neq only (complex has no total order, so
      ::  gth/gte/lth/lte crash); %conj is unary (see +trans-scalar).  Reductions
      ::  use the generic round-each-step path (no exact complex accumulator).
      =/  ord  |=([a=@ b=@] ^-(@ ~|('lagoon: %cplx has no total order; use abs/equ' !!)))
      ?+    `^bloq`bloq  !!
          %5
        ?+  fun  !!
          %add  ~(add ch:complex rnd)
          %sub  ~(sub ch:complex rnd)
          %mul  ~(mul ch:complex rnd)
          %div  ~(div ch:complex rnd)
          %equ  |=([a=@ b=@] ?:((~(equ ch:complex rnd) a b) ~(one ch:complex rnd) `@`0))
          %neq  |=([a=@ b=@] ?:((~(neq ch:complex rnd) a b) ~(one ch:complex rnd) `@`0))
          %gth  ord
          %gte  ord
          %lth  ord
          %lte  ord
        ==
          %6
        ?+  fun  !!
          %add  ~(add cs:complex rnd)
          %sub  ~(sub cs:complex rnd)
          %mul  ~(mul cs:complex rnd)
          %div  ~(div cs:complex rnd)
          %equ  |=([a=@ b=@] ?:((~(equ cs:complex rnd) a b) ~(one cs:complex rnd) `@`0))
          %neq  |=([a=@ b=@] ?:((~(neq cs:complex rnd) a b) ~(one cs:complex rnd) `@`0))
          %gth  ord
          %gte  ord
          %lth  ord
          %lte  ord
        ==
          %7
        ?+  fun  !!
          %add  ~(add cd:complex rnd)
          %sub  ~(sub cd:complex rnd)
          %mul  ~(mul cd:complex rnd)
          %div  ~(div cd:complex rnd)
          %equ  |=([a=@ b=@] ?:((~(equ cd:complex rnd) a b) ~(one cd:complex rnd) `@`0))
          %neq  |=([a=@ b=@] ?:((~(neq cd:complex rnd) a b) ~(one cd:complex rnd) `@`0))
          %gth  ord
          %gte  ord
          %lth  ord
          %lte  ord
        ==
          %8
        ?+  fun  !!
          %add  ~(add cq:complex rnd)
          %sub  ~(sub cq:complex rnd)
          %mul  ~(mul cq:complex rnd)
          %div  ~(div cq:complex rnd)
          %equ  |=([a=@ b=@] ?:((~(equ cq:complex rnd) a b) ~(one cq:complex rnd) `@`0))
          %neq  |=([a=@ b=@] ?:((~(neq cq:complex rnd) a b) ~(one cq:complex rnd) `@`0))
          %gth  ord
          %gte  ord
          %lth  ord
          %lte  ord
        ==
      ==  :: bloq cplx
      ::
        %fixp
      ::  fixed-point Q a.b (/lib/fixed), prec [a b] from meta.tail, element
      ::  width N = 2^bloq.  add/sub/mod and comparisons need only the width
      ::  and reuse /lib/twoc (identical to %int2 on the stored integer);
      ::  mul/div rescale by 2^b through /lib/fixed.  Comparisons return the
      ::  fixed value 1.0 (= 2^b) for true, 0.0 for false.  Divide-by-zero
      ::  crashes (fixed-point has no NaN).
      =/  prc  ;;([a=@ b=@] tail)
      =/  one  (bex b.prc)
      ?+  fun  !!
        %add  ~(add twoc:twoc bloq)
        %sub  ~(sub twoc:twoc bloq)
        %mul  |=([x=@ y=@] =/(p (mul:fixed x prc y prc) (scale:fixed -.p +.p prc)))
        %div  |=([x=@ y=@] -:(div:fixed x prc y prc))
        %mod  ~(rem twoc:twoc bloq)
        %gth  |=([x=@ y=@] ?:((~(gth twoc:twoc bloq) x y) one 0x0))
        %gte  |=([x=@ y=@] ?:((~(gte twoc:twoc bloq) x y) one 0x0))
        %lth  |=([x=@ y=@] ?:((~(lth twoc:twoc bloq) x y) one 0x0))
        %lte  |=([x=@ y=@] ?:((~(lte twoc:twoc bloq) x y) one 0x0))
        %equ  |=([x=@ y=@] ?:(.=(x y) one 0x0))
        %neq  |=([x=@ y=@] ?:(.=(x y) 0x0 one))
      ==
    ==  :: kind
  ::
  ++  trans-scalar
    |=  [=bloq =kind fun=ops]
    ^-  $-(@ @)
    ::  %conj is the complex conjugate; for every non-complex (real) kind it
    ::  is the identity, which is what makes +dotc reduce to +dot on reals.
    ?-    kind
        %uint
      ?+  fun  !!
        %abs   |=(b=@ b)
        %conj  |=(b=@ b)
      ==
      ::
        %int2
      ?+  fun  !!
        %abs   ~(abs twoc:twoc bloq)
        %conj  |=(b=@ b)
      ==
      ::
        %i754
      ?+    bloq  !!
          %7
        ?+  fun  !!
          %abs  |=(b=@ ?:((~(gte rq rnd) b .~~~0) b (~(mul rq rnd) b .~~~-1)))
          %conj  |=(b=@ b)
        ==
          %6
        ?+  fun  !!
          %abs  |=(b=@ ?:((~(gte rd rnd) b .~0) b (~(mul rd rnd) b .~-1)))
          %conj  |=(b=@ b)
        ==
          %5
        ?+  fun  !!
          %abs  |=(b=@ ?:((~(gte rs rnd) b .0) b (~(mul rs rnd) b .-1)))
          %conj  |=(b=@ b)
        ==
          %4
        ?+  fun  !!
          %abs  |=(b=@ ?:((~(gte rh rnd) b .~~0) b (~(mul rh rnd) b .~~-1)))
          %conj  |=(b=@ b)
        ==
      ==
      ::
        %unum
      ?+    bloq  !!
          %6
        ?+  fun  !!
          %abs   abs:rpd:unum
          %conj  |=(b=@ b)
        ==
          %5
        ?+  fun  !!
          %abs   abs:rps:unum
          %conj  |=(b=@ b)
        ==
          %4
        ?+  fun  !!
          %abs   abs:rph:unum
          %conj  |=(b=@ b)
        ==
          %3
        ?+  fun  !!
          %abs   abs:rpb:unum
          %conj  |=(b=@ b)
        ==
      ==
      ::
        %cplx
      ?+    bloq  !!
          %8
        ?+  fun  !!
          %abs   ~(abs cq:complex rnd)
          %conj  ~(conj cq:complex rnd)
        ==
          %7
        ?+  fun  !!
          %abs   ~(abs cd:complex rnd)
          %conj  ~(conj cd:complex rnd)
        ==
          %6
        ?+  fun  !!
          %abs   ~(abs cs:complex rnd)
          %conj  ~(conj cs:complex rnd)
        ==
          %5
        ?+  fun  !!
          %abs   ~(abs ch:complex rnd)
          %conj  ~(conj ch:complex rnd)
        ==
      ==
      ::  fixed-point abs is two's-complement abs at the element width; conj is
      ::  the identity (fixed-point is real), so +dotc reduces to +dot.
        %fixp
      ?+  fun  !!
        %abs   ~(abs twoc:twoc bloq)
        %conj  |=(b=@ b)
      ==
    ==
  ::
  ++  el-wise-op
    |=  [a=ray fun=$-(@ @)]
    ^-  ray
    ?>  (check a)
    %-  spac
    :-  meta.a
    =/  ali  (ravel a)
    %+  rep  bloq.meta.a
    %+  turn
      ali
    |=(e=@ (fun e))
 :: 
  ++  bin-op
    |=  [a=ray b=ray op=$-([@ @] @)]
    ^-  ray
    ?>  =(meta.a meta.b)
    ?>  (check a)
    ?>  (check b)
    %-  spac
    :-  meta.a
    =/  ali  (ravel a)
    =/  bob  (ravel b)
    %+  rep  bloq.meta.a
    %+  turn
      (gulf 0 (dec (lent ali)))
    |=  i=@
    (op (snag i ali) (snag i bob))
  ::
  ++  ter-op
    |=  [a=ray b=ray c=ray op=$-([@ @ @] @)]
    ^-  ray
    ?>  =(meta.a meta.b)
    ?>  =(meta.c meta.b)
    ?>  (check a)
    ?>  (check b)
    ?>  (check c)
    %-  spac:la
    :-  meta.a
    =/  ali  (ravel:la a)
    =/  bob  (ravel:la b)
    =/  car  (ravel:la c)
    %+  rep  bloq.meta.a
    %+  turn
      (gulf 0 (dec (lent ali)))
    |=  i=@
    (op (snag i ali) (snag i bob) (snag i car))
  --
--
