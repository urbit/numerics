  ::
::::  /lib/tinygrad
::
/-  ls=lagoon,
    ts=tinygrad
/+  *saloon
|%
::  +take: set +tg params
::
::    rnd: rounding mode
::    rtol: relative tolerance, use the correct bit width @r
::
++  take
  |=  [inrnd=rounding-mode inrtol=@r]
  %*(. tg rnd inrnd, rtol inrtol)
::
++  tg
  =+  [rnd=*rounding-mode rtol=`@r`0x1]
  ::
  ::  tinygrad opcodes
  ::
  |%
  ++  aux
    |%
    ++  const
      |=  [=meta:ls val=@]
      ^-  tensor:ts
      (fill:(lake rnd) meta val)
    ::
    ++  one
      |=  [=meta:ls]
      ^-  @
      ?-    kind.meta
          %uint
        1
        ::
          %real
        ?+    bloq.meta  !!
            %7  .~~~1
            %6  .~1
            %5  .1
            %4  .~~1
        ==  ::  bloq
      ==  ::  kind
    ::
    ++  pi-by-2
      |=  [=meta:ls]
      ^-  @
      ?+    kind.meta  !!
          %real
        ?+    bloq.meta  !!
            %7  ^~((~(div rq:math [rnd rtol]) pi:rq:math .~~~2))
            %6  ^~((~(div rd:math [rnd rtol]) pi:rd:math .~2))
            %5  ^~((~(div rs:math [rnd rtol]) pi:rs:math .2))
            %4  ^~((~(div rh:math [rnd rtol]) pi:rh:math .~~2))
        ==  ::  bloq
      ==  ::  kind
    ::
    ++  log2
      |=  [=meta:ls]
      ^-  @
      ?+  bloq.meta  !!
        %7  log2:rq:math
        %6  log2:rd:math
        %5  log2:rs:math
        %4  log2:rh:math
      ==
    ::
    ++  ilog2
      |=  [=meta:ls]
      ^-  @
      ?+  bloq.meta  !!
        %7  ^~((~(div rq:math [rnd rtol]) .~~~1 log2:rq:math))
        %6  ^~((~(div rd:math [rnd rtol]) .~1 log2:rd:math))
        %5  ^~((~(div rs:math [rnd rtol]) .1 log2:rs:math))
        %4  ^~((~(div rh:math [rnd rtol]) .~~1 log2:rh:math))
      ==
    --  ::  aux
  ::
  ++  ops
    |%
    ::
    ::  UnaryOps
    ::
    ::  EXP2  = hook_overflow(math.inf, lambda x: math.exp(x*math.log(2)))
    ++  exp2
      |=  a=tensor:ts
      ^-  tensor:ts
      (exp:(sake [rnd rtol]) (mul-scalar:la a (log2:aux meta.a)))
    ::  LOG2  = lambda x: math.log2(x) if x > 0 else -math.inf if x == 0 else math.nan,
    ++  log2
      |=  a=tensor:ts
      ^-  tensor:ts
      (log-2:(sake [rnd rtol]) a)
    ::  CAST
    ++  cast  el-wise-op:(lake rnd)
    ::  SIN   = math.sin
    ++  sin
      |=  a=tensor:ts
      ^-  tensor:ts
      (sin:(sake [rnd rtol]) a)
    ::  SQRT  = lambda x: math.sqrt(x) if x >= 0 else math.nan,
    ++  sqrt
      |=  a=tensor:ts
      ^-  tensor:ts
      (sqrt:(sake [rnd rtol]) a)
    ::  NEG   = lambda x: (not x) if isinstance(x, bool) else -x,
    ++  neg
      |=  a=tensor:ts
      ^-  tensor:ts
      (neg:(sake [rnd rtol]) a)
    ::
    :: BinaryOps
    ::
    ::  ADD   = operator.add
    ++  add  add:(sake [rnd rtol])
    ::  SUB   = operator.sub
    ++  sub  sub:(sake [rnd rtol])
    ::  MUL   = operator.mul
    ++  mul  mul:(sake [rnd rtol])
    ::  DIV   = lambda x,y: int(x/y) if isinstance(x, int) else (x/y if y != 0 else x*math.inf)
    ++  div  div:(sake [rnd rtol])
    ::  MAX   = operator.max
    ++  max
      |=  [a=tensor:ts b=tensor:ts]
      ^-  tensor:ts
      (bin-op:(lake rnd) a b (fun-scalar meta.a %max))
    ::  NOTE: MOD, CMPLT don't have to be implemented on vectors, just scalars
    ::  MOD   = lambda x,y: abs(int(x))%abs(int(y))*(1,-1)[x<0]
    ::    XX Lagoon is actually general modulo, check if an issue in reals
    ++  mod  mod:(lake rnd)
    ::  CMPLT = operator.lt
    ++  cmplt
      |=  [a=tensor:ts b=tensor:ts]
      ^-  tensor:ts
      (bin-op:(lake rnd) a b (fun-scalar meta.a %cmplt))
    ::  CMPEQ = operator.eq
    ++  cmpeq
      |=  [a=tensor:ts b=tensor:ts]
      ^-  tensor:ts
      (bin-op:(lake rnd) a b (fun-scalar meta.a %cmpeq))
    ::  XOR   = operator.xor
    ++  xor
      |=  [a=tensor:ts b=tensor:ts]
      ^-  tensor:ts
      ?>  =(meta.a meta.b)
      (spac:la meta.a (mix data:a data:b))
    ::
    :: TernaryOps
    ::
    ::  WHERE  = lambda x,y,z: y if x else z
    ++  where
      |=  [a=tensor:ts b=tensor:ts c=tensor:ts]
      ^-  tensor:ts
      |^
      (ter-op:(lake rnd) a b c where-op)
      ++  where-op
        |=  [a=@ b=@ c=@]
        ^-  @
        ?.  =(0 a)
          b
        c
      --
    ::  MULACC = lambda x,y,z,dtype: f"(({x}*{y})+{z})"
    ++  mulacc
      |=  [a=tensor:ts b=tensor:ts c=tensor:ts]
      ^-  tensor:ts
      (add:(sake [rnd rtol]) (mul:(sake [rnd rtol]) a b) c)
    ::
    :: ReduceOps
    ::
    ::  SUM   = auto()
    ++  sumred  cumsum:(lake rnd)
    ::  MAX   = auto()
    ++  maxred  max:(lake rnd)
    ::
    ::  AUX FNS
    ::
    ++  fun-scalar
      |=  [=meta fun=@tas]
      ^-  $-([@ @] @)
      =,  meta
      ?+    `^kind`kind  ~|(kind !!)
          %real
        ?+    `^bloq`bloq  !!
            %7
          ?+  fun  !!
            %cmplt  |=([a=@rq b=@rq] ?:((~(lth rq rnd) a b) 0x1 0x0))
            %cmpeq  |=([a=@rq b=@rq] ?:((~(equ rq rnd) a b) 0x1 0x0))
            %max    |=([a=@rq b=@rq] ?:((~(gth rq rnd) a b) a b))
          ==  ::  fn
            %6
          ?+  fun  !!
            %cmplt  |=([a=@rd b=@rd] ?:((~(lth rd rnd) a b) 0x1 0x0))
            %cmpeq  |=([a=@rd b=@rd] ?:((~(equ rd rnd) a b) 0x1 0x0))
            %max    |=([a=@rd b=@rd] ?:((~(gth rd rnd) a b) a b))
          ==  ::  fn
            %5
          ?+  fun  !!
            %cmplt  |=([a=@rs b=@rs] ?:((~(lth rs rnd) a b) 0x1 0x0))
            %cmpeq  |=([a=@rs b=@rs] ?:((~(equ rs rnd) a b) 0x1 0x0))
            %max    |=([a=@rs b=@rs] ?:((~(gth rs rnd) a b) a b))
          ==  ::  fn
            %4
          ?+  fun  !!
            %cmplt  |=([a=@rh b=@rh] ?:((~(lth rh rnd) a b) 0x1 0x0))
            %cmpeq  |=([a=@rh b=@rh] ?:((~(equ rh rnd) a b) 0x1 0x0))
            %max    |=([a=@rh b=@rh] ?:((~(gth rh rnd) a b) a b))
          ==  ::  fn
        ==  ::  bloq
      ==  ::  kind
    --  ::  ops
  
  ::
  ++  fns
    |%
    ::
    ::  tinygrad functions
    ::
    :: ++  contiguous
    :: ++  contiguous-backward
    :: ++  cast
    ++  neg
      |%
      ::  (checked)
      ++  forward   neg:ops
      ::  (checked)
      ++  backward  neg:ops
      --
    ::
    ++  reciprocal
      |%
      ::  1/x
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (div:ops (const:aux meta.a (one:aux meta.a)) a)
      :: -1/x^2
      ::  (checked)
      ++  backward  !!
      --
    ::
    ++  sin
      |%
      ::  (checked)
      ++  forward   sin:ops
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    ++  relu
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (where:ops (gth:(lake rnd) a (const:aux meta.a 0x0)) a (const:aux meta.a 0x0))
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    ++  log
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (log-2:(sake [rnd rtol]) (mul:ops a (const:aux meta.a (log2:aux meta.a))))
      ::  (checked)
      ++  backward  !!
      --
    ::
    ++  exp
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (mul:ops `tensor:ts`(const:aux meta.a (ilog2:aux meta.a)) `tensor:ts`(exp2:ops a))
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    ++  sqrt
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (sqrt:ops a)
      ::  (checked)
      ++  backward  !!
      --
    :: (checked against https://towardsdatascience.com/derivative-of-the-sigmoid-function-536880cf918e)
    ++  sigmoid
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        =/  one  (const:aux meta.a (one:aux meta.a))
        (div:ops one (add:ops one (exp2:ops (mul:ops (neg:ops (const:aux meta.a (ilog2:aux meta.a))) a))))
      ::  (checked)
      ++  backward  !!
      --
    ::
    :: Binary ops
    ::
    ++  less
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (cmplt:ops a b)
      ::  (checked)
      ++  backward  !!
      --
    ::
    ++  eq
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (cmpeq:ops a b)
      ::  (checked)
      ++  backward  !!
      --
    ::
    ++  xor
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (xor:ops a b)
      ::  (checked)
      ++  backward  !!
      --
    ::
    ++  add
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (add:ops a b)
      ::  (checked)
      ++  backward  !!
      --
    ::
    ++  sub
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (sub:ops a b)
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    ++  mul
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (mul:ops a b)
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    ++  div
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts]
        ^-  tensor:ts
        (div:ops a b)
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    :: Ternary ops
    ::
    ++  where
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts b=tensor:ts c=tensor:ts]
        ^-  tensor:ts
        (where:ops a b c)
      ::  (checked)
      ++  backward  !!
      --
    ::
    :: Reduce ops
    ::
    ++  sumred
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (sumred:ops a)
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    ++  maxred
      |%
      ::  (checked)
      ++  forward
        |=  a=tensor:ts
        ^-  tensor:ts
        (maxred:ops a)
      ::  (unchecked)
      ++  backward  !!
      --
    ::
    :: Movement ops
    ::
    ++  expand
      |%
      ++  forward
        |=  [a=tensor:ts shape=(list @)]
        ^-  tensor:ts
        |^
        =/  r  (zeros:(lake rnd) [shape bloq.meta.a kind.meta.a fxp.meta.a])
        ::  recursively set values into expanded tensor
        =/  ldx  (turn shape.meta.a dec)
        =/  dex  (reap (lent ldx) 0)
        |-
        ?:  =(ldx dex)  (set-item:(lake rnd) r dex (get-item:(lake rnd) a dex))
        $(dex (inc-index dex shape.meta.a), r (set-item:(lake rnd) r dex (get-item:(lake rnd) a dex)))
        ::  Using modular arithmetic, calculate the next index.
        ++  inc-index
          |=  [dex=(list @) str=(list @)]
          ^-  (list @)
          ?>  =((lent str) (lent dex))
          =/  i  0
          |-
          =/  j  (snag i dex)
          ?:  =(+(j) (snag i str))
            $(i +(i), dex (snap dex i 0))
          (snap dex i +(j))
        --
      ++  backward  !!
      --
    ++  reshape
      |%
      ::  (checked)
      ++  forward
        |=  [a=tensor:ts shape=(list @)]
        ^-  tensor:ts
        =.  shape.meta.a  shape
        a
      ::  (checked)
      ++  backward  !!
      --
    ++  permute  !!
    ++  shrink  !!
    ++  flip  !!
    --  ::  fns
  --  ::  tg
--
