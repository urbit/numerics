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
++  const
  |=  [=meta:ls val=@]
  ^-  ray:ls
  (fill:(lake rnd) meta val)
::
++  one
  |=  [=meta:ls]
  ^-  ray:ls
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
  ^-  ray:ls
  ?+    kind.meta  !!
      %real
    ?+    bloq.meta  !!
        %7  ^~((~(div rq:math rnd) pi:rq:math .~~~2))
        %6  ^~((~(div rd:math rnd) pi:rd:math .~2))
        %5  ^~((~(div rs:math rnd) pi:rs:math .2))
        %4  ^~((~(div rh:math rnd) pi:rh:math .~~2))
    ==  ::  bloq
  ==  ::  kind
::
++  log2
  |=  [=meta:ls]
  ^-  ray:ls
  =/  log2
    ?+  bloq.meta  !!
      %7  log2:rq:math
      %6  log2:rd:math
      %5  log2:rs:math
      %4  log2:rh:math
    ==
  ==  ::  bloq
::
++  ilog2
  |=  [=meta:ls]
  ^-  ray:ls
  =/  log2
    ?+  bloq.meta  !!
      %7  ^~((~(div rq:math rnd) .~~~1 log2:rq:math)
      %6  ^~((~(div rd:math rnd) .~1 log2:rd:math)
      %5  ^~((~(div rs:math rnd) .1 log2:rs:math)
      %4  ^~((~(div rh:math rnd) .~~1 log2:rh:math)
    ==
  ==  ::  bloq
::
++  tg
  =+  [rnd=*rounding-mode rtol=`@r`0x1]
  ::
  ::  tinygrad opcodes
  ::
  |%
  ++  ops
    |%
    ::
    ::  UnaryOps
    ::
    ::  EXP2  = hook_overflow(math.inf, lambda x: math.exp(x*math.log(2)))
    ++  exp2
      |=  a=ray:ls
      ^-  ray
      =/  log2
        ?+  bloq.meta.a  !!
          %7  log2:rq:math
          %6  log2:rd:math
          %5  log2:rs:math
          %4  log2:rh:math
        ==
      (exp:(sake [rnd rtol]) (mul-scalar:la a log2))
    ::  LOG2  = lambda x: math.log2(x) if x > 0 else -math.inf if x == 0 else math.nan,
    ++  log2
      |=  a=ray:ls
      ^-  ray
      (log-2:(sake [rnd rtol]) a)
    ::  CAST
    ++  cast  el-wise-op:(lake rnd)
    ::  SIN   = math.sin
    ++  sin
      |=  a=ray:ls
      ^-  ray
      (sin:(sake [rnd rtol]) a)
    ::  SQRT  = lambda x: math.sqrt(x) if x >= 0 else math.nan,
    ++  sqrt
      |=  a=ray:ls
      ^-  ray
      (sqrt:(sake [rnd rtol]) a)
    ::  NEG   = lambda x: (not x) if isinstance(x, bool) else -x,
    ++  neg
      |=  a=ray:ls
      ^-  ray
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
      |=  [a=ray:ls b=ray:ls]
      ^-  ray
      (bin-op:(lake rnd) a b (fun-scalar meta.a %max))
    ::  NOTE: MOD, CMPLT don't have to be implemented on vectors, just scalars
    ::  MOD   = lambda x,y: abs(int(x))%abs(int(y))*(1,-1)[x<0]
    ::    XX Lagoon is actually general modulo, check if an issue in reals
    ++  mod  mod:(lake rnd)
    ::  CMPLT = operator.lt
    ++  cmplt
      |=  [a=ray:ls b=ray:ls]
      ^-  ray
      (bin-op:(lake rnd) a b (fun-scalar meta.a %cmplt))
    ::  CMPEQ = operator.eq
    ++  cmpeq
      |=  [a=ray:ls b=ray:ls]
      ^-  ray
      (bin-op:(lake rnd) a b (fun-scalar meta.a %cmpeq))
    ::  XOR   = operator.xor
    ++  xor
      |=  [a=ray:ls b=ray:ls]
      ^-  ray
      ?>  =(meta.a meta.b)
      (spac:la meta.a (mix data:a data:b))
    ::
    :: TernaryOps
    ::
    ::  WHERE  = lambda x,y,z: y if x else z
    ++  where
      |=  [a=ray:ls b=ray:ls c=ray:ls]
      ^-  ray
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
      |=  [a=ray:ls b=ray:ls c=ray:ls]
      ^-  ray
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
    --
  ::
  ::  tinygrad functions
  ::
  +|  %functions
  :: ++  contiguous
  :: ++  contiguous-backward
  :: ++  cast
  ++  neg
    |%
    ++  forward   neg.ops
    ++  backward  neg.ops
    --
  ::
  ++  reciprocal
    |%
    ::  1/x
    ++  forward
      |=  a=ray
      ^-  ray
      (div:ops (const meta.a (one meta)) a)
    :: -1/x^2
    ++  backward
      |=  a=ray
      ^-  ray
      (div:ops (neg:ops (const meta.a (one meta))) (mul:ops a a))
    --
  ::
  ++  sin
    |%
    ::  sin(x)
    ++  forward   sin.ops
    ::  cos(x) = -sin(x - pi/2)
    ++  backward
      |=  a=ray
      ^-  ray
      (neg:ops (sin:ops (sub:ops a (const meta.a (pi-by-2 meta.a)))))
    --
  ::
  ++  relu
    |%
    ++  forward
      |=  a=ray
      ^-  ray
      (max:(lake rnd) a (const meta.a 0x0))
    :: return self.ret.const(0).e(BinaryOps.CMPLT, self.ret).cast(grad_output.dtype).e(BinaryOps.MUL, grad_output)
    ++  backward  !!
    --
  ::
  ++  log
    |%
    :: x.e(UnaryOps.LOG2).e(BinaryOps.MUL, x.const(math.log(2)))
    ++  forward
      |=  a=ray
      ^-  ray
      (log-2:(sake [rnd rtol]) (mul:ops a (const meta.a (log2 meta.a))))
    ++  backward
      |=  a=ray
      ^-  ray
      (div:ops (const meta.a (one meta)) a)
    --
  ::
  ++  exp
    |%
    :: x.e(BinaryOps.MUL, x.const(1/math.log(2))).e(UnaryOps.EXP2)
    ++  forward
      |=  a=ray
      ^-  ray
      (exp2:ops (mul:ops a (const meta.a (ilog2 meta.a))))
    ++  backward
      |=  a=ray
      ^-  ray
      (mul:ops (exp2:ops a) a)
    --
  ::
  ++  sqrt
    |%
    ++  forward
      |=  a=ray
      ^-  ray
      (sqrt:ops a)
    ++  backward
      |=  a=ray
      ^-  ray
      (div:ops a (mul:ops (const meta.a 2) (sqrt:ops a)))
    --
  ::
  ++  sigmoid
    |%
    ++  forward
      |=  a=ray
      ^-  ray
      (div:ops (const meta.a (one meta)) (add:ops (const meta.a (one meta)) (exp2:ops (mul:ops (neg:ops (const meta.a (ilog2 meta.a)))))))))))
    ++  backward
      |=  a=ray
      ^-  ray
      (mul:ops a (sub:ops (const meta.a 1) a))
    --
  ::
  :: Binary ops
  ::
  ++  less
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (cmplt:ops a b)
    ++  backward  !!
    --
  ::
  ++  eq
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (cmpeq:ops a b)
    ++  backward  !!
    --
  ::
  ++  xor
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (xor:ops a b)
    ++  backward  !!
    --
  ::
  ++  add
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (add:ops a b)
    ++  backward
      |=  [a=ray b=ray]
      ^-  ray
      [a b]
    --
  ::
  ++  sub
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (sub:ops a b)
    ++  backward
      |=  [a=ray b=ray]
      ^-  ray
      [a (neg:ops b)]
    --
  ::
  ++  mul
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (mul:ops a b)
    ++  backward
      |=  [a=ray b=ray]
      ^-  ray
      [b a]
    --
  ::
  ++  div
    |%
    ++  forward
      |=  [a=ray b=ray]
      ^-  ray
      (div:ops a b)
    ++  backward
      |=  [a=ray b=ray]
      ^-  ray
      :-  (div:ops a b)
      :: XXX ???
      (neg:ops (mul:ops a (div:ops y (mul:ops x y))))
    --
  ::
  :: Ternary ops
  ::
  ++  where
    |%
    ++  forward
      |=  [a=ray b=ray c=ray]
      ^-  ray
      (where:ops a b c)
    ++  backward  !!
    --
  ::
  :: Reduce ops
  ::
  ++  sumred
    |%
    ++  forward
      |=  a=ray
      ^-  ray
      (sumred:ops a)
    ++  backward
      |=  a=ray
      ^-  ray
      ::  XXX ??? not correct yet, just placeholder
      (expand a)
    --
  ::
  ++  maxred
    |%
    ++  forward
      |=  a=ray
      ^-  ray
      (maxred:ops a)
    ++  backward  !!
    --
  ::
  :: Movement ops
  ::
  ++  expand  !!
  ++  reshape  !!
  ++  permute  !!
  ++  shrink  !!
  ++  flip  !!
--
