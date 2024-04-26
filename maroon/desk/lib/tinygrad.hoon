  ::
::::  /lib/tinygrad
::
/-  ls=lagoon
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
  ++  add  add:(lake rnd)
  ::  SUB   = operator.sub
  ++  sub  sub:(lake rnd)
  ::  MUL   = operator.mul
  ++  mul  mul:(lake rnd)
  ::  DIV   = lambda x,y: int(x/y) if isinstance(x, int) else (x/y if y != 0 else x*math.inf)
  ++  div  div:(lake rnd)
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
    (add:(lake rnd) (mul:(lake rnd) a b) c)
  ::
  :: ReduceOps
  ::
  ::  SUM   = auto()
  ++  sumred  cumsum:(lake rnd)
  ::  MAX   = auto()
  ++  maxred  max:(lake rnd)
  ::
  ::  MovementOps
  ::


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
--
