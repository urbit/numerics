  ::
::::  /lib/tinygrad
::
/-  ls=lagoon
/+  *lagoon, math
|%
++  tgrm
  |=  [inrnd=rounding-mode inrtol=@r]
  %*(. tg rnd inrnd, rtol inrtol)
::
++  tg
  =+  [rnd=*rounding-mode rtol=`@r`0x1 flops=(list @)]
  |%
  ++  hc
    |%
    ++  exp
      |=  a=ray:ls
      ^-  ray
      (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %exp))
    ++  log-2
      |=  a=ray:ls
      ^-  ray
      (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %log-2))
    ++  sin
      |=  a=ray:ls
      ^-  ray
      (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %sin))
    ++  sqrt
      |=  a=ray:ls
      ^-  ray
      (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %sqrt))
    ++  neg
      |=  a=ray:ls
      ^-  ray
      (el-wise-op:la a (trans-scalar bloq.meta.a kind.meta.a %neg))
    ++  trans-scalar
      |=  [=bloq =kind fun=@ta]
      ^-  $-(@ @)
      ?+    kind  ~|(kind !!)
          %real
        ?+    bloq  !!
            %7
          ?+  fun  !!
            %exp  ~(exp rq:math [rnd rtol])
            %log-2  ~(log-2 rq:math [rnd rtol])
            %neg  ~(neg rq:math [rnd rtol])
            %sin  ~(sin rq:math [rnd rtol])
            %sqrt  ~(sqrt rq:math [rnd rtol])
          ==
            %6
          ?+  fun  !!
            %exp  ~(exp rd:math [rnd rtol])
            %log-2  ~(log-2 rd:math [rnd rtol])
            %neg  ~(neg rd:math [rnd rtol])
            %sin  ~(sin rd:math [rnd rtol])
            %sqrt  ~(sqrt rd:math [rnd rtol])
          ==
            %5
          ?+  fun  !!
            %exp  ~(exp rs:math [rnd rtol])
            %log-2  ~(log-2 rs:math [rnd rtol])
            %neg  ~(neg rs:math [rnd rtol])
            %sin  ~(sin rs:math [rnd rtol])
            %sqrt  ~(sqrt rs:math [rnd rtol])
          ==
            %4
          ?+  fun  !!
            %exp  ~(exp rh:math [rnd rtol])
            %log-2  ~(log-2 rh:math [rnd rtol])
            %neg  ~(neg rh:math [rnd rtol])
            %sin  ~(sin rh:math [rnd rtol])
            %sqrt  ~(sqrt rh:math [rnd rtol])
          ==
        ==
      ==
    --
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
    (exp:hc (mul-scalar:la a log2))
  ::  LOG2  = lambda x: math.log2(x) if x > 0 else -math.inf if x == 0 else math.nan,
  ++  log2
    |=  a=ray:ls
    ^-  ray
    (log-2:hc a)
  ::  CAST
  ++  cast  el-wise-op:(lake rnd:tg)
  ::  SIN   = math.sin
  ++  sin
    |=  a=ray:ls
    ^-  ray
    (sin:hc a)
  ::  SQRT  = lambda x: math.sqrt(x) if x >= 0 else math.nan,
  ++  sqrt
    |=  a=ray:ls
    ^-  ray
    (sqrt:hc a)
  ::  NEG   = lambda x: (not x) if isinstance(x, bool) else -x,
  ++  neg
    |=  a=ray:ls
    ^-  ray
    (neg:hc a)
  --
  ::
  :: BinaryOps
  ::
  ::  ADD   = operator.add
  ++  add  add:(lake rnd:tg)
  ::  SUB   = operator.sub
  ++  sub  sub:(lake rnd:tg)
  ::  MUL   = operator.mul
  ++  mul  mul:(lake rnd:tg)
  ::  DIV   = lambda x,y: int(x/y) if isinstance(x, int) else (x/y if y != 0 else x*math.inf)
  ++  div  div:(lake rnd:tg)
  ::  MAX   = operator.max
  ++  max
    |=  [a=ray:ls b=ray:ls]
    ^-  ray
    (bin-op:(lake rnd:tg) a b (fun-scalar meta.a %max))
  ::  NOTE: MOD, CMPLT don't have to be implemented on vectors, just scalars
  ::  MOD   = lambda x,y: abs(int(x))%abs(int(y))*(1,-1)[x<0]
  ::    XX Lagoon is actually general modulo, check if an issue in reals
  ++  mod  mod:(lake rnd:tg)
  ::  CMPLT = operator.lt
  ++  cmplt
    |=  [a=ray:ls b=ray:ls]
    ^-  ray
    (bin-op:(lake rnd:tg) a b (fun-scalar meta.a %cmplt))
  ::  CMPEQ = operator.eq
  ++  cmpeq
    |=  [a=ray:ls b=ray:ls]
    ^-  ray
    (bin-op:(lake rnd:tg) a b (fun-scalar meta.a %cmpeq))
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
    (ter-op a b c where-op)
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
    (add:(lake rnd:tg) (mul:(lake rnd:tg) a b) c)
  ::
  :: ReduceOps
  ::
  ::  SUM   = auto()
  ++  sumred  cumsum:(lake rnd:tg)
  ::  MAX   = auto()
  ++  maxred  max:(lake rnd:tg)
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
        %cmplt  |=([a=@rq b=@rq] ?:((~(lth rq rnd:tg) a b) 0x1 0x0))
        %cmpeq  |=([a=@rq b=@rq] ?:((~(equ rq rnd:tg) a b) 0x1 0x0))
        %max    |=([a=@rq b=@rq] ?:((~(gth rq rnd:tg) a b) a b))
      ==  ::  fn
        %6
      ?+  fun  !!
        %cmplt  |=([a=@rd b=@rd] ?:((~(lth rd rnd:tg) a b) 0x1 0x0))
        %cmpeq  |=([a=@rd b=@rd] ?:((~(equ rd rnd:tg) a b) 0x1 0x0))
        %max    |=([a=@rd b=@rd] ?:((~(gth rd rnd:tg) a b) a b))
      ==  ::  fn
        %5
      ?+  fun  !!
        %cmplt  |=([a=@rs b=@rs] ?:((~(lth rs rnd:tg) a b) 0x1 0x0))
        %cmpeq  |=([a=@rs b=@rs] ?:((~(equ rs rnd:tg) a b) 0x1 0x0))
        %max    |=([a=@rs b=@rs] ?:((~(gth rs rnd:tg) a b) a b))
      ==  ::  fn
        %4
      ?+  fun  !!
        %cmplt  |=([a=@rh b=@rh] ?:((~(lth rh rnd:tg) a b) 0x1 0x0))
        %cmpeq  |=([a=@rh b=@rh] ?:((~(equ rh rnd:tg) a b) 0x1 0x0))
        %max    |=([a=@rh b=@rh] ?:((~(gth rh rnd:tg) a b) a b))
      ==  ::  fn
    ==  ::  bloq
  ==  ::  kind
  ::
  ++  ter-op
    |=  [a=ray b=ray c=ray op=$-([@ @ @] @)]
    ^-  ray
    ?>  =(meta.a meta.b)
    ?>  =(meta.c meta.b)
    %-  spac:la
    :-  meta.a
    =/  ali  (ravel:la a)
    =/  bob  (ravel:la b)
    =/  car  (ravel:la c)
    %^  rev  bloq.meta.a  (lent ali)
    %+  rep  bloq.meta.a
    %+  turn
      (gulf 0 (dec (lent ali)))
    |=  i=@
    (op (snag i ali) (snag i bob) (snag i car))

--
