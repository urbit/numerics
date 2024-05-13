#   MAchine leaRning in hOON (Maroon)

**WIP ~2024.5.2 implementing tinygrad operations**

- Our current objective is to implement `++forward` inference.
- The `++backward` arms need to deal with `grad_output` arguments.

##  Code

### Layer Zero (Opcodes)

- UnaryOps
	- [x] EXP2  = hook_overflow(math.inf, lambda x: math.exp(x*math.log(2))),
	- [x] LOG2  = lambda x: math.log2(x) if x > 0 else -math.inf if x == 0 else math.nan,
	- [x] CAST  = lambda self,arg: FlopCounter(self.shape, self.consume_flops(), self.mem)
	- [x] SIN   = math.sin
	- [x] SQRT  = lambda x: math.sqrt(x) if x >= 0 else math.nan,
	- [x] NEG   = lambda x: (not x) if isinstance(x, bool) else -x,
- BinaryOps
	- [x] ADD   = operator.add
	- [x] SUB   = operator.sub
	- [x] MUL   = operator.mul
	- [x] DIV   = lambda x,y: int(x/y) if isinstance(x, int) else (x/y if y != 0 else x*math.inf)
	- [x] MAX   = operator.max
	- [x] MOD   = lambda x,y: abs(int(x))%abs(int(y))*(1,-1)[x<0]
	- [x] CMPLT = operator.lt
	- [x] CMPEQ = operator.eq
	- [x] XOR   = operator.xor
- TernaryOps
	- [x] WHERE  = lambda x,y,z: y if x else z
	- [x] MULACC = lambda x,y,z,dtype: f"(({x}*{y})+{z})" // TRITON.py
- ReduceOps
	- [x] SUM   = auto()
	- [x] MAX   = auto()
- BufferOps
	- [ ] LOAD  = lambda arg: FlopCounter(arg.st.shape, 0, {arg.idx: arg.dtype.itemsize * arg.st.real_size()}) → `fill:la`
	- [ ] CONST = lambda arg: FlopCounter(arg.st.shape, 0, {})
	- [ ] STORE =  lambda self,arg: FlopCounter(arg.st.shape, self.consume_flops(), {**self.mem, arg.idx: arg.dtype.itemsize * arg.st.real_size()}),
- LoadOps
	- [x] EMPTY = auto() → `zeros:la`
	- [ ] CONST = auto() → `fill:la`
	- [ ] COPY  = auto() → `return array`
	- [ ] CONTIGUOUS = auto()
	- [ ] CUSTOM     = auto()
	- [ ] ASSIGN     = auto() 

### Layer One (Functions)

Functions generally implement a `++forward` arm and a `++backward` (derivative) arm.

- Structure Ops
  - [ ] Contiguous
  - [ ] ContiguousBackward
  - [ ] Cast
- Unary Ops
  - [x] Neg
  - [x] Reciprocal
  - [x] Sin
  - [x] Relu
  - [x] Log
  - [x] Exp
  - [x] Sqrt
  - [x] Sigmoid
- Binary Ops
  - [x] Less
  - [x] Eq
  - [x] Xor
  - [x] Add
  - [x] Sub
  - [x] Mul
  - [x] Div
- Ternary Ops
  - [x] Where
- Reduce Ops
  - [x] Sum
  - [x] Max
- Movement Ops
  - [x] Expand
  - [x] Reshape
  - [ ] Permute
  - [ ] Pad - figure out the NumPy logic
  - [ ] Shrink - 
  - [ ] Flip - write fn to reverse dir along axis

shape=~[2 3 4 5]
axis=~[1 3]
shape=~[2 *3* 4 *5*]

```py
def forward(self, x:LazyBuffer, axis:Tuple[int, ...]) -> LazyBuffer:
    self.arg = tuple([-1 if i in set(axis) else 1 for i in range(len(x.shape))])
    return x.stride(self.arg)

  def forward(self, x:LazyBuffer, arg:Tuple[Tuple[int, int], ...]) -> LazyBuffer:
    self.narg = tuple([(p[0], s+p[0]) for s,p in zip(x.shape, arg)])
    return x.pad(arg)

def pad(self, arg:Tuple[Optional[Tuple[sint, sint]], ...], value:float=0.0) -> Tensor:
    if all(x is None or x == (0,0) for x in arg):
		return self
    ret = F.Pad.apply(self, arg=(narg:=tuple(x if x is not None else (0,0) for x in arg)))
    return ret if 0 == value else ret + F.Pad.apply(Tensor.ones_like(self), arg=narg).where(0, value)
```