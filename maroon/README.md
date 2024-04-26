#   MAchine leaRning in hOON (Maroon)

**WIP ~2024.4.26 implementing tinygrad operations**

```
# these are the llops your accelerator must implement, along with toCpu
# the Enum class doesn't work with mypy, this is static. sorry it's ugly
# NOTE: MOD, CMPLT don't have to be implemented on vectors, just scalars
# NOTE: many GPUs don't have DIV, but UnaryOps.RECIP doesn't work for integer division
- UnaryOps
	EXP2  = hook_overflow(math.inf, lambda x: math.exp(x*math.log(2))),
	LOG2  = lambda x: math.log2(x) if x > 0 else -math.inf if x == 0 else math.nan,
	CAST  = lambda self,arg: FlopCounter(self.shape, self.consume_flops(), self.mem)
	SIN   = math.sin
	SQRT  = lambda x: math.sqrt(x) if x >= 0 else math.nan, UnaryOps.SIN: math.sin,
	NEG   = lambda x: (not x) if isinstance(x, bool) else -x, 
- BinaryOps
	ADD   = operator.add
	SUB   = operator.sub
	MUL   = operator.mul
	DIV   = lambda x,y: int(x/y) if isinstance(x, int) else (x/y if y != 0 else x*math.inf)
	MAX   = operator.max
	MOD   = lambda x,y: abs(int(x))%abs(int(y))*(1,-1)[x<0]
	CMPLT = operator.lt
	CMPEQ = operator.eq
	XOR   = operator.xor 
- TernaryOps
	WHERE  = lambda x,y,z: y if x else z
	MULACC = lambda x,y,z,dtype: f"(({x}*{y})+{z})" // TRITON.py
- ReduceOps
	SUM   = auto()
	MAX   = auto()
- BufferOps
	LOAD  = lambda arg: FlopCounter(arg.st.shape, 0, {arg.idx: arg.dtype.itemsize * arg.st.real_size()}),
	CONST = lambda arg: FlopCounter(arg.st.shape, 0, {})
	STORE =  lambda self,arg: FlopCounter(arg.st.shape, self.consume_flops(), {**self.mem, arg.idx: arg.dtype.itemsize * arg.st.real_size()}),
- LoadOps
	EMPTY = auto()
	CONST = auto()
	COPY  = auto()
	CONTIGUOUS = auto()
	CUSTOM     = auto()
	ASSIGN     = auto() 
	
Op = Union[UnaryOps, BinaryOps, ReduceOps, LoadOps, TernaryOps, BufferOps]
OpType = Union[Type[UnaryOps], Type[BinaryOps], Type[ReduceOps], Type[LoadOps], Type[TernaryOps], Type[BufferOps]]
```

Operators:

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
	- [L] EMPTY = auto() → `zeros:la`
	- [ ] CONST = auto() → `fill:la`
	- [ ] COPY  = auto() → `return array`
	- [ ] CONTIGUOUS = auto()
	- [ ] CUSTOM     = auto()
	- [ ] ASSIGN     = auto() 

To implement:
- Lagoon: fix ++de-ray on leading zeros
- Saloon: adopt /lib/math
- Saloon: compatibility w/ Lagoon `+$ray`s
- Maroon:



```
;;((list (list @rs)) data:(de-ray:la:la (exp2:tg:tg (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.1 .2] ~[.3 .4] ~[.5 .6]]))))
;;((list (list @rs)) data:(de-ray:la:la (exp2:(tgrm %u .1e-5):tg (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.1 .2] ~[.3 .4] ~[.5 .6]]))))

(log2:(tgrm %u .1e-5):tg (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.1 .2] ~[.3 .4] ~[.5 .6]]))

(log2:(tgrm %u .1e-5):tg (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.1 .2] ~[.3 .4] ~[.5 .6]]))

(where:tg (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.1 .0] ~[.0 .0] ~[.0 .1]]) (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.2 .2] ~[.2 .2] ~[.2 .2]]) (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.3 .3] ~[.3 .3] ~[.3 .3]]))

(mulacc:tg (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.1 .0] ~[.0 .0] ~[.0 .1]]) (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.2 .2] ~[.2 .2] ~[.2 .2]]) (en-ray:la:la [~[3 2] 5 %real ~] ~[~[.3 .3] ~[.3 .3] ~[.3 .3]]))

```

