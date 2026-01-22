#   `/lib/math` for Urbit

An unjetted library for basic mathematical and special function support.

We support the following functions and special functions:

- `++add`, $+$ addition
- `++sub`, $-$ subtraction
- `++mul`, $\times$ multiplication
- `++div`, $/$ division
- `++mod`, modulo (remainder after division)
- `++fma`, $\text{fma}$ fused multiply-add
- `++sgn`, $\text{sgn}$ signum (also `++sig`)
- `++neg`, $-$ unary negation
- `++factorial`, $!$ factorial
- `++abs`, $\text{abs}$
- `++exp`, $\exp$
- `++expm1`, $\exp - 1$
- `++sin`, $\sin$
- `++cos`, $\cos$
- `++tan`, $\tan$
- `++pow-n`, $\text{pow}$ to integer power
- `++log`, $\log$ (natural logarithm)
- `++log-10`, $\log_{10}$ (log base-10)
- `++log-2`, $\log_{2}$ (log base-2)
- `++pow`, $\text{pow}$
- `++sqrt`, $\sqrt{}$ (also `++sqt`)
- `++cbrt`, $\sqrt[3]{}$ (also `++cbt`)
- `++arg` (alias for `++abs` in real values)

Logical functions:

- `++lth`, $<$
- `++lte` $\leq$ (also `++leq`)
- `++gth`, $>$
- `++gte`, $\geq$ (also `++geq`)
- `++equ`, $=$
- `++neq`, $\neq$
- `++is-close`
- `++all-close`
- `++is-int`

Constants (to machine accuracy):

- `++tau`, $\tau = 2 \pi$
- `++pi`, $\pi$
- `++e`, $e$ (base of natural logarithm)
- `++phi`, $\phi$, Euler's constant
- `++sqt2`, $\sqrt{2}$
- `++invsqt2`, $\frac{1}{\sqrt{2}}$
- `++log2`, $\log 2$
- `++invlog2`, $\frac{1}{\log 2}$
- `++log10`, $\log 10$
- `++huge`, largest valid number in `bloq` width
- `++tiny`, smallest valid number in `bloq` size

---

It would be nice to have the following special functions as well:

- $\arcsin$
- $\arccos$
- $\arctan$
- $\sinh$
- $\cosh$
- $\tanh$

We do not envision including the Bessel functions and other more abstruse functions.

We use naïve algorithms which are highly reproducible.  We special-case some arguments to make them tractable without catastrophic cancellation.

---

# `/lib/unum` for Urbit

A library for universal numbers (posits, quires, and valids) per the Gustafson (2019) formulation and the 2022 Posit Standard.

- $\text{posit}<8,0>$
- $\text{posit}<16,1>$
- $\text{posit}<32,2>$

##  Posit Standard Compliance

Two interfaces are maintained for `/lib/unum`:  one to hew to the `/lib/math` and `/lib/saloon` interface, and another to alias in the 2022 Posit Standard requirements.  The names are slightly modified from the standard to adhere to Hoon name requirements.

### Basic functions of one posit value argument

- `++negate` ← `(sub 0 val)`
- `++abs` ← no change
- `++sign` ← `++sgn`
- `++nearest-int`
- `++ceil`
- `++floor`
- `++next`
- `++prior`

### Comparison functions of two posit value arguments

- `++compare-equal` ← `++equ`
- `++compare-not-equal` ← `!(equ val)`
- `++compare-greater` ← `++gth`
- `++compare-greater-equal` ← `++gte`
- `++compare-less` ← `++lth`
- `++compare-less-equal` ← `++lte`

### Arithmetic functions of two posit value arguments

- `++addition` ← `++add`
- `++subtraction` ← `++sub`
- `++multiplication` ← `++mul`
- `++division` ← `++div`

### Elementary functions of one posit value argument

- `++sqrt` ← no change
- `++rsqrt` ← TODO
- `++exp` ← no change
- `++exp-minus-1` ← `++expm1`
- `++exp2` ← no change
- `++exp2-minus-1` ← `++exp2m1`
- `++exp10` ← no change
- `++exp10-minus-1` ← `++exp10m1`
- `++log` ← no change
- `++log-plus-1` ← `++logp1`
- `++log2` ← no change
- `++log2-plus-1` ← `++log2p1`
- `++log10` ← no change
- `++log10-plus-1` ← `++log10p1`
- `++sin` ← no change
- `++sin-pi` ← `(sin (mul pi val))`
- `++cos` ← no change
- `++cos-pi` ← `(cos (mul pi val))`
- `++tan` ← no change
- `++tan-pi` ← `(tan (mul pi val))`
- `++arcsin` ← no change
- `++arcsin-pi` ← `(arcsin (div val pi))`
- `++arccos` ← no change
- `++arccos-pi` ← `(arccos (div val pi))`
- `++arctan` ← no change
- `++arctan-pi` ← `(arctan (div val pi))`
- `++sinh` ← no change
- `++cosh` ← no change
- `++tanh` ← no change
- `++arcsinh` ← no change
- `++arccosh` ← no change
- `++arctanh` ← no change

### Elementary functions of two posit value arguments

- `++hypot` ← `(sqt (mul val-1 val-1) (mul val-2 val-2))`
- `++pow` ← no change
- `++arctan2` ← ??? complex
- `++arctan2-pi` ← `(div (arctan2 val) pi)`

### Functions of three posit value arguments

- `++fmm` ← `:(mul val-1 val-2 val-3)`

### Functions of one posit value argument and one integer argument

- `++compound` ← `(pow (add val 1) int)`
- `++root-n` ← `(pow val (div 1 int))`

### Functions involving quire value arguments

- `++p-to-q`
- `++q-negate`
- `++q-abs`
- `++q-add-p`
- `++q-sub-p`
- `++q-add-q`
- `++q-sub-q`
- `++q-mul-add`
- `++q-mul-sub`
- `++q-to-p`

### Conversions between different precisions

- `++p8-to-16`
- `++p8-to-32`
- `++p16-to-8`
- `++p16-to-32`
- `++p32-to-8`
- `++p32-to-16`

### Conversions between posit format and decimal format

### Conversions between posit format and IEEE 754 format
