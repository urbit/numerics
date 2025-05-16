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
- `++sin`, $\sin$
- `++cos`, $\cos$
- `++tan`, $\tan$
- `++pow-n`, $\text{pow}$ to integer power
- `++log`, $\log$ (natural logarithm)
- `++log-10`, $\log_{10}$ (log base-10)
- `++log-2`, $\log_{2}$ (log base-2)
- `++pow`, $\text{pow}$
- `++sqrt`, $\sqrt{}$ (also `++sqt`)
- `++cbrt`, $\sqrt[\uproot{3}]{}$ (also `++cbt`)
- `++arg` (alias for `++abs` in real values)

Logical functions:

- `++lth`, $<$
- `++lte` $\leq$ (also `++leq`)
- `++gth`, $>$
- `++gte`, $\geq$ (also `++geq`)
- `++equ`, $=$
- `++neq`, $\neq$
- `is-close`
- `all-close`
- `is-int`

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
