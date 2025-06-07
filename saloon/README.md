#   Scientific ALgorithms in hOON

Transcendental and algebraic functions for use with Lagoon `+$ray`s.

We support the following functions and special functions:

- `++add`, $+$ addition (pass-through from Lagoon)
- `++sub`, $-$ subtraction (pass-through from Lagoon)
- `++mul`, $\times$ multiplication (pass-through from Lagoon)
- `++div`, $/$ division (pass-through from Lagoon)
- `++fma`, $\text{fma}$ fused multiply-add
- `++neg`, $-$ unary negation
- `++factorial`, $!$ factorial
- `++abs`, $\text{abs}$ (pass-through from Lagoon)
- `++exp`, $\exp$
- `++sin`, $\sin$
- `++cos`, $\cos$
- `++tan`, $\tan$
- `++pow-n`, $\text{pow}$ to integer power
- `++log`, $\log$ (natural logarithm)
- `++log-10`, $\log_{10}$ (log base-10)
- `++log-2`, $\log_{2}$ (log base-2)
- `++pow`, $\text{pow}$
- `++sqrt`, $\sqrt$ (also `++sqt`)
- `++cbrt`, $\cbrt$ (also `++cbt`)

Logical functions:

- `++lth`, $<$ (pass-through from Lagoon)
- `++lte` $\leq$ (also `++leq`) (pass-through from Lagoon)
- `++gth`, $>$ (pass-through from Lagoon)
- `++gte`, $\geq$ (also `++geq`) (pass-through from Lagoon)
- `++equ`, $=$ (pass-through from Lagoon)
- `++neq`, $\neq$ (pass-through from Lagoon)
- `is-close` (pass-through from Lagoon)
- `all-close` (pass-through from Lagoon `++all`)
- `any-close` (pass-through from Lagoon `++any`)

##  References

- Milton Abramowitz & Irene Stegun, _Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables_.  1964–2010.
- Forman Acton, _Numerical Methods that (Usually) Work_, 1ed.  1997.
- [Bartosz Ciechanowski, “Float Exposed” (webapp)](https://float.exposed/0x00000001)
- [David Goldberg, “What Every Computer Scientist Should Know About Floating-Point Arithmetic”](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
- Parviz Moin, _Fundamentals of Engineering Numerical Analysis_. 2ed.  2001.
