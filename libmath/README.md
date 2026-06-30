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

A library for universal numbers (posits, quires, and valids) per the [2022 Posit Standard](https://posithub.org/docs/posit_standard-2.pdf) (Posit Working Group, John Gustafson chair).

We implement the three standard precisions ($\text{posit}\langle 8,2\rangle$, $\text{posit}\langle 16,2\rangle$, $\text{posit}\langle 32,2\rangle$) plus two extensions beyond the 2022 standard ($\text{posit}\langle 64,2\rangle$ and $\text{posit}\langle 128,2\rangle$), which are best-effort.

> **Standard, not legacy.**  The 2022 Posit Standard fixes the exponent size at **`es = 2` for every width** (§2 defines the exponent field as a 2-bit unsigned integer).  This differs from the original 2017 draft (and from SoftPosit's *fast* `p8`/`p16` types), which scaled `es` with width as `posit<8,0>`, `posit<16,1>`, `posit<32,2>`.  Only `posit32` coincides between the two conventions; `posit8` and `posit16` have different bit layouts.  We target the standard.  The reference oracle for 8- and 16-bit posits is therefore SoftPosit's `pX2` (es=2 at any width) path, *not* the fast `p8`/`p16` types; `p32` is used directly.

##  Format and encoding

A posit of precision $n$ has four bit-fields, most-significant first (§3.3):

| field | width | meaning |
|---|---|---|
| sign $S$ | 1 bit | integer $s \in \{0,1\}$; implicit value $1-3s$ (so $+1$ or $-2$) |
| regime $R$ | $k$ bits | run of bits equal to $R_0$, terminated by $\overline{R_0}$; signed integer $r = -k$ if $R_0=0$, else $r = k-1$ |
| exponent $E$ | 2 bits | unsigned $e \in \{0,1,2,3\}$ (`es = 2`); may be truncated past the LSB, then treated as 0 |
| fraction $F$ | up to $\max(0,n-5)$ bits | unsigned $f = F / 2^m$, $0 \le f < 1$; trailing truncated bits are 0 |

The represented value (with `useed` $= 2^{2^{es}} = 16$, so the regime contributes $4r$ to the power of two) is:

$$p = \bigl((1-3s) + f\bigr) \times 2^{(1-2s)\,(4r + e + s)}$$

**Exceptions** (all bits except $S$ are zero): $S=0 \Rightarrow$ posit $0$; $S=1 \Rightarrow$ **NaR** (Not a Real — the single non-real value, with bit pattern equal to the most-negative two's-complement integer, $1000\ldots0$).

Key constants per width (§3, Table 1; $n$ = bit width):

| property | posit8 | posit16 | posit32 | $\text{posit}n$ |
|---|---|---|---|---|
| fraction length | 0–3 | 0–11 | 0–27 | $0$–$\max(0,n-5)$ |
| $\textit{minPos}$ | $2^{-24}$ | $2^{-56}$ | $2^{-120}$ | $2^{-4n+8}$ |
| $\textit{maxPos}$ | $2^{24}$ | $2^{56}$ | $2^{120}$ | $2^{4n-8}$ |
| $\textit{pIntMax}$ | $16$ | $1024$ | $8388608$ | $2^{\lceil 4(n-3)/5\rceil}$ |
| quire bits | 128 | 256 | 512 | $16n$ |
| decimals to round-trip | 2 | 5 | 10 | — |

$\textit{minPos} = 1/\textit{maxPos}$, and every posit value is an integer multiple of $\textit{minPos}$.  $\textit{pIntMax}$ is the largest integer all of whose predecessors are exactly representable (so posit8 represents every integer only up to 16).

### Rounding

Posits have exactly **one** rounding mode (§4.1): round-to-nearest, **ties-to-even**, and **saturating** — a magnitude above $\textit{maxPos}$ rounds to $\pm\textit{maxPos}$ and a nonzero magnitude below $\textit{minPos}$ rounds to $\pm\textit{minPos}$.  Posits never round up to NaR nor down to zero.  Fused operations (`++fmm`, `++rsqrt`, `++exp-minus-1`, `++hypot`, …) are evaluated exactly and rounded once.

### Comparisons

Posit ordering is **identical to two's-complement signed-integer ordering of the raw bits** (§5.3), so comparisons need no decoding — a single signed compare suffices.  NaR shares the bit pattern of the most-negative integer, so it compares less than every real posit; `NaR == NaR` is true and `++sign` of NaR is NaR.

### Quire

The quire is a fixed-point exact accumulator of $16n$ bits (§3.4): sign (1) · carry guard (31) · integer part ($8n-16$) · fraction part ($8n-16$).  Its value is $2^{16-8n}$ times the two's-complement integer formed by all bits; $S=1$ with all other bits zero is NaR.  The quire makes dot products and sums of up to $\sim 2^{23+4n}$ terms exact before a single final rounding — the core reason posits beat floats for linear algebra.

### Auras

| aura | meaning | width |
|---|---|---|
| `@rpb` `@rph` `@rps` `@rpd` | posit | byte (8), half (16), single (32), double (64) |
| `@rqb` `@rqh` `@rqs` | quire | for posit8 / posit16 / posit32 |
| `@rvb` `@rvh` `@rvs` | valid (interval) | byte / half / single |

The `@rq*` quire auras nest under `@rq` (quad-precision float) by Hoon's prefix-nesting rule; this overlap is accepted intentionally.  There is currently **no literal syntax or pretty-printer** for posit auras — values are written as bit-casts, e.g. `` `@rpb`0b100.0000 `` for $1.0$.  A decimal printer (emitting the §6.3 minimum significant digits: 2/5/10/21 for posit8/16/32/64) would be an *optional local runtime patch*, decoupled from this library, with no upstream-timing commitment.

### Implemented

`lib/unum.hoon` implements the `es = 2` layout above as a single generic core `++pp` (keyed on `bloq`), specialized into `++rpb` / `++rph` / `++rps` / `++rpd` / `++rpq` (posit8/16/32/64/128).  The g-layer record `+$up` mirrors the standard library float `+$fn`.  Each width exposes:

- **decode / encode** — `++sea` (`@` → `+$up`), `++bit` (`+$up` → `@`, round-to-nearest-even, saturating)
- **constants** — `++zero` `++nar` `++one` `++maxpos` `++minpos` `++huge` `++tiny`, and `++pi` `++tau` `++e` `++phi` `++sqt2` `++invsqt2` `++log2` `++invlog2` `++log10` (correctly rounded at every width)
- **comparison / sign** — `++gth` `++lth` `++gte` `++lte` `++equ` `++neq`, `++neg` `++abs` `++sgn`
- **arithmetic** — `++add` `++sub` `++mul` `++div` (correctly rounded), `++fma` (fused multiply-add), `++sqt` (square root)
- **rounding** — `++round` with `++rnd` (nearest-even) / `++flr` (floor) / `++cel` (ceil)
- **integer conversion** — `++sun` (`@u`→), `++san` (`@s`→), `++toi` (→`(unit @s)`)
- **IEEE-754 conversion** — `++to-rh` `++to-rs` `++to-rd` `++to-rq` and `++from-rh` `++from-rs` `++from-rd` `++from-rq`, value-based across *any* posit/float width pair
- **quire** (the `16n`-bit exact accumulator) — `++p-to-q` `++q-to-p` `++q-mul-add` `++q-mul-sub` `++q-add-p` `++q-sub-p` `++q-add-q` `++q-sub-q` `++q-negate`, and `++fdp` (fused dot product, single rounding)
- **elementary** — `++exp` `++sin` `++cos` `++tan` `++pow-n` `++log` `++log-2` `++log-10` `++pow` `++is-close` (naive Taylor series in the `/lib/math` style: reproducible, accurate near 0, not range-reduced)

Arithmetic is verified against SoftPosit (the reference C implementation, via its Python wrapper): exhaustively over all posit8 pairs and sampled at posit16/32, with the offline harness in `tools/posit_check.py`; the on-ship suite is `tests/lib/unum-core` and `tests/lib/unum-fns`.

### Not yet implemented

- the standard-name alias interface (below): `++negate` `++addition` `++compare-less` `++sin-pi` `++compound` `++root-n` `++fmm` `++hypot` `++arctan2` etc., and the inverse / hyperbolic / `*-plus-1` / `*-minus-1` elementary functions
- `++next` / `++prior` / `++nearest-int` / `++ceil` / `++floor` under their standard names (the behaviors exist as `++rnd`/`++flr`/`++cel`)
- valids (the interval unum class, `@rv*`)
- jets (the library is pure Hoon; SoftPosit's C is the spec for a future jet)

##  Posit Standard Compliance (planned alias layer)

The following is the *planned* second interface for `/lib/unum`: standard-named aliases over the implemented operations above, plus the not-yet-written functions.  One interface hews to the `/lib/math` and `/lib/saloon` conventions; this one aliases in the 2022 Posit Standard names, slightly modified to adhere to Hoon name requirements.

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
