# `/lib/twoc` jets — vere reference

Reference (master) copy of the hand-maintained C jet source for the
two's-complement integer (`/lib/twoc`) jets, mirrored by hand into the vere
runtime.  Unlike `/lib/unum` (SoftUnum) or Lagoon (SoftBLAS), `/lib/twoc` needs
**no vendored library** — the arithmetic is plain integer C plus GNU MP, and
GMP is already vendored (`ext/gmp`) and linked into `pkg/noun`.

## What is jetted

`/lib/twoc` keeps all its logic in ONE width-keyed door `++twid` (`|_ wid=@`);
the bloq-keyed `++twoc` facade just re-exports `twid` gates at width
`(bex bloq)`.  So we register a SINGLE `%twid` core and **every** caller fires
through it:

- `/lib/fixed`            — `twid:twoc` directly, at arbitrary widths `a+b+1`.
- `/lib/unum`, Lagoon `%int2` — via the `twoc` facade (bloq → `bex bloq`).

First pass jets the 12 hot-path arms Lagoon `%int2` drives element-wise:
`add sub neg mul div rem pow abs gth lth lte gte`.  (`msb`, `s-to-twoc`,
`twoc-to-s` are a planned second pass.)

`wid` is read from the door sample at gate axis 30 (gate axis 7 = door, door
axis 6 = `wid`), the same shape `unum.c` uses for `bloq`.

## Dispatch: native c3 types vs GNU MP

| `wid`        | path                              |
|--------------|-----------------------------------|
| `wid <= 64`  | native `c3_d` (uint64), masked    |
| `wid <= 128` | native `unsigned __int128`, masked|
| `wid  > 128` | GNU MP (`mpz`), masked to `wid`   |

The power-of-two bloq widths (8/16/32/64/128 for bloq 3..7) all land in the
native paths; GMP is the exact backstop for genuinely arbitrary precision. All
native arithmetic is unsigned with explicit masking, so there is no
signed-overflow UB (`div(min,-1)` wraps to `min`, matching the Hoon).

## Declining to Hoon (bit-exactness for free)

The jet returns `u3_none` (pure-Hoon arm runs) when:

- an operand does not fit in `wid` bits (`u3r_met(0,x) > wid`) — the Hoon
  derives sign from raw bit-width (`xeb`), which diverges from masking only
  off-contract; declining keeps us bit-exact without replicating `xeb` edges.
- a divisor is zero in `div`/`rem` — the Hoon `div`/`rem` crash; we let them
  (and avoid a native `SIGFPE` in the serf).

## Validation

Standalone, before any ship build (`libmath/tools/twoc/`): the jet's arithmetic
is copied verbatim and checked against an independent two's-complement reference
on hardware `__int128`/`mpz` — `twoc_native_test.c` (native paths, 8280 checks,
widths 8/16/17/32/64/100, all sign edges) and `twoc_gmp_test.c` (GMP
cross-checked vs the hardware reference at widths 8..100 plus hand-verified
200-bit cases, 5692 checks).  Both pass.

On a live hoon-135 fakezod the `twid` jet FIRES bit-exactly through both call
paths, proven two ways: (1) a sentinel `+1` planted in the native `add` showed
`1` for jetted `add(0,0)` via BOTH `twid:twoc` (direct) and `twoc:twoc` (facade,
bloq 3 & 6), `0` interpreted — so a single `%twid` registration fires for both
callers; (2) timing the identical loop with vs without the hint: `twid`-direct
1M add-64 = 2.85 s interpreted → 0.87 s jetted (~3.3x); `twoc`-facade is ~1.75x
(the facade re-instantiates the door per call).  The modest margin is expected:
the interpreted arm is `(mod (^add a b) (bex wid))` and `^add`/`bex`/`mod` are
already jetted, so the jet collapses 3 dispatches + glue into 1 + native rather
than replacing interpreted bignum math.  Bigger wins need wide widths or
array-level jetting (one C call per `%int2` array op).

## Deltas to apply in vere (not full copies — see the vere branch/PR)

- `pkg/noun/jets/i/twoc.c` — this file (master copy lives in numerics).
- `pkg/noun/build.zig` — add `jets/i/twoc.c` to the noun sources (next to
  `jets/i/unum.c`).  No new `linkLibrary` — GMP is already linked.
- `pkg/noun/jets/w.h`, `q.h` — declare `u3wi_twid_*` / `u3qi_twid_*` (binops &
  comparisons `(c3_d,u3_atom,u3_atom)`, unary `(c3_d,u3_atom)`).
- `pkg/noun/jets/135/tree.c` — register a `twid` core under `non` (sibling of
  `unum`/`math`/`lagoon`): `{ "twid", 7, 0, _135_non__twid_d, no_hashes }`,
  with `_135_non__twid_<arm>_a[] = {{".2", u3wi_twid_<arm>}, {}}` and
  `_135_non__twid_d[]` listing the 12 arms at axis 7.  hoon-135 only (lowest
  shipped kelvin), as with unum.

## Hoon side (`libmath/desk/lib/twoc.hoon`)

Jet-hinted to match: `~% %non ..part ~` on the file core, `~/ %twid` on
`++twid`, and `~/ %<arm>` on each of the 12 jetted arms.  The `++twoc` facade is
left unhinted — it re-exports `twid` gates, which already carry the
registration, so a single `%twid` jet fires for both call paths.  (Verify by
FIRING, not by jet count.)
