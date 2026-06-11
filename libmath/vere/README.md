# libmath C jets (vere mirror)

C jet implementations for the `math.hoon` transcendental library, **vendored
from the vere runtime** (`pkg/noun`). The jets do not live in this repo's build
— this directory is a maintained mirror so the C and the Hoon evolve together.
There is no automatic cross-repo sync; apply changes to vere by hand.

## What's here
- `jets/i/math.c` — the jet bodies. Each runs the IDENTICAL algorithm as its
  Hoon arm in Berkeley SoftFloat, so jet output is **bit-exact** to the pure-Hoon
  reference (not merely faithful). Marshalling is chub-based (word-size-agnostic;
  the divergence that broke the @rq sub jet). Currently: `@rd exp`.
- `jets/tree-registration.snippet.c` — the `tree.c` registration block to add to
  each per-kelvin jet tree, plus the q.h/w.h/build.zig one-liners.

## How to apply to vere (`pkg/noun`)
1. Copy `jets/i/math.c` → `pkg/noun/jets/i/math.c`.
2. Add the registration block from `jets/tree-registration.snippet.c` to
   **each** of `pkg/noun/jets/13{5,6,7}/tree.c` (the active tree is picked by the
   ship's `hoon-version`; numerics ships are 136, stock test ships 135). `math`
   registers under the `non` chapter as a sibling of `lagoon`.
3. Declare in `pkg/noun/jets/q.h` (`u3qi_rd_exp`) and `w.h` (`u3wi_rd_exp`).
4. Add `"jets/i/math.c"` to the jet source list in `pkg/noun/build.zig`.
5. Rebuild; on a 64-bit build use `-Dvere64=true`.

## CRITICAL: the Hoon jet structure (the gotcha that cost a day)
`math.hoon` must mirror `/lib/lagoon`'s jet structure exactly:

```
=<  math                       :: export the engine so `/+ math` + `rd:math` work
~%  %non  ..part  ~            :: file |% anchored to the shared %non chapter
|%
++  math
  ^|
  =+  [rnd=*?(%n %u %d %z)]     :: <-- REQUIRED: the engine must be a DOOR
  ~/  %math                    ::     (sample-bearing), not a bare |%
  |%
  ++  rd
    ~/  %rd
    ^|
    |_  [r=... rtol=...]
    ++  exp  ~/  %exp  |=(x=@rd ...)
```

A **sample-less `~/`-hinted core anchors its parent to the dashboard ROOT**, so
`%math` (and everything under it) misses the hot dashboard — the core SPOTs in
`_cj_mile` but `_cj_hot_mean` returns `jax=0` and no C jet ever attaches. The
`^| =+ [rnd]` sample makes `~/` resolve the parent to `%non`, exactly as
lagoon's `++ la ^| =+ [rnd] ~/ %lagoon |%` does. Verified via a `_cj_hot_mean`
trace: `math par=0 MISS` (bare |%) -> `math par=non FOUND` (door).

The chapter chain is `pen -> hex -> non -> math -> <width door> -> <fn>`, e.g.
`non/math/rd/exp`.

## Verifying a jet fires
`~>(%bout (exp:rd:math .~1))` should report a few µs (jetted) vs hundreds of µs
(interpreted), and `` `@ux`(exp:rd:math .~1) `` must equal `0x4005.bf0a.8b14.576a`
(bit-exact). Use `tmux` for the ship so the dojo stays interactive; jet
registration only warms on first computation of each core (not at boot).
