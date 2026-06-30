# `/lib/unum` jets — vere reference

Reference (master) copies of the hand-maintained C jet sources for the posit
(`/lib/unum`) jets, mirrored by hand into the vere runtime (`urbit/vere`,
branch `sigilante/unum-jets`).  SoftUnum itself is *vendored* into vere
(`ext/softunum`, pinned tarball of `sigilante/SoftUnum`), so the only
hand-synced jet file is `unum.c`.

## Files here

- `noun/jets/i/unum.c` — the jet wrappers.  One `%unum` core, each op
  dispatched on `bloq` read from the `pp` door sample (gate axis 30), calling
  SoftUnum `p8_*`/`p16_*`/`p32_*`.  posit64/128 (bloq 6/7) return `u3_none`
  (fall back to the pure-Hoon arm) until SoftUnum covers them.

## Deltas applied in vere (not full copies — see the vere branch/PR)

- `ext/softunum/{build.zig,build.zig.zon}` — vendor SoftUnum (mirror
  `ext/softblas`); wired into `pkg/noun/build.zig{,.zon}`.
- `pkg/noun/build.zig` — add `jets/i/unum.c` to the noun sources.
- `pkg/noun/jets/w.h`, `q.h` — declare `u3wi_unum_*` / `u3qi_unum_*`.
- `pkg/noun/jets/135/tree.c` — register `non/unum/<arm>` in the **hoon-135**
  dashboard.  We ship to the lowest (current) kelvin only; the 408k pill boots
  hoon-135, so the 135 dashboard is the one consulted.  (Local testing against
  newer kelvins may also touch `136/137/tree.c`, but those are NOT part of the
  shipped change — only 135.)

## Hoon side (`libmath/desk/lib/unum.hoon`)

The library is jet-hinted to match: `~% %non ..part ~` on the file core,
`~/ %unum` on `++pp`, and `~/ %<arm>` on each jetted arm (e.g. `~/ %add`).

## Status

`++add:pp` **fires bit-exact and verified on a live hoon-135 fakezod**:
`(add:rpb:unum 0x40 0x40)` → `0x48`, `add:rph` → `0x4800`, `add:rps` →
`0x4800.0000`, each dispatching the right width (bloq 3/4/5 read from the door
sample at gate axis 30); posit64 (bloq 6) declines to the pure-Hoon arm and
still returns the correct value.  The remaining scalar/quire/conversion arms
follow the same one-jet-per-op, bloq-dispatch pattern.

# libmath C jets — 32-bit vere mirror

Vendored from the 32-bit `urbit/vere` runtime (`pkg/noun/jets`). This is the
32-bit twin of `libmath/vere64/`; the two `math.c` copies differ by exactly two
lines (the `u3r_mean` two-arg form). See **`../vere64/README.md`** for the full
documentation, apply steps, the Hoon jet-structure gotcha, and verification.
