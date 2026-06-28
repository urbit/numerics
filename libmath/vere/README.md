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
- `pkg/noun/jets/{135,136,137}/tree.c` — register `non/unum/<arm>` in **all
  three** dashboards (a jet registered only in 135 will not fire on a
  hoon-136/137 ship).

## Hoon side (`libmath/desk/lib/unum.hoon`)

The library is jet-hinted to match: `~% %non ..part ~` on the file core,
`~/ %unum` on `++pp`, and `~/ %<arm>` on each jetted arm (e.g. `~/ %add`).

## Status

`++add:pp` (posit8/16/32) implemented and building; on-ship firing
verification pending.  The remaining scalar/quire/conversion arms follow the
same one-jet-per-op, bloq-dispatch pattern.
