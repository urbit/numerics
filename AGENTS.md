# Agent merge guide: numerics → urbit/urbit and urbit/vere

How to propagate a change across the three repositories that hold the Urbit
numerical stack. Read this before syncing `/lib/math`, `/lib/unum`, `/lib/twoc`,
`/lib/lagoon`, `/lib/complex`, `/lib/fixed`, saloon, or maroon.

## The three repos and their roles

| repo | holds | path |
| --- | --- | --- |
| **urbit/numerics** | **canonical master** — the Hoon libs AND a C jet reference mirror | `<fam>/desk/…`, `<fam>/vere{,64}/…` |
| **urbit/urbit** | the Arvo desk that ships the Hoon | `pkg/arvo/lib/*.hoon`, `pkg/arvo/sur/*.hoon` |
| **urbit/vere** | the runtime that compiles and registers the jets | `pkg/noun/jets/…`, `pkg/noun/build.zig` |

**numerics is the source of truth.** The Hoon in `<fam>/desk/lib/*.hoon` is the
normative spec; a jet must match its Hoon arm bit-for-bit. Changes flow
numerics → urbit (Hoon) and numerics → vere (jets, applied by hand).

Exception in practice: a fix sometimes lands in a **vere PR first** (typically
refcount-linter or build fixes). When it does, **mirror it back to numerics** so
there is no delta. The invariant is that all three stay in sync, not that edits
only ever start in numerics.

## numerics layout (per family, per loom)

Each family directory with jets (`lagoon/`, `libmath/`, `maroon/`, `saloon/`)
contains:

- `<fam>/desk/` — the Hoon: `lib/*.hoon`, `sur/*.hoon`, `tests/`. This ships to
  urbit/urbit.
- `<fam>/vere/noun/jets/` — the C jet **reference mirror**: `i/*.c`, `i/*.h`,
  `q.h`, `w.h`, `<kelvin>/tree.c`. This is applied by hand into urbit/vere.
- `<fam>/vere64/…` — the post-vere64 twin, where the split still exists.

These jet mirrors are **per-family subsets**, not full vere trees. `libmath/`
carries only the math/unum/twoc jets; `lagoon/` only lagoon; and so on.

## The vere64 split (read this before editing any jet)

vere `develop` merged the 64-bit loom (`vere64`, urbit/vere#970): **one** source
now builds both looms via `-Dvere64`, with uniform names (`c3_h`, `u3i_half`)
instead of `#ifdef`. Direct atoms are 31-bit on the 32-bit loom, 63-bit on the
64-bit loom.

numerics still keeps two mirrors per family:

| | `<fam>/vere/` (pre-vere64, 32-bit twin) | `<fam>/vere64/` (post-vere64) |
| --- | --- | --- |
| half-word type | `c3_w` | `c3_h` |
| `u3r_mean` form | variadic: `u3r_mean(cor, ax, &a, …, 0)` | braced: `u3r_mean(cor, {ax, &a}, …)` |

- A change must go to **both** `vere/` and `vere64/` where the split exists.
- **Loom-agnostic** changes (a `u3_weak` retype, a comment/`@Refcount`
  annotation, a new arm's logic) are **byte-identical** in both mirrors.
- **Loom-specific** lines differ only in the two rows above.
- `lagoon/` is being collapsed into a single both-loom `lagoon/vere/` (mirrors
  unified vere develop); once merged, treat lagoon as one mirror. Check whether
  `<fam>/vere64/` still exists before assuming a split.

## Cross-repo diffs are noisy — trust the `.c`, not the headers

numerics carries **per-family** jet mirrors; the vere branches are **full
trees**. So a whole-file diff of the shared headers — `q.h`, `w.h`,
`<kelvin>/tree.c` — is dominated by **other jet families** and is not a real
delta. Do not read those diffs wholesale.

- **Use the family's `.c` file as the unit of truth.** It is self-contained;
  verify a sync by **byte-identity** (`diff -q`), not by staring at a header diff.
- To isolate a real (non-typing) delta, normalize the typedef on both sides
  first: `diff <(sed 's/u3_weak/u3_noun/g' A) <(sed 's/u3_weak/u3_noun/g' B)`.
- For `q.h`/`w.h`/`tree.c`, compare **only the family's own lines**
  (`grep 'u3qi_<fam>_'`, `grep '<fam>'`) or diff against the **same family's**
  mirror — never against the full vere tree.

## The refcount linter (dozreg, urbit/vere#1059)

Jets must pass dozreg's checker with **0 findings, including `--strict-weak`**.

- Any wrapper/helper that can `return u3_none` (a jet that punts to Hoon) must be
  typed **`u3_weak`** — in the `.c` definition **and** the `q.h`/`w.h`
  declaration. `u3_weak` is a `typedef` of `u3_noun` (no runtime effect); the
  type is what the linter reads.
- A local holding a possibly-none product (`u3r_at`, another `u3_weak` call) must
  be `u3_weak`.
- Pointer out-params carry an `@Refcount:` annotation, e.g.
  `//  @Refcount: fills transferred `val`` for a callee that always writes a
  transferred noun to `*val`; add `` on `c3y` `` for a conditional (u3r_cell-style)
  fill.
- Reads that assume a direct atom go through `u3r_cat`, not `u3x_atom`-as-integer.
- Run the linter locally before every jet PR; see the vere-side run notes.

## Standard sync workflow (change originates in numerics)

1. Edit the Hoon in `<fam>/desk/lib/*.hoon` (+ `sur/`, + `tests/`). This is the
   spec.
2. Edit the jet reference in `<fam>/vere/` **and** `<fam>/vere64/` (both looms).
   Keep the `.c` bit-exact to the Hoon arm.
3. Mirror the Hoon to **urbit/urbit** `pkg/arvo/lib/` (+ `sur/`, + tests).
4. Mirror the jets to **urbit/vere** `pkg/noun/jets/{i/*.c, q.h, w.h,
   <kelvin>/tree.c}` and add the `.c` to `pkg/noun/build.zig` if new.
5. Verify (all three): build both looms `-Werror`; refcount linter 0 under
   `--strict-weak`; run the family's desk test suite **jetted** on a fakezod
   (tmux driver), matching unjetted Hoon.

## Back-port workflow (fix landed in a vere PR first)

1. Pick the matching numerics mirror by loom form (braced `u3r_mean` → `vere64/`;
   variadic → `vere/`).
2. The family's `.c` can usually be **copied wholesale** from the vere branch
   (numerics = "vere minus the other families" for one family's `.c`). For
   `q.h`/`w.h`/`<kelvin>/tree.c`, merge in **only that family's lines** — never
   clobber the other families the numerics header already carries.
3. Verify the family's `.c` is byte-identical to the vere source; PR to numerics.

## Practical mechanics

- Work in **worktrees**, one per target branch; never disturb scratch branches.
- **Worktree builds** fail on `zig build` because `build.zig` reads
  `.git/logs/HEAD` and a worktree's `.git` is a file. Temporarily
  `mv .git .git-file-bak && ln -s <repo>/.git/worktrees/<name> .git`, build,
  then restore the file before committing.
- If a repo is on its default branch, **branch first**. One focused PR per
  concern (don't bundle unrelated families or a build fix with a jet fix).
- Registration ships to the **current kelvin** only (`<kelvin>/tree.c`, 135 for
  the 408k pill); local testing may touch 136/137 but those are not shipped.

## Provenance rule

Before trusting a "the spec says …" comment, or any C behavior that diverges from
the Hoon, **trace it in git**. numerics Hoon is ground truth. (Cautionary case: a
truncating C-`fmod` `%mod` was once fabricated against the Hoon, which always used
`+toi`; it was reverted across desk Hoon, both jets, and tests.)

## Pointers

- Per-family details: `<fam>/README.md` (e.g. `libmath/vere/README.md`,
  `libmath/vere64/README.md`, `lagoon/README.md`).
- Review/plan state: `WORKPLAN.md`; reference PDFs and audits in `doc/`.
