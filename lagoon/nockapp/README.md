# lagoon/nockapp

`/lib/lagoon` on NockVM: the Hoon, compiled with `hoonc`, plus Rust jets for
the `%i754` array arms built on [`sdfloat`/`sdblas`](https://github.com/sigilante/RustFloat),
so a NockApp computes exactly the bits the Hoon (and the vere C jets) produce.

```
crates/lagoon-jets/     the lagoon jets (`LAGOON_HOT`); `hot_state()` = hoon-138 built-ins
                        + the float-door jets + these
  src/bin/lagoon-kick   test driver: fires a hoonc --arbitrary trap under the jets
crates/hoon-float-jets/ jets for hoon-138's own rh/rs/rd/rq doors (`HOON_FLOAT_HOT`)
hoon/                   generated: lib/ sur/ copies of the desk with `..part` -> `..ut`,
                        the test files, run-tests.hoon (the lagoon suite as lazy traps),
                        float-tests.hoon and float-diff.hoon (the float doors)
scripts/sync-hoon.sh    regenerates hoon/lib, hoon/sur, hoon/run-tests.hoon
scripts/gen-float-*.py  regenerate the float-door tests; float-diff-report.py classifies
PORTING.md              rules for porting an arm; read before editing jets/
```

## Why a copy of the Hoon

Vere libraries register jets under `..part` (Arvo). hoonc compiles against
bare `hoon-138`, where the parent is `..ut`; the jet path becomes
`[k.138 one two tri qua pen non lagoon <arm>]`. That one token is the only
difference from `../desk`, so the copy is generated, never edited.

## Build and test

Toolchain: the pinned nightly in `rust-toolchain.toml` (nockvm needs it);
`hoonc` 0.2.0 on `PATH` for the Hoon.

```sh
scripts/sync-hoon.sh                                   # after any desk change
(cd hoon && hoonc --arbitrary --output run-tests.jam run-tests.hoon .)   # ~1-7 min
cargo build --release
./target/release/lagoon-kick hoon/run-tests.jam        # 153 tests, ~1.5 s
```

Differential check, the way the jets are verified: NockVM runs both the jet
and the raw Nock for every jet named in `NOCK_TEST_JETS` and bails with
`%jest` on any mismatch.

```sh
NOCK_TEST_JETS=k.138/one/two/tri/qua/pen/non/lagoon/add-rays,k.138/one/two/tri/qua/pen/non/lagoon/mmul \
  ./target/release/lagoon-kick hoon/run-tests.jam
LAGOON_TESTS=arithmetic/,linalg/ ...                   # filter rows by substring
LAGOON_JET_TRACE=1 ...                                 # one stderr line per jet call
LAGOON_JET_SABOTAGE=1 ...                              # harness self-check: add-rays is made wrong
```

## The float doors

`crates/hoon-float-jets` jets `add sub mul div sqt lth lte equ gte gth` of
hoon-138's `++rh`/`++rs`/`++rd`/`++rq` with sdfloat. Those doors already
carry `~%`/`~/` hints, so the jets register at `k.138/one/two/tri/<door>/<arm>`
with no Hoon change; `fma` stays Nock (sdfloat has no fused multiply-add).
Every float operation on NockVM, not just lagoon's, becomes bit-exact and fast.

```sh
scripts/gen-float-tests.py && (cd hoon && hoonc --arbitrary --output float-tests.jam float-tests.hoon .)
NOCK_TEST_JETS=k.138/one/two/tri/rs/add,...  ./target/release/lagoon-kick hoon/float-tests.jam
```

`float-tests.hoon` is 320 rows, one per (door, arm, value set, mode). Under
test mode the only rows that mismatch are add/sub/mul/div on the
overflow-capable value set in `%u`/`%d`/`%z`: the Hoon `++fl` overflows to
infinity in every mode where IEEE 754 (and sdfloat) saturates to the largest
finite value (urbit/urbit#7426). `scripts/float-diff-report.py` proves that
is the whole difference: over 6400 input pairs the jets and the pure Hoon
disagree 592 times, each time as ±MAX against ±inf, never otherwise. Once the
Hoon is fixed those rows pass unchanged. `HOON_FLOAT_JET_DISABLE=1` leaves
the doors to the Nock; `HOON_FLOAT_JET_TRACE=1` prints each call.

## Using the jets in a NockApp

```rust
let hot = lagoon_jets::hot_state();          // URBIT_HOT_STATE ++ HOON_FLOAT_HOT ++ LAGOON_HOT
boot::setup(&kernel_jam, cli, &hot, "my-app", None).await?;
```

`Hot::init` registers only what it is given: passing `LAGOON_HOT` alone
leaves `add`, `dec` and the bit operations running as raw Nock.

## Jets vs Hoon under NockVM

Any input a jet does not handle punts to the Nock, so jets are never a
correctness risk, only a speed one. Where the vere C jets and the Hoon
disagree, the Hoon wins here (test mode enforces it); see `PORTING.md` for
the cases the port surfaced.
