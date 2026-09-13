# lagoon/nockapp

`/lib/lagoon` on NockVM: the Hoon, compiled with `hoonc`, plus Rust jets for
the `%i754` array arms built on [`sdfloat`/`sdblas`](https://github.com/sigilante/RustFloat),
so a NockApp computes exactly the bits the Hoon (and the vere C jets) produce.

```
crates/lagoon-jets/     the jets (`LAGOON_HOT`), `hot_state()` = hoon-138 built-ins + these
  src/bin/lagoon-kick   test driver: fires a hoonc --arbitrary trap under the jets
hoon/                   generated: lib/ sur/ copies of the desk with `..part` -> `..ut`,
                        the test files, and run-tests.hoon (the suite as lazy traps)
scripts/sync-hoon.sh    regenerates hoon/ from ../desk and ../../libmath/desk
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

## Using the jets in a NockApp

```rust
let hot = lagoon_jets::hot_state();          // URBIT_HOT_STATE ++ LAGOON_HOT
boot::setup(&kernel_jam, cli, &hot, "my-app", None).await?;
```

`Hot::init` registers only what it is given: passing `LAGOON_HOT` alone
leaves `add`, `dec` and the bit operations running as raw Nock.

## Jets vs Hoon under NockVM

Any input a jet does not handle punts to the Nock, so jets are never a
correctness risk, only a speed one. Where the vere C jets and the Hoon
disagree, the Hoon wins here (test mode enforces it); see `PORTING.md` for
the cases the port surfaced.
