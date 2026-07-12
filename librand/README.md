# `/lib/rand` for Urbit

Deterministic pseudo-random number generation for numerical work (Monte
Carlo, ML init, shuffles, simulation). **Not cryptographic** — entropy
acquisition is Arvo's job (`eny`), cryptographic randomness is Zuse's job.

Full design is in `rand-spec.md` (repo root). Hoon reference implementation
first, jets follow (see `rand-spec.md` section 11).

## Status (milestone 1 of 9, `rand-spec.md` section 13)

Done:

- `++split-mix` — SplitMix64 (`+next`, `+split`, `+finalize`), KAT-checked
  against Vigna's reference C.
- `++seed` — `+from-atom`, `+from-eny`, `+mix` (the two-word compression
  primitive `+fork` will use in milestone 2).

Not yet implemented: `++philox` (primary engine), `++pcg`, `++fork`,
`++gen` facade, `++uni`, `++dist`, `++sample`, the Saloon `+rand-ray`
extension. See `NEXT-STEPS.md` and `rand-spec.md` section 13 for the full
milestone order.

## Layout

```
librand/
  README.md
  NEXT-STEPS.md
  desk/lib/rand.hoon         :: the library
  desk/tests/lib/rand.hoon   :: -test %/tests/lib/rand ~
```
