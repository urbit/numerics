//! NockVM jets for `/lib/lagoon`.
//!
//! Each jet mirrors one arm of the `+la` door for `%i754` rays, computing with
//! [`sdfloat`]/[`sdblas`] so the bits equal what the Hoon produces. Any input
//! a jet does not handle (another kind, an unexpected shape, an invalid ray)
//! punts: the interpreter runs the Nock instead, which is always correct.
//! Jets never bail where the Hoon would succeed.
//!
//! Registration: the Hoon is `~%  %non  ..ut  ~` around `/lib/lagoon`, the
//! door is `~/  %lagoon`, and each arm is `~/  %<name>`; see [`hot::LAGOON_HOT`].

pub mod hot;
pub mod jets;
pub mod ray;

pub use hot::{hot_state, LAGOON_HOT};
