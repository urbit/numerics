//! Hot-state registration for the lagoon jets.
//!
//! Path: `[k.138 one two tri qua pen non lagoon <arm>]`. `non` is the
//! `~%  %non  ..ut  ~` chapter wrapping `/lib/lagoon` (parent `..ut`, so it
//! nests under hoon-138's `%pen`), `lagoon` is the `~/  %lagoon` door, and
//! each arm is its `~/  %<name>` hint.

use nockvm::jets::hot::{HotEntry, K_138, URBIT_HOT_STATE};
use nockvm::jets::Jet;
use either::Either::Left;

use crate::jets;

macro_rules! lagoon_arm {
    ($name:literal, $jet:expr) => {
        (
            &[
                K_138,
                Left(b"one"),
                Left(b"two"),
                Left(b"tri"),
                Left(b"qua"),
                Left(b"pen"),
                Left(b"non"),
                Left(b"lagoon"),
                Left($name),
            ],
            1,
            $jet as Jet,
        )
    };
}

/// Every lagoon jet, for `boot::setup(.., &LAGOON_HOT, ..)` or
/// `nockapp::utils::create_context`.
pub const LAGOON_HOT: &[HotEntry] = &[
    lagoon_arm!(b"add-rays", jets::add_rays),
    lagoon_arm!(b"mmul", jets::mmul),
    lagoon_arm!(b"sub-rays", jets::elem::sub),
    lagoon_arm!(b"mul-rays", jets::elem::mul),
    lagoon_arm!(b"div-rays", jets::elem::div),
    lagoon_arm!(b"mod-rays", jets::elem::mod_rays),
    lagoon_arm!(b"add-scal", jets::elem::add_scalar),
    lagoon_arm!(b"sub-scal", jets::elem::sub_scalar),
    lagoon_arm!(b"mul-scal", jets::elem::mul_scalar),
    lagoon_arm!(b"div-scal", jets::elem::div_scalar),
    lagoon_arm!(b"mod-scal", jets::elem::mod_scalar),
    lagoon_arm!(b"gth", jets::cmp::gth),
    lagoon_arm!(b"gte", jets::cmp::gte),
    lagoon_arm!(b"lth", jets::cmp::lth),
    lagoon_arm!(b"lte", jets::cmp::lte),
    lagoon_arm!(b"abs", jets::cmp::abs),
    lagoon_arm!(b"max", jets::reduce::max),
    lagoon_arm!(b"min", jets::reduce::min),
    lagoon_arm!(b"argmax", jets::reduce::argmax),
    lagoon_arm!(b"argmin", jets::reduce::argmin),
    lagoon_arm!(b"cumsum", jets::reduce::cumsum),
    lagoon_arm!(b"dot", jets::reduce::dot),
    lagoon_arm!(b"trace", jets::reduce::trace),
    lagoon_arm!(b"transpose", jets::shape::transpose),
    lagoon_arm!(b"diag", jets::shape::diag),
    lagoon_arm!(b"ravel", jets::shape::ravel),
    lagoon_arm!(b"range", jets::shape::range),
    lagoon_arm!(b"linspace", jets::shape::linspace),
];

/// The hoon-138 built-in jets, the `rh`/`rs`/`rd`/`rq` door jets, and the
/// lagoon jets: what a NockApp should hand to `boot::setup` /
/// `create_context`. `Hot::init` registers only what it is given, so passing
/// `LAGOON_HOT` alone would leave `add`, `dec`, and the bit operations
/// running as raw Nock.
pub fn hot_state() -> Vec<HotEntry> {
    let mut v = Vec::with_capacity(URBIT_HOT_STATE.len() + hoon_float_jets::HOON_FLOAT_HOT.len() + LAGOON_HOT.len());
    v.extend_from_slice(URBIT_HOT_STATE);
    // Diagnostic: leave the float doors to the Nock, to diff jet vs Hoon.
    if std::env::var_os("HOON_FLOAT_JET_DISABLE").is_none() {
        v.extend_from_slice(hoon_float_jets::HOON_FLOAT_HOT);
    }
    v.extend_from_slice(LAGOON_HOT);
    v
}
