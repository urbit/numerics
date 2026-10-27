//! Elementwise comparisons and abs.
//! C: `u3qi_la_{gth,gte,lth,lte,abs}_i754`.
//!
//! Semantics come from the Hoon, which for `%i754` is `+bin-op` over
//! `?:((~(gth rs rnd) a b) .1 .0)` per width, and for `+abs` is `+el-wise-op`
//! over `?:((~(gte rs rnd) b .0) b (~(mul rs rnd) b .-1))`. The `++ff`
//! comparisons are `(fall (lth:fl ..) |)`: any NaN operand compares false,
//! `-0 == +0`, and the ordering is IEEE, i.e. exactly sdfloat's `lt`/`le`.
//! The result element is the kind's `1.0` or `+0.0`.
//!
//! `+abs` therefore differs from the C `f32_abs` (a sign-bit clear) in two
//! places that the Hoon decides: `-0` satisfies `gte b 0` and comes back as
//! `-0` unchanged, and a NaN fails it and comes back as the canonical NaN of
//! `mul`. Both are reproduced here rather than sdfloat's `abs`.

use nockvm::interpreter::Context;
use nockvm::jets::{JetErr, Result};
use nockvm::noun::Noun;
use sdblas::{BlasFloat, Round, F128, F16, F32, F64};

use super::{from_raw, same_meta, to_raw, two_rays};
use crate::by_bloq;
use crate::ray::{self, KIND_I754};

/// One of the four orderings, as `(a, b) -> bool` in IEEE terms with NaN
/// false, matching `++ff`.
#[derive(Clone, Copy)]
enum Ord {
    /// `a > b` = `lth:fl b a`.
    Gth,
    /// `a >= b` = `lte:fl b a` = `!(a < b)`, false on NaN.
    Gte,
    /// `a < b`.
    Lth,
    /// `a <= b` = `!(b < a)`, false on NaN.
    Lte,
}

fn cmp_w<F: BlasFloat>(x: &[u128], y: &[u128], o: Ord) -> Vec<u128> {
    x.iter()
        .zip(y)
        .map(|(&a, &b)| {
            let (a, b): (F, F) = (from_raw(a), from_raw(b));
            let t = match o {
                Ord::Gth => b.lt(a),
                Ord::Gte => b.le(a),
                Ord::Lth => a.lt(b),
                Ord::Lte => a.le(b),
            };
            to_raw(if t { F::ONE } else { F::ZERO })
        })
        .collect()
}

fn compare(context: &mut Context, subject: Noun, o: Ord) -> Result {
    let (a, b, _rnd) = two_rays(context, subject)?;
    let space = context.stack.noun_space();
    if a.kind != KIND_I754 || !same_meta(&a, &b, &space) {
        return Err(JetErr::Punt);
    }
    let x = ray::elems(&a, &space);
    let y = ray::elems(&b, &space);
    let out = by_bloq!(a.bloq, cmp_w(&x, &y, o));
    let data = ray::pack(context, a.width(), &out);
    Ok(ray::build(context, a.meta, data))
}

/// `+gth`: `~/  %gth`, `[a=ray b=ray] -> ray`.
pub fn gth(context: &mut Context, subject: Noun) -> Result {
    ray::trace("gth");
    compare(context, subject, Ord::Gth)
}

/// `+gte`: `~/  %gte`, `[a=ray b=ray] -> ray`.
pub fn gte(context: &mut Context, subject: Noun) -> Result {
    ray::trace("gte");
    compare(context, subject, Ord::Gte)
}

/// `+lth`: `~/  %lth`, `[a=ray b=ray] -> ray`.
pub fn lth(context: &mut Context, subject: Noun) -> Result {
    ray::trace("lth");
    compare(context, subject, Ord::Lth)
}

/// `+lte`: `~/  %lte`, `[a=ray b=ray] -> ray`.
pub fn lte(context: &mut Context, subject: Noun) -> Result {
    ray::trace("lte");
    compare(context, subject, Ord::Lte)
}

fn abs_w<F: BlasFloat>(x: &[u128], r: Round) -> Vec<u128> {
    x.iter()
        .map(|&v| {
            let b: F = from_raw(v);
            // `?:((gte b .0) b (mul b .-1))`
            to_raw(if F::ZERO.le(b) { b } else { b.mul(F::NEG_ONE, r) })
        })
        .collect()
}

/// `+abs`: `~/  %abs`, `a=ray -> ray`.
pub fn abs(context: &mut Context, subject: Noun) -> Result {
    ray::trace("abs");
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let a = ray::parse(sam, &space)?;
    let rnd = ray::rounding(subject, &space)?;
    if a.kind != KIND_I754 {
        return Err(JetErr::Punt);
    }
    let x = ray::elems(&a, &space);
    let out = by_bloq!(a.bloq, abs_w(&x, rnd));
    let data = ray::pack(context, a.width(), &out);
    Ok(ray::build(context, a.meta, data))
}
