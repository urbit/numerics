//! Reductions: max/min/argmax/argmin/cumsum/dot/trace.
//! C: `u3qi_la_{max,min,argmax,argmin,cumsum,dot,trace}_i754`.
//!
//! The Hoon reduces the ravel with `reel` (a right fold). For `max`/`min` the
//! fold's initial accumulator is the *first* element (the `|:` sample
//! default), and the step keeps the element when it compares strictly
//! greater (less) than the accumulator, so ties and NaNs keep the
//! accumulator. For `cumsum` the initial accumulator is `+0` and each step is
//! `add(x_i, acc)` from the last element to the first, which is also the
//! order the C kernel uses. `dot` is `cumsum` of the elementwise product,
//! `trace` is `cumsum` of the diagonal; both box the scalar with
//! `+scalar-to-ray` (all-ones shape of the same rank, input bloq/kind/tail).

use nockvm::interpreter::Context;
use nockvm::jets::{JetErr, Result};
use nockvm::noun::{Noun, D};
use sdblas::{BlasFloat, Round, F128, F16, F32, F64};

use crate::by_bloq;
use crate::jets::{from_raw, same_meta, to_raw, two_rays};
use crate::ray::{self, Ray, KIND_I754};

/// `LAGOON_JET_DISABLE=max,argmax,..` makes the named jets punt, to compare
/// the pure Hoon against the jet on the same run.
fn disabled(name: &str) -> bool {
    std::env::var("LAGOON_JET_DISABLE")
        .map(|v| v.split(',').any(|n| n == name))
        .unwrap_or(false)
}

fn one_ray(context: &mut Context, subject: Noun) -> std::result::Result<(Ray, Round, usize), JetErr> {
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let a = ray::parse(sam, &space)?;
    let rnd = ray::rounding(subject, &space)?;
    let rank = ray::list_atoms(a.shape_noun, &space)?.len();
    if a.kind != KIND_I754 {
        return Err(JetErr::Punt);
    }
    Ok((a, rnd, rank))
}

/// `+scalar-to-ray`: box one element with an all-ones shape of `rank`.
fn scalar_ray(context: &mut Context, a: &Ray, rank: usize, v: u128) -> Noun {
    let shape = vec![1u64; rank];
    let meta = ray::build_meta(context, &shape, a.bloq, a.kind, a.tail);
    let data = ray::pack(context, a.width(), &[v]);
    ray::build(context, meta, data)
}

/// `reel` with the first element as the initial accumulator, keeping `x_i`
/// when `gth(x_i, acc)` (`want_max`) or `lth(x_i, acc)` is true.
fn fold_extreme<F: BlasFloat>(x: &[u128], want_max: bool) -> u128 {
    let mut acc: F = from_raw(x[0]);
    for &v in x.iter().rev() {
        let xi: F = from_raw(v);
        // gth(a, b) is lth(b, a) in `++fl`; both are false on NaN.
        let take = if want_max { acc.lt(xi) } else { xi.lt(acc) };
        if take {
            acc = xi;
        }
    }
    to_raw(acc)
}

/// `reel` sum from `+0`, `add(x_i, acc)` from the last element down.
fn fold_sum<F: BlasFloat>(x: &[u128], r: Round) -> u128 {
    let mut acc = F::ZERO;
    for &v in x.iter().rev() {
        let xi: F = from_raw(v);
        acc = xi.add(acc, r);
    }
    to_raw(acc)
}

fn products<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
    x.iter()
        .zip(y)
        .map(|(&a, &b)| to_raw(from_raw::<F>(a).mul(from_raw::<F>(b), r)))
        .collect()
}

fn extreme(context: &mut Context, subject: Noun, want_max: bool) -> Result {
    let (a, _rnd, rank) = one_ray(context, subject)?;
    if a.len == 0 {
        return Err(JetErr::Punt); // the Hoon crashes on `-:(ravel a)`
    }
    let space = context.stack.noun_space();
    let x = ray::elems(&a, &space);
    let v = by_bloq!(a.bloq, fold_extreme(&x, want_max));
    Ok(scalar_ray(context, &a, rank, v))
}

fn arg_extreme(context: &mut Context, subject: Noun, want_max: bool) -> Result {
    let (a, _rnd, _rank) = one_ray(context, subject)?;
    if a.len == 0 {
        return Err(JetErr::Punt);
    }
    let space = context.stack.noun_space();
    let x = ray::elems(&a, &space);
    let v = by_bloq!(a.bloq, fold_extreme(&x, want_max));
    // `find` on the ravel: first index whose raw bits equal the extreme.
    match x.iter().position(|&e| e == v) {
        Some(i) => Ok(D(i as u64)),
        None => Err(JetErr::Punt),
    }
}

/// `+max`: `~/  %max`, `a=ray -> ray`.
pub fn max(context: &mut Context, subject: Noun) -> Result {
    ray::trace("max");
    if disabled("max") {
        return Err(JetErr::Punt);
    }
    extreme(context, subject, true)
}

/// `+min`: `~/  %min`, `a=ray -> ray`.
pub fn min(context: &mut Context, subject: Noun) -> Result {
    ray::trace("min");
    if disabled("min") {
        return Err(JetErr::Punt);
    }
    extreme(context, subject, false)
}

/// `+argmax`: `~/  %argmax`, `a=ray -> @ud`.
pub fn argmax(context: &mut Context, subject: Noun) -> Result {
    ray::trace("argmax");
    if disabled("argmax") {
        return Err(JetErr::Punt);
    }
    arg_extreme(context, subject, true)
}

/// `+argmin`: `~/  %argmin`, `a=ray -> @ud`.
pub fn argmin(context: &mut Context, subject: Noun) -> Result {
    ray::trace("argmin");
    if disabled("argmin") {
        return Err(JetErr::Punt);
    }
    arg_extreme(context, subject, false)
}

/// `+cumsum`: `~/  %cumsum`, `a=ray -> ray` (a scalar ray, see module doc).
pub fn cumsum(context: &mut Context, subject: Noun) -> Result {
    ray::trace("cumsum");
    if disabled("cumsum") {
        return Err(JetErr::Punt);
    }
    let (a, rnd, rank) = one_ray(context, subject)?;
    let space = context.stack.noun_space();
    let x = ray::elems(&a, &space);
    let v = by_bloq!(a.bloq, fold_sum(&x, rnd));
    Ok(scalar_ray(context, &a, rank, v))
}

/// `+dot`: `~/  %dot`, `[a=ray b=ray] -> ray`: `(cumsum (mul a b))`.
pub fn dot(context: &mut Context, subject: Noun) -> Result {
    ray::trace("dot");
    if disabled("dot") {
        return Err(JetErr::Punt);
    }
    let (a, b, rnd) = two_rays(context, subject)?;
    let space = context.stack.noun_space();
    if a.kind != KIND_I754 || !same_meta(&a, &b, &space) || a.len == 0 {
        return Err(JetErr::Punt);
    }
    let rank = ray::list_atoms(a.shape_noun, &space)?.len();
    let x = ray::elems(&a, &space);
    let y = ray::elems(&b, &space);
    let p = by_bloq!(a.bloq, products(&x, &y, rnd));
    let v = by_bloq!(a.bloq, fold_sum(&p, rnd));
    Ok(scalar_ray(context, &a, rank, v))
}

/// `+trace`: `~/  %trace`, `a=ray -> ray`: `(cumsum (diag a))`, square 2-D only.
pub fn trace(context: &mut Context, subject: Noun) -> Result {
    ray::trace("trace");
    if disabled("trace") {
        return Err(JetErr::Punt);
    }
    let (a, rnd, _rank) = one_ray(context, subject)?;
    let space = context.stack.noun_space();
    let dims = ray::list_atoms(a.shape_noun, &space)?;
    if dims.len() != 2 || dims[0] != dims[1] || dims[0] == 0 {
        return Err(JetErr::Punt);
    }
    let n = dims[0] as usize;
    let x = ray::elems(&a, &space);
    let d: Vec<u128> = (0..n).map(|i| x[i * n + i]).collect();
    let v = by_bloq!(a.bloq, fold_sum(&d, rnd));
    // `diag` yields meta `[~[n 1] bloq kind tail]`; `+scalar-to-ray` of that is
    // shape `~[1 1]` with the input's bloq/kind/tail.
    Ok(scalar_ray(context, &a, 2, v))
}
