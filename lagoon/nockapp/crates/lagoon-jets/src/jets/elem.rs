//! Elementwise arithmetic on two rays and ray-scalar arithmetic.
//! C: `u3qi_la_{sub,mul,div,mod,adds,subs,muls,divs,mods}_i754`.
//!
//! The Hoon is `+bin-op` over `+fun-scalar`: each element pair goes through
//! `~(op rX rnd)` for the width, independently; the scalar arms `+fill` a
//! ray with `n` and call the ray arm. So every arm here is one width-generic
//! elementwise loop with an explicit `Round`.

use nockvm::interpreter::Context;
use nockvm::jets::util::slot;
use nockvm::jets::{JetErr, Result};
use nockvm::noun::Noun;
use sdblas::{BlasFloat, Round, F128, F16, F32, F64};

use crate::by_bloq;
use crate::jets::{from_raw, same_meta, to_raw, two_rays};
use crate::ray::{self, KIND_I754};

/// The Hoon `%mod` for `%i754` (`+fun-scalar`, since numerics #78):
/// `a - b * san(need(toi_r(a / b)))`, every step rounded in the door mode
/// `r`, the quotient rounded to an integer in `r` too (`+toi`). A non-finite
/// quotient (zero divisor, NaN, inf) makes `toi` give `~` and `need` crash:
/// `None` here, and the jet punts so the Nock crashes. `san (toi q)` yields
/// `+0` for a quotient that rounds to zero, whatever its sign.
pub(crate) fn hoon_mod<F: BlasFloat>(a: F, b: F, r: Round) -> Option<F> {
    let q = a.div(b, r);
    if q.is_nan() || q.is_inf() {
        return None;
    }
    let mut t = q.round_to_int(r);
    if t.is_zero() {
        t = F::ZERO;
    }
    Some(a.sub(b.mul(t, r), r))
}

fn map2<F: BlasFloat>(x: &[u128], y: &[u128], r: Round, op: fn(F, F, Round) -> F) -> Vec<u128> {
    x.iter()
        .zip(y)
        .map(|(&p, &q)| to_raw(op(from_raw::<F>(p), from_raw::<F>(q), r)))
        .collect()
}

fn w_sub<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
    map2::<F>(x, y, r, |a, b, r| a.sub(b, r))
}
fn w_mul<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
    map2::<F>(x, y, r, |a, b, r| a.mul(b, r))
}
fn w_div<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
    map2::<F>(x, y, r, |a, b, r| a.div(b, r))
}
fn w_mod<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> std::result::Result<Vec<u128>, JetErr> {
    x.iter()
        .zip(y)
        .map(|(&p, &q)| hoon_mod(from_raw::<F>(p), from_raw::<F>(q), r).map(to_raw).ok_or(JetErr::Punt))
        .collect()
}

#[derive(Clone, Copy)]
enum Op {
    Sub,
    Mul,
    Div,
    Mod,
}

fn apply(bloq: u32, op: Op, x: &[u128], y: &[u128], r: Round) -> std::result::Result<Vec<u128>, JetErr> {
    match op {
        Op::Sub => Ok(by_bloq!(bloq, w_sub(x, y, r))),
        Op::Mul => Ok(by_bloq!(bloq, w_mul(x, y, r))),
        Op::Div => Ok(by_bloq!(bloq, w_div(x, y, r))),
        Op::Mod => by_bloq!(bloq, w_mod(x, y, r)),
    }
}

/// `[a=ray b=ray] -> ray` over `+bin-op`: same meta, elementwise.
fn rays(context: &mut Context, subject: Noun, op: Op) -> Result {
    let (a, b, rnd) = two_rays(context, subject)?;
    let space = context.stack.noun_space();
    if a.kind != KIND_I754 || !same_meta(&a, &b, &space) {
        return Err(JetErr::Punt);
    }
    let x = ray::elems(&a, &space);
    let y = ray::elems(&b, &space);
    let out = apply(a.bloq, op, &x, &y, rnd)?;
    let data = ray::pack(context, a.width(), &out);
    Ok(ray::build(context, a.meta, data))
}

/// `[a=ray n=@] -> ray`: the Hoon fills a ray of `a`'s meta with `n` (so `n`
/// is truncated to the element width by `+rep`) and calls the ray arm.
fn scalar(context: &mut Context, subject: Noun, op: Op) -> Result {
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let a = ray::parse(slot(sam, 2, &space)?, &space)?;
    let n = slot(sam, 3, &space)?
        .in_space(&space)
        .as_atom()
        .map_err(|_| JetErr::Punt)?
        .to_le_bytes();
    let rnd = ray::rounding(subject, &space)?;
    if a.kind != KIND_I754 {
        return Err(JetErr::Punt);
    }
    let w = a.width();
    let mut nb = [0u8; 16];
    for (i, byte) in n.iter().take(w).enumerate() {
        nb[i] = *byte;
    }
    let nv = u128::from_le_bytes(nb);
    let x = ray::elems(&a, &space);
    let y = vec![nv; x.len()];
    let out = apply(a.bloq, op, &x, &y, rnd)?;
    let data = ray::pack(context, w, &out);
    Ok(ray::build(context, a.meta, data))
}

/// `+sub`: `~/  %sub-rays`, `[a=ray b=ray] -> ray`.
pub fn sub(context: &mut Context, subject: Noun) -> Result {
    ray::trace("sub-rays");
    rays(context, subject, Op::Sub)
}

/// `+mul`: `~/  %mul-rays`, `[a=ray b=ray] -> ray`.
pub fn mul(context: &mut Context, subject: Noun) -> Result {
    ray::trace("mul-rays");
    rays(context, subject, Op::Mul)
}

/// `+div`: `~/  %div-rays`, `[a=ray b=ray] -> ray`.
pub fn div(context: &mut Context, subject: Noun) -> Result {
    ray::trace("div-rays");
    rays(context, subject, Op::Div)
}

/// `+mod`: `~/  %mod-rays`, `[a=ray b=ray] -> ray`. Per the Hoon: quotient
/// rounded in the door mode; punts (so the Nock crashes) on a zero divisor.
pub fn mod_rays(context: &mut Context, subject: Noun) -> Result {
    ray::trace("mod-rays");
    rays(context, subject, Op::Mod)
}

/// `+add-scalar`: `~/  %add-scal`, `[a=ray n=@] -> ray`.
pub fn add_scalar(context: &mut Context, subject: Noun) -> Result {
    ray::trace("add-scal");
    // `a + fill(n)`: reuse the add kernel's semantics (IEEE add is
    // commutative, so the operand order the C `?axpy` uses is immaterial).
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let a = ray::parse(slot(sam, 2, &space)?, &space)?;
    let n = slot(sam, 3, &space)?
        .in_space(&space)
        .as_atom()
        .map_err(|_| JetErr::Punt)?
        .to_le_bytes();
    let rnd = ray::rounding(subject, &space)?;
    if a.kind != KIND_I754 {
        return Err(JetErr::Punt);
    }
    let w = a.width();
    let mut nb = [0u8; 16];
    for (i, byte) in n.iter().take(w).enumerate() {
        nb[i] = *byte;
    }
    let nv = u128::from_le_bytes(nb);
    let x = ray::elems(&a, &space);
    let y = vec![nv; x.len()];
    fn w_add<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
        map2::<F>(x, y, r, |a, b, r| a.add(b, r))
    }
    let out = by_bloq!(a.bloq, w_add(&x, &y, rnd));
    let data = ray::pack(context, w, &out);
    Ok(ray::build(context, a.meta, data))
}

/// `+sub-scalar`: `~/  %sub-scal`, `[a=ray n=@] -> ray`.
pub fn sub_scalar(context: &mut Context, subject: Noun) -> Result {
    ray::trace("sub-scal");
    scalar(context, subject, Op::Sub)
}

/// `+mul-scalar`: `~/  %mul-scal`, `[a=ray n=@] -> ray`.
pub fn mul_scalar(context: &mut Context, subject: Noun) -> Result {
    ray::trace("mul-scal");
    scalar(context, subject, Op::Mul)
}

/// `+div-scalar`: `~/  %div-scal`, `[a=ray n=@] -> ray`.
pub fn div_scalar(context: &mut Context, subject: Noun) -> Result {
    ray::trace("div-scal");
    scalar(context, subject, Op::Div)
}

/// `+mod-scalar`: `~/  %mod-scal`, `[a=ray n=@] -> ray`.
pub fn mod_scalar(context: &mut Context, subject: Noun) -> Result {
    ray::trace("mod-scal");
    scalar(context, subject, Op::Mod)
}
