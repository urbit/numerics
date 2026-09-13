//! The jets. One function per `+la` arm; each dispatches on the ray's bloq
//! to a width-generic worker over [`sdblas::BlasFloat`].
//!
//! Layout: shared helpers and `add-rays`/`mmul` here; the other arms in
//! submodules by family. Every arm listed in `hot.rs` has a function here;
//! an unported one is a stub returning `JetErr::Punt`.

pub mod cmp;
pub mod elem;
pub mod reduce;
pub mod shape;

use nockvm::interpreter::Context;
use nockvm::jets::util::slot;
use nockvm::jets::{JetErr, Result};
use nockvm::noun::{Noun, D};
use sdblas::{BlasFloat, Round, Trans, F128, F16, F32, F64};
use sdfloat::Format;

use crate::ray::{self, Ray, KIND_I754};

/// Run `f::<F>()` for the width named by `bloq` (4..=7).
#[macro_export]
macro_rules! by_bloq {
    ($bloq:expr, $f:ident($($arg:expr),*)) => {
        match $bloq {
            4 => $f::<F16>($($arg),*),
            5 => $f::<F32>($($arg),*),
            6 => $f::<F64>($($arg),*),
            7 => $f::<F128>($($arg),*),
            _ => return Err(JetErr::Punt),
        }
    };
}

pub(crate) fn from_raw<F: Format>(v: u128) -> F {
    F::from_raw(v)
}
pub(crate) fn to_raw<F: Format>(v: F) -> u128 {
    v.raw()
}

pub(crate) fn two_rays(context: &mut Context, subject: Noun) -> std::result::Result<(Ray, Ray, Round), JetErr> {
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let a = ray::parse(slot(sam, 2, &space)?, &space)?;
    let b = ray::parse(slot(sam, 3, &space)?, &space)?;
    let rnd = ray::rounding(subject, &space)?;
    Ok((a, b, rnd))
}

pub(crate) fn same_meta(a: &Ray, b: &Ray, space: &nockvm::noun::NounSpace) -> bool {
    a.bloq == b.bloq
        && a.kind == b.kind
        && a.len == b.len
        && nockvm::ext::noun_equality(a.meta.in_space(space), b.meta.in_space(space))
}

// ---------------------------------------------------------------- add-rays

fn add_rays_w<F: BlasFloat>(x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
    // C: ?axpy(len, ONE, x, 1, y, 1): y <- 1*x + y.
    let xs: Vec<F> = x.iter().map(|&v| from_raw(v)).collect();
    let mut ys: Vec<F> = y.iter().map(|&v| from_raw(v)).collect();
    sdblas::axpy(xs.len() as u64, F::ONE, &xs, 1, &mut ys, 1, r);
    ys.into_iter().map(to_raw).collect()
}

/// `+add`: `~/  %add-rays`, `[a=ray b=ray] -> ray`.
pub fn add_rays(context: &mut Context, subject: Noun) -> Result {
    ray::trace("add-rays");
    let (a, b, rnd) = two_rays(context, subject)?;
    let space = context.stack.noun_space();
    if a.kind != KIND_I754 || !same_meta(&a, &b, &space) {
        return Err(JetErr::Punt);
    }
    let x = ray::elems(&a, &space);
    let y = ray::elems(&b, &space);
    let mut out = by_bloq!(a.bloq, add_rays_w(&x, &y, rnd));
    if std::env::var_os("LAGOON_JET_SABOTAGE").is_some() {
        // Harness self-check: a deliberately wrong result must be caught by
        // NOCK_TEST_JETS. Never set this outside that check.
        if let Some(e) = out.first_mut() {
            *e ^= 1;
        }
    }
    let data = ray::pack(context, a.width(), &out);
    Ok(ray::build(context, a.meta, data))
}

// -------------------------------------------------------------------- mmul

fn mmul_w<F: BlasFloat>(m: u64, n: u64, p: u64, x: &[u128], y: &[u128], r: Round) -> Vec<u128> {
    // C: ?gemm('N','N', M, N, P, ONE, x, N, y, P, ZERO, out(zeroed), P).
    let xs: Vec<F> = x.iter().map(|&v| from_raw(v)).collect();
    let ys: Vec<F> = y.iter().map(|&v| from_raw(v)).collect();
    let mut out = vec![F::ZERO; (m * p) as usize];
    sdblas::gemm(Trans::NoTrans, Trans::NoTrans, m, n, p, F::ONE, &xs, n, &ys, p, F::ZERO, &mut out, p, r);
    out.into_iter().map(to_raw).collect()
}

/// `+mmul`: `~/  %mmul`, `[a=ray b=ray] -> ray`; result meta is
/// `[~[M P] bloq kind ~]` (the Hoon builds it with `zeros`, tail `~`).
pub fn mmul(context: &mut Context, subject: Noun) -> Result {
    ray::trace("mmul");
    let (a, b, rnd) = two_rays(context, subject)?;
    let space = context.stack.noun_space();
    if a.kind != KIND_I754 || b.kind != KIND_I754 || a.bloq != b.bloq {
        return Err(JetErr::Punt);
    }
    let sa = ray::list_atoms(a.shape_noun, &space)?;
    let sb = ray::list_atoms(b.shape_noun, &space)?;
    if sa.len() != 2 || sb.len() != 2 || sa[1] != sb[0] {
        return Err(JetErr::Punt);
    }
    let (m, n, p) = (sa[0], sa[1], sb[1]);
    let x = ray::elems(&a, &space);
    let y = ray::elems(&b, &space);
    let out = by_bloq!(a.bloq, mmul_w(m, n, p, &x, &y, rnd));
    let data = ray::pack(context, a.width(), &out);
    let meta = ray::build_meta(context, &[m, p], a.bloq, a.kind, D(0));
    Ok(ray::build(context, meta, data))
}
