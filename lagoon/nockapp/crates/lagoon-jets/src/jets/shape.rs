//! Shape/construction: transpose/diag/ravel/range/linspace.
//! C: `u3qi_la_{transpose,diag,ravel_i754,range_i754,linspace_i754}`.
//!
//! `transpose`, `diag`, and `ravel` are pure bit movement and serve every
//! kind the codec accepts; `range` and `linspace` are `%i754` by definition
//! (the Hoon forces `kind.meta` to `%i754` and dispatches on bloq 4..=7).

use nockvm::ext::AtomExt;
use nockvm::interpreter::Context;
use nockvm::jets::util::slot;
use nockvm::jets::{JetErr, Result};
use nockvm::noun::{Atom, Noun, NounSpace, D, T};
use sdblas::{BlasFloat, Round, F128, F16, F32, F64};

use crate::by_bloq;
use crate::jets::{from_raw, to_raw};
use crate::ray::{self, KIND_I754};

/// Guard against the Hoon's unbounded `+range` loop (d = 0, NaN, or a step
/// too small to move): past this many elements the jet punts and lets the
/// Nock do whatever the Nock does.
const RANGE_CAP: usize = 1 << 24;

fn one_ray(context: &mut Context, subject: Noun) -> std::result::Result<ray::Ray, JetErr> {
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    ray::parse(sam, &space)
}

/// `[shape bloq kind tail]` of a bare `meta` sample (no data to check).
fn parse_meta(meta: Noun, space: &NounSpace) -> std::result::Result<(Vec<u64>, u32, Noun), JetErr> {
    let m = meta.in_space(space).as_cell().map_err(|_| JetErr::Punt)?;
    let shape = ray::list_atoms(m.head().noun(), space)?;
    let m2 = m.tail().as_cell().map_err(|_| JetErr::Punt)?;
    let bloq = m2.head().as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
    let m3 = m2.tail().as_cell().map_err(|_| JetErr::Punt)?;
    let tail = m3.tail().noun();
    if !(4..=7).contains(&bloq) {
        return Err(JetErr::Punt);
    }
    Ok((shape, bloq as u32, tail))
}

fn atom_u128(n: Noun, space: &NounSpace) -> std::result::Result<u128, JetErr> {
    let a = n.in_space(space).as_atom().map_err(|_| JetErr::Punt)?;
    if a.bit_size() > 128 {
        return Err(JetErr::Punt);
    }
    let mut b = [0u8; 16];
    let v = a.to_le_bytes();
    b[..v.len().min(16)].copy_from_slice(&v[..v.len().min(16)]);
    Ok(u128::from_le_bytes(b))
}

// --------------------------------------------------------------- transpose

/// `+transpose`: `~/  %transpose`, `a=ray -> ray`. 2-D only; the result is
/// `[~[c r] bloq kind ~]` (tail `~`, from the Hoon's `zeros`) with
/// `out[j][i] = a[i][j]`.
pub fn transpose(context: &mut Context, subject: Noun) -> Result {
    ray::trace("transpose");
    let a = one_ray(context, subject)?;
    let space = context.stack.noun_space();
    let shape = ray::list_atoms(a.shape_noun, &space)?;
    if shape.len() != 2 {
        return Err(JetErr::Punt);
    }
    let (r, c) = (shape[0] as usize, shape[1] as usize);
    let x = ray::elems(&a, &space);
    let mut out = vec![0u128; r * c];
    for i in 0..r {
        for j in 0..c {
            out[j * r + i] = x[i * c + j];
        }
    }
    let data = ray::pack(context, a.width(), &out);
    let meta = ray::build_meta(context, &[c as u64, r as u64], a.bloq, a.kind, D(0));
    Ok(ray::build(context, meta, data))
}

// -------------------------------------------------------------------- diag

/// `+diag`: `~/  %diag`, `a=ray -> ray`. Square 2-D only; the main diagonal
/// as an `n x 1` ray with meta `[~[n 1] bloq kind tail]` (tail kept).
pub fn diag(context: &mut Context, subject: Noun) -> Result {
    ray::trace("diag");
    let a = one_ray(context, subject)?;
    let space = context.stack.noun_space();
    let shape = ray::list_atoms(a.shape_noun, &space)?;
    if shape.len() != 2 || shape[0] != shape[1] || shape[0] == 0 {
        return Err(JetErr::Punt);
    }
    let n = shape[0] as usize;
    let x = ray::elems(&a, &space);
    let out: Vec<u128> = (0..n).map(|i| x[i * n + i]).collect();
    let data = ray::pack(context, a.width(), &out);
    let meta = ray::build_meta(context, &[n as u64, 1], a.bloq, a.kind, a.tail);
    Ok(ray::build(context, meta, data))
}

// ------------------------------------------------------------------- ravel

/// `+ravel`: `~/  %ravel`, `a=ray -> (list @)`: `(snip (rip bloq data))`.
/// Works on the data atom alone, exactly as the Hoon does (no `+check`):
/// every `2^bloq`-bit block up to the atom's top bit, minus the last one.
/// Serves any kind and any bloq up to 7.
pub fn ravel(context: &mut Context, subject: Noun) -> Result {
    ray::trace("ravel");
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let c = sam.in_space(&space).as_cell().map_err(|_| JetErr::Punt)?;
    let meta = c.head().noun();
    let data = c.tail().as_atom().map_err(|_| JetErr::Punt)?;
    let m = meta.in_space(&space).as_cell().map_err(|_| JetErr::Punt)?;
    let m2 = m.tail().as_cell().map_err(|_| JetErr::Punt)?;
    let bloq = m2.head().as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
    if bloq > 7 {
        return Err(JetErr::Punt);
    }
    let bits = data.bit_size();
    if bits == 0 {
        // (rip bloq 0) is ~ and (snip ~) crashes; let the Nock crash.
        return Err(JetErr::Punt);
    }
    let w = 1usize << bloq; // bits per element
    let blocks = bits.div_ceil(w); // (rip ...) length
    let n = blocks - 1; // after snip
    let bytes = data.to_le_bytes();
    // Element i occupies bits [i*w, (i+1)*w); w may be below 8 for bloq < 3.
    let elem = |i: usize| -> Vec<u8> {
        if w >= 8 {
            let wb = w / 8;
            let mut v = vec![0u8; wb];
            let start = i * wb;
            for (k, slot) in v.iter_mut().enumerate() {
                *slot = *bytes.get(start + k).unwrap_or(&0);
            }
            v
        } else {
            let bit = i * w;
            let byte = *bytes.get(bit / 8).unwrap_or(&0);
            vec![(byte >> (bit % 8)) & ((1u8 << w) - 1)]
        }
    };
    let mut list = D(0);
    for i in (0..n).rev() {
        let v = elem(i);
        let atom = Atom::from_bytes(&mut context.stack, &v).as_noun();
        list = T(&mut context.stack, &[atom, list]);
    }
    Ok(list)
}

// ------------------------------------------------------------------- range

fn range_w<F: BlasFloat>(a: u128, b: u128, d: u128, r: Round) -> Option<Vec<u128>> {
    let (a, b, d): (F, F, F) = (from_raw(a), from_raw(b), from_raw(d));
    // ba = b - a; descending iff ba < 0 (NaN is "not less", so ascending).
    let desc = b.sub(a, r).lt(F::ZERO);
    let mut out = vec![a];
    let mut cur = a;
    loop {
        let next = cur.add(d, r);
        let stop = if desc { next.le(b) } else { b.le(next) };
        if stop {
            return Some(out.into_iter().map(to_raw).collect());
        }
        if out.len() >= RANGE_CAP {
            return None;
        }
        out.push(next);
        cur = next;
    }
}

/// `+range`: `~/  %range`, `[=meta [a=@ b=@] d=@] -> ray`. Starts at `a`,
/// repeatedly adds `d` (rounded in the door's mode), stopping before the
/// first value that reaches `b` (`>=` ascending, `<=` descending, decided
/// by the sign of `b - a`). Result meta `[~[len] bloq %i754 tail]`.
pub fn range(context: &mut Context, subject: Noun) -> Result {
    ray::trace("range");
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let (_shape, bloq, tail) = parse_meta(slot(sam, 2, &space)?, &space)?;
    let a = atom_u128(slot(sam, 12, &space)?, &space)?;
    let b = atom_u128(slot(sam, 13, &space)?, &space)?;
    let d = atom_u128(slot(sam, 7, &space)?, &space)?;
    let rnd = ray::rounding(subject, &space)?;
    let out = by_bloq!(bloq, range_w(a, b, d, rnd)).ok_or(JetErr::Punt)?;
    let data = ray::pack(context, 1usize << (bloq - 3), &out);
    let meta = ray::build_meta(context, &[out.len() as u64], bloq, KIND_I754, tail);
    Ok(ray::build(context, meta, data))
}

// ---------------------------------------------------------------- linspace

fn linspace_w<F: BlasFloat>(a: u128, b: u128, n: u64, r: Round) -> Vec<u128> {
    let (a, b): (F, F) = (from_raw(a), from_raw(b));
    if n == 1 {
        return vec![to_raw(a)];
    }
    // d = (b - a) / (n - 1); element i = a + i * d for i < n - 1; last is b.
    let d = b.sub(a, r).div(F::from_u64(n - 1, r), r);
    let mut out = Vec::with_capacity(n as usize);
    for i in 0..n - 1 {
        out.push(to_raw(a.add(F::from_u64(i, r).mul(d, r), r)));
    }
    out.push(to_raw(b));
    out
}

/// `+linspace`: `~/  %linspace`, `[=meta [a=@ b=@] n=@ud] -> ray`. `n`
/// evenly spaced values from `a` to `b` inclusive; `n = 0` crashes in the
/// Hoon (punt), `n = 1` gives `[a]`. Result meta `[~[n] bloq %i754 tail]`.
pub fn linspace(context: &mut Context, subject: Noun) -> Result {
    ray::trace("linspace");
    let space = context.stack.noun_space();
    let sam = ray::sample(subject, &space)?;
    let (_shape, bloq, tail) = parse_meta(slot(sam, 2, &space)?, &space)?;
    let a = atom_u128(slot(sam, 12, &space)?, &space)?;
    let b = atom_u128(slot(sam, 13, &space)?, &space)?;
    let n = slot(sam, 7, &space)?.in_space(&space).as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
    if n == 0 || n > RANGE_CAP as u64 {
        return Err(JetErr::Punt);
    }
    let rnd = ray::rounding(subject, &space)?;
    let out = by_bloq!(bloq, linspace_w(a, b, n, rnd));
    let data = ray::pack(context, 1usize << (bloq - 3), &out);
    let meta = ray::build_meta(context, &[n], bloq, KIND_I754, tail);
    Ok(ray::build(context, meta, data))
}
