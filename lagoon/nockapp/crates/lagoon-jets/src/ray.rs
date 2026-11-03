//! Reading and writing `$ray` nouns.
//!
//! `/sur/lagoon`: `ray = [meta data=@ux]`, `meta = [shape=(list @) bloq kind tail]`.
//! `data` is one atom, row-major, element `i` at bit offset `i * 2^bloq` from
//! the LSB, with a single 1 pin bit at offset `len * 2^bloq` (the Hoon
//! `+check`: `(met bloq data) == (roll shape mul) + 1`).

use nockvm::interpreter::Context;
use nockvm::jets::util::slot;
use nockvm::jets::{JetErr, Result};
use nockvm::noun::{Atom, Noun, NounSpace, D, T};
use nockvm::ext::AtomExt;
use sdfloat::Round;

/// Little-endian `@tas` value of a short ASCII symbol.
pub const fn tas(s: &[u8]) -> u64 {
    let mut v = 0u64;
    let mut i = s.len();
    while i > 0 {
        i -= 1;
        v = (v << 8) | s[i] as u64;
    }
    v
}

pub const KIND_I754: u64 = tas(b"i754");

/// A parsed ray. `meta` is the original meta noun, reused verbatim in results.
#[derive(Clone, Copy)]
pub struct Ray {
    pub meta: Noun,
    pub shape_noun: Noun,
    pub bloq: u32,
    pub kind: u64,
    pub tail: Noun,
    pub data: Atom,
    /// Product of the shape.
    pub len: u64,
}

impl Ray {
    /// Element width in bytes.
    pub fn width(&self) -> usize {
        1usize << (self.bloq - 3)
    }
}

/// Walk a Hoon list of atoms.
pub fn list_atoms(mut n: Noun, space: &NounSpace) -> std::result::Result<Vec<u64>, JetErr> {
    let mut out = Vec::new();
    loop {
        let h = n.in_space(space);
        if let Ok(a) = h.as_atom() {
            if a.as_u64().map_err(|_| JetErr::Punt)? == 0 {
                return Ok(out);
            }
            return Err(JetErr::Punt);
        }
        let c = h.as_cell().map_err(|_| JetErr::Punt)?;
        let v = c.head().as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
        out.push(v);
        n = c.tail().noun();
    }
}

/// Parse `[meta data]`, validating as `+check` does. Punts on anything odd.
pub fn parse(ray: Noun, space: &NounSpace) -> std::result::Result<Ray, JetErr> {
    let c = ray.in_space(space).as_cell().map_err(|_| JetErr::Punt)?;
    let meta = c.head().noun();
    let data = c.tail().as_atom().map_err(|_| JetErr::Punt)?;
    let m = meta.in_space(space).as_cell().map_err(|_| JetErr::Punt)?;
    let shape_noun = m.head().noun();
    let m2 = m.tail().as_cell().map_err(|_| JetErr::Punt)?;
    let bloq = m2.head().as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
    let m3 = m2.tail().as_cell().map_err(|_| JetErr::Punt)?;
    let kind = m3.head().as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
    let tail = m3.tail().noun();
    if !(4..=7).contains(&bloq) {
        return Err(JetErr::Punt);
    }
    let dims = list_atoms(shape_noun, space)?;
    let mut len: u64 = 1;
    for d in &dims {
        len = len.checked_mul(*d).ok_or(JetErr::Punt)?;
    }
    // +check: (met bloq data) == len + 1, i.e. bit_size in (w*len, w*(len+1)].
    let w = 1u64 << bloq;
    let bits = data.bit_size() as u64;
    let lo = w.checked_mul(len).ok_or(JetErr::Punt)?;
    let hi = lo.checked_add(w).ok_or(JetErr::Punt)?;
    if !(bits > lo && bits <= hi) {
        return Err(JetErr::Punt);
    }
    Ok(Ray { meta, shape_noun, bloq: bloq as u32, kind, tail, data: data.atom(), len })
}

/// The raw element bytes (little-endian, pin excluded), exactly `len * width`.
pub fn bytes(r: &Ray, space: &NounSpace) -> Vec<u8> {
    let n = (r.len as usize) * r.width();
    let mut v = r.data.in_space(space).to_le_bytes();
    v.resize(n.max(v.len()), 0);
    v.truncate(n);
    v
}

/// Elements as raw bit patterns.
pub fn elems(r: &Ray, space: &NounSpace) -> Vec<u128> {
    let w = r.width();
    bytes(r, space)
        .chunks(w)
        .map(|c| {
            let mut b = [0u8; 16];
            b[..w].copy_from_slice(c);
            u128::from_le_bytes(b)
        })
        .collect()
}

/// Pack elements (raw bit patterns) into a pinned data atom.
pub fn pack(context: &mut Context, width: usize, elems: &[u128]) -> Atom {
    let mut v = Vec::with_capacity(elems.len() * width + 1);
    for e in elems {
        v.extend_from_slice(&e.to_le_bytes()[..width]);
    }
    v.push(1);
    Atom::from_bytes(&mut context.stack, &v)
}

/// Build `[meta data]` reusing `meta`.
pub fn build(context: &mut Context, meta: Noun, data: Atom) -> Noun {
    T(&mut context.stack, &[meta, data.as_noun()])
}

/// Build a meta `[shape bloq kind tail]` from a fresh shape.
pub fn build_meta(context: &mut Context, shape: &[u64], bloq: u32, kind: u64, tail: Noun) -> Noun {
    let mut list = D(0);
    for d in shape.iter().rev() {
        let a = Atom::new(&mut context.stack, *d).as_noun();
        list = T(&mut context.stack, &[a, list]);
    }
    T(&mut context.stack, &[list, D(bloq as u64), D(kind), tail])
}

/// The door's rounding mode: `rnd` sits at axis 30 of a gate core inside `+la`.
pub fn rounding(subject: Noun, space: &NounSpace) -> std::result::Result<Round, JetErr> {
    let r = slot(subject, 30, space).map_err(|_| JetErr::Punt)?;
    let v = r.in_space(space).as_atom().map_err(|_| JetErr::Punt)?.as_u64().map_err(|_| JetErr::Punt)?;
    Ok(match v {
        x if x == tas(b"n") => Round::NearEven,
        x if x == tas(b"u") => Round::Up,
        x if x == tas(b"d") => Round::Down,
        x if x == tas(b"z") => Round::Zero,
        x if x == tas(b"a") => Round::Away,
        _ => return Err(JetErr::Punt),
    })
}

/// The gate sample.
pub fn sample(subject: Noun, space: &NounSpace) -> Result {
    slot(subject, 6, space)
}

/// Print `msg` to stderr when `LAGOON_JET_TRACE` is set (for checking that a
/// jet actually fires; test mode alone is silent when a jet never matches).
pub fn trace(msg: &str) {
    if std::env::var_os("LAGOON_JET_TRACE").is_some() {
        eprintln!("lagoon-jet: {msg}");
    }
}
