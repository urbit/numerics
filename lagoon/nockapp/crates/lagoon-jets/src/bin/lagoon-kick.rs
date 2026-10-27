//! lagoon-kick: fire a hoonc `--arbitrary` trap under the lagoon jets.
//!
//! The build trap's value must itself be a trap (`hoon/run-tests.hoon`), so
//! that the tests run here, under these jets, rather than inside hoonc. That
//! inner trap is kicked too and must produce `(list [name=@t lines=(list @t)])`:
//! one row per test, an empty `lines` meaning pass. Set `NOCK_TEST_JETS` to a
//! comma-separated list of jet paths (`k.138/one/two/tri/qua/pen/non/lagoon/mmul`)
//! to have the interpreter run each such jet against the raw Nock and bail
//! on any mismatch.
//!
//!     lagoon-kick <trap.jam>

use std::{env, fs, process};

use lagoon_jets::hot_state;
use nockapp::kernel::boot::parse_test_jets;
use nockapp::utils::create_context;
use nockvm::ext::NounExt;
use nockvm::interpreter::{interpret, Context};
use nockvm::jets::cold::Cold;
use nockvm::jets::JetDispatchMode;
use nockvm::mem::{NockStack, NOCK_STACK_SIZE_MEDIUM};
use nockvm::noun::{Noun, D, T};

fn cord(n: Noun, space: &nockvm::noun::NounSpace) -> String {
    match n.in_space(space).as_atom() {
        Ok(a) => a.into_string().unwrap_or_else(|_| "<non-utf8>".into()),
        Err(_) => "<cell>".into(),
    }
}

fn list(mut n: Noun, space: &nockvm::noun::NounSpace) -> Vec<Noun> {
    let mut out = Vec::new();
    while let Ok(c) = n.in_space(space).as_cell() {
        out.push(c.head().noun());
        n = c.tail().noun();
    }
    out
}

fn main() {
    let path = match env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: lagoon-kick <trap.jam>");
            process::exit(2);
        }
    };
    let jam = fs::read(&path).unwrap_or_else(|e| {
        eprintln!("cannot read {path}: {e}");
        process::exit(1);
    });
    let test_jets = parse_test_jets(&env::var("NOCK_TEST_JETS").unwrap_or_default());
    let mut stack = NockStack::new(NOCK_STACK_SIZE_MEDIUM, 0);
    let cold = Cold::new(&mut stack);
    let hot = hot_state();
    let mut context = create_context(stack, &hot, cold, None, test_jets, JetDispatchMode::Exact);

    let trap = <Noun as NounExt>::cue_bytes_slice(&mut context.stack, &jam).unwrap_or_else(|e| {
        eprintln!("{path} is not a valid jam: {e:?}");
        process::exit(1);
    });
    let raw = env::args().any(|a| a == "--raw");
    let t0 = std::time::Instant::now();
    let built = unsafe {
        context.with_stack_frame(0, |context: &mut Context| {
            let kick = T(&mut context.stack, &[D(9), D(2), D(0), D(1)]);
            interpret(context, trap, kick)
        })
    };
    eprintln!("build kick: {} ms", t0.elapsed().as_millis());
    let built = match built {
        Ok(p) => p,
        Err(e) => {
            eprintln!("build trap crashed: {e:?}");
            process::exit(1);
        }
    };
    if raw {
        let space = context.stack.noun_space();
        match built.in_space(&space).as_atom() {
            Ok(a) => println!("atom: {:?}", a.as_u64().ok()),
            Err(_) => println!("cell"),
        }
        return;
    }
    let t0 = std::time::Instant::now();
    let product = unsafe {
        context.with_stack_frame(0, |context: &mut Context| {
            let kick = T(&mut context.stack, &[D(9), D(2), D(0), D(1)]);
            interpret(context, built, kick)
        })
    };
    eprintln!("suite kick: {} ms", t0.elapsed().as_millis());
    let product = match product {
        Ok(p) => p,
        Err(e) => {
            eprintln!("trap crashed: {e:?}");
            process::exit(1);
        }
    };
    let filter: Vec<String> = env::var("LAGOON_TESTS")
        .map(|s| s.split(',').filter(|p| !p.is_empty()).map(String::from).collect())
        .unwrap_or_default();
    let rows = {
        let space = context.stack.noun_space();
        list(product, &space)
            .into_iter()
            .filter_map(|row| {
                let c = row.in_space(&space).as_cell().ok()?;
                Some((cord(c.head().noun(), &space), c.tail().noun()))
            })
            .filter(|(name, _)| filter.is_empty() || filter.iter().any(|p| name.contains(p.as_str())))
            .collect::<Vec<_>>()
    };
    let mut failed = 0usize;
    for (name, trap) in &rows {
        let t0 = std::time::Instant::now();
        let res = unsafe {
            context.with_stack_frame(0, |context: &mut Context| {
                let kick = T(&mut context.stack, &[D(9), D(2), D(0), D(1)]);
                interpret(context, *trap, kick)
            })
        };
        let ms = t0.elapsed().as_millis();
        let space = context.stack.noun_space();
        match res {
            Ok(lines_noun) => {
                let lines = list(lines_noun, &space);
                if lines.is_empty() {
                    println!("ok    {name}  ({ms} ms)");
                } else {
                    failed += 1;
                    println!("FAIL  {name}  ({ms} ms)");
                    for l in lines {
                        println!("      {}", cord(l, &space));
                    }
                }
            }
            Err(e) => {
                failed += 1;
                println!("CRASH {name}  ({ms} ms): {e:?}");
            }
        }
    }
    println!("{} tests, {} failed", rows.len(), failed);
    if failed > 0 {
        process::exit(1);
    }
}
