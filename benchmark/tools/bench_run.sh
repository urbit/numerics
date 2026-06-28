#!/usr/bin/env bash
#
# bench_run.sh — drive the math timing grid on one fresh 4.4 ship and emit CSV.
#
# One jetted binary per word-size; the jet is toggled in Hoon by swapping
# math.hoon (NOT separate nojet binaries):
#   jetted     : math.hoon hints intact            -> +bench-grid %cheb ...
#   cheb-interp: comment `~% %non` AND `~/ %math`   -> +bench-grid %cheb ...
#   taylor     : bare math-taylor lib (no hints)    -> +bench-grid %taylor ...
#
# The WHOLE grid for one mode is run by a single dojo command (+bench-grid), which
# slogs `[%cell door arm set]` then the cell's `took ..` for all doors x arms x
# sets, ending with `[%grid-done]`.  Per-command dojo input is unreliable, so we
# issue only ~3 commands total and scrape the streamed slogs.
#
# Inputs are precomputed inside the gen OUTSIDE ~>(%bout); we scrape only `took`,
# so precompute is excluded.  %bout prints "took <label>/<dotted>"; the dotted
# digits ARE microseconds (e.g. s/93.308.561 = 93308561 us) — we strip the dots.
#
# Usage:
#   tools/bench_run.sh 32 ./bin/urbit-jet-32 results/$(date +%F)/bench-timing.csv
#   BENCH_SMOKE=1 tools/bench_run.sh 32 ./bin/urbit-jet-32 results/smoke/bench-timing.csv
#
set -uo pipefail

WS="${1:?word size 32|64}"
BIN="${2:?path to jet binary}"
OUT="${3:?output csv path}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
SESSION="bench${WS}"
PIER="piers/bench${WS}"
DESK="desk"
RAW="$(dirname "$OUT")/raw-${WS}"

N_JET=100000     # jetted reps (all doors)
N_INT=10000      # cheb-interp reps for @rd/@rs (base ops jetted)
N_RH=100         # cheb-interp reps for @rh (base ops NOT jetted -> ~1000x slower)
N_RQ=1000        # cheb-interp reps for @rq (f128 base ops likely un-jetted too)
# Taylor (iterative Sum x^i/i!) converges in thousands of terms for asin/acos/log
# -> ~12 min/cell at 10k.  Tiny counts; higher variance but enough to show the gap.
N_TRD=10; N_TRS=10; N_TRH=10; N_TRQ=10
SETS=5
GRID_TIMEOUT=5400    # 90-min cap per +bench-grid (a slow door bails instead of hanging)

# SKIP_JETTED=1 reuses an already-banked jetted run (interpreted-only re-run).
if [ "${BENCH_SMOKE:-0}" = "1" ]; then
  N_JET=1000; N_INT=1000; N_RH=50; N_RQ=200; N_TRD=30; N_TRS=30; N_TRH=15; N_TRQ=15; SETS=2; GRID_TIMEOUT=1800
fi

mkdir -p "$(dirname "$OUT")"
echo "word_size,impl,door,arm,n,set,elapsed_us,status" > "$OUT"

clean_prompt() {
  local i
  for i in 1 2 3; do
    tmux send-keys -t "$SESSION" C-u;   sleep 0.4
    tmux send-keys -t "$SESSION" Enter; sleep 0.8
  done
  sleep 1
}
# IMPORTANT: combined form (text + Enter in one call).  The -l flag + a separate
# Enter does NOT submit `+gen` commands in dojo; this form does.
send_cmd() { tmux send-keys -t "$SESSION" "$1" Enter; }

# dojo only reliably submits @ud args in DOTTED form (e.g. 100.000, not 100000) —
# plain digits leave the +gen line unsubmitted.  Group digits in threes.
dots() { perl -e '$_=shift; s/(?<=\d)(?=(\d{3})+$)/./g; print' "$1"; }

boot_ship() {
  # NB: macOS pgrep has no \b — match the plain pier path.
  if tmux has-session -t "$SESSION" 2>/dev/null \
     && pgrep -f "piers/bench${WS}" >/dev/null \
     && tmux capture-pane -t "$SESSION" -p 2>/dev/null | grep -q 'dojo>'; then
    echo "[ship up — reusing]"; return
  fi
  pkill -TERM -f "piers/bench${WS}" 2>/dev/null
  pkill -TERM -f "$(basename "$BIN") work" 2>/dev/null; sleep 2
  /bin/rm -rf "$PIER"
  tmux kill-session -t "$SESSION" 2>/dev/null
  tmux new-session -d -s "$SESSION" -x 220 -y 50
  tmux set-option -t "$SESSION" history-limit 200000 2>/dev/null
  tmux send-keys -t "$SESSION" "cd $ROOT && $BIN -F zod -c $PIER" Enter
  echo "[booting $SESSION — downloads urbit-v4.4.pill]"
  local w=0
  until tmux capture-pane -t "$SESSION" -p | grep -q 'dojo>'; do
    sleep 5; w=$((w+5)); [ $w -gt 360 ] && { echo "boot timeout"; exit 1; }
  done
  sleep 5
}

sync_desk() {   # $1 = math.hoon variant to install
  clean_prompt
  if ! ls "$PIER/base/gen" >/dev/null 2>&1; then
    send_cmd '|mount %base'; sleep 6
    local w=0
    until ls "$PIER/base/gen" >/dev/null 2>&1; do sleep 3; w=$((w+3)); [ $w -gt 60 ] && break; done
  fi
  /bin/cp -f "$DESK"/lib/bench-core.hoon "$DESK"/lib/bench-domains.hoon \
             "$DESK"/lib/bench-cells.hoon "$DESK"/lib/math-taylor.hoon "$PIER"/base/lib/
  /bin/cp -f "$1" "$PIER"/base/lib/math.hoon
  /bin/cp -f "$DESK"/gen/bench-math.hoon "$DESK"/gen/bench-grid.hoon "$PIER"/base/gen/
  clean_prompt
  send_cmd '|commit %base'
  local w=0
  until tmux capture-pane -t "$SESSION" -p -S -25 | grep -qE 'bench-grid|math/hoon'; do
    sleep 2; w=$((w+2)); [ $w -gt 60 ] && break; done
  sleep 2
}

# run one whole grid -> append parsed rows to $OUT
run_grid() {
  local tag="$1" gimpl="$2" sets="$3" nrd="$4" nrs="$5" nrh="$6" nrq="$7" scope="${8:-all}"
  echo "[grid] impl=$tag gen=$gimpl sets=$sets scope=$scope n=(rd:$nrd rs:$nrs rh:$nrh rq:$nrq)"
  clean_prompt
  # count existing done/fail markers; wait for a NEW one (stale markers from prior
  # runs must not short-circuit the wait or trip a false build-fail).
  local g0 f0
  g0=$(tmux capture-pane -t "$SESSION" -p -S -200000 | grep -c 'grid-done')
  f0=$(tmux capture-pane -t "$SESSION" -p -S -200000 | grep -c 'generator-build-fail')
  send_cmd "+bench-grid %${gimpl} $(dots "$sets") $(dots "$nrd") $(dots "$nrs") $(dots "$nrh") $(dots "$nrq") %${scope}"
  local w=0
  until [ "$(tmux capture-pane -t "$SESSION" -p -S -200000 | grep -c 'grid-done')" -gt "$g0" ]; do
    if [ "$(tmux capture-pane -t "$SESSION" -p -S -200000 | grep -c 'generator-build-fail')" -gt "$f0" ]; then
      echo "  [grid BUILD FAILED — see session $SESSION]"; return 1
    fi
    sleep 10; w=$((w+10))
    [ $((w % 60)) -eq 0 ] && echo "  ${w}s, cells so far: $(tmux capture-pane -t "$SESSION" -p -S -200000 | grep -c '%cell')"
    [ $w -gt "$GRID_TIMEOUT" ] && { echo "  [grid TIMEOUT]"; break; }
  done
  tmux capture-pane -t "$SESSION" -p -S -200000 > "${RAW}-${tag}.txt"
  # NB: parse via the standalone script (reads the raw as an ARG).  An inline
  # `python3 - … < raw <<'PY'` collides stdin (heredoc vs `< raw`) -> 0 rows.
  local before; before=$(wc -l < "$OUT")
  python3 "$ROOT/tools/parse_raw.py" "$tag" "$WS" "${RAW}-${tag}.txt" >> "$OUT"
  echo "  parsed $(( $(wc -l < "$OUT") - before )) rows"
}

# --- math.hoon variants ------------------------------------------------------
HINTED="$DESK/lib/math.hoon"
NOJET="/tmp/math-nojet-${WS}.hoon"
sed -e 's|^~%  %non  ..part  ~|:: ~%  %non  ..part  ~  ::NOJET|' \
    -e 's|^  ~/  %math|  :: ~/  %math  ::NOJET|' \
    "$HINTED" > "$NOJET"

# --- run --------------------------------------------------------------------
boot_ship

if [ "${SKIP_JETTED:-0}" != "1" ]; then
  echo "[mode: jetted — hinted math.hoon]"
  sync_desk "$HINTED"
  run_grid jetted cheb "$SETS" "$N_JET" "$N_JET" "$N_JET" "$N_JET" all
fi

echo "[mode: taylor — hinted math.hoon, tiny n (iterative is slow)]"
sync_desk "$HINTED"
run_grid taylor taylor "$SETS" "$N_TRD" "$N_TRS" "$N_TRH" "$N_TRQ" two

# cheb-interp must run on the NOJET BINARY (math.hoon hints intact, no jet
# registered).  Commenting ~% on the jet binary floods `fund: parent not found`.
# Run cheb separately on the nojet binary; SKIP_CHEB=1 omits it here.
if [ "${SKIP_CHEB:-0}" != "1" ]; then
  echo "[mode: cheb-interp — commented math.hoon (NB: prefer nojet binary)]"
  sync_desk "$NOJET"
  run_grid cheb cheb "$SETS" "$N_INT" "$N_INT" "$N_RH" "$N_RQ" two
fi

echo "[done] -> $OUT  ($(($(wc -l < "$OUT")-1)) rows)"
