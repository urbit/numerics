#!/bin/sh
# verify.sh <fn> [fn2 ...] -- diff harness <fn> outputs against the live Hoon
# (<fn>:rd:math) in the tmux 'dev' ship, bit-for-bit.  Run ./build.sh first to
# refresh /tmp/rd_harness.txt; this script reuses it.
#
#   ./build.sh > /tmp/rd_harness.txt 2>/dev/null
#   ./verify.sh exp log sin cos
norm() { printf '%s' "$1" | tr -d '0xX.' | sed 's/^0*//'; }   # canonical hex (no 0x/dots/leading-0)
# minimal dot-grouped @rd literal: strip leading 0s, group every 4 from the right
# (so 0x0000000000000000 -> 0x0, not a wrapping 0x0000.0000.0000.0000)
dotg() {
  h=$(printf '%s' "$1" | sed 's/^0x//;s/^0*//'); [ -z "$h" ] && h=0
  printf '0x%s' "$(printf '%s' "$h" | rev | sed 's/.\{4\}/&./g;s/\.$//' | rev)"
}
# warmup: wait until the dojo answers a known query (clears stale pane state)
tmux clear-history -t dev 2>/dev/null
i=0
while [ $i -lt 20 ]; do
  tmux send-keys -t dev C-c 2>/dev/null; sleep 0.3; tmux send-keys -t dev C-u 2>/dev/null; sleep 0.3
  tmux send-keys -t dev -l '`@ux`(exp:rd:math `@rd`0x3ff0.0000.0000.0000)' 2>/dev/null; sleep 0.2
  tmux send-keys -t dev Enter 2>/dev/null; sleep 1.5
  tmux capture-pane -t dev -p 2>/dev/null | grep -qE '^0x4005.bf0a' && break
  i=$((i+1))
done

for FN in "$@"; do
  ok=0; bad=0
  grep "^$FN " /tmp/rd_harness.txt | while read -r f in out _; do
    din=$(dotg "$in")
    tmux send-keys -t dev C-c 2>/dev/null; sleep 0.4
    tmux send-keys -t dev C-u 2>/dev/null; sleep 0.4
    tmux send-keys -t dev -l "\`@ux\`($FN:rd:math \`@rd\`$din)" 2>/dev/null; sleep 0.3
    tmux send-keys -t dev Enter 2>/dev/null; sleep 1.6
    got=$(tmux capture-pane -t dev -p 2>/dev/null | grep -E '^0x[0-9a-f]' | tail -1)
    if [ "$(norm "$got")" = "$(norm "$out")" ]; then
      ok=$((ok+1))
    else
      bad=$((bad+1)); printf 'DIFF %s %s  hoon=%s  C=%s\n' "$FN" "$in" "$got" "$out"
    fi
    printf '%s\n' "$ok $bad" > /tmp/rd_tally.$FN
  done
  read o b < /tmp/rd_tally.$FN 2>/dev/null
  printf '== %s: %s OK, %s DIFF ==\n' "$FN" "${o:-0}" "${b:-0}"
done
