#!/bin/sh
# verify_clean.sh <session> <fn> [fn2 ...]
# One-at-a-time, NO C-c (SIGINT corrupts the loom mid-event).  Clears the dojo
# line with C-u only, generous spacing, checks ship liveness after every query.
# Compares harness /tmp/rd_harness.txt outputs against live (<fn>:rd:math) bit-
# for-bit.  Aborts immediately and names the killer input if the serf dies.
SES="$1"; shift
norm() { printf '%s' "$1" | tr -d '0xX.' | sed 's/^0*//'; }
dotg() {
  h=$(printf '%s' "$1" | sed 's/^0x//;s/^0*//'); [ -z "$h" ] && h=0
  printf '0x%s' "$(printf '%s' "$h" | rev | sed 's/.\{4\}/&./g;s/\.$//' | rev)"
}
alive() { pgrep -fq "urbit -F $SES" && return 0 || return 1; }
ask() {  # $1 = dojo expression; echoes first 0x line of the response
  tmux send-keys -t "$SES" C-u 2>/dev/null; sleep 0.3
  tmux send-keys -t "$SES" -l "$1" 2>/dev/null; sleep 0.3
  tmux send-keys -t "$SES" Enter 2>/dev/null; sleep 1.6
  tmux capture-pane -t "$SES" -p 2>/dev/null | grep -E '^0x[0-9a-f]' | tail -1
}

for FN in "$@"; do
  ok=0; bad=0; dead=0
  while read -r f in out _; do
    [ "$f" = "$FN" ] || continue
    alive || { echo "ABORT: ship dead before $FN $in"; dead=1; break; }
    din=$(dotg "$in")
    got=$(ask "\`@ux\`($FN:rd:math \`@rd\`$din)")
    if ! alive; then echo "KILLER: $FN $din  (ship died on this input)"; dead=1; break; fi
    if [ "$(norm "$got")" = "$(norm "$out")" ]; then ok=$((ok+1))
    else bad=$((bad+1)); printf 'DIFF %s %s  hoon=%s  C=%s\n' "$FN" "$din" "$got" "$out"; fi
  done < /tmp/rd_harness.txt
  printf '== %s: %s OK, %s DIFF%s ==\n' "$FN" "$ok" "$bad" "$([ $dead = 1 ] && echo ' (ABORTED)')"
  [ $dead = 1 ] && break
done
