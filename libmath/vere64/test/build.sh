#!/bin/sh
# Build + run the @rd bit-exact harness.  Uses the SAME Berkeley SoftFloat the
# vere jets use.  The zig-archived .a is not 8-byte aligned for Apple ld, so we
# re-archive it with libtool (chmod the extracted .o's; they come out read-only).
set -e

SFINC="$(dirname "$(find "$HOME/urbit" -name softfloat.h 2>/dev/null | grep -i softfloat/source/include | head -1)")"
SFLIB="$(find "$HOME/urbit/vere-ml64" -name libsoftfloat.a 2>/dev/null | head -1)"
[ -n "$SFINC" ] && [ -n "$SFLIB" ] || { echo "softfloat not found (build vere-ml64 first)"; exit 1; }

WORK="$(mktemp -d)"
( cd "$WORK" && ar x "$SFLIB" && chmod u+rw ./*.o && libtool -static -o ./libsoftfloat.a ./*.o )
cc -O2 -I"$SFINC" "$(dirname "$0")/rd_check.c" "$WORK/libsoftfloat.a" -o "$WORK/rd_check"
"$WORK/rd_check"
rm -rf "$WORK"

cat <<'EOF'

# To compare against the Hoon (on a ship with /+ math), the dojo expression is:
#   =rd  rd:math
#   :: then for each input X above, check `@ux`(exp:rd:math X) / (log:rd:math X)
# matches the harness <out> bit-for-bit.
EOF
