#!/bin/bash
# Round-trip test for compact_db: builds a synthetic DB, compacts it and
# checks that compact_db's verification passes and that classify gives
# identical results (also in chunked mode) for both formats.
# Usage: tests/test_compact_db.sh   (after 'make' in src/)

set -eu
SRC=$(cd "$(dirname "$0")/../src" && pwd)
TESTS=$(cd "$(dirname "$0")" && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

run_case() {
  local K=$1 NT=$2
  local D=$TMP/k$K-m$NT
  echo "=== k=$K minimizer=$NT"
  python3 "$TESTS/gen_compact_test_data.py" "$D" $K
  "$SRC/db_sort" -M -d "$D/unsorted.kdb" -o "$D/db.kdb" -i "$D/db.idx" -n $NT 2>/dev/null
  "$SRC/build_taxdb" "$D/names.dmp" "$D/nodes.dmp" > "$D/taxDB"
  "$SRC/compact_db" -v -d "$D/db.kdb" -i "$D/db.idx" -o "$D/db.kdbc" -x "$D/db.idxc"

  classify() { # name db idx extra-args...
    local name=$1 db=$2 idx=$3; shift 3
    "$SRC/classify" -d "$db" -i "$idx" -a "$D/taxDB" -o "$D/$name.out" -r "$D/$name.report" "$@" "$D/reads.fa" 2>/dev/null
  }
  classify orig "$D/db.kdb" "$D/db.idx"
  classify compact "$D/db.kdbc" "$D/db.idxc"
  classify compact_mem "$D/db.kdbc" "$D/db.idxc" -M
  cmp "$D/orig.out" "$D/compact.out"
  cmp "$D/orig.out" "$D/compact_mem.out"
  cmp "$D/orig.report" "$D/compact.report"
  # chunked mode, small chunks (forces several chunks)
  local CH=$(( $(stat -c %s "$D/db.kdbc") / 5 ))
  classify orig_chunk "$D/db.kdb" "$D/db.idx" -x $(( $(stat -c %s "$D/db.kdb") / 5 + $(stat -c %s "$D/db.idx") ))
  classify compact_chunk "$D/db.kdbc" "$D/db.idxc" -x $CH
  cmp "$D/orig.out" "$D/compact_chunk.out"
  cmp "$D/orig.out" "$D/orig_chunk.out"
  echo "classified: $(grep -c '^C' "$D/orig.out") of $(wc -l < "$D/orig.out") reads; outputs identical"
}

# taxids needing more than 26 bits must be rejected
echo "=== oversized taxid"
D=$TMP/big
python3 "$TESTS/gen_compact_test_data.py" "$D" 21 > /dev/null
python3 - "$D/unsorted.kdb" <<'PYEOF'
import struct, sys
with open(sys.argv[1], 'r+b') as f:
    f.seek(72 + 2 * (4 + 8 * 42) + 6)   # first record: 6 key bytes, then taxid
    f.write(struct.pack('<I', 1 << 26))
PYEOF
"$SRC/db_sort" -M -d "$D/unsorted.kdb" -o "$D/db.kdb" -i "$D/db.idx" -n 7 2>/dev/null
if "$SRC/compact_db" -d "$D/db.kdb" -i "$D/db.idx" -o "$D/db.kdbc" -x "$D/db.idxc" 2> "$D/err"; then
  echo "compact_db accepted a 27-bit taxid"; exit 1
fi
grep -q "needs more than 26 bits" "$D/err"
echo "oversized taxid rejected"

run_case 31 10
run_case 31 12
run_case 25 9
run_case 21 7
echo "ALL OK"
