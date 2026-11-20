#!/bin/bash
# Round-trip test for compact_db: builds a synthetic DB, compacts it and
# checks that compact_db's verification passes and that classify gives
# identical results (also in chunked mode) for both formats.
# Usage: tests/test_compact_db.sh   (after 'make' in src/)
# KRAKEN_TEST_WIDE_INDEX=1 also tests a legacy 17-nt index, using a sparse
# 128 GiB file and mapping that much virtual address space.

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
  if [ "$K" = 31 ] && [ "$NT" = 10 ]; then
    # A second database has the same keys but conflicting taxids, making
    # first-database precedence observable.
    cp "$D/db.kdb" "$D/conflict.kdb"
    python3 - "$D/conflict.kdb" <<'PYEOF'
import struct, sys
with open(sys.argv[1], 'r+b') as f:
    data = bytearray(f.read())
    header = 72 + 2 * (4 + 8 * 62)
    for pos in range(header + 8, len(data), 12):
        struct.pack_into('<I', data, pos, 100)
    f.seek(0)
    f.write(data)
PYEOF
    cp "$D/reads.fa" "$D/batch.fa"
    printf '>short\nACGT\n>ambiguous\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n' >> "$D/batch.fa"
    awk 'NR == 2 {for (i = 0; i < 30; i++) print ">duplicate" i "\n" $0; exit}' \
      "$D/reads.fa" >> "$D/batch.fa"
    awk 'NR % 2 {id=substr($0,2); next} {print "@" id "\n" $0 "\n+"; gsub(/./,"I"); print}' \
      "$D/batch.fa" > "$D/batch.fq"
    printf '>short\nACGT\n>ambiguous\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n' > "$D/no-pending.fa"
    threads=(1)
    if [ "$(nproc)" -ge 2 ]; then threads+=(2); fi
    for input in "$D/batch.fa" "$D/batch.fq" "$D/no-pending.fa"; do
      for order in plain compact; do
        if [ "$order" = plain ]; then
          db1="$D/db.kdb"; idx1="$D/db.idx"
          db2="$D/db.kdbc"; idx2="$D/db.idxc"
        else
          db1="$D/db.kdbc"; idx1="$D/db.idxc"
          db2="$D/conflict.kdb"; idx2="$D/db.idx"
        fi
        common=(-d "$db1" -i "$idx1" -d "$db2" -i "$idx2"
                -d "$D/conflict.kdb" -i "$D/db.idx" -a "$D/taxDB" -t 1 -u 5000)
        "$SRC/classify" "${common[@]}" -B 0 -o "$D/baseline.out" -r "$D/baseline.report" \
          -C "$D/baseline.classified" -U "$D/baseline.unclassified" "$input" 2>/dev/null
        for t in "${threads[@]}"; do
          for size in 6000 1M; do
            for mode in sorted unsorted; do
              hint=()
              if [ "$mode" = unsorted ]; then hint=(-H); fi
              for strategy in batched outer; do
                outer=()
                if [ "$strategy" = outer ]; then outer=(-D); fi
                KRAKEN_BATCH_FORCE_HINTS=1 "$SRC/classify" "${common[@]}" -t "$t" -B "$size" \
                  "${hint[@]}" "${outer[@]}" -o "$D/batch.out" -r "$D/batch.report" \
                  -C "$D/batch.classified" -U "$D/batch.unclassified" "$input" 2>/dev/null
                cmp "$D/baseline.out" "$D/batch.out"
                cmp "$D/baseline.report" "$D/batch.report"
                cmp "$D/baseline.classified" "$D/batch.classified"
                cmp "$D/baseline.unclassified" "$D/batch.unclassified"
                rm "$D/batch.report"
              done
            done
          done
        done
        rm "$D/baseline.report"
      done
    done
    # First database has no keys, so later databases must resolve queries.
    python3 - "$D" <<'PYEOF'
import struct, sys
from pathlib import Path
d = Path(sys.argv[1])
header_size = 72 + 2 * (4 + 8 * 62)
header = bytearray((d / 'db.kdb').read_bytes()[:header_size])
struct.pack_into('<Q', header, 48, 0)
(d / 'empty.kdb').write_bytes(header)
(d / 'empty.idx').write_bytes(b'KRAKIX2' + bytes([10]) + bytes(8 * ((1 << 20) + 1)))
PYEOF
    gzip -c "$D/batch.fq" > "$D/batch.fq.gz"
    mkdir "$D/tmp"
    common=(-d "$D/empty.kdb" -i "$D/empty.idx" -d "$D/db.kdbc" -i "$D/db.idxc"
            -d "$D/conflict.kdb" -i "$D/db.idx" -a "$D/taxDB" -t 1 -u 5000 -B 6000)
    for executable in classify classifyExact; do
      [ -x "$SRC/$executable" ] || continue
      "$SRC/$executable" "${common[@]}" -o "$D/multi-baseline.out" -r "$D/multi-baseline.report" \
        -C "$D/multi-baseline.classified" -U "$D/multi-baseline.unclassified" \
        "$D/batch.fa" "$D/batch.fq.gz" "$D/no-pending.fa" 2>/dev/null
      TMPDIR="$D/tmp" OMP_DYNAMIC=TRUE KRAKEN_BATCH_FORCE_HINTS=1 \
        "$SRC/$executable" "${common[@]}" -t "${threads[-1]}" -D \
        -o "$D/multi-outer.out" -r "$D/multi-outer.report" -C "$D/multi-outer.classified" \
        -U "$D/multi-outer.unclassified" "$D/batch.fa" "$D/batch.fq.gz" "$D/no-pending.fa" 2>"$D/multi-outer.log"
      [ "$(grep -c '^Database-outer pass ' "$D/multi-outer.log")" -eq 3 ]
      for result in out report classified unclassified; do
        cmp "$D/multi-baseline.$result" "$D/multi-outer.$result"
      done
      rm "$D/multi-baseline.report" "$D/multi-outer.report"
      [ -z "$(ls -A "$D/tmp")" ]
    done
    # A found zero taxid must not fall through to a later conflicting database.
    python3 - "$D/conflict.kdb" "$D/zero.kdb" <<'PYEOF'
import struct, sys
data = bytearray(open(sys.argv[1], 'rb').read())
header = 72 + 2 * (4 + 8 * 62)
for pos in range(header + 8, len(data), 12):
    struct.pack_into('<I', data, pos, 0)
open(sys.argv[2], 'wb').write(data)
PYEOF
    zero=(-d "$D/zero.kdb" -i "$D/db.idx" -d "$D/conflict.kdb" -i "$D/db.idx"
          -a "$D/taxDB" -t 1 -u 5000 -B 6000 -r off)
    "$SRC/classify" "${zero[@]}" -o "$D/zero-batch.out" "$D/batch.fa" 2>/dev/null
    "$SRC/classify" "${zero[@]}" -D -o "$D/zero-outer.out" "$D/batch.fa" 2>/dev/null
    cmp "$D/zero-batch.out" "$D/zero-outer.out"
    "$SRC/classify" "${common[@]}" -D -r off -o "$D/no-pending.out" \
      "$D/no-pending.fa" 2>"$D/no-pending.log"
    if grep -q 'Database-outer pass 2' "$D/no-pending.log"; then
      echo "database-outer mode continued after all queries were resolved"; exit 1
    fi
    for invalid in '-q' '-x 1M' '-x 0' '-B 0'; do
      if "$SRC/classify" "${common[@]}" -D $invalid "$D/batch.fa" 2>"$D/outer-invalid.err"; then
        echo "database-outer mode accepted incompatible $invalid"; exit 1
      fi
      grep -q "database-outer mode requires batching" "$D/outer-invalid.err"
    done
    if TMPDIR="$D/no-such-directory" "$SRC/classify" "${common[@]}" -D "$D/batch.fa" 2>"$D/outer-tmp.err"; then
      echo "database-outer mode accepted an invalid TMPDIR"; exit 1
    fi
    grep -q "creating database-outer temporary file" "$D/outer-tmp.err"
    mkdir "$D/install" "$D/wrapper-db"
    cp "$TESTS/../scripts/krakenuniq" "$TESTS/../scripts/krakenlib.pm" \
      "$TESTS/../scripts/read_merger.pl" "$D/install/"
    ln -s "$SRC/classify" "$D/install/classify"
    ln -s "$D/db.kdbc" "$D/wrapper-db/database.kdb"
    ln -s "$D/db.idxc" "$D/wrapper-db/database.idx"
    ln -s "$D/taxDB" "$D/wrapper-db/taxDB"
    for strategy in batch outer; do
      outer=()
      if [ "$strategy" = outer ]; then outer=(--database-outer); fi
      "$D/install/krakenuniq" --db "$D/wrapper-db" --batch-size 6000 "${outer[@]}" \
        --output "$D/wrapper-$strategy.out" --report-file "$D/wrapper-$strategy.report" \
        "$D/batch.fa" 2>/dev/null
      grep -v '^#' "$D/wrapper-$strategy.report" > "$D/wrapper-$strategy.data"
    done
    cmp "$D/wrapper-batch.out" "$D/wrapper-outer.out"
    cmp "$D/wrapper-batch.data" "$D/wrapper-outer.data"
    "$D/install/krakenuniq" --db "$D/wrapper-db" --database-outer \
      --output "$D/wrapper-default.out" --report-file "$D/wrapper-default.report" \
      "$D/batch.fa" 2>/dev/null
    cmp "$D/wrapper-batch.out" "$D/wrapper-default.out"
    for strategy in batch outer; do
      outer=()
      if [ "$strategy" = outer ]; then outer=(--database-outer); fi
      "$D/install/krakenuniq" --db "$D/wrapper-db" --batch-size 6000 "${outer[@]}" \
        --paired --output "$D/paired-$strategy.out" --report-file "$D/paired-$strategy.report" \
        "$D/batch.fq" "$D/batch.fq" 2>/dev/null
      grep -v '^#' "$D/paired-$strategy.report" > "$D/paired-$strategy.data"
    done
    cmp "$D/paired-batch.out" "$D/paired-outer.out"
    cmp "$D/paired-batch.data" "$D/paired-outer.data"
    echo "database-outer output matches batching for multiple input files"
    if "$SRC/classify" -B nonsense 2>"$D/invalid.err"; then
      echo "invalid batch size accepted"; exit 1
    fi
    grep -q "invalid batch size" "$D/invalid.err"
    echo "hierarchical batch output matches the per-read path"
    if [ "${KRAKEN_TEST_WIDE_INDEX:-0}" = 1 ]; then
      # Only the queried offsets need physical storage in this sparse fixture.
      python3 - "$D" <<'PYEOF'
import random, struct, sys
from pathlib import Path
d = Path(sys.argv[1])
k, nt = 31, 17
mask = (1 << (2 * nt)) - 1
xor = 0xe37e28c4271b5a2d & mask
rng = random.Random(391)
def revcomp(n, length):
    result = 0
    for _ in range(length):
        result = (result << 2) | (3 - (n & 3))
        n >>= 2
    return result
while True:
    key = rng.getrandbits(2 * k)
    key = min(key, revcomp(key, k))
    windows = [(key >> (2 * i)) & mask for i in range(k - nt + 1)]
    bin_key = min(min(w, revcomp(w, nt)) ^ xor for w in windows)
    if bin_key > (1 << 32):
        break
header_size = 72 + 2 * (4 + 8 * 2 * k)
header = bytearray((d / "db.kdb").read_bytes()[:header_size])
struct.pack_into("<Q", header, 48, 1)
(d / "wide.kdb").write_bytes(header + struct.pack("<QI", key, 100))
with open(d / "wide.idx", "wb") as f:
    f.write(b"KRAKIX2" + bytes([nt]))
    f.truncate(8 + 8 * ((1 << (2 * nt)) + 1))
    f.seek(8 + 8 * (bin_key + 1))
    f.write(struct.pack("<Q", 1))
seq = "".join("ACGT"[(key >> (2 * i)) & 3] for i in reversed(range(k)))
(d / "wide.fa").write_text("".join(f">wide{i}\n{seq}\n" for i in range(100)) + ">ambig\n" + "N" * k + "\n")
(d / "wide.expected").write_text("".join(f"C\twide{i}\t100\t31\t100:1\n" for i in range(100)) + "U\tambig\t0\t31\tA:1\n")
PYEOF
      for mode in sorted unsorted; do
        hint=()
        if [ "$mode" = unsorted ]; then hint=(-H); fi
        for strategy in batched outer; do
          outer=()
          if [ "$strategy" = outer ]; then outer=(-D); fi
          KRAKEN_BATCH_FORCE_HINTS=1 "$SRC/classify" -d "$D/db.kdbc" -i "$D/db.idxc" \
            -d "$D/wide.kdb" -i "$D/wide.idx" -d "$D/db.kdbc" -i "$D/db.idxc" \
            -a "$D/taxDB" -t "${threads[-1]}" -B 1M "${hint[@]}" "${outer[@]}" -r off \
            -o "$D/wide.out" "$D/wide.fa" 2>/dev/null
          cmp "$D/wide.expected" "$D/wide.out"
        done
      done
      echo "wide legacy index and mixed-width databases match expected output"
    fi
  fi
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
