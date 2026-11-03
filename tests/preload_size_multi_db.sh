#!/bin/bash
# Checks that --preload-size (classify -x) with several databases (hierarchical search)
# produces output identical to classification without --preload-size.
# Usage: tests/preload_size_multi_db.sh KRAKENUNIQ_INSTALL_DIR
set -euo pipefail
BIN=$(cd "${1:?Usage: $0 KRAKENUNIQ_INSTALL_DIR}" && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

python3 - <<'PY'
import random, os
random.seed(1)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
def mutate(s):
    return ''.join(random.choice('ACGT') if random.random() < 0.01 else c for c in s)
g = {n: rnd(20000) for n in ('g1', 'g2', 'g3')}

for db in ('DB1', 'DB2'):
    os.makedirs(db + '/taxonomy')
    os.makedirs(db + '/library')
    with open(db + '/taxonomy/nodes.dmp', 'w') as f:
        for tid, parent, rank in ((1, 1, 'no rank'), (2, 1, 'superkingdom'), (10, 2, 'species'), (20, 2, 'species'), (30, 2, 'species')):
            f.write('%d\t|\t%d\t|\t%s\t|\n' % (tid, parent, rank))
    with open(db + '/taxonomy/names.dmp', 'w') as f:
        for tid, name in ((1, 'root'), (2, 'Bacteria'), (10, 'Species one'), (20, 'Species two'), (30, 'Species three')):
            f.write('%d\t|\t%s\t|\t\t|\tscientific name\t|\n' % (tid, name))

# g1 is in both databases, with different taxa, to test that DB1 takes precedence
lib = {'DB1': [('g1', 10, g['g1']), ('g3', 30, g['g3'])],
       'DB2': [('g2', 20, g['g2']), ('g1b', 20, g['g1']), ('g3b', 30, g['g3'][:10000])]}
for db, entries in lib.items():
    with open(db + '/library/seqs.fa', 'w') as fa, open(db + '/library/seqs.map', 'w') as mp:
        for name, tid, seq in entries:
            fa.write('>%s\n%s\n' % (name, seq))
            mp.write('%s\t%d\t%s\n' % (name, tid, name))

reads = []
names = list(g)
for i in range(4000):
    kind = i % 4
    if kind == 3:
        s = rnd(100)
    elif kind == 2:  # chimeric read
        a, b = random.sample(names, 2)
        s = g[a][random.randrange(19950):][:50]
        s += g[b][random.randrange(19950):][:50]
    else:
        src = g[random.choice(names)]
        p = random.randrange(len(src) - 100)
        s = mutate(src[p:p + 100])
        if random.random() < 0.5:
            s = rc(s)
    if random.random() < 0.02:
        s = s[:40] + 'N' + s[41:]
    reads.append(s)
with open('reads.fq', 'w') as f:
    for i, s in enumerate(reads):
        f.write('@r%d\n%s\n+\n%s\n' % (i, s, 'I' * len(s)))
PY

for db in DB1 DB2; do
  "$BIN/krakenuniq-build" --db $db --kmer-len 21 --minimizer-len 8 --threads 2 --jellyfish-hash-size 10M >$db.build.log 2>&1 || { cat $db.build.log; exit 1; }
done
ls -l DB1/database.* DB2/database.*

run() {
  local tag=$1; shift
  "$BIN/krakenuniq" --db DB1 --db DB2 --threads 2 --report-file $tag.report --output $tag.out "$@" reads.fq 2>$tag.log \
    || { cat $tag.log; return 1; }
  sort $tag.out -o $tag.out.sorted
}

run ref
echo "Calls without --preload-size (taxon: reads):"
cut -f1,3 ref.out | sort | uniq -c

status=0
for size in 600K 800K 2M 100M; do
  run x$size --preload-size $size
  if diff -q ref.out.sorted x$size.out.sorted >/dev/null && diff -q ref.report x$size.report >/dev/null; then
    echo "OK: --preload-size $size gives identical output"
  else
    echo "FAIL: --preload-size $size differs"
    diff ref.out.sorted x$size.out.sorted | head -n 10 || true
    diff ref.report x$size.report | head -n 10 || true
    status=1
  fi
done

# a buffer too small for even a single minimizer bin must fail with an error message
if run tiny --preload-size 1K; then
  echo "FAIL: --preload-size 1K should have failed"; status=1
else
  grep -q 'too small' tiny.log && echo "OK: too small --preload-size is rejected" || { echo "FAIL: no error message"; cat tiny.log; status=1; }
fi
exit $status
