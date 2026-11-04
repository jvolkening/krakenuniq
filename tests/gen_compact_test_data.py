#!/usr/bin/env python3
"""Generate a small synthetic KrakenUniq database (unsorted), a matching
taxonomy (names.dmp/nodes.dmp) and a FASTA file of reads for testing
compact_db. Usage: gen_compact_test_data.py OUTDIR [K] [SEED]"""
import os, random, struct, sys

outdir = sys.argv[1]
K = int(sys.argv[2]) if len(sys.argv) > 2 else 31
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
os.makedirs(outdir, exist_ok=True)
GENOME_LEN = int(os.environ.get('GENOME_LEN', 120000))
N_READS = int(os.environ.get('N_READS', 150))

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
MASK = (1 << (2 * K)) - 1


def revcomp(s):
    return ''.join(COMP[c] for c in reversed(s))


def canon(s):
    a = b = 0
    for c in s:
        a = (a << 2) | CODE[c]
    for c in revcomp(s):
        b = (b << 2) | CODE[c]
    return min(a, b)


# taxonomy: 1 -> genera 10..12 -> species; one species has the largest
# supported taxid (2^26 - 1)
BIG = (1 << 26) - 1
species = [100, 101, 102, 103, 104, BIG]
genus_of = {100: 10, 101: 10, 102: 11, 103: 11, 104: 12, BIG: 12}
genomes = {}
for sp in species:
    genomes[sp] = ''.join(random.choice('ACGT') for _ in range(GENOME_LEN))
# shared segments produce k-mers with LCA = genus or root
seg = genomes[100][1000:21000]
genomes[101] = genomes[101][:5000] + seg + genomes[101][25000:]
seg = genomes[100][30000:40000]
genomes[102] = genomes[102][:5000] + seg + genomes[102][15000:]

owners = {}
for sp, g in genomes.items():
    for i in range(len(g) - K + 1):
        owners.setdefault(canon(g[i:i + K]), set()).add(sp)


def lca(sps):
    if len(sps) == 1:
        return next(iter(sps))
    gs = {genus_of[s] for s in sps}
    return gs.pop() if len(gs) == 1 else 1


records = [(kmer, lca(sps)) for kmer, sps in owners.items()]
random.shuffle(records)

key_bits = 2 * K
hdr_size = 72 + 2 * (4 + 8 * key_bits)
hdr = bytearray(hdr_size)
hdr[0:8] = b'JFLISTDN'
struct.pack_into('<Q', hdr, 8, key_bits)
struct.pack_into('<Q', hdr, 16, 4)
struct.pack_into('<Q', hdr, 48, len(records))
key_len = (key_bits + 7) // 8
with open(os.path.join(outdir, 'unsorted.kdb'), 'wb') as f:
    f.write(hdr)
    for kmer, tax in records:
        f.write(kmer.to_bytes(key_len, 'little') + struct.pack('<I', tax))

with open(os.path.join(outdir, 'nodes.dmp'), 'w') as f:
    f.write('1\t|\t1\t|\tno rank\t|\n')
    for g in (10, 11, 12):
        f.write('%d\t|\t1\t|\tgenus\t|\n' % g)
    for sp in species:
        f.write('%d\t|\t%d\t|\tspecies\t|\n' % (sp, genus_of[sp]))
with open(os.path.join(outdir, 'names.dmp'), 'w') as f:
    f.write('1\t|\troot\t|\t\t|\tscientific name\t|\n')
    for t in [10, 11, 12] + species:
        f.write('%d\t|\ttaxon%d\t|\t\t|\tscientific name\t|\n' % (t, t))

with open(os.path.join(outdir, 'reads.fa'), 'w') as f:
    n = 0
    for sp, g in genomes.items():
        for _ in range(N_READS):
            p = random.randrange(len(g) - 150)
            r = list(g[p:p + 150])
            for j in range(len(r)):
                if random.random() < 0.01:
                    r[j] = random.choice('ACGT')
            r = ''.join(r)
            if random.random() < 0.5:
                r = revcomp(r)
            f.write('>read%d\n%s\n' % (n, r))
            n += 1
    for _ in range(100):
        f.write('>rand%d\n%s\n' % (n, ''.join(random.choice('ACGTN') for _ in range(150))))
        n += 1
print('%d k-mers, k=%d' % (len(records), K))
