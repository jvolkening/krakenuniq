#! /bin/sh

cd tests
. ./compat.sh

K=9

# Deterministic sequences: a long one (spans several parser buffers) with
# a few Ns, and a short one (shorter than K).
awk -v n=20000 'BEGIN {
  s = 12345; split("A C G T", b, " ")
  print ">long"
  line = ""
  for(i = 0; i < n; i++) {
    s = (s * 1103515245 + 12345) % 2147483648
    c = b[int(s / 65536) % 4 + 1]
    if(i % 3001 == 1500) c = "N"
    line = line c
    if(length(line) == 60) { print line; line = "" }
  }
  if(line != "") print line
  print ">short"
  print "ACGTN"
}' > ${pref}.fa

# Reference: count every step-th window of each run of valid bases
expected () {
    awk -v k=$K -v step=$1 '
      /^>/ { run = 0; seq = ""; next }
      {
        for(i = 1; i <= length($0); i++) {
          c = substr($0, i, 1); seq = seq c
          if(c !~ /[ACGT]/) { run = 0; seq = ""; continue }
          run++
          if(run >= k && (run - k) % step == 0) cnt[substr(seq, length(seq) - k + 1)]++
        }
      }
      END { for(m in cnt) print m, cnt[m] }' ${pref}.fa | sort
}

for step in 1 2 3 7; do
    echo "Count with --step $step" && \
        $JF count -m $K -s 1M -t $nCPUs --step $step --buffer-size 4096 -o ${pref}_$step -O ${pref}.fa && \
        $JF dump -c ${pref}_$step | sort > ${pref}_$step.got && \
        expected $step > ${pref}_$step.exp && \
        diff ${pref}_$step.exp ${pref}_$step.got
done

echo "Default step is 1" && \
    $JF count -m $K -s 1M -t $nCPUs --buffer-size 4096 -o ${pref}_def -O ${pref}.fa && \
    $JF dump -c ${pref}_def | sort | diff ${pref}_1.exp -

echo "Invalid step values are rejected"
for bad in 0 -1 x 1.5; do
    if $JF count -m $K -s 1M --step $bad -o ${pref}_bad -O ${pref}.fa > /dev/null 2>&1; then
        echo "--step $bad should fail"; exit 1
    fi
done
