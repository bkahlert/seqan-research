#!/bin/bash
if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` {fastqc_data.txt}"
  exit
fi
csplit --silent $1 /'^>>.[a-z]'/ {*}
grep -v '>>END_MODULE' xx01 > summary.tsv
grep -v '>>END_MODULE' xx02 > per_pos_qual.tsv
grep -v '>>END_MODULE' xx03 > seq_qual.tsv
grep -v '>>END_MODULE' xx04 > pos_nuc_dist.tsv
grep -v '>>END_MODULE' xx05 > pos_gc_content.tsv
grep -v '>>END_MODULE' xx06 > seq_gc_content.tsv
grep -v '>>END_MODULE' xx07 > pos_n_content.tsv
grep -v '>>END_MODULE' xx08 > readlength.tsv
grep -v '>>END_MODULE' xx09 > seq_dupl.tsv
grep -v '>>END_MODULE' xx11 > kmer_content.tsv
rm xx*
