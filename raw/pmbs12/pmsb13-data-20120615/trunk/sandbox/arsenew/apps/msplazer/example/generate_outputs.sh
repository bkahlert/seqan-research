#!/bin/sh
# Output generation script for msplazer

MSPLAZER=../msplazer
#MSPLAZER=../../../../build/RelDebug/sandbox/ktrappe/apps/msplazer/msplazer

# ============================================================
# Calling MSplazer
# ============================================================

# Calling first Stellar (externally) and then call MSplazer with option --matchfile (-m)
#./stellar -d chr10.fa -q 20_contigs_chr10_sim_uniqIds.fa -o 20_contigs_chr10_sim_uniqIds.gff -l 30 -k 16 -e 0.03 > 20_contigs_chr10_sim_uniqIds.out
${MSPLAZER} -d chr10.fa -q 20_contigs_chr10_sim_uniqIds.fa -m 20_contigs_chr10_sim_uniqIds.gff -st 1 -v -i externalStellar/ -j 20_contigs_chr10_sim_uniqIds_ > externalStellar/20_contigs_chr10_sim_uniqIds.stdout

# Calling only MSplazer which runs Stellar internally using the Stellar options -l 30 (minimal match length), -k 16 (q-gram length), -e 0.03 (error rate)
${MSPLAZER} -d chr10.fa -q 20_contigs_chr10_sim_uniqIds.fa -l 30 -k 16 -e 0.03 -st 1 -v -i externalStellar/ -j 20_contigs_chr10_sim_uniqIds_internalStellar_ > internalStellar/20_contigs_chr10_sim_uniqIds_internalStellar.stdout
