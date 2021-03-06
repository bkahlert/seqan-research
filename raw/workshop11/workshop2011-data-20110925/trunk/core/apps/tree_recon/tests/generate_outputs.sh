#!/bin/sh
#
# The current directory is assumed to contain the files example{1,2,3}.dist
# with Phylip distance matrices.  The script will then run the tree_recon
# program configured below with different parameters on all files and generate
# output files for them.
#
# They can be used to compare against the output of new tree_recon programs.

# seqan-apps was a checkout of September 25, 2009 as advocated by T. Rausch.
TREE_RECON=../../../../../seqan-apps/projects/library/cmake/apps/tree_recon

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${TREE_RECON} -m example${i}.dist -o example${i}.out
    ${TREE_RECON} -m example${i}.dist -o example${i}.out
done

# Run with all building method options.
for i in 1 2 3; do
    for b in nj min max avg wavg; do
        echo ${TREE_RECON} -m example${i}.dist -b ${b} -o example${i}.${b}.out
        ${TREE_RECON} -m example${i}.dist -b ${b} -o example${i}.${b}.out
    done
done

# Run with all output formats
for i in 1 2 3; do
    for f in dot newick; do
        echo ${TREE_RECON} -m example${i}.dist -f ${f} -o example${i}.${f}.out
        ${TREE_RECON} -m example${i}.dist -f ${f} -o example${i}.${f}.out
    done
done
