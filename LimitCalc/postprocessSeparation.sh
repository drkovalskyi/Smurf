#!/bin/bash

#
# a script to combine batch job output
# and evaluate hypothesis separation between M0 and M1
#

if [ ! $# -eq 1 ]; then
    echo "
USAGE: ./postprocessSeparation.sh TASK M0 M1
    TASK  - The name of the task directory
"
    exit 1
fi

TASK=$1
M0NAME=SMHiggs
M1NAME=Graviton
WORKDIR=`pwd`
MACRODIR=${WORKDIR}/../../


# combine output
cd ${TASK}
hadd LL_toyM0_fitM0_maxlls_tree.root output/LL_toy${M0NAME}_fit${M0NAME}_seed*_maxlls_tree.root
hadd LL_toyM0_fitM1_maxlls_tree.root output/LL_toy${M0NAME}_fit${M1NAME}_seed*_maxlls_tree.root
hadd LL_toyM1_fitM0_maxlls_tree.root output/LL_toy${M1NAME}_fit${M0NAME}_seed*_maxlls_tree.root
hadd LL_toyM1_fitM1_maxlls_tree.root output/LL_toy${M1NAME}_fit${M1NAME}_seed*_maxlls_tree.root

# do hypoSeparation
export NBINS=200
export XMIN=-30
export XMAX=30
root -l -b -q ${MACRODIR}/hypoSeparation.C\(\"${M0NAME}\",\"LL_toyM0_fitM0_maxlls_tree.root\",\"LL_toyM0_fitM1_maxlls_tree.root\",\"${M1NAME}\",\"LL_toyM1_fitM0_maxlls_tree.root\",\"LL_toyM1_fitM1_maxlls_tree.root\",$NBINS,$XMIN,$XMAX\)

cd ${WORKDIR}

