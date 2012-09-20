#!/bin/bash

#
# a script to combine batch job output
# and evaluate hypothesis separation between M0 and M1
#

if [ ! $# -eq 3 ]; then
    echo "
USAGE: ./postprocessSeparation.sh TASK M0 M1
    TASK  - The name of the task directory
    M0    - The card for model0
    M1    - The card for model1
"
    exit 1
fi

TASK=$1
M0=$2
M1=$3
M0NAME=`echo ${M0} | sed 's/\.txt//'`
M1NAME=`echo ${M1} | sed 's/\.txt//'`
WORKDIR=`pwd`
LANDSDIR=${WORKDIR}/../LandS

# combine output
cd ${TASK}
hadd LL_toyM0_fitM0_maxlls_tree.root output/LL_toy${M0NAME}_fit${M0NAME}_seed*_maxlls_tree.root
hadd LL_toyM0_fitM1_maxlls_tree.root output/LL_toy${M0NAME}_fit${M1NAME}_seed*_maxlls_tree.root
hadd LL_toyM1_fitM0_maxlls_tree.root output/LL_toy${M1NAME}_fit${M0NAME}_seed*_maxlls_tree.root
hadd LL_toyM1_fitM1_maxlls_tree.root output/LL_toy${M1NAME}_fit${M1NAME}_seed*_maxlls_tree.root

# do hypoSeparation
root -l ${LANDSDIR}/test/hypoSeparation.C\(\"${M0NAME}\",\"LL_toyM0_fitM0_maxlls_tree.root\",\"LL_toyM0_fitM1_maxlls_tree.root\",\"${M1NAME}\",\"LL_toyM1_fitM0_maxlls_tree.root\",\"LL_toyM1_fitM1_maxlls_tree.root\"\)

cd ${WORKDIR}

