#!/bin/bash

for PROC in data hww115 hww120 hww130 hww140 hww150 hww160 hww170 hww180 hww190 hww200 hww250 hww300 qqww ggww wjets wjets_data wjets_PassFail ttbar tw dyee dymm dytt wgamma wz zz ww2l_pythia; do
    echo "Checking " $PROC >> checkBatchjob.log 
    ./check.sh $PROC
done
