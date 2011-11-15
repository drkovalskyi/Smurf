#!/bin/bash

# You need to do kinit -A $USER@FNAL.GOV first to excute this script

export OUTDIR=/smurf/yygao/data/LP2011/WW/0j/ME/
#export OUTDIR=./

#for PROC in data hww115 hww120 hww130 hww140 hww150 hww160 hww170 hww180 hww190 hww200 hww250 hww300 hww350 hww400 hww450 hww500 hww550 hww600 qqww ggww wjets wjets_data wjets_PassFail ttbar tw dyee dymm dytt wgamma wz zz; do
#for PROC in data; do
for PROC in hww180 hww190 hww200 hww250 hww300; do
    scp ygao@cmslpc-sl5.fnal.gov:/uscms_data/d2/ygao/SmurfV6/WW/LP2011/{$PROC}.root $OUTDIR/${PROC}.root
#    cp $OUTDIR/${PROC}.root $OUTDIR/${PROC}_dXsec.root
done