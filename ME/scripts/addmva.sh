#!/bin/bash

# this scripts add the mva output to the existing smurfntuples
# USAGE: ./addmva.sh njets inputdir outputdir

export NJETS=$1;
export INPUTDIR=$2
export OUTPUTDIR=$3


mkdir -p log/

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./addmva.sh njets inputdir outputdir
        njet - njet bin e.g. 0, 1 or 2
        inputdir - input directories.. /smurf/yygao/data/EPS/WW/
        outputdir - input directories.. /smurf/yygao/data/EPS/WW/"
    exit 1
fi

#for MASS in 115 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600; do
for MASS in 130; do
    echo doing $MASS
    cd ../../MVA
    nohup ./add_bdt.sh $NJETS $MASS $INPUTDIR $OUTPUTDIR >& ../ME/scripts/log/addbdt_${MASS}_${NJETS}j.log 
    cd ../ME/scripts/
done
