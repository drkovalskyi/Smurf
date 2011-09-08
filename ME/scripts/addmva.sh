#!/bin/bash

# this scripts add the mva output to the existing smurfntuples
# USAGE: ./addmva.sh njets inputdir outputdir

export NJETS=$1;
export INPUTDIR=$2
export OUTPUTDIR=$3
export ADDME=$4

mkdir -p log/

if [ ! $# -eq 4 ]; then
    echo "USAGE: ./addmva.sh njets inputdir outputdir
        njet - njet bin e.g. 0, 1 or 2
        inputdir - input directories.. /smurf/yygao/data/LP2011/WW/
        outputdir - input directories.. /smurf/yygao/data/LP2011/WW/
	addme - set to 1 if you want to add me output before the bdt"
    exit 1
fi


#  self-protection, the ME code is invoked only for the 0-jet bin
if [ ${NJETS} == "0" ] && [ ${ADDME} == "1" ]; then
    export MEFLAG=1
else
    export MEFLAG=0
fi

## add me if specified
if [ ${MEFLAG} == "1" ]; then
    echo "adding the LR of ME..."
    ./addLR.sh ${INPUTDIR}/${NJETS}j/ ${INPUTDIR}/${NJETS}j/ME/
fi


## add bdt
echo "adding the BDT..."
for MASS in 115 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600; do
    echo doing $MASS
    cd ../../MVA
    nohup ./add_bdt.sh $NJETS $MASS $INPUTDIR $OUTPUTDIR $MEFLAG >& ../ME/scripts/log/addbdt_${MASS}_${NJETS}j.log 
    cd ../ME/scripts/
done
