#!/bin/bash

#
# script to create new inputs to test the shapeStat option in LandS
#

INPUTDIR=$1

if [ ! $# -eq 1 ]; then
    echo "USAGE: ./shapestat.sh inputdir
    inputdir - directory of the cards, e.g. ana_v6_1500pb_LP_FINAL"
    exit 1
fi

# check that directories exist

if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir ${INPUTCARD} containing the orignial cards doesnt exist!
        exit 2
fi


rm temp.txt

#for MH in 115 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600; do
for MH in 140; do
    export INPUTCARD=$INPUTDIR/$MH
    # create the split card directory for each mass point
    export SHAPESTAT=$INPUTCARD/shapeStat/
    echo creating directory "$SHAPESTAT"
    mkdir -p $SHAPESTAT
    # create the new cards according to njet and flavor
    # for NJETS in 0 1; do
    for NJETS in 0; do
	# for FLAVOR in sf of; do
	for FLAVOR in of; do
	    cp $INPUTCARD/hww${FLAVOR}_${NJETS}j.input.root $SHAPESTAT/
	    echo "Doing root -l -b -q shapeStat.C\($MH,$NJETS,\"${FLAVOR}\",\"${INPUTDIR}\"\);"
	    root -l -b -q shapeStat.C\($MH,$NJETS,\"${FLAVOR}\",\"${INPUTDIR}\"\)
	done
    done
    
    echo -e "$MH $MH/shapeStat/hwwsf_0j_shape.txt $MH/shapeStat/hwwof_0j_shape.txt $MH/shapeStat/hwwsf_1j_shape.txt $MH/shapeStat/hwwof_1j_shape.txt $MH/hww_2j_cut.txt" >> temp.txt
done

mv temp.txt $INPUTDIR/limits_nj_shape_shapeStat.txt



