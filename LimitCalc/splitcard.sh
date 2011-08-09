#!/bin/bash

#
# script to split the shape analysis input into cards
#

if [ ! $# -eq 1 ]; then
    echo "USAGE: ./splitcard.sh inputdir
    inputdir - directory of the cards, e.g. EPS2011_StatTest/"
    exit 1
fi

INPUTDIR=$1


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
    export SPLITCARDS=$INPUTCARD/splitcards/
    echo "$SPLITCARDS"
    mkdir -p $SPLITCARDS
    
    # split the cards according to njet and flavor
    for NJETS in 0 1; do
	for FLAVOR in sf of; do
	    echo "Doing root -l -b -q splitcard.C\($MH,$NJETS,\"${FLAVOR}\",\"${INPUTDIR}\"\);"
	    root -l -b -q splitcard.C\($MH,$NJETS,\"${FLAVOR}\",\"${INPUTDIR}\"\)
	done
    done
    
    echo -e "$MH $MH/splitcards/hwwsf_0j_shape_Bin*.txt $MH/splitcards/hwwof_0j_shape_Bin*.txt $MH/splitcards/hwwsf_1j_shape_Bin*.txt $MH/splitcards/hwwof_1j_shape_Bin*.txt $MH/hww_2j_cut.txt" >> temp.txt
done

mv temp.txt $INPUTDIR/limits_nj_shape_split.txt



