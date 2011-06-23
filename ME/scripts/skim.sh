#!/bin/bash

#
# script to skim smurf ntuples
#

if [ ! $# -eq 2 ]; then
    echo "USAGE: ./skim.sh   INPUTDIR OUTPUTDIR
        INPUTDIR - location of smurf ntuples to skim (e.g. /smurf/data/Run2011_Spring11_SmurfV3/mitf-alljets/)
        OUTPUTDIR - location to output skimmed ntuples"
    exit 1
fi

INPUTDIR=$1
OUTPUTDIR=$2

# check that directories exist

if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir doesnt exist!
        exit 2
fi

if [ ! -d $OUTPUTDIR ]; then
        echo Error: Output dir doesnt exist!
        exit 3
fi

# loop over root files in input dir
# and do the skim root script

for FILE in `ls $INPUTDIR | grep lfake.root`
do
        root -l -b -q smurfproducer.C\(\"$INPUTDIR\",\"$FILE\",\"$OUTPUTDIR\"\)
done

