#!/bin/bash

#
# script to skim smurf ntuples
#

INPUTDIR=$1
OUTPUTDIR=$2
SELECTION=$3

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./skim.sh   INPUTDIR OUTPUTDIR
        INPUTDIR - location of smurf ntuples to skim (e.g. /smurf/data/Run2011_Spring11_SmurfV3/mitf-alljets/)
        OUTPUTDIR - location to output skimmed ntuples
        SELECTION - selection, choose from ZZ, WW or PassFail "
    exit 1
fi


# check that directories exist

if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir doesnt exist!
        exit 2
fi

if [ ! -d $OUTPUTDIR ]; then
    echo Error: Output dir doesnt exist!
    exit 3
fi

mkdir -p $OUTPUTDIR/$SELECTION

# loop over root files in input dir
# and do the skim root script
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data-met20-1092ipb.root
hww115.root
hww120.root
hww130.root
hww140.root
hww150.root
hww160.root
hww170.root
hww180.root
hww190.root
hww200.root
hww250.root
hww300.root
zz.root
wz.root
ttbar.root
tw.root
qqww.root
ggww.root
wjets.root
wgamma.root
dyee.root
dymm.root
dytt.root
EOF

#hzz200.root
#hzz250.root
#hzz300.root
#hzz400.root

#for FILE in `ls $INPUTDIR | grep lfake.root`
for FILE in `cat list_samples.txt` ; do
    outputdir=$OUTPUTDIR/$SELECTION/
    echo doing "root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\"\);"
    root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\"\);
done

