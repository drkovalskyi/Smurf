
#!/bin/bash

#
# script to skim smurf ntuples
# choose from WW preselection, and PassFail Selections
# 

INPUTDIR=$1
OUTPUTDIR=$2
SELECTION=$3

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./skim.sh   INPUTDIR OUTPUTDIR
        INPUTDIR - location of smurf ntuples to skim (e.g. /smurf/data/Run2011_Spring11_SmurfV3/mitf-alljets/)
        OUTPUTDIR - location to output skimmed ntuples
        SELECTION - selection, choose from WW, PassFail"
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

# loop over root files in input dir
# and do the skim root script

if [ "$SELECTION" == 'WW' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
hww150.root
hww155.root
EOF
fi
data.root
data_3l.root
zz.root
wz.root
ttbar.root
tw.root
qqww.root
ggww.root
wjets.root
dyll.root
wgamma.root
wglll.root
hww110.root
hww115.root
hww120.root
hww125.root
hww130.root
hww135.root
hww140.root
hww145.root
hww150.root
hww155.root
hww160.root
hww170.root
hww180.root
hww190.root
hww200.root
hww250.root
hww300.root
hww350.root
hww400.root
hww450.root
hww500.root
hww550.root
hww600.root
hww700.root
hww800.root
hww900.root
hww1000.root
wwmcnlo.root
wwmcnloup.root
wwmcnlodown.root
ttbar_powheg.root



if [ "$SELECTION" == 'PassFail' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data.root
zz.root
wz.root
ttbar.root
tw.root
qqww.root
ggww.root
wgamma.root
wglll.root
dyll.root
wjets.root
EOF
fi


# Do the skimming...
for FILE in `cat list_samples.txt` ; do
    for JETBIN in 0 1 2 ; do 
	outputdir=$OUTPUTDIR/$SELECTION/${JETBIN}j/
	if [ "$SELECTION" == 'PassFail' ]; then
	    outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	fi
	mkdir -p $outputdir
	echo doing "root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",$JETBIN\);"
	root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",$JETBIN\);
    done
done

 
