#!/bin/bash

# 1. IMPORTANT: this assumes that the data file already inclueds the differential cross-sections for ME method
#    ln -s /smurf/ceballos/tmva/weights weights

# 2. To add the BDT output to the smurfntuples
#    Use 5.28 version root, this is currently available at lxplus

export NJETS=$1;
export MH=$2;
export INPUTDIR=$3
export OUTPUTDIR=$4
export MEFLAG=$5
 
if [ ! $# -eq 5 ]; then
	echo "USAGE: ./add_bdt.sh njets mH inputdir outputdir meflag
        njet - njet bin e.g. 0, 1 or 2
        mH   - SM Higgs mass hypothesis e.g. 120
        inputdir - input directories.. /smurf/yygao/data/EPS/WW/
        outputdir - input directories.. /smurf/yygao/data/EPS/WW/
	meflag - set to 1 if analyzing the inputs from ME code"
	exit 1
fi

# append the inputdir by the number of jets
if [ ${NJETS} == "0" ] && [ ${MEFLAG} == "1" ]; then
	INPUTDIR=${INPUTDIR}/${NJETS}j/ME/
else
    INPUTDIR=${INPUTDIR}/${NJETS}j/
fi

### samples must be in "data" folder
rm -f data
ln -s $INPUTDIR data

# check that input directories exist
if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir doesnt exist!
        exit 2
fi

# make output directory
OUTPUTDIR=${OUTPUTDIR}/mva/$MH/
mkdir -p $OUTPUTDIR
rm -f output
ln -s $OUTPUTDIR output

# this is the prefix added to the BDT added ntuples
export TAG=ntuples_${MH}train_${NJETS}jets;
#export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;
export METHODS=BDTG;

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
# define a list of the files to analyze
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data
wgamma
hww${MH}
qqww
ggww
wjets
wjets_data
wjets_PassFail
ttbar
tw
wz
zz
dyee
dymm
dytt
EOF

# ===========================================
# Fill the smurfntuples with the BDT
# ===========================================

export evaluateMVAFile=evaluateMVA_smurf.C+;
if [ ${NJETS} == "1" ]; then
  export evaluateMVAFile=evaluateMVA_smurf1.C+;
fi

for i in `cat list_samples.txt` ; do
    dataset=${i%%,*};
    echo "filling MVA information in sample: "  $sample
    ./root-5.28.sh -q -b ${evaluateMVAFile}\(\"data/${dataset}.root\",${MH},\"${METHODS}\",\"${TAG}\"\);
    mv $OUTPUTDIR/${TAG}_${dataset}.root $OUTPUTDIR/${dataset}_${NJETS}j.root
done

rm -f list_samples.txt;

