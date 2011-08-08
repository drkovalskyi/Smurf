#!/bin/bash

# Example:
# ./run_mva.sh 0 130 projlumi datalumi
#
# 1. input files should be available in "data" folder
#  IMPORTANT: this assumes that the data file already inclueds the differential cross-sections for ME method
#    ln -s /smurf/ceballos/tmva/weights weights
# 

# 2. To add the BDT output to the smurfntuples
#    Use 5.28 version root, this is currently available at lxplus

export NJETS=$1;
export MH=$2;
export INPUTDIR=$3
export OUTPUTDIR=$4


if [ ! $# -eq 4 ]; then
	echo "USAGE: ./run_limit.sh njets mH lumi
        njet - njet bin e.g. 0, 1 or 2
        mH   - SM Higgs mass hypothesis e.g. 120
        inputdir - input directories.. /smurf/yygao/data/EPS/WW/
        outputdir - input directories.. /smurf/yygao/data/EPS/WW/"
	exit 1
fi

# append the inputdir by the number of jets
INPUTDIR=${INPUTDIR}/${NJETS}j/

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

rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data-met20-1092ipb
wgamma
hww${MH}
qqww
ggww
wjets
wjets_data
ttbar
tw
stop
ttop
wz
zz
dyee
dymm
dytt
EOF


# ===========================================
# Fill the smurfntuples with the BDT
# ===========================================

export evaluateMVAFile=evaluateMVA_smurf.C++;
if [ ${NJETS} == "1" ]; then
  export evaluateMVAFile=evaluateMVA_smurf1.C++;
fi

for i in `cat list_samples.txt` ; do
    dataset=${i%%,*};
    echo "filling MVA information in sample: "  $sample
    ./root-5.28.sh -q -b ${evaluateMVAFile}\(\"data/${dataset}.root\",${MH},\"${METHODS}\",\"${TAG}\"\);
    mv $OUTPUTDIR/${TAG}_${dataset}.root $OUTPUTDIR/${dataset}_${NJETS}j.root
done

#tidy up
rm -f list_samples.txt;
