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
if [ ${MEFLAG} == "1" ]; then
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
export TAG=ntuples2012_PostICHEP_${MH}train_${NJETS}jets;
if [ ${NJETS} == "2" ]; then  
    TAG=ntuples2012_PostICHEP_125train_2jets;
fi
export METHODS=BDTG;


# If do the training...
export DO_TRAINING=0;
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${DO_TRAINING} == "1" ]; then
    export SIG_TRAIN=data/hww${MH}.root;
    export BKG_TRAIN=data/qqww${MH}.root;
    if [ ${NJETS} == "2" ]; then
	export SIG_TRAIN=data/hww${MH}.root;
	export BKG_TRAIN=data/top.root;
    fi  
    if [ ! -d $weights ]; then
	mkdir -p weights;
    fi
    ./root-5.28.sh -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi


### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
# define a list of the files to analyze
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data
qqww
ggww
ttbar_powheg
tw
wz
zz
dyll
wgamma
wglll
wjets
www
zgamma
www_PassFail
zgamma_PassFail
data_PassFail
qqww_PassFail
ggww_PassFail
wjets_PassFail
tw_PassFail
wz_PassFail
zz_PassFail
dyll_PassFail
wgamma_PassFail
wglll_PassFail
ttbar_powheg_PassFail
wwmcnlo
wwmcnloup
wwmcnlodown
ttbar
EOF


# ===========================================
# Fill the smurfntuples with the BDT
# ===========================================
export evaluateMVAFile=evaluateMVA_smurf_hww.C+;

# set running perid setting these depend on the lepton efficiency, fakerate and pu histograms 
# check the script directly to make sure all the histograms are read correctly

export PERIOD=2 

# add mva for the signal  
./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/hww${MH}.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",$PERIOD\);
mv $OUTPUTDIR/${TAG}_hww${MH}.root $OUTPUTDIR/hww${MH}_${NJETS}j.root

# add mva for all background
# add this for 0/1 Jet all mass points and only add 125 for the 2-jet bin
if [ ${NJETS} -lt 2 ] || [ $MH -eq 125 ]; then
    echo njets = $NJETS
    for i in `cat list_samples.txt` ; do
	dataset=${i%%,*};
	echo "filling MVA information in sample: "  $dataset
	PERIOD=2
	echo "period = $PERIOD"
	./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/${dataset}.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",$PERIOD\);
	mv $OUTPUTDIR/${TAG}_${dataset}.root $OUTPUTDIR/${dataset}_${NJETS}j.root
    done
fi

rm -f list_samples.txt;

