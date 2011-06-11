#!/bin/bash

# input files should be available in "data" folder
#
#ln -s /data/smurf/ceballos/data_41x data
#mkdir -p /data/smurf/ceballos/tmva/weights
#mkdir -p /data/smurf/ceballos/tmva/output
#ln -s /data/smurf/ceballos/tmva/weights
#ln -s /data/smurf/ceballos/tmva/output

export NJETS=$1;
export MH=$2;

export TAG=ntuples_${MH}train_${NJETS}jets;
export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${NJETS} == "1" ]; then
  export trainMVA_smurfFile=trainMVA_smurf1.C+;
fi
export DO_TRAINING=0;
if [ ${DO_TRAINING} == "1" ]; then
  export SIG_TRAIN=data/hww${MH}.root;
  export BKG_TRAIN=data/histo_hww_train_pu.root;
  #export BKG_TRAIN=data/ww2l_pythia.root;
  mkdir -p weights;
  root -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
### samples must be in "data" folder
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data/hww${MH}.root
data/data.root
data/histo_hww_std_data_pu.root
data/histo_hww_std_pu.root
EOF

export evaluateMVAFile=evaluateMVA_smurf.C+;
if [ ${NJETS} == "1" ]; then
  export evaluateMVAFile=evaluateMVA_smurf1.C+;
fi

for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"${METHODS}\",\"${TAG}\"\);
done
rm -f list_samples.txt;
