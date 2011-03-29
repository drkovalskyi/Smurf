#!/bin/bash

# input files should be available in "data" folder
# ln -s /smurf/ceballos/data_39x data

export NJETS=$1;
export MH=$2;

export TAG=ntuples_${MH}train_${NJETS}jets;
export METHODS=Likelihood,KNN,BDT,BDTD,MLPBNN;

### Training: change done hand made, it's an expert option
export DO_TRAINING=1;
if [ ${DO_TRAINING} == "1" ]; then
  export SIG_TRAIN=data/histo_hww${MH}_std_pu11_randomized.all.root;
  export BKG_TRAIN=data/histo_w10-vv-mg-v8-pu11_all_noskim_normalized.root;
  mkdir -p weights;
  root -l -q -b trainMVA_smurf.C+\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data/histo_hww${MH}_std_pu11_randomized.all.root
data/histo_hww_std_pu_randomized.all.root
data/histo_data_all_noskim_normalized.root
EOF

mkdir -p output;
for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  root -l -q -b evaluateMVA_smurf.C+\(\"${dataset}\",${MH},\"${METHODS}\",\"${TAG}\"\);
done
rm -f list_samples.txt;
