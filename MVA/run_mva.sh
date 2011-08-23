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
export WEIGHTSONLY=$3;

export TAG=ntuples_${MH}train_${NJETS}jets;
export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;
#export TAG=TEST_ntuples_${MH}train_${NJETS}jets;
#export METHODS=KNN,Likelihood,LikelihoodD,MLPBNN,BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${NJETS} == "1" ]; then
  export trainMVA_smurfFile=trainMVA_smurf1.C+;
elif [ ${NJETS} == "hzz" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_hzz.C+;
  export NJETS=999;
fi

export DO_TRAINING=0;
if [ ${DO_TRAINING} == "1" ]; then
  #export SIG_TRAIN=data/hzz${MH}.root;
  #export BKG_TRAIN=data/zz.root;
  #export BKG_TRAIN=data/background_training.root;
  #export BKG_TRAIN=data/background_training.root;
  export SIG_TRAIN=data/hww${MH}.root;
  export BKG_TRAIN=data/ww2l_pythia.root;
  #export BKG_TRAIN=data/background42x.root;
  mkdir -p weights;
  root -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
### samples must be in "data" folder
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data/backgroundA_3l.root
data/backgroundA.root
data/backgroundA_skim1.root
data/backgroundA_skim2.root
data/backgroundB.root
data/backgroundB_skim1.root
data/backgroundB_skim2.root
data/data_2fake.root
data/data_2l.root
data/data_2l_skim1.root
data/data_2l_skim2.root
data/data_lfake.root
data/dyee.root
data/dymm.root
data/dytt.root
data/ggww.root
data/mc_l_fake_pu.root
data/qqww.root
data/ttbar_mg.root
data/ttbar.root
data/tw_ds.root
data/tw.root
data/wgamma.root
data/wjets.root
data/ww2l_pythia.root
data/wz_py.root
data/wz.root
data/zz.root
data/hww${MH}.root
data/hzz${MH}.root
EOF

export evaluateMVAFile=evaluateMVA_smurf_hww.C+;
if [ ${NJETS} == "999" ]; then
  export evaluateMVAFile=evaluateMVA_smurf_hzz.C+;
fi

for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  if [ ${WEIGHTSONLY} == "1" ]; then
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"\",\"\",\"\",${NJETS},1,0\);
  else
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,0\);
  fi
  
done
rm -f list_samples.txt;
