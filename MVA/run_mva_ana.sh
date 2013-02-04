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
export doMultiClass=1;

#export TAG=TEST_${MH}train_${NJETS}jets;
#export METHODS=Likelihood,BDT,BDTD,BDTG;
export TAG=ntuples2012_PostICHEP_${MH}train_${NJETS}jets;
#export TAG=TEST_ntuples_${MH}train_${NJETS}jets;
#export METHODS=BDT,Likelihood,BDTG,BDTD,BDTB,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,BoostedFisher,LikelihoodD,LikelihoodPCA,FDA_GA,RuleFit;
export METHODS=BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${doMultiClass} == "1" ]; then
  export TAG=ntuples2012_MultiClass_${MH}train_${NJETS}jets;
  export trainMVA_smurfFile=trainMVA_smurf_MultiClass.C+;
  export METHODS=BDTG;
fi

if [ ${NJETS} == "hzz" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_hzz.C+;
  export TAG=ntuples_hzz_${MH}train
elif [ ${NJETS} == "wh3l" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_wh3l.C+;
  export TAG=ntuples_wh3l_${MH}train
elif [ ${NJETS} == "3" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_MultiClass.C+;
  export NJETS=2;
  export TAG=ntuples2012_PostICHEP_${MH}train_${NJETS}jets;
fi

export DO_TRAINING=0;
if [ ${DO_TRAINING} == "1" ]; then
  export SIG_TRAIN=data/hww${MH}.root;
  export BKG_TRAIN=data/training.root;
  if [ ${NJETS} == "hzz" ]; then
    export NJETS=999;
    export SIG_TRAIN=data/hzz${MH}.root;
    export BKG_TRAIN=data/zz.root;
  elif [ ${NJETS} == "wh3l" ]; then
    export NJETS=999;
    export SIG_TRAIN=data/hww_training.root;
    export BKG_TRAIN=data/backgroundD_3l.root;
 fi
  mkdir -p weights;
  root -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
  exit;
fi

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
### samples must be in "data" folder
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
output/hww110.root
output/hww115.root
output/hww120.root
output/hww125.root
output/hww130.root
output/hww135.root
output/hww140.root
output/hww145.root
output/hww150.root
output/hww155.root
output/hww160.root
output/hww170.root
output/hww180.root
output/hww190.root
output/hww200.root
output/hww250.root
output/hww300.root
output/hww350.root
output/hww400.root
output/hww450.root
output/hww500.root
output/hww550.root
output/hww600.root
data/hww_syst_skim6.root
data/data_skim6.root
data/backgroundA_skim6.root
EOF

export evaluateMVAFile=evaluateMVA_smurf_hww.C+;
if [ ${NJETS} == "hzz" ]; then
  export NJETS=999;
  export evaluateMVAFile=evaluateMVA_smurf_hzz.C+;
elif [ ${NJETS} == "wh3l" ]; then
  export METHODS=BDT,BDTG,Likelihood;
  export TAG=ntuples_wh3l_120train;
  export NJETS=999;
  export evaluateMVAFile=evaluateMVA_smurf_wh3l.C+;
elif [ ${NJETS} == "wh2l" ]; then
  export METHODS=BDTG;
  export TAG=ntuples_wh2l_125train;
  export NJETS=999;
  export evaluateMVAFile=evaluateMVA_smurf_hww.C+;
fi

for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  if [ ${WEIGHTSONLY} == "1" ]; then
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"\",\"\",\"\",${NJETS},1,0,\"/data\",${doMultiClass},3\);
  else
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"/data\",${doMultiClass},3\);
  fi
  
done
rm -f list_samples.txt;
