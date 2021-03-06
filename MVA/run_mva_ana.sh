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
elif [ ${NJETS} == "wh2l" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_wh2l.C+;
  export TAG=ntuples_wh2l_${MH}train
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
  elif [ ${NJETS} == "wh2l" ]; then
    export NJETS=999;
    export SIG_TRAIN=data/hww_training.root;
    export BKG_TRAIN=data/backgroundA_skim8.root;
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
data/xww125p6_s12-g125ww4l-2hp-v19.root
data/xww125p6_s12-gqq125ww4l-2hp-v19.root
data/xww125p6_s12-g125ww4l-2ph7-v19.root
data/xww125p6_s12-gqq125ww4l-2ph7-v19.root
data/xww125p6_s12-g125ww4l-2hm-v19.root
data/xww125p6_s12-gqq125ww4l-2hm-v19.root
data/xww125p6_s12-g125ww4l-2ph2-v19.root
data/xww125p6_s12-gqq125ww4l-2ph2-v19.root
data/xww125p6_s12-gqq125ww4l-2ph3-v19.root
data/xww125p6_s12-g125ww4l-2ph3-v19.root
data/xww125p6_s12-g125ww4l-2ph6-v19.root
data/xww125p6_s12-gqq125ww4l-2ph6-v19.root
data/xww125p6_s12-g125ww4l-2mh9-v19.root
data/xww125p6_s12-gqq125ww4l-2mh9-v19.root
data/xww125p6_s12-g125ww4l-2mh10-v19.root
data/xww125p6_s12-gqq125ww4l-2mh10-v19.root
data/xww125p6_s12-g125ww4l-2bp-v19.root
data/xww125p6_s12-gqq125ww4l-2bp-v19.root
data/xww125p6_s12-v125ww4l-mix-v19.root
data/xww125p6_s12-x125ww4l-0mt-v19.root
data/xww125p6_s12-x125ww4l-0phf05ph0-v19.root
data/xww125p6_s12-x125ww4l-0phf05-v19.root
data/xww125p6_s12-x125ww4l-0pm-v19.root
data/xww125p6_s12-x125ww4l-0l1-v19.root
data/xww125p6_s12-x125ww4l-0l1f05ph180-v19.root
data/xww125p6_s12-v125ww4l-1m-v19.root
data/xww125p6_s12-v125ww4l-1p-v19.root
data/xww125p6_s12-x125ww4l-0ph-v19.root
data/xww125p6_s12-x125ww4l-0mf05ph0-v19.root
data/xww125p6_s12-x125ww4l-0m-jjh-v19.root
data/xww125p6_s12-x125ww4l-0p-jjh-v19.root
data/xww125p6_s12-x125ww4l-0ph-vbf-v19.root
data/xww125p6_s12-x125ww4l-0mf05ph0-vbf-v19.root
data/xww125p6_s12-x125ww4l-0m-vbf-v19.root
data/xww125p6_s12-x125ww4l-0p-vbf-v19.root
data/xww125p6_s12-h125ww4l-0ph-wh-v19.root
data/xww125p6_s12-x125ww4l-0mf05ph0-wh-v19.root
data/xww125p6_s12-h125ww4l-0m-wh-v19.root
data/xww125p6_s12-h125ww4l-0p-wh-v19.root
data/xww125p6_s12-x125ww4l-0ph-zh-v19.root
data/xww125p6_s12-x125ww4l-0mf05ph0-zh-v19.root
data/xww125p6_s12-x125ww4l-0p-zh-v19.root
data/xww125p6_s12-x125ww4l-0m-zh-v19.root
data/hww110.root
data/hww115.root
data/hww120.root
data/hww125.root
data/hww130.root
data/hww135.root
data/hww140.root
data/hww145.root
data/hww150.root
data/hww155.root
data/hww160.root
data/hww170.root
data/hww180.root
data/hww190.root
data/hww200.root
data/hww250.root
data/hww300.root
data/hww350.root
data/hww400.root
data/hww450.root
data/hww500.root
data/hww550.root
data/hww600.root
data/xww1m125.root
data/xww1p125.root
data/xww1mix125.root
data/xww0p125.root
data/xww0m125.root
data/xww2p125.root
data/xww2pqq125.root
data/hww_syst_skim6.root
data/hww_syst_sherpa_skim6.root
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
