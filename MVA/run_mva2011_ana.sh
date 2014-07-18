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

#export TAG=TEST_${MH}train_${NJETS}jets;
#export METHODS=BDTG;
export TAG=ntuples_${MH}train_${NJETS}jets;
#export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;
#export TAG=TEST_ntuples_${MH}train_${NJETS}jets;
#export METHODS=BDT,Likelihood,BDTG,BDTD,BDTB,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,BoostedFisher,LikelihoodD,LikelihoodPCA,FDA_GA,RuleFit;
export METHODS=BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${NJETS} == "hzz" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_hzz.C+;
  export TAG=ntuples_hzz_${MH}train
elif [ ${NJETS} == "wh3l" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_wh3l.C+;
  export TAG=ntuples_wh3l_${MH}train
elif [ ${NJETS} == "wh2l" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_wh2l.C+;
  export TAG=ntuples_wh2l_${MH}train
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
    export BKG_TRAIN=data/backgroundD_3l.root
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
data2011/xww125p6_g125ww4l-2bp.root
data2011/xww125p6_g125ww4l-2mh10.root
data2011/xww125p6_g125ww4l-2mh9.root
data2011/xww125p6_g125ww4l-2ph3.root
data2011/xww125p6_g125ww4l-2ph6.root
data2011/xww125p6_gqq125ww4l-2hp.root
data2011/xww125p6_gqq125ww4l-2mh9.root
data2011/xww125p6_gqq125ww4l-2ph7.root
data2011/xww125p6_h125ww4l-0m-wh.root
data2011/xww125p6_v125ww4l-1m.root
data2011/xww125p6_v125ww4l-1p.root
data2011/xww125p6_v125ww4l-mix.root
data2011/xww125p6_x125ww4l-0l1f05ph180.root
data2011/xww125p6_x125ww4l-0mf05ph0.root
data2011/xww125p6_x125ww4l-0mf05ph0-vbf.root
data2011/xww125p6_x125ww4l-0phf05ph0.root
data2011/xww125p6_x125ww4l-0ph.root
data2011/xww125p6_x125ww4l-0p-jjh.root
data2011/xww125p6_x125ww4l-0pm.root
data2011/hww110.root
data2011/hww115.root
data2011/hww118.root
data2011/hww120.root
data2011/hww122.root
data2011/hww124.root
data2011/hww126.root
data2011/hww128.root
data2011/hww130.root
data2011/hww135.root
data2011/hww140.root
data2011/hww150.root
data2011/hww160.root
data2011/hww170.root
data2011/hww180.root
data2011/hww190.root
data2011/hww200.root
data2011/hww250.root
data2011/hww300.root
data2011/hww350.root
data2011/hww400.root
data2011/hww450.root
data2011/hww500.root
data2011/hww550.root
data2011/hww600.root
data2011/xww0p125.root
data2011/xww0m125.root
data2011/xww2p125.root
data2011/xww2pqq125.root
data2011/hww_syst_skim6.root
data2011/data_skim6.root
data2011/backgroundA_skim6.root
EOF
#data2011/hww${MH}.root
#data2011/xww0m125.root
#data2011/xww0p125.root
#data2011/xww2p125.root
#data2011/xww2pqq125.root

export evaluateMVAFile=evaluateMVA_smurf_hww2011.C+;
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
  export evaluateMVAFile=evaluateMVA_smurf_hww2011.C+;
fi

for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  if [ ${WEIGHTSONLY} == "1" ]; then
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"\",\"\",\"\",${NJETS},1,0,\"/data\",0,4\);
  else
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"/data\",0,4\);
  fi
  
done
rm -f list_samples.txt;
