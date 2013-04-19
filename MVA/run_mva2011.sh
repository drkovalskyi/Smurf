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
data2011/hww${MH}.root
data2011/h${MH}tt-vtth.root
data2011/wh3l${MH}.root
data2011/xww0m${MH}.root
data2011/xww0p${MH}.root
data2011/xww2p${MH}.root
data2011/xww2pqq${MH}.root
data2011/zhww${MH}.root
data2011/zh105inv.root
data2011/zh115inv.root
data2011/zh125inv.root
data2011/zh135inv.root
data2011/zh145inv.root
data2011/backgroundA.root
data2011/data_mit_2fake.root
data2011/data.root
data2011/dyee.root
data2011/dymm.root
data2011/dytt.root
data2011/ggww.root
data2011/hww_syst.root
data2011/qqww.root
data2011/training.root
data2011/ttbar_mg.root
data2011/ttbar.root
data2011/tw_ds.root
data2011/tw.root
data2011/wg3l.root
data2011/wgamma_lgamma.root
data2011/wgamma.root
data2011/wjets.root
data2011/ww_mcnlo_down.root
data2011/ww_mcnlo.root
data2011/ww_mcnlo_up.root
data2011/www.root
data2011/wz_py.root
data2011/wz.root
data2011/zz_mg.root
data2011/zz_py.root
EOF

export evaluateMVAFile=evaluateMVA_smurf_hww2011.C+;
if [ ${NJETS} == "hzz" ]; then
  export NJETS=999;
  export evaluateMVAFile=evaluateMVA_smurf_hzz.C+;
elif [ ${NJETS} == "wh3l" ]; then
  export METHODS=BDT,BDTG,Likelihood;
  export TAG=ntuples_wh3l_120train;
  export NJETS=999;
  export evaluateMVAFile=evaluateMVA_smurf_wh3l.C+;
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
