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

#export TAG=ntuplesNoMZCut_${MH}train_${NJETS}jets;
#export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;
export TAG=ntuples_${MH}train_${NJETS}jets;
export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;
#export TAG=TEST_ntuples_${MH}train_${NJETS}jets;
#export METHODS=BDTD,Likelihood,LikelihoodD,MLPBNN,BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${NJETS} == "hzz" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_hzz.C+;
  export NJETS=999;
fi

export DO_TRAINING=0;
if [ ${DO_TRAINING} == "1" ]; then
  export SIG_TRAIN=data/hww${MH}.root;
  export BKG_TRAIN=data/training.root;
  if [ ${NJETS} == "hzz" ]; then
    export SIG_TRAIN=data/hzz${MH}.root;
    export BKG_TRAIN=data/zz.root;
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
data/backgroundA.root
data/backgroundC_3l.root
data/backgroundC.root
data/backgroundC_skim1.root
data/backgroundC_skim2.root
data/backgroundD_3l.root
data/backgroundD.root
data/backgroundD_skim1.root
data/data_2l_3l.root
data/data_2l.root
data/data_2l_skim1.root
data/data_2l_skim2.root
data/data-emb-tau123.root
data/data_lfake.root
data/data_mit_2l_3l.root
data/data_mit_2l.root
data/data_mit_2l_skim1.root
data/data_mit_2l_skim2.root
data/data_mit_lfake.root
data/dyee.root
data/dymm.root
data/dytt.root
data/ggww.root
data/hww_syst.root
data/hww_syst_skim3.root
data/mc_l_fake_pu.root
data/qqww_py.root
data/qqww.root
data/training.root
data/ttbar_mg.root
data/ttbar.root
data/tw_ds.root
data/tw.root
data/wgamma_lgamma.root
data/wgamma.root
data/wjets.root
data/ww_mcnlo_down.root
data/ww_mcnlo.root
data/ww_mcnlo_up.root
data/wz_py.root
data/wz.root
data/zz_mg.root
data/zz_py.root
data/hww${MH}.root
EOF

export evaluateMVAFile=evaluateMVA_smurf_hww.C+;
if [ ${NJETS} == "hzz" ]; then
  export evaluateMVAFile=evaluateMVA_smurf_hzz.C+;
fi

for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  if [ ${WEIGHTSONLY} == "1" ]; then
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"\",\"\",\"\",${NJETS},1,0,\"/data\",3\);
  else
    root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"/data\",3\);
  fi
  
done
rm -f list_samples.txt;
