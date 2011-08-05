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
  #export BKG_TRAIN=data/ww2l_pythia.root;
  export BKG_TRAIN=data/background42x.root;
  mkdir -p weights;
  root -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi
#exit;

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
### samples must be in "data" folder
#data/data_2l.root
#data/hww${MH}.root
#data/background42x_spring11dy.root
#data/background42x.root
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data/background42x.root
data/background42x_spring11dy.root
data/background42x_wwpythia_met20_zveto.root
data/data_2l.root
data/data-met20-1092ipb.root
data/dataToy_2l_met20_zveto_1.root
data/dataToy_2l_met20_zveto_2.root
data/dataToy_2l_met20_zveto_3.root
data/dyee.root
data/dymm.root
data/dytt.root
data/ggww.root
data/qqww.root
data/stop.root
data/ttbar.root
data/ttop.root
data/tw.root
data/wgamma.root
data/wjets.root
data/ww2l_mcnlo.root
data/wz.root
data/zz.root
data/hww${MH}.root
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
