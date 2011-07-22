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
export LABEL=$3
export INCLUDE_PHOTONS=$4
export DO_TRAINING=$5
#
#export TAG=ntuples_${MH}train_${NJETS}jets;
export TAG=ntuples_hzz_${MH}train_${LABEL};
#export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;
export METHODS=MLPBNN,BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf_hzz.C+;
if [ ${DO_TRAINING} == "1" ]; then
  export SIG_TRAIN=data/hzz${MH}.root;
#  export BKG_TRAIN=data/background_vv.root;
  export BKG_TRAIN=data/zz.root;
  if [ ${INCLUDE_PHOTONS} == "1" ]; then
  export PHO_TRAIN=data/data_photons.goodlumi.root;
  fi
  mkdir -p weights;
  root -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${PHO_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi
#exit;

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
### samples must be in "data" folder
rm -f list_samples.txt;
#data/data_2l.root
#data/hww${MH}.root
#data/background_mconly.root
#data/background.root
#data/zz.root
#data/hzz250.root
#data/background_mconly.root
#data/background.root
#data/zz.root
cat > list_samples.txt <<EOF
data/background.root
data/hzz${MH}.root
EOF

export evaluateMVAFile=evaluateMVA_smurf_hzz.C+;
  
for i in `cat list_samples.txt` ; do
  dataset=${i%%,*};
  echo "filling MVA information in sample: "  $dataset
  root -l -q -b ${evaluateMVAFile}\(\"${dataset}\",${MH},\"${LABEL}\",\"${METHODS}\",\"${TAG}\"\);
done
rm -f list_samples.txt;

root -l -q -b evaluateMVA_smurf_hzz_photons.C+\(\"data/data_photons.goodlumi.root\",${MH},\"${LABEL}\",\"${METHODS}\",\"${TAG}\"\);




#My test samples

# root -l evaluateMVA_smurf.C+\(\"PowhegPythia_ggHToWW_mH130.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l evaluateMVA_smurf.C+\(\"PowhegPythia_ggHToWW_mH130_RScaleUp_FScaleUp.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l evaluateMVA_smurf.C+\(\"PowhegPythia_ggHToWW_mH130_RScaleDown_FScaleDown.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)

# root -l evaluateMVA_smurf.C+\(\"qqww.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l evaluateMVA_smurf.C+\(\"histo_p11-ww2l-mcnlo-v1g1-pu_all_noskim.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l evaluateMVA_smurf.C+\(\"histo_p11-ww2lup-mcnlo-v1g1-pu_all_noskim.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l evaluateMVA_smurf.C+\(\"histo_p11-ww2ldown-mcnlo-v1g1-pu_all_noskim.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)

# root -l -q evaluateMVA_smurf.C+\(\"HwwAcceptance_p11-vvj-v1g1-pu_1_noskim.smurf.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l -q evaluateMVA_smurf.C+\(\"HwwAcceptance_p11-ww2l-mcnlo-v1g1-pu_1_noskim_0000.root_copy.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l -q evaluateMVA_smurf.C+\(\"HwwAcceptance_p11-ww2lup-mcnlo-v1g1-pu_1_noskim_0000.root_copy.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l -q evaluateMVA_smurf.C+\(\"HwwAcceptance_p11-ww2ldown-mcnlo-v1g1-pu_1_noskim_0000.root_copy.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_0jets\"\)
# root -l -q evaluateMVA_smurf1.C+\(\"HwwAcceptance_p11-vvj-v1g1-pu_1_noskim.smurf.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_1jets\"\)
# root -l -q evaluateMVA_smurf1.C+\(\"HwwAcceptance_p11-ww2l-mcnlo-v1g1-pu_1_noskim_0000.root_copy.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_1jets\"\)
# root -l -q evaluateMVA_smurf1.C+\(\"HwwAcceptance_p11-ww2lup-mcnlo-v1g1-pu_1_noskim_0000.root_copy.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_1jets\"\)
# root -l -q evaluateMVA_smurf1.C+\(\"HwwAcceptance_p11-ww2ldown-mcnlo-v1g1-pu_1_noskim_0000.root_copy.root\",130,\"KNN,BDT,BDTD,MLPBNN,BDTG\",\"ntuples_130train_1jets\"\)
