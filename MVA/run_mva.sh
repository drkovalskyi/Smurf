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
export doMultiClass=0;

#export TAG=TEST_${MH}train_${NJETS}jets;
#export METHODS=Likelihood,BDT,BDTD,BDTG;
export TAG=ntuples2012_PostICHEP_${MH}train_${NJETS}jets;
#export TAG=TEST_ntuples_${MH}train_${NJETS}jets;
export METHODS=Fisher,BDT,Likelihood,BDTG,BDTD,BDTB,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,BoostedFisher,LikelihoodD,LikelihoodPCA,RuleFit;
#export METHODS=BDTG;

### Training: change done hand made, it's an expert option
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${doMultiClass} == "1" ]; then
  export TAG=ntuples2012_MultiClass_${MH}train_${NJETS}jets;
  export trainMVA_smurfFile=trainMVA_smurf_MultiClass.C+;
  export METHODS=BDTG;
fi

if [ ${NJETS} == "zhinv" ]; then
  export trainMVA_smurfFile=trainMVA_smurf_zhinv.C+;
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
  if [ ${NJETS} == "zhinv" ]; then
    export NJETS=999;
    export SIG_TRAIN=output/zh125inv.root;
    export BKG_TRAIN=output/backgroundA_skim10.root;
  elif [ ${NJETS} == "wh3l" ]; then
    export NJETS=999;
    export SIG_TRAIN=data/hww_training.root;
    export BKG_TRAIN=data/backgroundD_3l.root;
  elif [ ${NJETS} == "wh2l" ]; then
    export NJETS=900;
    export SIG_TRAIN=data/hww_training.root;
    export BKG_TRAIN=data/backgroundA_skim8.root;
  elif [ ${NJETS} == "2" ]; then
    export NJETS=2;
    export SIG_TRAIN=data/hww_training.root;
    export BKG_TRAIN=data/training.root;
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
data/hww${MH}.root
data/vhtt${MH}.root
data/xww1m${MH}.root
data/xww1p${MH}.root
data/xww1mix${MH}.root
data/xww0m${MH}.root
data/xww0p${MH}.root
data/xww2p${MH}.root
data/xww2pqq${MH}.root
data/zhww${MH}.root
data/hzz4l${MH}.root
data/gghww${MH}_minlo.root
data/zh125_hgoldgoldgamma.root
data/dyllpt100.root
data/bbww_8tev.root
data/bbww_14tev.root
data/hpp_mh200.root
data/hpp_mh800.root
data/wwss_qed_2_qcd_99_matching.root
data/wwss_qed_2_qcd_99_no_matching.root
data/wwss_qed_4_qcd_99_lt012.root
data/wwss_qed_4_qcd_99_lm0123.root
data/wwss_qed_4_qcd_99_ls0ls1.root
data/wzgamma3l_ewk_sm.root
data/wz3l_ewk_sm.root
data/wz3l_ewk_lt012.root
data/wz3l_ewk_lm0123.root
data/wz3l_ewk_ls0ls1.root
data/wwss_qed_4_qcd_99_ls_lm_lt.root
data/wwss_qed_4_qcd_99_sm.root
data/wwss_qed_2_qcd_99_sm.root
data/wwss_qed_4_qcd_0_sm.root
data/wzgamma_qed_3_qcd_99_sm.root
data/wzgamma_qed_5_qcd_99_sm.root
data/wzgamma_qed_5_qcd_0_sm.root
data/wzgamma_qed_5_qcd_0_lt1_1p25.root
data/wwss_ph_wh.root
data/wwss_ph_noh.root
data/zh105inv.root
data/zh115inv.root
data/zh125inv.root
data/zh135inv.root
data/zh145inv.root
data/zh175inv.root
data/zh200inv.root
data/zh300inv.root
data/backgroundA.root
data/backgroundB.root
data/backgroundC.root
data/backgroundD.root
data/backgroundEWK_skim8.root
data/backgroundEWK_skim14.root
data/hww_syst.root
data/hww_syst_sherpa.root
data/data_2fake.root
data/data_llg.root
data/lgamma_3l.root
data/data.root
data/data_ztt.root
data/dyll.root
data/gamma50.root
data/gamma75.root
data/gamma90.root
data/ggww.root
data/hzz4l125.root
data/lgamma50.root
data/lgamma75.root
data/lgamma90.root
data/qqww_py.root
data/qqww_powheg.root
data/qqww.root
data/training.root
data/ttbar_sherpa.root
data/ttbar_powheg.root
data/ttbar_mg.root
data/ttbar_madspin.root
data/ttbar_mcnlo.root
data/ttbar_amcnlo.root
data/ttw_amcnlo.root
data/ttw_mg.root
data/tw.root
data/wgammafo.root
data/wgamma.root
data/wglee.root
data/wglll.root
data/wglmm.root
data/wjets.root
data/ww2j_mg_ewk.root
data/ww2j_mg_qcd.root
data/ww2j_ph.root
data/ww_dps.root
data/wwmcnlodown.root
data/wwmcnlo.root
data/wwmcnloup.root
data/wwss.root
data/www_amcnlo.root
data/www_mg.root
data/www.root
data/wz.root
data/wz_py.root
data/wzz_amcnlo.root
data/wzz_mg.root
data/zgammafo.root
data/zgamma.root
data/zz4l_powheg.root
data/zz.root
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
