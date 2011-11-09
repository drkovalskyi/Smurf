#!/bin/csh

# ntuples should be on "data" folder

setenv NJETS      $1;
setenv MH         $2;
setenv MAKEINPUTS $3;

#setenv TAG       winterWW_ntuples_${MH}train_${NJETS}jets_winterWjets_ntuples_${MH}train_${NJETS}jets;
setenv TAG       ntuples_${MH}train_${NJETS}jets;
if ($MH == 0) then
  setenv TAG       ntuples_115train_${NJETS}jets;
endif

if ($NJETS == 2) then
  setenv TAG       ntuples_${MH}train_1jets;
endif

#setenv SIG_TEST data_LP_42x/${TAG}_hww${MH}.root
#setenv BKG_TEST data_LP_42x/${TAG}_backgroundC_skim2.root
#setenv DAT_TEST data_LP_42x/${TAG}_data_2l_skim2.root;
#setenv SYS_TEST data_LP_42x/${TAG}_hww_syst_skim3.root;
setenv SIG_TEST data/${TAG}_hww${MH}.root
setenv BKG_TEST data/${TAG}_backgroundC_skim2.root
setenv DAT_TEST data/${TAG}_data_2l_skim2.root;
setenv SYS_TEST data/${TAG}_hww_syst_skim3.root;
#setenv SIG_TEST data_LP/${TAG}_hww${MH}.root
#setenv BKG_TEST data_LP/${TAG}_background42x_METVeto_ZVeto.root
#setenv DAT_TEST data_LP/${TAG}_data_2l_METVeto_ZVeto.root
#setenv SYS_TEST data_LP/${TAG}_hww_syst_skim3.root;

### Perform analysis
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",0,1,0,\"$SYS_TEST\"\); --> mm
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",1,1,0,\"$SYS_TEST\"\); --> me
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",2,1,0,\"$SYS_TEST\"\); --> em
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",3,1,0,\"$SYS_TEST\"\); --> ee
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,0,\"$SYS_TEST\"\); --> ll
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,0,\"$SYS_TEST\"\); --> sf
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,0,\"$SYS_TEST\"\); --> of

if ($MAKEINPUTS == "") then
  setenv INPUTDIR /data/smurf/ceballos/inputLimits/ana_v7_test
  mkdir -p ${INPUTDIR}/${MH}
  if ($NJETS == 0) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\"\);
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\"\);
   mv hwwsf_0j.input.root ${INPUTDIR}/${MH}/hwwsf_0j.input.root
   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwsf_0j_shape.txt
   mv hwwof_0j.input.root ${INPUTDIR}/${MH}/hwwof_0j.input.root
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwof_0j_shape.txt

   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwsf_0j_cut.txt
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwof_0j_cut.txt

  else if ($NJETS == 1) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\"\);
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\"\);
   mv hwwsf_1j.input.root ${INPUTDIR}/${MH}/hwwsf_1j.input.root
   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwsf_1j_shape.txt
   mv hwwof_1j.input.root ${INPUTDIR}/${MH}/hwwof_1j.input.root
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwof_1j_shape.txt

   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwsf_1j_cut.txt
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwof_1j_cut.txt

  else if ($NJETS == 2) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,\"$SYS_TEST\"\);
   mv output/histo_limits_${TAG}_2j_chan4_mh${MH}_cut.txt ${INPUTDIR}/${MH}/hww_2j_cut.txt

  endif

else
  if ($MAKEINPUTS == "") then
    setenv MAKEINPUTS 4;
  endif
  root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",$MAKEINPUTS,1,\"$SYS_TEST\"\);

endif
