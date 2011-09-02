#!/bin/csh

# ntuples should be on "data" folder

setenv NJETS      $1;
setenv MH         $2;
setenv MAKEINPUTS $3;

setenv TAG       ntuples_${MH}train_${NJETS}jets;
if ($NJETS == 2) then
  setenv TAG       ntuples_${MH}train_1jets;
endif

setenv SIG_TEST data/${TAG}_hww${MH}.root
setenv BKG_TEST data/${TAG}_backgroundA_skim1.root
setenv DAT_TEST data/${TAG}_data_2l_skim1.root;
#setenv SYS_TEST data/${TAG}_hww_syst_skim2.root;
setenv SYS_TEST "";

### Perform analysis
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",0,1,0,\"$SYS_TEST\"\); --> mm
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",1,1,0,\"$SYS_TEST\"\); --> me
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",2,1,0,\"$SYS_TEST\"\); --> em
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",3,1,0,\"$SYS_TEST\"\); --> ee
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,0,\"$SYS_TEST\"\); --> ll
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,0,\"$SYS_TEST\"\); --> sf
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,0,\"$SYS_TEST\"\); --> of

if ($MAKEINPUTS == "") then
  setenv INPUTDIR /data/smurf/ceballos/inputLimits/ana_v6_test;
  mkdir -p ${INPUTDIR}/${MH}
  if ($NJETS == 0) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\"\);
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\"\);
   mv hwwsf_0j.input.root ${INPUTDIR}/${MH}/hwwsf_0j.input.root
   mv output/histo_limits_ntuples_${MH}train_0jets_0j_chan5_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwsf_0j_shape.txt
   mv hwwof_0j.input.root ${INPUTDIR}/${MH}/hwwof_0j.input.root
   mv output/histo_limits_ntuples_${MH}train_0jets_0j_chan6_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwof_0j_shape.txt

   mv output/histo_limits_ntuples_${MH}train_0jets_0j_chan5_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwsf_0j_cut.txt
   mv output/histo_limits_ntuples_${MH}train_0jets_0j_chan6_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwof_0j_cut.txt

  else if ($NJETS == 1) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\"\);
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\"\);
   mv hwwsf_1j.input.root ${INPUTDIR}/${MH}/hwwsf_1j.input.root
   mv output/histo_limits_ntuples_${MH}train_1jets_1j_chan5_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwsf_1j_shape.txt
   mv hwwof_1j.input.root ${INPUTDIR}/${MH}/hwwof_1j.input.root
   mv output/histo_limits_ntuples_${MH}train_1jets_1j_chan6_mh${MH}_shape.txt ${INPUTDIR}/${MH}/hwwof_1j_shape.txt

   mv output/histo_limits_ntuples_${MH}train_1jets_1j_chan5_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwsf_1j_cut.txt
   mv output/histo_limits_ntuples_${MH}train_1jets_1j_chan6_mh${MH}_cut.txt   ${INPUTDIR}/${MH}/hwwof_1j_cut.txt

  else if ($NJETS == 2) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,\"$SYS_TEST\"\);
   mv output/histo_limits_ntuples_${MH}train_1jets_2j_chan4_mh${MH}_cut.txt ${INPUTDIR}/${MH}/hww_2j_cut.txt

  endif

else
  if ($MAKEINPUTS == "") then
    setenv MAKEINPUTS 4;
  endif
  root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",$MAKEINPUTS,1,\"$SYS_TEST\"\);

endif
