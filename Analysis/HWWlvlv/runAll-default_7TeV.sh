#!/bin/csh

setenv NJETS      $1;
setenv MH         $2;
setenv MAKEINPUTS $3;
setenv CAT         4;

setenv TAG       ntuples_${MH}train_${NJETS}jets;
if ($MH == 0) then
  setenv TAG       ntuples_160train_${NJETS}jets;
endif

if ($NJETS == 2) then
  setenv TAG       ntuples_${MH}train_1jets;
  if ($MH == 0) then
    setenv TAG       ntuples_160train_1jets;
  endif
endif

### data_Summer11
setenv SIG_TEST data2011/${TAG}_hww${MH}.root
setenv BKG_TEST data2011/${TAG}_backgroundA_skim6.root
setenv DAT_TEST data2011/${TAG}_data_skim6.root;
setenv SYS_TEST data2011/${TAG}_hww_syst_skim6.root;

### Perform analysis
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",0,1,0,\"$SYS_TEST\",$CAT\); --> mm
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",1,1,0,\"$SYS_TEST\",$CAT\); --> me
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",2,1,0,\"$SYS_TEST\",$CAT\); --> em
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",3,1,0,\"$SYS_TEST\",$CAT\); --> ee
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,0,\"$SYS_TEST\",$CAT\); --> ll
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,0,\"$SYS_TEST\",$CAT\); --> sf
 #root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,0,\"$SYS_TEST\",$CAT\); --> of

setenv STRCATA ""
setenv STRCATB ""
if ($CAT == "13") then
  setenv STRCATA "_lt"
  setenv STRCATB  "lt"
endif
if ($MAKEINPUTS == "") then
  setenv INPUTDIR /data/smurf/ceballos/inputLimits/ana_test
  mkdir -p ${INPUTDIR}/${MH}
  if ($NJETS == 0) then
   root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_0j.input_7TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_shape${STRCATA}_7TeV.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_0j_shape_7TeV.txt
   mv hwwof${STRCATB}_0j.input_7TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_shape${STRCATA}_7TeV.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_0j_shape_7TeV.txt

   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_cut${STRCATA}_7TeV.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_0j_cut_7TeV.txt
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_cut${STRCATA}_7TeV.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_0j_cut_7TeV.txt

  else if ($NJETS == 1) then
   root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_1j.input_7TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_shape${STRCATA}_7TeV.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_1j_shape_7TeV.txt
   mv hwwof${STRCATB}_1j.input_7TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_shape${STRCATA}_7TeV.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_1j_shape_7TeV.txt

   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_cut${STRCATA}_7TeV.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_1j_cut_7TeV.txt
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_cut${STRCATA}_7TeV.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_1j_cut_7TeV.txt

  else if ($NJETS == 2) then
   root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_2j.input_7TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_2j_chan5_mh${MH}_shape${STRCATA}_7TeV.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_2j_shape_7TeV.txt
   mv hwwof${STRCATB}_2j.input_7TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_2j_chan6_mh${MH}_shape${STRCATA}_7TeV.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_2j_shape_7TeV.txt
   mv output/histo_limits_${TAG}_2j_chan5_mh${MH}_cut${STRCATA}_7TeV.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_2j_cut_7TeV.txt
   mv output/histo_limits_${TAG}_2j_chan6_mh${MH}_cut${STRCATA}_7TeV.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_2j_cut_7TeV.txt

  endif

else
  if ($MAKEINPUTS == "") then
    setenv MAKEINPUTS 4;
  endif
  root -l -q -b PlotHiggsRes2011.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",$MAKEINPUTS,1,\"$SYS_TEST\",$CAT\);

endif
