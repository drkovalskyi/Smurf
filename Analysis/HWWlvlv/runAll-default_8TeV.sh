#!/bin/csh

setenv NJETS      $1;
setenv MH         $2;
setenv MAKEINPUTS $3;
setenv CAT         1;

setenv TAG       ntuples2012_${MH}train_${NJETS}jets;

if ($NJETS == 2) then
  #setenv TAG       ntuples2012_125train_${NJETS}jets;
  setenv TAG       ntuples2012_${MH}train_1jets;
endif

### data_Summer12
setenv SIG_TEST data2012/${TAG}_hww${MH}.root
#setenv SIG_TEST data2012/${TAG}_hww125.root
setenv BKG_TEST data2012/${TAG}_backgroundA_skim2.root
setenv DAT_TEST data2012/${TAG}_data_skim2.root;
#setenv DAT_TEST data2012/${TAG}_hww125.root;
setenv SYS_TEST data2012/${TAG}_hww_syst_skim6.root;

### Perform analysis
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",0,1,0,\"$SYS_TEST\",$CAT\); --> mm
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",1,1,0,\"$SYS_TEST\",$CAT\); --> me
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",2,1,0,\"$SYS_TEST\",$CAT\); --> em
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",3,1,0,\"$SYS_TEST\",$CAT\); --> ee
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,0,\"$SYS_TEST\",$CAT\); --> ll
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,0,\"$SYS_TEST\",$CAT\); --> sf
 #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,0,\"$SYS_TEST\",$CAT\); --> of

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
   root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_0j.input_8TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_shape${STRCATA}_8TeV.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_0j_shape_8TeV.txt
   mv hwwof${STRCATB}_0j.input_8TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_shape${STRCATA}_8TeV.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_0j_shape_8TeV.txt

   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_0j_cut_8TeV.txt
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_0j_cut_8TeV.txt

  else if ($NJETS == 1) then
   root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_1j.input_8TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_shape${STRCATA}_8TeV.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_1j_shape_8TeV.txt
   mv hwwof${STRCATB}_1j.input_8TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_shape${STRCATA}_8TeV.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_1j_shape_8TeV.txt

   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_1j_cut_8TeV.txt
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_1j_cut_8TeV.txt

  else if ($NJETS == 2) then
   root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_2j.input_8TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_2j_chan5_mh${MH}_shape${STRCATA}_8TeV.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_2j_shape_8TeV.txt
   mv hwwof${STRCATB}_2j.input_8TeV.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_2j_chan6_mh${MH}_shape${STRCATA}_8TeV.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_2j_shape_8TeV.txt
   mv output/histo_limits_${TAG}_2j_chan5_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_2j_cut_8TeV.txt
   mv output/histo_limits_${TAG}_2j_chan6_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_2j_cut_8TeV.txt

   #root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,\"$SYS_TEST\",$CAT\);
   #mv output/histo_limits_${TAG}_2j_chan4_mh${MH}_cut${STRCATA}_8TeV.txt   ${INPUTDIR}/${MH}/hww${STRCATB}_2j_cut_8TeV.txt

  endif

else
  if ($MAKEINPUTS == "") then
    setenv MAKEINPUTS 4;
  endif
  root -l -q -b PlotHiggsRes2012.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",$MAKEINPUTS,1,\"$SYS_TEST\",$CAT\);

endif
