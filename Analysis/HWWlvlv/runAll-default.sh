#!/bin/csh

# ntuples should be on "data" folder

setenv NJETS      $1;
setenv MH         $2;
setenv MAKEINPUTS $3;
setenv CAT         2;

#setenv TAG       winterWW_ntuples_${MH}train_${NJETS}jets_winterWjets_ntuples_${MH}train_${NJETS}jets;
setenv TAG       ntuples_${MH}train_${NJETS}jets;
if ($MH == 0) then
  setenv TAG       ntuples_110train_${NJETS}jets;
endif

if ($NJETS == 2) then
  setenv TAG       ntuples_${MH}train_1jets;
  if ($MH == 0) then
    setenv TAG       ntuples_115train_1jets;
  endif
endif

### data_Summer11
setenv SIG_TEST data/${TAG}_hww${MH}.root
setenv BKG_TEST data/${TAG}_backgroundC_skim2.root
setenv DAT_TEST data/${TAG}_data_2l_skim2.root;
#setenv DAT_TEST data/${TAG}_hww124.root
setenv SYS_TEST data/${TAG}_hww_syst_skim3.root;
### data_Fall11
#setenv SIG_TEST data/${TAG}_hww${MH}.root
#setenv BKG_TEST data/${TAG}_backgroundD_skim2.root
#setenv DAT_TEST data/${TAG}_data_skim2.root;
#setenv SYS_TEST data/${TAG}_hww_syst_skim3.root;

### Perform analysis
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",0,1,0,\"$SYS_TEST\",$CAT\); --> mm
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",1,1,0,\"$SYS_TEST\",$CAT\); --> me
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",2,1,0,\"$SYS_TEST\",$CAT\); --> em
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",3,1,0,\"$SYS_TEST\",$CAT\); --> ee
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,0,\"$SYS_TEST\",$CAT\); --> ll
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,0,\"$SYS_TEST\",$CAT\); --> sf
 #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,0,\"$SYS_TEST\",$CAT\); --> of

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
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_0j.input.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_shape${STRCATA}.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_0j_shape.txt
   mv hwwof${STRCATB}_0j.input.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_shape${STRCATA}.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_0j_shape.txt

   mv output/histo_limits_${TAG}_0j_chan5_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_0j_cut.txt
   mv output/histo_limits_${TAG}_0j_chan6_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_0j_cut.txt

  else if ($NJETS == 1) then
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   mv hwwsf${STRCATB}_1j.input.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_shape${STRCATA}.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_1j_shape.txt
   mv hwwof${STRCATB}_1j.input.root ${INPUTDIR}/${MH}/.
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_shape${STRCATA}.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_1j_shape.txt

   mv output/histo_limits_${TAG}_1j_chan5_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_1j_cut.txt
   mv output/histo_limits_${TAG}_1j_chan6_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_1j_cut.txt

  else if ($NJETS == 2) then
   #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",5,1,\"$SYS_TEST\",$CAT\);
   #root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",6,1,\"$SYS_TEST\",$CAT\);
   #mv hwwsf${STRCATB}_2j.input.root ${INPUTDIR}/${MH}/.
   #mv output/histo_limits_${TAG}_2j_chan5_mh${MH}_shape${STRCATA}.txt ${INPUTDIR}/${MH}/hwwsf${STRCATB}_2j_shape.txt
   #mv hwwof${STRCATB}_2j.input.root ${INPUTDIR}/${MH}/.
   #mv output/histo_limits_${TAG}_2j_chan6_mh${MH}_shape${STRCATA}.txt ${INPUTDIR}/${MH}/hwwof${STRCATB}_2j_shape.txt
   mv output/histo_limits_${TAG}_2j_chan5_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hwwsf${STRCATB}_2j_cut.txt
   mv output/histo_limits_${TAG}_2j_chan6_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hwwof${STRCATB}_2j_cut.txt

   root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",4,1,\"$SYS_TEST\",$CAT\);
   mv output/histo_limits_${TAG}_2j_chan4_mh${MH}_cut${STRCATA}.txt   ${INPUTDIR}/${MH}/hww${STRCATB}_2j_cut.txt

  endif

else
  if ($MAKEINPUTS == "") then
    setenv MAKEINPUTS 4;
  endif
  root -l -q -b PlotHiggsRes.C+\($NJETS,$MH,\"$TAG\",\"$SIG_TEST\",\"$BKG_TEST\"\,\"$DAT_TEST\",$MAKEINPUTS,1,\"$SYS_TEST\",$CAT\);

endif
