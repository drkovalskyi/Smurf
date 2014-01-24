#!/bin/bash

source ~/EVAL_SH64 4_2_8_patch4

cd $CMSSW_BASE/src;

echo THEDATA: $THEDATA;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(2,"/data/smurf/data/Run2011_Summer11_SmurfV7_42X/mitf-alljets/backgroundC_skim1.root","$THEDATA")';
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/;
root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(2,"/data/smurf/data/Run2011_Summer11_SmurfV7_42X/mitf-alljets/backgroundC_skim1.root","$THEDATA")';
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/;
