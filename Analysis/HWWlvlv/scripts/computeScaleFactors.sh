root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(0)' & tee HWWAnalysis_DYBkgEstimate0.log
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(0,"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011A/backgroundC_skim1.root","/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011A/data_2l_skim1.root")' & tee HWWAnalysis_TopBkgEstimate0.log
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(0,"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011A/backgroundC_skim1.root","/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011A/data_2l_skim1.root")' & tee HWWAnalysis_WWBkgEstimate0.log
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors0.h 
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors0.h 
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors0.h 

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(1)' & tee HWWAnalysis_DYBkgEstimate1.log
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(1,"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011B/backgroundC_skim1.root","/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011B/data_2l_skim1.root")' & tee HWWAnalysis_TopBkgEstimate1.log
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(1,"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011B/backgroundC_skim1.root","/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Run2011B/data_2l_skim1.root")' & tee HWWAnalysis_WWBkgEstimate1.log
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors1.h 
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors1.h 
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors1.h 

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(2)' & tee HWWAnalysis_DYBkgEstimate2.log
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(2,"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/backgroundC_skim1.root","/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/data_2l_skim1.root")' & tee HWWAnalysis_TopBkgEstimate2.log
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(2,"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/backgroundC_skim1.root","/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/data_2l_skim1.root")' & tee HWWAnalysis_WWBkgEstimate2.log
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors2.h 
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors2.h 
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors2.h
