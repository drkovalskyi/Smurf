root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+ & tee HWWAnalysis_DYBkgEstimate.log
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'("/data/smurf/sixie/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/backgroundC_skim1.root","/data/smurf/sixie/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/data_2l_skim1.root")' & tee HWWAnalysis_TopBkgEstimate.log
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'("/data/smurf/sixie/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/backgroundC_skim1.root","/data/smurf/sixie/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/data_2l_skim1.root")' & tee HWWAnalysis_WWBkgEstimate.log
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/
