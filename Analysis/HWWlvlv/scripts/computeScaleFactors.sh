
#grep -e data: -e S: Smurf/Analysis/HWWlvlv/headers/log_ww_sel.txt| awk '{if(NR%2==1)printf("%d ",$2);else printf("%f %f\n",$2,$4);}'| awk '{if(NR%7==1)a="ll";else if(NR%7==2)a="SF";else if(NR%7==3)a="OF";else if(NR%7==4)a="mm";else if(NR%7==5)a="em";else if(NR%7==6)a="me";else if(NR%7==0)a="ee";if(NR<=7)b="0j";else if(NR<=14)b="1j";else if(NR<=21)b="2j";else if(NR<=28)b="0j";else if(NR<=35)b="1j";else if(NR<=42)b="2j";printf("chan(%2s-%2s) ==> data: %4d S: %7.2f B: %7.2f -> %5.3f +/- %5.3f\n",b,a,$1,$2,$3,($1-$3)/$2,sqrt($1)/$2);}'

export NSEL=$1;

if [ $1 == 1 ]; then
cd $CMSSW_BASE/src/Ana/nt_scripts/;

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_20_met.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_20_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' optimalCuts_52x.C;
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,4,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,5,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,6,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,0,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,1,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,2,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,3,0)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,14,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,15,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,16,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,10,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,11,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,12,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,13,0)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,24,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,25,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,26,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,20,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,21,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,22,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,23,0)';

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' optimalCuts_52x.C;
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,4,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,5,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,6,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,0,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,1,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,2,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,3,0)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,14,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,15,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,16,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,10,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,11,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,12,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,13,0)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,24,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,25,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,26,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,20,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,21,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,22,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,23,0)';

sed -ie 's/) useDYMVA = false/) useDYMVA = true/' optimalCuts_52x.C;
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,4,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,5,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,6,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,0,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,1,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,2,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,3,0)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,14,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,15,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,16,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,10,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,11,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,12,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,13,0)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,24,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,25,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,26,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,20,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,21,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,22,0)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundA_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",0,23,0)';
sed -ie 's/) useDYMVA = true/) useDYMVA = false/' optimalCuts_52x.C;

rm optimalCuts_52x.Ce;

else

sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(0)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_20_met.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_20_met.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(0)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_20_met.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;

sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(0)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_10_met.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(0)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
cp TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopVBFBkgScaleFactors_20_10_met.h;
mv TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopVBFBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(0)';
cp WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/WWBkgScaleFactors.h;
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h;

sed -ie 's/bool forBDTAna = false/bool forBDTAna = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(0)';
cat $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h DYBkgScaleFactorsBDT.h > DYBkgScaleFactors.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
sed -ie 's/bool forBDTAna = true/bool forBDTAna = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
mv DYEstimateTableBDT.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTableBDT_20_10_met.txt;
mv DYBkgScaleFactorsBDT.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactorsBDT_20_10_met.h;

rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.Ce;
rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.Ce;

fi
