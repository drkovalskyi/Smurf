
#grep -e data: -e S: Smurf/Analysis/HWWlvlv/headers/log_ww_sel.txt| awk '{if(NR%2==1)printf("%d ",$2);else printf("%f %f\n",$2,$4);}'| awk '{if(NR%7==1)a="ll";else if(NR%7==2)a="SF";else if(NR%7==3)a="OF";else if(NR%7==4)a="mm";else if(NR%7==5)a="em";else if(NR%7==6)a="me";else if(NR%7==0)a="ee";if(NR<=7)b="0j";else if(NR<=14)b="1j";else if(NR<=21)b="2j";else if(NR<=28)b="0j";else if(NR<=35)b="1j";else if(NR<=42)b="2j";else if(NR<=49)b="0j";else if(NR<=56)b="1j";else if(NR<=63)b="2j";printf("chan(%2s-%2s) ==> data: %4d S: %7.2f B: %7.2f -> %5.3f +/- %5.3f\n",b,a,$1,$2,$3,($1-$3)/$2,sqrt($1)/$2);}'

#same-sign, 5 final states, data, bkg, qqww, ggww, top, wjets, vv, zjets, vg
#grep -e data -e "bdg(27)" -e "bdg(28)" -e "bdg(29)" -e "bdg(30)" -e "bdg(xW)" -e "bdg(xZ)" -e "bdg(tX)" -e "bdg(Vg)" ../../Smurf/Analysis/HWWlvlv/headers/log_ss_ww_20_20.txt|awk '{if($1!="data:")printf("%8.3f %8.3f ",$3,$5);else printf("%4d ",$2);if(NR%9==0)printf("\n");}' > kkk
#awk '{printf("    & %4d & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f\n",$1,($2+$4+$6+$8+$10+$12+$14+$16),sqrt($3*$3+$5*$5+$7*$7+$9*$9+$11*$11+$13*$13+$15*$15*$17*$17),$6,$7,$8,$9,$14,$15,$10,$11,($2+$4),sqrt($3*$3+$5*$5),$12,$13,$16,$17)}' kkk


export NSEL=$1;

if [ $1 == 1 ]; then
cd $CMSSW_BASE/src/Ana/nt_scripts/;

sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' optimalCuts_42x.C;
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,4,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,5,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,6,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,0,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,1,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,2,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,3,2)';

root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,14,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,15,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,16,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,10,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,11,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,12,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,13,2)';

root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,24,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,25,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,26,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,20,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,21,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,22,2)';
root -l -b -q optimalCuts_42x.C+'(29,"","",0,"ntuples_42x_v9/backgroundA_smik6.root","ntuples_42x_v9/hww160.root","ntuples_42x_v9/data_smik6.root",1,23,2)';

else

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor2011.C+'(4)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_10_met_7TeV.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met_7TeV.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_7TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors2011.C+'(4)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met_7TeV.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_7TeV.h;
cp TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopVBFBkgScaleFactors_20_10_met_7TeV.h;
mv TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopVBFBkgScaleFactors_7TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor2011.C+'(4)';
cp WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/WWBkgScaleFactors_7TeV.h;
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_7TeV.h;

sed -ie 's/bool forBDTAna = false/bool forBDTAna = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor2011.C;
root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor2011.C+'(4)';
cat $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_7TeV.h DYBkgScaleFactorsBDT.h > DYBkgScaleFactors.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_7TeV.h;
sed -ie 's/bool forBDTAna = true/bool forBDTAna = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor2011.C;
mv DYEstimateTableBDT.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTableBDT_20_10_met_7TeV.txt;
mv DYBkgScaleFactorsBDT.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactorsBDT_20_10_met_7TeV.h;

cat $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met_7TeV.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactorsBDT_20_10_met_7TeV.h > temp_7TeV.h;
mv temp_7TeV.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met_7TeV.h;

rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.Ce;
rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.Ce;

fi
