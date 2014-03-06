
#grep -e data: -e S: Smurf/Analysis/HWWlvlv/headers/log_ww_sel.txt| awk '{if(NR%2==1)printf("%d ",$2);else printf("%f %f\n",$2,$4);}'| awk '{if(NR%7==1)a="ll";else if(NR%7==2)a="SF";else if(NR%7==3)a="OF";else if(NR%7==4)a="mm";else if(NR%7==5)a="em";else if(NR%7==6)a="me";else if(NR%7==0)a="ee";if(NR<=7)b="0j";else if(NR<=14)b="1j";else if(NR<=21)b="2j";else if(NR<=28)b="0j";else if(NR<=35)b="1j";else if(NR<=42)b="2j";else if(NR<=49)b="0j";else if(NR<=56)b="1j";else if(NR<=63)b="2j";printf("chan(%2s-%2s) ==> data: %4d S: %7.2f B: %7.2f -> %5.3f +/- %5.3f\n",b,a,$1,$2,$3,($1-$3)/$2,sqrt($1)/$2);}'

#same-sign, 5 final states, data, bkg, qqww, ggww, top, wjets, vv, zjets, vg
#grep -e data -e "bdg(27)" -e "bdg(28)" -e "bdg(29)" -e "bdg(30)" -e "bdg(xW)" -e "bdg(xZ)" -e "bdg(tX)" -e "bdg(19)" -e "bdg(20)" -e "bdg(21)" log|grep -v root|awk '{if($1!="data:")printf("%8.3f %8.3f ",$3,$5);else printf("%4d ",$2);if(NR%11==0)printf("\n");}' > kkk
#awk '{printf("    & %4d & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f\n",$1,($2+$4+$6+$8+$10+$16+$18+$20),sqrt($3*$3+$5*$5+$7*$7+$9*$9+$11*$11+$17*$17+$19*$19+$21*$21),sqrt($2*$2*0.30*0.30+$4*$4*0.30*0.30+$6*$6*0.50*0.50+$8*$8*0.10*0.10+$10*$10+$16*$16*0.36*0.36+$18*$18*0.30*0.30+$20*$20*0.10*0.10),$12+$14,sqrt($13*$13+$15*$15),$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$16,$17,$18,$19,$20,$21)}' kkk

#grep -A1 data log|grep pm > aaa
#awk '{printf("%7.2f",$3);if(NR%4!=0)printf(",");else printf("\n");}' aaa

source ~/EVAL_SH65 5_3_14;

export NSEL=$1;

if [ $1 == 1 ]; then
cd $CMSSW_BASE/src/Ana/nt_scripts/;

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_20_met.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_20_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' optimalCuts_53x.C;
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,5,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,6,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,0,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,1,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,2,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,3,3)';

root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,15,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,16,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,10,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,11,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,12,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,13,3)';

root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,24,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,25,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,26,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,20,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,21,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,22,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,23,3)';

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' optimalCuts_53x.C;
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,5,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,6,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,0,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,1,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,2,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,3,3)';

root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,15,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,16,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,10,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,11,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,12,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,13,3)';

root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,24,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,25,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,26,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,20,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,21,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,22,3)';
root -l -b -q optimalCuts_53x.C+'(29,0,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,23,3)';

rm optimalCuts_53x.Ce;

elif [ $1 == 2 ]; then
echo "making nice plots";
cd $CMSSW_BASE/src/Ana/nt_scripts/;

root -l -b -q optimalCuts_53x.C+'(29,1,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_ptmax.root;
root -l -b -q optimalCuts_53x.C+'(29,2,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_ptmin.root;
root -l -b -q optimalCuts_53x.C+'(29,7,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_massll.root;
root -l -b -q optimalCuts_53x.C+'(29,8,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_mt.root;
root -l -b -q optimalCuts_53x.C+'(29,20,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_deltaphill.root;
root -l -b -q optimalCuts_53x.C+'(29,57,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,4,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_deltarll.root;

root -l -b -q optimalCuts_53x.C+'(29,1,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_ptmax.root;
root -l -b -q optimalCuts_53x.C+'(29,2,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_ptmin.root;
root -l -b -q optimalCuts_53x.C+'(29,7,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_massll.root;
root -l -b -q optimalCuts_53x.C+'(29,8,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_mt.root;
root -l -b -q optimalCuts_53x.C+'(29,20,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_deltaphill.root;
root -l -b -q optimalCuts_53x.C+'(29,57,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_deltarll.root;
root -l -b -q optimalCuts_53x.C+'(29,55,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_dphidilepmet.root;
root -l -b -q optimalCuts_53x.C+'(29,53,"ntuples_53x/backgroundA_skim6.root","ntuples_53x/hww125.root","ntuples_53x/data_skim6.root",1,14,3)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_dphidilepjet.root;

else

sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(3)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_20_met.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_20_met.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(3)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_20_met.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;

sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(3)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_10_met.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(3)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
cp TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopVBFBkgScaleFactors_20_10_met.h;
mv TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopVBFBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(3)';
cp WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/WWBkgScaleFactors.h;
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h;

sed -ie 's/bool forBDTAna = false/bool forBDTAna = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(3)';
cat $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h DYBkgScaleFactorsBDT.h > DYBkgScaleFactors.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
sed -ie 's/bool forBDTAna = true/bool forBDTAna = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
mv DYEstimateTableBDT.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTableBDT_20_10_met.txt;
mv DYBkgScaleFactorsBDT.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactorsBDT_20_10_met.h;

cat $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactorsBDT_20_10_met.h > temp.h;
mv temp.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h;

rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.Ce;
rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.Ce;

fi
