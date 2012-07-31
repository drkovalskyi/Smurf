
#grep -e data: -e S: Smurf/Analysis/HWWlvlv/headers/log_ww_sel.txt| awk '{if(NR%2==1)printf("%d ",$2);else printf("%f %f\n",$2,$4);}'| awk '{if(NR%7==1)a="ll";else if(NR%7==2)a="SF";else if(NR%7==3)a="OF";else if(NR%7==4)a="mm";else if(NR%7==5)a="em";else if(NR%7==6)a="me";else if(NR%7==0)a="ee";if(NR<=7)b="0j";else if(NR<=14)b="1j";else if(NR<=21)b="2j";else if(NR<=28)b="0j";else if(NR<=35)b="1j";else if(NR<=42)b="2j";else if(NR<=49)b="0j";else if(NR<=56)b="1j";else if(NR<=63)b="2j";printf("chan(%2s-%2s) ==> data: %4d S: %7.2f B: %7.2f -> %5.3f +/- %5.3f\n",b,a,$1,$2,$3,($1-$3)/$2,sqrt($1)/$2);}'

#same-sign, 5 final states, data, bkg, qqww, ggww, top, wjets, vv, zjets, vg
#grep -e data -e "bdg(27)" -e "bdg(28)" -e "bdg(29)" -e "bdg(30)" -e "bdg(xW)" -e "bdg(xZ)" -e "bdg(tX)" -e "bdg(Vg)" ../../Smurf/Analysis/HWWlvlv/headers/log_ss_ww_20_20.txt|awk '{if($1!="data:")printf("%8.3f %8.3f ",$3,$5);else printf("%4d ",$2);if(NR%9==0)printf("\n");}' > kkk
#awk '{printf("    & %4d & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f & %6.1f $\\pm$ %5.1f\n",$1,($2+$4+$6+$8+$10+$12+$14+$16),sqrt($3*$3+$5*$5+$7*$7+$9*$9+$11*$11+$13*$13+$15*$15*$17*$17),$6,$7,$8,$9,$14,$15,$10,$11,($2+$4),sqrt($3*$3+$5*$5),$12,$13,$16,$17)}' kkk


export NSEL=$1;

if [ $1 == 1 ]; then
cd $CMSSW_BASE/src/Ana/nt_scripts/;

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_20_met.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_20_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' optimalCuts_52x.C;
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,4,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,5,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,6,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,0,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,1,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,2,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,3,2)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,14,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,15,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,16,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,10,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,11,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,12,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,13,2)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,24,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,25,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,26,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,20,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,21,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,22,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,23,2)';

cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h  $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
cp $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' optimalCuts_52x.C;
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,4,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,5,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,6,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,0,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,1,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,2,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,3,2)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,14,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,15,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,16,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,10,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,11,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,12,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,13,2)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,24,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,25,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,26,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,20,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,21,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,22,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,23,2)';

sed -ie 's/) useDYMVA = false/) useDYMVA = true/' optimalCuts_52x.C;
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,4,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,5,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,6,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,0,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,1,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,2,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,3,2)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,14,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,15,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,16,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,10,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,11,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,12,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,13,2)';

root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,24,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,25,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,26,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,20,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,21,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,22,2)';
root -l -b -q optimalCuts_52x.C+'(29,"","",0,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww160.root","ntuples_52x/data_skim2.root",1,23,2)';
sed -ie 's/) useDYMVA = true/) useDYMVA = false/' optimalCuts_52x.C;

rm optimalCuts_52x.Ce;

elif [ $1 == 2 ]; then
echo "making nice plots";
cd $CMSSW_BASE/src/Ana/nt_scripts/;

root -l -b -q optimalCuts_52x.C+'(29,"","",1,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,4,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_ptmax.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",2,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,4,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_ptmin.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",7,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,4,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_massll.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",8,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,4,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_mt.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",20,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,4,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_deltaphill.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",57,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,4,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_0j_mh125_deltarll.root;

root -l -b -q optimalCuts_52x.C+'(29,"","",1,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_ptmax.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",2,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_ptmin.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",7,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_massll.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",8,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_mt.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",20,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_deltaphill.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",57,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_deltarll.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",55,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_dphidilepmet.root;
root -l -b -q optimalCuts_52x.C+'(29,"","",53,"ntuples_52x/backgroundB_skim2.root","ntuples_52x/hww125.root","ntuples_52x/data_skim2.root",1,14,2)';
mv histo_nice.root /data/smurf/ceballos/distributions/note_hww8tev_ichep2012/dist/wwpresel_1j_mh125_dphidilepjet.root;

else

sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
sed -ie 's/bool WWXSSel = false/bool WWXSSel = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(2)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_20_met.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_20_met.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(2)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_20_met.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;

sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
sed -ie 's/bool WWXSSel = true/bool WWXSSel = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(2)';
mv DYEstimateTable.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTable_20_10_met.txt;
cp DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactors_20_10_met.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.C+'(2)';
cp TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopBkgScaleFactors_20_10_met.h;
mv TopBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h;
cp TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/TopVBFBkgScaleFactors_20_10_met.h;
mv TopVBFBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopVBFBkgScaleFactors_8TeV.h;

root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkg/ComputeWWBkgScaleFactor.C+'(2)';
cp WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/WWBkgScaleFactors.h;
mv WWBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h;

sed -ie 's/bool forBDTAna = false/bool forBDTAna = true/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
root -l -b -q $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C+'(2)';
cat $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h DYBkgScaleFactorsBDT.h > DYBkgScaleFactors.h;
mv DYBkgScaleFactors.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h;
sed -ie 's/bool forBDTAna = true/bool forBDTAna = false/' $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.C;
mv DYEstimateTableBDT.txt $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYEstimateTableBDT_20_10_met.txt;
mv DYBkgScaleFactorsBDT.h $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/headers/DYBkgScaleFactorsBDT_20_10_met.h;

rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/DYBkg/ComputeDYBkgScaleFactor.Ce;
rm -f $CMSSW_BASE/src/Smurf/Analysis/HWWlvlv/TopBkg/ComputeTopScaleFactors.Ce;

fi
