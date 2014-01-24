#!/bin/sh

export YEAR=$1;

if [ $YEAR == 2012 ]; then
root -l -q -b Smurf/Analysis/Zll/GJScaleFactorA.C+'("Ana/nt_scripts/ntuples_53x/data_ph50.root","Ana/nt_scripts/ntuples_53x/data_ph75.root","Ana/nt_scripts/ntuples_53x/data_ph.root")';
hadd -f outputPH0.root output1.root output2.root output3.root;
rm -f output1.root output2.root output3.root;
root -l -q -b Smurf/Analysis/Zll/GJScaleFactorB.C+'("Ana/nt_scripts/ntuples_53x/outputPH0.root","Ana/nt_scripts/ntuples_53x/data_skim10.root",2012)';
mv outputG2012.root Ana/nt_scripts/ntuples_53x/;
rm -f outputPH0.root;
hadd -f Ana/nt_scripts/ntuples_53x/hww_syst_skim10_gj.root Ana/nt_scripts/ntuples_53x/hww_syst_skim10.root Ana/nt_scripts/ntuples_53x/outputG2012.root;

else
root -l -q -b Smurf/Analysis/Zll/GJScaleFactorA.C+'("Ana/nt_scripts/ntuples_42x_v9/data_ph50.root","Ana/nt_scripts/ntuples_42x_v9/data_ph75.root","Ana/nt_scripts/ntuples_42x_v9/data_ph.root")';
hadd -f outputPH0.root output1.root output2.root output3.root;
rm -f output1.root output2.root output3.root;
root -l -q -b Smurf/Analysis/Zll/GJScaleFactorB.C+'("Ana/nt_scripts/ntuples_42x_v9/outputPH0.root","Ana/nt_scripts/ntuples_42x_v9/data_skim10.root",2011)'
mv outputG2011.root Ana/nt_scripts/ntuples_42x_v9/;
rm -f outputPH0.root;
hadd -f Ana/nt_scripts/ntuples_42x_v9/hww_syst_skim10_gj.root Ana/nt_scripts/ntuples_42x_v9/hww_syst_skim10.root Ana/nt_scripts/ntuples_42x_v9/outputG2011.root;

fi
