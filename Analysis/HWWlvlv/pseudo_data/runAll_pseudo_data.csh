#!/bin/tcsh

cd /home/ceballos/releases/CMSSW_4_2_8_patch4/src/Smurf/Analysis/HWWlvlv;

@ count = 1;
while ($count < 23)
 echo $count
 setenv THEDATA /data/smurf/ceballos/ntuples/pseudo-data_4700ipb/pseudo_data_$count.root;
 /home/ceballos/releases/CMSSW_4_2_8_patch4/src/Smurf/Analysis/HWWlvlv/scripts/computeScaleFactors.sh; 
 ./runAllJobsOneJetBin.sh 0 $count;
 ./runAllJobsOneJetBin.sh 1 $count;
 
 @ count = $count + 1;
 
end
