#!/bin/bash

# Example:
# ./run_mva.sh 0 130
#
# 1. input files should be available in "data" folder
#  IMPORTANT: this assumes that the data file already inclueds the differential cross-sections for ME method
#  ln -s /smurf/yygao/SmurfV3/datame/ data
#  ln -s /smurf/ceballos/tmva/weights .
# 

# 2. To add the BDT output to the smurfntuples
#    Use 5.28 version root, this is currently available at lxplus

export NJETS=$1;
export MH=$2;

mkdir -p output


# this is the prefix added to the BDT added ntuples
export TAG=ntuples_${MH}train_${NJETS}jets;
export METHODS=BDT;
#export METHODS=KNN,BDT,BDTD,MLPBNN,BDTG;

### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
### samples must be in "data" folder
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
hww${MH}_ME.root
ww_ME.root
ggww_ME.root
wjets_ME.root
ttbar_ME.root
tw_ME.root
wz_ME.root
zz_ME.root
EOF

# ===========================================
# Fill the smurfntuples with the BDT and LR
# ===========================================

for i in `cat list_samples.txt` ; do
  sample=${i%%,*};
  echo "filling MVA information in sample: "  $sample
  ./root-5.28.sh -q -b evaluateMVA_smurf.C+\(\"data/$sample\",${MH},\"${METHODS}\",\"${TAG}\"\);
  echo "filling ME LR in sample: "  $sample
  root -q -b LR.C+\(${MH},\"${TAG}\",\"$sample\",\"data/\",\"output/\",-1\);

done


# ===========================================
# output LR overlay distributions and limit calculation inputs
# ===========================================
echo "make the final plots for ${MH} point"
mkdir -p output/plots;
mkdir -p output/limits;
root -l -q -b makePlots.C\(${MH},\"output/\"\);
    
# ===========================================
# run the limit calculation for both BDT/ME
# http://www.t2.ucsd.edu/tastwiki/bin/view/Smurf/MVAShape
#   cvs co -d LandS UserCode/mschen/LandS
#   cd LandS/
#   make all
# ===========================================

export LANDS=/tas03/home/yygao/WWAnalysis/CMSSW_3_11_3_ME_test/src/LandS/test/lands.exe
cd output/limits/;
echo "Calculating expected limits using BDT output for mH = " ${MH}
${LANDS} -d hww${MH}_BDT.card -tB 1000 -tPB 30 -M Bayesian --doExpectation 1 -t 10 > hww${MH}_BDT.res; 
echo "Calculating expected limits using BDT output for mH = " ${MH}
${LANDS} -d hww${MH}_ME.card -tB 1000 -tPB 30 -M Bayesian --doExpectation 1 -t 10 > hww${MH}_ME.res; 
cd ../../

rm -f list_samples.txt;
