#!/bin/sh 

TOTALTOY=$1
TESTNAME=$2

# set up fit output directory
[ ! -d "fitoutput_${TESTNAME}" ] && echo Create "fitoutput_${TESTNAME}" && mkdir fitoutput_${TESTNAME}

#do fit 
for ((x=0; x<$TOTALTOY; x++))
do 
    echo "... fitting toy $x"

    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_$x.txt > fitoutput_${TESTNAME}/log_$x.txt 2>&1 
    rm fitoutput_${TESTNAME}/fit_${x}_nominalShape.root 
    rm fitoutput_${TESTNAME}/fit_${x}_maxllfit.root 
    rm fitoutput_${TESTNAME}/fit_${x}_fittedShape_mu0.root 
    
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_newdefault_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_newdefault_$x.txt > fitoutput_${TESTNAME}/log_newdefault_$x.txt 2>&1 
    rm fitoutput_${TESTNAME}/fit_newdefault_${x}_nominalShape.root 
    rm fitoutput_${TESTNAME}/fit_newdefault_${x}_maxllfit.root 
    rm fitoutput_${TESTNAME}/fit_newdefault_${x}_fittedShape_mu0.root 

    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_CR1_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_CR1_$x.txt > fitoutput_${TESTNAME}/log_CR1_$x.txt 2>&1 
    rm fitoutput_${TESTNAME}/fit_CR1_${x}_nominalShape.root 
    rm fitoutput_${TESTNAME}/fit_CR1_${x}_maxllfit.root 
    rm fitoutput_${TESTNAME}/fit_CR1_${x}_fittedShape_mu0.root 

    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_CR2_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_CR2_$x.txt > fitoutput_${TESTNAME}/log_CR2_$x.txt 2>&1 
    rm fitoutput_${TESTNAME}/fit_CR2_${x}_nominalShape.root 
    rm fitoutput_${TESTNAME}/fit_CR2_${x}_maxllfit.root 
    rm fitoutput_${TESTNAME}/fit_CR2_${x}_fittedShape_mu0.root 

done 
