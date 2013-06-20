#!/bin/sh 

TOTALTOY=$1
TESTNAME=$2

# set up fit output directory
[ ! -d "fitoutput_${TESTNAME}" ] && echo Create "fitoutput_${TESTNAME}" && mkdir fitoutput_${TESTNAME}

#do fit 
for ((x=0; x<$TOTALTOY; x++))
do 
    echo "... fitting toy $x"
     
    #  nominal fit for uncertainties to other backgrounds 
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_$x.txt > fitoutput_${TESTNAME}/log_$x.txt 2>&1 
   
    # fit newdefault
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_newdefault_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_newdefault_$x.txt > fitoutput_${TESTNAME}/log_newdefault_$x.txt 2>&1 

    #  fit CR1 to get nuisances  
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_CR1_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_CR1_$x.txt > fitoutput_${TESTNAME}/log_CR1_$x.txt 2>&1 

    #  fit CR2 to get nuisances  
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_CR2_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_CR2_$x.txt > fitoutput_${TESTNAME}/log_CR2_$x.txt 2>&1 
    
    #  fit newdefault to using CR1 nuisances  
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_newdefault_usingCR1_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_newdefault_$x.txt --ManualNuisanceFeeding fitoutput_${TESTNAME}/log_CR1_$x.txt > fitoutput_${TESTNAME}/log_newdefault_usingCR1_$x.txt 2>&1 
   
    #  fit newdefault to using CR2 nuisances  
    ../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name fitoutput_${TESTNAME}/fit_newdefault_usingCR2_$x -d 125/toycards_${TESTNAME}/hwwof_0j_shape_8TeV_newdefault_$x.txt --ManualNuisanceFeeding fitoutput_${TESTNAME}/log_CR2_$x.txt > fitoutput_${TESTNAME}/log_newdefault_usingCR2_$x.txt 2>&1 
    
    # clean-up 
    rm fitoutput_${TESTNAME}/*_nominalShape.root 
    rm fitoutput_${TESTNAME}/*_maxllfit.root 
    rm fitoutput_${TESTNAME}/*_fittedShape_mu0.root 

done 
