#!/bin/sh  

####################################################################################
# Step 1: make cards with prefit qqWW + postfit others and only qqWW syst 
#
# make 3 sets of cards : 
#    - newdefault : entire mT mll region 
#    - CR1(remove CR2)
#    - CR2(remove CR1)  
####################################################################################

#
# get postfit shape from nominal card  
# 
../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name nominal -d 125/hwwof_0j_shape_8TeV.txt

rm nominal_nominalShape.root
rm nominal_maxllfit.root
rm nominal_fittedShape_mu0.root

#
# make cards with only qqWW syst 
#
root -b -q remakecards.C\(1\) # newdefault
root -b -q remakecards.C\(2\) # CR1
root -b -q remakecards.C\(3\) # CR2


#
# Fit using cards made above  
#
../fixPath.pl 125 

../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name newdefault -d 125/hwwof_0j_shape_8TeV_newdefault.txt
rm newdefault_nominalShape.root
rm newdefault_maxllfit.root
rm newdefault_fittedShape_mu0.root
../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name CR1 -d 125/hwwof_0j_shape_8TeV_CR1.txt
rm CR1_CR1Shape.root
rm CR1_maxllfit.root
rm CR2_fittedShape_mu0.root
../../../LandS/test/lands.exe --bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 --name CR2 -d 125/hwwof_0j_shape_8TeV_CR2.txt
rm CR2_nominalShape.root
rm CR2_maxllfit.root
rm CR2_fittedShape_mu0.root

####################################################################################
# Step 2: generate toys 
#
# IMPORTANT 
#   - the card to generate toys should be made by hands
#   - put the the card name after -d
#   - it's more convenient if the card name contains TESTNAME below 
#
####################################################################################

TESTNAME="FullqqWWSyst"

# generate toys 

../../../LandS/test/lands.exe -d 125/hwwof_0j_shape_8TeV_foytoygeneration_FullqqWWSyst.txt -M Hybrid  -m 125 --minuitSTRATEGY 0  --bWriteToys 1 -n $TESTNAME --nToysForCLsb 1000 --nToysForCLb 1 --singlePoint 1 --seed 12345  -rMin 0 -rMax 5  --freq


####################################################################################
# Step 3: get postfit uncertainty
####################################################################################

TOTALTOY=1000

# make cards using psuedo-data
# In the script, make sure to use correct toys 
[ ! -d "125/toycards_${TESTNAME}" ] && echo Create "125/toycards_${TESTNAME}" && mkdir 125/toycards_${TESTNAME}
root -b -q makecardfortoy_runAll.C\($TOTALTOY,\"$TESTNAME\"\)

# do fit 
./fit.sh $TOTALTOY $TESTNAME

# get postfit histograms  
root -b -q getposthisto_runAll.C\($TOTALTOY,\"$TESTNAME\"\) 

# get postfit uncertainties  
root -b -q getPostUncert_runAll.C\($TOTALTOY,\"$TESTNAME\"\)


####################################################################################
# Step 4: draw plots
####################################################################################
root -b -q compareshapes_runAll.C
