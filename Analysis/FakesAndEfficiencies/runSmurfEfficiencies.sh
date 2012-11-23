#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

#RUNLIST=Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON
#TAG=V00-02-03
#ELE2012=/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/merged_${RUNLIST}.root
#MU2012=/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/merged_${RUNLIST}.root
#MCEE=/smurf/dlevans/LeptonTree/V00-02-02/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged_ee.root
#MCMM=/smurf/dlevans/LeptonTree/V00-02-02/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged_mm.root
#OUTNAME="2012May25V0"

RUNLIST="HCP"
OUTNAME="_V00-02-07_HCP_V0"
ELE2012="/smurf/dlevans/LeptonTree/V00-02-07/SingleElectron/merged_${RUNLIST}.root"
#MU2012="/smurf/dlevans/LeptonTree/V00-02-07/SingleMu/merged_${RUNLIST}.root"

OUTNAME="_V00-02-07_trigNameFix_HCP_V1"
MU2012="/smurf/dlevans/LeptonTree/V00-02-07_trigNameFix/SingleMu/merged_${RUNLIST}.root"

MCEE=/smurf/dlevans/LeptonTree/V00-02-07/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/merged_ee.root
MCMM=/smurf/dlevans/LeptonTree/V00-02-07/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/merged_mm.root
PU52X=8
PU53X=9


#
# options
#

MuonTagAndProbe=1
MuonFakeRate=2
ElectronTagAndProbe=3
ElectronFakeRate=4
MuonTagAndProbeMC=5
ElectronTagAndProbeMC=6

#
# electrons
#

ElectronIsoDenominator2012=30
ElectronIDDenominator2012=31
ElectronFO=32
ElectronIDIsoDenominator=33
ElectronIDIsoRndTagDenominator=34
ElectronIDIsoRndTagDblNoDzDenominator=35

ElectronIsoNumerator2012=40
ElectronIDNumerator2012=41
ElectronIDIsoNumerator2012=42
ElectronTrigSglNumerator2012=43
ElectronTrigDblLeadNumerator2012=44
ElectronTrigDblTrailNumerator2012=45
ElectronTrigDblNumerator2012=46

# ID wrt Iso
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,0,2,2,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_ID\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIsoDenominator2012,$ElectronIDNumerator2012\,60,120,\"$MCEE\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MCEE}\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_ID\",\"png\",1,${PU53X},0,$ElectronTagAndProbeMC,$ElectronIsoDenominator2012,$ElectronIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_ID/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_ID/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_ID/eff.root\"\)

# Iso wrt ID
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,0,2,2,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_Iso\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDDenominator2012,$ElectronIsoNumerator2012\,60,120,\"$MCEE\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MCEE}\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_Iso\",\"png\",1,${PU53X},0,$ElectronTagAndProbeMC,$ElectronIDDenominator2012,$ElectronIsoNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_Iso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_Iso/eff.root\"\)

# Single Trigger wrt ID+Iso
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoDenominator,$ElectronTrigSglNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblLeadNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblTrailNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)

# Double Trigger (last step = dz cut) wrt ID+Iso + trailing leg (before dz cut) (random tag)
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDblNoDzDenominator,$ElectronTrigDblNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)


#
# muons
#

MuonIsoDenominator2012=10
MuonIDDenominator2012=11
MuonFO=12
MuonIDIsoDenominator=13
MuonIDIsoRndTagDenominator=14
MuonIDIsoRndTagDblNoDzDenominator=15

MuonIsoNumerator2012=20
MuonIDNumerator2012=21
MuonIDIsoNumerator2012=22
MuonTrigSglNumerator2012=23
MuonTrigDblLeadNumerator2012=24
MuonTrigDblTrailNumerator2012=25
MuonTrigDblNumerator2012=26

# ID wrt Iso
#root -b -q plotEff.C+\(\"MuonOffline.bins\",2,0,2,2,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_ID\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIsoDenominator2012,$MuonIDNumerator2012\,60,120,\"$MCMM\"\)
#root -b -q plotEff.C+\(\"MuonOffline.bins\",0,0,0,0,\"${MCMM}\",\"HWW_MuonMC${OUTNAME}_NM1Eff_ID\",\"png\",1,${PU53X},0,$MuonTagAndProbeMC,$MuonIsoDenominator2012,$MuonIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_ID/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_ID/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_ID/eff.root\"\)

# Iso wrt ID
#root -b -q plotEff.C+\(\"MuonOffline.bins\",2,0,2,2,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_Iso\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDDenominator2012,$MuonIsoNumerator2012\,60,120,\"$MCMM\"\)
#root -b -q plotEff.C+\(\"MuonOffline.bins\",0,0,0,0,\"${MCMM}\",\"HWW_MuonMC${OUTNAME}_NM1Eff_Iso\",\"png\",1,${PU53X},0,$MuonTagAndProbeMC,$MuonIDDenominator2012,$MuonIsoNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_Iso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_Iso/eff.root\"\)

# Single Trigger wrt ID+Iso
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoDenominator,$MuonTrigSglNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblLeadNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)

#DATASETS="DoubleMu_Run2012A-13Jul2012-v1_AOD_190456_193621 \
#DoubleMu_Run2012A-recover-06Aug2012-v1_AOD_190782_190949 \
#DoubleMu_Run2012B-13Jul2012-v4_AOD_193834_196531 \
#DoubleMu_Run2012C-PromptReco-v2_AOD_198934_202950 \
#DoubleMu_Run2012C-24Aug2012-v1_AOD_198022_198523"
#for DATASET in $DATASETS; do
#    DATAPATH=/smurf/dlevans/LeptonTree/V00-02-07/$DATASET/merged_${RUNLIST}.root
#    TESTOUTNAME=${OUTNAME}_${DATASET}
#    # default syst jet 15
#    root -b -q plotEff.C+\(\"MuonFR.bins\",0,0,0,0,\"${DATAPATH}\",\"HWW_Muon${TESTOUTNAME}_FRJet15_IDISO\",\"png\",1,0,0,$MuonFakeRate,$MuonFO15,$MuonIDIsoNumerator2012\)
#    root -b -q plotDataMC.C+\(\"HWW_Muon${TESTOUTNAME}_FRJet15_IDISO/extra\",\"null\",\"HWW_Muon${TESTOUTNAME}_FRJet15_IDISO/eff.root\"\)
#    root -b -q printEff.C+\(\"HWW_Muon${TESTOUTNAME}_FRJet15_IDISO/extra\",\"HWW_Muon${TESTOUTNAME}_FRJet15_IDISO/eff.root\",\"HWW_Muon${TESTOUTNAME}_FRJet15_IDISO/eff.root\"\)
#done

#MU2012="/smurf/dlevans/LeptonTree/V00-02-07/SingleMu_Run2012A-13Jul2012-v1_AOD_190456_193621/merged_${RUNLIST}.root"
#OUTNAME="_V00-02-07_Run2012A_V0"
#MU2012="/smurf/dlevans/LeptonTree/V00-02-07/SingleMu_Run2012B-13Jul2012-v1_AOD_193834_196531/merged_${RUNLIST}.root"
#OUTNAME="_V00-02-07_Run2012B_V0"
#MU2012="/smurf/dlevans/LeptonTree/V00-02-07/SingleMu_Run2012C-24Aug2012-v1_AOD_198022_198523/merged_${RUNLIST}.root"
#OUTNAME="_V00-02-07_Run2012Cv1_V0"
#MU2012="/smurf/dlevans/LeptonTree/V00-02-07/SingleMu_Run2012C-24Aug2012-v1_AOD_198022_198523/merged_${RUNLIST}.root"
#OUTNAME="_V00-02-07_Run2012Cv2_V0"


root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblTrailNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)

# Double Trigger (last step = dz cut) wrt ID+Iso + trailing leg (before dz cut) (random tag)
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDblNoDzDenominator,$MuonTrigDblNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)

#rm *.so *.d

