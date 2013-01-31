#!/bin/bash

#TAG="V00-02-06"
#root -b -q skimMCLeptonTree.C+\(0,\"/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_ee.root\"\)
#root -b -q skimMCLeptonTree.C+\(1,\"/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_mm.root\"\)

#TAG="V00-02-07"
#DATASET="DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM"
#root -b -q skimMCLeptonTree.C+\(0,\"/smurf/dlevans/LeptonTree/${TAG}/${DATASET}\",\"merged_ee.root\"\)
#root -b -q skimMCLeptonTree.C+\(1,\"/smurf/dlevans/LeptonTree/${TAG}/${DATASET}\",\"merged_mm.root\"\)

#root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/FR_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_e.root\",37509,10000000\)
#root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/FR_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_m.root\",37509,10000000\)
#root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/FR_DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/\",\"merged_m.root\",1930.9940,4900000\)
#root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/FR_DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/\",\"merged_e.root\",1930.9940,5000000\)

root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-5to15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/\",\"merged_processed_e.root\",4.2639499E10,689184\)
root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_e.root\",9.8828742E8,361524\)
root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_e.root\",6.609E7,6000000\)
root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_e.root\",8148778.0,5998860\)
root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM/\",\"merged_processed_e.root\",1033680.0,5694864\)
root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM/\",\"merged_processed_e.root\",156293.3,5785732\)
root -b -q skimMCLeptonTree.C+\(4,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_e.root\",34138.15,5714398\)

root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-5to15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/\",\"merged_processed_m.root\",4.2639499E10,689184\)
root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_m.root\",9.8828742E8,361524\)
root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_m.root\",6.609E7,6000000\)
root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_m.root\",8148778.0,5998860\)
root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM/\",\"merged_processed_m.root\",1033680.0,5694864\)
root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM/\",\"merged_processed_m.root\",156293.3,5785732\)
root -b -q skimMCLeptonTree.C+\(5,\"/smurf/dlevans/LeptonTree/V00-02-09/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM/\",\"merged_processed_m.root\",34138.15,5714398\)

