
#
# muons
#

This="HWW_Muon_V00-02-09_Moriond_V0_NM1Eff_ID"
HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Muon_V00-02-07_HCP_V0_NM1Eff_ID"
root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

This="HWW_Muon_V00-02-09_Moriond_V0_NM1Eff_Iso"
HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Muon_V00-02-07_HCP_V0_NM1Eff_Iso"
root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

This="HWW_Muon_V00-02-09_Moriond_V0_NM1Eff_TrigSgl"
HCP="/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/HWW_Muon_V00-02-07_trigNameFix_HCP_V1_NM1Eff_TrigSgl"
root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

This="HWW_Muon_V00-02-09_Moriond_V0_NM1Eff_TrigLeadDbl"
HCP="/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/HWW_Muon_V00-02-07_trigNameFix_HCP_V1_NM1Eff_TrigLeadDbl"
root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

This="HWW_Muon_V00-02-09_Moriond_V0_NM1Eff_TrigTrailDbl"
HCP="/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/HWW_Muon_V00-02-07_trigNameFix_HCP_V1_NM1Eff_TrigTrailDbl"
root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

#
# electrons
#

#This="HWW_Electron_V00-02-09_Moriond_V0_NM1Eff_ID"
#HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Electron_V00-02-07_HCP_V0_NM1Eff_ID"
#root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

#This="HWW_Electron_V00-02-09_Moriond_V0_NM1Eff_Iso"
#HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Electron_V00-02-07_HCP_V0_NM1Eff_Iso"
#root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

#This="HWW_Electron_V00-02-09_Moriond_V0_NM1Eff_TrigSgl"
#HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Electron_V00-02-07_HCP_V0_NM1Eff_TrigSgl"
#root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

#This="HWW_Electron_V00-02-09_Moriond_V0_NM1Eff_TrigLeadDbl"
#HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Electron_V00-02-07_HCP_V0_NM1Eff_TrigLeadDbl"
#root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

#This="HWW_Electron_V00-02-09_Moriond_V0_NM1Eff_TrigTrailDbl"
#HCP="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Electron_V00-02-07_HCP_V0_NM1Eff_TrigTrailDbl"
#root -b -q printEff.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${This}/compare\",\"${HCP}/eff.root\",\"${This}/eff.root\"\)

#rm *.so *.d

