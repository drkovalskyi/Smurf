
#
# muons
#

#ICHEP_ID="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Muon2012June22V1_NM1Eff_ID"
#HCP_ID="HWW_Muon_V00-02-07_HCP_V0_NM1Eff_ID"
#root -b -q printEff.C+\(\"${HCP_ID}/compare_ICHEP2012\",\"${ICHEP_ID}/eff.root\",\"${HCP_ID}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_ID}/compare_ICHEP2012\",\"${ICHEP_ID}/eff.root\",\"${HCP_ID}/eff.root\"\)

#ICHEP_Iso="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Muon2012June22V1_NM1Eff_Iso"
#HCP_Iso="HWW_Muon_V00-02-07_HCP_V0_NM1Eff_Iso"
#root -b -q printEff.C+\(\"${HCP_Iso}/compare_ICHEP2012\",\"${ICHEP_Iso}/eff.root\",\"${HCP_Iso}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_Iso}/compare_ICHEP2012\",\"${ICHEP_Iso}/eff.root\",\"${HCP_Iso}/eff.root\"\)

#ICHEP_TrigSgl="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Muon2012June22V1_NM1Eff_TrigSgl"
#HCP_TrigSgl="HWW_Muon_V00-02-07_HCP_V0_NM1Eff_TrigSgl"
#root -b -q printEff.C+\(\"${HCP_TrigSgl}/compare_ICHEP2012\",\"${ICHEP_TrigSgl}/eff.root\",\"${HCP_TrigSgl}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_TrigSgl}/compare_ICHEP2012\",\"${ICHEP_TrigSgl}/eff.root\",\"${HCP_TrigSgl}/eff.root\"\)

#ICHEP_TrigLeadDbl="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Muon2012June22V1_NM1Eff_TrigLeadDbl"
#HCP_TrigLeadDbl="HWW_Muon_V00-02-07_HCP_V0_NM1Eff_TrigLeadDbl"
#root -b -q printEff.C+\(\"${HCP_TrigLeadDbl}/compare_ICHEP2012\",\"${ICHEP_TrigLeadDbl}/eff.root\",\"${HCP_TrigLeadDbl}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_TrigLeadDbl}/compare_ICHEP2012\",\"${ICHEP_TrigLeadDbl}/eff.root\",\"${HCP_TrigLeadDbl}/eff.root\"\)

#ICHEP_TrigTrailDbl="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Muon2012June22V1_NM1Eff_TrigTrailDbl"
#HCP_TrigTrailDbl="HWW_Muon_V00-02-07_HCP_V0_NM1Eff_TrigTrailDbl"
#root -b -q printEff.C+\(\"${HCP_TrigTrailDbl}/compare_ICHEP2012\",\"${ICHEP_TrigTrailDbl}/eff.root\",\"${HCP_TrigTrailDbl}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_TrigTrailDbl}/compare_ICHEP2012\",\"${ICHEP_TrigTrailDbl}/eff.root\",\"${HCP_TrigTrailDbl}/eff.root\"\)


#
# electrons
#

#ICHEP_ID="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Electron2012June22V1_NM1Eff_ID"
#HCP_ID="HWW_Electron_V00-02-07_HCP_V0_NM1Eff_ID"
#root -b -q printEff.C+\(\"${HCP_ID}/compare_ICHEP2012\",\"${ICHEP_ID}/eff.root\",\"${HCP_ID}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_ID}/compare_ICHEP2012\",\"${ICHEP_ID}/eff.root\",\"${HCP_ID}/eff.root\"\)

#ICHEP_Iso="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Electron2012June22V1_NM1Eff_Iso"
#HCP_Iso="HWW_Electron_V00-02-07_HCP_V0_NM1Eff_Iso"
#root -b -q printEff.C+\(\"${HCP_Iso}/compare_ICHEP2012\",\"${ICHEP_Iso}/eff.root\",\"${HCP_Iso}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_Iso}/compare_ICHEP2012\",\"${ICHEP_Iso}/eff.root\",\"${HCP_Iso}/eff.root\"\)

This="HWW_Electron_V00-02-09_Moriond_V0_NM1Eff_TrigSgl"
HCP_TrigSgl="/smurf/dlevans/Efficiencies/V00-02-07_HCP_V0/HWW_Electron_V00-02-07_HCP_V0_NM1Eff_TrigSgl"
root -b -q printEff.C+\(\"${This}/compare_ICHEP2012\",\"${HCP_TrigSgl}/eff.root\",\"${This}/eff.root\"\)
root -b -q plotDataMC.C+\(\"${This}/compare_ICHEP2012\",\"${HCP_TrigSgl}/eff.root\",\"${This}/eff.root\"\)

#ICHEP_TrigLeadDbl="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Electron2012June22V1_NM1Eff_TrigLeadDbl"
#HCP_TrigLeadDbl="HWW_Electron_V00-02-07_HCP_V0_NM1Eff_TrigLeadDbl"
#root -b -q printEff.C+\(\"${HCP_TrigLeadDbl}/compare_ICHEP2012\",\"${ICHEP_TrigLeadDbl}/eff.root\",\"${HCP_TrigLeadDbl}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_TrigLeadDbl}/compare_ICHEP2012\",\"${ICHEP_TrigLeadDbl}/eff.root\",\"${HCP_TrigLeadDbl}/eff.root\"\)

#ICHEP_TrigTrailDbl="/smurf/dlevans/Efficiencies/V00-02-06_V1/HWW_Electron2012June22V1_NM1Eff_TrigTrailDbl"
#HCP_TrigTrailDbl="HWW_Electron_V00-02-07_HCP_V0_NM1Eff_TrigTrailDbl"
#root -b -q printEff.C+\(\"${HCP_TrigTrailDbl}/compare_ICHEP2012\",\"${ICHEP_TrigTrailDbl}/eff.root\",\"${HCP_TrigTrailDbl}/eff.root\"\)
#root -b -q plotDataMC.C+\(\"${HCP_TrigTrailDbl}/compare_ICHEP2012\",\"${ICHEP_TrigTrailDbl}/eff.root\",\"${HCP_TrigTrailDbl}/eff.root\"\)

#rm *.so *.d

