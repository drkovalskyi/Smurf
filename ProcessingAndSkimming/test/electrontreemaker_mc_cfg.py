import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")


#
# general config
#
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_42_V14::All"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#
# lepton maker
#

process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

#
# input
#
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'file:/smurf/dlevans/488D4763-00B1-E011-B110-BCAEC532971C.root'
      'file:/smurf/cerati/Run2011B_Photon_AOD_PromptReco-v1_000_176_799_0099D394-BDE5-E011-82EC-BCAEC5329703.root',
      'file:/smurf/cerati/Run2011A_Photon_AOD_PromptReco-v4_000_167_913_24574F08-84A3-E011-AE85-003048F0258C.root'
    )
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.leptonTreeMakerSequenceForElectronMC)

