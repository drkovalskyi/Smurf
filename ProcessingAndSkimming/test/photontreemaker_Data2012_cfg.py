import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")


#
# general config
#
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_52_V7::All"
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#
# lepton maker
#

process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

#
# PFIso
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')

#
# input
#
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/smurf/cerati/Run2012A_DoubleElectron_AOD_PromptReco-v1_000_191_247_04825687-3588-E111-82CE-BCAEC518FF63.root'
    )
)

process.leptonTreeMaker2012.pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3Residual')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.photonFilters * process.leptonTreeMakerSequenceData2012 * process.pfParticleSelectionSequence * process.eleIsoSequence * process.leptonTreeMaker2012)

