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
        #'file:/smurf/cerati/Run2012A_DoubleElectron_AOD_PromptReco-v1_000_191_247_04825687-3588-E111-82CE-BCAEC518FF63.root'
        'root://xrootd.unl.edu//store/data/Run2012C/DoubleElectron/AOD/PromptReco-v2/000/200/075/C85ECB28-46DE-E111-A416-0025901D62A0.root'
    )
)

process.leptonTreeMaker2012.runFR = cms.bool(False)
process.leptonTreeMaker2012.runGamma = cms.bool(False)
process.leptonTreeMaker2012.pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3Residual')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.electronPath = cms.Path(process.electronFilters * process.leptonTreeMakerSequenceData2012 * process.leptonTreeMaker2012)
process.photonPath = cms.Path(process.photonFilters * process.leptonTreeMakerSequenceData2012 * process.leptonTreeMaker2012)

