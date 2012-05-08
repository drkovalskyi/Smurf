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
# input
#
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/smurf/dlevans/LeptonTree/singlemu_190702_66_37385741.root'
        #'file:/smurf/dlevans/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/659/5E3B8EB2-BF82-E111-AD9D-001D09F2932B.root'
    )
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.muonFilters * process.leptonTreeMakerSequenceData2012 * process.leptonTreeMaker2012)

