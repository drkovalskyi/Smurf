import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")

#
# general config
#
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START52_V5::All"
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
        'file:/smurf/dlevans/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/659/5E3B8EB2-BF82-E111-AD9D-001D09F2932B.root'
    )
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.electronFilters * process.leptonTreeMakerSequenceData2012 * process.pfParticleSelectionSequence * process.eleIsoSequence * process.leptonTreeMaker2012)

