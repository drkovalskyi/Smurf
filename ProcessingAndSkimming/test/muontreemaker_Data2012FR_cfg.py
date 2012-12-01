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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#
# lepton maker
#

process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

#
# input
#
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://xrootd.unl.edu//store/data/Run2012C/SingleMu/AOD/PromptReco-v2/000/202/237/DA293EDC-27F9-E111-AFA3-001D09F24682.root'
        'root://xrootd.unl.edu//store/data/Run2012C/DoubleMu/AOD/PromptReco-v2/000/200/075/5A5F1D8A-D3DD-E111-8F14-BCAEC5329709.root'
        #'file:/smurf/dlevans/LeptonTree/singlemu_190702_66_37385741.root'
        #'file:/smurf/dlevans/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/659/5E3B8EB2-BF82-E111-AD9D-001D09F2932B.root'
    )
)

process.leptonTreeMaker2012.runTP = cms.bool(False)
process.leptonTreeMaker2012.runGamma = cms.bool(False)
process.leptonTreeMaker2012.pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3Residual')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.muonFilters * process.leptonTreeMakerSequenceData2012FR * process.leptonTreeMaker2012)

