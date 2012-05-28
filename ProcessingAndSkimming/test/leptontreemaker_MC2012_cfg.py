import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")

#
# general config
#
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START52_V9::All"
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
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/2EB0D0B7-C275-E111-8CFB-001A64789E04.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/32AAFDB8-C875-E111-AD40-0025901248FA.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/4208A68B-BC75-E111-ADF0-001A64789498.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/6A89213C-BF75-E111-A798-0025B3E063F0.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/B4EADB63-BE75-E111-A87D-003048D479D6.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/BA7F6AC5-C475-E111-B9A7-001A6478935C.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/DA10D9B6-BD75-E111-BB07-003048673F12.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/DA2D50ED-C075-E111-ACCA-002481E15522.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/F2F8B140-1476-E111-851C-00259020081C.root',
#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/F63CF568-0B76-E111-9385-00E08178C0E3.root'
'file:/smurf/dlevans/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v1/0000/2EB0D0B7-C275-E111-8CFB-001A64789E04.root'
#'rfio:/castor/cern.ch/user/e/emanuele/AODSummer12/00514868-B47A-E111-803E-001D0967DDC3.root'

    )
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.leptonTreeMakerSequenceMC * process.leptonTreeMaker2012)

