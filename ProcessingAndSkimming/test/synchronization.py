import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.GlobalTag.globaltag = "GR_R_42_V14::All"
process.GlobalTag.globaltag = "GR_R_52_V7::All"
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

process.hww = cms.EDFilter('Smurf::SmurfTreeProducer',
    inputMuonCollection       =  cms.InputTag('muons'),
    inputElectronCollection   =  cms.InputTag('gsfElectrons'),
    inputJetCollection        =  cms.InputTag('ak5PFJets'),
    inputPVCollection         =  cms.InputTag('offlinePrimaryVertices'),
    inputPFCandCollection     =  cms.InputTag('particleFlow'),
    rhoIsoAll                 =  cms.InputTag('kt6PFJets', 'rho'),
    rhoIsoNeutral2011         =  cms.InputTag('kt6PFJetsCentralNeutral2011', 'rho'),
    conversions               =  cms.InputTag('allConversions'),
    beamspot                  =  cms.InputTag('offlineBeamSpot'),

    smurfTree               =  cms.string('test.root'),

#     rhoIsoAllInputTag           =  cms.InputTag('kt6PFJets',   'rho', 'RECO'),
#     rhoIsoAllCentralInputTag    =  cms.InputTag('kt6PFJetsDeterministicIso',   'rho'),
#     rhoIsoNeutralInputTag       =  cms.InputTag('kt6PFJetsCentralNeutral',     'rho'),
#     pfNoPUCandsInputTag         =  cms.InputTag('pfNoPileUp'),

#     metInputTag                 =  cms.InputTag('pfMet'),


#     pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3'),

)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      # 'file:/afs/cern.ch/cms/higgs/higgsww/dmytro/samples/A6DE4085-8191-E111-BF4E-001E67396D5B.root'
      # 'file:/tas/dmytro/samples/MC/A6DE4085-8191-E111-BF4E-001E67396D5B.root'
      'file:/home/users/dmytro/projects/HWW-53X/src/A6DE4085-8191-E111-BF4E-001E67396D5B.root'
    )
)

#process.leptonTreeMaker2011.pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3Residual')

# process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.hww)

