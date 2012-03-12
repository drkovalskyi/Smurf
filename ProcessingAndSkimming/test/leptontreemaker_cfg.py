import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")


process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/smurf/dlevans/488D4763-00B1-E011-B110-BCAEC532971C.root'
    )
)

#
# rho
#

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsDeterministicIso = kt4PFJets.clone(
    rParam = 0.6,
    doAreaFastjet = True,
    doRhoFastjet = True,
    voronoiRfact = 0.9,
    Rho_EtaMax = 2.5,
    Ghost_EtaMax = 2.5
)

process.leptonTreeMaker = cms.EDProducer('LeptonTreeMaker',

    #
    # input tags
    #

    electronsInputTag      =  cms.InputTag('gsfElectrons'),
    muonsInputTag          =  cms.InputTag('muons'),
    photonsInputTag        =  cms.InputTag('photons'),
    primaryVertexInputTag  =  cms.InputTag('offlinePrimaryVertices'),
    rhoInputTag            =  cms.InputTag('kt6PFJetsDeterministicIso', 'rho'),
    pfCandsInputTag        =  cms.InputTag('particleFlow'),
    metInputTag            =  cms.InputTag('pfMet'),
    jetsInputTag           =  cms.InputTag('ak5PFJets'),

    pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3'),

    #
    # define triggers
    #

    electronTPTriggerNames = cms.untracked.vstring( 
        "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*",
        "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v*",
        "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v*",
        "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v*"),

    muonTPTriggerNames     = cms.untracked.vstring(
        "HLT_Mu15_v*",
        "HLT_Mu24_v*",
        "HLT_Mu30_v*",
        "HLT_Mu40_v*",
        "HLT_Mu40_eta2p1_v*"
        "HLT_IsoMu17_v*",
        "HLT_IsoMu17_eta2p1_v*",
        "HLT_IsoMu24_v*",
        "HLT_IsoMu30_eta2p1_v*"),

    electronFRTriggerNames = cms.untracked.vstring(
        "HLT_Ele8_CaloIdL_CaloIsoVL_v*",
        "HLT_Ele17_CaloIdL_CaloIsoVL_v8",
        "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v8",
        "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*"),

    muonFRTriggerNames     = cms.untracked.vstring(
        "HLT_Mu8_v*",
        "HLT_Mu15_v*"),

    photonTriggerNames     = cms.untracked.vstring(
        "HLT_Photon20_CaloIdVL_IsoL_v*",
        "HLT_Photon30_CaloIdVL_IsoL_v*",
        "HLT_Photon50_CaloIdVL_IsoL_v*",
        "HLT_Photon75_CaloIdVL_IsoL_v*",
        "HLT_Photon90_CaloIdVL_IsoL_v*"),

)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)
  
process.p = cms.Path(process.kt6PFJetsDeterministicIso * process.leptonTreeMaker)
process.e = cms.EndPath(process.out)

