import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff') # Import the Jet RECO modules

#
# general config
#

process.GlobalTag.globaltag = "GR_R_42_V14::All"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.kt6PFJets.doRhoFastjet  = True                # Turn-on the FastJet density calculation

#
# input
#

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/smurf/dlevans/488D4763-00B1-E011-B110-BCAEC532971C.root'
      #'file:/smurf/cerati/Run2011B_Photon_AOD_PromptReco-v1_000_176_799_0099D394-BDE5-E011-82EC-BCAEC5329703.root',
      #'file:/smurf/cerati/Run2011A_Photon_AOD_PromptReco-v4_000_167_913_24574F08-84A3-E011-AE85-003048F0258C.root'
    )
)

#
# filters
#

process.hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),

  # provide list of HLT paths (or patterns) you want
    HLTPaths = cms.vstring(
        "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*",
        "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v*",
        "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v*",
        "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v*",
        "HLT_Mu15_v*",
        "HLT_Mu24_v*",
        "HLT_Mu30_v*",
        "HLT_Mu40_v*",
        "HLT_Mu40_eta2p1_v*"
        "HLT_IsoMu17_v*",
        "HLT_IsoMu17_eta2p1_v*",
        "HLT_IsoMu24_v*",
        "HLT_IsoMu30_eta2p1_v*",
        "HLT_Ele8_CaloIdL_CaloIsoVL_v*",
        "HLT_Ele17_CaloIdL_CaloIsoVL_v8",
        "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v8",
        "HLT_Mu8_v*",
        "HLT_Mu15_v*",
        "HLT_Photon20_CaloIdVL_IsoL_v*",
        "HLT_Photon30_CaloIdVL_IsoL_v*",
        "HLT_Photon50_CaloIdVL_IsoL_v*",
        "HLT_Photon75_CaloIdVL_IsoL_v*",
        "HLT_Photon90_CaloIdVL_IsoL_v*"),

    eventSetupPathsKey = cms.string(''),    # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),                 # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False)                  # throw exception on unknown path names
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

process.kt6PFJetsDeterministicJEC = kt4PFJets.clone(
    rParam = 0.6,
    doAreaFastjet = True,
    doRhoFastjet = True,
    voronoiRfact = 0.9,
    Rho_EtaMax = 5.0,
    Ghost_EtaMax = 5.0
)

#
# lepton maker
#

process.leptonTreeMaker = cms.EDProducer('LeptonTreeMaker',

    #
    # input tags
    #

    electronsInputTag      =  cms.InputTag('gsfElectrons'),
    muonsInputTag          =  cms.InputTag('muons'),
    photonsInputTag        =  cms.InputTag('photons'),
    primaryVertexInputTag  =  cms.InputTag('offlinePrimaryVertices'),
    rhoIsoInputTag         =  cms.InputTag('kt6PFJetsDeterministicJEC', 'rho'),
    rhoJECInputTag         =  cms.InputTag('kt6PFJetsDeterministicIso', 'rho'),
    pfCandsInputTag        =  cms.InputTag('particleFlow'),
    metInputTag            =  cms.InputTag('pfMet'),
    jetsInputTag           =  cms.InputTag('ak5PFJets'),

    conversionsInputTag      =  cms.InputTag('allConversions'),
    beamSpotInputTag      =  cms.InputTag('offlineBeamSpot'),

    pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3'),

    #this works for crab jobs, can be changed for local submission
    pathToBDTWeights = cms.string('src/Smurf/ProcessingAndSkimming/data'),
    #pathToBDTWeights = cms.string('data'),

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
        "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v8"),

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

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('myOutputFile.root')
#)
#process.out.outputCommands = cms.untracked.vstring( 'drop *' )
#process.out.outputCommands.extend(cms.untracked.vstring('keep *_*_*_SMURF'))
#process.e = cms.EndPath(process.out)

fastJetSequence = cms.Sequence(process.kt6PFJets * process.kt6PFJetsDeterministicIso * process.kt6PFJetsDeterministicJEC)
#process.p = cms.Path(process.hltHighLevel * fastJetSequence * process.leptonTreeMaker)
process.p = cms.Path(fastJetSequence * process.leptonTreeMaker)


