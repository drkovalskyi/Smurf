import FWCore.ParameterSet.Config as cms

leptonTreeMaker = cms.EDProducer('LeptonTreeMaker',

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
    #pathToBDTWeights = cms.string('src/Smurf/ProcessingAndSkimming/data'),
    pathToBDTWeights = cms.string('data'),

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
