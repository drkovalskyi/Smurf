import FWCore.ParameterSet.Config as cms

leptonTreeMaker = cms.EDProducer('LeptonTreeMaker',

    #
    # input tags
    #

    electronsInputTag           =  cms.InputTag('gsfElectrons'),
    muonsInputTag               =  cms.InputTag('muons'),
    photonsInputTag             =  cms.InputTag('photons'),
    genParticlesInputTag        =  cms.InputTag('genParticles'),
    primaryVertexInputTag       =  cms.InputTag('offlinePrimaryVertices'),
    rhoIsoAllInputTag           =  cms.InputTag('kt6PFJets',   'rho', 'RECO'),
    rhoIsoAllCentralInputTag    =  cms.InputTag('kt6PFJetsDeterministicIso',   'rho'),
    rhoIsoNeutralInputTag       =  cms.InputTag('kt6PFJetsCentralNeutral',     'rho'),
    rhoIsoNeutral2011InputTag   =  cms.InputTag('kt6PFJetsCentralNeutral2011', 'rho'),
    pfCandsInputTag             =  cms.InputTag('particleFlow'),
    pfNoPUCandsInputTag         =  cms.InputTag('pfNoPileUp'),

    metInputTag                 =  cms.InputTag('pfMet'),
    jetsInputTag                =  cms.InputTag('ak5PFJets'),

    conversionsInputTag      =  cms.InputTag('allConversions'),
    beamSpotInputTag      =  cms.InputTag('offlineBeamSpot'),

    pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3'),

    #
    # define triggers
    # to record branches for
    #

    photonTriggers  = cms.untracked.VInputTag(
        cms.InputTag('HLT_Photon20_CaloIdVL_IsoL_v*::HLT_Photon20_CaloIdVL_IsoL'),
        cms.InputTag('HLT_Photon20_CaloIdVL_IsoL_v*::HLT_Photon30_CaloIdVL_IsoL'),
        cms.InputTag('HLT_Photon20_CaloIdVL_IsoL_v*::HLT_Photon50_CaloIdVL_IsoL'),
        cms.InputTag('HLT_Photon20_CaloIdVL_IsoL_v*::HLT_Photon75_CaloIdVL_IsoL'),
        cms.InputTag('HLT_Photon20_CaloIdVL_IsoL_v*::HLT_Photon90_CaloIdVL_IsoL')),

    muTriggers      = cms.untracked.VInputTag(
        cms.InputTag('HLT_Mu8_v*::HLT_Mu8'),
        cms.InputTag('HLT_Mu15_v*::HLT_Mu15')),

    eleTriggers     = cms.untracked.VInputTag(
        cms.InputTag('HLT_Ele8_CaloIdL_CaloIsoVL_v*::HLT_Ele8_CaloIdL_CaloIsoVL'),
        cms.InputTag('HLT_Ele17_CaloIdL_CaloIsoVL_v*::HLT_Ele17_CaloIdL_CaloIsoVL'),
        cms.InputTag('HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*::HLT_Ele8_CaloIdL_CaloIsoVL_Jet40'))

)
