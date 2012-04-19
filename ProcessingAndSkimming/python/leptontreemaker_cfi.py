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

    isoValInputTags         = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),

    conversionsInputTag      =  cms.InputTag('allConversions'),
    beamSpotInputTag      =  cms.InputTag('offlineBeamSpot'),

    pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3'),

    #this works for crab jobs, can be changed for local submission
    pathToBDTWeights = cms.string('src/Smurf/ProcessingAndSkimming/data/'),
    #pathToBDTWeights = cms.string('./data/'),

    #
    # define triggers
    # used to select events
    #

    electronFRTriggerNames = cms.untracked.VInputTag(
        cms.InputTag("HLT_Ele8_CaloIdL_CaloIsoVL_v*"),
        cms.InputTag("HLT_Ele17_CaloIdL_CaloIsoVL_v*"),
        cms.InputTag("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*")),

    muonFRTriggerNames     = cms.untracked.VInputTag(
        cms.InputTag("HLT_Mu8_v*"),
        cms.InputTag("HLT_Mu15_v*")),

    photonTriggerNames     = cms.untracked.VInputTag(
        cms.InputTag("HLT_Photon20_CaloIdVL_IsoL_v*"),
        cms.InputTag("HLT_Photon30_CaloIdVL_IsoL_v*"),
        cms.InputTag("HLT_Photon50_CaloIdVL_IsoL_v*"),
        cms.InputTag("HLT_Photon75_CaloIdVL_IsoL_v*"),
        cms.InputTag("HLT_Photon90_CaloIdVL_IsoL_v*")),

    #
    # define triggers
    # to measure in tag and probe 
    #

    muTriggers      = cms.untracked.VInputTag(),
    eleTriggers     = cms.untracked.VInputTag(),

)
