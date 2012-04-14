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

    electronTPTriggerNames = cms.untracked.VInputTag( 
        cms.InputTag("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter"),
        cms.InputTag("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter"),
        cms.InputTag("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter")),

    muonTPTriggerNames     = cms.untracked.VInputTag(
        cms.InputTag("HLT_Mu15_v*"),
        cms.InputTag("HLT_Mu24_v*"),
        cms.InputTag("HLT_Mu30_v*"),
        cms.InputTag("HLT_Mu40_v*"),
        cms.InputTag("HLT_Mu40_eta2p1_v*"),
        cms.InputTag("HLT_IsoMu17_v*"),
        cms.InputTag("HLT_IsoMu17_eta2p1_v*"),
        cms.InputTag("HLT_IsoMu24_v*"),
        cms.InputTag("HLT_IsoMu30_eta2p1_v*")),

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

    measureSingleEle    = cms.untracked.VInputTag(),
    measureLeadingDoubleEle = cms.untracked.VInputTag(),
    measureTrailingDoubleEle = cms.untracked.VInputTag(),
    measureDoubleEleDZ = cms.untracked.VInputTag(),

    measureSingleMu24       = cms.untracked.VInputTag(),
    measureSingleMu30       = cms.untracked.VInputTag(),
    measureLeadingDoubleMu  = cms.untracked.VInputTag(),
    measureTrailingDoubleMu = cms.untracked.VInputTag(),
    measureDoubleMuDZ       = cms.untracked.VInputTag(),

)
