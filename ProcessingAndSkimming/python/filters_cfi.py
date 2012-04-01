import FWCore.ParameterSet.Config as cms

#
# photon filters
#

highPtPhotons = cms.EDFilter("CandViewSelector",
src = cms.InputTag("photons"),
cut = cms.string("pt > 40.0")
)

highPtPhotonFilter = cms.EDFilter("CandViewCountFilter",
src = cms.InputTag("highPtPhotons"),
minNumber = cms.uint32(1),
)

centralPhotons = cms.EDFilter("CandViewSelector",
src = cms.InputTag("highPtPhotons"),
cut = cms.string("abs(eta)<1.479")
)

centralPhotonFilter = cms.EDFilter("CandViewCountFilter",
src = cms.InputTag("centralPhotons"),
minNumber = cms.uint32(1),
)

#
# electron filters
#

highPtElectrons = cms.EDFilter("CandViewSelector",
srs = cms.InputTag("gsfElectrons"),
cut = cms.string("pt > 20.0")
)

highPtElectronFilter = cms.EDFilter("CandViewCountFilter",
src = cms.InputTag("highPtElectrons"),
minNumber = cms.uint32(1)
)

#
# muon filters
#

highPtMuons = cms.EDFilter("CandViewSelector",
srs = cms.InputTag("muons"),
cut = cms.string("pt > 20.0")
)

highPtMuonFilter = cms.EDFilter("CandViewCountFilter",
src = cms.InputTag("highPtMuons"),
minNumber = cms.uint32(1)
)

#
# hlt filters
#

hltHighLevel = cms.EDFilter("HLTHighLevel",
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
