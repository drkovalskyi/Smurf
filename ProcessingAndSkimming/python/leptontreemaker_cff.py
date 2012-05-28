import FWCore.ParameterSet.Config as cms

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoJets.Configuration.RecoPFJets_cff import * # Import the Jet RECO modules
kt6PFJets.doRhoFastjet  = True                # Turn-on the FastJet density calculation

#
# rho
#

from RecoJets.JetProducers.kt4PFJets_cfi import *

# for 2011 muon definition
kt6PFJetsCentralNeutral2011 = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True)
kt6PFJetsCentralNeutral2011.Ghost_EtaMax = cms.double(2.5)
kt6PFJetsCentralNeutral2011.Rho_EtaMax = cms.double(2.5)

# for Egamma
kt6PFJetsDeterministicIso = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
kt6PFJetsDeterministicIso.Rho_EtaMax = cms.double(2.5)

# for JEC
kt6PFJetsDeterministicJEC = kt4PFJets.clone(
    rParam = 0.6,
    doAreaFastjet = True,
    doRhoFastjet = True,
    voronoiRfact = 0.9,
    Rho_EtaMax = 5.0,
    Ghost_EtaMax = 5.0
)

fastJetSequence = cms.Sequence(kt6PFJets * kt6PFJetsCentralNeutral2011 * kt6PFJetsDeterministicIso * kt6PFJetsDeterministicJEC)

#
# pf no pileup collection
#

from CommonTools.ParticleFlow.pfNoPileUp_cff import *
pfPileUp.PFCandidates = "particleFlow"
pfNoPileUp.bottomCollection = "particleFlow"

#
# lepton maker
#

from Smurf.ProcessingAndSkimming.leptontreemaker_cfi import *

# For 2011
leptonTreeMaker2011 = leptonTreeMaker.clone()
leptonTreeMaker2011.rhoIsoNeutralInputTag       =  cms.InputTag('kt6PFJetsCentralNeutral2011',     'rho')

# For 2012
leptonTreeMaker2012 = leptonTreeMaker.clone()

leptonTreeMaker2012.photonTriggers     = cms.untracked.VInputTag(
        cms.InputTag("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v*::HLT_Photon22_R9Id90_HE10_Iso40_EBOnly"),
        cms.InputTag("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*::HLT_Photon36_R9Id90_HE10_Iso40_EBOnly"),
        cms.InputTag("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*::HLT_Photon50_R9Id90_HE10_Iso40_EBOnly"),
        cms.InputTag("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*::HLT_Photon75_R9Id90_HE10_Iso40_EBOnly"),
        cms.InputTag("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*::HLT_Photon90_R9Id90_HE10_Iso40_EBOnly"))

leptonTreeMaker2012.muTriggers     = cms.untracked.VInputTag(

        # HLT_Mu17_Mu8_v*
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL1sL1DoubleMu10MuOpen:HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen'),
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8:HLT_Mu17_Mu8_TrailingLeg'),
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:HLT_Mu17_Mu8_LeadingLeg'),
        cms.InputTag('HLT_Mu17_Mu8_v*::HLT_Mu17_Mu8'),

        # HLT_Mu17_TkMu8_v*
        cms.InputTag('HLT_Mu17_TkMu8_v*:hltL2fL1sDoubleMu10MuOpenL1f0L2Filtered10:HLT_Mu17_TkMu8_TrailingLeg'),
        cms.InputTag('HLT_Mu17_TkMu8_v*:hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17:HLT_Mu17_TkMu8_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TkMu8_v*:hltDiMuonGlbFiltered17TrkFiltered8:HLT_Mu17_TkMu8_TrailingLegTrkFiltered'),
        cms.InputTag('HLT_Mu17_TkMu8_v*::HLT_Mu17_TkMu8'),

        # HLT_IsoMu24_eta2p1_v*
        cms.InputTag('HLT_IsoMu24_eta2p1_v*:hltL1sMu16Eta2p1:HLT_IsoMu24_eta2p1_L1sMu16Eta2p1'),
        cms.InputTag('HLT_IsoMu24_eta2p1_v*::HLT_IsoMu24_eta2p1'),

        # for FR
        cms.InputTag("HLT_Mu8_v*::HLT_Mu8"),
        cms.InputTag("HLT_Mu17_v*::HLT_Mu17")

)

leptonTreeMaker2012.eleTriggers     = cms.untracked.VInputTag(

        # HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltL1sL1DoubleEG137:HLT_Ele17_Ele8_L1sL1DoubleEG137'),
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter:HLT_Ele17_Ele8_LeadingLeg'),
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter:HLT_Ele17_Ele8_TrailingLeg'),
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Ele17_Ele8'),

        # HLT_Ele27_WP80_v*
        cms.InputTag('HLT_Ele27_WP80_v*:hltL1sL1SingleEG20ORL1SingleEG22:HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22'),
        cms.InputTag('HLT_Ele27_WP80_v*::HLT_Ele27_WP80'),

        # HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*
        cms.InputTag('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter:HLT_Ele17_Ele8_Mass50_LeadingLeg'),
        cms.InputTag('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter:HLT_Ele17_Ele8_Mass50_TrailingLeg'),

        # HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*
        cms.InputTag('HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter:HLT_Ele20_SC4_Mass50_LeadingLeg'),
        cms.InputTag('HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter:HLT_Ele20_SC4_Mass50_TrailingLeg'),

        # HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3
        cms.InputTag('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter:HLT_Ele32_SC17_Mass50_LeadingLeg'),
        cms.InputTag('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter:HLT_Ele32_SC17_Mass50_TrailingLeg'),

        # for FR
        cms.InputTag("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*::HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30"),
        cms.InputTag("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"),
        cms.InputTag("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*::HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30"),
        cms.InputTag("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"),
        cms.InputTag("HLT_Ele8_CaloIdT_TrkIdVL_v*::HLT_Ele8_CaloIdT_TrkIdVL")
)

#
# JetID
#


from CMGTools.External.pujetidsequence_cff import *
puJetId.jets = cms.InputTag("ak5PFJets")
puJetId.applyJec = cms.bool(True)
puJetId.inputIsCorrected = cms.bool(False)
puJetMva.jets = cms.InputTag("ak5PFJets")
puJetMva.applyJec = cms.bool(True)
puJetMva.inputIsCorrected = cms.bool(False)


#
# filters
#

from Smurf.ProcessingAndSkimming.filters_cfi import *

# sequences
leptonTreeMakerSequenceMC   = cms.Sequence(fastJetSequence * pfNoPileUpSequence * puJetIdSqeuence)
leptonTreeMakerSequenceData2011 = cms.Sequence(hltHighLevel2011 * fastJetSequence * pfNoPileUpSequence * puJetIdSqeuence)
leptonTreeMakerSequenceData2012 = cms.Sequence(hltHighLevel2012 * fastJetSequence * pfNoPileUpSequence * puJetIdSqeuence)

# electrons and muons
electronFilters = cms.Sequence(highPtElectrons * highPtElectronFilter)
muonFilters = cms.Sequence(highPtMuons * highPtMuonFilter)

# filter sequences
# photons
photonFilters = cms.Sequence(highPtPhotons * highPtPhotonFilter * centralPhotons * centralPhotonFilter)
leptonTreeMakerSequenceForPhotonMC   = cms.Sequence(photonFilters * leptonTreeMakerSequenceMC)
leptonTreeMakerSequenceForPhotonData = cms.Sequence(photonFilters * leptonTreeMakerSequenceData2011)

#
# random numbers
#

RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    leptonTreeMaker2012 = cms.PSet(
        initialSeed = cms.untracked.uint32(764575546),
         engineName = cms.untracked.string('TRandom3')
    ),
    leptonTreeMaker2011 = cms.PSet(
        initialSeed = cms.untracked.uint32(764575546),
         engineName = cms.untracked.string('TRandom3')
    )
)

