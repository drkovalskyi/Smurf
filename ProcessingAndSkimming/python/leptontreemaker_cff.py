import FWCore.ParameterSet.Config as cms

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoJets.Configuration.RecoPFJets_cff import * # Import the Jet RECO modules
kt6PFJets.doRhoFastjet  = True                # Turn-on the FastJet density calculation

#
# rho
#

from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJetsDeterministicIso = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
kt6PFJetsDeterministicIso.Rho_EtaMax = cms.double(2.5)

kt6PFJetsDeterministicJEC = kt4PFJets.clone(
    rParam = 0.6,
    doAreaFastjet = True,
    doRhoFastjet = True,
    voronoiRfact = 0.9,
    Rho_EtaMax = 5.0,
    Ghost_EtaMax = 5.0
)

fastJetSequence = cms.Sequence(kt6PFJets * kt6PFJetsDeterministicIso * kt6PFJetsDeterministicJEC)

#
# lepton maker
#

from Smurf.ProcessingAndSkimming.leptontreemaker_cfi import *

# For 2011
leptonTreeMaker2011 = leptonTreeMaker.clone()

# For 2012
leptonTreeMaker2012 = leptonTreeMaker.clone()
leptonTreeMaker2012.electronTPTriggerNames = cms.untracked.VInputTag(
        cms.InputTag("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"),
        cms.InputTag("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter"),
        cms.InputTag("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"),
        cms.InputTag("HLT_Ele27_WP80_v*"))
leptonTreeMaker2012.muonTPTriggerNames     = cms.untracked.VInputTag(
        cms.InputTag("HLT_IsoMu17_Mu8_v*"),
        cms.InputTag("HLT_IsoMu24_eta2p1_v*"))
leptonTreeMaker2012.electronFRTriggerNames = cms.untracked.VInputTag(       
        cms.InputTag("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*"),
        cms.InputTag("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"),
        cms.InputTag("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*"),
        cms.InputTag("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"))
leptonTreeMaker2012.muonFRTriggerNames     = cms.untracked.VInputTag(
        cms.InputTag("HLT_Mu8_v*"),
        cms.InputTag("HLT_Mu17_v*"))
leptonTreeMaker2012.photonTriggerNames     = cms.untracked.VInputTag(
        cms.InputTag("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v*"),
        cms.InputTag("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*"),
        cms.InputTag("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*"),
        cms.InputTag("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*"),
        cms.InputTag("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*"))
leptonTreeMaker2012.measureSingleEle    = cms.untracked.VInputTag(
        cms.InputTag('HLT_Ele27_WP80_v*'))
leptonTreeMaker2012.measureLeadingDoubleEle = cms.untracked.VInputTag(
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter'))
leptonTreeMaker2012.measureTrailingDoubleEle = cms.untracked.VInputTag(
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter'))
leptonTreeMaker2012.measureDoubleEleDZ = cms.untracked.VInputTag(
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ'))
leptonTreeMaker2012.measureSingleMu24       = cms.untracked.VInputTag(
        cms.InputTag('HLT_IsoMu24_eta2p1_v*'))
leptonTreeMaker2012.measureSingleMu30       = cms.untracked.VInputTag(
        cms.InputTag('HLT_IsoMu30_eta2p1_v*'))
leptonTreeMaker2012.measureLeadingDoubleMu  = cms.untracked.VInputTag(
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17'))
leptonTreeMaker2012.measureTrailingDoubleMu = cms.untracked.VInputTag(
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8'))
leptonTreeMaker2012.measureDoubleMuDZ       = cms.untracked.VInputTag(
        cms.InputTag('HLT_Mu17_Mu8_v*:hltDiMuonMu17Mu8DzFiltered0p2'))


#
# filters
#

from Smurf.ProcessingAndSkimming.filters_cfi import *

# sequences
leptonTreeMakerSequenceMC   = cms.Sequence(fastJetSequence)
leptonTreeMakerSequenceData2011 = cms.Sequence(hltHighLevel2011 * fastJetSequence)
leptonTreeMakerSequenceData2012 = cms.Sequence(hltHighLevel2012 * fastJetSequence)

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
    )
)

