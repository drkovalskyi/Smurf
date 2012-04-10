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
leptonTreeMaker2012.electronTPTriggerNames = cms.untracked.vstring(
        "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*")
leptonTreeMaker2012.muonTPTriggerNames     = cms.untracked.vstring(
        "HLT_IsoMu17_v*")
leptonTreeMaker2012.electronFRTriggerNames = cms.untracked.vstring(        
        "HLT_Ele8_CaloIdL_CaloIsoVL_v*",
        "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*",
        "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        "HLT_Ele8_CaloIdT_TrkIdVL_v*")
leptonTreeMaker2012.muonFRTriggerNames     = cms.untracked.vstring(
        "HLT_Mu8_v*")
leptonTreeMaker2012.photonTriggerNames     = cms.untracked.vstring(
        "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v*",
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*",
        "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*",
        "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*")

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

