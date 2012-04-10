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

# lepton maker
from Smurf.ProcessingAndSkimming.leptontreemaker_cfi import *
# filters
from Smurf.ProcessingAndSkimming.filters_cfi import *

# sequences
leptonTreeMakerSequenceMC   = cms.Sequence(fastJetSequence * leptonTreeMaker)
leptonTreeMakerSequenceData = cms.Sequence(hltHighLevel * leptonTreeMakerSequenceMC)

# filter sequences
# photons
photonFilters = cms.Sequence(highPtPhotons * highPtPhotonFilter * centralPhotons * centralPhotonFilter)
leptonTreeMakerSequenceForPhotonMC   = cms.Sequence(photonFilters * leptonTreeMakerSequenceMC)
leptonTreeMakerSequenceForPhotonData = cms.Sequence(photonFilters * leptonTreeMakerSequenceData)

# electrons
electronFilters = cms.Sequence(highPtElectrons * highPtElectronFilter)
leptonTreeMakerSequenceForElectronMC   = cms.Sequence(electronFilters * leptonTreeMakerSequenceMC)
leptonTreeMakerSequenceForElectronData = cms.Sequence(electronFilters * leptonTreeMakerSequenceData)

# muons
muonFilters = cms.Sequence(highPtMuons * highPtMuonFilter)
leptonTreeMakerSequenceForMuonMC   = cms.Sequence(muonFilters * leptonTreeMakerSequenceMC)
leptonTreeMakerSequenceForMuonData = cms.Sequence(muonFilters * leptonTreeMakerSequenceData)

