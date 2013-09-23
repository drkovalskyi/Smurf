// -*- C++ -*-
// $Id$

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// #include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"

#include "Smurf/ProcessingAndSkimming/interface/SmurfTreeProducer.h"
#include "Smurf/ProcessingAndSkimming/interface/Selections.h"
#include "TFile.h"
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

namespace Smurf {
  bool isGoodVertex(const reco::Vertex& vertex){
    if (!vertex.isValid() || vertex.isFake()) return false;
    if (vertex.ndof() <= 4.) return false;
    if (vertex.position().Rho() > 2.0) return false;
    if (fabs(vertex.position().Z()) > 24.0) return false;
    return true;
  }

  const reco::Vertex* SmurfTreeProducer::primaryVertex( edm::Event& iEvent) const
  {
    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByLabel(m_primaryVertex, vtxs);
    if (vtxs->empty()) return 0;
    if (isGoodVertex(vtxs->front()))  return &vtxs->front();
    return 0;
  }
 
  bool SmurfTreeProducer::passedElectronSelection( edm::Event& iEvent, const reco::GsfElectron& ele, ele::SelectionType selection ) const
  {
    const reco::Vertex* pv = primaryVertex(iEvent);
    if (!pv) return false;
    double pt = ele.pt();
    double d0 = fabs(ele.gsfTrack()->dxy(pv->position()));
    double dz = fabs(ele.gsfTrack()->dz(pv->position()));

    if ( selection==ele::All || selection==ele::Acceptance ){
      if (fabs(ele.eta()) > 2.5)    return false;
    }
    if ( selection==ele::All || selection==ele::d0 )
      if (d0 > 0.02) return false;
    if ( selection==ele::All || selection==ele::dz )
      if (dz > 0.1) return false;
    if ( selection==ele::All || selection==ele::Conversion ){
      edm::Handle<reco::ConversionCollection> conversions;
      iEvent.getByLabel(m_conversionsInputTag, conversions);
      edm::Handle<reco::BeamSpot> beamspot;
      iEvent.getByLabel(m_beamSpotInputTag, beamspot);
      unsigned int mhits = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      bool conv = ConversionTools::hasMatchedConversion(ele, conversions, beamspot->position());
      if (mhits > 0)      return false;
      if (conv)           return false;
    }
    if ( selection==ele::All || selection==ele::Iso ){
      edm::Handle<reco::PFCandidateCollection> pfs;
      iEvent.getByLabel(m_pfCands, pfs);
      edm::Handle<reco::VertexCollection> vtxs;
      iEvent.getByLabel(m_primaryVertex, vtxs);
      edm::Handle<double> rho;
      iEvent.getByLabel(m_rhoIsoAllInputTag, rho);
      // std::cout << "rho: " << *rho << std::endl;
      float chiso     = 0.0; float emiso = 0.0; float nhiso = 0.0; float dbeta = 0.0;
      double ea04      = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04,
									 ele.eta(), ElectronEffectiveArea::kEleEAData2012);
      // std::cout << "ea04: " << ea04 << std::endl;
      smurfselections::PFIsolation2012(ele, *pfs, *vtxs, m_pfPileUpAlgo, 0, 0.4,
				       chiso, emiso, nhiso, dbeta, false, true, true, false);
      float rhoPrime  = std::max(0.0, *rho);
      float neutral   = std::max(0.0, emiso + nhiso - ea04 * rhoPrime);
      // std::cout << "neutral: " << neutral << std::endl;
      // std::cout << "chiso + neutral)/ele.pt(): " << (chiso + neutral)/ele.pt() << std::endl;
      if ((chiso + neutral)/ele.pt() > 0.15) return false;
    }
    if ( selection==ele::All || selection==ele::Id ){
      float eta = fabs(ele.superCluster()->eta());
      float mvaValue  = m_reader_egammaPOG2012MVA->mvaValue(ele, *pv, *m_ttBuilder, *m_clusterTools, false);

      if (ele.pt() > 20.0) {
        if (eta < 0.8    && mvaValue < 0.94) return false;
	if (eta >= 0.8   && eta < 1.479 && mvaValue < 0.85) return false;
        if (eta >= 1.479 && mvaValue < 0.92) return false;
      } else {
        if (eta < 0.8    && mvaValue < 0.00) return false;
	if (eta >= 0.8   && eta < 1.479 && mvaValue < 0.10) return false;
        if (eta >= 1.479 && mvaValue < 0.62) return false;
      }
    }
    if ( selection==ele::All || selection==ele::FO ){
//     if (fabs(etaSC)<1.479) {
// 	if (cms2.els_sigmaIEtaIEta().at(i)>0.01 ||
// 	    fabs(cms2.els_dEtaIn().at(i))>0.007 ||
// 	    fabs(cms2.els_dPhiIn().at(i))>0.15 ||
// 	    cms2.els_hOverE().at(i)>0.12 ||
// 	    cms2.els_tkIso().at(i)/pt>0.2 ||
// 	    TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
// 	    //cms2.els_ecalIso().at(i)/pt>0.20 ||////FIXME                                                                                                                            
// 	    cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
//       } else {
// 	if (cms2.els_sigmaIEtaIEta().at(i)>0.03 ||
// 	    fabs(cms2.els_dEtaIn().at(i))>0.009 ||
// 	    fabs(cms2.els_dPhiIn().at(i))>0.10 ||
// 	    cms2.els_hOverE().at(i)>0.10 ||
// 	    cms2.els_tkIso().at(i)/pt>0.2 ||
// 	    //TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||                                                                                                              
// 	    cms2.els_ecalIso().at(i)/pt>0.20 ||
// 	    cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
//       }

      
      if (ele.dr03TkSumPt()/pt          > 0.2)  return false;
      if (ele.dr03HcalTowerSumEt()/pt   > 0.2)  return false;
      if (ele.isEB()) {
        if (ele.sigmaIetaIeta()                               > 0.01)  return false;
        if (fabs(ele.deltaEtaSuperClusterTrackAtVtx())        > 0.007) return false;
        if (fabs(ele.deltaPhiSuperClusterTrackAtVtx())        > 0.15)  return false;
        if (ele.hadronicOverEm()                              > 0.12)  return false;
        if (std::max(ele.dr03EcalRecHitSumEt() - 1.0, 0.0)/pt > 0.2)   return false;
      } else {
        if (ele.sigmaIetaIeta()                               > 0.03)  return false;
        if (fabs(ele.deltaEtaSuperClusterTrackAtVtx())        > 0.009) return false;
        if (fabs(ele.deltaPhiSuperClusterTrackAtVtx())        > 0.10)  return false;
        if (ele.hadronicOverEm()                              > 0.10)  return false;
        if (ele.dr03EcalRecHitSumEt()/pt                      > 0.2)   return false;
      }
    }
    return true;
  }

  bool SmurfTreeProducer::passedMuonSelection( edm::Event& iEvent, const reco::Muon& muon, muon::SelectionType selection ) const
  {
    const reco::Vertex* pv = primaryVertex(iEvent);
    if (!pv) return false;

    const reco::TrackRef siTrack = muon.innerTrack();
    double d0 = siTrack.isNonnull() ? muon.innerTrack()->dxy(pv->position())    : 999.9;
    double dz = siTrack.isNonnull() ? muon.innerTrack()->dz(pv->position())     : 999.9;
    unsigned int nLayers = siTrack.isNonnull()      ?  muon.track()->hitPattern().trackerLayersWithMeasurement() : 0;
    unsigned int nValidPixelHits = siTrack.isNonnull() ?  muon.innerTrack()->hitPattern().numberOfValidPixelHits() : 0;
    double ptErr = siTrack.isNonnull()      ?  muon.innerTrack()->ptError() : 99999.9;
    double kink  = muon.combinedQuality().trkKink;
    
    if ( selection==muon::All || selection==muon::Acceptance ){
      if (fabs(muon.eta()) > 2.4)    return false;
      if (!muon.isTrackerMuon()&&!muon.isGlobalMuon()) return false;
    }
    if ( selection==muon::All || selection==muon::d0 ){
      if (muon.pt() < 20.0 && fabs(d0) >= 0.01)    return false;           
      else if (fabs(d0) >= 0.02)                   return false;
    }
    if ( selection==muon::All || selection==muon::dz )
      if (fabs(dz) >= 0.1)           return false;
    if ( selection==muon::All || selection==muon::Id ){
      if (nLayers < 6)               return false;
      if (nValidPixelHits == 0)      return false;
      if (kink > 20.0)                return false;
      if (!muon.isPFMuon())          return false;   
      bool goodTightMuon = true;
      if (!muon.isGlobalMuon() || !muon.isTrackerMuon()) goodTightMuon = false;
      if (goodTightMuon) {
        const reco::TrackRef globalTrack = muon.globalTrack();
        double nValidMuonHits = globalTrack.isNonnull()   ? globalTrack->hitPattern().numberOfValidMuonHits() : 0;
        double chi2ndof       = globalTrack.isNonnull()   ? globalTrack->chi2()/globalTrack->ndof()           : 999.9;
        if (muon.numberOfMatches() < 2)    goodTightMuon = false;
        if (chi2ndof >= 10)                 goodTightMuon = false;
        if (nValidMuonHits == 0)            goodTightMuon = false;
      }
      bool goodTrackerMuon = true;
      if (muon.isTrackerMuon()) {
        if (!::muon::isGoodMuon(muon, ::muon::TMLastStationTight)) goodTrackerMuon = false;
      } else goodTrackerMuon = false;                 
      if (!goodTrackerMuon && !goodTightMuon)   return false;
      if (ptErr / muon.innerTrack()->pt() > 0.1)   return false;
    }
    if ( selection==muon::All || selection==muon::Iso ){
      static int counter(0);
      edm::Handle<reco::PFCandidateCollection> pfs;
      iEvent.getByLabel(m_pfCands, pfs);
      edm::Handle<double> rho;
      iEvent.getByLabel(m_rhoIsoAllInputTag, rho);
      // iEvent.getByLabel(m_rhoIsoNeutral2011InputTag, rho);
      const reco::GsfElectronCollection nullEls;
      const reco::MuonCollection nullMus;
      
      double mvaValue = 
	m_reader_muonHZZ2012IsoRingsMVA->mvaValue(muon, *pv, *pfs, *rho, 
						  MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
      if ( counter>0 ){
	std::cout << "Run : " << iEvent.id().run() << " \tEvent : " << iEvent.id().event() << "Rho : " << *rho << "\n" 
		  << "Pt/Eta : " << muon.pt() << " / " << muon.eta() << " \tMVA : " << mvaValue << std::endl;
	counter--;
      }
      double eta = fabs(muon.eta());
      if (muon.pt() > 20.0) {
        if (eta < 1.479) {
	  if (mvaValue < 0.82)   return false;
	} else {
	  if (mvaValue < 0.86)   return false;
	}
      } else {        
        if (eta < 1.479) {
	  if (mvaValue < 0.86)   return false;
	} else {
	  if (mvaValue < 0.82)   return false;
	}
      }               
    }
    return true;
  }

// bool SmurfTreeProducer::muAccept( const reco::Muon& muon ) const
// {
//   return fabs(muon.eta())<2.4;
// }
// bool SmurfTreeProducer::elAccept( const reco::GsfElectron& el ) const
// {
//   return fabs(el.eta())<2.5;
// }

// bool SmurfTreeProducer::muonid(  const reco::Muon& muon ) const
// {
//   if ( ! muon.isGlobalMuon() ) return false;
//   if ( ! muon.isTrackerMuon() ) return false;
//   if ( muon.globalTrack()->normalizedChi2()>10 ) return false;
//   if ( muon.globalTrack()->hitPattern().numberOfValidMuonHits() == 0 ) return false;
//   if ( muon.innerTrack()->numberOfValidHits() < 11 ) return false;
//   return true;
// }
// double SmurfTreeProducer::muonisoValue(  const reco::Muon& muon ) const
// {
//   double sum = muon.isolationR03().sumPt + muon.isolationR03().emEt + muon.isolationR03().hadEt;
//   return sum/muon.pt();
// }
// bool SmurfTreeProducer::muoniso(  const reco::Muon& muon ) const
// {
//   return muonisoValue(muon) < 0.15;
// }

// bool SmurfTreeProducer::vbtfId( const edm::Event& iEvent, 
// 				 std::string name, const reco::GsfElectronRef& el ) const
// {
//   bool passedTrivialSelector = false;
//   if ( name == "simpleEleId70relIso" ){
//     if (el->isEB())
//       passedTrivialSelector = el->sigmaIetaIeta()<0.01 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.03 &&
// 	fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.004 && el->hadronicOverEm()<0.025;
//     else
//       passedTrivialSelector = el->sigmaIetaIeta()<0.03 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.02 &&
// 	fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.005 && el->hadronicOverEm()<0.025;
//   }
//   if ( name == "simpleEleId80relIso" ){
//     if (el->isEB())
//       passedTrivialSelector = el->sigmaIetaIeta()<0.01 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.06 &&
// 	fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.004 && el->hadronicOverEm()<0.04;
//     else
//       passedTrivialSelector = el->sigmaIetaIeta()<0.03 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.03 &&
// 	fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.007 && el->hadronicOverEm()<0.025;
//   }
//   if ( name == "simpleEleId90relIso" ){
//     if (el->isEB())
//       passedTrivialSelector = el->sigmaIetaIeta()<0.01 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.8 &&
// 	fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.007 && el->hadronicOverEm()<0.12;
//     else
//       passedTrivialSelector = el->sigmaIetaIeta()<0.03 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.7 &&
// 	fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.009 && el->hadronicOverEm()<0.05;
//   }
//   return passedTrivialSelector;

//   /*
//   edm::Handle<edm::ValueMap<float> > map;
//   iEvent.getByLabel(name,map);
//   // get rid of floats to avoid rounding errors
//   unsigned int flag = abs(int((*map)[el]+0.01));
//   return (flag & 0x01);
//   */
//   /*
//   if (strcmp(name,"simpleEleId80relIso")==0){
//     edm::Handle<reco::TrackCollection> ctfTracks;
//     iEvent.getByLabel("generalTracks",ctfTracks);
//     SimpleCutBasedElectronIDSelectionFunctor selector(SimpleCutBasedElectronIDSelectionFunctor::relIso95, 3.8112, ctfTracks);
//     pat::strbitset ret;
//     selector(pat::Electron(*el),ret);
//     if ( ret.test("sihih_EB") && ret.test("dphi_EB") && ret.test("deta_EB") && ret.test("hoe_EB") &&
// 	 ret.test("sihih_EE") && ret.test("dphi_EE") && ret.test("deta_EE") && ret.test("hoe_EE") )
//       return true;
//   }
//   return false;
//   */
//   /*
//   if (el->isEB())
//     return el->sigmaIetaIeta()<0.01 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.06 &&
//       fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.004 && el->hadronicOverEm()<0.04;
//   else
//     return el->sigmaIetaIeta()<0.03 && fabs(el->deltaPhiSuperClusterTrackAtVtx()<0.03) &&
//       fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.007 && el->hadronicOverEm()<0.025;
//   */
// }

// double SmurfTreeProducer::elisoValue( const reco::GsfElectron& el) const
// {
//   double sum = el.dr03TkSumPt() + std::max(el.dr03EcalRecHitSumEt()-1.0,0.0) + 
//     el.dr03HcalTowerSumEt();
//   return sum/el.pt();
// }
// bool SmurfTreeProducer::eliso( const reco::GsfElectron& el) const
// {
//   return elisoValue(el) < 0.1;
// }

// bool SmurfTreeProducer::convRejection_partner( const edm::Event& iEvent, const reco::GsfElectron& el ) const
// {
//   edm::Handle<reco::TrackCollection> tracks;
//   iEvent.getByLabel("generalTracks", tracks);
//   // double bfield = 3.8112;
//   // ConversionFinder convFinder;
//   // return !ConversionFinder::isElFromConversion(el, tracks, bfield, 0.02, 0.02, 0.45);
//   ConversionFinder convFinder;
//   ConversionInfo convInfo = convFinder.getConversionInfo(el, tracks, 3.8112);
//   return fabs(convInfo.dist()) > 0.02 || fabs(convInfo.dcot()) > 0.02;
// }

// bool SmurfTreeProducer::convRejection_hits( const reco::GsfElectron& el ) const
// {
//   return el.gsfTrack()->trackerExpectedHitsInner().numberOfHits() == 0;
// }


// bool SmurfTreeProducer::metPreselection( const WW& ww ) const
// {
//   return ww.met>20;
// }

// bool SmurfTreeProducer::metFull( const WW& ww ) const
// {
//   double tightDPhi = fabs(deltaPhi(ww.lt_p4.phi(), ww.metPhi)); 
//   double looseDPhi = fabs(deltaPhi(ww.ll_p4.phi(), ww.metPhi)); 
//   double deltaPhi = std::min(tightDPhi, looseDPhi);

//   double pMet = ww.met;

//   if (deltaPhi < TMath::Pi()/2) pMet = pMet*TMath::Sin(deltaPhi);

//   if ( pMet < 20 ) return false;
//   if (ww.type() == WW::EE || ww.type() == WW::MM) {
//     if ( ww.met < 45 ) return false;
//     // if ( !metBalance(i_hyp) ) return false;
//   }
//   return true;
// }

// bool SmurfTreeProducer::metBad( const WW& ww ) const
// {
//   double tightDPhi = deltaPhi(ww.lt_p4.phi(), ww.metPhi); 
//   double looseDPhi = deltaPhi(ww.ll_p4.phi(), ww.metPhi); 
//   double deltaPhi = std::min(tightDPhi, looseDPhi);

//   double pMet = ww.met;

//   if (deltaPhi < TMath::Pi()/2) pMet = pMet*TMath::Sin(deltaPhi);

//   if ( pMet < 20 ) return false;
//   if (ww.type() == WW::EE || ww.type() == WW::MM) {
//     if ( pMet < 45 ) return false;
//     // if ( !metBalance(i_hyp) ) return false;
//   }
//   return true;
// }

// template<typename T> T
// SmurfTreeProducer::mostEnergeticJet( const std::vector<T>& collection, const WW& ww,
// 				      double maxEta ) const
// {
//   T ojet;
//   for ( typename std::vector<T>::const_iterator jet = collection.begin();
// 	jet != collection.end(); ++jet ){
//     if ( reco::deltaR(ww.lt_p4,*jet) < 0.3 ) continue;
//     if ( reco::deltaR(ww.ll_p4,*jet) < 0.3 ) continue;
//     if ( jet->pt() < ojet.pt() || fabs(jet->eta()) > maxEta ) continue;
//     ojet = *jet;
//   }
//   return ojet;
// }

// template<typename T> bool
// SmurfTreeProducer::jetVeto( const T& collection, const WW& ww,
// 			     double maxPt, double maxEta ) const
// {
//   return mostEnergeticJet(collection,ww,maxEta).pt() < maxPt;
// }


// bool SmurfTreeProducer::chargeMatched( const reco::GsfElectron& el ) const
// {
//   reco::TrackRef ctf = el.closestCtfTrackRef();
//   if ( !ctf.isNull() && ctf->charge()*el.charge() < 0 ) return false;
//   return true;
// }
  
  SmurfTreeProducer::SmurfTreeProducer(const edm::ParameterSet& iConfig):
    m_muons(              iConfig.getParameter<edm::InputTag>("inputMuonCollection" ) ),
    m_electrons(          iConfig.getParameter<edm::InputTag>("inputElectronCollection" ) ),
    m_jets(               iConfig.getParameter<edm::InputTag>("inputJetCollection" ) ),
    m_primaryVertex(      iConfig.getParameter<edm::InputTag>("inputPVCollection" ) ),
    m_pfCands(            iConfig.getParameter<edm::InputTag>("inputPFCandCollection" ) ),
    m_rhoIsoAllInputTag(  iConfig.getParameter<edm::InputTag>("rhoIsoAll" ) ),
    m_rhoIsoNeutral2011InputTag(  iConfig.getParameter<edm::InputTag>("rhoIsoNeutral2011" ) ),
    m_conversionsInputTag(iConfig.getParameter<edm::InputTag>("conversions" ) ),
    m_beamSpotInputTag(   iConfig.getParameter<edm::InputTag>("beamspot" ) ),
    m_monitor(true)
  {
    m_file = TFile::Open(iConfig.getParameter<std::string>("smurfTree").c_str(), "RECREATE");
    m_tree = new SmurfTree();
    m_tree->CreateTree();
    m_tree->tree_->SetDirectory(m_file);
    m_l1MinPt = 20;
    m_l2MinPt = 10;
    m_matchCharge = true;

    std::string cmssw_base = getenv("CMSSW_BASE");
    // Iso rings
    std::vector<std::string> muonHZZ2012IsoRingsWeights;
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
    m_reader_muonHZZ2012IsoRingsMVA = new MuonMVAEstimator();
    m_reader_muonHZZ2012IsoRingsMVA->initialize("muonHZZ2012IsoRingsMVA", 
						MuonMVAEstimator::kIsoRings, true, muonHZZ2012IsoRingsWeights);
    m_pfPileUpAlgo = new PFPileUpAlgo();
    
    m_reader_egammaPOG2012MVA = new EGammaMvaEleEstimator();
    std::vector<std::string> myManualCatWeigthsTrig;
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat1.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat2.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat3.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat4.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat5.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat6.weights.xml");
    m_reader_egammaPOG2012MVA->initialize("BDT", EGammaMvaEleEstimator::kTrig, true, myManualCatWeigthsTrig);

  }
  
  SmurfTreeProducer::~SmurfTreeProducer(){}

  void SmurfTreeProducer::endJob(){ m_monitor.print(); }



  void SmurfTreeProducer::fillCommonVariables(const edm::Event& iEvent)
  {
    
    m_tree->run_       = iEvent.id().run()        ;
    m_tree->event_     = iEvent.id().event()      ;
    m_tree->lumi_      = iEvent.luminosityBlock() ;
    // m_tree->nvtx_      = smurfselections::CountGoodPV(vtx_h_);
    // m_tree->scale1fb_  = 1.0;
  }

  
  bool SmurfTreeProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    // if (iEvent.id().event()!=1764047) return false; //REMOVE ME
    edm::Handle<reco::MuonCollection> muons;
    iEvent.getByLabel(m_muons,muons);
    edm::Handle<reco::GsfElectronCollection> els;
    iEvent.getByLabel(m_electrons, els);
    edm::Handle<reco::PFJetCollection> jets;
    iEvent.getByLabel(m_jets,jets);

    // unbiased revertexing
    edm::ESHandle<TransientTrackBuilder> ttBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
    m_ttBuilder = ttBuilder.product();
    // ecal tools
    EcalClusterLazyTools clusterTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));
    m_clusterTools = &clusterTools;
    
    m_tree->InitVariables();
    fillCommonVariables(iEvent);
  
    m_monitor.count(*m_tree, SmurfTree::mm, "All events");

    std::vector<std::pair<unsigned int, const reco::Candidate*> > leptons;
    for (unsigned int iMuon = 0; iMuon < muons->size(); ++iMuon){
      if (muons->at(iMuon).pt()<m_l1MinPt && muons->at(iMuon).pt()<m_l2MinPt) continue;
      leptons.push_back( std::pair<unsigned int, const reco::Candidate*>(iMuon,&muons->at(iMuon)) );
    }
    for (unsigned int iEl = 0; iEl < els->size(); ++iEl){
      if (els->at(iEl).pt()<m_l1MinPt && els->at(iEl).pt()<m_l2MinPt) continue;
      leptons.push_back( std::pair<unsigned int, const reco::Candidate*>(iEl,&els->at(iEl)) );
    }
    if ( leptons.size() < 2 ){
      return false;
    }
    m_monitor.count(*m_tree, SmurfTree::mm, "two leptons");
    // edm::Handle<reco::METCollection> mets;
    // iEvent.getByLabel("tcMet",mets);

    for (unsigned int iLep1 = 0; iLep1 < leptons.size()-1; ++iLep1)
      for (unsigned int iLep2 = iLep1+1; iLep2 < leptons.size(); ++iLep2){
	if ( (leptons.at(iLep1).second->pt()<m_l1MinPt || leptons.at(iLep2).second->pt()<m_l2MinPt) &&
	     (leptons.at(iLep1).second->pt()<m_l2MinPt || leptons.at(iLep2).second->pt()<m_l1MinPt) ) continue;
	unsigned int ilt = leptons.at(iLep1).first;
	unsigned int ill = leptons.at(iLep2).first;
	const reco::Candidate* lt = leptons.at(iLep1).second;
	const reco::Candidate* ll = leptons.at(iLep2).second;
	if (lt->pt()<ll->pt()){
	  ill = leptons.at(iLep1).first;
	  ilt = leptons.at(iLep2).first;
	  ll = leptons.at(iLep1).second;
	  lt = leptons.at(iLep2).second;
	}
	SmurfTree::Type type = SmurfTree::mm;
	if ( abs(lt->pdgId()) == 13 && abs(ll->pdgId()) == 11 ) type = SmurfTree::me;
	if ( abs(lt->pdgId()) == 11 && abs(ll->pdgId()) == 11 ) type = SmurfTree::ee;
	if ( abs(lt->pdgId()) == 11 && abs(ll->pdgId()) == 13 ) type = SmurfTree::em;
	
	// if (type != SmurfTree::ee) continue; // REMOVE ME
	m_monitor.count(*m_tree, type, "pt 20/10");
	
	if ( m_matchCharge && leptons.at(iLep1).second->charge()*leptons.at(iLep2).second->charge() > 0 ) continue;
	
	m_monitor.count(*m_tree, type, "charge + pt 20/10");

	if ( abs(lt->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)lt, muon::Acceptance)) continue;
	if ( abs(ll->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)ll, muon::Acceptance)) continue;
	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::Acceptance)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::Acceptance)) continue;

	m_monitor.count(*m_tree, type, "charge + pt 20/10 + acceptance");

	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::Iso)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::Iso)) continue;

	m_monitor.count(*m_tree, type, "+ele isolation");

	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::Conversion)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::Conversion)) continue;

	m_monitor.count(*m_tree, type, "+conversion");

	if ( abs(lt->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)lt, muon::d0)) continue;
	if ( abs(ll->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)ll, muon::d0)) continue;
	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::d0)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::d0)) continue;

	m_monitor.count(*m_tree, type, "+ d0");

	if ( abs(lt->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)lt, muon::dz)) continue;
	if ( abs(ll->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)ll, muon::dz)) continue;
	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::dz)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::dz)) continue;

	m_monitor.count(*m_tree, type, "+ dz");

	if ( abs(lt->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)lt, muon::Id)) continue;
	if ( abs(ll->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)ll, muon::Id)) continue;
	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::Id)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::Id)) continue;

	m_monitor.count(*m_tree, type, "+ id");

	if ( abs(lt->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)lt, muon::All)) continue;
	if ( abs(ll->pdgId()) == 13 && !passedMuonSelection(iEvent,*(const reco::Muon*)ll, muon::All)) continue;
	if ( abs(lt->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)lt, ele::All)) continue;
	if ( abs(ll->pdgId()) == 11 && !passedElectronSelection(iEvent,*(const reco::GsfElectron*)ll, ele::All)) continue;
	
	m_monitor.count(*m_tree, type, "lepton id/iso pt20/10");

//       m_tree->InitVariables();
//       fillCommonVariables(iEvent);
//       ww.metPhi = mets->front().phi();
//       ww.met    = mets->front().pt();

//       if ( abs(lt->pdgId()) == 13 ){
// 	ww.lt_mu = reco::MuonRefVector::value_type(muons,ilt);
// 	ww.lt_relIso = muonisoValue(*ww.lt_mu);
//       }
//       if ( abs(ll->pdgId()) == 13 ){
// 	ww.ll_mu = reco::MuonRefVector::value_type(muons,ill);
// 	ww.ll_relIso = muonisoValue(*ww.ll_mu);
//       }
//       if ( abs(lt->pdgId()) == 11 ){
// 	ww.lt_el = reco::GsfElectronRefVector::value_type(els,ilt);
// 	ww.lt_relIso = elisoValue(*ww.lt_el);
//       }
//       if ( abs(ll->pdgId()) == 11 ){
// 	ww.ll_el = reco::GsfElectronRefVector::value_type(els,ill);
// 	ww.ll_relIso = elisoValue(*ww.ll_el);
//       }
//       ww.lt_p4 = lt->p4();
//       ww.lt_id = lt->pdgId();
//       ww.ll_p4 = ll->p4();
//       ww.ll_id = ll->pdgId();
//       ww.cuts.setType(lt->pdgId(),ll->pdgId());

//       ww.p4 = ww.lt_p4+ww.ll_p4;
//       if ( ww.p4.mass() < m_minMass ) continue;
      
//       ww.lt_d0BS = dxy(lt,beamSpot->position());
//       ww.ll_d0BS = dxy(ll,beamSpot->position());
//       if ( m_primaryVertex ){
// 	ww.lt_d0PV = dxy(lt, m_primaryVertex->position());
// 	ww.ll_d0PV = dxy(ll, m_primaryVertex->position());
// 	ww.lt_dzPV = dz(lt, m_primaryVertex->position());
// 	ww.ll_dzPV = dz(ll, m_primaryVertex->position());
//       }
//       if ( !ww.lt_mu.isNull() ){
// 	if ( muonid( *ww.lt_mu ) )  ww.cuts.add_lt_mu_Cut( WW::MuonId );
// 	if ( muoniso( *ww.lt_mu ) ) ww.cuts.add_lt_mu_Cut( WW::MuonIso );
// 	if ( muond0( ww.lt_d0PV ) ) ww.cuts.add_lt_mu_Cut( WW::MuonImpPV );
// 	// if ( muond0( ww.lt_d0BS ) ) ww.cuts.add_lt_mu_Cut( WW::MuonImpBS );
// 	if ( muAccept(*ww.lt_mu) )  ww.cuts.add_lt_mu_Cut( WW::MuonAccept );
//       }
//       if ( !ww.ll_mu.isNull() ){
// 	if ( muonid( *ww.ll_mu ) )  ww.cuts.add_ll_mu_Cut( WW::MuonId );
// 	if ( muoniso( *ww.ll_mu ) ) ww.cuts.add_ll_mu_Cut( WW::MuonIso );
// 	if ( muond0( ww.ll_d0PV ) ) ww.cuts.add_ll_mu_Cut( WW::MuonImpPV );
// 	// if ( muond0( ww.ll_d0BS ) ) ww.cuts.add_ll_mu_Cut( WW::MuonImpBS );
// 	if ( muAccept(*ww.ll_mu) )  ww.cuts.add_ll_mu_Cut( WW::MuonAccept );
//       }
//       if ( !ww.lt_el.isNull() ){
// 	if ( vbtfId(iEvent,"simpleEleId70relIso",ww.lt_el) )  ww.cuts.add_lt_el_Cut( WW::EleVBTF70Id );
// 	if ( vbtfId(iEvent,"simpleEleId80relIso",ww.lt_el) )  ww.cuts.add_lt_el_Cut( WW::EleVBTF80Id );
// 	if ( vbtfId(iEvent,"simpleEleId90relIso",ww.lt_el) )  ww.cuts.add_lt_el_Cut( WW::EleVBTF90Id );
// 	if ( eliso( *ww.lt_el ) )                             ww.cuts.add_lt_el_Cut( WW::EleIso );
// 	if ( eld0( ww.lt_d0PV ) )                             ww.cuts.add_lt_el_Cut( WW::EleImpPV );
// 	// if ( eld0( ww.lt_d0BS ) )                             ww.cuts.add_lt_el_Cut( WW::EleImpBS );
// 	if ( convRejection_partner(iEvent, *ww.lt_el) )       ww.cuts.add_lt_el_Cut( WW::ConvPartner );
// 	if ( convRejection_hits(*ww.lt_el) )                  ww.cuts.add_lt_el_Cut( WW::ConvNoMiss );
// 	// if ( chargeMatched(*ww.lt_el) )                       ww.cuts.add_lt_el_Cut( WW::ElChargeMatch );
// 	if ( elAccept( *ww.lt_el ) )                          ww.cuts.add_lt_el_Cut( WW::EleAccept );
//       }
//       if ( !ww.ll_el.isNull() ){
// 	if ( vbtfId(iEvent,"simpleEleId70relIso",ww.ll_el) )  ww.cuts.add_ll_el_Cut( WW::EleVBTF70Id );
// 	if ( vbtfId(iEvent,"simpleEleId80relIso",ww.ll_el) )  ww.cuts.add_ll_el_Cut( WW::EleVBTF80Id );
// 	if ( vbtfId(iEvent,"simpleEleId90relIso",ww.ll_el) )  ww.cuts.add_ll_el_Cut( WW::EleVBTF90Id );
// 	if ( eliso( *ww.ll_el ) )                             ww.cuts.add_ll_el_Cut( WW::EleIso );
// 	if ( eld0( ww.ll_d0PV ) )                             ww.cuts.add_ll_el_Cut( WW::EleImpPV );
// 	// if ( eld0( ww.ll_d0BS ) )                             ww.cuts.add_ll_el_Cut( WW::EleImpBS );
// 	if ( convRejection_partner(iEvent, *ww.ll_el) )       ww.cuts.add_ll_el_Cut( WW::ConvPartner );
// 	if ( convRejection_hits(*ww.ll_el) )                  ww.cuts.add_ll_el_Cut( WW::ConvNoMiss );
// 	// if ( chargeMatched(*ww.ll_el) )                       ww.cuts.add_ll_el_Cut( WW::ElChargeMatch );
// 	if ( elAccept( *ww.ll_el ) )                          ww.cuts.add_ll_el_Cut( WW::EleAccept );
//       }

//       if ( metPreselection(ww) )                    ww.cuts.add_hypo_Cut( WW::METBase );
//       if ( metFull(ww) )                            ww.cuts.add_hypo_Cut( WW::METFull );
//       if ( (ww.lt_p4.pt()>20 && ww.ll_p4.pt()>10) ||
// 	   (ww.ll_p4.pt()>20 && ww.lt_p4.pt()>10) ) ww.cuts.add_hypo_Cut( WW::Pt2010 );
//       if ( ww.lt_p4.pt()>20 && ww.ll_p4.pt()>20 )   ww.cuts.add_hypo_Cut( WW::Pt2020 );
//       if ( ww.lt_id * ww.ll_id < 0 )                ww.cuts.add_hypo_Cut( WW::ChargeMatch );
//       if ( jetVeto(*pfjets, ww, 25, 5 ) )           ww.cuts.add_hypo_Cut( WW::JetVeto25_5 );
//       if ( jetVeto(*pfjets, ww, 20, 3 ) )           ww.cuts.add_hypo_Cut( WW::JetVeto20_3 );
//       if ( jetVeto(*pfjets, ww, 30, 5 ) )           ww.cuts.add_hypo_Cut( WW::JetVeto30_5 );
//       if ( softmuon(iEvent,ww) )                    ww.cuts.add_hypo_Cut( WW::TopSoftMuon );
//       if ( extralepton(iEvent,ww) )                 ww.cuts.add_hypo_Cut( WW::ExtraLepton );

//       ww.jet_p4  = mostEnergeticJet(*pfjets, ww, 5).p4();

//       prod->push_back(ww);
    }
  // bool keepEvent = !prod->empty();
  // iEvent.put(prod);
  return true;
}

// double SmurfTreeProducer::dxy(const reco::Candidate* cand, const math::XYZPoint& vertex) const
// {
//   if ( const reco::Muon* muon = dynamic_cast<const reco::Muon*>(cand) )
//     if ( !muon->innerTrack().isNull() )
//       return muon->innerTrack()->dxy(vertex);
//   if ( const reco::GsfElectron* el = dynamic_cast<const reco::GsfElectron*>(cand) )
//     return el->gsfTrack()->dxy(vertex);
//   return -9999.;
// }

// double SmurfTreeProducer::dz(const reco::Candidate* cand, const math::XYZPoint& vertex) const
// {
//   if ( const reco::Muon* muon = dynamic_cast<const reco::Muon*>(cand) )
//     if ( !muon->innerTrack().isNull() )
//       return muon->innerTrack()->dz(vertex);
//   if ( const reco::GsfElectron* el = dynamic_cast<const reco::GsfElectron*>(cand) )
//     return el->gsfTrack()->dz(vertex);
//   return -9999.;
// }

// bool 
// SmurfTreeProducer::softmuon( const edm::Event& iEvent, const WW& ww ) const
// {
//   edm::Handle<reco::MuonCollection> muons;
//   iEvent.getByLabel(m_muons,muons);
//   for ( reco::MuonCollection::const_iterator muon = muons->begin(); 
// 	muon != muons->end(); ++muon )
//     {
//       if ( muon->pt()<3 ) continue;
//       if ( ww.lt_mu.isNonnull() && (&*muon == ww.lt_mu.get()) ) continue;
//       if ( ww.ll_mu.isNonnull() && (&*muon == ww.ll_mu.get()) ) continue;
//       if ( ! muon->isTrackerMuon() ) continue;
//       if ( ! muon::isGoodMuon(*muon, muon::TMLastStationAngTight) ) continue;
//       if ( muon->innerTrack()->numberOfValidHits() < 11 ) continue;
//       if ( !m_primaryVertex || 
// 	   fabs(dxy(&*muon, m_primaryVertex->position()))>0.2 ) continue;
//       return false;
//     }
//   return true;
// }

// bool 
// SmurfTreeProducer::extralepton( const edm::Event& iEvent, const WW& ww ) const
// {
//   edm::Handle<reco::MuonCollection> muons;
//   iEvent.getByLabel(m_muons,muons);
//   for ( reco::MuonCollection::const_iterator muon = muons->begin(); 
// 	muon != muons->end(); ++muon )
//     {
//       if ( muon->pt()<10 ) continue;
//       if ( ww.lt_mu.isNonnull() && (&*muon == ww.lt_mu.get()) ) continue;
//       if ( ww.ll_mu.isNonnull() && (&*muon == ww.ll_mu.get()) ) continue;
//       if ( ! muonid(*muon) ) continue;
//       if ( ! m_primaryVertex || ! muond0(dxy(&*muon, m_primaryVertex->position())) ) continue;
//       if ( ! muoniso(*muon) ) continue;
//       return false;
//     }

//   edm::Handle<reco::GsfElectronCollection> els;
//   iEvent.getByLabel(m_electrons, els);
//   for ( unsigned int i=0; i<els->size(); ++i)
//     {
//       if ( els->at(i).pt()<10 ) continue;
//       if ( reco::deltaR(ww.lt_p4, els->at(i).p4()) < 0.1 ) continue;
//       if ( reco::deltaR(ww.ll_p4, els->at(i).p4()) < 0.1 ) continue;
//       if ( ! vbtfId(iEvent,"simpleEleId80relIso",
// 		    reco::GsfElectronRefVector::value_type(els,i)) ) continue;
//       if ( ! m_primaryVertex || ! eld0(dxy(&els->at(i), m_primaryVertex->position())) ) continue;
//       if ( ! eliso(els->at(i)) ) continue;
//       return false;
//     }
//   return true;
// }
}

//Define this as a plug-in
DEFINE_FWK_MODULE(Smurf::SmurfTreeProducer);
