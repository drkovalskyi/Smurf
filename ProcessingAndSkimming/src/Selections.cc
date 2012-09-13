#include "Smurf/ProcessingAndSkimming/interface/Selections.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include <algorithm>
#include <utility>

//
// utilities
//

bool smurfselections::compareJetPt(std::pair<reco::PFJet, float> lv1, std::pair<reco::PFJet, float> lv2) {
  return lv1.first.pt()*lv1.second > lv2.first.pt()*lv2.second;
}

std::vector<std::pair<reco::PFJet, float> > smurfselections::goodJets(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const edm::Handle<edm::View<reco::PFJet> > &jets_h, const reco::Candidate &cand1, const reco::Candidate &cand2,
	const JetCorrector *corrector, float ptCut)

{

    edm::Handle<edm::ValueMap<int> > puJetIdFlag;
    iEvent.getByLabel("puJetMva","fullId",puJetIdFlag);

    std::vector<std::pair<reco::PFJet, float> > goodJets;
    edm::View<reco::PFJet> jetsCollection = *(jets_h.product());
    edm::View<reco::PFJet>::const_iterator jet;
    for (jet = jetsCollection.begin(); jet != jetsCollection.end(); ++jet) {

        // jec
        int idx = jet - jetsCollection.begin();
        edm::RefToBase<reco::Jet> jetRef(edm::Ref<edm::View<reco::PFJet> >(jets_h, idx));

        #ifdef RELEASE_4XY
        float jec = corrector->correction(*jet, jetRef, iEvent, iSetup);
        #else
        float jec = corrector->correction(*jet, iEvent, iSetup);
        #endif

	int    idflag = (*puJetIdFlag)[jetRef];
	if (PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )==0) continue;

        // cuts
        if (jet->pt() * jec <= ptCut) continue;
        if (fabs(jet->eta()) >= 5.0) continue;
        if (reco::deltaR(jet->eta(), jet->phi(), cand1.eta(), cand1.phi()) < 0.3) continue;
        if (reco::deltaR(jet->eta(), jet->phi(), cand2.eta(), cand2.phi()) < 0.3) continue;

        // push back
        goodJets.push_back(std::make_pair(*jet, jec));

    }

    std::sort(goodJets.begin(), goodJets.end(), compareJetPt);
    return goodJets;
}

std::vector<std::pair<reco::PFJet, float> > smurfselections::goodJets(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const edm::Handle<edm::View<reco::PFJet> > &jets_h, const reco::Candidate &cand, 
	const JetCorrector *corrector, float ptCut)
{


    edm::Handle<edm::ValueMap<int> > puJetIdFlag;
    iEvent.getByLabel("puJetMva","fullId",puJetIdFlag);

    std::vector<std::pair<reco::PFJet, float> > goodJets;
    edm::View<reco::PFJet> jetsCollection = *(jets_h.product());
    edm::View<reco::PFJet>::const_iterator jet;
    for (jet = jetsCollection.begin(); jet != jetsCollection.end(); ++jet) {

        // jec
        int idx = jet - jetsCollection.begin();
        edm::RefToBase<reco::Jet> jetRef(edm::Ref<edm::View<reco::PFJet> >(jets_h, idx));

        #ifdef RELEASE_4XY
        float jec = corrector->correction(*jet, jetRef, iEvent, iSetup);
        #else
        float jec = corrector->correction(*jet, iEvent, iSetup);
        #endif

	int    idflag = (*puJetIdFlag)[jetRef];
	if (PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )==0) continue;

        // cuts
        if (jet->pt() * jec <= ptCut) continue;
        if (fabs(jet->eta()) >= 5.0) continue;
        if (reco::deltaR(jet->eta(), jet->phi(), cand.eta(), cand.phi()) < 0.3) continue;

        // push back
        goodJets.push_back(std::make_pair(*jet, jec));

    }

    std::sort(goodJets.begin(), goodJets.end(), compareJetPt);
    return goodJets;
}

float smurfselections::electronIsoValuePF(const reco::PFCandidateCollection &pfCandCollection,
        const reco::GsfElectron& el, const reco::Vertex& vtx, float coner, float minptn, float dzcut,
        float footprintdr, float gammastripveto, float elestripveto, int filterId)
{

    float pfciso = 0.;
    float pfniso = 0.;
    float pffootprint = 0.;
    float pfjurveto = 0.;
    float pfjurvetoq = 0.;

    reco::TrackRef siTrack  = el.closestCtfTrackRef();
    reco::GsfTrackRef gsfTrack = el.gsfTrack();

    if (gsfTrack.isNull() && siTrack.isNull()) return -9999.;

    float eldz = gsfTrack.isNonnull() ? gsfTrack->dz(vtx.position()) : siTrack->dz(vtx.position());
    float eleta = el.eta();

    reco::PFCandidateCollection::const_iterator pf;
    for (pf = pfCandCollection.begin(); pf != pfCandCollection.end(); ++pf) {

        float pfeta = pf->eta();    
        float dR = reco::deltaR(pfeta, pf->phi(), eleta, el.phi());
        if (dR>coner) continue;

        float deta = fabs(pfeta - eleta);
        int pfid = abs(pf->pdgId());
        float pfpt = pf->pt();

        if (filterId!=0 && filterId!=pfid) continue;

        if (pf->charge()==0) {
            //neutrals
            if (pfpt>minptn) {
                pfniso+=pfpt;
                if (dR<footprintdr && pfid==130) pffootprint+=pfpt;
                if (deta<gammastripveto && pfid==22)  pfjurveto+=pfpt;
            }
        } else {
            //charged  
            //avoid double counting of electron itself
            //if either the gsf or the ctf track are shared with the candidate, skip it
            const reco::TrackRef pfTrack  = pf->trackRef();
            if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
            //if (pfid==11 && pf->gsfTrackRef().isNonnull()) {  // SYNC 18/04/12
            if (pf->gsfTrackRef().isNonnull()) {                // SYNC 18/04/12
                if (gsfTrack.isNonnull() && gsfTrack.key()==pf->gsfTrackRef().key()) continue;
            } 
            //check electrons with gsf track
            if (pfid==11 && pf->gsfTrackRef().isNonnull()) {
                if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
                continue;//and avoid double counting
            }
            //then check anything that has a ctf track
            if (pfTrack.isNonnull()) {//charged (with a ctf track)
                if(fabs( pfTrack->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
            }
        } 
    }
    return (pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq)/el.pt();
}

float smurfselections::muonIsoValuePF(const reco::PFCandidateCollection &pfCandCollection,
        const reco::Muon& mu, const reco::Vertex& vtx, 
        float coner, float minptn, float dzcut, int filterId)
{

    float pfciso = 0;
    float pfniso = 0;

    const reco::TrackRef siTrack  = mu.innerTrack();

    float mudz = siTrack.isNonnull() ? siTrack->dz(vtx.position()) : mu.standAloneMuon()->dz(vtx.position());

    reco::PFCandidateCollection::const_iterator pf;
    for (pf = pfCandCollection.begin(); pf != pfCandCollection.end(); ++pf) {

        float dR = reco::deltaR(pf->eta(), pf->phi(), mu.eta(), mu.phi());
        if (dR>coner) continue;

        int pfid = abs(pf->pdgId());
        if (filterId!=0 && filterId!=pfid) continue;

        float pfpt = pf->pt();
        if (pf->charge()==0) {
            //neutrals
            if (pfpt>minptn) pfniso+=pfpt;
        } else {
            //charged
            //avoid double counting of muon itself
            const reco::TrackRef pfTrack  = pf->trackRef();
            if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
            //first check electrons with gsf track
            if (abs(pf->pdgId())==11 && pf->gsfTrackRef().isNonnull()) {
                if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
                continue;//and avoid double counting
            }
            //then check anything that has a ctf track
            if (pfTrack.isNonnull()) {//charged (with a ctf track)
                if( fabs(pfTrack->dz(vtx.position()) - mudz) <dzcut) {//dz cut
                    pfciso+=pfpt;
                }
            } 
        }
    } 

    return (pfciso+pfniso)/mu.pt();
}

unsigned int smurfselections::CountGoodPV(const edm::Handle<reco::VertexCollection> &pvCollection)
{
    unsigned int nvtx = 0; 
    for (unsigned int ivtx = 0; ivtx < pvCollection->size(); ++ivtx)
    { 
        if (pvCollection->at(ivtx).isFake())                    continue; 
        if (pvCollection->at(ivtx).ndof() <= 4.)                continue; 
        if (pvCollection->at(ivtx).position().Rho() > 2.0)      continue; 
        if (fabs(pvCollection->at(ivtx).position().Z()) > 24.0) continue; 
        nvtx++; 
    }
    return nvtx;
}

//
// 2011 selections
//

bool smurfselections::passMuonFO2011(const edm::View<reco::Muon>::const_iterator &muon, 
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex)
{

    const reco::TrackRef siTrack = muon->innerTrack();
    float d0 = siTrack.isNonnull() ? fabs(muon->innerTrack()->dxy(vertex.position()))    : 999.9;
    float dz = siTrack.isNonnull() ? fabs(muon->innerTrack()->dz(vertex.position()))     : 999.9;
    unsigned int nValidHits = siTrack.isNonnull()      ?  muon->innerTrack()->numberOfValidHits()                   : 0;
    unsigned int nValidPixelHits = siTrack.isNonnull() ?  muon->innerTrack()->hitPattern().numberOfValidPixelHits() : 0;
    float ptErr = siTrack.isNonnull()      ?  muon->innerTrack()->ptError() : 99999.9;
    float kink  = muon->combinedQuality().trkKink;

    if (fabs(muon->eta()) > 2.4)    return false;
    if (nValidHits < 11)            return false;
    if (dz >= 0.1)                  return false;
    if (d0 >= 0.2)                  return false;
    if (nValidPixelHits == 0)       return false;
    if (ptErr / muon->pt() > 0.1)   return false;
    if (kink > 20.0)                return false;

    bool goodGlobalMuon = true;
    if (muon->isGlobalMuon()) {
        const reco::TrackRef globalTrack = muon->globalTrack();
        float nValidMuonHits = globalTrack.isNonnull()   ? globalTrack->hitPattern().numberOfValidMuonHits() : 0;
        float chi2ndof       = globalTrack.isNonnull()   ? globalTrack->chi2()/globalTrack->ndof()           : 999.9;
        if (muon->numberOfMatches() < 2)    goodGlobalMuon = false;
        if (chi2ndof >= 10)                 goodGlobalMuon = false;
        if (nValidMuonHits == 0)            goodGlobalMuon = false;
    } else goodGlobalMuon = false;

    bool goodTrackerMuon = true;
    if (muon->isTrackerMuon()) {
        if (!muon::isGoodMuon(*muon, muon::TMLastStationTight))       goodTrackerMuon = false;
    } else goodTrackerMuon = false;

    if (!(goodTrackerMuon || goodGlobalMuon))   return false;
    float isoVal = muonIsoValuePF(pfCandCollection, *muon, vertex, 0.3, 1.0, 0.1, 0);
    if (isoVal / muon->pt() > 0.40) return false;

    // pass muon fo
    return true;

}

bool smurfselections::passMuonID2011(const edm::View<reco::Muon>::const_iterator &muon, 
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex)
{

    // muon must pass FO
    if (!passMuonFO2011(muon, pfCandCollection, vertex))      return false;

    // only additional ID is tighter d0 cut
    const reco::TrackRef siTrack = muon->innerTrack();
    float d0 = siTrack.isNonnull() ? fabs(muon->innerTrack()->dxy(vertex.position())) : 999.9;
    if (muon->pt() > 20.0) {
        if (d0 >= 0.02)    return false;
    } else {
        if (d0 >= 0.01)    return false;
    }

    // pass muon id
    return true;

}

bool smurfselections::passMuonIso2011(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex)
{
    float isoVal = muonIsoValuePF(pfCandCollection, *muon, vertex, 0.3, 1.0, 0.1, 0);
    float pt = muon->pt();
    float eta = fabs(muon->eta());
    if (pt > 20.0) {
        if      (eta < 1.479  && isoVal >= 0.13)              return false;
        if      (eta >= 1.479 && isoVal >= 0.09)              return false;
    } else {
        if      (eta < 1.479  && isoVal >= 0.06)              return false;
        if      (eta >= 1.479 && isoVal >= 0.05)              return false;
    }
    return true;
}

bool smurfselections::passElectronFO2011(const reco::GsfElectronRef &electron, 
        const reco::Vertex &vertex, const Point &beamspot, 
        const edm::Handle<reco::ConversionCollection> &conversions)
{

    float pt = electron->pt();
    float d0 = fabs(electron->gsfTrack()->dxy(vertex.position()));
    float dz = fabs(electron->gsfTrack()->dz(vertex.position()));
    if (electron->dr03TkSumPt()/pt      > 0.2)      return false;
    if (electron->dr03HcalTowerSumEt()/pt > 0.2)    return false;
    if (d0 > 0.02)                                  return false;
    if (dz > 0.1)                                   return false;

    unsigned int mhits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    bool conv = ConversionTools::hasMatchedConversion(*electron, conversions, beamspot);
    if (mhits > 0)      return false;
    if (conv)           return false;

    if (electron->isEB()) {
        if (electron->sigmaIetaIeta()                               > 0.01)  return false;
        if (fabs(electron->deltaEtaSuperClusterTrackAtVtx())        > 0.007) return false;
        if (fabs(electron->deltaPhiSuperClusterTrackAtVtx())        > 0.15)  return false;
        if (electron->hadronicOverEm()                              > 0.12)  return false; 
        if (std::max(electron->dr03EcalRecHitSumEt() - 1.0, 0.0)/pt > 0.2)   return false;
    } else {
        if (electron->sigmaIetaIeta()                               > 0.03)  return false;
        if (fabs(electron->deltaEtaSuperClusterTrackAtVtx())        > 0.009) return false;
        if (fabs(electron->deltaPhiSuperClusterTrackAtVtx())        > 0.10)  return false;
        if (electron->hadronicOverEm()                              > 0.10)  return false;
        if (electron->dr03EcalRecHitSumEt()/pt                      > 0.2)   return false;
    }

    return true;
}

bool smurfselections::passElectronID2011(const reco::GsfElectronRef &electron,
        const reco::Vertex &vertex, const Point &beamspot,  
        const edm::Handle<reco::ConversionCollection> &conversions,
        const float &mvaValue)
{

    float pt = electron->pt();
    float eta = fabs(electron->superCluster()->eta());

    if (pt < 20.0) {
        if (eta >= 0.0 && eta < 1.0 &&      mvaValue <= 0.139)     return false;
        if (eta >= 1.0 && eta < 1.479 &&    mvaValue <= 0.525)     return false;
        if (eta >= 1.479 && eta < 2.5 &&    mvaValue <= 0.543)     return false;
    } else {
        if (eta >= 0.0 && eta < 1.0 &&      mvaValue <= 0.947)     return false;
        if (eta >= 1.0 && eta < 1.479 &&    mvaValue <= 0.950)     return false;
        if (eta >= 1.479 && eta < 2.5 &&    mvaValue <= 0.884)     return false;
    }

    if (electron->isEB()) {
        if (electron->sigmaIetaIeta()                             > 0.01)  return false;
    } else {
        if (electron->sigmaIetaIeta()                             > 0.03)  return false;
    }

    if (!passElectronFO2011(electron, vertex, beamspot, conversions))   return false;
    return true;
}

bool smurfselections::passElectronIso2011(const reco::GsfElectronRef &electron,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex)
{
    float isoVal = electronIsoValuePF(pfCandCollection, *electron, vertex, 0.4, 1.0, 0.1, 0.07, 0.025, -999., 0);
    if      (fabs(electron->superCluster()->eta()) < 1.479  && isoVal >= 0.13)    return false;
    else if (fabs(electron->superCluster()->eta()) >= 1.479 && isoVal >= 0.09)    return false;
    return true;
}

std::pair<double,double> smurfselections::trackerMET(std::vector<const reco::Candidate*>& objs, double deltaZCut, 
						     const reco::PFCandidateCollection &pfCandCollection,
						     const reco::Vertex &vertex) 
{
  double pX(0), pY(0);
  for (unsigned int it_o=0;it_o<objs.size();it_o++) {
    pX -= objs[it_o]->px();
    pY -= objs[it_o]->py();
  }
  reco::PFCandidateCollection::const_iterator pf;
  for (pf = pfCandCollection.begin(); pf != pfCandCollection.end(); ++pf) {
    //reject neutrals
    if (pf->charge()==0) continue;
    //avoid overlap with object under test
    bool overlap = false;
    for (unsigned int it_o=0;it_o<objs.size();it_o++) {
      float dR = reco::deltaR(pf->eta(), pf->phi(), objs[it_o]->eta(), objs[it_o]->phi());
      if ( fabs(dR)<0.1 ) overlap = true;
    }
    if (overlap) continue;
    //dz cut
    const reco::TrackRef pfTrack  = pf->trackRef();
    if (pfTrack.isNonnull()==false) continue;
    if( fabs(pfTrack->dz(vertex.position()))>deltaZCut ) continue;
    pX -= pf->px();
    pY -= pf->py();
  }

  double met    = sqrt(pX * pX + pY * pY);
  double metphi = atan2(pY, pX);
  return std::make_pair(met, metphi);

}

bool smurfselections::threeChargesAgree(const reco::GsfElectron &ele)
{
    #ifdef RELEASE_4XY
    const reco::TrackRef     ctfTkRef = ele.closestCtfTrackRef();
    #else 
    const reco::TrackRef     ctfTkRef = ele.closestTrack();
    #endif
    const reco::GsfTrackRef  gsfTkRef = ele.gsfTrack();
    if (ctfTkRef.isNonnull() && gsfTkRef.isNonnull()) {
        if ((ctfTkRef->charge() == gsfTkRef->charge()) && (ctfTkRef->charge() == ele.scPixCharge())) return true;
    }
    return false;
}

bool smurfselections::passPhotonSelection2011(const edm::View<reco::Photon>::const_iterator &photon, const double &rhoIso)
{
        float eta = fabs(photon->eta());
        float pt = photon->pt();
        if (pt < 30.0)                          return false;
        if (eta > 3.0)                          return false;
        if (photon->hasPixelSeed())             return false;
        if (photon->hadronicOverEm() >= 0.05)   return false;
        //r9 cut
        if (photon->r9() < 0.9)   return false;

        if (eta <= 1.479) {
            if (photon->sigmaIetaIeta() >= 0.011)   return false;
            if (photon->sigmaIetaIeta() <  0.001)   return false;
            if (photon->trkSumPtHollowConeDR04()  >= 2.0 + 0.001*pt + 0.0167*rhoIso)     return false;
            if (photon->ecalRecHitSumEtConeDR04() >= 4.2 + 0.006*pt + 0.183*rhoIso)      return false;
            if (photon->hcalTowerSumEtConeDR04()  >= 2.2 + 0.0025*pt + 0.062*rhoIso)     return false;
        } else {
            if (photon->sigmaIetaIeta() >= 0.03)    return false;
            if (photon->trkSumPtHollowConeDR04()  >= 2.0 + 0.001*pt + 0.032*rhoIso)      return false;
            if (photon->ecalRecHitSumEtConeDR04() >= 4.2 + 0.006*pt + 0.090*rhoIso)      return false;
            if (photon->hcalTowerSumEtConeDR04()  >= 2.2 + 0.0025*pt + 0.180*rhoIso)     return false;
        }

    return true;

}

//
// 2012 selections
//

bool smurfselections::passElectronFO2012(const reco::GsfElectronRef &electron,
        const reco::Vertex &vertex, const Point &beamspot,
        const edm::Handle<reco::ConversionCollection> &conversions)
{

    float pt = electron->pt();
    float d0 = fabs(electron->gsfTrack()->dxy(vertex.position()));
    float dz = fabs(electron->gsfTrack()->dz(vertex.position()));
    if (electron->dr03TkSumPt()/pt          > 0.2)  return false;
    if (electron->dr03HcalTowerSumEt()/pt   > 0.2)  return false;
    if (d0 > 0.02)                                  return false;
    if (dz > 0.1)                                   return false;

    unsigned int mhits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    bool conv = ConversionTools::hasMatchedConversion(*electron, conversions, beamspot);
    if (mhits > 0)      return false;
    if (conv)           return false;

    if (electron->isEB()) {
        if (electron->sigmaIetaIeta()                               > 0.01)  return false;
        if (fabs(electron->deltaEtaSuperClusterTrackAtVtx())        > 0.007) return false;
        if (fabs(electron->deltaPhiSuperClusterTrackAtVtx())        > 0.15)  return false;
        if (electron->hadronicOverEm()                              > 0.12)  return false;
        if ((electron->dr03EcalRecHitSumEt() - 1.0)/pt > 0.2)   return false;
    } else {
        if (electron->sigmaIetaIeta()                               > 0.03)  return false;
        if (fabs(electron->deltaEtaSuperClusterTrackAtVtx())        > 0.009) return false;
        if (fabs(electron->deltaPhiSuperClusterTrackAtVtx())        > 0.10)  return false;
        if (electron->hadronicOverEm()                              > 0.10)  return false;
        if (electron->dr03EcalRecHitSumEt()/pt                      > 0.2)   return false;
    }

    return true;
}

bool smurfselections::passElectronID2012(const reco::GsfElectronRef &electron,
        const reco::Vertex &vertex, const Point &beamspot,
        const edm::Handle<reco::ConversionCollection> &conversions,
        const float &mvaValue)
{

    float eta = fabs(electron->superCluster()->eta());
    if (electron->pt() > 20.0) {
        if (eta < 0.8                 && mvaValue <= 0.94) return false;
        if (eta >= 0.8 && eta < 1.479 && mvaValue <= 0.85) return false;
        if (eta >= 1.479              && mvaValue <= 0.92) return false;
    } else {
        if (eta < 0.8                 && mvaValue <= 0.00) return false;
        if (eta >= 0.8 && eta < 1.479 && mvaValue <= 0.10) return false;
        if (eta >= 1.479              && mvaValue <= 0.62) return false;
    }

    if (!passElectronFO2012(electron, vertex, beamspot, conversions))   return false;
    return true;
}

bool smurfselections::passElectronIso2012(const reco::GsfElectronRef &electron,
        const float &pfiso_ch, const float &pfiso_em, const float &pfiso_nh,
        const float &ea, const float &rho)
{
    float rhoPrime  = std::max(float(0.0), rho);
    float neutral   = std::max(float(0.0), pfiso_em + pfiso_nh - ea * rhoPrime);
    if ((pfiso_ch + neutral)/electron->pt() > 0.15) return false;
    return true;
}

bool smurfselections::passMuonFO2012(const edm::View<reco::Muon>::const_iterator &muon,
        const float &mvaValue,
        const reco::Vertex &vertex)
{
    
    const reco::TrackRef siTrack = muon->innerTrack();
    float d0 = siTrack.isNonnull() ? fabs(muon->innerTrack()->dxy(vertex.position()))    : 999.9;
    float dz = siTrack.isNonnull() ? fabs(muon->innerTrack()->dz(vertex.position()))     : 999.9;
    unsigned int nLayers = siTrack.isNonnull()      ?  muon->track()->hitPattern().trackerLayersWithMeasurement() : 0;
    unsigned int nValidPixelHits = siTrack.isNonnull() ?  muon->innerTrack()->hitPattern().numberOfValidPixelHits() : 0;
    float ptErr = siTrack.isNonnull()      ?  muon->innerTrack()->ptError() : 99999.9;
    float kink  = muon->combinedQuality().trkKink;

    if (fabs(muon->eta()) > 2.4)    return false;
    if (nLayers <= 5)               return false;
    if (dz >= 0.1)                  return false;
    if (d0 >= 0.2)                  return false;
    if (nValidPixelHits == 0)       return false;
    if (ptErr / muon->pt() > 0.1)   return false;
    if (kink > 20.0)                return false;
    if (mvaValue <= -0.6)           return false;
    #ifdef RELEASE_52X
    if (!muon->isPFMuon())          return false;
    #endif

    bool goodGlobalMuon = true;
    if (muon->isGlobalMuon()) {
        const reco::TrackRef globalTrack = muon->globalTrack();
        float nValidMuonHits = globalTrack.isNonnull()   ? globalTrack->hitPattern().numberOfValidMuonHits() : 0;
        float chi2ndof       = globalTrack.isNonnull()   ? globalTrack->chi2()/globalTrack->ndof()           : 999.9;
        if (muon->numberOfMatches() < 2)    goodGlobalMuon = false;
        if (chi2ndof >= 10)                 goodGlobalMuon = false;
        if (nValidMuonHits == 0)            goodGlobalMuon = false;
    } else goodGlobalMuon = false;
    
    bool goodTrackerMuon = true;
    if (muon->isTrackerMuon()) {
        if (!muon::isGoodMuon(*muon, muon::TMLastStationTight))       goodTrackerMuon = false;
    } else goodTrackerMuon = false;                 

    if (!(goodTrackerMuon || goodGlobalMuon))   return false;
    
    // pass muon fo     
    return true;
}

bool smurfselections::passMuonID2012(const edm::View<reco::Muon>::const_iterator &muon,
        const float &mvaValue,
        const reco::Vertex &vertex)
{

    // muon must pass FO
    if (!passMuonFO2012(muon, mvaValue, vertex))      return false;
        
    // only additional ID is tighter d0 cut
    const reco::TrackRef siTrack = muon->innerTrack();              
    float d0 = siTrack.isNonnull() ? fabs(muon->innerTrack()->dxy(vertex.position())) : 999.9;
    if (muon->pt() > 20.0) { 
        if (d0 >= 0.02)    return false;
    } else {
        if (d0 >= 0.01)    return false;
    }

    // pass muon id
    return true;

}

bool smurfselections::passMuonIso2012(const edm::View<reco::Muon>::const_iterator &muon,
        const float &mvaValue)
{
    float eta = fabs(muon->eta());
    if (muon->pt() > 20.0) {
        if (eta < 1.479  && mvaValue <= 0.82)   return false;
        if (eta >= 1.479 && mvaValue <= 0.86)   return false;
    } else {        
        if (eta < 1.479  && mvaValue <= 0.86)   return false;
        if (eta >= 1.479 && mvaValue <= 0.82)   return false;
    }               
    return true;
}

void smurfselections::PFIsolation2012(const reco::GsfElectron& el, const reco::PFCandidateCollection &pfCands, 
        const reco::VertexCollection &vertexCollection, const PFPileUpAlgo *pfPileUpAlgo,
        const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh, float &dbeta,
        bool applyEBVeto, bool applyEEVeto, bool emulatePFNoPileup, bool removeElectronTracks)
{   

    // isolation sums
    pfiso_ch = 0.0;
    pfiso_em = 0.0; 
    pfiso_nh = 0.0;
    dbeta = 0.0;

    // loop on pfcandidates
    reco::PFCandidateCollection::const_iterator pf = pfCands.begin();
    for (pf = pfCands.begin(); pf != pfCands.end(); ++pf) {
            
        // skip electrons and muons
        if (pf->particleId() == reco::PFCandidate::e)     continue;
        if (pf->particleId() == reco::PFCandidate::mu)    continue;
    
        // deltaR between electron and cadidate
        const float dR      = deltaR(pf->eta(), pf->phi(), el.eta(), el.phi());
        const float dEta    = fabs(pf->eta() - el.eta());
        if (dR > R)                             continue;

        if (pf->particleId() == reco::PFCandidate::h) {

            // vertex matching
            if (emulatePFNoPileup) {
                if (pfPileUpAlgo->chargedHadronVertex(vertexCollection, *pf) != vertexIndex) {
                    dbeta += pf->pt();
                    continue;
                }
            } else {
                const math::XYZPoint &pv = vertexCollection.at(vertexIndex).position();
                float eldz = el.gsfTrack().isNonnull() ? el.gsfTrack()->dz(pv) : el.closestCtfTrackRef()->dz(pv);
                if (fabs(pf->trackRef()->dz(pv) - eldz) > 0.1) {
                    dbeta += pf->pt();
                    continue;
                }
            }

            // remove electron tracks
            if (removeElectronTracks) {
                const reco::TrackRef    pfTrack     = pf->trackRef();
                const reco::GsfTrackRef pfGsfTrack  = pf->gsfTrackRef();
                const reco::TrackRef    elTrack     = el.closestCtfTrackRef();
                const reco::GsfTrackRef elGsfTrack  = el.gsfTrack();
                if (pfGsfTrack.isNonnull() && elGsfTrack.isNonnull()
                    && pfGsfTrack.key() ==  elGsfTrack.key()) continue;
                if (pfTrack.isNonnull() && elTrack.isNonnull() 
                    && pfTrack.key() ==  elTrack.key()) continue;
            }

        }

        // geometric vetoes
        if (!el.isEB() && applyEEVeto) {
            if (pf->particleId() == reco::PFCandidate::h      && dR <= 0.015)   continue;
            if (pf->particleId() == reco::PFCandidate::gamma  && dR <= 0.08)    continue;
        }
        if (el.isEB() && applyEBVeto) {
            if (pf->particleId() == reco::PFCandidate::h      && dR <= 0.015)   continue;
            if (pf->particleId() == reco::PFCandidate::gamma  && dEta <= 0.025)    continue;
        }

        // remove mis-identified electron super-clusters
        if (el.gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 0 
                && pf->mva_nothing_gamma() > 0.99 && el.superCluster().isNonnull() 
                && pf->superClusterRef().isNonnull() 
                && el.superCluster() == pf->superClusterRef())   continue;

        // add to isolation sum
        if (pf->particleId() == reco::PFCandidate::h)       pfiso_ch += pf->pt();
        if (pf->particleId() == reco::PFCandidate::gamma)   pfiso_em += pf->pt();
        if (pf->particleId() == reco::PFCandidate::h0)      pfiso_nh += pf->pt();

    }

}

double smurfselections::getElectronRadialIsolation(const reco::GsfElectron &el, const reco::PFCandidateCollection &PFCandidates, 
        double cone_size, double neutral_et_threshold, bool barrel_veto)
{   
    double radial_iso = 0;

    for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); iP != PFCandidates.end(); ++iP) {

        //************************************************************
        // Veto any PFmuon, or PFEle
        if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
        //************************************************************

        //************************************************************
        // New Isolation Calculations
        //************************************************************
        double dr = reco::deltaR(iP->p4().eta(), iP->p4().phi(), el.p4().eta(), el.p4().phi());

        if (dr > cone_size) continue;

        if (!el.isEB()) {
            if (iP->particleId() == reco::PFCandidate::h && dr <= 0.015) continue;
            if (iP->particleId() == reco::PFCandidate::gamma && dr <=0.08) continue;
        }
        else if (barrel_veto && el.mvaOutput().mva < -0.1) {
            if (iP->particleId() == reco::PFCandidate::h && dr <= 0.015) continue;
            if (iP->particleId() == reco::PFCandidate::gamma && dr <=0.08) continue;            
        }
    
        //Charged
        if(iP->trackRef().isNonnull()) {
            radial_iso += iP->pt() * (1 - 3*dr) / el.pt();
        }
        else if (iP->pt() > neutral_et_threshold) {
            radial_iso += iP->pt() * (1 - 3*dr) / el.pt();
        } 
    } //loop over PF candidates
        
    return radial_iso;        
}

//
// muons
//

bool smurfselections::passMuonHPASameSign(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex)
{               
    if (!passMuonIsPOGTight(muon, vertex))                                  return false;
    if (muon->isEnergyValid() && muon->isolationR03().emVetoEt > 4.0)       return false;
    if (muon->isEnergyValid() && muon->isolationR03().hadVetoEt > 6.0)      return false;
    if (fabs(muon->innerTrack()->dxy(vertex.position())) >= 0.02)           return false;
    if (fabs(muon->innerTrack()->dz(vertex.position())) >= 0.1)             return false;
    return true;
}    

bool smurfselections::passMuonIsPOGTight(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex)
{
    if (!muon->isGlobalMuon())                                              return false;
    #ifdef RELEASE_52X
    if (!muon->isPFMuon())                                                  return false;
    #endif
    if (!muon->globalTrack().isNonnull())                                   return false;
    if (muon->globalTrack()->normalizedChi2() >= 10.)                       return false;
    if (muon->globalTrack()->hitPattern().numberOfValidMuonHits() == 0)     return false;
    if (muon->numberOfMatchedStations() <= 1)                               return false;
    if (!muon->innerTrack().isNonnull())                                    return false;
    if (fabs(muon->innerTrack()->dxy(vertex.position())) >= 0.2)            return false;
    if (muon->innerTrack()->hitPattern().numberOfValidPixelHits() == 0)     return false;
    if (!muon->track().isNonnull())                                         return false;
    if (muon->track()->hitPattern().trackerLayersWithMeasurement() <= 5)    return false;
    return true;
}

bool smurfselections::passMuonIsPOGSoft(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex)
{

    if (!muon::isGoodMuon(*muon, muon::TMOneStationTight))                  return false;
    if (!muon->innerTrack().isNonnull())                                    return false;
    if (muon->innerTrack()->hitPattern().numberOfValidTrackerHits() <= 10)  return false;
    if (muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1) return false;
    if (muon->innerTrack()->normalizedChi2() >= 1.8)                        return false;
    if (fabs(muon->innerTrack()->dxy(vertex.position())) >= 3.0)            return false;
    if (fabs(muon->innerTrack()->dz(vertex.position())) >= 30.0)            return false;
    return true;
}

// original implementation of radial PF muon isolation taken from here
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/folgueras/MuonIsolation/Tools/src/RadialIsolation.cc?view=markup
double smurfselections::getMuonRadialIsolation(const reco::Muon &mu, const reco::PFCandidateCollection &PFCandidates, 
        double cone_size, double neutral_et_threshold)
{
    double radial_iso = 0;

    reco::TrackRef muTrk = mu.track();
    if (muTrk.isNull())
        muTrk = mu.standAloneMuon();
    if (muTrk.isNull())
        return -9999;

    for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); iP != PFCandidates.end(); ++iP) {
        if(iP->trackRef().isNonnull() && mu.innerTrack().isNonnull() && edm::refToPtr(iP->trackRef()) == edm::refToPtr(mu.innerTrack())) continue;  // exclude the muon itself

        //************************************************************
        // New Isolation Calculations
        //************************************************************
        double dr = reco::deltaR(iP->p4().eta(), iP->p4().phi(), mu.p4().eta(), mu.p4().phi());

        if (dr > cone_size) continue;
        if (dr < 0.01) continue;  // inner veto cone

        //Charged
        if(iP->trackRef().isNonnull()) {
            //************************************************************
            // Veto any PFmuon, or PFEle
            if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
            //************************************************************
            radial_iso += iP->pt() * (1 - 3*dr) / mu.pt();
        }
        else if (iP->pt() > neutral_et_threshold) {
            radial_iso += iP->pt() * (1 - 3*dr) / mu.pt();
        }
    } //loop over PF candidates

    return radial_iso;
}



