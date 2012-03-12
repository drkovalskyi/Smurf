
#include "Smurf/ProcessingAndSkimming/interface/Selections.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include <algorithm>

//
// utilities
//

float smurfselections::d0Vertex(const float &d0, const float &phi, const Point &vertex)
{
    return (d0 - vertex.x() * sin(phi) + vertex.y() * cos(phi));
}

float smurfselections::dzVertex(const reco::Candidate &cand, const Point &vertex)
{
    float dz = (cand.vz() - vertex.z()) - 
        ((cand.vx() - vertex.x()) * cand.px() + (cand.vy() - vertex.y()) * cand.py()) / cand.pt() * cand.pz()/cand.pt();
    return dz;
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
            if (pfid==11 && pf->gsfTrackRef().isNonnull()) {
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
    return pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq;
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
                if(fabs( pfTrack->dz(vtx.position()) - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
            } 
        }
    } 
    return pfciso+pfniso;
}

//
// 2011 selections
//

bool smurfselections::passMuonFO2011(const edm::View<reco::Muon>::const_iterator &muon, const Point &vertex)
{
    float d0 = fabs(d0Vertex(muon->innerTrack()->d0(), muon->innerTrack()->phi(), vertex));
    float dz = fabs(dzVertex(*muon, vertex));

    if (fabs(muon->eta()) > 2.4)                                        return false;
    if (muon->innerTrack()->numberOfValidHits() < 11)                   return false;
    if (dz >= 0.1)                                                      return false;
    if (d0 >= 0.02)                                                     return false;
    if (muon->innerTrack()->hitPattern().numberOfValidPixelHits() == 0) return false;

    bool goodGlobalMuon = true;
    if (muon->isGlobalMuon()) {
        if (muon->numberOfMatches() < 2)                                    goodGlobalMuon = false;
        if (muon->globalTrack()->chi2()/muon->globalTrack()->ndof() >= 10)  goodGlobalMuon = false;
        if (muon->globalTrack()->hitPattern().numberOfValidMuonHits() == 0) goodGlobalMuon = false;
    } else goodGlobalMuon = false;

    bool goodTrackerMuon = true;
    if (muon->isTrackerMuon()) {
        if (!muon::isGoodMuon(*muon, muon::TMLastStationTight))       goodTrackerMuon = false;
    } else goodTrackerMuon = false;

    if (!(goodTrackerMuon || goodGlobalMuon))   return false;

    // pass muon fo
    return true;

}

bool smurfselections::passMuonID2011(const edm::View<reco::Muon>::const_iterator &muon, const Point &vertex)
{

    // muon must pass FO
    if (!passMuonFO2011(muon, vertex))      return false;

    // only additional ID is tighter d0 cut
    float d0 = fabs(d0Vertex(muon->innerTrack()->d0(), muon->innerTrack()->phi(), vertex));
    if (muon->pt() > 20.0 && d0 >= 0.01)    return false;
    else if (d0 >= 0.02)                    return false;

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
        if (eta < 1.479 && isoVal >= 0.13)    return false;
        else if (isoVal >= 0.09)              return false;
    } else {
        if (eta < 1.479 && isoVal >= 0.06)    return false;
        else if (isoVal >= 0.05)              return false;
    }
    return true;
}

bool smurfselections::passElectronFO2011(const edm::View<reco::GsfElectron>::const_iterator &electron, 
        const Point &vertex, const Point &beamspot, 
        const edm::Handle<reco::ConversionCollection> &conversions)
{

    float pt = electron->pt();
    float d0 = fabs(d0Vertex(electron->gsfTrack()->d0(), electron->gsfTrack()->phi(), vertex));
    float dz = fabs(dzVertex(*electron, vertex));
    if (electron->dr03TkSumPt()/pt      > 0.2)      return false;
    if (electron->dr03HcalTowerSumEt()  > 0.2)      return false;
    if (d0 > 0.02)                                  return false;
    if (dz > 0.1)                                   return false;

    unsigned int mhits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    bool conv = ConversionTools::hasMatchedConversion(*electron, conversions, beamspot);
    if (mhits > 0)      return false;
    if (conv)           return false;

    if (electron->isEB()) {
        if (electron->scSigmaIEtaIEta()                             > 0.01)  return false;
        if (electron->deltaEtaSuperClusterTrackAtVtx()              > 0.007) return false;
        if (electron->deltaPhiSuperClusterTrackAtVtx()              > 0.15)  return false;
        if (electron->hadronicOverEm()                              > 0.12)  return false; 
        if (std::max(electron->dr03EcalRecHitSumEt() - 1.0, 0.0)/pt > 0.2)   return false;
    } else {
        if (electron->scSigmaIEtaIEta()                             > 0.03)  return false;
        if (electron->deltaEtaSuperClusterTrackAtVtx()              > 0.009) return false;
        if (electron->deltaPhiSuperClusterTrackAtVtx()              > 0.10)  return false;
        if (electron->hadronicOverEm()                              > 0.10)  return false;
        if (electron->dr03EcalRecHitSumEt()/pt                      > 0.2)   return false;
    }

    return true;
}

bool smurfselections::passElectronID2011(const edm::View<reco::GsfElectron>::const_iterator &electron,
        const Point &vertex, const Point &beamspot,  
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

    if (!passElectronFO2011(electron, vertex, beamspot, conversions))   return false;
    return true;
}

bool smurfselections::passElectronIso2011(const edm::View<reco::GsfElectron>::const_iterator &electron,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex)
{
    float isoVal = electronIsoValuePF(pfCandCollection, *electron, vertex, 0.4, 1.0, 0.1, 0.07, 0.025, -999., 0);
    if (fabs(electron->superCluster()->eta()) < 1.479 && isoVal >= 0.13)    return false;
    else if (isoVal >= 0.09)                                                return false;
    return true;
}

