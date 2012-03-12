#ifndef SMURF_PROCESSINGANDSKIMMING_SELECTIONS_H
#define SMURF_PROCESSINGANDSKIMMING_SELECTIONS_H

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

namespace smurfselections {

typedef math::XYZPoint Point;

//
// utilities
//

float d0Vertex(const float &d0, const float &phi, const Point &vertex);
float dzVertex(const reco::Candidate &cand, const Point &vertex);
float electronIsoValuePF(const reco::PFCandidateCollection &pfCandCollection,
        const reco::GsfElectron& el, const reco::Vertex& vtx, 
        float coner, float minptn, float dzcut,
        float footprintdr, float gammastripveto, float elestripveto, int filterId);
float muonIsoValuePF(const reco::PFCandidateCollection &pfCandCollection,
        const reco::Muon& mu, const reco::Vertex& vtx, 
        float coner, float minptn, float dzcut, int filterId);

//
// 2011 selections
//

bool passMuonFO2011(const edm::View<reco::Muon>::const_iterator &muon, 
        const Point &vertex);
bool passMuonID2011(const edm::View<reco::Muon>::const_iterator &muon,
        const Point &vertex);
bool passMuonIso2011(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex);

bool passElectronFO2011(const edm::View<reco::GsfElectron>::const_iterator &electron, 
        const Point &vertex, const Point &beamspot,  
        const edm::Handle<reco::ConversionCollection> &conversions);
bool passElectronID2011(const edm::View<reco::GsfElectron>::const_iterator &electron,
        const Point &vertex, const Point &beamspot,  
        const edm::Handle<reco::ConversionCollection> &conversions, 
        const float &mvaValue);
bool passElectronIso2011(const edm::View<reco::GsfElectron>::const_iterator &electron,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex);

//
// 2012 selections
//

}

#endif

