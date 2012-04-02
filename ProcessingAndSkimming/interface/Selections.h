#ifndef SMURF_PROCESSINGANDSKIMMING_SELECTIONS_H
#define SMURF_PROCESSINGANDSKIMMING_SELECTIONS_H

//
// Get around backward incompatibility in
// interface for JetCorrector::correction
//#define RELEASE_52X
#define RELEASE_4XY
//
//

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

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
 bool compareJetPt(std::pair<reco::PFJet, float> lv1, std::pair<reco::PFJet, float> lv2);
std::vector<std::pair<reco::PFJet, float> > goodJets(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const edm::Handle<edm::View<reco::PFJet> > &jets_h, const reco::Candidate &cand1, const reco::Candidate &cand2,
        const JetCorrector *corrector, float ptCut);
std::vector<std::pair<reco::PFJet, float> > goodJets(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const edm::Handle<edm::View<reco::PFJet> > &jets_h, const reco::Candidate &cand,
	const JetCorrector *corrector, float ptCut);
//trackerMET: the objs vector should contain the two leptons or the photon or whatever you do not want to be taken from the PFCandidateCollection
std::pair<double,double> trackerMET(std::vector<const reco::Candidate*>& objs, 
				    double deltaZCut, 
				    const reco::PFCandidateCollection &pfCandCollection,
				    const reco::Vertex &vertex);

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

