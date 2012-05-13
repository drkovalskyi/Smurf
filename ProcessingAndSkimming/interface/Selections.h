#ifndef SMURF_PROCESSINGANDSKIMMING_SELECTIONS_H
#define SMURF_PROCESSINGANDSKIMMING_SELECTIONS_H

//
// Get around backward incompatibility in
// interface for JetCorrector::correction
#define RELEASE_52X
//#define RELEASE_4XY
//
//

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"

namespace smurfselections {

typedef math::XYZPoint Point;

//
// utilities
//

float d0Vertex(const float &d0, const float &phi, const reco::Vertex &vertex);
float dzVertex(const reco::Candidate &cand, const reco::Vertex &vertex);
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
unsigned int CountGoodPV(const edm::Handle<reco::VertexCollection> &pvCollection);

//
// 2011 selections
//

bool passMuonFO2011(const edm::View<reco::Muon>::const_iterator &muon, 
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex);
bool passMuonID2011(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex);
bool passMuonIso2011(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex);

bool passElectronFO2011(const reco::GsfElectronRef &electron, 
        const reco::Vertex &vertex, const Point &beamspot,  
        const edm::Handle<reco::ConversionCollection> &conversions);
bool passElectronID2011(const reco::GsfElectronRef &electron,
        const reco::Vertex &vertex, const Point &beamspot,  
        const edm::Handle<reco::ConversionCollection> &conversions, 
        const float &mvaValue);
bool passElectronIso2011(const reco::GsfElectronRef &electron,
        const reco::PFCandidateCollection &pfCandCollection,
        const reco::Vertex &vertex);

bool threeChargesAgree(const reco::GsfElectron &ele);

bool passPhotonSelection2011(const edm::View<reco::Photon>::const_iterator &photon,
        const double &rhoIso);

//
// 2012 selections
//

bool passElectronFO2012(const reco::GsfElectronRef &electron,
        const reco::Vertex &vertex, const Point &beamspot,
        const edm::Handle<reco::ConversionCollection> &conversions);
bool passElectronID2012(const reco::GsfElectronRef &electron,
        const reco::Vertex &vertex, const Point &beamspot,
        const edm::Handle<reco::ConversionCollection> &conversions,
        const float &mvaValue);
bool passElectronIso2012(const reco::GsfElectronRef &electron,
        const float &pfiso_ch, const float &pfiso_em, const float &pfiso_nh,
        const float &ea, const float &rho);

bool passMuonFO2012(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex);
bool passMuonID2012(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex);
bool passMuonIso2012(const edm::View<reco::Muon>::const_iterator &muon,
        const float &mvaValue);


// for 2012 pf isolation
void PFIsolation2012(const reco::GsfElectron& el, const reco::PFCandidateCollection &pfCands,
        const reco::VertexCollection &vertexCollection, const PFPileUpAlgo *pfPileUpAlgo,
        const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh, float &dbeta,
        bool applyEBVeto, bool applyEEVeto, bool emulatePFNoPileup, bool removeElectronTracks);

double getElectronRadialIsolation(const reco::GsfElectron &ele, const reco::PFCandidateCollection &PFCandidates, 
    double cone_size = 0.3, double neutral_et_threshold = 1., bool barrel_veto = false);

//
// for muons
//

bool passMuonHPASameSign(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex);

bool passMuonIsPOGTight(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex);
bool passMuonIsPOGSoft(const edm::View<reco::Muon>::const_iterator &muon,
        const reco::Vertex &vertex);
double getMuonRadialIsolation(const reco::Muon &mu, const reco::PFCandidateCollection &PFCandidates, 
    double cone_size = 0.3, double neutral_et_threshold = 1.);

}

#endif

