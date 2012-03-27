// -*- C++ -*-
//
// Package:    LeptonTreeMaker
// Class:      LeptonTreeMaker
// 
/**\class LeptonTreeMaker LeptonTreeMaker.cc Smurf/LeptonTreeMaker/src/LeptonTreeMaker.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Dave Evans,510 1-015,+41227679496,
//         Created:  Thu Mar  8 11:43:50 CET 2012
// $Id: LeptonTreeMaker.cc,v 1.7 2012/03/22 17:19:56 cerati Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <vector>

// fw include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// objects
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// user include files
#include "Smurf/Core/LeptonTree.h"
#include "Smurf/ProcessingAndSkimming/interface/Selections.h"

// root include files
#include "TFile.h"
#include "TTree.h"
#include "TRegexp.h"
#include "TString.h"

//
// class declaration
//

class LeptonTreeMaker : public edm::EDProducer {
    public:
        explicit LeptonTreeMaker(const edm::ParameterSet&);
        ~LeptonTreeMaker();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

        // ----------functions to fill tree  ---------------

        // efficiencies
        void fillElectronTagAndProbeTree(const edm::Event& iEvent, 
                const TransientTrackBuilder *ttBuilder, 
                EcalClusterLazyTools *clusterTools);
        void fillMuonTagAndProbeTree(const edm::Event& iEvent);

        // fake rates
        void fillElectronFakeRateTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
                const TransientTrackBuilder *ttBuilder, 
                const EcalClusterLazyTools *clusterTools);

        void fillMuonFakeRateTree(const edm::Event& iEvent, const edm::EventSetup &iSetup);

        // photons
        void fillPhotonTree(const edm::Event& iEvent, const edm::EventSetup &iSetup);

        // common variables
        void fillCommonVariables(const edm::Event& iEvent);

        // did any trigger pass
        bool eventPassTrigger(const std::vector<std::string> &trigNames);
        bool eventPassTrigger(const std::string &trigName);
        double getTriggerPrescale(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string &trigName);
        void triggerObjects();

        // ----------member data ---------------------------

        // lepton tree
        TFile       *leptonFile_;
        LeptonTree  *leptonTree_;

        // input tags
        edm::InputTag electronsInputTag_;
        edm::InputTag muonsInputTag_;
        edm::InputTag photonsInputTag_;
        edm::InputTag primaryVertexInputTag_;
        edm::InputTag rhoJECInputTag;
        edm::InputTag rhoIsoInputTag;
        edm::InputTag pfCandsInputTag_;
        edm::InputTag metInputTag_;
        edm::InputTag jetsInputTag_;

        edm::InputTag conversionsInputTag_;
        edm::InputTag beamSpotInputTag_;

        // jet corrections
        std::string pfJetCorrectorL1FastL2L3_;

        std::string pathToBDTWeights_;

        // trigger names
        std::vector<std::string> electronTPTriggerNames_;
        std::vector<std::string> muonTPTriggerNames_;
        std::vector<std::string> electronFRTriggerNames_;
        std::vector<std::string> muonFRTriggerNames_;
        std::vector<std::string> photonTriggerNames_;

        // trigger related
        HLTConfigProvider          hltConfig_;
        std::string                processName_;
        const edm::TriggerResults* triggerResults_;
  
        // electron id related
        ElectronIDMVA *electronIDMVA_;

        // common products
        double rhoJEC_;
        double rhoIso_;
        reco::PFCandidateCollection pfCandCollection_;
        edm::View<reco::Vertex> vertexCollection_;
        reco::Vertex pv_;

};

//
// constants, enums and typedefs
//

typedef math::XYZTLorentzVectorD LorentzVector;
typedef math::XYZPoint Point;

//
// static data member definitions
//

//
// constructors and destructor
//
LeptonTreeMaker::LeptonTreeMaker(const edm::ParameterSet& iConfig)
{

    //
    // get input tags
    //

    electronsInputTag_      =  iConfig.getParameter<edm::InputTag>("electronsInputTag");
    muonsInputTag_          =  iConfig.getParameter<edm::InputTag>("muonsInputTag");
    photonsInputTag_        =  iConfig.getParameter<edm::InputTag>("photonsInputTag");
    primaryVertexInputTag_  =  iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    rhoJECInputTag          =  iConfig.getParameter<edm::InputTag>("rhoJECInputTag");
    rhoIsoInputTag          =  iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
    pfCandsInputTag_        =  iConfig.getParameter<edm::InputTag>("pfCandsInputTag");
    metInputTag_            =  iConfig.getParameter<edm::InputTag>("metInputTag");
    jetsInputTag_           =  iConfig.getParameter<edm::InputTag>("jetsInputTag");

    conversionsInputTag_    =  iConfig.getParameter<edm::InputTag>("conversionsInputTag");
    beamSpotInputTag_       =  iConfig.getParameter<edm::InputTag>("beamSpotInputTag");

    pfJetCorrectorL1FastL2L3_ = iConfig.getParameter<std::string>("pfJetCorrectorL1FastL2L3");

    pathToBDTWeights_       = iConfig.getParameter<std::string>("pathToBDTWeights");

    electronTPTriggerNames_ =  iConfig.getUntrackedParameter<std::vector<std::string> >("electronTPTriggerNames");
    muonTPTriggerNames_     =  iConfig.getUntrackedParameter<std::vector<std::string> >("muonTPTriggerNames"); 
    electronFRTriggerNames_ =  iConfig.getUntrackedParameter<std::vector<std::string> >("electronFRTriggerNames");
    muonFRTriggerNames_     =  iConfig.getUntrackedParameter<std::vector<std::string> >("muonFRTriggerNames");
    photonTriggerNames_     =  iConfig.getUntrackedParameter<std::vector<std::string> >("photonTriggerNames");

    //
    // set up lepton tree
    //

    leptonFile_ = TFile::Open("leptonTree.root", "RECREATE");
    assert(leptonFile_);
    leptonTree_ = new LeptonTree();
    leptonTree_->CreateTree();
    leptonTree_->tree_->SetDirectory(leptonFile_);

    //
    // set up trigger
    //

    processName_ = "";

    //
    // set up electron id
    //

    electronIDMVA_ = new ElectronIDMVA();
    electronIDMVA_->Initialize("BDTG method",
            pathToBDTWeights_+"/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/Subdet2HighPt_WithIPInfo_BDTG.weights.xml" ,                
            ElectronIDMVA::kWithIPInfo);

}


LeptonTreeMaker::~LeptonTreeMaker()
{

    //
    // tidy up
    //

    if (electronIDMVA_) delete electronIDMVA_;

    //
    // save and close lepton tree
    //

    leptonFile_->cd();
    leptonTree_->tree_->Write();
    leptonFile_->Close();
}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
LeptonTreeMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //
    // get trigger information
    //

    // trigger event handle
    if (processName_ == "") {
        edm::Handle<trigger::TriggerEvent> triggerEvent_h;
        iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEvent_h);
        if (!triggerEvent_h.isValid()) {
            throw cms::Exception("[LeptonTreeMaker][produce] Error getting TriggerEvent product from Event!");
        }
        processName_ = triggerEvent_h.provenance()->processName();
        bool changed(true);
        hltConfig_.init(iEvent.getRun(), iSetup, processName_, changed);
    }

    // trigger results
    edm::Handle<edm::TriggerResults> triggerResults_h;
    iEvent.getByLabel(edm::InputTag("TriggerResults", "", processName_), triggerResults_h);
    if (!triggerResults_h.isValid()) {
        throw cms::Exception("LeptonTreeMaker][produce] Error getting TriggerResults product from Event!");
    }
    triggerResults_=triggerResults_h.product();

    // sanity check
    assert(triggerResults_->size() == hltConfig_.size());

    //  
    // set up tools
    //

    // unbiased revertexing
    edm::ESHandle<TransientTrackBuilder> ttBuilder_h;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder_h);
    const TransientTrackBuilder *ttBuilder = ttBuilder_h.product();

    // ecal tools
    EcalClusterLazyTools *clusterTools = new EcalClusterLazyTools(iEvent, iSetup, 
        edm::InputTag("reducedEcalRecHitsEB"), 
        edm::InputTag("reducedEcalRecHitsEE"));

    //
    // common products  
    //

    // rho for jec
    edm::Handle<double> rhoJEC_h;
    iEvent.getByLabel(rhoJECInputTag, rhoJEC_h);
    rhoJEC_ = *(rhoJEC_h.product());

    // rho for isolation
    edm::Handle<double> rhoIso_h;
    iEvent.getByLabel(rhoIsoInputTag, rhoIso_h);
    rhoIso_ = *(rhoIso_h.product());

    // pf candidates
    edm::Handle<reco::PFCandidateCollection> pfCand_h;
    iEvent.getByLabel(pfCandsInputTag_, pfCand_h);
    pfCandCollection_ = *(pfCand_h.product());

    // vertices
    edm::Handle<edm::View<reco::Vertex> > vtx_h;
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h);
    vertexCollection_ = *(vtx_h.product());
    if (vertexCollection_.size() > 0) pv_ = vertexCollection_.at(0);
    else return;

    //triggerObjects();

    //
    // decide what to do based on trigger information
    //

    bool passElectronTPTrigger  = eventPassTrigger(electronTPTriggerNames_);
    bool passMuonTPTrigger      = eventPassTrigger(muonTPTriggerNames_);
    bool passElectronFRTrigger  = eventPassTrigger(electronFRTriggerNames_);
    bool passMuonFRTrigger      = eventPassTrigger(muonFRTriggerNames_);
    bool passPhotonTrigger      = eventPassTrigger(photonTriggerNames_);

    // fake rates
    if (passElectronFRTrigger)  fillElectronFakeRateTree(iEvent, iSetup, ttBuilder, clusterTools);
    if (passMuonFRTrigger)      fillMuonFakeRateTree(iEvent, iSetup);

    // efficiency
    if (passElectronTPTrigger)  fillElectronTagAndProbeTree(iEvent, ttBuilder, clusterTools);
    if (passMuonTPTrigger)      fillMuonTagAndProbeTree(iEvent);

    // photons
    if (passPhotonTrigger)      fillPhotonTree(iEvent, iSetup);

    //
    // tidy up
    //

    if (clusterTools) delete clusterTools;

}

// ------------ method called once each job just before starting event loop  ------------
    void 
LeptonTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeptonTreeMaker::endJob() {
}

// ------------ method called when starting to processes a run  ------------
    void 
LeptonTreeMaker::beginRun(edm::Run &r, edm::EventSetup const &c)
{
    if (processName_ != "") {
        bool changed(true);
        if (!hltConfig_.init(r, c, processName_, changed)) {
            throw cms::Exception("[LeptonTreeMaker][beginRun] Config extraction failure with process name " + processName_);
        }
    }
}

// ------------ method called when ending the processing of a run  ------------
    void 
LeptonTreeMaker::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
    void 
LeptonTreeMaker::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
    void 
LeptonTreeMaker::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

void LeptonTreeMaker::fillElectronTagAndProbeTree(const edm::Event& iEvent,
                const TransientTrackBuilder *ttBuilder,
                EcalClusterLazyTools *clusterTools)
{

    // electrons
    edm::Handle<edm::View<reco::GsfElectron> > els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    edm::View<reco::GsfElectron> electronCollection = *(els_h.product());
    if (electronCollection.size() < 2) return;

    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    iEvent.getByLabel(conversionsInputTag_, conversions_h);

    // beamspot
    edm::Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
    const reco::BeamSpot &thebs = *(beamspot_h.product());

    // look for tag and probe
    edm::View<reco::GsfElectron>::const_iterator tag;
    edm::View<reco::GsfElectron>::const_iterator probe;

    for (tag = electronCollection.begin(); tag != electronCollection.end(); ++tag) {

        // look for good tag
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.5) continue;
        float mvaValue = electronIDMVA_->MVAValue(&*tag, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
        if (!smurfselections::passElectronID2011(tag, pv_.position(), thebs.position(), conversions_h, mvaValue)) continue;
        if (!smurfselections::passElectronIso2011(tag, pfCandCollection_, pv_)) continue;

        for (probe = electronCollection.begin(); probe != electronCollection.end(); ++probe) {

            // look for good probe
            if (tag == probe)             continue;
            if (probe->pt() < 10.0)       continue;
            if (fabs(probe->eta()) > 2.5) continue;

            // construct tag and probe mass
            LorentzVector p4 = tag->p4() + probe->p4();

            // fill the tree
            fillCommonVariables(iEvent);
            leptonTree_->eventSelection_     = LeptonTree::ZeeTagAndProbe;
            leptonTree_->probe_              = probe->p4();
            leptonTree_->qProbe_             = probe->charge();
            leptonTree_->tag_                = tag->p4();
            leptonTree_->qTag_               = tag->charge();
            leptonTree_->tagAndProbeMass_    = p4.M();

            mvaValue = electronIDMVA_->MVAValue(&*probe, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
            if (smurfselections::passElectronFO2011(probe, pv_.position(), thebs.position(), conversions_h))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
            if (smurfselections::passElectronID2011(probe, pv_.position(), thebs.position(), conversions_h, mvaValue))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
            if (smurfselections::passElectronIso2011(probe, pfCandCollection_, pv_)) 
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);

            leptonTree_->tree_->Fill(); 

        }

    }

}

void LeptonTreeMaker::fillMuonTagAndProbeTree(const edm::Event& iEvent)
{

    // muons
    edm::Handle<edm::View<reco::Muon> > mus_h;
    iEvent.getByLabel(muonsInputTag_, mus_h);
    edm::View<reco::Muon> muonCollection = *(mus_h.product());
    if (muonCollection.size() < 2) return;

    // look for tag and probe
    edm::View<reco::Muon>::const_iterator tag;
    edm::View<reco::Muon>::const_iterator probe;

    for (tag = muonCollection.begin(); tag != muonCollection.end(); ++tag) {

        // look for good tag
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.4) continue;
        if (!smurfselections::passMuonID2011(tag, pv_.position()))         continue;
        if (!smurfselections::passMuonIso2011(tag, pfCandCollection_, pv_)) continue;

        for (probe = muonCollection.begin(); probe != muonCollection.end(); ++probe) {

            // look for good probe
            if (tag == probe)            continue;
            if (probe->pt() < 10.0)      continue;
            if (fabs(probe->eta()) > 2.4) continue;

            // construct tag and probe mass
            LorentzVector p4 = tag->p4() + probe->p4();

            // fill the tree
            fillCommonVariables(iEvent);
            leptonTree_->eventSelection_     = LeptonTree::ZmmTagAndProbe;
            leptonTree_->probe_              = probe->p4();
            leptonTree_->qProbe_             = probe->charge();
            leptonTree_->tag_                = tag->p4();
            leptonTree_->qTag_               = tag->charge();
            leptonTree_->tagAndProbeMass_    = p4.M();

            if (smurfselections::passMuonFO2011(probe, pv_.position()))         leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
            if (smurfselections::passMuonID2011(probe, pv_.position()))         leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
            if (smurfselections::passMuonIso2011(probe, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);

            leptonTree_->tree_->Fill();

        }

    }

}

void LeptonTreeMaker::fillElectronFakeRateTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                const TransientTrackBuilder *ttBuilder,
                const EcalClusterLazyTools *clusterTools)
{

    // electrons
    edm::Handle<edm::View<reco::GsfElectron> > els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    edm::View<reco::GsfElectron> electronCollection = *(els_h.product());
    if (electronCollection.size() < 1) return;
    
    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    iEvent.getByLabel(conversionsInputTag_, conversions_h);

    // beamspot
    edm::Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
    const reco::BeamSpot &thebs = *(beamspot_h.product());
    
    // met
    edm::Handle<edm::View<reco::PFMET> > met_h;
    iEvent.getByLabel(metInputTag_, met_h);
    float met = met_h->front().et();
    if (met > 20.0) return;

    // jets
    edm::Handle<edm::View<reco::PFJet> > jets_h;
    iEvent.getByLabel(jetsInputTag_, jets_h);
    const JetCorrector* corrector = JetCorrector::getJetCorrector(pfJetCorrectorL1FastL2L3_, iSetup);

    // look for a good FO
    edm::View<reco::GsfElectron>::const_iterator it;
    edm::View<reco::GsfElectron>::const_iterator fo;
    unsigned int nFO = 0;
    for (it = electronCollection.begin(); it != electronCollection.end(); ++it) {
        if (it->pt()        < 10.0) continue;
        if (fabs(it->eta()) > 2.4)  continue;
        if (!smurfselections::passElectronFO2011(it, pv_.position(), thebs.position(), conversions_h)) continue;
        ++nFO;
        fo = it;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {
        fillCommonVariables(iEvent);

        if (eventPassTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v*"))       leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle8;
        if (eventPassTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v*"))      leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle17VL;
        if (eventPassTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*")) leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle8VLJet40;
        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();

        float mvaValue = electronIDMVA_->MVAValue(&*fo, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
        leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
        if (smurfselections::passElectronID2011(fo, pv_.position(), thebs.position(), conversions_h, mvaValue))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
        if (smurfselections::passElectronIso2011(fo, pfCandCollection_, pv_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);

        // jets
        leptonTree_->jet1_ = LorentzVector(0,0,0,0);
        leptonTree_->jet2_ = LorentzVector(0,0,0,0);
        leptonTree_->jet3_ = LorentzVector(0,0,0,0);
        std::vector<std::pair<reco::PFJet, float> > jets30 = smurfselections::goodJets(iEvent, iSetup, jets_h, *fo, corrector, 30.);
        leptonTree_->njets_ = jets30.size();
        std::vector<std::pair<reco::PFJet, float> > jets15 = smurfselections::goodJets(iEvent, iSetup, jets_h, *fo, corrector, 15.);
        if (jets15.size() > 0) leptonTree_->jet1_ = jets15.at(0).first.p4() * jets15.at(0).second;
        if (jets15.size() > 1) leptonTree_->jet2_ = jets15.at(1).first.p4() * jets15.at(1).second;
        if (jets15.size() > 2) leptonTree_->jet3_ = jets15.at(2).first.p4() * jets15.at(2).second;

        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillMuonFakeRateTree(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // muons
    edm::Handle<edm::View<reco::Muon> > mus_h;
    iEvent.getByLabel(muonsInputTag_, mus_h);
    edm::View<reco::Muon> muonCollection = *(mus_h.product());
    if (muonCollection.size() < 1) return;

    // met
    edm::Handle<edm::View<reco::PFMET> > met_h;
    iEvent.getByLabel(metInputTag_, met_h);
    float met = met_h->front().et();
    if (met > 20.0) return;

    // jets
    edm::Handle<edm::View<reco::PFJet> > jets_h;
    iEvent.getByLabel(jetsInputTag_, jets_h);
    const JetCorrector* corrector = JetCorrector::getJetCorrector(pfJetCorrectorL1FastL2L3_, iSetup);

    // look for a good FO
    edm::View<reco::Muon>::const_iterator it;
    edm::View<reco::Muon>::const_iterator fo;
    unsigned int nFO = 0;
    for (it = muonCollection.begin(); it != muonCollection.end(); ++it) {
        if (it->pt()        < 10.0) continue;
        if (fabs(it->eta()) > 2.4)  continue;
        if (!smurfselections::passMuonFO2011(it, pv_.position())) continue;
        ++nFO;
        fo = it;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {
        fillCommonVariables(iEvent);    

        if (eventPassTrigger("HLT_Mu8_v*"))  leptonTree_->eventSelection_ |= LeptonTree::QCDFakeMu8;
        if (eventPassTrigger("HLT_Mu15_v*")) leptonTree_->eventSelection_ |= LeptonTree::QCDFakeMu15;
        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();

        leptonTree_->leptonSelection_                                                                 |= (LeptonTree::PassMuFO);
        if (smurfselections::passMuonID2011(fo, pv_.position()))         leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
        if (smurfselections::passMuonIso2011(fo, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);

        // jets
        leptonTree_->jet1_ = LorentzVector(0,0,0,0);
        leptonTree_->jet2_ = LorentzVector(0,0,0,0);
        leptonTree_->jet3_ = LorentzVector(0,0,0,0);
        std::vector<std::pair<reco::PFJet, float> > jets30 = smurfselections::goodJets(iEvent, iSetup, jets_h, *fo, corrector, 30.);
        leptonTree_->njets_ = jets30.size();
        std::vector<std::pair<reco::PFJet, float> > jets15 = smurfselections::goodJets(iEvent, iSetup, jets_h, *fo, corrector, 15.);
        if (jets15.size() > 0) leptonTree_->jet1_ = jets15.at(0).first.p4() * jets15.at(0).second;
        if (jets15.size() > 1) leptonTree_->jet2_ = jets15.at(1).first.p4() * jets15.at(1).second;
        if (jets15.size() > 2) leptonTree_->jet3_ = jets15.at(2).first.p4() * jets15.at(2).second;

        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillPhotonTree(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // photons
    edm::Handle<edm::View<reco::Photon> > pho_h;
    iEvent.getByLabel(photonsInputTag_, pho_h);
    edm::View<reco::Photon> photonCollection = *(pho_h.product());
    if (photonCollection.size() < 1) return;

    // met
    edm::Handle<edm::View<reco::PFMET> > met_h;
    iEvent.getByLabel(metInputTag_, met_h);
    float met = met_h->front().et();
    float metPhi = met_h->front().phi();
    
    // jets
    edm::Handle<edm::View<reco::PFJet> > jets_h;
    iEvent.getByLabel(jetsInputTag_, jets_h);
    const JetCorrector* corrector = JetCorrector::getJetCorrector(pfJetCorrectorL1FastL2L3_, iSetup);

    // look for a good photon
    edm::View<reco::Photon>::const_iterator it;
    edm::View<reco::Photon>::const_iterator photon;
    unsigned int nPhotons = 0;
    for (it = photonCollection.begin(); it != photonCollection.end(); ++it) {
        // see https://twiki.cern.ch/twiki/bin/view/CMS/Vgamma2011PhotonID
        float eta = fabs(it->eta());
        float pt = it->pt();
        if (pt < 30.0)                      continue;
        if (eta > 3.0)                      continue;
        if (it->hasPixelSeed())             continue;
        if (it->hadronicOverEm() >= 0.05)   continue;
	//r9 cut
        if (it->r9() < 0.9)   continue;

        if (eta <= 1.479) {
            if (it->sigmaIetaIeta() >= 0.011)   continue;
            if (it->sigmaIetaIeta() <  0.001)   continue;
            if (it->trkSumPtHollowConeDR04()  >= 2.0 + 0.001*pt + 0.0167*rhoIso_)     continue;
            if (it->ecalRecHitSumEtConeDR04() >= 4.2 + 0.006*pt + 0.183*rhoIso_)      continue;
            if (it->hcalTowerSumEtConeDR04()  >= 2.2 + 0.0025*pt + 0.062*rhoIso_)     continue;
        } else {
            if (it->sigmaIetaIeta() >= 0.03)    continue;
            if (it->trkSumPtHollowConeDR04()  >= 2.0 + 0.001*pt + 0.032*rhoIso_)      continue;
            if (it->ecalRecHitSumEtConeDR04() >= 4.2 + 0.006*pt + 0.090*rhoIso_)      continue;
            if (it->hcalTowerSumEtConeDR04()  >= 2.2 + 0.0025*pt + 0.180*rhoIso_)     continue;
        }

        ++nPhotons;
        photon = it;
    }

    // if exactly one good photon
    // then fill the photon tree
    if (nPhotons == 1) {
        fillCommonVariables(iEvent);
        leptonTree_->eventSelection_     = LeptonTree::PhotonSelection;
        leptonTree_->probe_              = photon->p4();
        leptonTree_->qProbe_             = 0;
        leptonTree_->met_                = met;
        leptonTree_->metPhi_             = metPhi;

	const reco::Candidate* phocand = &(*photon);
	std::vector<const reco::Candidate*> phos;
 	phos.push_back(phocand);
	std::pair<double,double> tkmet = smurfselections::trackerMET(phos, 0.1, pfCandCollection_, pv_);
	leptonTree_->trackMet_           = tkmet.first;
	leptonTree_->trackMetPhi_        = tkmet.second;

        // jets
        leptonTree_->jet1_ = LorentzVector(0,0,0,0);
        leptonTree_->jet2_ = LorentzVector(0,0,0,0);
        leptonTree_->jet3_ = LorentzVector(0,0,0,0);
        std::vector<std::pair<reco::PFJet, float> > jets30 = smurfselections::goodJets(iEvent, iSetup, jets_h, *photon, corrector, 30.);
        leptonTree_->njets_ = jets30.size();
        std::vector<std::pair<reco::PFJet, float> > jets15 = smurfselections::goodJets(iEvent, iSetup, jets_h, *photon, corrector, 15.);
        if (jets15.size() > 0) leptonTree_->jet1_ = jets15.at(0).first.p4() * jets15.at(0).second;
        if (jets15.size() > 1) leptonTree_->jet2_ = jets15.at(1).first.p4() * jets15.at(1).second;
        if (jets15.size() > 2) leptonTree_->jet3_ = jets15.at(2).first.p4() * jets15.at(2).second;

	//trigger info
	//here the order matters: need hltPrescale to store the lowest prescale of fired triggers
        if (eventPassTrigger("HLT_Photon20_CaloIdVL_IsoL_v*")) {
	  leptonTree_->eventSelection_ |= LeptonTree::Photon20CaloIdVLIsoL;
	  leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon20_CaloIdVL_IsoL_v*");
	}
	if (eventPassTrigger("HLT_Photon30_CaloIdVL_IsoL_v*")) {
	  leptonTree_->eventSelection_ |= LeptonTree::Photon30CaloIdVLIsoL;
	  leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon30_CaloIdVL_IsoL_v*");
	}
	if (eventPassTrigger("HLT_Photon50_CaloIdVL_IsoL_v*")) {
	  leptonTree_->eventSelection_ |= LeptonTree::Photon50CaloIdVLIsoL;
	  leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon50_CaloIdVL_IsoL_v*");
	}
	if (eventPassTrigger("HLT_Photon75_CaloIdVL_IsoL_v*")) {
	  leptonTree_->eventSelection_ |= LeptonTree::Photon75CaloIdVLIsoL;
	  leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon75_CaloIdVL_IsoL_v*");
	}
	if (eventPassTrigger("HLT_Photon90_CaloIdVL_IsoL_v*")) {
	  leptonTree_->eventSelection_ |= LeptonTree::Photon90CaloIdVLIsoL;
	  leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon90_CaloIdVL_IsoL_v*");
	}

        leptonTree_->tree_->Fill();

    }

}

void LeptonTreeMaker::fillCommonVariables(const edm::Event& iEvent)
{

    // fill
    leptonTree_->run_       = iEvent.id().run()        ;
    leptonTree_->event_     = iEvent.id().event()      ;
    leptonTree_->lumi_      = iEvent.luminosityBlock() ;
    leptonTree_->rho_       = rhoIso_                  ;
    leptonTree_->nvtx_      = vertexCollection_.size()  ;

    edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
    bool bPuInfo=iEvent.getByLabel("addPileupInfo", puInfoH);
    if(bPuInfo) {
      for (std::vector<PileupSummaryInfo>::const_iterator itr = puInfoH->begin(); itr != puInfoH->end(); ++itr ){
	if (itr->getBunchCrossing()== 0) leptonTree_->npu_         = itr->getPU_NumInteractions();
	if (itr->getBunchCrossing()==+1) leptonTree_->npuPlusOne_  = itr->getPU_NumInteractions();
	if (itr->getBunchCrossing()==-1) leptonTree_->npuMinusOne_ = itr->getPU_NumInteractions();
      }
    }

}

bool LeptonTreeMaker::eventPassTrigger(const std::vector<std::string> &trigNames)
{

    for(unsigned int j = 0; j < trigNames.size(); j++) {
        if (eventPassTrigger(trigNames.at(j))) return true;
    }
    return false;
}

bool LeptonTreeMaker::eventPassTrigger(const std::string &trigName)
{

  for(unsigned int i = 0; i<hltConfig_.size(); i++) {
      bool result = triggerResults_->accept(i);
      if (!result) continue;
      TString hltTrigName(hltConfig_.triggerName(i));
      TString pattern(trigName);
      hltTrigName.ToLower();
      pattern.ToLower();
      TRegexp reg(Form("%s", pattern.Data()), true);
      if (hltTrigName.Index(reg) >= 0) return true;
  }
  return false;
}

double LeptonTreeMaker::getTriggerPrescale(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string &trigName)
{

  for(unsigned int i = 0; i<hltConfig_.size(); i++) {
      bool result = triggerResults_->accept(i);
      if (!result) continue;
      TString hltTrigName(hltConfig_.triggerName(i));
      TString pattern(trigName);
      hltTrigName.ToLower();
      pattern.ToLower();
      TRegexp reg(Form("%s", pattern.Data()), true);
      if (hltTrigName.Index(reg) >= 0) return hltConfig_.prescaleValue(iEvent, iSetup, hltConfig_.triggerName(i));
  }
  return 0.;
}

void LeptonTreeMaker::triggerObjects()
{

    int index = -1;
    std::string trigName = "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_L1_SingleEG12_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*";
    for(unsigned int i = 0; i<hltConfig_.size(); i++) {
        TString hltTrigName(hltConfig_.triggerName(i));
        TString pattern(trigName);
        hltTrigName.ToLower();
        pattern.ToLower();
        TRegexp reg(Form("%s", pattern.Data()), true);
        if (hltTrigName.Index(reg) >= 0) {
            index = i;
            break;
        }
    }
    
    std::cout << index << std::endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonTreeMaker);
