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
// $Id$
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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"

// user include files
#include "Smurf/Core/LeptonTree.h"

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
        void fillElectronTagAndProbeTree(const edm::Event& iEvent);
        void fillMuonTagAndProbeTree(const edm::Event& iEvent);

        // fake rates
        void fillElectronFakeRateTree(const edm::Event& iEvent);
        void fillMuonFakeRateTree(const edm::Event& iEvent);

        // photons
        void fillPhotonTree(const edm::Event& iEvent);

        // common variables
        void fillCommonVariables(const edm::Event& iEvent);

        // did any trigger pass
        bool eventPassTrigger(const std::vector<std::string> &trigNames);

        // ----------member data ---------------------------

        // lepton tree
        TFile       *leptonFile_;
        LeptonTree  *leptonTree_;

        // input tags
        edm::InputTag electronsInputTag_;
        edm::InputTag muonsInputTag_;
        edm::InputTag photonsInputTag_;
        edm::InputTag primaryVertexInputTag_;
        edm::InputTag rhoInputTag_;

        // trigger names
        std::vector<std::string> electronTPTriggerNames_;
        std::vector<std::string> muonTPTriggerNames_;
        std::vector<std::string> electronFRTriggerNames_;
        std::vector<std::string> muonFRTriggerNames_;
        std::vector<std::string> photonTriggerNames_;

        // trigger related
        HLTConfigProvider       hltConfig_;
        std::string             processName_;

        // electron id related
        ElectronIDMVA *electronIDMVA_;

};

//
// constants, enums and typedefs
//

typedef math::XYZTLorentzVectorD LorentzVector;

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
    rhoInputTag_            =  iConfig.getParameter<edm::InputTag>("rhoInputTag");

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

    // sanity check
    assert(triggerResults_h->size() == hltConfig_.size());

    //
    // decide what to do based on trigger information
    //

    bool passElectronTPTrigger  = eventPassTrigger(electronTPTriggerNames_);
    bool passMuonTPTrigger      = eventPassTrigger(muonTPTriggerNames_);
    bool passElectronFRTrigger  = eventPassTrigger(electronFRTriggerNames_);
    bool passMuonFRTrigger      = eventPassTrigger(muonFRTriggerNames_);
    bool passPhotonTrigger      = eventPassTrigger(photonTriggerNames_);

    // fake rates
    if (passElectronFRTrigger)  fillElectronFakeRateTree(iEvent);
    if (passMuonFRTrigger)      fillElectronFakeRateTree(iEvent);

    // efficiency
    if (passElectronTPTrigger)  fillElectronTagAndProbeTree(iEvent);
    if (passMuonTPTrigger)      fillElectronTagAndProbeTree(iEvent);

    // photons
    if (passPhotonTrigger)      fillPhotonTree(iEvent);

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

void LeptonTreeMaker::fillElectronTagAndProbeTree(const edm::Event& iEvent)
{

    // electrons
    edm::Handle<edm::View<reco::GsfElectron> > els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    edm::View<reco::GsfElectron> electronCollection = *(els_h.product());
    if (electronCollection.size() < 2) return;

    // look for tag and probe
    edm::View<reco::GsfElectron>::const_iterator tag;
    edm::View<reco::GsfElectron>::const_iterator probe;

    for (tag = electronCollection.begin(); tag != electronCollection.end(); ++tag) {

        // look for good tag
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.4) continue;
        


        for (probe = electronCollection.begin(); probe != electronCollection.end(); ++probe) {

            // look for good probe
            if (tag == probe)            continue;
            if (probe->pt() < 10.0)      continue;
            if (fabs(probe->eta()) >2.4) continue;

            // construct tag and probe mass
            LorentzVector p4 = tag->p4() + probe->p4();

            // fill the tree
            fillCommonVariables(iEvent);
            leptonTree_->eventSelection_     = LeptonTree::ZeeTagAndProbe;
            leptonTree_->pt_                 = probe->pt();
            leptonTree_->eta_                = probe->eta();
            leptonTree_->phi_                = probe->phi();
            leptonTree_->tagAndProbeMass_    = p4.M();
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

        for (probe = muonCollection.begin(); probe != muonCollection.end(); ++probe) {

            // look for good probe
            if (tag == probe)            continue;
            if (probe->pt() < 10.0)      continue;
            if (fabs(probe->eta()) >2.4) continue;

            // construct tag and probe mass
            LorentzVector p4 = tag->p4() + probe->p4();

            // fill the tree
            fillCommonVariables(iEvent);
            leptonTree_->eventSelection_     = LeptonTree::ZmmTagAndProbe;
            leptonTree_->pt_                 = probe->pt();
            leptonTree_->eta_                = probe->eta();
            leptonTree_->phi_                = probe->phi();
            leptonTree_->tagAndProbeMass_    = p4.M();
            leptonTree_->tree_->Fill();

        }

    }

}

void LeptonTreeMaker::fillElectronFakeRateTree(const edm::Event& iEvent)
{

    // electrons
    edm::Handle<edm::View<reco::GsfElectron> > els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    edm::View<reco::GsfElectron> electronCollection = *(els_h.product());
    if (electronCollection.size() < 1) return;

    // look for a good FO
    edm::View<reco::GsfElectron>::const_iterator it;
    edm::View<reco::GsfElectron>::const_iterator fo;
    unsigned int nFO = 0;
    for (it = electronCollection.begin(); it != electronCollection.end(); ++it) {
        if (it->pt()        < 10.0) continue;
        if (fabs(it->eta()) > 2.4)  continue;
        ++nFO;
        fo = it;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {
        fillCommonVariables(iEvent);
        leptonTree_->eventSelection_     = LeptonTree::QCDFakeEle8;   // FIXME
        leptonTree_->pt_                 = fo->pt();
        leptonTree_->eta_                = fo->eta();
        leptonTree_->phi_                = fo->phi();
        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillMuonFakeRateTree(const edm::Event& iEvent)
{

    // muons
    edm::Handle<edm::View<reco::Muon> > mus_h;
    iEvent.getByLabel(muonsInputTag_, mus_h);
    edm::View<reco::Muon> muonCollection = *(mus_h.product());
    if (muonCollection.size() < 1) return;

    // look for a good FO
    edm::View<reco::Muon>::const_iterator it;
    edm::View<reco::Muon>::const_iterator fo;
    unsigned int nFO = 0;
    for (it = muonCollection.begin(); it != muonCollection.end(); ++it) {
        if (it->pt()        < 10.0) continue;
        if (fabs(it->eta()) > 2.4)  continue;
        ++nFO;
        fo = it;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {
        fillCommonVariables(iEvent);    
        leptonTree_->eventSelection_     = LeptonTree::QCDFakeMu8;   // FIXME
        leptonTree_->pt_                 = fo->pt();
        leptonTree_->eta_                = fo->eta();
        leptonTree_->phi_                = fo->phi();
        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillPhotonTree(const edm::Event& iEvent)
{

    // photons
    edm::Handle<edm::View<reco::Photon> > pho_h;
    iEvent.getByLabel(photonsInputTag_, pho_h);
    edm::View<reco::Photon> photonCollection = *(pho_h.product());
    if (photonCollection.size() < 1) return;

    // look for a good photon
    edm::View<reco::Photon>::const_iterator it;
    edm::View<reco::Photon>::const_iterator photon;
    unsigned int nPhotons = 0;
    for (it = photonCollection.begin(); it != photonCollection.end(); ++it) {
        if (it->pt()        < 30.0) continue;
        if (fabs(it->eta()) > 3.0)  continue;
        ++nPhotons;
        photon = it;
    }

    // if exactly one good photon
    // then fill the photon tree
    if (nPhotons == 1) {
        fillCommonVariables(iEvent);
        leptonTree_->eventSelection_     = LeptonTree::PhotonSelection;
        leptonTree_->pt_                 = photon->pt();
        leptonTree_->eta_                = photon->eta();
        leptonTree_->phi_                = photon->phi();
        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillCommonVariables(const edm::Event& iEvent)
{

    // vertices
    edm::Handle<edm::View<reco::Vertex> > vtx_h;
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h);
    edm::View<reco::Vertex> vertexCollection = *(vtx_h.product());

    // rho
    edm::Handle<double> rho_h;
    iEvent.getByLabel(rhoInputTag_, rho_h);
    double rho = *(rho_h.product());

    // fill
    leptonTree_->run_       = iEvent.id().run()        ;
    leptonTree_->event_     = iEvent.id().event()      ;
    leptonTree_->lumi_      = iEvent.luminosityBlock() ;
    leptonTree_->rho_       = rho                      ;
    leptonTree_->nvtx_      = vertexCollection.size()  ;

}

bool LeptonTreeMaker::eventPassTrigger(const std::vector<std::string> &trigNames)
{

  for(unsigned int i = 0; i<hltConfig_.size(); i++) {
    for(unsigned int j = 0; j < trigNames.size(); j++) {
      TString hltTrigName(hltConfig_.triggerName(i));
      TString pattern(trigNames.at(j));
      hltTrigName.ToLower();
      pattern.ToLower();
      TRegexp reg(Form("%s", pattern.Data()), true);
      if (hltTrigName.Index(reg) >= 0) return true;
    }
  }

    return false;

}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonTreeMaker);
