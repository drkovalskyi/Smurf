//ryAOD -*- C++ -*-
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
// $Id: LeptonTreeMaker.cc,v 1.21 2012/04/16 13:02:45 dlevans Exp $
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

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
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

// MVAs...
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/MuonIDMVA.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronMVAEstimator.h"

// user include files
#include "Smurf/Core/LeptonTree.h"
#include "Smurf/ProcessingAndSkimming/interface/Selections.h"
#include "Smurf/ProcessingAndSkimming/interface/Utilities.h"

// root include files
#include "TFile.h"
#include "TTree.h"
#include "TRegexp.h"
#include "TString.h"
#include "TPRegexp.h"

//
// constants, enums and typedefs
//

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;
typedef math::XYZTLorentzVectorD LorentzVector;
typedef math::XYZPoint Point;

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
        void fillElectronTagAndProbeTree(const edm::Event& iEvent, const edm::EventSetup &setup,
                const TransientTrackBuilder *ttBuilder, 
                EcalClusterLazyTools *clusterTools);
        void fillMuonTagAndProbeTree(const edm::Event& iEvent, const edm::EventSetup &setup,
                const TransientTrackBuilder *ttBuilder);

        // fake rates
        void fillElectronFakeRateTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
                const TransientTrackBuilder *ttBuilder, 
                EcalClusterLazyTools *clusterTools,
                const JetCorrector *jetCorrector);

        void fillMuonFakeRateTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
                const TransientTrackBuilder *ttBuilder, const JetCorrector *jetCorrector);

        // photons
        void fillPhotonTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
                const JetCorrector *jetCorrector);

        // init dynamic trigger branch values
        void initTriggerBranchValues();

        // common variables
        void fillCommonVariables(const edm::Event& iEvent);
        void fillJets(const edm::Event& iEvent, const edm::EventSetup &iSetup, 
                const reco::Candidate &cand1, const JetCorrector *jetCorrector);

        // find the trigger versions
        // from their names...
        void getTriggerVersions(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                const std::vector<edm::InputTag> &trigNames, std::vector<unsigned int> &versions);

        // did any trigger pass
        void objectMatchTrigger(const edm::Event &iEvent, const edm::EventSetup &iSetup,
                const std::vector<edm::InputTag> &trigNames,
                const trigger::TriggerObjectCollection &allObjects,
                const LorentzVector &obj, std::vector<unsigned int> &prescale);

        bool eventPassTrigger(const std::vector<edm::InputTag> &trigNames);
        bool eventPassTrigger(const std::string &trigName);
        double getTriggerPrescale(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string &trigName);

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

        // pf isolation related
        std::vector<edm::InputTag> isoValInputTags_;

        // jet corrections
        std::string pfJetCorrectorL1FastL2L3_;

        std::string pathToBDTWeights_;

        // trigger names
        std::vector<edm::InputTag> electronFRTriggerNames_;
        std::vector<edm::InputTag> muonFRTriggerNames_;
        std::vector<edm::InputTag> photonTriggerNames_;

        // trigger related
        const trigger::TriggerEvent   *triggerEvent_;
        HLTConfigProvider          hltConfig_;
        std::string                processName_;
        const edm::TriggerResults* triggerResults_;

        std::vector<edm::InputTag> muTriggers_;    
        std::vector<unsigned int> muTriggerPrescalesTag_;
        std::vector<unsigned int> muTriggerPrescalesProbe_;
        std::vector<unsigned int> muTriggerVersions_;
        std::vector<edm::InputTag> eleTriggers_;       
        std::vector<unsigned int> eleTriggerPrescalesTag_;
        std::vector<unsigned int> eleTriggerPrescalesProbe_;
        std::vector<unsigned int> eleTriggerVersions_;

        // electron id related
        ElectronIDMVA           *electronIDMVA_;
        ElectronIDMVA           *electronIDIsoMVA_;
        MuonIDMVA               *muonIDMVA_;
        ElectronMVAEstimator    *egammaIDMVA_;

        // common products
        double rhoJEC_;
        double rhoIso_;
        reco::PFCandidateCollection pfCandCollection_;
        edm::Handle<reco::VertexCollection> vtx_h_;
        edm::Handle<edm::View<reco::PFJet> > jets_h_;
        reco::Vertex pv_;
        JetCorrector *jetCorrector_;

        // random number generators
        float rndm_;

};

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

    isoValInputTags_        =  iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

    pfJetCorrectorL1FastL2L3_   = iConfig.getParameter<std::string>("pfJetCorrectorL1FastL2L3");
    pathToBDTWeights_           = iConfig.getParameter<std::string>("pathToBDTWeights");

    electronFRTriggerNames_ =  iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("electronFRTriggerNames");
    muonFRTriggerNames_     =  iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("muonFRTriggerNames");
    photonTriggerNames_     =  iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("photonTriggerNames");
    muTriggers_             = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("muTriggers");
    eleTriggers_            = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("eleTriggers");

    //
    // set up lepton tree
    //

    leptonFile_ = TFile::Open("leptonTree.root", "RECREATE");
    assert(leptonFile_);
    leptonTree_ = new LeptonTree();
    leptonTree_->CreateTree();
    leptonTree_->tree_->SetDirectory(leptonFile_);

    //
    // add dynamic branches
    // to store trigger
    // --- first branch stores prescale (0 trigger fails, >= 1 prescale of passing trigger)
    // --- second branch stores trigger version number (regexp for _vXX)
    //

    muTriggerPrescalesTag_.reserve(muTriggers_.size());
    muTriggerPrescalesProbe_.reserve(muTriggers_.size());
    muTriggerVersions_.reserve(muTriggers_.size());
    for (unsigned int i = 0; i < muTriggers_.size(); ++i) {
        const char *trigName = muTriggers_[i].process().c_str();
        leptonTree_->tree_->Branch(Form("%s_tag", trigName), &muTriggerPrescalesTag_[i], Form("%s_tag/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_probe", trigName), &muTriggerPrescalesProbe_[i], Form("%s_probe/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_version", trigName), &muTriggerVersions_[i], Form("%s_version/i", trigName));
    }

    eleTriggerPrescalesTag_.reserve(eleTriggers_.size());
    eleTriggerPrescalesProbe_.reserve(eleTriggers_.size());
    eleTriggerVersions_.reserve(eleTriggers_.size());
    for (unsigned int i = 0; i < eleTriggers_.size(); ++i) {
        const char *trigName = eleTriggers_[i].process().c_str();
        leptonTree_->tree_->Branch(Form("%s_tag", trigName), &eleTriggerPrescalesTag_[i], Form("%s_tag/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_probe", trigName), &eleTriggerPrescalesProbe_[i], Form("%s_probe/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_version", trigName), &eleTriggerVersions_[i], Form("%s_version/i", trigName));
    }

    //
    // set up trigger
    //

    processName_ = "";

    //
    // set up electron id
    //

    electronIDMVA_ = new ElectronIDMVA();
    electronIDMVA_->Initialize("BDTG method",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml",                
            ElectronIDMVA::kWithIPInfo);

    electronIDIsoMVA_ = new ElectronIDMVA();
    electronIDIsoMVA_->Initialize("BDTG method",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml",
            ElectronIDMVA::kIDIsoCombined);

    muonIDMVA_ = new MuonIDMVA();
    muonIDMVA_->Initialize("BDTG method", 
            pathToBDTWeights_+"/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml",
            pathToBDTWeights_+"/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml",
            MuonIDMVA::kIDIsoCombinedDetIso);

    std::vector<std::string> myManualCatWeigthsTrig;
    myManualCatWeigthsTrig.push_back(pathToBDTWeights_+"/Electrons_BDTG_TrigV0_Cat1.weights.xml");
    myManualCatWeigthsTrig.push_back(pathToBDTWeights_+"/Electrons_BDTG_TrigV0_Cat2.weights.xml");
    myManualCatWeigthsTrig.push_back(pathToBDTWeights_+"/Electrons_BDTG_TrigV0_Cat3.weights.xml");
    myManualCatWeigthsTrig.push_back(pathToBDTWeights_+"/Electrons_BDTG_TrigV0_Cat4.weights.xml");
    myManualCatWeigthsTrig.push_back(pathToBDTWeights_+"/Electrons_BDTG_TrigV0_Cat5.weights.xml");
    myManualCatWeigthsTrig.push_back(pathToBDTWeights_+"/Electrons_BDTG_TrigV0_Cat6.weights.xml");
    egammaIDMVA_ = new ElectronMVAEstimator();
    egammaIDMVA_->initialize("BDT",
            ElectronMVAEstimator::kTrig,
            true,
            myManualCatWeigthsTrig);

}

LeptonTreeMaker::~LeptonTreeMaker()
{

    //
    // tidy up
    //

    if (electronIDMVA_)     delete electronIDMVA_;
    if (electronIDIsoMVA_)  delete electronIDIsoMVA_;
    if (egammaIDMVA_)       delete egammaIDMVA_;
    if (muonIDMVA_)         delete muonIDMVA_;

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
    edm::Handle<trigger::TriggerEvent> triggerEvent_h_;
    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEvent_h_);
    if (!triggerEvent_h_.isValid()) {
        throw cms::Exception("[LeptonTreeMaker][produce] Error getting TriggerEvent product from Event!");
    }
    triggerEvent_ = triggerEvent_h_.product();

    // set process name if not set
    if (processName_ == "") {

        processName_ = triggerEvent_h_.provenance()->processName();
        bool changed(true);
        hltConfig_.init(iEvent.getRun(), iSetup, processName_, changed);

        // dump save tags for specific paths
        smurfutilities::DumpSaveTags("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Ele27_WP80_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_IsoMu17_Mu8_v", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_IsoMu24_eta2p1_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu17_Mu8_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu17_TkMu8_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_IsoMu30_eta2p1_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Ele27_WP80_v*", hltConfig_);

    }

    // trigger results
    edm::Handle<edm::TriggerResults> triggerResults_h;
    iEvent.getByLabel(edm::InputTag("TriggerResults", "", processName_), triggerResults_h);
    if (!triggerResults_h.isValid()) {
        throw cms::Exception("LeptonTreeMaker][produce] Error getting TriggerResults product from Event!");
    }
    triggerResults_ = triggerResults_h.product();

    // sanity check
    assert(triggerResults_->size() == hltConfig_.size());

    //  
    // set up tools
    //

    // random numbers
   edm::Service<edm::RandomNumberGenerator> rng;
   if (!rng.isAvailable()) {
     throw cms::Exception("Configuration")
       << "LeptonTreeMaker requires the RandomNumberGeneratorService\n"
          "which is not present in the configuration file.  You must add the service\n"
          "in the configuration file or remove the modules that require it.";
   }
   CLHEP::RandFlat rndgen(rng->getEngine());
    rndm_ = rndgen.fire();

    // unbiased revertexing
    edm::ESHandle<TransientTrackBuilder> ttBuilder_h;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder_h);
    const TransientTrackBuilder *ttBuilder = ttBuilder_h.product();

    // ecal tools
    EcalClusterLazyTools *clusterTools = new EcalClusterLazyTools(iEvent, iSetup, 
            edm::InputTag("reducedEcalRecHitsEB"), 
            edm::InputTag("reducedEcalRecHitsEE"));

    // rho for jec
    edm::Handle<double> rhoJEC_h;
    iEvent.getByLabel(rhoJECInputTag, rhoJEC_h);
    rhoJEC_ = *(rhoJEC_h.product());

    // jets
    iEvent.getByLabel(jetsInputTag_, jets_h_);
    const JetCorrector *jetCorrector = JetCorrector::getJetCorrector(pfJetCorrectorL1FastL2L3_, iSetup);

    // rho for isolation
    edm::Handle<double> rhoIso_h;
    iEvent.getByLabel(rhoIsoInputTag, rhoIso_h);
    rhoIso_ = *(rhoIso_h.product());

    // pf candidates
    edm::Handle<reco::PFCandidateCollection> pfCand_h;
    iEvent.getByLabel(pfCandsInputTag_, pfCand_h);
    pfCandCollection_ = *(pfCand_h.product());

    // vertices
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h_);
    const reco::VertexCollection vertexCollection = *(vtx_h_.product());
    if (vertexCollection.size() > 0) pv_ = vertexCollection.at(0);
    else return;

    //
    // decide what to do based on trigger information
    //

    bool passElectronFRTrigger  = eventPassTrigger(electronFRTriggerNames_);
    bool passMuonFRTrigger      = eventPassTrigger(muonFRTriggerNames_);
    bool passPhotonTrigger      = eventPassTrigger(photonTriggerNames_);

    // fake rates
    if (passElectronFRTrigger)  fillElectronFakeRateTree(iEvent, iSetup, ttBuilder, clusterTools, jetCorrector);
    if (passMuonFRTrigger)      fillMuonFakeRateTree(iEvent, iSetup, ttBuilder, jetCorrector);

    // efficiency
    fillElectronTagAndProbeTree(iEvent, iSetup, ttBuilder, clusterTools);
    fillMuonTagAndProbeTree(iEvent, iSetup, ttBuilder);

    // photons
    if (passPhotonTrigger)      fillPhotonTree(iEvent, iSetup, jetCorrector);

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

    //
    // configure HLT
    //

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

void LeptonTreeMaker::fillElectronTagAndProbeTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
        const TransientTrackBuilder *ttBuilder,
        EcalClusterLazyTools *clusterTools)
{

    leptonTree_->InitVariables();
    initTriggerBranchValues();

    // electrons
    edm::Handle<reco::GsfElectronCollection> els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    reco::GsfElectronCollection electronCollection = *(els_h.product());
    if (electronCollection.size() < 2) return;

    // iso deposits
    // for pf iso
    IsoDepositVals isoVals(isoValInputTags_.size());
    for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
        iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
    }

    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    iEvent.getByLabel(conversionsInputTag_, conversions_h);

    // beamspot
    edm::Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
    const reco::BeamSpot &thebs = *(beamspot_h.product());

    // triggers
    const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();

    // look for tag and probe
    unsigned int nEle = els_h->size();
    for (unsigned int itag = 0; itag < nEle; ++itag) {

        // look for good tag
        reco::GsfElectronRef tag(els_h, itag);
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.5) continue;
        float mvaValue = electronIDMVA_->MVAValue(&*tag, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
        if (!smurfselections::passElectronID2011(tag, pv_, thebs.position(), conversions_h, mvaValue)) continue;
        if (!smurfselections::passElectronIso2011(tag, pfCandCollection_, pv_)) continue;

        // if real data
        // then store tag trigger matching
        if (iEvent.isRealData()) 
            objectMatchTrigger(iEvent, iSetup, eleTriggers_, allObjects, tag->p4(), eleTriggerPrescalesTag_);

        for (unsigned int iprobe = 0; iprobe < nEle; ++iprobe) {

            if (itag == iprobe)             continue;
            reco::GsfElectronRef probe(els_h, iprobe);

            // look for good probe
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

            mvaValue                  = electronIDMVA_->MVAValue(&*probe, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
            float mvaIdIsoValue       = electronIDIsoMVA_->MVAValue(&*probe, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
            float egammaMVAValue      = egammaIDMVA_->mvaValue(*probe, pv_, *ttBuilder, *clusterTools);
            leptonTree_->electronMVA_       = mvaValue;
            leptonTree_->electronIDIsoMVA_  = mvaIdIsoValue;
            leptonTree_->egammaMVA_         = egammaMVAValue;

            if (smurfselections::passElectronFO2011(probe, pv_, thebs.position(), conversions_h))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
            if (smurfselections::passElectronID2011(probe, pv_, thebs.position(), conversions_h, mvaValue))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
            if (smurfselections::passElectronIso2011(probe, pfCandCollection_, pv_)) 
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);

            // cut based electron id
            double iso_ch =  (*(isoVals)[0])[probe];
            double iso_em = (*(isoVals)[1])[probe];
            double iso_nh = (*(isoVals)[2])[probe];
            leptonTree_->vetoId_    = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::VETO, probe, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);
            leptonTree_->looseId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE, probe, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);
            leptonTree_->mediumId_  = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::MEDIUM, probe, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);
            leptonTree_->tightId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::TIGHT, probe, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);

            // input variables
            leptonTree_->pfmva_     = probe->mva();
            leptonTree_->sceta_     = probe->superCluster()->eta();
            leptonTree_->scenergy_  = probe->superCluster()->energy();
            leptonTree_->chargesAgree_ = smurfselections::threeChargesAgree(*probe);
            leptonTree_->detain_    = probe->deltaEtaSuperClusterTrackAtVtx();
            leptonTree_->dphiin_    = probe->deltaPhiSuperClusterTrackAtVtx();
            leptonTree_->sieie_     = probe->sigmaIetaIeta();
            leptonTree_->hoe_       = probe->hadronicOverEm();
            #ifdef RELEASE_52X
            leptonTree_->hoetow_     = probe->hcalOverEcalBc();
            #endif
            leptonTree_->ooemoop_   = (1.0/probe->ecalEnergy() - probe->eSuperClusterOverP()/probe->ecalEnergy());
            leptonTree_->d0vtx_     = probe->gsfTrack()->dxy(pv_.position());
            leptonTree_->dzvtx_     = probe->gsfTrack()->dz(pv_.position());
            leptonTree_->vfitprob_  = ConversionTools::hasMatchedConversion(*probe, conversions_h, thebs.position());
            leptonTree_->mhit_      = probe->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
            leptonTree_->ecaliso_   = probe->dr03EcalRecHitSumEt();
            leptonTree_->hcaliso_   = probe->dr03HcalTowerSumEt();
            leptonTree_->trkiso_    = probe->dr03TkSumPt();
            leptonTree_->pfchiso_   = iso_ch;
            leptonTree_->pfemiso_   = iso_em;
            leptonTree_->pfnhiso_   = iso_nh;
            leptonTree_->pfaeff_    = smurfselections::GetEGammaEffectiveArea(probe->superCluster()->eta());

            // probe trigger matching
            if (iEvent.isRealData()) {
                objectMatchTrigger(iEvent, iSetup, eleTriggers_, allObjects, probe->p4(), eleTriggerPrescalesProbe_);
                getTriggerVersions(iEvent, iSetup, eleTriggers_, eleTriggerVersions_);
            }

            // fill the tree
            leptonTree_->tree_->Fill(); 

        }

    }

}

void LeptonTreeMaker::fillMuonTagAndProbeTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
        const TransientTrackBuilder *ttBuilder)
{

    leptonTree_->InitVariables();
    initTriggerBranchValues();

    // muons
    edm::Handle<edm::View<reco::Muon> > mus_h;
    iEvent.getByLabel(muonsInputTag_, mus_h);
    edm::View<reco::Muon> muonCollection = *(mus_h.product());
    if (muonCollection.size() < 2) return;

    // triggers
    const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();

    // look for tag and probe
    edm::View<reco::Muon>::const_iterator tag;
    edm::View<reco::Muon>::const_iterator probe;

    for (tag = muonCollection.begin(); tag != muonCollection.end(); ++tag) {

        // look for good tag
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.4) continue;
        if (!smurfselections::passMuonID2011(tag, pv_))                     continue;
        if (!smurfselections::passMuonIso2011(tag, pfCandCollection_, pv_)) continue;

        // if real data
        // then store tag trigger matching
        if (iEvent.isRealData()) 
            objectMatchTrigger(iEvent, iSetup, muTriggers_, allObjects, tag->p4(), muTriggerPrescalesTag_);

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

            float mvaValue          = muonIDMVA_->MVAValue(&*probe, pv_, ttBuilder, rhoIso_, false); 
            leptonTree_->muonMVA_   = mvaValue;

            if (smurfselections::passMuonFO2011(probe, pv_))                     leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
            if (smurfselections::passMuonID2011(probe, pv_))                     leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
            if (smurfselections::passMuonIso2011(probe, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);

            // probe trigger matching
            if (iEvent.isRealData()) {
                objectMatchTrigger(iEvent, iSetup, muTriggers_, allObjects, probe->p4(), muTriggerPrescalesProbe_);
                getTriggerVersions(iEvent, iSetup, muTriggers_, muTriggerVersions_);
            }

            leptonTree_->tree_->Fill();

        }

    }

}

void LeptonTreeMaker::fillElectronFakeRateTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const TransientTrackBuilder *ttBuilder,
        EcalClusterLazyTools *clusterTools,
        const JetCorrector *jetCorrector)
{

    leptonTree_->InitVariables();

    // electrons
    edm::Handle<reco::GsfElectronCollection> els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    reco::GsfElectronCollection electronCollection = *(els_h.product());
    if (electronCollection.size() < 1) return;

    // iso deposits
    // for pf iso
    IsoDepositVals isoVals(isoValInputTags_.size());
    for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
        iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
    }

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

    // look for a good FO
    unsigned int nEle = els_h->size();
    unsigned int nFO = 0;
    reco::GsfElectronRef fo;
    for (unsigned int iele = 0; iele < nEle; ++iele) {
        reco::GsfElectronRef ele(els_h, iele);
        if (ele->pt()        < 10.0) continue;
        if (fabs(ele->eta()) > 2.4)  continue;
        if (!smurfselections::passElectronFO2011(ele, pv_, thebs.position(), conversions_h)) continue;
        ++nFO;
        fo = ele;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {

        fillCommonVariables(iEvent);
        fillJets(iEvent, iSetup, *fo, jetCorrector);

        // 2011 triggers
        if (eventPassTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v*"))       leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle8;
        if (eventPassTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v*"))      leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle17;
        if (eventPassTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*")) leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle8Jet40;

        // additional triggers used in 2012
        // not covered by the above
        if (eventPassTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"))         leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle8;
        if (eventPassTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*"))   leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle8Jet30;
        if (eventPassTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"))        leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle17;
        if (eventPassTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*"))  leptonTree_->eventSelection_ |= LeptonTree::QCDFakeEle17Jet30;

        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();

        float mvaValue            = electronIDMVA_->MVAValue(&*fo, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
        float mvaIdIsoValue       = electronIDIsoMVA_->MVAValue(&*fo, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIso_);
        float egammaMVAValue      = egammaIDMVA_->mvaValue(*fo, pv_, *ttBuilder, *clusterTools);
        leptonTree_->electronMVA_       = mvaValue;
        leptonTree_->electronIDIsoMVA_  = mvaIdIsoValue;
        leptonTree_->egammaMVA_         = egammaMVAValue;

        leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
        if (smurfselections::passElectronID2011(fo, pv_, thebs.position(), conversions_h, mvaValue))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
        if (smurfselections::passElectronIso2011(fo, pfCandCollection_, pv_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);

        // cut based electron id
        double iso_ch = (*(isoVals)[0])[fo];
        double iso_em = (*(isoVals)[1])[fo];
        double iso_nh = (*(isoVals)[2])[fo];
        leptonTree_->vetoId_    = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::VETO, fo, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);
        leptonTree_->looseId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE, fo, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);
        leptonTree_->mediumId_  = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::MEDIUM, fo, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);
        leptonTree_->tightId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::TIGHT, fo, conversions_h, thebs, vtx_h_, iso_ch, iso_em, iso_nh, rhoIso_);

        leptonTree_->pfmva_     = fo->mva();
        leptonTree_->sceta_     = fo->superCluster()->eta();
        leptonTree_->scenergy_  = fo->superCluster()->energy();
        leptonTree_->chargesAgree_ = smurfselections::threeChargesAgree(*fo);
        leptonTree_->detain_    = fo->deltaEtaSuperClusterTrackAtVtx();
        leptonTree_->dphiin_    = fo->deltaPhiSuperClusterTrackAtVtx();
        leptonTree_->sieie_     = fo->sigmaIetaIeta();
        leptonTree_->hoe_       = fo->hadronicOverEm();
        #ifdef RELEASE_52X
        leptonTree_->hoetow_     = fo->hcalOverEcalBc();
        #endif
        leptonTree_->ooemoop_   = (1.0/fo->ecalEnergy() - fo->eSuperClusterOverP()/fo->ecalEnergy());
        leptonTree_->d0vtx_     = fo->gsfTrack()->dxy(pv_.position());
        leptonTree_->dzvtx_     = fo->gsfTrack()->dz(pv_.position());
        leptonTree_->vfitprob_  = ConversionTools::hasMatchedConversion(*fo, conversions_h, thebs.position());
        leptonTree_->mhit_      = fo->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
        leptonTree_->ecaliso_   = fo->dr03EcalRecHitSumEt();
        leptonTree_->hcaliso_   = fo->dr03HcalTowerSumEt();
        leptonTree_->trkiso_    = fo->dr03TkSumPt();
        leptonTree_->pfchiso_   = iso_ch;
        leptonTree_->pfemiso_   = iso_nh;
        leptonTree_->pfnhiso_   = iso_em;
        leptonTree_->pfaeff_    = smurfselections::GetEGammaEffectiveArea(fo->superCluster()->eta());

        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillMuonFakeRateTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const TransientTrackBuilder *ttBuilder, const JetCorrector *jetCorrector)
{

    leptonTree_->InitVariables();

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

    // look for a good FO
    edm::View<reco::Muon>::const_iterator it;
    edm::View<reco::Muon>::const_iterator fo;
    unsigned int nFO = 0;
    for (it = muonCollection.begin(); it != muonCollection.end(); ++it) {
        if (it->pt()        < 10.0) continue;
        if (fabs(it->eta()) > 2.4)  continue;
        if (!smurfselections::passMuonFO2011(it, pv_)) continue;
        ++nFO;
        fo = it;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {

        fillCommonVariables(iEvent);    
        fillJets(iEvent, iSetup, *fo, jetCorrector);

        // 2011 triggers
        if (eventPassTrigger("HLT_Mu8_v*"))  leptonTree_->eventSelection_ |= LeptonTree::QCDFakeMu8;
        if (eventPassTrigger("HLT_Mu15_v*")) leptonTree_->eventSelection_ |= LeptonTree::QCDFakeMu15;

        // additional triggers used in 2012
        // not covered by the above
        if (eventPassTrigger("HLT_Mu17_v*")) leptonTree_->eventSelection_ |= LeptonTree::QCDFakeMu17;

        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();

        float mvaValue          = muonIDMVA_->MVAValue(&*fo, pv_, ttBuilder, rhoIso_, false);
        leptonTree_->muonMVA_   = mvaValue;

        leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
        if (smurfselections::passMuonID2011(fo, pv_))                     leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
        if (smurfselections::passMuonIso2011(fo, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);

        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillPhotonTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const JetCorrector *jetCorrector)
{

    leptonTree_->InitVariables();

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

    // look for a good photon
    edm::View<reco::Photon>::const_iterator it;
    edm::View<reco::Photon>::const_iterator photon;
    unsigned int nPhotons = 0;
    for (it = photonCollection.begin(); it != photonCollection.end(); ++it) {

        if (!smurfselections::passPhotonSelection2011(it, rhoIso_)) continue;
        ++nPhotons;
        photon = it;

    }

    // if exactly one good photon
    // then fill the photon tree
    if (nPhotons == 1) {

        fillCommonVariables(iEvent);
        fillJets(iEvent, iSetup, *photon, jetCorrector);

        leptonTree_->eventSelection_     = LeptonTree::PhotonSelection;
        leptonTree_->probe_              = photon->p4();
        leptonTree_->qProbe_             = 0;
        leptonTree_->met_                = met;
        leptonTree_->metPhi_             = metPhi;
        leptonTree_->sumet_              = met_h->front().sumEt();
        leptonTree_->metSig_             = met_h->front().significance();

        const reco::Candidate* phocand = &(*photon);
        std::vector<const reco::Candidate*> phos;
        phos.push_back(phocand);
        std::pair<double,double> tkmet = smurfselections::trackerMET(phos, 0.1, pfCandCollection_, pv_);
        leptonTree_->trackMet_           = tkmet.first;
        leptonTree_->trackMetPhi_        = tkmet.second;

        // trigger info 2011
        // here the order matters: need hltPrescale to store the lowest prescale of fired triggers
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

        // trigger info 2012
        if (eventPassTrigger("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v*")) {
            leptonTree_->eventSelection_ |= LeptonTree::Photon22R9Id90HE10Iso40EB;
            leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v*");
        }
        if (eventPassTrigger("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*")) {
            leptonTree_->eventSelection_ |= LeptonTree::Photon36R9Id90HE10Iso40EB;
            leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*");
        }
        if (eventPassTrigger("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*")) {
            leptonTree_->eventSelection_ |= LeptonTree::Photon50R9Id90HE10Iso40EB;
            leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*");
        }
        if (eventPassTrigger("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*")) {
            leptonTree_->eventSelection_ |= LeptonTree::Photon75R9Id90HE10Iso40EB;
            leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*");
        }
        if (eventPassTrigger("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*")) {
            leptonTree_->eventSelection_ |= LeptonTree::Photon90R9Id90HE10Iso40EB;
            leptonTree_->hltPrescale_ = getTriggerPrescale(iEvent, iSetup, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*");
        }

        leptonTree_->tree_->Fill();

    }

}

void LeptonTreeMaker::initTriggerBranchValues()
{
    // init the dynamic trigger branches
    for (unsigned int i = 0; i < muTriggers_.size(); ++i) {
        muTriggerPrescalesTag_[i] = 0;
        muTriggerPrescalesProbe_[i] = 0;
        muTriggerVersions_[i] = 0;
    }

    for (unsigned int i = 0; i < eleTriggers_.size(); ++i) {
        eleTriggerPrescalesTag_[i] = 0;
        eleTriggerPrescalesProbe_[i] = 0;
        eleTriggerVersions_[i] = 0;
    }
}

void LeptonTreeMaker::fillCommonVariables(const edm::Event& iEvent)
{

    // fill
    leptonTree_->run_       = iEvent.id().run()        ;
    leptonTree_->event_     = iEvent.id().event()      ;
    leptonTree_->lumi_      = iEvent.luminosityBlock() ;
    leptonTree_->rnd_       = rndm_                    ;
    leptonTree_->rho_       = rhoIso_                  ;
    leptonTree_->nvtx_      = vtx_h_->size()  ;
    leptonTree_->scale1fb_  = 1.0;

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

void LeptonTreeMaker::fillJets(const edm::Event& iEvent, const edm::EventSetup &iSetup, 
            const reco::Candidate &cand1, const JetCorrector *jetCorrector)
{

    leptonTree_->jet1_ = LorentzVector(0,0,0,0);
    leptonTree_->jet2_ = LorentzVector(0,0,0,0);
    leptonTree_->jet3_ = LorentzVector(0,0,0,0);
    std::vector<std::pair<reco::PFJet, float> > jets30 = smurfselections::goodJets(iEvent, iSetup, jets_h_, cand1, jetCorrector, 30.);
    leptonTree_->njets_ = jets30.size();
    std::vector<std::pair<reco::PFJet, float> > jets15 = smurfselections::goodJets(iEvent, iSetup, jets_h_, cand1, jetCorrector, 15.);
    if (jets15.size() > 0) leptonTree_->jet1_ = jets15.at(0).first.p4() * jets15.at(0).second;
    if (jets15.size() > 1) leptonTree_->jet2_ = jets15.at(1).first.p4() * jets15.at(1).second;
    if (jets15.size() > 2) leptonTree_->jet3_ = jets15.at(2).first.p4() * jets15.at(2).second;

}

void LeptonTreeMaker::getTriggerVersions(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
        const std::vector<edm::InputTag> &trigNames, std::vector<unsigned int> &versions)
{

    TPRegexp re("._v(.*)");

    // loop on trigger names
    for (unsigned int t = 0; t < trigNames.size(); ++t) {
        TString trigName = trigNames[t].label();

        // loop on triggers in menu
        for (unsigned int i = 0; i < hltConfig_.size(); i++) {

            // get name of ith trigger
            TString hltTrigName(hltConfig_.triggerName(i));
            hltTrigName.ToLower(); 

            // test if it matches this trigger name
            // with any version
            TString pattern(trigName);
            pattern.ToLower();
            TRegexp reg(Form("%s", pattern.Data()), true);

            // if trigger matches
            // then extract version number
            if (hltTrigName.Index(reg) >= 0) {

                TObjArray *substrArr = re.MatchS(hltTrigName);
                if (substrArr->GetLast() == 1) {
                    versions[t] = ((TObjString*)substrArr->At(1))->GetString().Atoi();
                } else {
                    versions[t] = 0;
                }

            }
        }
            
    }

}

void LeptonTreeMaker::objectMatchTrigger(const edm::Event &iEvent, const edm::EventSetup &iSetup,
                const std::vector<edm::InputTag> &trigNames,
                const trigger::TriggerObjectCollection &allObjects,
                const LorentzVector &obj, std::vector<unsigned int> &prescale)
{

    for (unsigned int t = 0; t < trigNames.size(); ++t) {
        prescale[t] = smurfutilities::MatchTriggerObject(iEvent, iSetup,
            trigNames[t].label(), trigNames[t].instance(),
            processName_, hltConfig_, triggerResults_, triggerEvent_, allObjects, obj);
    } 

}

bool LeptonTreeMaker::eventPassTrigger(const std::vector<edm::InputTag> &trigNames)
{

    for(unsigned int j = 0; j < trigNames.size(); j++) {
        if (eventPassTrigger(trigNames[j].label())) return true;
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

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonTreeMaker);
