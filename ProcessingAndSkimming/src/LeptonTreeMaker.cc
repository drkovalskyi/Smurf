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
// $Id: LeptonTreeMaker.cc,v 1.46 2012/09/13 14:19:35 dlevans Exp $
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
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"

// MVAs...
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/MuonIDMVA.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"

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

        // synchronisation
        void testElectronObjects(const edm::Event& iEvent, const edm::EventSetup &iSetup,
                const TransientTrackBuilder *ttBuilder,
                EcalClusterLazyTools *clusterTools);

        void testMuonObjects(const edm::Event& iEvent, const edm::EventSetup &setup,
                const TransientTrackBuilder *ttBuilder);

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

        // get trigger prescales without object matching
        void getTriggerPrescales(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
                const std::vector<edm::InputTag> &trigNames, std::vector<unsigned int> &prescale);

        // ----------member data ---------------------------

        // lepton tree
        TFile       *leptonFile_;
        LeptonTree  *leptonTree_;

        // input tags
        edm::InputTag electronsInputTag_;
        edm::InputTag muonsInputTag_;
        edm::InputTag photonsInputTag_;
        edm::InputTag primaryVertexInputTag_;
        edm::InputTag rhoIsoAllInputTag;
        edm::InputTag rhoIsoAllCentralInputTag;
        edm::InputTag rhoIsoNeutralInputTag;
        edm::InputTag rhoIsoNeutral2011InputTag;

        edm::InputTag pfCandsInputTag_;
        edm::InputTag pfNoPUCandsInputTag_;
        edm::InputTag metInputTag_;
        edm::InputTag jetsInputTag_;
        edm::InputTag conversionsInputTag_;
        edm::InputTag beamSpotInputTag_;
        edm::InputTag genParticlesInputTag_;

        // jet corrections
        std::string pfJetCorrectorL1FastL2L3_;

        // trigger names
        std::vector<edm::InputTag> electronFRTriggerNames_;
        std::vector<edm::InputTag> muonFRTriggerNames_;

        // trigger related
        const trigger::TriggerEvent   *triggerEvent_;
        HLTConfigProvider          hltConfig_;
        std::string                processName_;
        const edm::TriggerResults* triggerResults_;

        std::vector<edm::InputTag> photonTriggers_;
        std::vector<unsigned int> photonTriggerPrescales_;
        std::vector<unsigned int> photonTriggerVersions_;
        std::vector<edm::InputTag> muTriggers_;    
        std::vector<unsigned int> muTriggerPrescalesTag_;
        std::vector<unsigned int> muTriggerPrescalesProbe_;
        std::vector<unsigned int> muTriggerVersions_;
        std::vector<edm::InputTag> eleTriggers_;       
        std::vector<unsigned int> eleTriggerPrescalesTag_;
        std::vector<unsigned int> eleTriggerPrescalesProbe_;
        std::vector<unsigned int> eleTriggerVersions_;

        // electron id related
        ElectronIDMVA           *reader_electronHWW2011MVA_;
        EGammaMvaEleEstimator   *reader_egammaTrigMVA_;
        //EGammaMvaEleEstimator   *reader_egammaNonTrigMVA_;

        // muon id related
        MuonMVAEstimator        *reader_muonHZZ2012IsoRingsMVA_;

        // common products
        double rhoIsoAll_;
        double rhoIsoAllCentral_;
        double rhoIsoNeutral_;
        double rhoIsoNeutral2011_;

        reco::PFCandidateCollection pfCandCollection_;
        reco::PFCandidateCollection pfNoPUCandCollection_;
        edm::Handle<reco::VertexCollection> vtx_h_;
        reco::VertexCollection vertexCollection_;
        edm::Handle<edm::View<reco::PFJet> > jets_h_;
        reco::Vertex pv_;
        JetCorrector *jetCorrector_;
        reco::GenParticleCollection genParticleCollection_;

        // random number generators
        float rndm_;
        PFPileUpAlgo *pfPileUpAlgo_;

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
    genParticlesInputTag_   =  iConfig.getParameter<edm::InputTag>("genParticlesInputTag");
    primaryVertexInputTag_  =  iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    rhoIsoAllInputTag       =  iConfig.getParameter<edm::InputTag>("rhoIsoAllInputTag");
    rhoIsoAllCentralInputTag =  iConfig.getParameter<edm::InputTag>("rhoIsoAllCentralInputTag");
    rhoIsoNeutralInputTag   =  iConfig.getParameter<edm::InputTag>("rhoIsoNeutralInputTag");
    rhoIsoNeutral2011InputTag   =  iConfig.getParameter<edm::InputTag>("rhoIsoNeutral2011InputTag");

    pfNoPUCandsInputTag_    =  iConfig.getParameter<edm::InputTag>("pfNoPUCandsInputTag");
    pfCandsInputTag_        =  iConfig.getParameter<edm::InputTag>("pfCandsInputTag");

    metInputTag_            =  iConfig.getParameter<edm::InputTag>("metInputTag");
    jetsInputTag_           =  iConfig.getParameter<edm::InputTag>("jetsInputTag");
    conversionsInputTag_    =  iConfig.getParameter<edm::InputTag>("conversionsInputTag");
    beamSpotInputTag_       =  iConfig.getParameter<edm::InputTag>("beamSpotInputTag");

    pfJetCorrectorL1FastL2L3_   = iConfig.getParameter<std::string>("pfJetCorrectorL1FastL2L3");

    photonTriggers_         = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("photonTriggers");
    muTriggers_             = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("muTriggers");
    eleTriggers_            = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("eleTriggers");

    // for pf no pileup emulation
    pfPileUpAlgo_ = new PFPileUpAlgo();

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

    printf("struct TriggerResults {\n");
    for (unsigned int i = 0; i < photonTriggers_.size(); ++i) {
        printf("unsigned int %s_;\n", photonTriggers_[i].process().c_str());
    }
    for (unsigned int i = 0; i < muTriggers_.size(); ++i) {
        printf("unsigned int %s_tag_;\n", muTriggers_[i].process().c_str());
        printf("unsigned int %s_probe_;\n", muTriggers_[i].process().c_str());
    }
    for (unsigned int i = 0; i < eleTriggers_.size(); ++i) {
        printf("unsigned int %s_tag_;\n", eleTriggers_[i].process().c_str());
        printf("unsigned int %s_probe_;\n", eleTriggers_[i].process().c_str());
    }
    printf("};\n");

    photonTriggerPrescales_.reserve(photonTriggers_.size());
    photonTriggerVersions_.reserve(photonTriggers_.size());
    for (unsigned int i = 0; i < photonTriggers_.size(); ++i) {
        const char *trigName = photonTriggers_[i].process().c_str();
        leptonTree_->tree_->Branch(Form("%s", trigName), &photonTriggerPrescales_[i], Form("%s/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_version", trigName), &photonTriggerVersions_[i], Form("%s_version/i", trigName));
        printf("eventTree->SetBranchAddress(\"%s\", &triggerResults.%s);\n", trigName, trigName);
    }

    muTriggerPrescalesTag_.reserve(muTriggers_.size());
    muTriggerPrescalesProbe_.reserve(muTriggers_.size());
    muTriggerVersions_.reserve(muTriggers_.size());
    for (unsigned int i = 0; i < muTriggers_.size(); ++i) {
        const char *trigName = muTriggers_[i].process().c_str();
        leptonTree_->tree_->Branch(Form("%s_tag", trigName), &muTriggerPrescalesTag_[i], Form("%s_tag/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_probe", trigName), &muTriggerPrescalesProbe_[i], Form("%s_probe/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_version", trigName), &muTriggerVersions_[i], Form("%s_version/i", trigName));
        printf("eventTree->SetBranchAddress(\"%s_tag\", &triggerResults.%s_tag_);\n", trigName, trigName);
        printf("eventTree->SetBranchAddress(\"%s_probe\", &triggerResults.%s_probe_);\n", trigName, trigName);
    }

    eleTriggerPrescalesTag_.reserve(eleTriggers_.size());
    eleTriggerPrescalesProbe_.reserve(eleTriggers_.size());
    eleTriggerVersions_.reserve(eleTriggers_.size());
    for (unsigned int i = 0; i < eleTriggers_.size(); ++i) {
        const char *trigName = eleTriggers_[i].process().c_str();
        leptonTree_->tree_->Branch(Form("%s_tag", trigName), &eleTriggerPrescalesTag_[i], Form("%s_tag/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_probe", trigName), &eleTriggerPrescalesProbe_[i], Form("%s_probe/i", trigName));
        leptonTree_->tree_->Branch(Form("%s_version", trigName), &eleTriggerVersions_[i], Form("%s_version/i", trigName));
        printf("eventTree->SetBranchAddress(\"%s_tag\", &triggerResults.%s_tag_);\n", trigName, trigName);
        printf("eventTree->SetBranchAddress(\"%s_probe\", &triggerResults.%s_probe_);\n", trigName, trigName);
    }

    //
    // set up trigger
    //

    processName_ = "";

    //
    // set up electron id
    //

    std::string cmssw_base = getenv("CMSSW_BASE");

    reader_electronHWW2011MVA_ = new ElectronIDMVA();
    reader_electronHWW2011MVA_->Initialize("BDTG method",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml",                
            ElectronIDMVA::kWithIPInfo);

    std::vector<std::string> myManualCatWeigthsTrig;
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat1.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat2.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat3.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat4.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat5.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat6.weights.xml");
    reader_egammaTrigMVA_ = new EGammaMvaEleEstimator();
    reader_egammaTrigMVA_->initialize("BDT",
            EGammaMvaEleEstimator::kTrig,
            true,
            myManualCatWeigthsTrig);

    //std::vector<std::string> myManualCatWeigthsNonTrig;
    //myManualCatWeigthsNonTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml");
    //myManualCatWeigthsNonTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml");
    //myManualCatWeigthsNonTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml");
    //myManualCatWeigthsNonTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml");
    //myManualCatWeigthsNonTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml");
    //myManualCatWeigthsNonTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml");
    //reader_egammaNonTrigMVA_ = new EGammaMvaEleEstimator();
    //reader_egammaNonTrigMVA_->initialize("BDT",
    //        EGammaMvaEleEstimator::kNonTrig,
    //        true,
    //        myManualCatWeigthsNonTrig);

    // si/andrew Iso rings
    std::vector<std::string> muonHZZ2012IsoRingsWeights;
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
    reader_muonHZZ2012IsoRingsMVA_ = new MuonMVAEstimator();
    reader_muonHZZ2012IsoRingsMVA_->initialize("muonHZZ2012IsoRingsMVA", MuonMVAEstimator::kIsoRings, true, muonHZZ2012IsoRingsWeights);

}

LeptonTreeMaker::~LeptonTreeMaker()
{

    //
    // tidy up
    //

    if (reader_electronHWW2011MVA_)         delete reader_electronHWW2011MVA_;
    if (reader_egammaTrigMVA_)           delete reader_egammaTrigMVA_;
    //if (reader_egammaNonTrigMVA_)           delete reader_egammaNonTrigMVA_;
    if (reader_muonHZZ2012IsoRingsMVA_)     delete reader_muonHZZ2012IsoRingsMVA_;
    if (pfPileUpAlgo_)                      delete pfPileUpAlgo_;

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

        // for SNT
        smurfutilities::DumpSaveTags("HLT_DoubleEle14_CaloIdT_TrkIdVL_Mass8_PFMET40_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleEle14_CaloIdT_TrkIdVL_Mass8_PFMET50_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT175_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT225_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleMu14_Mass8_PFMET40_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleMu14_Mass8_PFMET50_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleMu8_Mass8_PFHT175_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleMu8_Mass8_PFHT225_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleRelIso1p0Mu5_Mass8_PFHT175_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_DoubleRelIso1p0Mu5_Mass8_PFHT225_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT175_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT225_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu14_Ele14_CaloIdT_TrkIdVL_Mass8_PFMET40_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_Mu14_Ele14_CaloIdT_TrkIdVL_Mass8_PFMET50_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_RelIso1p0Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT175_v*", hltConfig_);
        smurfutilities::DumpSaveTags("HLT_RelIso1p0Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT225_v*", hltConfig_);

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

    // jets
    iEvent.getByLabel(jetsInputTag_, jets_h_);
    const JetCorrector *jetCorrector = JetCorrector::getJetCorrector(pfJetCorrectorL1FastL2L3_, iSetup);

    // rho for isolation
    // all candidate types
    // used for electrons
    edm::Handle<double> rhoIsoAll_h;
    iEvent.getByLabel(rhoIsoAllInputTag, rhoIsoAll_h);
    rhoIsoAll_ = *(rhoIsoAll_h.product());

    // central candidates - used by egamma
    edm::Handle<double> rhoIsoAllCentral_h;
    iEvent.getByLabel(rhoIsoAllCentralInputTag, rhoIsoAllCentral_h);
    rhoIsoAllCentral_ = *(rhoIsoAllCentral_h.product());

    // neutrals only
    // used for muons
    edm::Handle<double> rhoIsoNeutral_h;
    iEvent.getByLabel(rhoIsoNeutralInputTag, rhoIsoNeutral_h);
    rhoIsoNeutral_ = *(rhoIsoNeutral_h.product());

    // netrals only
    // the definition used by muons in 2011
    // that is not stored in the event
    edm::Handle<double> rhoIsoNeutral2011_h;
    iEvent.getByLabel(rhoIsoNeutral2011InputTag, rhoIsoNeutral2011_h);
    rhoIsoNeutral2011_ = *(rhoIsoNeutral2011_h.product());

    // pf candidates
    edm::Handle<reco::PFCandidateCollection> pfCand_h;
    iEvent.getByLabel(pfCandsInputTag_, pfCand_h);
    pfCandCollection_ = *(pfCand_h.product());

    // pf no pileup candidates
    edm::Handle<reco::PFCandidateCollection> pfNoPUCand_h;
    iEvent.getByLabel(pfNoPUCandsInputTag_, pfNoPUCand_h);
    pfNoPUCandCollection_ = *(pfNoPUCand_h.product());

    // gen particles
    if (!iEvent.isRealData()) {
        edm::Handle<reco::GenParticleCollection> genParticles_h;
        iEvent.getByLabel(genParticlesInputTag_, genParticles_h);
        genParticleCollection_ = *(genParticles_h.product());
    }

    // vertices
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h_);
    vertexCollection_ = *(vtx_h_.product());
    if (vertexCollection_.size() > 0) pv_ = vertexCollection_.at(0);
    else return;

    //
    // decide what to do based on trigger information
    //

    bool debug = false;

    if (debug) {
        //testMuonObjects(iEvent, iSetup, ttBuilder);
        testElectronObjects(iEvent, iSetup, ttBuilder, clusterTools);
    }

    if (!debug) {

        // fake rates
        fillElectronFakeRateTree(iEvent, iSetup, ttBuilder, clusterTools, jetCorrector);
        fillMuonFakeRateTree(iEvent, iSetup, ttBuilder, jetCorrector);

        // efficiency
        fillElectronTagAndProbeTree(iEvent, iSetup, ttBuilder, clusterTools);
        fillMuonTagAndProbeTree(iEvent, iSetup, ttBuilder);

        // photons
        fillPhotonTree(iEvent, iSetup, jetCorrector);

    }

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

//
// debug
//

void LeptonTreeMaker::testElectronObjects(const edm::Event& iEvent, const edm::EventSetup &iSetup,
        const TransientTrackBuilder *ttBuilder,
        EcalClusterLazyTools *clusterTools)
{
    
    // electrons
    edm::Handle<reco::GsfElectronCollection> els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    reco::GsfElectronCollection electronCollection = *(els_h.product());

    // beamspot
    edm::Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
    const reco::BeamSpot &thebs = *(beamspot_h.product());
    
    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    iEvent.getByLabel(conversionsInputTag_, conversions_h);
 
    // look for tag and probe
    unsigned int nEle = els_h->size();
    for (unsigned int itag = 0; itag < nEle; ++itag) {
        
        // look for good tag
        reco::GsfElectronRef ele(els_h, itag);
        if (ele->pt() < 20.0)       continue;
        if (fabs(ele->superCluster()->eta()) > 2.5) continue;
   
        bool passFO = false;
        if (smurfselections::passElectronFO2012(ele, pv_, thebs.position(), conversions_h)) {
            passFO = true;
        }

        // pf isolation
        float pfchiso04_ = 0.0;  float pfemiso04_ = 0.0; float pfnhiso04_ = 0.0; float dbeta = 0.0;
        smurfselections::PFIsolation2012(*ele, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, pfchiso04_, pfemiso04_, pfnhiso04_, dbeta, false, true, true, false);
        float ea04_        = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, ele->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        float egammaPOGNeutral = TMath::Max(float(0.0), pfnhiso04_ + pfemiso04_ - ea04_ * TMath::Max(float(0.0), float(rhoIsoAll_)));
        float egammaPOGIso = (pfchiso04_ + egammaPOGNeutral) / ele->pt();

        // mva
        //std::cout << "==== computing MVA value ===" << std::endl;
        float egammaTrigMVA_         = reader_egammaTrigMVA_->mvaValue(*ele, pv_, *ttBuilder, *clusterTools, false);
        //std::cout << "=== done ===" << std::endl;       

        printf("%u %u %u : ", iEvent.id().run(), iEvent.luminosityBlock(), iEvent.id().event());
        printf("pt, eta = %4.4f, %4.4f :\t", ele->pt(), ele->eta());
        printf("passFO = %u : ", passFO);
        printf("rho = %4.4f : ", rhoIsoAll_);
        printf("Iso(ch, em, nh), EA, RelIso = (%4.4f, %4.4f, %4.4f), %4.4f, %4.4f : ", 
            pfchiso04_, pfemiso04_, pfnhiso04_, ea04_, egammaPOGIso);
        printf("MVA = %4.4f\n", egammaTrigMVA_);
 
    }

}


void LeptonTreeMaker::testMuonObjects(const edm::Event& iEvent, const edm::EventSetup &setup,
    const TransientTrackBuilder *ttBuilder)
{

    // muons
    edm::Handle<edm::View<reco::Muon> > mus_h;
    iEvent.getByLabel(muonsInputTag_, mus_h);
    edm::View<reco::Muon> muonCollection = *(mus_h.product());
    edm::View<reco::Muon>::const_iterator mu;
    const reco::GsfElectronCollection nullEls;
    const reco::MuonCollection nullMus;

    // init tree
    initTriggerBranchValues();  
    leptonTree_->InitVariables();

    // loop and print debug tests
    //reader_muonHZZ2012IsoRingsMVA_->SetPrintMVADebug(true);
    for (mu = muonCollection.begin(); mu != muonCollection.end(); ++mu) {

        if (mu->pt() < 10.0)        continue;
        if (fabs(mu->eta()) > 2.4)  continue;

        float ringsMVA = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*mu, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
        bool id        = smurfselections::passMuonID2012(mu, ringsMVA, pv_);

        printf("%u %u %u : ", iEvent.id().run(), iEvent.luminosityBlock(), iEvent.id().event());
        printf("pt, eta = %4.4f, %4.4f :\t", mu->pt(), mu->eta());
        printf("passID = %u : ", id);
        printf("rho = %4.4f : ", rhoIsoAll_);
        printf("RingsMVA = %4.4f\n", ringsMVA);

        leptonTree_->probe_ = mu->p4();
        leptonTree_->muonHZZ2012IsoRingsMVA_ = ringsMVA;
        leptonTree_->tree_->Fill();
    }

    //reader_muonHZZ2012IsoRingsMVA_->SetPrintMVADebug(false);

}

//
// tag and probe
//

void LeptonTreeMaker::fillElectronTagAndProbeTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
        const TransientTrackBuilder *ttBuilder,
        EcalClusterLazyTools *clusterTools)
{

    // electrons
    edm::Handle<reco::GsfElectronCollection> els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    reco::GsfElectronCollection electronCollection = *(els_h.product());
    if (electronCollection.size() < 2) return;

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
        if (fabs(tag->superCluster()->eta()) > 2.5) continue;

        // 2012 ID selection
        float mvaValue  = reader_egammaTrigMVA_->mvaValue(*tag, pv_, *ttBuilder, *clusterTools, false);
        if (!smurfselections::passElectronID2012(tag, pv_, thebs.position(), conversions_h, mvaValue)) continue;

        // 2012 Iso selection
        float chiso     = 0.0; float emiso = 0.0; float nhiso = 0.0; float dbeta = 0.0;
        float ea04      = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04,
                tag->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        smurfselections::PFIsolation2012(*tag, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4,
                chiso, emiso, nhiso, dbeta, false, true, true, false);
        if (!smurfselections::passElectronIso2012(tag, chiso, emiso, nhiso, ea04, rhoIsoAll_)) continue;

        for (unsigned int iprobe = 0; iprobe < nEle; ++iprobe) {

            if (itag == iprobe)             continue;
            reco::GsfElectronRef probe(els_h, iprobe);

            // look for good probe
            if (probe->pt() < 10.0)       continue;
            if (fabs(probe->superCluster()->eta()) > 2.5) continue;

            // construct tag and probe mass
            LorentzVector p4 = tag->p4() + probe->p4();

            initTriggerBranchValues();
            leptonTree_->InitVariables();

            // this is per event
            fillCommonVariables(iEvent);

            // fill the tree
            leptonTree_->eventSelection_     = LeptonTree::ZeeTagAndProbe;
            leptonTree_->probe_              = probe->p4();
            leptonTree_->qProbe_             = probe->charge();
            leptonTree_->tag_                = tag->p4();
            leptonTree_->qTag_               = tag->charge();
            leptonTree_->tagAndProbeMass_    = p4.M();

            leptonTree_->electronHWW2011MVA_       = reader_electronHWW2011MVA_->MVAValue(&*probe, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);
            leptonTree_->egammaPOG2012MVA_         = reader_egammaTrigMVA_->mvaValue(*probe, pv_, *ttBuilder, *clusterTools, false);
            //leptonTree_->el_egammaNonTrigMVA_      = reader_egammaNonTrigMVA_->mvaValue(*probe, pv_, *ttBuilder, *clusterTools, false);

            // pf isolation
            smurfselections::PFIsolation2012(*probe, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.3, 
                leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, leptonTree_->dbeta03_, false, true, true, false);
            smurfselections::PFIsolation2012(*probe, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
                leptonTree_->pfchiso04_, leptonTree_->pfemiso04_, leptonTree_->pfnhiso04_, leptonTree_->dbeta04_, false, true, true, false);
            leptonTree_->radiso03_    = smurfselections::getElectronRadialIsolation(*probe, pfNoPUCandCollection_, 0.3);
            leptonTree_->radiso04_    = smurfselections::getElectronRadialIsolation(*probe, pfNoPUCandCollection_, 0.4);
            leptonTree_->iso2011_     = smurfselections::electronIsoValuePF(pfCandCollection_, *probe, pv_, 0.4, 1.0, 0.1, 0.07, 0.025, -999., 0);
            leptonTree_->ea04_        = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, 
                probe->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
            leptonTree_->ea03_        = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, 
                probe->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
            // particle flow variations
            float chiso = 0.0; float emiso = 0.0; float nhiso = 0.0; float dbeta = 0.0;
            // remove shared tracks
            smurfselections::PFIsolation2012(*probe, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
                chiso, emiso, nhiso, dbeta, false, true, true, true);
            leptonTree_->el_test_pfchiso04_trkveto_ = chiso;
            // apply dz cut instead of pf no pileup
            smurfselections::PFIsolation2012(*probe, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
                chiso, emiso, nhiso, dbeta, false, true, false, false);
            leptonTree_->el_test_pfchiso04_dzcut_ = chiso;
            // apply vetoes in the barrel
            smurfselections::PFIsolation2012(*probe, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
                chiso, emiso, nhiso, dbeta, true, true, true, false);
            leptonTree_->el_test_pfchiso04_ebveto_ = chiso;
            leptonTree_->el_test_pfemiso04_ebveto_ = emiso;
            
            // cut based ele id
            leptonTree_->vetoId_    = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::VETO, probe, conversions_h, thebs, vtx_h_, 
                leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);
            leptonTree_->looseId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE, probe, conversions_h, thebs, vtx_h_, 
                leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);
            leptonTree_->mediumId_  = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::MEDIUM, probe, conversions_h, thebs, vtx_h_, 
                leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);
            leptonTree_->tightId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::TIGHT, probe, conversions_h, thebs, vtx_h_, 
                leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);

            // higgs ids
            if (smurfselections::passElectronFO2011(probe, pv_, thebs.position(), conversions_h)) 
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
            if (smurfselections::passElectronID2011(probe, pv_, thebs.position(), conversions_h, leptonTree_->electronHWW2011MVA_))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
            if (smurfselections::passElectronIso2011(probe, pfCandCollection_, pv_))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);
            if (smurfselections::passElectronFO2012(probe, pv_, thebs.position(), conversions_h))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFOICHEP2012);
            if (smurfselections::passElectronID2012(probe, pv_, thebs.position(), conversions_h, leptonTree_->egammaPOG2012MVA_))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIDICHEP2012);
            if (smurfselections::passElectronIso2012(probe, 
                leptonTree_->pfchiso04_, leptonTree_->pfemiso04_, 
                leptonTree_->pfnhiso04_, leptonTree_->ea04_, rhoIsoAll_))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIsoICHEP2012);

            // input variables
            leptonTree_->pfmva_     = probe->mva();
            leptonTree_->sceta_     = probe->superCluster()->eta();
            leptonTree_->scenergy_  = probe->superCluster()->energy();
            leptonTree_->chargesAgree_ = smurfselections::threeChargesAgree(*probe);
            leptonTree_->detain_    = probe->deltaEtaSuperClusterTrackAtVtx();
            leptonTree_->dphiin_    = probe->deltaPhiSuperClusterTrackAtVtx();
            leptonTree_->sieie_     = probe->sigmaIetaIeta();
            leptonTree_->fbrem_     = probe->fbrem();
            leptonTree_->eopin_     = probe->eSuperClusterOverP();
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

            // trigger matching
            if (iEvent.isRealData()) {
                // leptons
                objectMatchTrigger(iEvent, iSetup, eleTriggers_, allObjects, tag->p4(), eleTriggerPrescalesTag_);
                objectMatchTrigger(iEvent, iSetup, eleTriggers_, allObjects, probe->p4(), eleTriggerPrescalesProbe_);
                getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
                getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
                getTriggerVersions(iEvent, iSetup, eleTriggers_, eleTriggerVersions_);
            }

            // probe gen matching
            if (!iEvent.isRealData()) {
                leptonTree_->gen_drs1_ = smurfutilities::MatchGenParticle(genParticleCollection_, probe->p4(), 11, 1);
                leptonTree_->gen_drs3_ = smurfutilities::MatchGenParticle(genParticleCollection_, probe->p4(), 11, 3);
            }

            // fill the tree
            leptonTree_->tree_->Fill(); 

        }

    }

}

void LeptonTreeMaker::fillMuonTagAndProbeTree(const edm::Event& iEvent, const edm::EventSetup &iSetup,
        const TransientTrackBuilder *ttBuilder)
{

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
    const reco::GsfElectronCollection nullEls;
    const reco::MuonCollection nullMus;

    for (tag = muonCollection.begin(); tag != muonCollection.end(); ++tag) {

        // look for good tag
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.4) continue;

        // 2012 ID selection
        float mvaValue = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*tag, pv_, pfCandCollection_,
            rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
        if (!smurfselections::passMuonID2012(tag, mvaValue, pv_))   continue;

        // 2012 Iso selection
        if (!smurfselections::passMuonIso2012(tag, mvaValue))  continue;

        for (probe = muonCollection.begin(); probe != muonCollection.end(); ++probe) {

            // look for good probe
            if (tag == probe)            continue;
            if (probe->pt() < 10.0)      continue;
            if (fabs(probe->eta()) > 2.4) continue;

            // construct tag and probe mass
            LorentzVector p4 = tag->p4() + probe->p4();

            initTriggerBranchValues();
            leptonTree_->InitVariables();

            // this is per event
            fillCommonVariables(iEvent);

            // fill the tree
            leptonTree_->eventSelection_     = LeptonTree::ZmmTagAndProbe;
            leptonTree_->probe_              = probe->p4();
            leptonTree_->qProbe_             = probe->charge();
            leptonTree_->tag_                = tag->p4();
            leptonTree_->qTag_               = tag->charge();
            leptonTree_->tagAndProbeMass_    = p4.M();

            #ifdef RELEASE_52X
            leptonTree_->pfemiso03_      = probe->pfIsolationR03().sumPhotonEt;
            leptonTree_->pfchiso03_      = probe->pfIsolationR03().sumChargedHadronPt;
            leptonTree_->pfnhiso03_      = probe->pfIsolationR03().sumNeutralHadronEt;
            leptonTree_->pfemiso04_      = probe->pfIsolationR04().sumPhotonEt;
            leptonTree_->pfchiso04_      = probe->pfIsolationR04().sumChargedHadronPt;
            leptonTree_->pfnhiso04_      = probe->pfIsolationR04().sumNeutralHadronEt;
            leptonTree_->dbeta03_        = probe->pfIsolationR03().sumPUPt;
            leptonTree_->dbeta04_        = probe->pfIsolationR04().sumPUPt;
            #endif
            leptonTree_->radiso03_       = smurfselections::getMuonRadialIsolation(*probe, pfNoPUCandCollection_, 0.3);
            leptonTree_->radiso04_       = smurfselections::getMuonRadialIsolation(*probe, pfNoPUCandCollection_, 0.4);
            leptonTree_->iso2011_        = smurfselections::muonIsoValuePF(pfCandCollection_, *probe, pv_, 0.3, 1.0, 0.1, 0);
            leptonTree_->eaem04_         = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIso04, probe->p4().eta(), MuonEffectiveArea::kMuEAData2012);
            leptonTree_->eanh04_         = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIso04, probe->p4().eta(), MuonEffectiveArea::kMuEAData2012);

            leptonTree_->muonHZZ2012IsoRingsMVA_    = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*probe, pv_, pfCandCollection_, rhoIsoAll_, 
                MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);

            if (smurfselections::passMuonFO2011(probe, pfCandCollection_, pv_))  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
            if (smurfselections::passMuonID2011(probe, pfCandCollection_, pv_))  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
            if (smurfselections::passMuonIso2011(probe, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);
            if (smurfselections::passMuonFO2012(probe, 
                leptonTree_->muonHZZ2012IsoRingsMVA_,pv_))                       leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFOICHEP2012);
            if (smurfselections::passMuonID2012(probe, 
                leptonTree_->muonHZZ2012IsoRingsMVA_, pv_))                      leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIDICHEP2012);
            if (smurfselections::passMuonIso2012(probe, 
                leptonTree_->muonHZZ2012IsoRingsMVA_))                           leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsoICHEP2012);
            if (smurfselections::passMuonIsPOGTight(probe, pv_))                 leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGTight);
            if (smurfselections::passMuonIsPOGSoft(probe, pv_))                  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGSoft);
            if (smurfselections::passMuonHPASameSign(probe, pv_))                leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsHPASS);
            if (probe->isTrackerMuon())                                          leptonTree_->leptonSelection_ |= (LeptonTree::PassMuTrackerMuon);
            if (probe->isGlobalMuon())                                           leptonTree_->leptonSelection_ |= (LeptonTree::PassMuGlobalMuon);
            #ifdef RELEASE_52X
            if (probe->isPFMuon())                                               leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPF);
            #endif

            // probe trigger matching
            if (iEvent.isRealData()) {
                objectMatchTrigger(iEvent, iSetup, muTriggers_, allObjects, probe->p4(), muTriggerPrescalesProbe_);
                objectMatchTrigger(iEvent, iSetup, muTriggers_, allObjects, tag->p4(), muTriggerPrescalesTag_);
                getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
                getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
                getTriggerVersions(iEvent, iSetup, muTriggers_, muTriggerVersions_);
            }

            // probe gen matching
            if (!iEvent.isRealData()) {
                leptonTree_->gen_drs1_ = smurfutilities::MatchGenParticle(genParticleCollection_, probe->p4(), 13, 1);
                leptonTree_->gen_drs3_ = smurfutilities::MatchGenParticle(genParticleCollection_, probe->p4(), 13, 3);
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
    initTriggerBranchValues();

    // electrons
    edm::Handle<reco::GsfElectronCollection> els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    reco::GsfElectronCollection electronCollection = *(els_h.product());
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

    // look for a good FO
    unsigned int nEle = els_h->size();
    unsigned int nFO = 0;
    reco::GsfElectronRef fo;
    for (unsigned int iele = 0; iele < nEle; ++iele) {
        reco::GsfElectronRef ele(els_h, iele);
        if (ele->pt()        < 10.0) continue;
        if (fabs(ele->superCluster()->eta()) > 2.5)  continue;
        if (!smurfselections::passElectronFO2012(ele, pv_, thebs.position(), conversions_h)) continue;
        ++nFO;
        fo = ele;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {

        fillCommonVariables(iEvent);
        fillJets(iEvent, iSetup, *fo, jetCorrector);

        // probe trigger matching
        if (iEvent.isRealData()) {
            // leptons
            getTriggerPrescales(iEvent, iSetup, eleTriggers_, eleTriggerPrescalesProbe_);
            getTriggerVersions(iEvent, iSetup, eleTriggers_, eleTriggerVersions_);
            // photons
            getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
            getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
        }

        leptonTree_->eventSelection_     |= LeptonTree::QCDFakeEle;
        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();
        leptonTree_->mt_                 = smurfutilities::Mt(met, fo->p4().Pt(), reco::deltaPhi(met_h->front().phi(), fo->p4().Phi()));

        leptonTree_->electronHWW2011MVA_       = reader_electronHWW2011MVA_->MVAValue(&*fo, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);;
        leptonTree_->egammaPOG2012MVA_         = reader_egammaTrigMVA_->mvaValue(*fo, pv_, *ttBuilder, *clusterTools);
        //leptonTree_->el_egammaNonTrigMVA_      = reader_egammaNonTrigMVA_->mvaValue(*fo, pv_, *ttBuilder, *clusterTools);

        // pf isolation
        smurfselections::PFIsolation2012(*fo, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.3, 
            leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, leptonTree_->dbeta03_, false, true, true, false);
        smurfselections::PFIsolation2012(*fo, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
            leptonTree_->pfchiso04_, leptonTree_->pfemiso04_, leptonTree_->pfnhiso04_, leptonTree_->dbeta04_, false, true, true, false);
        leptonTree_->radiso03_    = smurfselections::getElectronRadialIsolation(*fo, pfNoPUCandCollection_, 0.3);
        leptonTree_->radiso04_    = smurfselections::getElectronRadialIsolation(*fo, pfNoPUCandCollection_, 0.4);
        leptonTree_->iso2011_     = smurfselections::electronIsoValuePF(pfCandCollection_, *fo, pv_, 0.4, 1.0, 0.1, 0.07, 0.025, -999., 0);
        leptonTree_->ea04_        = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, 
            fo->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        leptonTree_->ea03_        = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, 
            fo->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);

        // particle flow variations
        float chiso = 0.0; float emiso = 0.0; float nhiso = 0.0; float dbeta = 0.0;
        // remove shared tracks
        smurfselections::PFIsolation2012(*fo, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
            chiso, emiso, nhiso, dbeta, false, true, true, true);
        leptonTree_->el_test_pfchiso04_trkveto_ = chiso;
        // apply dz cut instead of pf no pileup
        smurfselections::PFIsolation2012(*fo, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
            chiso, emiso, nhiso, dbeta, false, true, false, false);
        leptonTree_->el_test_pfchiso04_dzcut_ = chiso;
        // apply vetoes in the barrel
        smurfselections::PFIsolation2012(*fo, pfCandCollection_, vertexCollection_, pfPileUpAlgo_, 0, 0.4, 
            chiso, emiso, nhiso, dbeta, true, true, true, false);
        leptonTree_->el_test_pfchiso04_ebveto_ = chiso;
        leptonTree_->el_test_pfemiso04_ebveto_ = emiso;

        // cut based ele id
        leptonTree_->vetoId_    = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::VETO, fo, conversions_h, thebs, vtx_h_, 
            leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);
        leptonTree_->looseId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE, fo, conversions_h, thebs, vtx_h_, 
            leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);
        leptonTree_->mediumId_  = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::MEDIUM, fo, conversions_h, thebs, vtx_h_, 
            leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);
        leptonTree_->tightId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::TIGHT, fo, conversions_h, thebs, vtx_h_, 
            leptonTree_->pfchiso03_, leptonTree_->pfemiso03_, leptonTree_->pfnhiso03_, rhoIsoAllCentral_);

        // higgs ids
        if (smurfselections::passElectronFO2011(fo, pv_, thebs.position(), conversions_h))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
        if (smurfselections::passElectronID2011(fo, pv_, thebs.position(), conversions_h, leptonTree_->electronHWW2011MVA_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
        if (smurfselections::passElectronIso2011(fo, pfCandCollection_, pv_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);
        if (smurfselections::passElectronFO2012(fo, pv_, thebs.position(), conversions_h))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFOICHEP2012);
        if (smurfselections::passElectronID2012(fo, pv_, thebs.position(), conversions_h, leptonTree_->egammaPOG2012MVA_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIDICHEP2012);
        if (smurfselections::passElectronIso2012(fo,
            leptonTree_->pfchiso04_, leptonTree_->pfemiso04_,
            leptonTree_->pfnhiso04_, leptonTree_->ea04_, rhoIsoAll_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIsoICHEP2012);

        leptonTree_->pfmva_     = fo->mva();
        leptonTree_->sceta_     = fo->superCluster()->eta();
        leptonTree_->scenergy_  = fo->superCluster()->energy();
        leptonTree_->chargesAgree_ = smurfselections::threeChargesAgree(*fo);
        leptonTree_->detain_    = fo->deltaEtaSuperClusterTrackAtVtx();
        leptonTree_->dphiin_    = fo->deltaPhiSuperClusterTrackAtVtx();
        leptonTree_->sieie_     = fo->sigmaIetaIeta();
        leptonTree_->fbrem_     = fo->fbrem();
        leptonTree_->eopin_     = fo->eSuperClusterOverP();
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

        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillMuonFakeRateTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const TransientTrackBuilder *ttBuilder, const JetCorrector *jetCorrector)
{

    leptonTree_->InitVariables();
    initTriggerBranchValues();

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
    const reco::GsfElectronCollection nullEls;
    const reco::MuonCollection nullMus;
    unsigned int nFO = 0;
    for (it = muonCollection.begin(); it != muonCollection.end(); ++it) {
        if (it->pt()        < 10.0) continue;
        if (fabs(it->eta()) > 2.4)  continue;
        float mvaValue = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*it, pv_, pfCandCollection_, rhoIsoAll_, 
            MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
        if (!smurfselections::passMuonFO2012(it, mvaValue, pv_)) continue;
        ++nFO;
        fo = it;
    }

    // if exactly one good FO
    // then fill the photon tree
    if (nFO == 1) {

        fillCommonVariables(iEvent);    
        fillJets(iEvent, iSetup, *fo, jetCorrector);

        // probe trigger matching
        if (iEvent.isRealData()) {
            // leptons
            getTriggerPrescales(iEvent, iSetup, muTriggers_, muTriggerPrescalesProbe_);
            getTriggerVersions(iEvent, iSetup, muTriggers_, muTriggerVersions_);
            // photons
            getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
            getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
        }

        leptonTree_->eventSelection_    |= LeptonTree::QCDFakeMu;
        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();
        leptonTree_->mt_                 = smurfutilities::Mt(met, fo->p4().Pt(), reco::deltaPhi(met_h->front().phi(), fo->p4().Phi()));

        #ifdef RELEASE_52X
        leptonTree_->pfemiso03_      = fo->pfIsolationR03().sumPhotonEt;
        leptonTree_->pfchiso03_      = fo->pfIsolationR03().sumChargedHadronPt;
        leptonTree_->pfnhiso03_      = fo->pfIsolationR03().sumNeutralHadronEt;
        leptonTree_->pfemiso04_      = fo->pfIsolationR04().sumPhotonEt;
        leptonTree_->pfchiso04_      = fo->pfIsolationR04().sumChargedHadronPt;
        leptonTree_->pfnhiso04_      = fo->pfIsolationR04().sumNeutralHadronEt;
        leptonTree_->dbeta03_        = fo->pfIsolationR03().sumPUPt;
        leptonTree_->dbeta04_        = fo->pfIsolationR04().sumPUPt;
        #endif
        leptonTree_->radiso03_       = smurfselections::getMuonRadialIsolation(*fo, pfNoPUCandCollection_, 0.3);
        leptonTree_->radiso04_       = smurfselections::getMuonRadialIsolation(*fo, pfNoPUCandCollection_, 0.4);
        leptonTree_->iso2011_        = smurfselections::muonIsoValuePF(pfCandCollection_, *fo, pv_, 0.3, 1.0, 0.1, 0);
        leptonTree_->eaem04_         = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIso04, fo->p4().eta(), MuonEffectiveArea::kMuEAData2012);
        leptonTree_->eanh04_         = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIso04, fo->p4().eta(), MuonEffectiveArea::kMuEAData2012);

        leptonTree_->muonHZZ2012IsoRingsMVA_    = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*fo, pv_, pfCandCollection_, rhoIsoAll_, 
            MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);

        leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
        if (smurfselections::passMuonID2011(fo, pfCandCollection_, pv_))  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
        if (smurfselections::passMuonIso2011(fo, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);
        if (smurfselections::passMuonFO2012(fo, 
            leptonTree_->muonHZZ2012IsoRingsMVA_, pv_))                   leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFOICHEP2012);
        if (smurfselections::passMuonID2012(fo, 
            leptonTree_->muonHZZ2012IsoRingsMVA_, pv_))                   leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIDICHEP2012);
        if (smurfselections::passMuonIso2012(fo, 
            leptonTree_->muonHZZ2012IsoRingsMVA_))                        leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsoICHEP2012);
        if (smurfselections::passMuonIsPOGTight(fo, pv_))                 leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGTight);
        if (smurfselections::passMuonIsPOGSoft(fo, pv_))                  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGSoft);
        if (smurfselections::passMuonHPASameSign(fo, pv_))                leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsHPASS);
        if (fo->isTrackerMuon())                                          leptonTree_->leptonSelection_ |= (LeptonTree::PassMuTrackerMuon);
        if (fo->isGlobalMuon())                                           leptonTree_->leptonSelection_ |= (LeptonTree::PassMuGlobalMuon);
        #ifdef RELEASE_52X
        if (fo->isPFMuon())                                               leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPF);
        #endif

        leptonTree_->tree_->Fill();
    }

}

void LeptonTreeMaker::fillPhotonTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
        const JetCorrector *jetCorrector)
{

    leptonTree_->InitVariables();
    initTriggerBranchValues();

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

        if (!smurfselections::passPhotonSelection2011(it, rhoIsoAll_)) continue;
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
        leptonTree_->mt_                 = smurfutilities::Mt(met, photon->p4().Pt(), reco::deltaPhi(metPhi, photon->p4().Phi()));

        const reco::Candidate* phocand = &(*photon);
        std::vector<const reco::Candidate*> phos;
        phos.push_back(phocand);
        std::pair<double,double> tkmet = smurfselections::trackerMET(phos, 0.1, pfCandCollection_, pv_);
        leptonTree_->trackMet_           = tkmet.first;
        leptonTree_->trackMetPhi_        = tkmet.second;

        // probe trigger matching
        if (iEvent.isRealData()) {
            getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
            getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
            unsigned int lowestPrescale = 0;
            for (unsigned int i = 0; i < photonTriggerPrescales_.size(); ++i) {
                if (photonTriggerPrescales_[i] != 0) {
                    if (photonTriggerPrescales_[i] < lowestPrescale || lowestPrescale == 0) 
                        lowestPrescale = photonTriggerPrescales_[i];
                }
            }
            leptonTree_->hltPrescale_ = lowestPrescale;
        }

        leptonTree_->tree_->Fill();

    }

}

void LeptonTreeMaker::initTriggerBranchValues()
{
    // init the dynamic trigger branches

    photonTriggerPrescales_.clear();
    photonTriggerVersions_.clear();
    for (unsigned int i = 0; i < photonTriggers_.size(); ++i) {
        photonTriggerPrescales_.push_back(0);
        photonTriggerVersions_.push_back(0);
    }

    muTriggerPrescalesTag_.clear();
    muTriggerPrescalesProbe_.clear();
    for (unsigned int i = 0; i < muTriggers_.size(); ++i) {
        muTriggerPrescalesTag_.push_back(0);
        muTriggerPrescalesProbe_.push_back(0);
        muTriggerVersions_.push_back(0);
    }

    eleTriggerPrescalesTag_.clear();
    eleTriggerPrescalesProbe_.clear();
    eleTriggerVersions_.clear();
    for (unsigned int i = 0; i < eleTriggers_.size(); ++i) {
        eleTriggerPrescalesTag_.push_back(0);
        eleTriggerPrescalesProbe_.push_back(0);
        eleTriggerVersions_.push_back(0);
    }
}

void LeptonTreeMaker::fillCommonVariables(const edm::Event& iEvent)
{

    // fill
    leptonTree_->run_       = iEvent.id().run()        ;
    leptonTree_->event_     = iEvent.id().event()      ;
    leptonTree_->lumi_      = iEvent.luminosityBlock() ;
    leptonTree_->rnd_       = rndm_                    ;
    leptonTree_->rhoIsoAll_     = rhoIsoAll_;
    leptonTree_->rhoIsoAllCentral_ = rhoIsoAllCentral_;
    leptonTree_->rhoIsoNeutral_ = rhoIsoNeutral_;
    leptonTree_->nvtx_      = smurfselections::CountGoodPV(vtx_h_);
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
    std::vector<std::pair<reco::PFJet, float> > jets5 = smurfselections::goodJets(iEvent, iSetup, jets_h_, cand1, jetCorrector, 5.);
    if (jets5.size() > 0) {
        leptonTree_->jet1_          = jets5.at(0).first.p4() * jets5.at(0).second;
        leptonTree_->dPhiProbeJet1_ = reco::deltaPhi(cand1.p4().Phi(), jets5.at(0).first.p4().Phi());
    }
    if (jets5.size() > 1) leptonTree_->jet2_ = jets5.at(1).first.p4() * jets5.at(1).second;
    if (jets5.size() > 2) leptonTree_->jet3_ = jets5.at(2).first.p4() * jets5.at(2).second;


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
                    unsigned int version = ((TObjString*)substrArr->At(1))->GetString().Atoi();
                    versions[t] = version;
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

void LeptonTreeMaker::getTriggerPrescales(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
        const std::vector<edm::InputTag> &trigNames, std::vector<unsigned int> &prescale)
{

    for (unsigned int t = 0; t < trigNames.size(); ++t) {

        for (unsigned int i = 0; i < hltConfig_.size(); i++) {

            // did trigger pass
            bool result = triggerResults_->accept(i);
            if (!result) continue;

            // does trigger have the right name
            TString hltTrigName(hltConfig_.triggerName(i));
            TString pattern(trigNames[t].label());
            hltTrigName.ToLower();
            pattern.ToLower();
            TRegexp reg(Form("%s", pattern.Data()), true);

            // if so, get prescale value
            if (hltTrigName.Index(reg) >= 0) 
                prescale[t] = hltConfig_.prescaleValue(iEvent, iSetup, hltConfig_.triggerName(i));

        }

    }

}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonTreeMaker);
