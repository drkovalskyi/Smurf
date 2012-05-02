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
// $Id: LeptonTreeMaker.cc,v 1.33 2012/05/02 10:41:56 dlevans Exp $
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
#include "DataFormats/Math/interface/deltaPhi.h"

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
        edm::InputTag rhoJECInputTag;
        edm::InputTag rhoIsoAllInputTag;
        edm::InputTag rhoIsoNeutralInputTag;
        edm::InputTag rhoIsoNeutral2011InputTag;

        edm::InputTag pfCandsInputTag_;
        edm::InputTag pfNoPUCandsInputTag_;
        edm::InputTag metInputTag_;
        edm::InputTag jetsInputTag_;
        edm::InputTag conversionsInputTag_;
        edm::InputTag beamSpotInputTag_;

        // pf isolation related
        std::vector<edm::InputTag> eleIsoVal03InputTags_;
        std::vector<edm::InputTag> eleIsoVal04InputTags_;
        IsoDepositVals *eleIsoVals03_;
        IsoDepositVals *eleIsoVals04_; 

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
        ElectronIDMVA           *reader_electronHWW2011IDIsoMVA_;
        EGammaMvaEleEstimator   *reader_egammaPOG2012MVA_;

        MuonIDMVA               *reader_muonHWW2011IDIsoMVA_;
        MuonMVAEstimator        *reader_muonHZZ2012IDMVA_;
        MuonMVAEstimator        *reader_muonHZZ2012IDIsoRingsMVA_;
        MuonMVAEstimator        *reader_muonHZZ2012IsoRingsMVA_;
        MuonMVAEstimator        *reader_muonHZZ2012IsoDRMVA_;

        // common products
        double rhoJEC_;
        double rhoIsoAll_;
        double rhoIsoNeutral_;
        double rhoIsoNeutral2011_;

        reco::PFCandidateCollection pfCandCollection_;
        reco::PFCandidateCollection pfNoPUCandCollection_;
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
    rhoIsoAllInputTag       =  iConfig.getParameter<edm::InputTag>("rhoIsoAllInputTag");
    rhoIsoNeutralInputTag   =  iConfig.getParameter<edm::InputTag>("rhoIsoNeutralInputTag");
    rhoIsoNeutral2011InputTag   =  iConfig.getParameter<edm::InputTag>("rhoIsoNeutral2011InputTag");

    pfNoPUCandsInputTag_        =  iConfig.getParameter<edm::InputTag>("pfNoPUCandsInputTag");
    pfCandsInputTag_        =  iConfig.getParameter<edm::InputTag>("pfCandsInputTag");

    metInputTag_            =  iConfig.getParameter<edm::InputTag>("metInputTag");
    jetsInputTag_           =  iConfig.getParameter<edm::InputTag>("jetsInputTag");
    conversionsInputTag_    =  iConfig.getParameter<edm::InputTag>("conversionsInputTag");
    beamSpotInputTag_       =  iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
    eleIsoVal03InputTags_   =  iConfig.getParameter<std::vector<edm::InputTag> >("eleIsoVal03InputTags");
    eleIsoVal04InputTags_   =  iConfig.getParameter<std::vector<edm::InputTag> >("eleIsoVal04InputTags");

    pfJetCorrectorL1FastL2L3_   = iConfig.getParameter<std::string>("pfJetCorrectorL1FastL2L3");

    photonTriggers_         = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("photonTriggers");
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

    reader_electronHWW2011IDIsoMVA_ = new ElectronIDMVA();
    reader_electronHWW2011IDIsoMVA_->Initialize("BDTG method",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml",
            ElectronIDMVA::kIDIsoCombined);

    std::vector<std::string> myManualCatWeigthsTrig;
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat1.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat2.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat3.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat4.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat5.weights.xml");
    myManualCatWeigthsTrig.push_back(cmssw_base+"/src/Smurf/ProcessingAndSkimming/data/Electrons_BDTG_TrigV0_Cat6.weights.xml");
    reader_egammaPOG2012MVA_ = new EGammaMvaEleEstimator();
    reader_egammaPOG2012MVA_->initialize("BDT",
            EGammaMvaEleEstimator::kTrig,
            true,
            myManualCatWeigthsTrig);

    reader_muonHWW2011IDIsoMVA_ = new MuonIDMVA();
    reader_muonHWW2011IDIsoMVA_->Initialize("BDTG method",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml",
            cmssw_base+"/src/HiggsAnalysis/HiggsToWW2Leptons/data/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml",
            MuonIDMVA::kIDIsoCombinedDetIso);

    // si/andrew ID
    std::vector<std::string> muonHZZ2012IDWeights;
    //muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    //muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    //muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    //muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    //muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml");
    //muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml");
    muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml");
    muonHZZ2012IDWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml");

    reader_muonHZZ2012IDMVA_ = new MuonMVAEstimator();
    reader_muonHZZ2012IDMVA_->initialize("muonHZZ2012IDMVA", MuonMVAEstimator::kID, true, muonHZZ2012IDWeights);

    // si/andrew Iso rings
    std::vector<std::string> muonHZZ2012IsoRingsWeights;
    //muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    //muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    //muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    //muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    //muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
    //muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
    muonHZZ2012IsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
    reader_muonHZZ2012IsoRingsMVA_ = new MuonMVAEstimator();
    reader_muonHZZ2012IsoRingsMVA_->initialize("muonHZZ2012IsoRingsMVA", MuonMVAEstimator::kIsoRings, true, muonHZZ2012IsoRingsWeights);

    // si/andrew ID+Iso rings combined
    std::vector<std::string> muonHZZ2012IDIsoRingsWeights;
    //muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_barrel_lowpt.weights.xml");
    //muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_barrel_highpt.weights.xml");
    //muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_endcap_lowpt.weights.xml");
    //muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_endcap_highpt.weights.xml");
    //muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_tracker.weights.xml");
    muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_barrel_lowpt.weights.xml");
    muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_endcap_lowpt.weights.xml");
    muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_barrel_highpt.weights.xml");
    muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_endcap_highpt.weights.xml");
    muonHZZ2012IDIsoRingsWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_tracker.weights.xml");

    reader_muonHZZ2012IDIsoRingsMVA_ = new MuonMVAEstimator();
    reader_muonHZZ2012IDIsoRingsMVA_->initialize("muonHZZ2012IDIsoRingsMVA", MuonMVAEstimator::kIDIsoRingsCombined, true, muonHZZ2012IDIsoRingsWeights);

    // Santiago Iso dR
    std::vector<std::string> muonHZZ2012IsoDRWeights;
    muonHZZ2012IsoDRWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_santi-V1_LB_BDT.weights.xml");
    muonHZZ2012IsoDRWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_santi-V1_LE_BDT.weights.xml");
    muonHZZ2012IsoDRWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_santi-V1_HB_BDT.weights.xml");
    muonHZZ2012IsoDRWeights.push_back(cmssw_base+"/src/Muon/MuonAnalysisTools/data/MuonIsoMVA_santi-V1_HE_BDT.weights.xml");
    reader_muonHZZ2012IsoDRMVA_ = new MuonMVAEstimator();
    reader_muonHZZ2012IsoDRMVA_->initialize("muonHZZ2012IsoDRMVA", MuonMVAEstimator::kIsoDeltaR, true, muonHZZ2012IsoDRWeights);

}

LeptonTreeMaker::~LeptonTreeMaker()
{

    //
    // tidy up
    //

    if (reader_electronHWW2011MVA_)         delete reader_electronHWW2011MVA_;
    if (reader_electronHWW2011IDIsoMVA_)    delete reader_electronHWW2011IDIsoMVA_;
    if (reader_egammaPOG2012MVA_)           delete reader_egammaPOG2012MVA_;
    if (reader_muonHWW2011IDIsoMVA_)        delete reader_muonHWW2011IDIsoMVA_;
    if (reader_muonHZZ2012IDMVA_)           delete reader_muonHZZ2012IDMVA_;
    if (reader_muonHZZ2012IDIsoRingsMVA_)   delete reader_muonHZZ2012IDIsoRingsMVA_;
    if (reader_muonHZZ2012IsoRingsMVA_)     delete reader_muonHZZ2012IsoRingsMVA_;
    if (reader_muonHZZ2012IsoDRMVA_)        delete reader_muonHZZ2012IsoDRMVA_;

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

    // rho for jec
    edm::Handle<double> rhoJEC_h;
    iEvent.getByLabel(rhoJECInputTag, rhoJEC_h);
    rhoJEC_ = *(rhoJEC_h.product());

    // jets
    iEvent.getByLabel(jetsInputTag_, jets_h_);
    const JetCorrector *jetCorrector = JetCorrector::getJetCorrector(pfJetCorrectorL1FastL2L3_, iSetup);

    // rho for isolation
    // all candidate types
    // used for electrons
    edm::Handle<double> rhoIsoAll_h;
    iEvent.getByLabel(rhoIsoAllInputTag, rhoIsoAll_h);
    rhoIsoAll_ = *(rhoIsoAll_h.product());

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

    // isolation for electrons
    eleIsoVals03_ = new IsoDepositVals(eleIsoVal03InputTags_.size());
    for (size_t j = 0; j < eleIsoVal03InputTags_.size(); ++j) {
        iEvent.getByLabel(eleIsoVal03InputTags_[j], (*eleIsoVals03_)[j]);
    }
    eleIsoVals04_ = new IsoDepositVals(eleIsoVal04InputTags_.size());
    for (size_t j = 0; j < eleIsoVal04InputTags_.size(); ++j) {
        iEvent.getByLabel(eleIsoVal04InputTags_[j], (*eleIsoVals04_)[j]);
    }

    // vertices
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h_);
    const reco::VertexCollection vertexCollection = *(vtx_h_.product());
    if (vertexCollection.size() > 0) pv_ = vertexCollection.at(0);
    else return;

    //
    // decide what to do based on trigger information
    //

    // fake rates
    fillElectronFakeRateTree(iEvent, iSetup, ttBuilder, clusterTools, jetCorrector);
    fillMuonFakeRateTree(iEvent, iSetup, ttBuilder, jetCorrector);

    // efficiency
    fillElectronTagAndProbeTree(iEvent, iSetup, ttBuilder, clusterTools);
    fillMuonTagAndProbeTree(iEvent, iSetup, ttBuilder);

    // photons
    fillPhotonTree(iEvent, iSetup, jetCorrector);

    //
    // tidy up
    //

    if (clusterTools) delete clusterTools;
    if (eleIsoVals03_) delete eleIsoVals03_;
    if (eleIsoVals04_) delete eleIsoVals04_;

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

    initTriggerBranchValues();  
    leptonTree_->InitVariables(); 

    // this is per event
    fillCommonVariables(iEvent);

    // triggers
    const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();
    getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
    getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);

    // look for tag and probe
    unsigned int nEle = els_h->size();
    for (unsigned int itag = 0; itag < nEle; ++itag) {

        // look for good tag
        reco::GsfElectronRef tag(els_h, itag);
        if (tag->pt() < 20.0)       continue;
        if (fabs(tag->eta()) > 2.5) continue;
        float mvaValue = reader_electronHWW2011MVA_->MVAValue(&*tag, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);
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
            leptonTree_->eventSelection_     = LeptonTree::ZeeTagAndProbe;
            leptonTree_->probe_              = probe->p4();
            leptonTree_->qProbe_             = probe->charge();
            leptonTree_->tag_                = tag->p4();
            leptonTree_->qTag_               = tag->charge();
            leptonTree_->tagAndProbeMass_    = p4.M();

            mvaValue                  = reader_electronHWW2011MVA_->MVAValue(&*probe, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);
            leptonTree_->electronHWW2011MVA_       = mvaValue;
            leptonTree_->electronHWW2011IDIsoMVA_  = reader_electronHWW2011IDIsoMVA_->MVAValue(&*probe, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);

            //std::cout << "event: " << leptonTree_->event_ << " ";
            leptonTree_->egammaPOG2012MVA_         = reader_egammaPOG2012MVA_->mvaValue(*probe, pv_, *ttBuilder, *clusterTools, false);

            if (smurfselections::passElectronFO2011(probe, pv_, thebs.position(), conversions_h))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
            if (smurfselections::passElectronID2011(probe, pv_, thebs.position(), conversions_h, mvaValue))
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
            if (smurfselections::passElectronIso2011(probe, pfCandCollection_, pv_)) 
                leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);

            // pf isolation
            leptonTree_->el_pfchiso03_   = (*(*eleIsoVals03_)[0])[probe];
            leptonTree_->el_pfemiso03_   = (*(*eleIsoVals03_)[1])[probe];
            leptonTree_->el_pfnhiso03_   = (*(*eleIsoVals03_)[2])[probe];
            leptonTree_->el_pfchiso04_   = (*(*eleIsoVals04_)[0])[probe];
            leptonTree_->el_pfemiso04_   = (*(*eleIsoVals04_)[1])[probe];
            leptonTree_->el_pfnhiso04_   = (*(*eleIsoVals04_)[2])[probe];
            leptonTree_->el_radiso03_    = smurfselections::getElectronRadialIsolation(*probe, pfNoPUCandCollection_, 0.3);
            leptonTree_->el_radiso04_    = smurfselections::getElectronRadialIsolation(*probe, pfNoPUCandCollection_, 0.4);
            leptonTree_->el_iso_         = smurfselections::electronIsoValuePF(pfCandCollection_, *probe, pv_, 0.4, 1.0, 0.1, 0.07, 0.025, -999., 0);
            leptonTree_->el_ea04_        = smurfselections::getEGammaEffectiveArea(probe->superCluster()->eta(), smurfselections::EGAMMA2012_04);

            // cut based ele id
            leptonTree_->vetoId_    = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::VETO, probe, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);
            leptonTree_->looseId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE, probe, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);
            leptonTree_->mediumId_  = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::MEDIUM, probe, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);
            leptonTree_->tightId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::TIGHT, probe, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);

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

            // probe trigger matching
            if (iEvent.isRealData()) {
                // leptons
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

    // muons
    edm::Handle<edm::View<reco::Muon> > mus_h;
    iEvent.getByLabel(muonsInputTag_, mus_h);
    edm::View<reco::Muon> muonCollection = *(mus_h.product());
    if (muonCollection.size() < 2) return;

    initTriggerBranchValues();
    leptonTree_->InitVariables();    

    // this is per event
    fillCommonVariables(iEvent);

    // triggers
    const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();
    getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
    getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);

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
            leptonTree_->eventSelection_     = LeptonTree::ZmmTagAndProbe;
            leptonTree_->probe_              = probe->p4();
            leptonTree_->qProbe_             = probe->charge();
            leptonTree_->tag_                = tag->p4();
            leptonTree_->qTag_               = tag->charge();
            leptonTree_->tagAndProbeMass_    = p4.M();

            leptonTree_->mu_pfemiso04_      = probe->pfIsolationR04().sumPhotonEt;
            leptonTree_->mu_pfchiso04_      = probe->pfIsolationR04().sumChargedHadronPt;
            leptonTree_->mu_pfnhiso04_      = probe->pfIsolationR04().sumNeutralHadronEt;
            leptonTree_->mu_radiso03_       = smurfselections::getMuonRadialIsolation(*probe, pfNoPUCandCollection_, 0.3);
            leptonTree_->mu_radiso04_       = smurfselections::getMuonRadialIsolation(*probe, pfNoPUCandCollection_, 0.4);
            leptonTree_->mu_iso_            = smurfselections::muonIsoValuePF(pfCandCollection_, *probe, pv_, 0.3, 1.0, 0.1, 0);
            leptonTree_->mu_eaem04_         = smurfselections::getMuonEffectiveArea(probe->p4().eta(), smurfselections::MUON2012_EM04);
            leptonTree_->mu_eanh04_         = smurfselections::getMuonEffectiveArea(probe->p4().eta(), smurfselections::MUON2012_NH04);

            const reco::GsfElectronCollection nullEls;
            const reco::MuonCollection nullMus;
            leptonTree_->muonHWW2011IDIsoMVA_       = reader_muonHWW2011IDIsoMVA_->MVAValue(&*probe, pv_, ttBuilder, rhoIsoAll_, false);
            leptonTree_->muonHZZ2012IDMVA_          = reader_muonHZZ2012IDMVA_->mvaValue(*probe, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
            leptonTree_->muonHZZ2012IDIsoRingsMVA_  = reader_muonHZZ2012IDIsoRingsMVA_->mvaValue(*probe, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);

//reader_muonHZZ2012IsoRingsMVA_->SetPrintMVADebug(true);
//            std::cout << leptonTree_->run_ << " " << leptonTree_->lumi_ << " " << leptonTree_->event_ << std::endl;
            leptonTree_->muonHZZ2012IsoRingsMVA_    = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*probe, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
//reader_muonHZZ2012IsoRingsMVA_->SetPrintMVADebug(false);

            leptonTree_->muonHZZ2012IsoDRMVA_       = reader_muonHZZ2012IsoDRMVA_->mvaValue(*probe, pfNoPUCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC);

            if (smurfselections::passMuonFO2011(probe, pv_))                     leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
            if (smurfselections::passMuonID2011(probe, pv_))                     leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
            if (smurfselections::passMuonIso2011(probe, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);
            if (smurfselections::passMuonIsPOGTight(probe, pv_))                 leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGTight);
            if (smurfselections::passMuonIsPOGSoft(probe, pv_))                  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGSoft);
            #ifdef RELEASE_52X
            if (probe->isPFMuon())                                               leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPF);
            #endif

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

        // probe trigger matching
        if (iEvent.isRealData()) {
            // leptons
            const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();
            objectMatchTrigger(iEvent, iSetup, eleTriggers_, allObjects, fo->p4(), eleTriggerPrescalesProbe_);
            getTriggerVersions(iEvent, iSetup, eleTriggers_, eleTriggerVersions_);
            // photons
            getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
            getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
        }

        leptonTree_->eventSelection_     |= LeptonTree::QCDFakeEle;
        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();
        leptonTree_->mt_                 = smurfutilities::Mt(met, fo->p4().Pt(), reco::deltaPhi(met_h->front().phi(), fo->p4().Phi()));

        float mvaValue            = reader_electronHWW2011MVA_->MVAValue(&*fo, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);
        leptonTree_->electronHWW2011MVA_       = mvaValue;
        leptonTree_->electronHWW2011IDIsoMVA_  = reader_electronHWW2011IDIsoMVA_->MVAValue(&*fo, pv_, *clusterTools, ttBuilder, pfCandCollection_, rhoIsoAll_);
        leptonTree_->egammaPOG2012MVA_         = reader_egammaPOG2012MVA_->mvaValue(*fo, pv_, *ttBuilder, *clusterTools);

        leptonTree_->leptonSelection_ |= (LeptonTree::PassEleFO);
        if (smurfselections::passElectronID2011(fo, pv_, thebs.position(), conversions_h, mvaValue))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleID);
        if (smurfselections::passElectronIso2011(fo, pfCandCollection_, pv_))
            leptonTree_->leptonSelection_ |= (LeptonTree::PassEleIso);

        // pf isolation
        leptonTree_->el_pfchiso03_   = (*(*eleIsoVals03_)[0])[fo];
        leptonTree_->el_pfemiso03_   = (*(*eleIsoVals03_)[1])[fo];
        leptonTree_->el_pfnhiso03_   = (*(*eleIsoVals03_)[2])[fo];
        leptonTree_->el_pfchiso04_   = (*(*eleIsoVals04_)[0])[fo];
        leptonTree_->el_pfemiso04_   = (*(*eleIsoVals04_)[1])[fo];
        leptonTree_->el_pfnhiso04_   = (*(*eleIsoVals04_)[2])[fo];
        leptonTree_->el_radiso03_    = smurfselections::getElectronRadialIsolation(*fo, pfNoPUCandCollection_, 0.3);
        leptonTree_->el_radiso04_    = smurfselections::getElectronRadialIsolation(*fo, pfNoPUCandCollection_, 0.4);
        leptonTree_->el_iso_         = smurfselections::electronIsoValuePF(pfCandCollection_, *fo, pv_, 0.4, 1.0, 0.1, 0.07, 0.025, -999., 0);
        leptonTree_->el_ea04_        = smurfselections::getEGammaEffectiveArea(fo->superCluster()->eta(), smurfselections::EGAMMA2012_04);

        // cut based ele id
        leptonTree_->vetoId_    = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::VETO, fo, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);
        leptonTree_->looseId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE, fo, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);
        leptonTree_->mediumId_  = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::MEDIUM, fo, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);
        leptonTree_->tightId_   = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::TIGHT, fo, conversions_h, thebs, vtx_h_, leptonTree_->el_pfchiso03_, leptonTree_->el_pfemiso03_, leptonTree_->el_pfnhiso03_, rhoIsoAll_);

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

        // probe trigger matching
        if (iEvent.isRealData()) {
            // leptons
            const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();
            objectMatchTrigger(iEvent, iSetup, muTriggers_, allObjects, fo->p4(), muTriggerPrescalesProbe_);
            getTriggerVersions(iEvent, iSetup, muTriggers_, muTriggerVersions_);
            // photons
            getTriggerPrescales(iEvent, iSetup, photonTriggers_, photonTriggerPrescales_);
            getTriggerVersions(iEvent, iSetup, photonTriggers_, photonTriggerVersions_);
        }

        leptonTree_->eventSelection_    |= LeptonTree::QCDFakeMu;
        leptonTree_->probe_              = fo->p4();
        leptonTree_->qProbe_             = fo->charge();
        leptonTree_->mt_                 = smurfutilities::Mt(met, fo->p4().Pt(), reco::deltaPhi(met_h->front().phi(), fo->p4().Phi()));

        leptonTree_->mu_pfemiso04_      = fo->pfIsolationR04().sumPhotonEt;
        leptonTree_->mu_pfchiso04_      = fo->pfIsolationR04().sumChargedHadronPt;
        leptonTree_->mu_pfnhiso04_      = fo->pfIsolationR04().sumNeutralHadronEt;
        leptonTree_->mu_radiso03_       = smurfselections::getMuonRadialIsolation(*fo, pfNoPUCandCollection_, 0.3);
        leptonTree_->mu_radiso04_       = smurfselections::getMuonRadialIsolation(*fo, pfNoPUCandCollection_, 0.4);
        leptonTree_->mu_iso_            = smurfselections::muonIsoValuePF(pfCandCollection_, *fo, pv_, 0.3, 1.0, 0.1, 0);
        leptonTree_->mu_eaem04_         = smurfselections::getMuonEffectiveArea(fo->p4().eta(), smurfselections::MUON2012_EM04);
        leptonTree_->mu_eanh04_         = smurfselections::getMuonEffectiveArea(fo->p4().eta(), smurfselections::MUON2012_NH04);

        const reco::GsfElectronCollection nullEls;
        const reco::MuonCollection nullMus;
        leptonTree_->muonHWW2011IDIsoMVA_       = reader_muonHWW2011IDIsoMVA_->MVAValue(&*fo, pv_, ttBuilder, rhoIsoAll_, false);
        leptonTree_->muonHZZ2012IDMVA_          = reader_muonHZZ2012IDMVA_->mvaValue(*fo, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
        leptonTree_->muonHZZ2012IDIsoRingsMVA_  = reader_muonHZZ2012IDIsoRingsMVA_->mvaValue(*fo, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
        leptonTree_->muonHZZ2012IsoRingsMVA_    = reader_muonHZZ2012IsoRingsMVA_->mvaValue(*fo, pv_, pfCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC, nullEls, nullMus);
        leptonTree_->muonHZZ2012IsoDRMVA_       = reader_muonHZZ2012IsoDRMVA_->mvaValue(*fo, pfNoPUCandCollection_, rhoIsoAll_, MuonEffectiveArea::kMuEAFall11MC);

        leptonTree_->leptonSelection_ |= (LeptonTree::PassMuFO);
        if (smurfselections::passMuonID2011(fo, pv_))                     leptonTree_->leptonSelection_ |= (LeptonTree::PassMuID);
        if (smurfselections::passMuonIso2011(fo, pfCandCollection_, pv_)) leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIso);
        if (smurfselections::passMuonIsPOGTight(fo, pv_))                 leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGTight);
        if (smurfselections::passMuonIsPOGSoft(fo, pv_))                  leptonTree_->leptonSelection_ |= (LeptonTree::PassMuIsPOGSoft);
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
    std::vector<std::pair<reco::PFJet, float> > jets15 = smurfselections::goodJets(iEvent, iSetup, jets_h_, cand1, jetCorrector, 15.);
    if (jets15.size() > 0) {
        leptonTree_->jet1_          = jets15.at(0).first.p4() * jets15.at(0).second;
        leptonTree_->dPhiProbeJet1_ = reco::deltaPhi(cand1.p4().Phi(), jets15.at(0).first.p4().Phi());
    }
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
