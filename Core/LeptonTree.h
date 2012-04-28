#ifndef LeptonTree_H
#define LeptonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

#include "Math/LorentzVector.h"

#include <cmath>
#include "assert.h"

//
// Ntuple structure:
//
// Ntuple content:



class LeptonTree {

    public:
        /// float doesn't have dictionary by default, so use double
        typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

        /// bit map
        /// DON'T CHANGE ORDER
        enum LeptonSelection {
            PassEleSC                       = 1UL<<0,  // 
            PassEleReco                     = 1UL<<1,  // 
            PassEleFO                       = 1UL<<2,  // 
            PassEleID                       = 1UL<<3,  // 
            PassEleIso                      = 1UL<<4,  //
            PassMuCTFTrack                  = 1UL<<5,  // 
            PassMuGlobalOrTrackerMuon       = 1UL<<6,  // 
            PassMuIsPF                      = 1UL<<7,  //
            PassMuIsPOGTight                = 1UL<<8,  //
            PassMuIsPOGSoft                 = 1UL<<9,  //
            PassMuFO                        = 1UL<<10, //
            PassMuID                        = 1UL<<11, //
            PassMuIso                       = 1UL<<12, //
            PassPhotonID                    = 1UL<<13, //
            PassPhotonIso                   = 1UL<<14  //
        };

        enum EventSelection {
            ZeeTagAndProbe                  = 1UL<<0,  // 
            ZmmTagAndProbe                  = 1UL<<1,  // 
            OniaEETagAndProbe               = 1UL<<2,  // 
            OniaMuMuTagAndProbe             = 1UL<<3,  // 
            QCDFakeEle                      = 1UL<<4,  //
            QCDFakeMu                       = 1UL<<5,  // 
            ZJetsFakeEleSelection           = 1UL<<6,  // 
            ZJetsFakeMuSelection            = 1UL<<7,  // 
            PhotonSelection                 = 1UL<<8   //
        };

        /// variables
        unsigned int   event_;
        unsigned int   run_;
        unsigned int   lumi_;
        float          rnd_;
        unsigned int   nvtx_;
        unsigned int   npu_;
        unsigned int   npuPlusOne_;
        unsigned int   npuMinusOne_;
        unsigned int   leptonSelection_;
        unsigned int   eventSelection_;
        float          rhoIsoAll_;
        float          rhoIsoNeutral_;
        float          tagAndProbeMass_;
        LorentzVector  tag_;
        LorentzVector  probe_;
        int   qTag_;
        int   qProbe_;
        float          scale1fb_;
        LorentzVector  jet1_;
        LorentzVector  jet2_;
        LorentzVector  jet3_;
        float          met_;
        float          metPhi_;
        float          trackMet_;
        float          trackMetPhi_;
        unsigned int   njets_;
        float          hltPrescale_;
        float          sumet_;
        float          metSig_;
        float          mt_;
        float          dPhiProbeJet1_;

        //
        // MVA IDs
        //

        // HWW MVAs
        float electronHWW2011MVA_;
        float electronHWW2011IDIsoMVA_;
        float muonHWW2011IDIsoMVA_;

        // POG MVAs
        float egammaPOG2012MVA_;

        // other MVAs
        float muonHZZ2012MVA_;

        //
        // for electron studies with new data
        //    

        // masks for 2012 cut based ID
        unsigned int vetoId_;
        unsigned int looseId_;
        unsigned int mediumId_;
        unsigned int tightId_;

        // input variables
        float sceta_;
        float scenergy_;
        bool chargesAgree_;
        float pfmva_;
        float ooemoop_;
        float eopin_;
        float fbrem_;
        float detain_;
        float dphiin_;
        float hoe_;
        float hoetow_;
        float sieie_;
        float d0vtx_;
        float dzvtx_;
        bool vfitprob_;
        float mhit_;
        float ecaliso_;
        float hcaliso_;
        float trkiso_;

        float el_pfemiso03_;
        float el_pfchiso03_;
        float el_pfnhiso03_;
        float el_pfemiso04_;
        float el_pfchiso04_;
        float el_pfnhiso04_;
        float el_ea04_;

        float mu_pfemiso04_;
        float mu_pfchiso04_;
        float mu_pfnhiso04_;
        float mu_eaem04_;
        float mu_eanh04_;

    public:
        /// this is the main element
        TTree *tree_;
        TFile *f_;

        /// hold the names of variables to facilitate things (filled during Init)
        std::vector<std::string> variables_;

        /// default constructor  
        LeptonTree() : tagPtr_(&tag_), probePtr_(&probe_) {}
        /// default destructor
        ~LeptonTree(){ 
            if (f_) f_->Close();  
        };

        /// initialize varibles and fill list of available variables
        void InitVariables();

        /// load a LeptonTree
        void LoadTree(const char* file){
            f_ = TFile::Open(file);
            assert(f_);
            tree_ = dynamic_cast<TTree*>(f_->Get("leptons"));
            assert(tree_);
        }

        /// create a LeptonTree
        void CreateTree(){
            tree_ = new TTree("leptons","Lepton ntuples");
            f_ = 0;
            InitVariables();
            //book the branches
            tree_->Branch("event"            , &event_            ,   "event/i");
            tree_->Branch("run"              , &run_              ,   "run/i");
            tree_->Branch("lumi"             , &lumi_             ,   "lumi/i");
            tree_->Branch("rnd"              , &rnd_              ,   "rnd/F");
            tree_->Branch("nvtx"             , &nvtx_             ,   "nvtx/i");
            tree_->Branch("npu"              , &npu_              ,   "npu/i");
            tree_->Branch("npuPlusOne"       , &npuPlusOne_       ,   "npuPlusOne/i");
            tree_->Branch("npuMinusOne"      , &npuMinusOne_      ,   "npuMinusOne/i");
            tree_->Branch("leptonSelection"  , &leptonSelection_  ,   "leptonSelection/i");
            tree_->Branch("eventSelection"   , &eventSelection_   ,   "eventSelection/i");
            tree_->Branch("rhoIsoAll"              , &rhoIsoAll_              ,   "rhoIsoAll/F");
            tree_->Branch("rhoIsoNeutral"              , &rhoIsoNeutral_              ,   "rhoIsoNeutral/F");
            tree_->Branch("tagAndProbeMass"  , &tagAndProbeMass_  ,   "tagAndProbeMass/F");
            tree_->Branch("tag"              , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tagPtr_);
            tree_->Branch("probe"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &probePtr_);
            tree_->Branch("qTag"             , &qTag_             ,   "qTag/I");
            tree_->Branch("qProbe"           , &qProbe_           ,   "qProbe/I");
            tree_->Branch("scale1fb"         , &scale1fb_         ,   "scale1fb/F");
            tree_->Branch("jet1"             , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet1_);
            tree_->Branch("jet2"             , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet2_);
            tree_->Branch("jet3"             , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet3_);
            tree_->Branch("met"              , &met_              ,   "met/F");
            tree_->Branch("metPhi"           , &metPhi_           ,   "metPhi/F");
            tree_->Branch("trackMet"         , &trackMet_         ,   "trackMet/F");
            tree_->Branch("trackMetPhi"      , &trackMetPhi_      ,   "trackMetPhi/F");
            tree_->Branch("njets"            , &njets_            ,   "njets/i");
            tree_->Branch("hltPrescale"      , &hltPrescale_      ,   "hltPrescale/F");
            tree_->Branch("sumet"            , &sumet_            ,   "sumet/F");
            tree_->Branch("metSig"           , &metSig_           ,   "metSig/F");
            tree_->Branch("mt"            , &mt_            ,   "mt/F");
            tree_->Branch("dPhiProbeJet1"            , &dPhiProbeJet1_            ,   "dPhiProbeJet1/F");

            tree_->Branch("electronHWW2011MVA"      , &electronHWW2011MVA_      ,   "electronHWW2011MVA/F");
            tree_->Branch("electronHWW2011IDIsoMVA"      , &electronHWW2011IDIsoMVA_      ,   "electronHWW2011IDIsoMVA/F");
            tree_->Branch("muonHWW2011MVA"          , &muonHWW2011IDIsoMVA_          ,   "muonHWW2011MVA/F");
            tree_->Branch("egammaPOG2012MVA"        , &egammaPOG2012MVA_        ,   "egammaPOG2012MVA/F");
            tree_->Branch("muonHZZ2012MVA"        , &muonHZZ2012MVA_        ,   "muonHZZ2012MVA/F");

            // extra branches for electron ID studies
            tree_->Branch("vetoId"           , &vetoId_           ,   "vetoId/i");
            tree_->Branch("looseId"          , &looseId_          ,   "looseId/i");
            tree_->Branch("mediumId"         , &mediumId_         ,   "mediumId/i");
            tree_->Branch("tightId"          , &tightId_          ,   "tightId/i");
            tree_->Branch("pfmva"               , &pfmva_               ,   "pfmva/F");
            tree_->Branch("sceta"               , &sceta_               ,   "sceta/F");
            tree_->Branch("scenergy"               , &scenergy_               ,   "scenergy/F");
            tree_->Branch("chargesAgree"               , &chargesAgree_               ,   "chargesAgree/O");
            tree_->Branch("eopin"               , &eopin_               ,   "eopin/F");
            tree_->Branch("ooemoop"               , &ooemoop_               ,   "ooemoop/F");
            tree_->Branch("fbrem"               , &fbrem_               ,   "fbrem/F");
            tree_->Branch("detain"               , &detain_               ,   "detain/F");
            tree_->Branch("dphiin"               , &dphiin_               ,   "dphiin/F");
            tree_->Branch("hoe"               , &hoe_               ,   "hoe/F");
            tree_->Branch("hoetow"               , &hoetow_               ,   "hoetow/F");
            tree_->Branch("sieie"               , &sieie_               ,   "sieie/F");
            tree_->Branch("d0vtx"               , &d0vtx_               ,   "d0vtx/F");
            tree_->Branch("dzvtx"               , &dzvtx_               ,   "dzvtx/F");
            tree_->Branch("vfitprob"               , &vfitprob_               ,   "vfitprob/O");
            tree_->Branch("mhit"               , &mhit_               ,   "mhit/F");
            tree_->Branch("ecaliso"               , &ecaliso_               ,   "ecaliso/F");
            tree_->Branch("hcaliso"               , &hcaliso_               ,   "hcaliso/F");
            tree_->Branch("trkiso"               , &trkiso_               ,   "trkiso/F");
            tree_->Branch("el_pfemiso03"               , &el_pfemiso03_               ,   "el_pfemiso03/F");
            tree_->Branch("el_pfchiso03"               , &el_pfchiso03_               ,   "el_pfchiso03/F");
            tree_->Branch("el_pfnhiso03"               , &el_pfnhiso03_               ,   "el_pfnhiso03/F");
            tree_->Branch("el_pfemiso04"               , &el_pfemiso04_               ,   "el_pfemiso04/F");
            tree_->Branch("el_pfchiso04"               , &el_pfchiso04_               ,   "el_pfchiso04/F");
            tree_->Branch("el_pfnhiso04"               , &el_pfnhiso04_               ,   "el_pfnhiso04/F");
            tree_->Branch("el_ea04"               , &el_ea04_               ,   "mu_ea04/F");

            tree_->Branch("mu_pfemiso04"               , &mu_pfemiso04_               ,   "mu_pfemiso04/F");
            tree_->Branch("mu_pfchiso04"               , &mu_pfchiso04_               ,   "mu_pfchiso04/F");
            tree_->Branch("mu_pfnhiso04"               , &mu_pfnhiso04_               ,   "mu_pfnhiso04/F");
            tree_->Branch("mu_eaem04"               , &mu_eaem04_               ,   "mu_eaem04/F");
            tree_->Branch("mu_eanh04"               , &mu_eanh04_               ,   "mu_eanh04/F");

        }

        // initialze a LeptonTree
        void InitTree(){
            assert(tree_);
            // don't forget to set pointers to zero before you set address
            // or you will fully appreciate that "ROOT sucks" :)
            InitVariables();
            //Set branch address
            Int_t currentState = gErrorIgnoreLevel;
            // gErrorIgnoreLevel = kError;
            gErrorIgnoreLevel = kBreak;
            tree_->SetBranchAddress("event",            &event_);
            tree_->SetBranchAddress("run",              &run_);
            tree_->SetBranchAddress("lumi",             &lumi_);
            tree_->SetBranchAddress("rnd",              &rnd_);
            tree_->SetBranchAddress("nvtx",             &nvtx_);
            tree_->SetBranchAddress("npu",              &npu_);
            tree_->SetBranchAddress("npuPlusOne",       &npuPlusOne_);
            tree_->SetBranchAddress("npuMinusOne",      &npuMinusOne_);
            tree_->SetBranchAddress("leptonSelection",  &leptonSelection_);
            tree_->SetBranchAddress("eventSelection",   &eventSelection_);
            tree_->SetBranchAddress("rhoIsoAll",              &rhoIsoAll_);
            tree_->SetBranchAddress("rhoIsoNeutral",              &rhoIsoNeutral_);
            tree_->SetBranchAddress("tagAndProbeMass",  &tagAndProbeMass_);
            tree_->SetBranchAddress("tag",              &tagPtr_);
            tree_->SetBranchAddress("probe",            &probePtr_);
            tree_->SetBranchAddress("qTag",             &qTag_);
            tree_->SetBranchAddress("qProbe",           &qProbe_);
            tree_->SetBranchAddress("scale1fb",         &scale1fb_);
            tree_->SetBranchAddress("jet1",             &jet1_);
            tree_->SetBranchAddress("jet2",             &jet2_);
            tree_->SetBranchAddress("jet3",             &jet3_);
            tree_->SetBranchAddress("met",              &met_);
            tree_->SetBranchAddress("metPhi",           &metPhi_);
            tree_->SetBranchAddress("trackMet",         &trackMet_);
            tree_->SetBranchAddress("trackMetPhi",      &trackMetPhi_);
            tree_->SetBranchAddress("njets",            &njets_);
            tree_->SetBranchAddress("hltPrescale",      &hltPrescale_);
            tree_->SetBranchAddress("sumet",            &sumet_);
            tree_->SetBranchAddress("metSig",           &metSig_);
            tree_->SetBranchAddress("mt",            &mt_);
            tree_->SetBranchAddress("dPhiProbeJet1",            &dPhiProbeJet1_);

            tree_->SetBranchAddress("electronHWW2011MVA",      &electronHWW2011MVA_);
            tree_->SetBranchAddress("electronHWW2011IDIsoMVA",      &electronHWW2011IDIsoMVA_);
            tree_->SetBranchAddress("muonHWW2011MVA",          &muonHWW2011IDIsoMVA_);
            tree_->SetBranchAddress("egammaPOG2012MVA",        &egammaPOG2012MVA_);
            tree_->SetBranchAddress("muonHZZ2012MVA",   &muonHZZ2012MVA_);

            tree_->SetBranchAddress("vetoId",            &vetoId_);
            tree_->SetBranchAddress("looseId",            &looseId_);
            tree_->SetBranchAddress("mediumId",            &mediumId_);
            tree_->SetBranchAddress("tightId",            &tightId_);
            tree_->SetBranchAddress("pfmva",          &pfmva_);
            tree_->SetBranchAddress("sceta",          &sceta_);
            tree_->SetBranchAddress("scenergy",          &scenergy_);
            tree_->SetBranchAddress("chargesAgree",          &chargesAgree_);
            tree_->SetBranchAddress("eopin",          &eopin_);
            tree_->SetBranchAddress("ooemoop",          &ooemoop_);
            tree_->SetBranchAddress("fbrem",          &fbrem_);
            tree_->SetBranchAddress("detain",          &detain_);
            tree_->SetBranchAddress("dphiin",          &dphiin_);
            tree_->SetBranchAddress("hoe",          &hoe_);
            tree_->SetBranchAddress("hoetow",          &hoetow_);
            tree_->SetBranchAddress("sieie",          &sieie_);
            tree_->SetBranchAddress("d0vtx",          &d0vtx_);
            tree_->SetBranchAddress("dzvtx",          &dzvtx_);
            tree_->SetBranchAddress("vfitprob",          &vfitprob_);
            tree_->SetBranchAddress("mhit",          &mhit_);
            tree_->SetBranchAddress("ecaliso",          &ecaliso_);
            tree_->SetBranchAddress("hcaliso",          &hcaliso_);
            tree_->SetBranchAddress("trkiso",          &trkiso_);
            tree_->SetBranchAddress("el_pfemiso03",          &el_pfemiso03_);
            tree_->SetBranchAddress("el_pfchiso03",          &el_pfchiso03_);
            tree_->SetBranchAddress("el_pfnhiso03",          &el_pfnhiso03_);
            tree_->SetBranchAddress("el_pfemiso04",          &el_pfemiso04_);
            tree_->SetBranchAddress("el_pfchiso04",          &el_pfchiso04_);
            tree_->SetBranchAddress("el_pfnhiso04",          &el_pfnhiso04_);
            tree_->SetBranchAddress("el_ea04",          &el_ea04_);

            tree_->SetBranchAddress("mu_pfemiso04",          &mu_pfemiso04_);
            tree_->SetBranchAddress("mu_pfchiso04",          &mu_pfchiso04_);
            tree_->SetBranchAddress("mu_pfnhiso04",          &mu_pfnhiso04_);
            tree_->SetBranchAddress("mu_eaem04",          &mu_eaem04_);
            tree_->SetBranchAddress("mu_eanh04",          &mu_eanh04_);

            gErrorIgnoreLevel = currentState;
        }

        /// get a built in type variable by name
        double Get(std::string value);
        /// compare two LeptonTrees for a given event on a given level of precision; 
        /// returns the variables that failed the comparison 
        std::vector<std::string> Compare(LeptonTree* value, double prec=0.005);

    private:

        LorentzVector *tagPtr_;
        LorentzVector *probePtr_;

}; 

inline void 
LeptonTree::InitVariables(){
    // create list of available variables
    if(variables_.empty()){
        //make sure that this is only done once
        variables_.push_back(std::string("event"            ));
        variables_.push_back(std::string("run"              ));
        variables_.push_back(std::string("lumi"             ));
        variables_.push_back(std::string("rnd"             ));
        variables_.push_back(std::string("nvtx"             ));
        variables_.push_back(std::string("npu"              ));
        variables_.push_back(std::string("npuPlusOne"       ));
        variables_.push_back(std::string("npuMinusOne"      ));
        variables_.push_back(std::string("leptonSelection"  ));
        variables_.push_back(std::string("eventSelection"   ));
        variables_.push_back(std::string("rhoIsoAll"              ));
        variables_.push_back(std::string("rhoIsoNeutral"              ));
        variables_.push_back(std::string("tagAndProbeMass"  ));
        variables_.push_back(std::string("tag"              ));
        variables_.push_back(std::string("probe"            ));
        variables_.push_back(std::string("qTag"             ));
        variables_.push_back(std::string("qProbe"           ));
        variables_.push_back(std::string("scale1fb"         ));
        variables_.push_back(std::string("jet1"             ));
        variables_.push_back(std::string("jet2"             ));
        variables_.push_back(std::string("jet3"             ));
        variables_.push_back(std::string("met"              ));
        variables_.push_back(std::string("metPhi"           ));
        variables_.push_back(std::string("trackMet"         ));
        variables_.push_back(std::string("trackMetPhi"      ));
        variables_.push_back(std::string("njets"            ));
        variables_.push_back(std::string("hltPrescale"      ));
        variables_.push_back(std::string("sumet"            ));
        variables_.push_back(std::string("metSig"           ));
        variables_.push_back(std::string("mt"            ));
        variables_.push_back(std::string("dPhiProbeJet1"            ));;
        variables_.push_back(std::string("electronHWW2011MVA"      ));
        variables_.push_back(std::string("electronHWW2011IDIsoMVA" ));    
        variables_.push_back(std::string("muonHWW2011MVA"          ));
        variables_.push_back(std::string("egammaPOG2012MVA"        ));
        variables_.push_back(std::string("muonHZZ2012MVA"        ));

        variables_.push_back(std::string("vetoId"));
        variables_.push_back(std::string("looseId"));
        variables_.push_back(std::string("mediumId"));
        variables_.push_back(std::string("tightId"));
        variables_.push_back(std::string("pfmva"));
        variables_.push_back(std::string("sceta"));
        variables_.push_back(std::string("scenergy"));
        variables_.push_back(std::string("chargesAgree"));
        variables_.push_back(std::string("eopin"));  
        variables_.push_back(std::string("ooemoop"));
        variables_.push_back(std::string("fbrem"));
        variables_.push_back(std::string("detain"));
        variables_.push_back(std::string("dphiin"));
        variables_.push_back(std::string("hoe"));
        variables_.push_back(std::string("hoetow"));
        variables_.push_back(std::string("sieie")); 
        variables_.push_back(std::string("d0vtx")); 
        variables_.push_back(std::string("dzvtx")); 
        variables_.push_back(std::string("vfitprob"));
        variables_.push_back(std::string("mhit"));      
        variables_.push_back(std::string("ecaliso"));     
        variables_.push_back(std::string("hcaliso"));    
        variables_.push_back(std::string("trkiso"));
        variables_.push_back(std::string("el_pfemiso03")); 
        variables_.push_back(std::string("el_pfchiso03")); 
        variables_.push_back(std::string("el_pfnhiso03"));    
        variables_.push_back(std::string("el_pfemiso04"));
        variables_.push_back(std::string("el_pfchiso04"));
        variables_.push_back(std::string("el_pfnhiso04"));
        variables_.push_back(std::string("el_ea04"));

        variables_.push_back(std::string("mu_pfemiso04"));
        variables_.push_back(std::string("mu_pfchiso04"));
        variables_.push_back(std::string("mu_pfnhiso04"));
        variables_.push_back(std::string("mu_eaem04"));
        variables_.push_back(std::string("mu_eanh04"));

    }
    // inizialize variables
    event_                = 0;
    run_                  = 0;
    lumi_                 = 0;
    rnd_                  = 0;
    nvtx_                 = 0;
    npu_                  = 0;
    npuPlusOne_           = 0;
    npuMinusOne_          = 0;
    leptonSelection_      = 0;
    eventSelection_       = 0;
    rhoIsoAll_  = -999;
    rhoIsoNeutral_    = -999;
    tagAndProbeMass_      = -999;
    tag_                  = LorentzVector();
    probe_                = LorentzVector();
    qTag_                 = 0;
    qProbe_               = 0;
    scale1fb_             = 0;
    jet1_                 = LorentzVector();
    jet2_                 = LorentzVector();
    jet3_                 = LorentzVector();
    met_                  = -999;
    metPhi_               = -999;
    trackMet_             = -999;
    trackMetPhi_          = -999;
    njets_                = 0;
    hltPrescale_          = 1;
    sumet_                = -999;
    metSig_               = -999;
    mt_                   = -999;
    dPhiProbeJet1_        = -999;

    electronHWW2011MVA_          = -999;
    electronHWW2011IDIsoMVA_     = -999;
    muonHWW2011IDIsoMVA_         = -999;
    egammaPOG2012MVA_            = -999;
    muonHZZ2012MVA_              = -999;

    vetoId_ = 0;
    looseId_ = 0;
    mediumId_ = 0;
    tightId_ = 0;
    pfmva_ = 0;
    sceta_ = 0;
    scenergy_ = 0;
    chargesAgree_ = 0;
    eopin_ = 0;
    ooemoop_ = 0;
    fbrem_ = 0;
    detain_ = 0;
    dphiin_ = 0;
    hoe_ = 0;
    hoetow_ = 0;
    sieie_ = 0;
    d0vtx_ = 0;
    dzvtx_ = 0;
    vfitprob_ = 0;
    mhit_ = 0;
    ecaliso_ = 0;
    hcaliso_ = 0;
    trkiso_ = 0;

    el_pfemiso03_ = 0.;
    el_pfchiso03_ = 0.;
    el_pfnhiso03_ = 0.;
    el_pfemiso04_ = 0.;
    el_pfchiso04_ = 0.;
    el_pfnhiso04_ = 0.;
    el_ea04_      = 0.;

    mu_pfemiso04_ = 0.;
    mu_pfchiso04_ = 0.;
    mu_pfnhiso04_ = 0.;
    mu_eaem04_    = 0.;
    mu_eanh04_    = 0.;

}

    inline double
LeptonTree::Get(std::string value)
{
    if(value=="event"            ) { return this->event_;	           }
    if(value=="run"              ) { return this->run_;	           }
    if(value=="lumi"             ) { return this->lumi_;	           }
    if(value=="rnd"              ) { return this->rnd_;               }
    if(value=="nvtx"             ) { return this->nvtx_;	           }
    if(value=="npu"              ) { return this->npu_;	           }
    if(value=="npuPlusOne"       ) { return this->npuPlusOne_;	   }
    if(value=="npuMinusOne"      ) { return this->npuMinusOne_;	   }
    if(value=="leptonSelection"  ) { return this->leptonSelection_;  }
    if(value=="eventSelection"   ) { return this->eventSelection_;   }
    if(value=="rhoIsoAll"              ) { return this->rhoIsoAll_;	           }
    if(value=="rhoIsoNeutral"              ) { return this->rhoIsoNeutral_;            }
    if(value=="tagAndProbeMass"  ) { return this->tagAndProbeMass_;  }
    if(value=="qTag"             ) { return this->qTag_;	           }
    if(value=="qProbe"           ) { return this->qProbe_;           }
    if(value=="scale1fb"         ) { return this->scale1fb_;         }
    if(value=="met"              ) { return this->met_;              }
    if(value=="metPhi"           ) { return this->metPhi_;           }
    if(value=="trackMet"         ) { return this->trackMet_;         }
    if(value=="trackMetPhi"      ) { return this->trackMetPhi_;      }
    if(value=="njets"            ) { return this->njets_;            }
    if(value=="hltPrescale"      ) { return this->hltPrescale_;      }
    if(value=="sumet"            ) { return this->sumet_;            }
    if(value=="metSig"           ) { return this->metSig_;           }
    if(value=="mt"            ) { return this->mt_;            }
    if(value=="dPhiProbeJet1"            ) { return this->dPhiProbeJet1_;            }

    if(value=="vetoId"              ) { return this->vetoId_;           }
    if(value=="looseId"              ) { return this->looseId_;           }
    if(value=="mediumId"              ) { return this->mediumId_;           }
    if(value=="tightId"              ) { return this->tightId_;           }
    if(value=="pfmva"              ) { return this->pfmva_;           }
    if(value=="sceta"              ) { return this->sceta_;           }
    if(value=="scenergy"              ) { return this->scenergy_;           }
    if(value=="chargesAgree"              ) { return this->chargesAgree_;           }
    if(value=="eopin"              ) { return this->eopin_;           }
    if(value=="ooemoop"              ) { return this->ooemoop_;           }
    if(value=="fbrem"              ) { return this->fbrem_;           }
    if(value=="detain"              ) { return this->detain_;           }
    if(value=="dphiin"              ) { return this->dphiin_;           }
    if(value=="hoe"              ) { return this->hoe_;           }
    if(value=="hoetow"              ) { return this->hoetow_;           }
    if(value=="sieie"              ) { return this->sieie_;           }
    if(value=="d0vtx"              ) { return this->d0vtx_;           }
    if(value=="dzvtx"              ) { return this->dzvtx_;           }
    if(value=="vfitprob"              ) { return this->vfitprob_;           }
    if(value=="mhit"              ) { return this->mhit_;           }
    if(value=="ecaliso"              ) { return this->ecaliso_;           }
    if(value=="hcaliso"              ) { return this->hcaliso_;           }
    if(value=="trkiso"              ) { return this->trkiso_;           }
    if(value=="el_pfemiso03"              ) { return this->el_pfemiso03_;           }
    if(value=="el_pfchiso03"              ) { return this->el_pfchiso03_;           }
    if(value=="el_pfnhiso03"              ) { return this->el_pfnhiso03_;           }
    if(value=="el_pfemiso04"              ) { return this->el_pfemiso04_;           }
    if(value=="el_pfchiso04"              ) { return this->el_pfchiso04_;           }
    if(value=="el_pfnhiso04"              ) { return this->el_pfnhiso04_;           }
    if(value=="el_ea04"              ) { return this->el_ea04_;           }

    if(value=="mu_pfemiso04"              ) { return this->mu_pfemiso04_;           }
    if(value=="mu_pfchiso04"              ) { return this->mu_pfchiso04_;           }
    if(value=="mu_pfnhiso04"              ) { return this->mu_pfnhiso04_;           }
    if(value=="mu_eaem04"              ) { return this->mu_eaem04_;           }
    if(value=="mu_eanh04"              ) { return this->mu_eanh04_;           }

    return -9999.; 
}

inline std::vector<std::string> 
LeptonTree::Compare(LeptonTree* value, double prec){
    std::vector<std::string> fails;
    // this should alway fit with ultimate precision
    if( this->event_ != value->event_ ){ fails.push_back( "event" ); }
    if( this->run_   != value->run_   ){ fails.push_back( "run"   ); }
    if( this->lumi_  != value->lumi_  ){ fails.push_back( "lumi"  ); }

    // check within  (relative) precision
    for(std::vector<std::string>::const_iterator var=variables_.begin(); var!=variables_.end(); ++var){
        if( fabs(Get(*var)-value->Get(*var))/Get(*var)>prec ) fails.push_back(*var);
    }
    return fails;
}


#endif
