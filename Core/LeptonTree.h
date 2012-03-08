#ifndef LeptonTree_H
#define LeptonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

#include <cmath>
#include "assert.h"

//
// Ntuple structure:
//
// Ntuple content:



class LeptonTree {
 public:

  /// bit map
  /// DON'T CHANGE ORDER
  enum LeptonSelection {

    PassEleSC                       = 1UL<<0,  // 
    PassEleReco                     = 1UL<<1,  // 

    PassEleFO                       = 1UL<<2,  // 
    PassEleID                       = 1UL<<3,  // 
    PassEleIso                      = 1UL<<4,  //

    PassEleTrigDoubleEleLeadingLeg  = 1UL<<5,  // 
    PassEleTrigDoubleEleTrailingLeg = 1UL<<6,  // 
    PassEleTrigSingleEle            = 1UL<<7,  // 
    PassEleTrigMuEGEleLeadingLeg    = 1UL<<8,  // 
    PassEleTrigMuEGEleTrailingLeg   = 1UL<<9,  // 

    PassMuCTFTrack                  = 1UL<<10, // 
    PassMuGlobalOrTrackerMuon       = 1UL<<11, // 

    PassMuFO                        = 1UL<<12, //
    PassMuID                        = 1UL<<13, //
    PassMuIso                       = 1UL<<14, //

    PassMuTrigDoubleMuLeadingLeg    = 1UL<<15,  // 
    PassMuTrigDoubleMuTrailingLeg   = 1UL<<16,  // 
    PassMuTrigSingleMu              = 1UL<<17,  // 
    PassMuTrigMuEGMuLeadingLeg      = 1UL<<18,  // 
    PassMuTrigMuEGMuTrailingLeg     = 1UL<<19,  // 
    

  };

  enum EventSelection {
    ZeeTagAndProbe                  = 1UL<<0,  // 
    ZmmTagAndProbe                  = 1UL<<1,  // 
    OniaEETagAndProbe               = 1UL<<2,  // 
    OniaMuMuTagAndProbe             = 1UL<<3,  // 
    QCDFakeEle8                     = 1UL<<4,  // 
    QCDFakeEle8VL                   = 1UL<<5,  // 
    QCDFakeEle17VL                  = 1UL<<6,  // 
    QCDFakeEle8VLJet40              = 1UL<<7,  // 
    QCDFakeEle8WithTrkIDIso         = 1UL<<8,  // 
    QCDFakeMu8                      = 1UL<<9,  // 
    QCDFakeMu15                     = 1UL<<10,  // 
    ZJetsFakeEleSelection           = 1UL<<11,  // 
    ZJetsFakeMuSelection            = 1UL<<12,  // 
  };


  /// variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   nvtx_;
  unsigned int   leptonSelection_;
  unsigned int   eventSelection_;
  float          rho_;
  float          tagAndProbeMass_;
  float          pt_;
  float          eta_;
  float          phi_;
  float          q_;
  float          scale1fb_;
  float          leadingAwayJetPt_;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  LeptonTree(){}
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
    tree_->Branch("nvtx"             , &nvtx_             ,   "nvtx/i");
    tree_->Branch("leptonSelection"  , &leptonSelection_  ,   "leptonSelection/i");
    tree_->Branch("eventSelection"   , &eventSelection_   ,   "eventSelection/i");
    tree_->Branch("rho"              , &rho_              ,   "rho/F");
    tree_->Branch("tagAndProbeMass"  , &tagAndProbeMass_  ,   "tagAndProbeMass/F");
    tree_->Branch("pt"               , &pt_               ,   "pt/F");
    tree_->Branch("eta"              , &eta_              ,   "eta/F");
    tree_->Branch("phi"              , &phi_              ,   "phi/F");
    tree_->Branch("q"                , &q_                ,   "q/F");
    tree_->Branch("scale1fb"           , &scale1fb_           ,   "scale1fb/F");
    tree_->Branch("leadingAwayJetPt"        , &leadingAwayJetPt_        ,   "leadingAwayJetPt/F");
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
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("leptonSelection",          &leptonSelection_);
    tree_->SetBranchAddress("eventSelection",          &eventSelection_);
    tree_->SetBranchAddress("rho",          &rho_);
    tree_->SetBranchAddress("tagAndProbeMass",          &tagAndProbeMass_);
    tree_->SetBranchAddress("pt",          &pt_);
    tree_->SetBranchAddress("eta",          &eta_);
    tree_->SetBranchAddress("phi",          &phi_);
    tree_->SetBranchAddress("q",          &q_);
    tree_->SetBranchAddress("scale1fb",          &scale1fb_);
    tree_->SetBranchAddress("leadingAwayJetPt",          &leadingAwayJetPt_);

    gErrorIgnoreLevel = currentState;
  }

  /// get a built in type variable by name
  double Get(std::string value);
  /// compare two LeptonTrees for a given event on a given level of precision; 
  /// returns the variables that failed the comparison 
  std::vector<std::string> Compare(LeptonTree* value, double prec=0.005);

  private:

}; 

inline void 
LeptonTree::InitVariables(){
  // create list of available variables
  if(variables_.empty()){
    //make sure that this is only done once
    variables_.push_back(std::string("event" ));
    variables_.push_back(std::string("run"   ));
    variables_.push_back(std::string("lumi"  ));
    variables_.push_back(std::string("nvtx"  ));
    variables_.push_back(std::string("leptonSelection" ));
    variables_.push_back(std::string("eventSelection"  ));
    variables_.push_back(std::string("rho" ));
    variables_.push_back(std::string("tagAndProbeMass" ));
    variables_.push_back(std::string("pt"  ));
    variables_.push_back(std::string("eta" ));
    variables_.push_back(std::string("phi" ));
    variables_.push_back(std::string("q"   ));
    variables_.push_back(std::string("scale1fb"  ));
    variables_.push_back(std::string("leadingAwayJetPt" ));

  }
  // inizialize variables
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  nvtx_          = 0;
  leptonSelection_  = 0;
  eventSelection_   = 0;
  rho_              = 0;
  tagAndProbeMass_  = 0;
  pt_               = 0;
  eta_              = -999;
  phi_              = -999;
  q_                = 0;
  scale1fb_           = 0;
  leadingAwayJetPt_        = 0;

}

inline double
LeptonTree::Get(std::string value)
{
  if(value=="event"           ) { return this->event_;	          }
  if(value=="run"             ) { return this->run_;	          }
  if(value=="lumi"            ) { return this->lumi_;	          }
  if(value=="nvtx"            ) { return this->nvtx_;	          }
  if(value=="leptonSelection" ) { return this->leptonSelection_;  }
  if(value=="eventSelection"  ) { return this->eventSelection_;   }
  if(value=="rho"             ) { return this->rho_;	          }
  if(value=="tagAndProbeMass" ) { return this->tagAndProbeMass_;  }
  if(value=="pt"              ) { return this->pt_;	          }
  if(value=="eta"             ) { return this->eta_;	          }
  if(value=="phi"             ) { return this->phi_;              }
  if(value=="q"               ) { return this->q_;	          }
  if(value=="scale1fb"          ) { return this->scale1fb_;	          }
  if(value=="leadingAwayJetPt" ) { return this->leadingAwayJetPt_; }

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
