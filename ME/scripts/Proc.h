#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "../TVar.hh"
#include "TString.h"

const int kNDilep=4;

enum Analysis {
  HWWANALYSIS = 0,
  HZZANALYSIS = 1,
};

class Proc {
  
 public:
  
  Proc(TVar::Process proc, float lumi, float massCut, TString inputDir, Analysis analysis, float metCut, float dphiCut, float mtCut)
    : proc_(proc), lumi_(lumi), massCut_(massCut), inputDir_(inputDir), analysis_(analysis), metCut_(metCut), dphiCut_(dphiCut), mtCut_(mtCut)
  {
    // call function to set total cross-sections
    initTotXsec();
    // call function to set BR
    initBR();
    // call function to init yields
    initYields();
    // call function to calculate acceptance
    CalculateAcceptance();
  }
  
  ~Proc() {};
  
  const float & GetAcceptance(unsigned int i) { return acceptance_[i]; }
  const float & GetYield(unsigned int i)      {return yield_[i]; }
  const float & GetMCFMXsec() { return MCFMXsec_; }
  
 private:
  
  TVar::Process proc_;
  float lumi_;
  float massCut_; // only for hww analysis
  TString inputDir_;
  float NLOXsec_;
  float MCFMXsec_;
  float BR_[kNDilep];
  float yield_[kNDilep];
  float acceptance_[kNDilep];
  int analysis_;
  float metCut_; // only for hzz analysis
  float dphiCut_; // only for hzz analysis
  float mtCut_;  // only for hzz analysis

  // function to set NLO and MCFM toal cross-sections
  void initTotXsec();  
  // function to set yields
  void initYields();
  // function to set BR
  void initBR();
  // function to calculate acceptance from yields
  void CalculateAcceptance();
};
