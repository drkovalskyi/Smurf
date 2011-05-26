#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"

#include <vector>
#include <math.h>
#include <iostream>

// Matrix Element related
#include "../TVar.hh"

const unsigned int kNDilep = 4;


class Sample {

public:
  
  Sample(TVar::Process LRProcess, TString datasample, TChain *chain, Color_t color, bool stack, float lumi);
  ~Sample();
    
  const TString name() { return datasample_; }
  const TVar::Process GetLRProcess() { return LRProcess_;}
  const bool Stack() {return stack_; }
  const TChain *GetChain() {return chain_; }
  TH1F* GetHistogram(TString var, Int_t nBins, Float_t binMin, Float_t binMax);
  
private:
  
  //
  // member functions
  //
  
  //
  // member data
  //
  
  TVar::Process LRProcess_;
  TString datasample_;
  TChain *chain_; 
  bool stack_;
  Color_t color_;
  Float_t lumi_;
};

