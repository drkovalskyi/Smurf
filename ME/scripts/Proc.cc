#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Proc.h"


// function to initialize BR 
void Proc::initBR() {
  if (proc_ == TVar::WW || (proc_ <= TVar::HWW300 && proc_ >= TVar::HWW115 )) {
    BR_[0] = (1.0+0.1736)*(1.0+0.1736)/9;
    BR_[1] = (1.0+0.1736)*(1.0+0.1784)/9;
    BR_[2] = (1.0+0.1736)*(1.0+0.1784)/9;
    BR_[3] = (1.0+0.1784)*(1.0+0.1784)/9;
  }
  else if (proc_ == TVar::Wp_1jet || proc_ == TVar::Wm_1jet) {
    BR_[0] = (1.0+0.1736)/3;
    BR_[1] = (2.0+0.1736+0.1784)/6;
    BR_[2] = (2.0+0.1736+0.1784)/6;
    BR_[3] = (1.0+0.1784)/3;
  }
  else if (proc_ == TVar::ZZ || ( proc_ >= TVar::HZZ200 && proc_ <= TVar::HZZ400) || proc_==TVar::WZ) {
    BR_[0] = (1.0+0.1736*0.1736)/3;
    BR_[1] = (0.1736*0.1784)/3;
    BR_[2] = (0.1736*0.1784)/3;
    BR_[3] = (1.0+0.1784*0.1784)/3;
  }
  else {
    std::cout << "Initialize BR() for process " << TVar::SmurfProcessName(proc_) << " WARNING! Unsupported processes...setting BR to 0...\n";
    for (int j = 0; j < kNDilep; j++)
      BR_[j] = 0.0;
  }
  return;
}

// function to calculate the acceptance
void Proc::CalculateAcceptance(){
  for (int i = 0; i < kNDilep; i++) {
    float denominator = lumi_ * BR_[i] * NLOXsec_;
    if (denominator!=0) acceptance_[i]= yield_[i]/denominator; 
    else  
      acceptance_[i] = 1E+20;
  }
}

// function to calculate acceptance from yields
void Proc::initYields() {
  TString fileName = TVar::SmurfProcessName(proc_);
  if (proc_ == TVar::WW) fileName = "qqww"; 
  TFile *f = TFile::Open(TString(inputDir_ + fileName.Data() + ".root"));
      
  if(f==0x0) {
    std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) << "WARNING!" << TString(inputDir_ + fileName + ".root") << " doesn't exist\n";
    for (int j = 0; j < kNDilep; j++)
      yield_[j] = 0.0;
    return;
  }
  
  std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) 
       <<" from " << TString(inputDir_ + fileName + ".root") << std::endl; 
  
  TTree *tree = (TTree*) f->Get("tree");
  double tot_yield = 0;
  std::cout<<"Yields from file "<< fileName << ".root: mm/me/em/ee:  ";
  
  for (int j=0; j<kNDilep; j++){ 
    TH1F *tmp = new TH1F("tmp", "tmp", 20, 0, 20);
    tree->Project("tmp", "dPhi", Form("scale1fb*(type==%i&&dilep.mass()<%f&&mt>80.)", j, massCut_));
    if (tmp!= 0x0) yield_[j] = tmp->Integral(0, 9999); 
    else yield_[j]=0;
    tot_yield += yield_[j];
    std::cout << Form("%.3f",yield_[j])<<" " ;
    delete tmp;
  } 
  std::cout << "; total yield = " << Form("%.3f", tot_yield) << "\n"; 
  
  f->Close();
}


void Proc::initTotXsec() {
  
  switch (proc_) {
    
  case (TVar::WW): 
    NLOXsec_ = 4.5;
    MCFMXsec_ = 0.333;
    break;

  case (TVar::Wp_1jet):
    NLOXsec_ = 31314.0;
    MCFMXsec_ = 18788.0/3.0;
    break;

  case (TVar::Wm_1jet):
    NLOXsec_ = 31314.0;
    MCFMXsec_ =  12525.0/3.0;
    break;

  case (TVar::HWW115):
    NLOXsec_ = 0.165009;
    MCFMXsec_ = 0.005534;
    break;

  case (TVar::HWW120):
    NLOXsec_ = 0.249642;
    MCFMXsec_ = 0.008370;
    break;

  case (TVar::HWW130):
    NLOXsec_ = 0.452090;
    MCFMXsec_ = 0.016125;
    break;

  case (TVar::HWW140):
    NLOXsec_ = 0.641773;
    MCFMXsec_ = 0.024054;
    break;

  case (TVar::HWW150):
    NLOXsec_ = 0.770471;
    MCFMXsec_ = 0.029674;
    break;

  case (TVar::HWW160):
    NLOXsec_ = 0.866443;
    MCFMXsec_ = 0.033716;
    break;

  case (TVar::HWW170):
    NLOXsec_ = 0.782962;
    MCFMXsec_ = 0.031889;
    break;

  case (TVar::HWW180):
    NLOXsec_ = 0.659328;
    MCFMXsec_ = 0.027389;
    break;

  case (TVar::HWW190):
    NLOXsec_ = 0.486486;
    MCFMXsec_ = 0.020545;
    break;

  case (TVar::HWW200):
    NLOXsec_ = 0.408305;
    MCFMXsec_ = 0.017441;
    break;

  case (TVar::HWW210):
    NLOXsec_ = 0.358465;
    MCFMXsec_ =  0.015340;
    break;
    
  case (TVar::HWW220):
    NLOXsec_ = 0.321398;
    MCFMXsec_ = 0.013686;
    break;

  case (TVar::HWW230):
    NLOXsec_ = 0.290454;
    MCFMXsec_ = 0.0123;
    break;

  case (TVar::HWW250):
    NLOXsec_ = 0.243724;
    MCFMXsec_ = 0.010063;
    break;
    
  case (TVar::HWW300):
    NLOXsec_ = 0.175652;
    MCFMXsec_ = 0.006398;
    break;

  case (TVar::HWW350):
    NLOXsec_ = 0.160052;
    MCFMXsec_ = 0.004251;
    break;

  case (TVar::HWW400):
    NLOXsec_ = 0.124330;
    MCFMXsec_ = 0.002464;
    break;

  case (TVar::HWW450):
    NLOXsec_ = 0.078433;
    MCFMXsec_ = 0.001637;
    break;

  case (TVar::HWW500):
    NLOXsec_ = 0.048702;
    MCFMXsec_ = 0.00156;
    break;

  case (TVar::HWW550):
    NLOXsec_ = 0.030364;
    MCFMXsec_ = 0.000848;
    break;

  case (TVar::HWW600):
    NLOXsec_ = 0.019184;
    MCFMXsec_ = 0.000633;
    break;
    
 case (TVar::HZZ250):
   NLOXsec_ = 0.03974;
   MCFMXsec_ = 0.004906;
   break;

 case (TVar::HZZ300):
   NLOXsec_ = 0.02999;
   MCFMXsec_ = 0.003274;
    break;

 case (TVar::HZZ400):
   NLOXsec_ = 0.0221156;
   MCFMXsec_ = 0.001323;
    break;

 case (TVar::ZZ):
    NLOXsec_ = 0.238;
    MCFMXsec_ = 0.06644;
    break;
    
 case (TVar::WZ):
    NLOXsec_ = 0.202;
    MCFMXsec_ = 0.060;
    break;

  default:
    NLOXsec_ = 0.0;
    MCFMXsec_ = 0.0;
    break;
  }
    
}
