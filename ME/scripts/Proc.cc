#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Proc.h"


// function to initialize BR 
void Proc::initBR() {
  if (proc_ == TVar::WW || (proc_ <= TVar::HWW300 && proc_ >= TVar::HWW120 )) {
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

  TFile *f = TFile::Open(TString(inputDir_ + TVar::SmurfProcessName(proc_).Data() + "_ME.root"));
  if(f==0x0) {
    std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) << " WARNING!" << TString(inputDir_ + TVar::SmurfProcessName(proc_).Data()+"_ME.root") << " doesn't exist\n";
    for (int j = 0; j < kNDilep; j++)
      yield_[j] = 0.0;
    return;
  }
  
  std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) 
       <<" from " << TString(inputDir_ + TVar::SmurfProcessName(proc_).Data()+"_ME.root") << std::endl; 
  
  TTree *tree = (TTree*) f->Get("tree");
  double tot_yield = 0;
  std::cout<<"Yields from file "<< TVar::SmurfProcessName(proc_)  <<"_ME.root: mm/me/em/ee:  ";
  
  for (int j=0; j<kNDilep; j++){ 
    TH1F *tmp = new TH1F("tmp", "tmp", 20, 0, 20);
    if(j==0 || j==2)
      tree->Project("tmp", "dPhi", Form("scale1fb*(type==%i&&dilep.mass()<%f)", j, massCut_));
    if(j==1 || j==3)
      tree->Project("tmp", "dPhi", Form("scale1fb*((type==%i&&lep2.pt()>15)&&dilep.mass()<%f)", j, massCut_));
    
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
    MCFMXsec_ = 2.983;
    break;

  case (TVar::Wp_1jet):
    NLOXsec_ = 31314.0;
    MCFMXsec_ = 18788.0;
    break;

  case (TVar::Wm_1jet):
    NLOXsec_ = 31314.0;
    MCFMXsec_ =  12525.0;
    break;

  case (TVar::HWW120):
    NLOXsec_ = 0.249642;
    MCFMXsec_ = 0.07533;
    break;

  case (TVar::HWW130):
    NLOXsec_ = 0.452090;
    MCFMXsec_ = 0.1451;
    break;

  case (TVar::HWW140):
    NLOXsec_ = 0.641773;
    MCFMXsec_ = 0.2165;
    break;

  case (TVar::HWW150):
    NLOXsec_ = 0.770471;
    MCFMXsec_ = 0.2671;
    break;

  case (TVar::HWW160):
    NLOXsec_ = 0.866443;
    MCFMXsec_ = 0.3034;
    break;

  case (TVar::HWW170):
    NLOXsec_ = 0.782962;
    MCFMXsec_ = 0.2870;
    break;

  case (TVar::HWW180):
    NLOXsec_ = 0.659328;
    MCFMXsec_ = 0.2465;
    break;

  case (TVar::HWW190):
    NLOXsec_ = 0.486486;
    MCFMXsec_ = 0.1849;
    break;

  case (TVar::HWW200):
    NLOXsec_ = 0.408305;
    MCFMXsec_ = 0.1570;
    break;

  case (TVar::HWW210):
    NLOXsec_ = 0.358465;
    MCFMXsec_ =  0.1380;
    break;
    
  case (TVar::HWW220):
    NLOXsec_ = 0.321398;
    MCFMXsec_ = 0.1232;
    break;

  case (TVar::HWW230):
    NLOXsec_ = 0.290454;
    MCFMXsec_ = 0.1107;
    break;

  case (TVar::HWW250):
    NLOXsec_ = 0.243724;
    MCFMXsec_ = 0.0906;
    break;
    
  case (TVar::HWW300):
    NLOXsec_ = 0.243724;
    MCFMXsec_ = 0.05758;
    break;
    
  default:
    NLOXsec_ = 0.0;
    MCFMXsec_ = 0.0;
    break;
  }
    
}
