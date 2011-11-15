#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Proc.h"
#include "TCut.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


// function to initialize BR 
void Proc::initBR() {
  if (proc_ == TVar::WW || (proc_ <= TVar::HWW600 && proc_ >= TVar::HWW115 )) {
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
  else if (proc_ == TVar::ZZ || ( proc_ >= TVar::HZZ200 && proc_ <= TVar::HZZ600) || proc_==TVar::WZ) {
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
  // hard-coded ntuple names.. fixme..
  TString fileName = TVar::SmurfProcessName(proc_);
  if ( proc_ == TVar::WW && analysis_ == HWWANALYSIS)  fileName = "qqww"; 
  if ( proc_ == TVar::WW && analysis_ == HZZANALYSIS)  fileName = "ww-madgraph";
  if ( proc_ == TVar::WZ && analysis_ == HZZANALYSIS)  fileName = "wz3l-madgraph";
  if ( proc_ == TVar::ZZ && analysis_ == HZZANALYSIS)  fileName = "zz2l-madgraph";
  if ( proc_ >= TVar::HZZ200 && proc_ <= TVar::HZZ600 && analysis_ == HZZANALYSIS) fileName = "gf"+fileName;
 
  TFile *f = TFile::Open(TString(inputDir_ + fileName.Data() + ".root"));
  
  if(f==0x0) {
    std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) << "WARNING!" << TString(inputDir_ + fileName + ".root") << " doesn't exist\n";
    for (int j = 0; j < kNDilep; j++)
      yield_[j] = 0.0;
    return;
  }
  
  TTree *tree = (TTree*) f->Get("tree");
  
  if(tree == 0x0) {
    std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) << "WARNING!" << TString(inputDir_ + fileName + ".root") << " doesn't exist\n";
    for (int j = 0; j < kNDilep; j++)
      yield_[j] = 0.0;
    return;
  }
  
  std::cout << "Proc:initYields for process " << TVar::SmurfProcessName(proc_) 
	    <<" from " << TString(inputDir_ + fileName + ".root") << std::endl; 


  // Initialize the branches to use to calculate LR
  int type_ = 0;
  float scale1fb_ = 1.;
  LorentzVector*  dilep_ = 0;
  float mt_ = 0.;
  float met_ = 0.;
  LorentzVector*  jet1_ = 0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  float metPhi_ = 0.;
  unsigned int nvtx_ = 0;
  float pmet_ = 0.;
  float pTrackMet_ = 0.;
  LorentzVector*  lep2_ = 0;
  
  tree->SetBranchAddress( "nvtx",  &nvtx_);
  tree->SetBranchAddress( "type",  &type_);
  tree->SetBranchAddress( "scale1fb",  &scale1fb_);
  tree->SetBranchAddress( "dilep", &dilep_); 
  tree->SetBranchAddress( "mt",  &mt_);
  tree->SetBranchAddress( "met",  &met_);
  tree->SetBranchAddress( "pmet",  &pmet_);
  tree->SetBranchAddress( "pTrackMet",  &pTrackMet_);
  tree->SetBranchAddress( "jet1", &jet1_);
  tree->SetBranchAddress( "jet2", &jet2_);
  tree->SetBranchAddress( "jet3", &jet3_);
  tree->SetBranchAddress( "metPhi",  &metPhi_); 
  tree->SetBranchAddress( "lep2", &lep2_); 

  // start filling the numbers...
  
  double tot_yield = 0;
  std::cout<<"Yields from file "<< fileName << ".root: mm/me/em/ee:  ";

  TH1F *tmp[4];
  
  for ( int i = 0; i < 4 ; i ++) {
    tmp[i] = new TH1F(Form("tmp_%i", i), Form("tmp_%i", i), 20, 0, 20);
    tmp[i]->Sumw2();
  }
  
  // loop over the events to fill histograms
  for(int ievt = 0; ievt< tree->GetEntries(); ievt++){
    tree->GetEntry(ievt);  
    bool passpresel = true;
    
    if ( analysis_ == HWWANALYSIS) {
      if ( dilep_->M() > massCut_ ) passpresel = false;
      if ( mt_ < 80.) passpresel = false;
      if ( dilep_->Pt() < 45) passpresel = false;
      if ( type_ == 0 || type_ == 3) {
	if ( dilep_->M() < 20) passpresel = false;
	if (TMath::Min(pmet_,pTrackMet_) < (37.+nvtx_/2.)  ) passpresel = false;
	if ( lep2_->Pt() < 15) passpresel = false;
      }
    }
    
    if ( analysis_ == HZZANALYSIS) {
      float dPhi1, dPhi2, dPhi3;
      jet1_->Pt() > 30 ? dPhi1 = acos(cos(metPhi_ - jet1_->Phi())) : dPhi1 = 999.9;
      jet2_->Pt() > 30 ? dPhi2 = acos(cos(metPhi_ - jet2_->Phi())) : dPhi2 = 999.9;
      jet3_->Pt() > 30 ? dPhi3 = acos(cos(metPhi_ - jet3_->Phi())) : dPhi3 = 999.9;
      
      // apply the dphi Cut
      if ( TMath::Min(dPhi1, TMath::Min(dPhi2, dPhi3)) < dphiCut_ ) passpresel = false;
      if ( met_ < metCut_) passpresel = false;

      float termA = sqrt( dilep_->Pt()*dilep_->Pt() + dilep_->M()*dilep_->M() );
      float termB = sqrt( met_*met_ + dilep_->M()*dilep_->M() );
      float newX = dilep_->Px() + met_ * cos(metPhi_);
      float newY = dilep_->Py() + met_ * sin(metPhi_);
      float termC = newX*newX + newY*newY;
      float mt = sqrt( pow(termA + termB, 2) - termC );
      
      if ( mt < mtCut_) passpresel = false;

      if ( dilep_->Pt() < 55) passpresel = false;
    }
    
    if ( passpresel ) 
      tmp[type_]->Fill(mt_, scale1fb_);
    }

  // get the yields from histograms
  
  for ( int i = 0; i < 4; i++) {
    if (tmp[i] != 0x0 ) 
      yield_[i] = tmp[i]->Integral(0, 9999); 
    else 
      yield_[i]=0;
    tot_yield += yield_[i];
    std::cout << Form("%.3f",yield_[i])<<" " ;
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
    MCFMXsec_ = 0.004908;
    break;
    
  case (TVar::HZZ300):
    NLOXsec_ = 0.0300396;
    MCFMXsec_ = 0.003274;
    break;
    
  case (TVar::HZZ350):
    NLOXsec_ = 0.0286009;
    MCFMXsec_ = 0.002241;
    break;
    
  case (TVar::HZZ400):
    NLOXsec_ = 0.0220830;
    MCFMXsec_ = 0.001325;
    break;
    
  case (TVar::HZZ500):
    NLOXsec_ = 0.0089522;
    MCFMXsec_ = 0.000636;
    break;
    
  case (TVar::HZZ600):
    NLOXsec_ = 0.0035900;
    MCFMXsec_ = 0.000353;
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
