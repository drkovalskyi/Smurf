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
  else if (proc_ == TVar::ZZ || ( proc_ >= TVar::HZZ200 && proc_ <= TVar::HZZ600) || proc_==TVar::WZ || proc_==TVar::Z_2l) {
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
  float pmet_ = 0.;
  float pTrackMet_ = 0.;
  LorentzVector*  lep2_ = 0;
  float dymva_ = 0.;
  unsigned int njets_ = 0;
  float dPhiDiLepJet1_ = 0.;

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
  tree->SetBranchAddress( "dymva",  &dymva_);
  tree->SetBranchAddress( "njets",  &njets_);
  tree->SetBranchAddress( "dPhiDiLepJet1",  &dPhiDiLepJet1_);

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
      if ( mt_ > mtCut_ ) passpresel = false;
      if ( dilep_->Pt() < 45) passpresel = false;
      if ( type_ == 0 || type_ == 3) {
	float metValue_ = TMath::Min(pmet_,pTrackMet_);
	if ( proc_ >= TVar::HWW110  && proc_ <= TVar::HWW140 ) {
	  if ( njets_ == 0 && dymva_ > 0.6 ) passpresel = false;
	  if ( njets_ == 1 && dymva_ > 0.3 ) passpresel = false;
	} else {
	  if (metValue_ < 45. ) passpresel = false;
	  if ( jet1_->Pt() > 15 && dPhiDiLepJet1_ > (165./180.*TMath::Pi()) ) 
	    passpresel = false;
	}
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

//
// Note:: The MCFM xsections are updated with the MCFM6.2
// The 7 TeV cross-sections are not done yet for mH>600.
//

void Proc::initTotXsec() {
  
  bool EightTeV = false;

  switch (proc_) {
    
  case (TVar::WW): 
    NLOXsec_ = 4.5;
    MCFMXsec_ = 0.333;
    if ( EightTeV ) {
      NLOXsec_ = 4.5;
      MCFMXsec_ = 0.4;
    }
    break;

  case (TVar::Wp_1jet):
    NLOXsec_ = 31314.0;
    MCFMXsec_ = 18788.0/3.0;
    if ( EightTeV ) {
      NLOXsec_ = 4.5;
      MCFMXsec_ = 21222./3.0;
    }
    break;

  case (TVar::Wm_1jet):
    NLOXsec_ = 31314.0;
    MCFMXsec_ =  12525.0/3.0;
    if ( EightTeV ) {
      NLOXsec_ = 31314.0;
      MCFMXsec_ =  14995.0/3.0;
    }
    break;

  case (TVar::HWW110):
    NLOXsec_ = 0.165009;
    MCFMXsec_ = 0.003447;
    if ( EightTeV ) {
      NLOXsec_ = 0.165009;
      MCFMXsec_ = 0.00479;
    }
    
    break;

  case (TVar::HWW115):
    NLOXsec_ = 0.165009;
    MCFMXsec_ = 0.005754;
    if ( EightTeV ) {
      NLOXsec_ = 0.165009;
      MCFMXsec_ = 0.007462;
    }
    break;

  case (TVar::HWW120):
    NLOXsec_ = 0.249642;
    MCFMXsec_ = 0.008743;
    if ( EightTeV ) {
      NLOXsec_ = 0.249642;
      MCFMXsec_ = 0.011372;
    }
    break;

  case (TVar::HWW125):
    NLOXsec_ = 0.249642;
    MCFMXsec_ = 0.012234;
    if ( EightTeV ) {
      NLOXsec_ = 0.249642;
      MCFMXsec_ = 0.015957;
    }
    break;
    
  case (TVar::HWW130):
    NLOXsec_ = 0.452090;
    MCFMXsec_ = 0.0159;
    if ( EightTeV ) {
      NLOXsec_ = 0.452090;
      MCFMXsec_ = 0.0208;
    }
    break;

 case (TVar::HWW135):
   NLOXsec_ = 0.452090;
   MCFMXsec_ = 0.01936;
    if ( EightTeV ) {
      NLOXsec_ = 0.452090;
      MCFMXsec_ = 0.0254;
    }
    break;

  case (TVar::HWW140):
    NLOXsec_ = 0.641773;
    MCFMXsec_ = 0.02234;
    if ( EightTeV ) {
      NLOXsec_ = 0.641773;
      MCFMXsec_ = 0.0294;
    }
    break;

  case (TVar::HWW145):
    NLOXsec_ = 0.641773;
    MCFMXsec_ = 0.025;
    if ( EightTeV ) {
      NLOXsec_ = 0.641773;
      MCFMXsec_ = 0.033;
    }
    break;

  case (TVar::HWW150):
    NLOXsec_ = 0.770471;
    MCFMXsec_ = 0.02674;
    if ( EightTeV ) {
      NLOXsec_ = 0.770471;
      MCFMXsec_ = 0.0354;
    }
    break;

  case (TVar::HWW155):
    NLOXsec_ = 0.770471;
    MCFMXsec_ = 0.0281;
    if ( EightTeV ) {
      NLOXsec_ = 0.770471;
      MCFMXsec_ = 0.0373;
    }
    break;


  case (TVar::HWW160):
    NLOXsec_ = 0.866443;
    MCFMXsec_ = 0.0293;
    if ( EightTeV ) {
      NLOXsec_ = 0.782962;
      MCFMXsec_ = 0.039;
    }
    break;

  case (TVar::HWW170):
    NLOXsec_ = 0.782962;
    MCFMXsec_ = 0.0269;
    if ( EightTeV ) {
      NLOXsec_ = 0.782962;
      MCFMXsec_ = 0.036;
    }
    break;

  case (TVar::HWW180):
    NLOXsec_ = 0.659328;
    MCFMXsec_ = 0.0225;
    if ( EightTeV ) {
      NLOXsec_ = 0.659328;
      MCFMXsec_ = 0.030;
    }
    break;

  case (TVar::HWW190):
    NLOXsec_ = 0.486486;
    MCFMXsec_ = 0.0164;
    if ( EightTeV ) {
      NLOXsec_ = 0.486486;
      MCFMXsec_ = 0.0222;
    }
    break;

  case (TVar::HWW200):
    NLOXsec_ = 0.408305;
    MCFMXsec_ = 0.0136;
    if ( EightTeV ) {
      NLOXsec_ = 0.408305;
      MCFMXsec_ = 0.0185;
    }
    break;

  case (TVar::HWW250):
    NLOXsec_ = 0.243724;
    MCFMXsec_ = 0.007;
    if ( EightTeV ) {
      NLOXsec_ = 0.243724;
      MCFMXsec_ = 0.0097;
    }
    break;
    
  case (TVar::HWW300):
    NLOXsec_ = 0.175652;
    MCFMXsec_ = 0.004024;
    if ( EightTeV ) {
      NLOXsec_ = 0.175652;
      MCFMXsec_ = 0.0057;
    }
    break;

  case (TVar::HWW350):
    NLOXsec_ = 0.160052;
    MCFMXsec_ = 0.002439;
    if ( EightTeV ) {
      NLOXsec_ = 0.160052;
      MCFMXsec_ = 0.00353;
    }
    break;

  case (TVar::HWW400):
    NLOXsec_ = 0.124330;
    MCFMXsec_ = 0.0013;
    if ( EightTeV ) {
      NLOXsec_ = 0.124330;
      MCFMXsec_ = 0.00192;
    }
    break;

  case (TVar::HWW450):
    NLOXsec_ = 0.078433;
    MCFMXsec_ = 0.0008;
    if ( EightTeV ) {
      NLOXsec_ = 0.124330;
      MCFMXsec_ = 0.0012;
    }
    break;

  case (TVar::HWW500):
    NLOXsec_ = 0.048702;
    MCFMXsec_ = 0.00052;
    if ( EightTeV ) {
      NLOXsec_ = 0.048702;
      MCFMXsec_ = 0.0008;
    }
    break;

  case (TVar::HWW550):
    NLOXsec_ = 0.030364;
    MCFMXsec_ = 0.000356;
    if ( EightTeV ) {
      NLOXsec_ = 0.030364;
      MCFMXsec_ = 0.000559;
    }
    break;

  case (TVar::HWW600):
    NLOXsec_ = 0.019184;
    MCFMXsec_ = 0.00025;
    if ( EightTeV ) {
      NLOXsec_ = 0.019184;
      MCFMXsec_ = 0.0004;
    }
    break;

  case (TVar::HWW700):
    NLOXsec_ = 0.019184;
    MCFMXsec_ = 0.000126;
    if ( EightTeV ) {
      NLOXsec_ = 0.019184;
      MCFMXsec_ = 0.0002;
    }
    break;

  case (TVar::HWW800):
    NLOXsec_ = 0.019184;
    MCFMXsec_ = 0.000065;
    if ( EightTeV ) {
      NLOXsec_ = 0.019184;
      MCFMXsec_ = 0.00011;
    }
    break;

  case (TVar::HWW900):
    NLOXsec_ = 0.019184;
    MCFMXsec_ = 0.000034;
    if ( EightTeV ) {
      NLOXsec_ = 0.019184;
      MCFMXsec_ = 0.00006;
    }
    break;

  case (TVar::HWW1000):
    NLOXsec_ = 0.019184;
    MCFMXsec_ = 0.000018;
    if ( EightTeV ) {
      NLOXsec_ = 0.019184;
      MCFMXsec_ = 0.000033;
    }
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

  case (TVar::Z_2l):
    NLOXsec_ = 2.0*1666.0;
    MCFMXsec_ = 6440.0;
    break;    

  default:
    NLOXsec_ = 0.0;
    MCFMXsec_ = 0.0;
    break;
  }
    
}
