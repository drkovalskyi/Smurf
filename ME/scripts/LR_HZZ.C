// This code calculates the LR of HWW hypothesis
// LR = Psig / (Psig +  Sum_i{fbkg_i * Pbkg_i})
// Psig and Pbkg_i are normalized to a constant individually
// 

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "../TVar.hh"
#include "Proc.cc"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

using namespace std;

void getProcess(int mH, TVar::Process & k, float & metCut, float & dphiCut, float & mtCut);
float weightME(TString fileName, TString inputSmurfFDir, TString meFDir);

//###################
//# main function
//###################
void LR_HZZ(int mH, TString fileName, TString inputSmurfFDir, TString meFDir, int nev, float lumi, TVar::VerbosityLevel verbosity = TVar::INFO)  
{
  
  TString inputFileName = meFDir + fileName + ".root";
  TFile* fin = new TFile(inputFileName);
  std::cout << "Opening " << inputFileName << "\n";
  
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) { 
    std::cout << "Error:" << inputFileName << " is not properly filled, exitting\n";
    return; 
  }

  TString outFileName = inputFileName;
  outFileName.ReplaceAll(".root", "_ME.root");
    
  TFile *newfile = new TFile(outFileName,"recreate");
  std::cout << "creating " << outFileName << "...\n";

  

  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

  // Add additional event weight to account for the possible difference 
  // between the total events before the differential cross-section calcualtion
  // and the input. This is to protect the cases where some of the jobs 
  // fail the differential cross-section calcuation in the lpc grid
  double scaleME = 1.0;
  evt_tree->Branch("scaleME",  &scaleME, "scaleME/D");
  double scalefactor = weightME(fileName, inputSmurfFDir, meFDir);

  // The variable we want to recalculate
  double LR[60];
  ch->SetBranchAddress("LR", LR);

  // Initialize the branches to use to calculate LR
  double dXsec_ [60];
  int type_ = 0;
  unsigned int run_ = 0;
  unsigned int event_ = 0; 
  LorentzVector*  dilep_ = 0;
  float mt_ = 0.;
  float met_ = 0.;
  unsigned int njets_ = 0;
  LorentzVector*  jet1_ = 0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  float metPhi_ = 0.;

  ch->SetBranchAddress( "run",   &run_);
  ch->SetBranchAddress( "event", &event_);
  ch->SetBranchAddress( "type",  &type_);
  ch->SetBranchAddress("dXsec",  &dXsec_);
  ch->SetBranchAddress( "dilep", &dilep_); 
  ch->SetBranchAddress( "mt",  &mt_);
  ch->SetBranchAddress( "njets"      , &njets_     );     
  ch->SetBranchAddress( "met",  &met_);
  ch->SetBranchAddress( "jet1", &jet1_);
  ch->SetBranchAddress( "jet2", &jet2_);
  ch->SetBranchAddress( "jet3", &jet3_);
  ch->SetBranchAddress( "metPhi",  &metPhi_); 

  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = ch->GetEntries();
  if(nev == -1) Ntot = ch->GetEntries();
  if(nev > 0 && nev < Ntot) Ntot = nev;
  
  printf("Total number of events = %d\n", Ntot);
  
  float metCut (0.0), dphiCut(-1.0), mtCut (-1.0);
  
  TVar::Process k;
  getProcess(mH, k, metCut, dphiCut, mtCut);


  float massCut = 99999.; // useless
  Proc *higgs = new Proc(k, lumi, massCut, meFDir, HZZANALYSIS, metCut, dphiCut, mtCut);
  Proc *ww    =  new Proc(TVar::WW, lumi, massCut, meFDir, HZZANALYSIS, metCut, dphiCut, mtCut);
  Proc *wz    =  new Proc(TVar::WZ, lumi, massCut, meFDir, HZZANALYSIS, metCut, dphiCut, mtCut);
  Proc *zz    = new Proc(TVar::ZZ, lumi, massCut, meFDir, HZZANALYSIS, metCut, dphiCut, mtCut);
  
  
if (verbosity >= TVar::DEBUG) {
    cout << "higgs->GetMCFMXsec() = " << higgs->GetMCFMXsec() << "\n";
    cout << "ww->GetMCFMXsec() = " << ww->GetMCFMXsec() << "\n";
    cout << "wz->GetMCFMXsec() = " << wz->GetMCFMXsec() << "\n";
    cout << "zz->GetMCFMXsec() = " << zz->GetMCFMXsec() << "\n";
  }
 
  for(int ievt=0;ievt<Ntot;ievt++){
    
    ch->GetEntry(ievt);  
    
    int diltype = type_;
    
    if ( type_ == 1 || type_ == 2) diltype = 0;
    
    if (verbosity >= TVar::DEBUG) {
      cout << "\n ** START LR Construction, run = " << run_ << "; event = " << event_ << " for Signal " << TVar::SmurfProcessName(k) << "\n"; 
      cout << "Event kinematics type == " << type_ << "\n";
    }
    float dPhi1, dPhi2, dPhi3;
    jet1_->Pt() > 30 ? dPhi1 = acos(cos(metPhi_ - jet1_->Phi())) : dPhi1 = 999.9;
    jet2_->Pt() > 30 ? dPhi2 = acos(cos(metPhi_ - jet2_->Phi())) : dPhi2 = 999.9;
    jet3_->Pt() > 30 ? dPhi3 = acos(cos(metPhi_ - jet3_->Phi())) : dPhi3 = 999.9;
    
    // apply the dphi and dilepton pT cuts
    if ( TMath::Min(dPhi1, TMath::Min(dPhi2, dPhi3)) < dphiCut || dilep_->Pt() < 55) {
      LR[k] = -1.; 
      LR[TVar::ZZ] = -1;
    }
    
    else {
      
      // ==== Construct HZZ LR
      // get the signal event probability
      double numer = 0.0;
      if (verbosity >= TVar::DEBUG) {
	std::cout << " higgs->GetMCFMXsec() = " << higgs->GetMCFMXsec() << "\n";
	std::cout << " higgs->GetAcceptance(" << diltype << ") = " << higgs->GetAcceptance(diltype) << "\n";
      } 
	
      if ( higgs->GetMCFMXsec() * higgs->GetAcceptance(diltype) > 0.)
	numer = dXsec_[k] / (higgs->GetMCFMXsec() * higgs->GetAcceptance(diltype)); 
      if (verbosity >= TVar::DEBUG)
	cout<< "PHZZ " << TVar::SmurfProcessName(k) << " "  <<numer<< "\t dXsec "<< dXsec_[k] << "\n";
      
      // get the background yields 
      double yield_bg = ww->GetYield(diltype) + wz->GetYield(diltype) +  zz->GetYield(diltype);
      if (verbosity >= TVar::DEBUG)
	cout<<"bg_yield="<<yield_bg<<"\n";
      
      // add WW background to the denominator
      double denom  = numer;
      float bkgconst = 4.0;
      
      if ( ww->GetMCFMXsec() * ww->GetAcceptance(diltype) > 0.) 
	denom += bkgconst*dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(diltype)) * ww->GetYield(diltype)/yield_bg;
      if (verbosity >= TVar::DEBUG)
	cout<<" PWW = "<<  dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(diltype)) * ww->GetYield(diltype)/yield_bg;
      
      // add WZ to the denominator
      if ( wz->GetMCFMXsec() * wz->GetAcceptance(diltype) > 0.) 
	denom += bkgconst*dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(diltype)) * wz->GetYield(diltype)/yield_bg;
      if (verbosity >= TVar::DEBUG)
	cout<<" PWZ = "<<  dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(diltype)) * wz->GetYield(diltype)/yield_bg;
      
      // add ZZ to the denominator
      if ( zz->GetMCFMXsec() * zz->GetAcceptance(diltype) > 0.) 
	denom += bkgconst*dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(diltype)) * zz->GetYield(diltype)/yield_bg;
      if (verbosity >= TVar::DEBUG)
	cout<<" PZZ = "<<  dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(diltype)) * zz->GetYield(diltype)/yield_bg;
      
      
      if(denom!=0)
	LR[k]=numer/denom;
      
      if (verbosity >= TVar::DEBUG) 
	cout<<"LR_HZZ["<<k<<"]= "<<LR[k]<<"\n";
      
      // ==== Construct WZ/ZZ LR
      // initialize the numerator
      numer = 0.0;
      // redefining the numerator and denumerator
      if ( zz->GetMCFMXsec() * zz->GetAcceptance(diltype) > 0.) 
	numer = dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(diltype)); 
      if (verbosity >= TVar::DEBUG)
	cout<<" PZZ = "<<  dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(diltype));
      if ( wz->GetMCFMXsec() * wz->GetAcceptance(diltype) > 0.) 
	numer += dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(diltype)); 
      if (verbosity >= TVar::DEBUG)
	cout<<" PWZ = "<<  dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(diltype));
      
      // get the background yields
      yield_bg = ww->GetYield(diltype);
      if (verbosity >= TVar::DEBUG)
	cout<<"bg_yield="<<yield_bg<<"\n";
      
      // add WW background to the denominator
      denom  = numer;
      if ( ww->GetMCFMXsec() * ww->GetAcceptance(diltype) > 0.) 
	denom += dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(diltype)) * ww->GetYield(diltype)/yield_bg;
      if (verbosity >= TVar::DEBUG)
	cout<<" PWW = "<<  dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(diltype)) * ww->GetYield(diltype)/yield_bg;
      
      if(denom!=0)
	LR[TVar::ZZ]=numer/denom;
      if (verbosity >= TVar::DEBUG) 
	cout<<"LR_WZ/ZZ["<<TVar::ZZ<<"]= "<<LR[TVar::ZZ]<<"\n";
    }
    
    scaleME = scalefactor;
    evt_tree->Fill();
    
    
  }//nevent
  
  delete higgs;
  delete ww;
  delete wz;
  delete zz;
  
  cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
  delete fin;
}  


void getProcess(int mH, TVar::Process & k, float & metCut, float & dphiCut, float & mtCut)
{
  switch (mH) {

  case (200):
    k = TVar::HZZ200;
    metCut = 50.0;
    dphiCut = 0.5;
    break;
  case (250):
    k = TVar::HZZ250;
    metCut = 70;
    dphiCut = 0.5;
    break;
  case (300):
    k = TVar::HZZ300;
    metCut = 80;
    dphiCut = 0.5;
    break;
  case (350):
    k = TVar::HZZ350;
    metCut = 90;
    dphiCut = 0.5;
    break;
  case (400):
    k = TVar::HZZ400;
    metCut = 90;
    dphiCut = 0.5;
    break;
  case (500):
    k = TVar::HZZ500;
    metCut = 100.;
    dphiCut = 0.5;
    break;
  case (600):
    k = TVar::HZZ600;
    metCut = 120.;
    dphiCut = 0.5;
    break;
  default:
    break;
  }
}




float weightME(TString fileName, TString inputSmurfFDir, TString meFDir)
{
  TString inputSmurfFileName = inputSmurfFDir + fileName + ".root";
  TFile *fin = TFile::Open(inputSmurfFileName, "READ");
  if (fin == 0x0 ) {
    std::cout << "ERROR: weightME() file " << inputSmurfFileName + " doesn't exsit, weight set to 1....";
    return 1.0;
  }

 
  TTree *treein = (TTree*)fin->Get("tree");
  float Nin = treein->GetEntries();
  
  TString meFileName = meFDir + fileName + ".root";
  TFile *fout = TFile::Open(meFileName, "READ");
  if (fout == 0x0 ) {
    std::cout << "ERROR: weightME() file " << meFileName + " doesn't exsit, weight set to 0....";
    return 1.0;
  }
  TTree *treeout = (TTree*)fout->Get("tree");
  float Nout = treeout->GetEntries();
  fin->Close();
  fout->Close();
  return Nin/Nout;
}
