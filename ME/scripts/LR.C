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

void getProcess(int mH, TVar::Process & k, float & massCut);
//###################
//# main function
//###################
void LR(int mH, TString bdtprefix, TString fileName, TString meFDir, TString mebdtFDir,  int nev, TVar::VerbosityLevel verbosity = TVar::INFO)  
{
  
  TVar::Process k;
  float massCut;
  getProcess(mH, k, massCut);
  
  cout << "Input File: "<< fileName << "\n"; 
  TFile* fin = new TFile(mebdtFDir + bdtprefix + "_" + fileName);
  cout << mebdtFDir + bdtprefix + "_" + fileName <<endl;
  TString outFileName = mebdtFDir + fileName;
  outFileName.ReplaceAll("_ME.root", Form("_LR_%s.root", TVar::SmurfProcessName(k).Data()));
    
  TFile *newfile = new TFile(outFileName,"recreate");
  std::cout << "creating" << outFileName << "...\n";

  
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

  // The variable we want to recalculate
  double LR[60];
  ch->SetBranchAddress("LR", LR);

  // Initialize the branches to use to calculate LR
  double dXsecList [60];
  int type_ = 0;
  unsigned int run_ = 0;
  unsigned int event_ = 0; 
  LorentzVector*  lep2_ = 0;
  LorentzVector*  dilep_ = 0;

  ch->SetBranchAddress( "run",   &run_);
  ch->SetBranchAddress( "event", &event_);
  ch->SetBranchAddress( "type",  &type_);
  ch->SetBranchAddress("dXsec",  &dXsecList);
  ch->SetBranchAddress( "lep2",  &lep2_      ); 
  ch->SetBranchAddress( "dilep", &dilep_      ); 

  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = ch->GetEntries();
  if(nev == -1) Ntot = ch->GetEntries();
  if(nev > 0 && nev < Ntot) Ntot = nev;
  
  printf("Total number of events = %d\n", Ntot);
  
  float lumi = 1000.0; 
  Proc *higgs = new Proc(k, lumi, massCut, meFDir);
  Proc *ww    =  new Proc(TVar::WW, lumi, massCut, meFDir);
  Proc *wpj    = new Proc(TVar::Wp_1jet, lumi, massCut, meFDir);
  Proc *wmj    = new Proc(TVar::Wm_1jet, lumi, massCut, meFDir);
  
if (verbosity >= TVar::DEBUG) {
    cout << "higgs->GetMCFMXsec() = " << higgs->GetMCFMXsec() << "\n";
    cout << "ww->GetMCFMXsec() = " << ww->GetMCFMXsec() << "\n";
    cout << "wpj->GetMCFMXsec() = " << wpj->GetMCFMXsec() << "\n";
    cout << "wmj->GetMCFMXsec() = " << wmj->GetMCFMXsec() << "\n";
  }
 
  for(int ievt=0;ievt<Ntot;ievt++){

    ch->GetEntry(ievt);  

    if(dilep_->mass() > massCut) continue;

    if (verbosity >= TVar::DEBUG) 
      cout << "\n ** START LR Construction, run = " << run_ << "; event = " << event_ << " for Signal " << TVar::SmurfProcessName(k) << "\n"; 
    
    // get the signal event probability
    double numer = dXsecList[k] / (higgs->GetMCFMXsec() * higgs->GetAcceptance(type_)); 
    if (verbosity >= TVar::DEBUG)
      cout<< "PHWW " << TVar::SmurfProcessName(k) << " "  <<numer<< "\t dXsec "<< dXsecList[k] << "\n";
    
    // get the background yields
    double yield_bg = ww->GetYield(type_) + wpj->GetYield(type_);
    if (verbosity >= TVar::DEBUG)
      cout<<"bg_yield="<<yield_bg<<"\n";
    
    // add WW background to the denominator
    double denom  = numer;
    denom += dXsecList[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWW = "<<  dXsecList[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;

    // add Wjet background combining W+jet and W-jet assume the same acceptance for Wpj and Wmj 
    denom +=  (dXsecList[TVar::Wp_1jet ] + dXsecList[TVar::Wm_1jet]) /( ( wpj->GetMCFMXsec() + wmj->GetMCFMXsec()) * wpj->GetAcceptance(type_)) * wpj->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWj = " << (dXsecList[TVar::Wp_1jet ] + dXsecList[TVar::Wm_1jet]) /( ( wpj->GetMCFMXsec() + wmj->GetMCFMXsec()) * wpj->GetAcceptance(type_)) * wpj->GetYield(type_)/yield_bg << "\n";
    
    if(denom!=0)
      LR[k]=numer/denom;
    if (verbosity >= TVar::DEBUG) 
      cout<<"LR_HWW["<<k<<"]= "<<LR[k]<<"\n";
    
    evt_tree->Fill();
    
    
  }//nevent
  
  delete higgs;
  delete ww;
  delete wpj;
  delete wmj;
  
  cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
  delete fin;
}  

void getProcess(int mH, TVar::Process & k, float& massCut)
{
  switch (mH) {
  case (120):
    k = TVar::HWW120;
    massCut = 70.0;
    break;
  case (130):
    k = TVar::HWW130;
    massCut = 80.0;
    break;
  case (140):
    k = TVar::HWW140;
    massCut = 90.0;
    break;
  case (150):
    k = TVar::HWW150;
    massCut = 100.0;
    break;
  case (160):
    k = TVar::HWW160;
    massCut = 100.0;
    break;
  case (170):
    k = TVar::HWW170;
    massCut = 100.0;
    break;
  case (180):
    k = TVar::HWW180;
    massCut = 110.0;
    break;
  case (190):
    k = TVar::HWW190;
    massCut = 120.0;
    break;
  case (200):
    k = TVar::HWW200;
    massCut = 130.0;
    break;
  case (210):
    k = TVar::HWW210;
    massCut = 140.0;
    break;
  case (220):
    k = TVar::HWW220;
    massCut = 150.0;
    break;
  case (230):
    k = TVar::HWW230;
    massCut = 230.0;
    break;
  case (250):
    k = TVar::HWW250;
    massCut = 250.0;
    break;
  case (300):
    k = TVar::HWW300;
    massCut = 300.0;
    break;
  default:
    break;
  }
  

}
