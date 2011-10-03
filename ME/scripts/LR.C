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
float weightME(TString fileName, TString inputSmurfFDir, TString meFDir);

//###################
//# main function
//###################
void LR(int mH, TString fileName, TString inputSmurfFDir, TString meFDir, int nev, float lumi, TVar::VerbosityLevel verbosity = TVar::INFO)  
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
  
  ch->SetBranchAddress( "run",   &run_);
  ch->SetBranchAddress( "event", &event_);
  ch->SetBranchAddress( "type",  &type_);
  ch->SetBranchAddress("dXsec",  &dXsec_);
  ch->SetBranchAddress( "dilep", &dilep_); 
  ch->SetBranchAddress( "mt",  &mt_);
  
  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = ch->GetEntries();
  if(nev == -1) Ntot = ch->GetEntries();
  if(nev > 0 && nev < Ntot) Ntot = nev;
  
  printf("Total number of events = %d\n", Ntot);
  
  
  TVar::Process k;
  float massCut;
  getProcess(mH, k, massCut);


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
    
    if (verbosity >= TVar::DEBUG) 
      cout << "\n ** START LR Construction, run = " << run_ << "; event = " << event_ << " for Signal " << TVar::SmurfProcessName(k) << "\n"; 
    
    // calculate the LR only for the events that pass the pre-selection
    if ( dilep_->mass() > massCut || mt_ < 80.)  
      LR[k] = -1.;
    else {
      // get the signal event probability
      double numer = dXsec_[k] / (higgs->GetMCFMXsec() * higgs->GetAcceptance(type_)); 
      if (verbosity >= TVar::DEBUG)
	cout<< "PHWW " << TVar::SmurfProcessName(k) << " "  <<numer<< "\t dXsec "<< dXsec_[k] << "\n";
      
      // get the background yields
      double yield_bg = ww->GetYield(type_) + wpj->GetYield(type_);
      if (verbosity >= TVar::DEBUG)
	cout<<"bg_yield="<<yield_bg<<"\n";
      
      // add WW background to the denominator
      double denom  = numer;
      denom += dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
      if (verbosity >= TVar::DEBUG)
	cout<<" PWW = "<<  dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
      
      // add Wjet background combining W+jet and W-jet assume the same acceptance for Wpj and Wmj 
      if ( mH < 150) {
	denom +=  15*(dXsec_[TVar::Wp_1jet ] + dXsec_[TVar::Wm_1jet]) /( ( wpj->GetMCFMXsec() + wmj->GetMCFMXsec()) * wpj->GetAcceptance(type_)) * wpj->GetYield(type_)/yield_bg;
      }
      else
	denom +=  (dXsec_[TVar::Wp_1jet ] + dXsec_[TVar::Wm_1jet]) /( ( wpj->GetMCFMXsec() + wmj->GetMCFMXsec()) * wpj->GetAcceptance(type_)) * wpj->GetYield(type_)/yield_bg;
      if (verbosity >= TVar::DEBUG)
	cout<<" PWj = " << (dXsec_[TVar::Wp_1jet ] + dXsec_[TVar::Wm_1jet]) /( ( wpj->GetMCFMXsec() + wmj->GetMCFMXsec()) * wpj->GetAcceptance(type_)) * wpj->GetYield(type_)/yield_bg << "\n";
      
      if(denom!=0)
	LR[k]=numer/denom;
    }
    
    if (verbosity >= TVar::DEBUG) 
      cout<<"LR_HWW["<<k<<"]= "<<LR[k]<<"\n";
    
    scaleME = scalefactor;
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
  case (115):
    k = TVar::HWW115;
    massCut = 70.0;
    break;
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
  case (350):
    k = TVar::HWW350;
    massCut = 350.0;
    break;
  case (400):
    k = TVar::HWW400;
    massCut = 400.0;
    break;
  case (450):
    k = TVar::HWW450;
    massCut = 450.0;
    break;
  case (500):
    k = TVar::HWW500;
    massCut = 500.0;
    break;
  case (550):
    k = TVar::HWW550;
    massCut = 550.0;
    break;
  case (600):
    k = TVar::HWW600;
    massCut = 600.0;
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
