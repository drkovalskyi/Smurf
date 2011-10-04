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

void getProcess(int mH, TVar::Process & k);
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
  unsigned int njets_ = 0;
  
  ch->SetBranchAddress( "run",   &run_);
  ch->SetBranchAddress( "event", &event_);
  ch->SetBranchAddress( "type",  &type_);
  ch->SetBranchAddress("dXsec",  &dXsec_);
  ch->SetBranchAddress( "dilep", &dilep_); 
  ch->SetBranchAddress( "mt",  &mt_);
  ch->SetBranchAddress( "njets"      , &njets_     );     


  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = ch->GetEntries();
  if(nev == -1) Ntot = ch->GetEntries();
  if(nev > 0 && nev < Ntot) Ntot = nev;
  
  printf("Total number of events = %d\n", Ntot);
  
  
  TVar::Process k;
  getProcess(mH, k);

  float massCut = 99999.; // useless
  Proc *higgs = new Proc(k, lumi, massCut, meFDir, HZZANALYSIS);
  Proc *ww    =  new Proc(TVar::WW, lumi, massCut, meFDir, HZZANALYSIS);
  Proc *wz    =  new Proc(TVar::WZ, lumi, massCut, meFDir, HZZANALYSIS);
  Proc *zz    = new Proc(TVar::ZZ, lumi, massCut, meFDir, HZZANALYSIS);
  
  
if (verbosity >= TVar::DEBUG) {
    cout << "higgs->GetMCFMXsec() = " << higgs->GetMCFMXsec() << "\n";
    cout << "ww->GetMCFMXsec() = " << ww->GetMCFMXsec() << "\n";
    cout << "wz->GetMCFMXsec() = " << wz->GetMCFMXsec() << "\n";
    cout << "zz->GetMCFMXsec() = " << zz->GetMCFMXsec() << "\n";
  }
 
  for(int ievt=0;ievt<Ntot;ievt++){
    
    ch->GetEntry(ievt);  
    
    if (verbosity >= TVar::DEBUG) 
      cout << "\n ** START LR Construction, run = " << run_ << "; event = " << event_ << " for Signal " << TVar::SmurfProcessName(k) << "\n"; 
    
    // if ( njets_ > 0 ) continue;

    // ==== Construct HZZ LR
    
    // get the signal event probability
    double numer = dXsec_[k] / (higgs->GetMCFMXsec() * higgs->GetAcceptance(type_)); 
    if (verbosity >= TVar::DEBUG)
      cout<< "PHWW " << TVar::SmurfProcessName(k) << " "  <<numer<< "\t dXsec "<< dXsec_[k] << "\n";
    
    // get the background yields 
    double yield_bg = ww->GetYield(type_) + wz->GetYield(type_) +  zz->GetYield(type_);
    if (verbosity >= TVar::DEBUG)
      cout<<"bg_yield="<<yield_bg<<"\n";
    
    // add WW background to the denominator
    double denom  = numer;
    denom += dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWW = "<<  dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
    
    // add WZ to the denominator
    denom += dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_)) * wz->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWZ = "<<  dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_)) * wz->GetYield(type_)/yield_bg;
    
    // add ZZ to the denominator
    denom += dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_)) * zz->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PZZ = "<<  dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_)) * zz->GetYield(type_)/yield_bg;
    
    
    if(denom!=0)
      LR[k]=numer/denom;
    
    if (verbosity >= TVar::DEBUG) 
      cout<<"LR_HZZ["<<k<<"]= "<<LR[k]<<"\n";
    
    // ==== Construct WZ/ZZ LR
    // redefining the numerator and denumerator
    numer = dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_)); 
    if (verbosity >= TVar::DEBUG)
      cout<<" PZZ = "<<  dXsec_[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_));

    numer += dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_)); 
    if (verbosity >= TVar::DEBUG)
      cout<<" PWZ = "<<  dXsec_[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_));
    
    // get the background yields
    yield_bg = ww->GetYield(type_);
    if (verbosity >= TVar::DEBUG)
      cout<<"bg_yield="<<yield_bg<<"\n";
    
    // add WW background to the denominator
    denom  = numer;
    denom += dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWW = "<<  dXsec_[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;

    if(denom!=0)
      LR[TVar::ZZ]=numer/denom;
    if (verbosity >= TVar::DEBUG) 
      cout<<"LR_WZ/ZZ["<<TVar::ZZ<<"]= "<<LR[TVar::ZZ]<<"\n";
    
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


void getProcess(int mH, TVar::Process & k)
{
  switch (mH) {

  case (200):
    k = TVar::HZZ200;
    break;
  case (250):
    k = TVar::HZZ250;
    break;
  case (300):
    k = TVar::HZZ300;
    break;
  case (350):
    k = TVar::HZZ350;
    break;
  case (400):
    k = TVar::HZZ400;
    break;
  case (500):
    k = TVar::HZZ500;
    break;
  case (600):
    k = TVar::HZZ600;
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
