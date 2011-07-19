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
float weightME(TString fileName, TString meFDir, TString mebdtFDir, TString bdtprefix);

//###################
//# main function
//###################
void LR_HZZ(int mH, TString bdtprefix, TString fileName, TString meFDir, TString mebdtFDir,  int nev, TVar::VerbosityLevel verbosity = TVar::INFO)  
{
  
  TVar::Process k;
  float massCut;
  getProcess(mH, k, massCut);
  
  cout << "Input File: "<< fileName << "\n"; 
  TFile* fin = new TFile(mebdtFDir + fileName + "_ME_merged.root");
  cout << mebdtFDir + bdtprefix + "_" + fileName + "_ME.root"<<endl;
  TString outFileName = mebdtFDir + fileName + "_ME.root" ;
  outFileName.ReplaceAll("_ME.root", Form("_LR_%s.root", TVar::SmurfProcessName(k).Data()));
    
  TFile *newfile = new TFile(outFileName,"recreate");
  std::cout << "creating " << outFileName << "...\n";

  
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");


  // Add additional event weight to account for the possible difference 
  // between the total events before the differential cross-section calcualtion
  // and the input. This is to protect the cases where some of the jobs 
  // fail the differential cross-section calcuation in the lpc grid
  double scaleME = 1.0;
  evt_tree->Branch("scaleME",  &scaleME, "scaleME/D");
  double scalefactor = weightME(fileName, meFDir, mebdtFDir, bdtprefix);

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
  
  float lumi = 1092.0; 
  Proc *higgs = new Proc(k, lumi, massCut, meFDir);
  Proc *ww    =  new Proc(TVar::WW, lumi, massCut, meFDir);
  Proc *wz    =  new Proc(TVar::WZ, lumi, massCut, meFDir);
  Proc *zz    = new Proc(TVar::ZZ, lumi, massCut, meFDir);
  
if (verbosity >= TVar::DEBUG) {
    cout << "higgs->GetMCFMXsec() = " << higgs->GetMCFMXsec() << "\n";
    cout << "ww->GetMCFMXsec() = " << ww->GetMCFMXsec() << "\n";
    cout << "wz->GetMCFMXsec() = " << wz->GetMCFMXsec() << "\n";
    cout << "zz->GetMCFMXsec() = " << zz->GetMCFMXsec() << "\n";
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
    double yield_bg = ww->GetYield(type_) + wz->GetYield(type_) +  zz->GetYield(type_);
    if (verbosity >= TVar::DEBUG)
      cout<<"bg_yield="<<yield_bg<<"\n";
    
    // add WW background to the denominator
    double denom  = numer;
    denom += dXsecList[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWW = "<<  dXsecList[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;

    denom += dXsecList[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_)) * wz->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWZ = "<<  dXsecList[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_)) * wz->GetYield(type_)/yield_bg;

    denom += dXsecList[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_)) * zz->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PZZ = "<<  dXsecList[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_)) * zz->GetYield(type_)/yield_bg;

    if(denom!=0)
      LR[k]=numer/denom;
    if (verbosity >= TVar::DEBUG) 
      cout<<"LR_HWW["<<k<<"]= "<<LR[k]<<"\n";



    // LR_ZZ
     // get the signal event probability
    numer = dXsecList[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_)); 
    if (verbosity >= TVar::DEBUG)
      cout<<" PZZ = "<<  dXsecList[TVar::ZZ] / (zz->GetMCFMXsec() * zz->GetAcceptance(type_));

    numer += dXsecList[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_)); 
    if (verbosity >= TVar::DEBUG)
      cout<<" PWZ = "<<  dXsecList[TVar::WZ] / (wz->GetMCFMXsec() * wz->GetAcceptance(type_));
    
    
    // get the background yields
    yield_bg = ww->GetYield(type_);
    if (verbosity >= TVar::DEBUG)
      cout<<"bg_yield="<<yield_bg<<"\n";
    
    // add WW background to the denominator
    denom  = numer;

    denom += dXsecList[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;
    if (verbosity >= TVar::DEBUG)
      cout<<" PWW = "<<  dXsecList[TVar::WW] / (ww->GetMCFMXsec() * ww->GetAcceptance(type_)) * ww->GetYield(type_)/yield_bg;

    if(denom!=0)
      LR[TVar::ZZ]=numer/denom;
    if (verbosity >= TVar::DEBUG) 
      cout<<"LR_HWW["<<TVar::ZZ<<"]= "<<LR[TVar::ZZ]<<"\n";

    scaleME = scalefactor;
    
    if (fileName=="zz") scaleME*=7.41/5.9;
    //    if ((fileName=="dyee")||(fileName=="dymm")) scaleME*=16.4/0.81;
    if ((fileName=="qqww")||(fileName=="ggww")||(fileName=="ttbar")||(fileName=="tw")) scaleME*=1.25;

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

void getProcess(int mH, TVar::Process & k, float& massCut)
{
  switch (mH) {
  case (200):
    k = TVar::HZZ200;
    massCut = 1000.0;
    break;
  case (250):
    k = TVar::HZZ250;
    massCut = 1000.0;
    break;
  case (300):
    k = TVar::HZZ300;
    massCut = 1000.0;
    break;
  case (400):
    k = TVar::HZZ400;
    massCut = 1000.0;
    break;
  default:
    break;
  }
}

float weightME(TString fileName, TString meFDir, TString mebdtFDir, TString bdtprefix)
{
  TFile *fin = TFile::Open(meFDir + fileName + ".root");
  if (fin == 0x0 ) {
    std::cout << "ERROR: weightME() file " << meFDir + fileName + ".root doesn't exsit, weight set to 1....";
    return 1.0;
  }

  TTree *treein = (TTree*)fin->Get("tree");
  float Nin = treein->GetEntries();
  
  TFile *fout = TFile::Open(mebdtFDir + bdtprefix + "_" + fileName + "_ME.root");
  if (fout == 0x0 ) {
    std::cout << "ERROR: weightME() file " << meFDir + fileName + ".root doesn't exsit, weight set to 0....";
    return 1.0;
  }
  TTree *treeout = (TTree*)fout->Get("tree");
  float Nout = treeout->GetEntries();
  fin->Close();
  fout->Close();
  return Nin/Nout;
}
