#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "commonFunction.h"
#include "../TVar.hh"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

using namespace std;

const int kNDilep=4;

float lumi=1000.;
float yield [kNProc][kNDilep];
float acceptance [kNProc][kNDilep];

float BR [kNProc][kNDilep] = {{ (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9,  (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)/3, (2.0+0.1736+0.1784)/6, (2.0+0.1736+0.1784)/6,   (1+0.1784)/3},
			      { (1.0+0.1736)/3, (2.0+0.1736+0.1784)/6, (2.0+0.1736+0.1784)/6,   (1+0.1784)/3},
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/6, (2.0+0.1736+0.1784)/6,   (1+0.1784)/3}, // Wpgamma placeholder
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/6, (2.0+0.1736+0.1784)/6,   (1+0.1784)/3}, // Wpgamma placeholder
			      { 2*0.2*(0.0363+0.0370*0.1736*0.1736), 0.0, 0.0, 2*0.2*(0.0363+0.0370*0.1784*0.1784) },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }, 
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1736)*(1.0+0.1784)/9, (1.0+0.1784)*(1.0+0.1784)/9 }
};

// The NLO xsec includes the W->l BR, where l includes e/mu/tau
// The values are obtained through http://ceballos.web.cern.ch/ceballos/hwwlnln/sigma_hww_2l2n_ec
float NLOXsec[kNProc] = {  4.5, 31314, 31314, 31314.0, 31314.0, 5.9, 0.249642, 0.452090, 0.641773, 0.770471, 0.866443, 0.782962, 0.659328, 0.486486, 0.408305, 0.358465, 0.321398, 0.290454, 0.243724, 0.175652};
// The MCFM total cross-sections 
float MCFMXsec[kNProc] = { 2.983, 18788.0, 12525.0,  11270, 11270.0, 4.3, 0.07533, 0.1451, 0.2165, 0.2671, 0.3034, 0.2870, 0.2465, 0.1849, 0.1570, 0.1380, 0.1232, 0.1107, 0.0906, 0.05758};

void CalculateAcceptance(){

    for (int i=0; i<kNProc; i++){
      for (int j=0; j<kNDilep; j++){
      float denominator = lumi* BR[i][j] * NLOXsec[i];
      if (denominator!=0) acceptance [i][j]= yield [i][j]/denominator; else acceptance [i][j]=1E+20;
      }
    }

}

void InitializeYields(TString inputDir, Float_t massCut, TVar::VerbosityLevel verbosity) 
{
  for (int i=0;i<kNProc;i++) {
    TFile *f = TFile::Open(TString(inputDir+name(i)+".root"));
    if(f==0x0) continue;
    if (verbosity >= TVar::INFO) cout << TString("../"+name(i)+".root") <<endl;
    
    TTree *tree = (TTree*) f->Get("tree");
    
    double tot_yield = 0;
    if (verbosity >= TVar::INFO) cout<<"Yields from file "<< name(i)<<".root: mm/me/em/ee:  ";    
    for (int j=0; j<kNDilep; j++){ 
      TH1F *tmp = new TH1F("tmp", "tmp", 20, 0, 20);
      if(j==0 || j==2)
	tree->Project("tmp", "dPhi", Form("scale1fb*(type==%i&&dilep.mass()<%f)", j, massCut));
      if(j==1 || j==3)
	tree->Project("tmp", "dPhi", Form("scale1fb*((type==%i&&lep2.pt()>15)&&dilep.mass()<%f)", j, massCut));
      
      if (tmp!= 0x0) yield[i][j] = tmp->Integral(0, 9999); 
      else yield[i][j]=0;
      tot_yield+=yield[i][j];
      if (verbosity >= TVar::INFO) cout << Form("%.3f",yield[i][j])<<" " ;
      delete tmp;
    }       
    if (verbosity >= TVar::INFO) cout << "; total yield = " << Form("%.3f", tot_yield) << "\n"; 
    f->Close();
  }
}

void CalculateLR(TString meFDir, TString fileName, TString outputDir, int maxevt, int k, float massCut, TVar::VerbosityLevel verbosity);

//###################
//# main function
//###################
void LR(TString smurfFDir, TString meFDir, TString fileName, TString outputDir, int nev, int k, float massCut, TVar::VerbosityLevel verbosity)  
{
  int maxevt=nev;
  InitializeYields(smurfFDir, massCut, verbosity);
  CalculateAcceptance();
  CalculateLR(meFDir, fileName, outputDir, maxevt, k, massCut, verbosity);
}

void CalculateLR(TString meFDir, TString fileName, TString outputDir, int maxevt, int k, float massCut, TVar::VerbosityLevel verbosity) {

  if (verbosity >= TVar::INFO) cout << "Input File: "<< fileName << "\n"; 
  TFile* fin = new TFile(meFDir + fileName);
  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll("_ME.root", Form("_LR_%s.root", name(k).c_str()));
			  
  TFile *newfile = new TFile(outFileName,"recreate");
  if (verbosity >= TVar::INFO) std::cout << "creating" << outFileName << "...\n";

  
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
  if(maxevt == -1) Ntot = ch->GetEntries();
  if(maxevt > 0 && maxevt < Ntot) Ntot = maxevt;
  
  printf("Total number of events = %d\n", Ntot);
  
  for(int ievt=0;ievt<Ntot;ievt++){
    
    ch->GetEntry(ievt);  

    if(dilep_->mass() > massCut) continue;
    
    if (verbosity >= TVar::INFO && (ievt % 100 == 0)) 
      std::cout << "Doing Event: " << ievt << std::endl;
    
    // 
    // Begin HWW likelihood 
    // LR = Psig / (Psig +  Sum_i{fbkg_i * Pbkg_i})
    // Psig and Pbkg_i are normalized to a constant individually
    // 
    if (verbosity > TVar::INFO)  cout << "\n ** START LR Construction, run = " << run_ << "; event = " << event_ << " for Signal " << name(k) << "\n"; 

    // get the signal event probability
    double numer = dXsecList[k] / (MCFMXsec[k]*acceptance[k][type_]);
    if (verbosity > TVar::INFO) cout<< "PHWW " << name(k) << " "  <<numer<< "\t dXsec "<< dXsecList[k] << "\n";
    
    // get the background yields
    double yield_bg=yield[proc_WW][type_]+yield[proc_Wpj][type_];
    if (verbosity > TVar::INFO) cout<<"bg_yield="<<yield_bg<<"\n";

    double denom  = numer;
    // add WW background to the denominator  
    denom +=  dXsecList[proc_WW ] /(MCFMXsec[proc_WW ]*acceptance[proc_WW ][type_]) * yield[proc_WW][type_]/yield_bg;
    if (verbosity > TVar::INFO) cout<<" PWW= "<< dXsecList[proc_WW ] / (MCFMXsec[proc_WW ]*acceptance[proc_WW ][type_]);
    
    // add Wjet background combining W+jet and W-jet, note that
    // 1. Same acceptance is used for Wpj and Wmj 
    // 2. yield[proc_Wpj] = yield[proc_Wmj], which are filled as the combined w+jet yield
    denom +=  (dXsecList[proc_Wpj ] + dXsecList[proc_Wmj]) /( ( MCFMXsec[proc_Wpj ] + MCFMXsec[proc_Wmj]) * acceptance[proc_Wpj ][type_]) * yield[proc_Wpj][type_]/yield_bg;
    if (verbosity > TVar::INFO) cout<<" PWj= "<< (dXsecList[proc_Wpj ] + dXsecList[proc_Wmj]) /( ( MCFMXsec[proc_Wpj ] + MCFMXsec[proc_Wmj]) * acceptance[proc_Wpj ][type_])  * yield[proc_Wpj][type_]/yield_bg << "\n";
    
    
    if(denom!=0)
      LR[k]=numer/denom;
    if (verbosity > TVar::INFO) cout<<"LR_HWW["<<k<<"]= "<<LR[k]<<"\n";
    
    evt_tree->Fill();
    
  }//nevent
  
  if (verbosity >= TVar::INFO)     cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
  delete fin;
}  
