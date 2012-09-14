 /*
   runME_test(smurfFDir, "ww.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, 0, TVar::INFO);
   0 - location of the smurf ntuples
   1 - the sample to run over
   2 - the location of the output me files
   3 - random number seed
   4 - smear level (0: no smearing, 1: include kx/ky, 2: include resolution which is not implemented so far)
   5 - number of integration steps
   6 - thresold for the dXsec/dXsecErr to print out information
   7 - maximum events
   8 - start event number
   9 - level of output
 */
 
#include "TVar.hh"
#include "TEvtProb.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double ERRORthreshold=1.0;
using namespace std;

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity);

//###################
//# main function
//###################
void runME_HZZ4l(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity=TVar::INFO){
  
  if (verbosity >= TVar::INFO) cout <<"=== Calculating differential cross-section ==========" <<endl;  
  xseccalc(inputDir, fileName, outputDir, maxevt, verbosity); 
 
}

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity){

  if (verbosity >= TVar::INFO) cout << "Input File: " << fileName << " \n";
  
  TFile* fin = new TFile(inputDir+fileName);
  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll(".root","_ME.root");
  cout << outFileName <<endl;
  TFile *newfile = new TFile(outFileName,"recreate");

  TTree* ch=(TTree*)fin->Get("passedEvents"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // Declare the matrix element related variables to be added to the existing ntuples
  double dXsecList;
  double dXsecErrList;
  
  evt_tree->Branch("dXsec"      ,&dXsecList           ,"dXsec/D");
  evt_tree->Branch("dXsecErr"   ,&dXsecErrList        ,"dXsecErr/D");
  
  // Initialize the branches to use to calculate the differential cross-sections
  Double_t EL1_ = 0.;
  Double_t pXL1_ = 0.;
  Double_t pYL1_ = 0.;
  Double_t pZL1_ = 0.;

  Double_t EL2_ = 0.;
  Double_t pXL2_ = 0.;
  Double_t pYL2_ = 0.;
  Double_t pZL2_ = 0.;
  
  Double_t EL3_ = 0.;
  Double_t pXL3_ = 0.;
  Double_t pYL3_ = 0.;
  Double_t pZL3_ = 0.;

  Double_t EL4_ = 0.;
  Double_t pXL4_ = 0.;
  Double_t pYL4_ = 0.;
  Double_t pZL4_ = 0.;
  
  ch->SetBranchAddress( "EL1"       , &EL1_      );   
  ch->SetBranchAddress( "pXL1"      , &pXL1_     );   
  ch->SetBranchAddress( "pYL1"      , &pYL1_     );   
  ch->SetBranchAddress( "pZL1"      , &pZL1_     );   

  ch->SetBranchAddress( "EL2"       , &EL2_      );   
  ch->SetBranchAddress( "pXL2"      , &pXL2_     );   
  ch->SetBranchAddress( "pYL2"      , &pYL2_     );   
  ch->SetBranchAddress( "pZL2"      , &pZL2_      );   

  ch->SetBranchAddress( "EL3"       , &EL3_      );   
  ch->SetBranchAddress( "pXL3"      , &pXL3_     );   
  ch->SetBranchAddress( "pYL3"      , &pYL3_     );   
  ch->SetBranchAddress( "pZL3"      , &pZL3_      );   

  ch->SetBranchAddress( "EL4"       , &EL4_      );   
  ch->SetBranchAddress( "pXL4"      , &pXL4_     );   
  ch->SetBranchAddress( "pYL4"      , &pYL4_     );   
  ch->SetBranchAddress( "pZL4"      , &pZL4_      );   
  
  // The reconstructed event level information  
  // Create the instance of TEvtProb to calculate the differential cross-section
  // and set the calculation global variables common for all processes

  TEvtProb Xcal2;  
  Xcal2.SetApplyFake(0);
  Xcal2.SetSmearLevel(2);
  Xcal2.SetSeed(1);
  Xcal2.SetMatrixElement(TVar::MCFM);
  
  if (verbosity >= TVar::INFO) cout << "Integration Seed= "<< Xcal2._seed << " SmearLevel= " << Xcal2._smearLevel << " Ncalls = " << Xcal2._ncalls <<  endl;  


  hzz4l_event_type hzz4l_event;
  
  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = maxevt; 
  
  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);
  
  for(int ievt = 0; ievt < Ntot; ievt++){
    if (verbosity >= TVar::INFO && (ievt % 5 == 0)) 
      std::cout << "Doing Event: " << ievt << std::endl;
    
    dXsecList = 0.;
    dXsecErrList = 0.;
    
    ch->GetEntry(ievt);           
    
    hzz4l_event.p[0].SetXYZM(pXL1_, pYL1_, pZL1_, 0.);
    hzz4l_event.p[1].SetXYZM(pXL2_, pYL2_, pZL2_, 0.);
    hzz4l_event.p[2].SetXYZM(pXL3_, pYL3_, pZL3_, 0.);
    hzz4l_event.p[3].SetXYZM(pXL4_, pYL4_, pZL4_, 0.);
    
    if (verbosity >= TVar::DEBUG) {
      cout << "\n=========================================================\n";
      cout << "Entry: " << ievt << "\n";
      cout << "Input: ==================================================" <<endl;
      printf("lep1 p3 = (%4.4f, %4.4f, %4.4f)  lep2 p3 = (%4.4f, %4.4f, %4.4f)\n",
	     pXL1_, pYL1_, pZL1_, pXL2_, pYL2_, pZL2_); 
      printf("lep3 p3 = (%4.4f, %4.4f, %4.4f)  lep4 p3 = (%4.4f, %4.4f, %4.4f)\n",
	     pXL3_, pYL3_, pZL3_, pXL4_, pYL4_, pZL4_); 
      cout << "=========================================================\n";
    } 
    // finish loading event information
    
    // ==== Begin the differential cross-section calculation
    TVar::Process Proc = TVar::ZZ_4l;
    double Xsec = 0.;
    double XsecErr = 0.;
    Xcal2.SetNcalls(1); // useless for ZZ->4l
    Xcal2.XsecCalc(Proc,hzz4l_event, Xsec, XsecErr, verbosity);
    
    // fill in the differential cross-section and errors
    dXsecList = Xsec;  
    dXsecErrList = XsecErr;
    
    evt_tree->Fill();
    
  }//nevent
  
  if (verbosity >= TVar::INFO) cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
}  
