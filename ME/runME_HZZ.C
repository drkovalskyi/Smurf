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
#include "SmurfTree.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double ERRORthreshold=1.0;
using namespace std;

void NeutrinoIntegration(int process, TString inputDir, TString fileName, TString outputDir, int seed, int SmearLevel,int ncalls,int maxevt, int evtstart, TVar::VerbosityLevel verbosity);

//###################
//# main function
//###################
void runME_HZZ(TString inputDir, TString fileName, TString outputDir, int seed,int SmearLevel,int ncalls,double Error, int nev, int evtstart = 0, TVar::VerbosityLevel verbosity=TVar::INFO){

  ERRORthreshold=Error;
  int process=TVar::HZZ;
  int maxevt=nev;
  
  if (verbosity >= TVar::INFO) cout <<"=== Neutrino Integration ==========" <<endl;  
  NeutrinoIntegration(process, inputDir, fileName, outputDir, seed, SmearLevel,ncalls,maxevt, evtstart, verbosity); 
 
}

void NeutrinoIntegration(int process,TString inputDir, TString fileName, TString outputDir, int seed, int SmearLevel,int ncalls,int maxevt, int evtstart, TVar::VerbosityLevel verbosity){

  if (verbosity >= TVar::INFO) cout << "Input File: " << fileName << " seed " << seed << " SmearLevel " << SmearLevel << " ncalls " << ncalls << endl;
  
  TFile* fin = new TFile(inputDir+fileName);
  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll(".root","_ME.root");
  cout << outFileName <<endl;
  TFile *newfile = new TFile(outFileName,"recreate");

  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // Declare the matrix element related variables to be added to the existing ntuples
  int ProcIdx=process; // what's this?
  int nProc=60;
  double dXsecList   [60];
  double dXsecErrList[60];
  double LR[60];
  
  evt_tree->Branch("ProcIdx"    ,&ProcIdx            ,"ProcIdx/I");
  evt_tree->Branch("nProc"      ,&nProc              ,"nProc/I");
  evt_tree->Branch("dXsec"      ,dXsecList           ,"dXsec[nProc]/D");
  evt_tree->Branch("dXsecErr"   ,dXsecErrList        ,"dXsecErr[nProc]/D");
  evt_tree->Branch("LR"   ,      LR                  ,"LR[nProc]/D");
  
  // Initialize the branches to use to calculate the differential cross-sections
  unsigned int run_ = 0;
  unsigned int event_ = 0;
  unsigned int njets_ = 0;
  LorentzVector*  lep1_ = 0;
  LorentzVector*  lep2_ = 0;
  int lq1_ = 0;
  int lq2_ = 0;
  int type_ = 0;
  float metPhi_ = 0;
  float met_ = 0;
  int dstype_ = 0;
  
  
  ch->SetBranchAddress( "run"        , &run_       );   
  ch->SetBranchAddress( "event"      , &event_     );   
  ch->SetBranchAddress( "lep1"       , &lep1_      );   
  ch->SetBranchAddress( "lep2"       , &lep2_      ); 
  ch->SetBranchAddress( "njets"      , &njets_     );     
  ch->SetBranchAddress( "lq1"        , &lq1_       );     
  ch->SetBranchAddress( "lq2"        , &lq2_       );     
  ch->SetBranchAddress( "type"       , &type_      );   
  ch->SetBranchAddress( "met"        , &met_       );       
  ch->SetBranchAddress( "metPhi"     , &metPhi_    );       
  ch->SetBranchAddress( "dstype"     , &dstype_    );   
  
  // The reconstructed event level information  
  // Create the instance of TEvtProb to calculate the differential cross-section
  // and set the calculation global variables common for all processes
  
  TEvtProb Xcal2;  
  Xcal2.SetApplyFake(1);
  Xcal2.SetSmearLevel(SmearLevel);
  Xcal2.SetSeed(seed);
  Xcal2.SetMatrixElement(TVar::MCFM);

  
  if (verbosity >= TVar::INFO) cout << "Integration Seed= "<< Xcal2._seed << " SmearLevel= " << Xcal2._smearLevel << " Ncalls = " << Xcal2._ncalls <<  endl;  
  cdf_event_type cms_event;
  
  //==========================================
  // Loop All Events
  //==========================================
 
  int Ntot = ch->GetEntries();
  if(evtstart+maxevt<Ntot) Ntot=evtstart+maxevt;
  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);
  
  for(int ievt=evtstart;ievt<Ntot;ievt++){
 
    if (verbosity >= TVar::INFO && (ievt % 5 == 0)) 
        std::cout << "Doing Event: " << ievt << std::endl;
 
    for(int idx=0;idx<nProc;idx++) {
      dXsecList   [idx] = 0;
      dXsecErrList[idx] = 0;
    }
    
    ch->GetEntry(ievt);            
    
    // impose trailing electron pT > 15 GeV (obselete for SmurfV5)
    // if( (type_ == 1||type_ == 3) && lep2_->Pt()<15 ) continue;

    // analyse only the 0-jet bin
    if (njets_ > 0) continue;
    
    // Initialize the reco-level cms_event, which contains the complete event
    // for the differential cross-section calculation
    
    int lep1_Type (0);
    int lep2_Type (0);
    
    switch (type_) {
    case 0: //mm
      lep1_Type = 13;
      lep2_Type = 13;
      break;
    case 1: // me
      lep1_Type = 13;
      lep2_Type = 11;
      break;
    case 2: // em
      lep1_Type = 11;
      lep2_Type = 13;
      break;
    case 3: // ee
      lep1_Type = 11;
      lep2_Type = 11;
      break;
    default:
      break;
    }
    
    cms_event.MetX = met_*cos(metPhi_);
    cms_event.MetY = met_*sin(metPhi_);
    
    if (lq1_>lq2_){
      cms_event.p[0].SetXYZM(lep1_->Px(), lep1_->Py(), lep1_->Pz(), 0.0);
      cms_event.p[1].SetXYZM(lep2_->Px(), lep2_->Py(), lep2_->Pz(), 0.0);
      cms_event.PdgCode[0]=lep1_Type;
      cms_event.PdgCode[1]=lep2_Type;    
    }
    else{
      cms_event.p[0].SetXYZM(lep2_->Px(), lep2_->Py(), lep2_->Pz(), 0.0);
      cms_event.p[1].SetXYZM(lep1_->Px(), lep1_->Py(), lep1_->Pz(), 0.0);
      cms_event.PdgCode[0]=lep2_Type;
      cms_event.PdgCode[1]=lep1_Type;    
    }
    
  if (verbosity >= TVar::DEBUG) {
    cout << "\n=========================================================\n";
    cout << "Entry: " << ievt << "   Run: " << run_ << "   Event: " << event_ <<endl;
    cout << "Input: ==================================================" <<endl;
    printf("lep1(%i) p3 = (%4.4f, %4.4f, %4.4f)  lep2(%i) p3 = (%4.4f, %4.4f, %4.4f)\n",
	   lep1_Type*lq1_, lep1_->Px(), lep1_->Py(), lep1_->Pz(), 
	   lep2_Type*lq2_, lep2_->Px(), lep2_->Py(), lep2_->Pz());
    cout << "=========================================================\n";
  }
  // finish loading event information
  
  // ==== Begin the differential cross-section calculation
  TVar::Process ProcInt;
  double Xsec=0;
  double XsecErr=0;
  double Ratio=0;
  int HiggsMASS[30]={0,0,0,0,0,115,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600,200,250,300,400};
  
  // The SM background processes
  int NProcessCalculate=0;
  TVar::Process processList[30];
  
  processList[ NProcessCalculate++]=TVar::WW;
  //  processList[ NProcessCalculate++]=TVar::Wp_1jet;
  //  processList[ NProcessCalculate++]=TVar::Wm_1jet;
  processList[ NProcessCalculate++]=TVar::ZZ;
  processList[ NProcessCalculate++]=TVar::WZ;
  
  // adjust the number of processes to run according to the datatype
  // run single higgs mass hypothesis for those hww MC
  // This is kind of dangerous, as we reply on the smurfntuples 
  // to be filled correctly...
  // processList[ NProcessCalculate++]=TVar::HZZ;

 switch ( dstype_) {
  case SmurfTree::hww200:
  case SmurfTree::vbfhww200:
    processList[ NProcessCalculate++]=TVar::HZZ200;
    break;
  case SmurfTree::hww250:
  case SmurfTree::vbfhww250:
    processList[ NProcessCalculate++]=TVar::HZZ250;
    break;
  case SmurfTree::hww300:
  case SmurfTree::vbfhww300:
    processList[ NProcessCalculate++]=TVar::HZZ300;
    break;
  case SmurfTree::hww400:
  case SmurfTree::vbfhww400:
    processList[ NProcessCalculate++]=TVar::HZZ400;
    break;
  default:
    processList[ NProcessCalculate++]=TVar::HZZ200;
    processList[ NProcessCalculate++]=TVar::HZZ250;
    processList[ NProcessCalculate++]=TVar::HZZ300;
    processList[ NProcessCalculate++]=TVar::HZZ400;
    break;
 }
 
  
    
  for(int iproc = 0; iproc < NProcessCalculate; iproc++) { 
    ProcInt=processList[iproc];
    
    // Load the MC based boost, efficiency and generator FR
    bool setMCFR = false;
    Xcal2.SetMCHist(ProcInt, "Util.root", setMCFR, verbosity);
    
    // Set the fakerate from data measurements
    if (!setMCFR && (iproc == TVar::Wp_1jet || iproc == TVar::Wm_1jet)) {
      TString elFRFile = "ElectronFakeRates_SmurfV5_V4.root";
      TString elFRHist = "ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta";
      TString muFRFile = "MuonFakeRate_SmurfV5_M1.root";
      TString muFRHist = "frEtaPt";
      Xcal2.SetFRHist(elFRFile, elFRHist, muFRFile, muFRHist, verbosity);
    }
    if (verbosity >= TVar::DEBUG) 
      printf(" Calculate Evt %4i Run %9i Evt %8i Proc %4i %s Lep %4i %4i\n", ievt, run_, event_, ProcInt, TVar::ProcessName(ProcInt).Data(), lq1_*lep1_Type, lq2_*lep2_Type);
    
    // -- correct the lepton fo pt for W + jet only
    // For W+ hypothesis, assume the l- is the FO
    double scale_fo = 1.0;
    if(ProcInt == TVar::Wp_1jet) {
      if( TMath::Abs(cms_event.PdgCode[1]) == 11) {
	scale_fo = Xcal2._FRhist.els_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.els_ptres->GetXaxis()->FindBin(cms_event.p[1].Pt()));
      }
      if( TMath::Abs(cms_event.PdgCode[1]) == 13) {
	scale_fo = Xcal2._FRhist.mus_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.mus_ptres->GetXaxis()->FindBin(cms_event.p[1].Pt()));
      }
      cms_event.p[1].SetXYZM( cms_event.p[1].Px()*scale_fo, cms_event.p[1].Py()*scale_fo, cms_event.p[1].Pz()*scale_fo, 0);
    }
    
    // For W- hypothesis, assume the l+ is the FO
    if(ProcInt == TVar::Wm_1jet) {
      if( TMath::Abs(cms_event.PdgCode[0]) == 11) {
	scale_fo = Xcal2._FRhist.els_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.els_ptres->GetXaxis()->FindBin(cms_event.p[0].Pt()));
      }
      if( TMath::Abs(cms_event.PdgCode[0]) == 13) {
	scale_fo = Xcal2._FRhist.mus_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.mus_ptres->GetXaxis()->FindBin(cms_event.p[0].Pt()));
      }	  
      cms_event.p[0].SetXYZM( cms_event.p[0].Px()*scale_fo, cms_event.p[0].Py()*scale_fo, cms_event.p[0].Pz()*scale_fo, 0);
    }
    
    // -- Do the integration over unknown parameters
    if (ProcInt >= TVar::HZZ200 && ProcInt <= TVar::HZZ400){
      Xcal2.SetProcess(TVar::HZZ);
      Xcal2.SetHiggsMass(HiggsMASS[ProcInt]);
      Xcal2.SetNcalls(ncalls);
      if (verbosity >= TVar::DEBUG) cout<< "Higgs Mass: " << HiggsMASS[ProcInt]<<"GeV \n";    
      Xcal2.NeutrinoIntegrate(TVar::HZZ,cms_event, &Xsec, &XsecErr, verbosity);
    }
    else {
      Xcal2.SetNcalls(ncalls);
      Xcal2.NeutrinoIntegrate(ProcInt,cms_event, &Xsec, &XsecErr, verbosity);
    }
    
    // fill in the differential cross-section and errors
    dXsecList[ProcInt]=Xsec;  
    dXsecErrList[ProcInt]=XsecErr;
    
    if (Xsec>0) Ratio = XsecErr/Xsec;
    if (Ratio > ERRORthreshold){
      if (verbosity >= TVar::DEBUG) 
	cout    << "IntegrateNTrials " << " Ncalls " << Xcal2._ncalls << " " << ProcInt << " " << TVar::ProcessName(ProcInt)
		<< " Xsec = " <<  Xsec << " +- " << XsecErr << " ( " << Ratio << " ) " << endl;
    }
    
    // setting back the lepton pT for the W + jet hypothesis
    if(ProcInt == TVar::Wp_1jet) 
      cms_event.p[1].SetXYZM( cms_event.p[1].Px()/scale_fo, cms_event.p[1].Py()/scale_fo, cms_event.p[1].Pz()/scale_fo, 0);
    
    if(ProcInt == TVar::Wm_1jet) 
      cms_event.p[0].SetXYZM( cms_event.p[0].Px()/scale_fo, cms_event.p[0].Py()/scale_fo, cms_event.p[0].Pz()/scale_fo, 0);
  } 
  
  // print some event based information for debugging purposes
  if (verbosity >= TVar::DEBUG) {
    double Mll  =(cms_event.p[0]+cms_event.p[1]).M();
    double Phill=TVector2::Phi_0_2pi(cms_event.p[0].DeltaPhi(cms_event.p[1]));
    if(Phill>TMath::Pi()) Phill=2*TMath::Pi()-Phill;			   
    
    cout << "START summary ====================" <<endl;
    printf(" Evt %4i/%4i Run %9i Evt %8i Proc %4i %s lep %4i %4i njets %d \n", ievt,Ntot, run_, event_, ProcIdx, TVar::ProcessName(ProcIdx).Data(), lq1_*lep1_Type, lq2_*lep2_Type,njets_);
    printf(" MetX %8.8f MetY %8.8f Mll %8.8f Phill %8.8f\n",cms_event.MetX, cms_event.MetY,Mll,Phill);
    printf(" L1 %d ( %8.8f , %8.8f , %8.8f , %8.8f)\n",cms_event.PdgCode[0],cms_event.p[0].Px(),cms_event.p[0].Py(),cms_event.p[0].Pz(),cms_event.p[0].Energy());
    printf(" L2 %d ( %8.8f , %8.8f , %8.8f , %8.8f)\n",cms_event.PdgCode[1],cms_event.p[1].Px(),cms_event.p[1].Py(),cms_event.p[1].Pz(),cms_event.p[1].Energy());
    
    for(int j=0;j<NProcessCalculate;j++){
      TVar::Process proc = processList[j];
      printf("%2i %8s  :", proc,TVar::ProcessName(proc).Data());
      double ratio = 0.0;
      if(dXsecList[processList[j]] > 0.0) ratio=dXsecErrList[processList[j]]/dXsecList[processList[j]];
      cout <<  dXsecList[processList[j]] << " +- " << dXsecErrList[processList[j]] << " ( " << ratio <<" ) "<<endl;
    }
    cout << "END summary =====================" <<endl;
  } 
  
  evt_tree->Fill();
  
  }//nevent
  
  if (verbosity >= TVar::INFO) cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
}  
