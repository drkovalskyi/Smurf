/*
  test version
  nohup root -l runME_test.C\+\(\"/smurf/data/Run2011_Spring11_SmurfV3/tas-zerojet/\"\,\"ww.root\"\,\"./\"\,10\,1\,100000\,1.0\,0.0\,500\) >& ww_ME_test.log &
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

enum process {
  proc_WW,
  proc_Wpj,
  proc_Wmj,
  proc_Wpg,
  proc_Wmg,
  proc_ZZ,
  proc_HWW120,
  proc_HWW130,
  proc_HWW140,
  proc_HWW150,
  proc_HWW160,
  proc_HWW170,
  proc_HWW180,
  proc_HWW190,
  proc_HWW200,
  proc_HWW210,
  proc_HWW220,
  proc_HWW230,
  proc_HWW250,
  proc_HWW300,
  kNProc
};

// this is needed to interface with the smurf trees
static std::string name(int  dstype){
  switch (dstype){
  case proc_WW:     return "ww";
  case proc_Wpj:    return "wjets";
  case proc_Wmj:    return "wjets";
  case proc_Wpg:    return "wgamma";
  case proc_Wmg:    return "wgamma";
  case proc_ZZ:     return "zz";
  case proc_HWW120: return "hww120";
  case proc_HWW130: return "hww130";
  case proc_HWW140: return "hww140";
  case proc_HWW150: return "hww150";
  case proc_HWW160: return "hww160";
  case proc_HWW170: return "hww170";
  case proc_HWW180: return "hww180";
  case proc_HWW190: return "hww190";
  case proc_HWW200: return "hww200";
  case proc_HWW210: return "hww210";
  case proc_HWW220: return "hww220";
  case proc_HWW230: return "hww230";
  case proc_HWW250: return "hww250";
  case proc_HWW300: return "hww300";
  default:     return "uknown";
  }
};

void NeutrinoIntegration(int process, TString inputDir, TString fileName, TString outputDir, int seed, int SmearLevel,int ncalls,int maxevt, int Mass, TVar::VerbosityLevel verbosity);

//###################
//# main function
//###################
void runME_test(TString inputDir, TString fileName, TString outputDir, int seed,int SmearLevel,int ncalls,double Error, int Mass, int nev, TVar::VerbosityLevel verbosity=TVar::INFO){

  ERRORthreshold=Error;
  int process=TVar::HWW;
  int maxevt=nev;
  
  if (verbosity >= TVar::INFO) cout <<"=== Neutrino Integration ==========" <<endl;  
  NeutrinoIntegration(process, inputDir, fileName, outputDir, seed, SmearLevel,ncalls,maxevt, Mass, verbosity); 
 
}


void NeutrinoIntegration(int process,TString inputDir, TString fileName, TString outputDir, int seed, int SmearLevel,int ncalls,int maxevt, int Mass, TVar::VerbosityLevel verbosity){

  if (verbosity >= TVar::INFO) cout << "Input File: " << fileName << " seed " << seed << " SmearLevel " << SmearLevel << " ncalls " << ncalls << endl;

  TFile* fin = new TFile(inputDir+fileName);
  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll(".root","_ME.root");
  cout << outFileName <<endl;
  TFile *newfile = new TFile(outFileName,"recreate");
  
  TEvtProb Xcal2;  
  
  Xcal2.SetApplyFake(1);
  Xcal2.SetSmearLevel(SmearLevel);
  Xcal2.SetSeed(seed);
  Xcal2.SetMatrixElement(TVar::MCFM);
  Xcal2.SetNcalls(ncalls);
  
  
  if (verbosity >= TVar::INFO) cout << "Integration Seed= "<< Xcal2._seed << " SmearLevel= " << Xcal2._smearLevel << " Ncalls = " << Xcal2._ncalls <<  endl;  
  
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  int ProcIdx=process;
  int nProc=60;
  double dXsecList   [60];
  double dXsecErrList[60];
  double LR[60];
  
  evt_tree->Branch("ProcIdx"    ,&ProcIdx            ,"ProcIdx/I");
  evt_tree->Branch("nProc"      ,&nProc              ,"nProc/I");
  evt_tree->Branch("dXsec"      ,dXsecList           ,"dXsec[nProc]/D");
  evt_tree->Branch("dXsecErr"   ,dXsecErrList        ,"dXsecErr[nProc]/D");
  evt_tree->Branch("LR"   ,      LR                  ,"LR[nProc]/D");
  
  cdf_event_type cms_event;
  
  // Initialize the branches to use
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
    
  
  ch->SetBranchAddress( "run"        , &run_       );   
  ch->SetBranchAddress( "event"      , &event_     );   
  ch->SetBranchAddress( "lep1"       , &lep1_      );   
  ch->SetBranchAddress( "lep2"       , &lep2_      ); 
  ch->SetBranchAddress( "njets"      , &njets_     );     
  ch->SetBranchAddress( "lq1"        , &lq1_       );     
  ch->SetBranchAddress( "lq2"        , &lq2_       );     
  ch->SetBranchAddress( "type"        , &type_       );   
  ch->SetBranchAddress( "met"        , &met_       );       
  ch->SetBranchAddress( "metPhi"        , &metPhi_       );       
  
  
  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = ch->GetEntries();
  if(maxevt<Ntot) Ntot=maxevt;
  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);
  
  for(int ievt=0;ievt<Ntot;ievt++){
   
    if (verbosity >= TVar::INFO && (ievt % 5 == 0)) 
        std::cout << "Doing Event: " << ievt << std::endl;
 
    for(int idx=0;idx<nProc;idx++) {
      dXsecList   [idx] = 0;
      dXsecErrList[idx] = 0;
    }
    
    ch->GetEntry(ievt);            
    // impose trailing electron pT > 15 GeV
    //    if (type_!=3) continue;
    if((type_==1||type_==3) && lep2_->Pt()<15 ) continue;

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
    
    double Mll  =(cms_event.p[0]+cms_event.p[1]).M();
    double Phill=TVector2::Phi_0_2pi(cms_event.p[0].DeltaPhi(cms_event.p[1]));
    if(Phill>TMath::Pi()) Phill=2*TMath::Pi()-Phill;			   
      

    if (verbosity >= TVar::DEBUG) {
        cout<<endl<<endl<<endl;
        cout << "=========================================================" <<endl;
        cout << "Entry: " << ievt << "   Run: " << run_ << "   Event: " << event_ <<endl;
        cout << "Input: ==================================================" <<endl;
        for(int lep=0;lep<2;lep++){
            printf("lep%d : %d %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	        lep, cms_event.PdgCode[lep], cms_event.p[lep].Px(),
	         cms_event.p[lep].Py(), cms_event.p[lep].Pz(),
	         cms_event.p[lep].Energy(), cms_event.p[lep].M());
        } 
        cout << "========================================================="<<endl;
    }
    
    // Processes considered
    int NProcessCalculate=0;
    TVar::Process processList[10];
    processList[ NProcessCalculate++]=TVar::WW;
    processList[ NProcessCalculate++]=TVar::HWW130;
    processList[ NProcessCalculate++]=TVar::HWW160;
    processList[ NProcessCalculate++]=TVar::HWW200;
    processList[ NProcessCalculate++]=TVar::Wp_1jet;
    processList[ NProcessCalculate++]=TVar::Wm_1jet;
    //processList[ NProcessCalculate++]=TVar::ZZ; 
    
    TVar::Process ProcInt;
    TVar::Process Vproc;
    
      double Xsec=0;
      double XsecErr=0;
      double Ratio=0;

      int HiggsMASS[20]={0,0,0,0,0,0,120,130,140,150,160,170,180,190,200,210,220,230,250,300};

      for(int iproc=0; iproc<NProcessCalculate; iproc++){ 
	
	// -- correct the lepton fo pt
	// For W+ hypothesis, assume the l- is the FO
	double scale_fo = 1.0;
	if(processList[iproc]==TVar::Wp_1jet) {
	  if( TMath::Abs(cms_event.PdgCode[1]) == 11) {
	    scale_fo = Xcal2._FRhist.els_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.els_ptres->GetXaxis()->FindBin(cms_event.p[1].Pt()));
	  }
	  if( TMath::Abs(cms_event.PdgCode[1]) == 13) {
	    scale_fo = Xcal2._FRhist.mus_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.mus_ptres->GetXaxis()->FindBin(cms_event.p[1].Pt()));
	  }
	  cms_event.p[1].SetXYZM( cms_event.p[1].Px()*scale_fo, cms_event.p[1].Py()*scale_fo, cms_event.p[1].Pz()*scale_fo, 0);
	}
	
	// For W- hypothesis, assume the l+ is the FO
	if(processList[iproc]==TVar::Wm_1jet) {
	  if( TMath::Abs(cms_event.PdgCode[0]) == 11) {
	    scale_fo = Xcal2._FRhist.els_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.els_ptres->GetXaxis()->FindBin(cms_event.p[0].Pt()));
	  }
	  if( TMath::Abs(cms_event.PdgCode[0]) == 13) {
	    scale_fo = Xcal2._FRhist.mus_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.mus_ptres->GetXaxis()->FindBin(cms_event.p[0].Pt()));
	  }	  
	  cms_event.p[0].SetXYZM( cms_event.p[0].Px()*scale_fo, cms_event.p[0].Py()*scale_fo, cms_event.p[0].Pz()*scale_fo, 0);
	}
	
	ProcInt=processList[iproc];
	Xcal2.SetNcalls(ncalls);
	Xcal2.SetMCHist(ProcInt, verbosity); 
	// Xcal2.SetFRHist(); 
	
	if (verbosity >= TVar::DEBUG) printf(" Calculate Evt %4i Run %9i Evt %8i Proc %4i %s Lep %4i %4i\n", ievt, run_, event_, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type);

	if ((processList[iproc]>=TVar::HWW120) && (processList[iproc]<=TVar::HWW300)){
	  Xcal2.SetProcess(TVar::HWW);
	  Xcal2.SetHiggsMass(HiggsMASS[processList[iproc]]);
	  if (verbosity >= TVar::DEBUG) cout<< "Higgs Mass: " << HiggsMASS[processList[iproc]]<<"GeV \n";
	  
	  if ( HiggsMASS[processList[iproc]]<160.8 )
	    Xcal2.SetHWWPhaseSpace(TVar::MH);
	  else
	    Xcal2.SetHWWPhaseSpace(TVar::MHMW);

	  Xcal2.NeutrinoIntegrate(TVar::HWW,cms_event, &Xsec, &XsecErr, verbosity);
	} 
	else Xcal2.NeutrinoIntegrate(ProcInt,cms_event, &Xsec, &XsecErr, verbosity);
	 
	dXsecList[ProcInt]=Xsec;  
	dXsecErrList[ProcInt]=XsecErr;
	if (Xsec>0) Ratio = XsecErr/Xsec;
	
	if (Ratio > ERRORthreshold){
	  if (verbosity >= TVar::ERROR) 
            cout    << "IntegrateNTrials " << " Ncalls " << Xcal2._ncalls << " " << Vproc << " " << TVar::ProcessName(Vproc)
	                << " Xsec = " <<  Xsec << " +- " << XsecErr << " ( " << Ratio << " ) " << endl;
	}

	// setting back the cms_event for other hypothesis
	if(processList[iproc]==TVar::Wp_1jet) 
	  cms_event.p[1].SetXYZM( cms_event.p[1].Px()/scale_fo, cms_event.p[1].Py()/scale_fo, cms_event.p[1].Pz()/scale_fo, 0);
	
	if(processList[iproc]==TVar::Wm_1jet) 
	  cms_event.p[0].SetXYZM( cms_event.p[0].Px()/scale_fo, cms_event.p[0].Py()/scale_fo, cms_event.p[0].Pz()/scale_fo, 0);
		  
      }
      

    if (verbosity >= TVar::DEBUG) {
        cout << "START summary ====================" <<endl;
        printf(" Evt %4i/%4i Run %9i Evt %8i Proc %4i %s lep %4i %4i njets %d \n", ievt,Ntot, run_, event_, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type,njets_);
        printf(" MetX %8.8f MetY %8.8f Mll %8.8f Phill %8.8f\n",cms_event.MetX, cms_event.MetY,Mll,Phill);
        printf(" L1 %d ( %8.8f , %8.8f , %8.8f , %8.8f)\n",cms_event.PdgCode[0],cms_event.p[0].Px(),cms_event.p[0].Py(),cms_event.p[0].Pz(),cms_event.p[0].Energy());
        printf(" L2 %d ( %8.8f , %8.8f , %8.8f , %8.8f)\n",cms_event.PdgCode[1],cms_event.p[1].Px(),cms_event.p[1].Py(),cms_event.p[1].Pz(),cms_event.p[1].Energy());
      
        for(int j=0;j<NProcessCalculate;j++){
	        TVar::Process proc = processList[j];
	        printf("%2i %8s  :", proc,TVar::ProcessName(proc).Data());
	        double ratio=0;
	        if(dXsecList[processList[j]]>0) ratio=dXsecErrList[processList[j]]/dXsecList[processList[j]];
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
