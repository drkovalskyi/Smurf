
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



const int kNDilep=4;

float lumi=1000.0;
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
float NLOXsec[kNProc] = {  4.5*0.919, 31314.0, 31314.0, 31314.0, 31314.0, 5.9, 0.249642, 0.452090, 0.641773, 0.770471, 0.866443, 0.782962, 0.659328, 0.486486, 0.408305, 0.358465, 0.321398, 0.290454, 0.243724, 0.175652};
// Note that MCFMXsec is obtained from a standalone MCFM calculations with no generator cuts applied,
// The following new numbers are including W->l branching fraction, ignoring currrently the Wgamma nubmers
// float MCFMXsec[kNProc] = { 2.983, 4032.1, 2694.8,  11270, 11270.0, 4.3, 0.07533, 0.1451, 0.2165, 0.2671, 0.3034, 0.2870, 0.2465, 0.1849, 0.1570, 0.1380, 0.1232, 0.1107, 0.0906, 0.05758};
float MCFMXsec[kNProc] = { 2.983, 16481.0, 10987.0,  11270, 11270.0, 4.3, 0.07533, 0.1451, 0.2165, 0.2671, 0.3034, 0.2870, 0.2465, 0.1849, 0.1570, 0.1380, 0.1232, 0.1107, 0.0906, 0.05758};

void CalculateAcceptance(){

    for (int i=0; i<kNProc; i++){
      for (int j=0; j<kNDilep; j++){
      float denominator = lumi* BR[i][j] * NLOXsec[i];
      if (denominator!=0) acceptance [i][j]= yield [i][j]/denominator; else acceptance [i][j]=1E+20;
      }
    }

}

void InitializeYields(TString inputDir){
  for (int i=0;i<kNProc;i++) {
    TFile *f = TFile::Open(TString(inputDir+name(i)+".root"));
    if(f==0x0) continue;
    cout << TString("../"+name(i)+".root") <<endl;
    TTree *tree = (TTree*) f->Get("tree");
    
    TH1F *hdphi[kNDilep];
    double tot_yield = 0;
    cout<<"Yields from file "<< name(i)<<".root: mm/me/em/ee:  ";    
    for (int j=0; j<kNDilep; j++){ 
      hdphi[j] = new TH1F(Form("hdphi_%i", j), Form("hdphi_%i", j), 20, 0, 4); 
      if(j==0 || j==2)
	 tree->Project(Form("hdphi_%i", j), "dPhi", Form("scale1fb*(type==%i)", j));
      if(j==1||j==3)
	tree->Project(Form("hdphi_%i", j), "dPhi", Form("scale1fb*(type==%i&&lep2.pt()>15)", j));
      if (hdphi!= 0x0) yield[i][j]=hdphi[j]->Integral(0, 9999); 
      else yield[i][j]=0;
      tot_yield+=yield[i][j];
      cout << Form("%.3f",yield[i][j])<<" " ;
    }       
    cout << "; total yield = " << Form("%.3f", tot_yield) << "\n"; 
  }
}


void NeutrinoIntegration(int process, TString inputDir, TString fileName, TString outputDir, int seed, int SmearLevel,int ncalls,int maxevt, int Mass);

//###################
//# main function
//###################
void runME_WWHWWWj(TString inputDir, TString fileName, TString outputDir, int seed,int SmearLevel,int ncalls,double Error, int Mass, int nev){

  ERRORthreshold=Error;
  int process=TVar::HWW;
  int maxevt=nev;

  InitializeYields(inputDir);
  CalculateAcceptance();

  cout <<"=== Neutrino Integration ==========" <<endl;  
  NeutrinoIntegration(process, inputDir, fileName, outputDir, seed, SmearLevel,ncalls,maxevt, Mass); 
 
}


void NeutrinoIntegration(int process,TString inputDir, TString fileName, TString outputDir, int seed, int SmearLevel,int ncalls,int maxevt, int Mass){

  cout << "Input File: "<< fileName <<" seed " << seed <<" SmearLevel " << SmearLevel <<" ncalls "<< ncalls<<endl;
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
  
  
  cout <<"Integration Seed= "<< Xcal2._seed << " SmearLevel= "<< Xcal2._smearLevel << " Ncalls = " << Xcal2._ncalls <<  endl;  
  
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
  printf("Total number of events = %d\n", Ntot);
  
  for(int ievt=0;ievt<Ntot;ievt++){
    
    for(int idx=0;idx<nProc;idx++) {
      dXsecList   [idx] = 0;
      dXsecErrList[idx] = 0;
    }
    
    ch->GetEntry(ievt);            
    // impose trailing electron pT > 15 GeV
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
      
    cout<<endl<<endl<<endl;
    cout <<"========================================================="<<endl;
    cout<<"Entry: "<<ievt<<"   Run: "<<run_<<"   Event: "<<event_<<endl;
    
    cout <<"Input: =================================================="<<endl;
    for(int lep=0;lep<2;lep++){
      printf("lep%d : %d %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	     lep, cms_event.PdgCode[lep], cms_event.p[lep].Px(),
	     cms_event.p[lep].Py(), cms_event.p[lep].Pz(),
	     cms_event.p[lep].Energy(), cms_event.p[lep].M());
    } 
    cout <<"========================================================="<<endl;
    
    
    // Processes considered
    int NProcessCalculate=0;
    TVar::Process processList[10];
    processList[ NProcessCalculate++]=TVar::WW;
    processList[ NProcessCalculate++]=TVar::HWW160;
    //processList[ NProcessCalculate++]=TVar::ZZ; 
    processList[ NProcessCalculate++]=TVar::Wp_1jet;
    processList[ NProcessCalculate++]=TVar::Wm_1jet;
    
    TVar::Process ProcInt;
    TVar::Process Vproc;
    
      double Xsec=0;
      double XsecErr=0;
      double Ratio=0;

      int HiggsMASS[20]={0,0,0,0,0,0,120,130,140,150,160,170,180,190,200,210,220,230,250,300};

      for(int iproc=0; iproc<NProcessCalculate; iproc++){ 
	
	// -- correct the lepton fo pt
	// For W+ hypothesis, assume the l- is the FO
	if(processList[iproc]==TVar::Wp_1jet) {
	  double scale_fo = 1.0;
	  if( TMath::Abs(cms_event.PdgCode[1]) == 11) {
	    scale_fo = Xcal2._FRhist.els_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.els_ptres->GetXaxis()->FindBin(cms_event.p[1].Pt()));
	    // cout << "TEvtProb::"<<__LINE__<< ": electron with pT "<< cms_event.p[1].Pt() << ", scale_fo = " << scale_fo<<endl;
	  }
	  if( TMath::Abs(cms_event.PdgCode[1]) == 13) {
	    scale_fo = Xcal2._FRhist.mus_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.mus_ptres->GetXaxis()->FindBin(cms_event.p[1].Pt()));
	    // cout << "TEvtProb::"<<__LINE__<< ": muon with pT "<< cms_event.p[1].Pt() << ", scale_fo = " << scale_fo<<endl;
	  }
	  cms_event.p[1].SetXYZM( cms_event.p[1].Px()*scale_fo, cms_event.p[1].Py()*scale_fo, cms_event.p[1].Py()*scale_fo, 0);
	}
	
	// For W- hypothesis, assume the l+ is the FO
	if(processList[iproc]==TVar::Wm_1jet) {
	  double scale_fo = 1.0;
	  if( TMath::Abs(cms_event.PdgCode[1]) == 11) {
	    scale_fo = Xcal2._FRhist.els_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.els_ptres->GetXaxis()->FindBin(cms_event.p[0].Pt()));
	    // cout << "TEvtProb::"<<__LINE__<< ": electron with pT "<< cms_event.p[0].Pt() << ", scale_fo = " << scale_fo<<endl;
	  }
	  if( TMath::Abs(cms_event.PdgCode[1]) == 13) {
	    scale_fo = Xcal2._FRhist.mus_ptres->ProfileX()->GetBinContent( Xcal2._FRhist.mus_ptres->GetXaxis()->FindBin(cms_event.p[0].Pt()));
	    // cout << "TEvtProb::"<<__LINE__<< ": muon with pT "<< cms_event.p[0].Pt() << ", scale_fo = " << scale_fo<<endl;
	  }
	  cms_event.p[0].SetXYZM( cms_event.p[0].Px()*scale_fo, cms_event.p[0].Py()*scale_fo, cms_event.p[0].Py()*scale_fo, 0);
	}
	

	ProcInt=processList[iproc];
	Xcal2.SetNcalls(ncalls);
	Xcal2.SetMCHist(ProcInt); 
	// Xcal2.SetFRHist(); 
	
	printf(" Calculate Evt %4i Run %9i Evt %8i Proc %4i %s Lep %4i %4i\n", ievt, run_, event_, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type);

	if ((processList[iproc]>=TVar::HWW120) && (processList[iproc]<=TVar::HWW300)){
	  Xcal2.SetProcess(TVar::HWW);
	  Xcal2.SetHiggsMass(HiggsMASS[processList[iproc]]);
	  cout<<HiggsMASS[processList[iproc]]<<"\n";
	  
	  if ( HiggsMASS[processList[iproc]]<160.8 )
	    Xcal2.SetHWWPhaseSpace(TVar::MH);
	  else
	    Xcal2.SetHWWPhaseSpace(TVar::MHMW);

	  Xcal2.NeutrinoIntegrate(TVar::HWW,cms_event, &Xsec, &XsecErr);
	} 
	else Xcal2.NeutrinoIntegrate(ProcInt,cms_event, &Xsec, &XsecErr);
	 
	dXsecList[ProcInt]=Xsec;  
	dXsecErrList[ProcInt]=XsecErr;
	if (Xsec>0) Ratio = XsecErr/Xsec;
	
	if (Ratio>ERRORthreshold){
	  cout <<"IntegrateNTrials "<<" Ncalls " << Xcal2._ncalls<<" "<<Vproc<<" "<< TVar::ProcessName(Vproc)
	       <<" Xsec = " <<  Xsec << " +- " << XsecErr << " ( " << Ratio << " ) " << endl;
	}
      }
            

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

      
      double denom=0.0, numer=0.0;
      double yield_bg=0.0;
      
      // Begin HWW likelihood

      cout << "START LR Construction ============" <<endl;

      for(int k=proc_HWW120; k<proc_HWW300; k++){
        yield_bg = yield[proc_WW][type_]+yield[proc_Wpj][type_];
	cout<<"bg_yield="<<yield_bg<<"\n";
	numer=1/(MCFMXsec[k]*acceptance[k][type_]) * dXsecList[k];
	denom  = numer;
	cout<< "PHWW " << name(k) << " "  <<numer<<"\n";

      	denom += 1/(MCFMXsec[proc_WW ]*acceptance[proc_WW ][type_]) * dXsecList[proc_WW ] * yield[proc_WW][type_]/yield_bg;
	cout<<" PWW= "<< 1/(MCFMXsec[proc_WW]*acceptance[proc_WW][type_]) * dXsecList[proc_WW] * yield[proc_WW][type_]/yield_bg <<"\n";	

	//denom += 1/(MCFMXsec[proc_Wpj ]*acceptance[proc_Wpj ][type_]) * dXsecList[proc_Wpj ] * 0.5*yield[proc_Wpj][type_]/yield_bg;
	// cout<<" PWpj= "<< 1/(MCFMXsec[proc_Wpj ]*acceptance[proc_Wpj ][type_]) * dXsecList[proc_Wpj ] * 0.5*yield[proc_Wpj][type_]/yield_bg <<"\n";
	
	denom += 1/(MCFMXsec[proc_Wpj ]*acceptance[proc_Wpj ][type_]) * dXsecList[proc_Wpj ] *  MCFMXsec[proc_Wpj]/(MCFMXsec[proc_Wpj]+MCFMXsec[proc_Wmj]) * yield[proc_Wpj][type_]/yield_bg;
	cout<<" PWpj= "<<  1/(MCFMXsec[proc_Wpj ]*acceptance[proc_Wpj ][type_]) * dXsecList[proc_Wpj ] *  MCFMXsec[proc_Wpj]/(MCFMXsec[proc_Wpj]+MCFMXsec[proc_Wmj]) * yield[proc_Wpj][type_]/yield_bg << "\n";


	//denom += 1/(MCFMXsec[proc_Wmj ]*acceptance[proc_Wmj ][type_]) * dXsecList[proc_Wmj ] * 0.5*yield[proc_Wpj][type_]/yield_bg;
	//cout<<" PWmj= "<< 1/(MCFMXsec[proc_Wmj ]*acceptance[proc_Wmj ][type_]) * dXsecList[proc_Wmj ] * 0.5*yield[proc_Wpj][type_]/yield_bg <<"\n";	

	denom += 1/(MCFMXsec[proc_Wmj ]*acceptance[proc_Wmj ][type_]) * dXsecList[proc_Wmj ] * MCFMXsec[proc_Wmj]/(MCFMXsec[proc_Wpj]+MCFMXsec[proc_Wmj]) * yield[proc_Wpj][type_]/yield_bg;
	cout<<" PWmj= "<< 1/(MCFMXsec[proc_Wmj ]*acceptance[proc_Wmj ][type_]) * dXsecList[proc_Wmj ] * MCFMXsec[proc_Wmj]/(MCFMXsec[proc_Wpj]+MCFMXsec[proc_Wmj]) * yield[proc_Wpj][type_]/yield_bg << "\n" <<endl;

	if(denom!=0)
	  LR[k]=numer/denom;
	  cout<<"LR_HWW["<<k<<"]= "<<LR[k]<<"\n";
      }      
      
      evt_tree->Fill();
      
    }//nevent
    
    cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
    
    newfile->cd(); 
    evt_tree->Write(); 
    newfile->Close();
    
}  
