#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include "TRandom.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include <iomanip>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "TH2D.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HZZllvv/factors.h"
#include "Smurf/HZZ/HiggsQCDScaleSystematics.h"
#include "Smurf/HZZ/PSUESystematics.h"
#include "Smurf/HZZ/LeptonScaleLookup.h"
#include "Smurf/Analysis/HZZllvv/PileupReweighting.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void plotHistsInPad(TH1D* h1, TH1D* h2);
void setHist(TH1D* h, int color, int style);
void setPair(TH1D* h1, TH1D* h2);
void setGraph(TGraphErrors* g, int color, int marker);
TGraphErrors* makeSignificanceCurve(TH1D* sig, TH1D* bgd, TH1D* dat, const char* name);
TGraphErrors* makeGraphFromHists   (TH1D* sig, TH1D* bgd, const char* name);

int    verboseLevel =   0;
const double sigmaB = 0.35;

//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void MakeHZZRes
(
 UInt_t  nJetsType   	 = 0,
 UInt_t  mH      	 = 250,
 TString outTag          = "default",
 TString sigInputFile    = "data/ntuples_250train_hzzjets_hzz250.root",
 TString bgdInputFile    = "data/ntuples_250train_hzzjets_zz.root",
 TString datInputFile    = "data/ntuples_250train_hzzjets_data_2l.root",
 Int_t   wwDecay         = 0,
 bool fillInfoNote       = false,
 int  TheVerboseLevel    = 0
 )
{

  //************************************************************************************************
  // Setup
  //************************************************************************************************
  verboseLevel = TheVerboseLevel;

  TString sigFile1 = sigInputFile;
  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  TString c1Name = TString("plot_") + outTag + TString(".gif");
  TString c2Name = TString("plot_") + outTag + TString("_effvsbkg.gif");

  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("certifiedUCSD.json"); 


  //************************************************************************************************
  // Load Trees
  //************************************************************************************************

  TChain *chsignal = new TChain("tree");
  chsignal->Add(sigFile1);

   TChain *chbackground = new TChain("tree");
//   TChain *chbackground = new TChain("HwwTree1");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
//   TChain *chdata = new TChain("HwwTree0");
  chdata->Add(datFile1);

  TTree *signal     = (TTree*) chsignal;
  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;





  //************************************************************************************************
  // DY Bkg Photon+Jets Template Reweight Factors
  //************************************************************************************************
  TFile *fReweightFactorFile = TFile::Open("PhotonJetsReweightFactors.root");
  string photonTemplateLabel = "Spring11MC_BarrelOnlyNoCuts";
  TH2D *fReweightFactor = (TH2D*)(fReweightFactorFile->Get(("DYReweightFactorPtNJet_"+photonTemplateLabel).c_str()));
  assert(fReweightFactor);
  fReweightFactor->SetDirectory(0);
  TH2D *fReweightFactor_ee = (TH2D*)(fReweightFactorFile->Get(("DYReweightFactorPtNJet_ee_"+photonTemplateLabel).c_str()));
  assert(fReweightFactor_ee);
  fReweightFactor_ee->SetDirectory(0);
  TH2D *fReweightFactor_mm = (TH2D*)(fReweightFactorFile->Get(("DYReweightFactorPtNJet_mm_"+photonTemplateLabel).c_str()));
  fReweightFactor_mm->SetDirectory(0);
  assert(fReweightFactor_mm);

  fReweightFactorFile->Close();
  delete fReweightFactorFile;

  //************************************************************************************************
  // Lepton Efficiencies
  //************************************************************************************************
  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/Run2011_Spring11_SmurfV6/lepton_eff/efficiency_results_v6_187pb_newhistnames.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  LeptonScaleLookup trigLookup("/data/smurf/data/Run2011_Spring11_SmurfV6/lepton_eff/efficiency_results_v6_187pb_newhistnames.root");

  TFile *fLeptonFRFileM = TFile::Open("/data/smurf/sixie/FakeRates/FakeRates_SmurfV6.PromptRecoV1V2.200ipb.root");
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open("/data/smurf/sixie/FakeRates/FakeRates_SmurfV6.PromptRecoV1V2.200ipb.root");
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  //************************************************************************************************
  // Higgs Pt Reweighting
  //************************************************************************************************
  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ggHWW_KFactors_PowhegToHQT.root");
  TH1D *HiggsPtKFactor;
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", mH);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) {
    HiggsPtKFactor->SetDirectory(0);
  }
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;

  //************************************************************************************************
  // ZZ Dilepton Pt Reweighting
  //************************************************************************************************
  TFile *fZZKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ZZ_KFactor.root");
  TH1D *ZZKFactor;
  ZZKFactor = (TH1D*)(fZZKFactorFile->Get("KFactorZZ_DileptonPt"));
  if (ZZKFactor) {
    ZZKFactor->SetDirectory(0);
  }
  assert(ZZKFactor);
  fZZKFactorFile->Close();
  delete fZZKFactorFile;

  //************************************************************************************************
  // MC -> Data ScaleFactors
  //************************************************************************************************
  double TopAndWWScaleFactor[3]            = { 1.615, 1.671, 1.000 };
  double TopAndWWScaleFactorKappa[3]       = { 1.359, 1.440, 1.000 };


  //************************************************************************************************
  // Luminosity
  //************************************************************************************************

  double scaleFactorLum = 0.187;


  //************************************************************************************************
  // higgs masses
  //************************************************************************************************
  Int_t mHIndex = -1;
  if (mH == 250) mHIndex = 0;
  if (mH == 275) mHIndex = 1;
  if (mH == 300) mHIndex = 2;
  if (mH == 325) mHIndex = 3;
  if (mH == 350) mHIndex = 4;
  if (mH == 375) mHIndex = 5;
  if (mH == 400) mHIndex = 6;
  if (mH == 425) mHIndex = 7;
  if (mH == 450) mHIndex = 8;
  if (mH == 475) mHIndex = 9;
  if (mH == 500) mHIndex = 10;
  if (mH == 525) mHIndex = 11;
  if (mH == 550) mHIndex = 12;
  if (mH == 575) mHIndex = 13;
  if (mH == 600) mHIndex = 14;
  if (mHIndex == -1 ) return;



  //************************************************************************************************
  // Yields Tables and Histograms
  //************************************************************************************************



  const int nMVAVars = 2; //number of different MVA variables
  int    nBinHis[nMVAVars] = { 200,  200};
  double minHis[nMVAVars]  = {-1.0, -1.0};
  double maxHis[nMVAVars]  = { 2.0,  1.0};

  //----------------------------------------------------------------------------
  const UInt_t nSignalChannels = 5; //total, ggH, VBF H, WH/ZH, ttH
  const int nBkgChannels = 8; //number of different background channels
  TH1D* sigMVA[nMVAVars][6];
  TH1D* bgdMVA[nMVAVars];
  TH1D* datMVA[nMVAVars];
  TH1D* bgdMVADecays[nMVAVars][nBkgChannels];
  for(int i=0; i<nMVAVars; i++) {
    for(int j=0; j<6; j++){
      sigMVA[i][j] = new TH1D(Form("sigMVA_%d_%d",i,j), Form("sigMVA_%d_%d",i,j), nBinHis[i], minHis[i], maxHis[i]);
      sigMVA[i][j]->Sumw2();
    }
    bgdMVA[i] = new TH1D(Form("bgdMVA_%d",i), Form("bgdMVA_%d",i), nBinHis[i], minHis[i], maxHis[i]);
    datMVA[i] = new TH1D(Form("datMVA_%d",i), Form("datMVA_%d",i), nBinHis[i], minHis[i], maxHis[i]);
    bgdMVA[i]->Sumw2();
    datMVA[i]->Sumw2();
    for(int j=0; j<nBkgChannels; j++) {
      bgdMVADecays[i][j] = new TH1D(Form("bgdMVADecays_%d_%d",i,j), Form("bgdMVADecays_%d_%d",i,j), nBinHis[i], minHis[i], maxHis[i]);
      bgdMVADecays[i][j]->Sumw2();
    }
  }

  TH1D* histoSignal = new TH1D("histoSignal", "histoSignal", nBinHis[1], minHis[1], maxHis[1]);
  histoSignal->Sumw2();
  TH1D* histoBkg[nBkgChannels];
  for(int k=0; k<nBkgChannels; k++) {
    histoBkg[k] = (TH1D*) histoSignal->Clone(Form("histoBkg%d",k));
    histoBkg[k]->Scale(0.0);
  }
  TH1D* histoData = (TH1D*)histoSignal->Clone("histoData");

  //Yields
  float nSigAcc[nSignalChannels];
  float nSigCut[nSignalChannels];
  float nSigMVA[nSignalChannels];
  float nSigEAcc[nSignalChannels];
  float nSigECut[nSignalChannels];
  float nSigEMVA[nSignalChannels];
  for(int k=0; k<nSignalChannels; k++) {
    nSigAcc[k]  = 0.0;
    nSigCut[k]  = 0.0;
    nSigMVA[k]  = 0.0;
    nSigEAcc[k] = 0.0;
    nSigECut[k] = 0.0;
    nSigEMVA[k] = 0.0;    
  }

  float nBgdAcc = 0.0;
  float nBgdCut = 0.0;
  float nBgdMVA = 0.0;
  float nBgdEAcc = 0.0;
  float nBgdECut = 0.0;
  float nBgdEMVA = 0.0;
  float nBgdAccDecays[nBkgChannels];
  float nBgdCutDecays[nBkgChannels];
  float nBgdMVADecays[nBkgChannels];
  float nBgdEAccDecays[nBkgChannels];
  float nBgdECutDecays[nBkgChannels];
  float nBgdEMVADecays[nBkgChannels];
  for(int k=0; k<nBkgChannels; k++) {
    nBgdAccDecays[k]  = 0.0;
    nBgdCutDecays[k]  = 0.0;
    nBgdMVADecays[k]  = 0.0;
    nBgdEAccDecays[k] = 0.0;
    nBgdECutDecays[k] = 0.0;
    nBgdEMVADecays[k] = 0.0;    
  }

  float nDatAcc = 0.0;
  float nDatCut = 0.0;
  float nDatMVA = 0.0;


  //************************************************************************************************
  // mH dependant Cuts
  // [NJETS][mH]
  //************************************************************************************************
  double cutMETLow[3][15]        = { { 50, 50, 70 , 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150},
                                     { 60, 60, 100, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150},
                                     { 80, 80, 100, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150} };
  double cutMTLow[3][15]         = { { 180,180,260, 75, 80, 80, 90,110,120,120,120,120,120,120,120},
                                     { 180,180,260, 75, 80, 80, 90,110,120,120,120,120,120,120,120},
                                     { 180,180,240, 75, 80, 80, 90,110,120,120,120,120,120,120,120} };
  double cutMTHigh[3][15]        = { { 220,220,320,125,130,150,160,170,180,190,200,210,220,230,250}, 
                                     { 200,200,320,125,130,150,160,170,180,190,200,210,220,230,250}, 
                                     { 200,200,320,125,130,150,160,170,180,190,200,210,220,230,250} };


  // TMVA cuts (known a posteriori)
  int useVar = 0; // which MVA to be used
  double useCut = -999.0;
  if    (nJetsType == 0){
    if     (mH == 250){
      useVar = 0;
      useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 275){
      useVar = 0;
      useCut =  0.675-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 300){
      useVar = 1;
      useCut =  0.09-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 325){
      useVar = 0;
      useCut = 0.735-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 350){
      useVar = 0;
      useCut = 0.585-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 375){
      useVar = 0;
      useCut = 0.875-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 400){
      useVar = 0;
      useCut = 0.855-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 425){
      useVar = 0;
      useCut = 0.855-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 450){
      useVar = 0;
      useCut = 0.615-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 475){
      useVar = 0;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 500){
      useVar = 0;
      useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 300){
      useVar = 0;
      useCut = 0.805-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 350){
      useVar = 0;
      useCut = 0.835-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 400){
      useVar = 0;
      useCut = 0.845-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 450){
      useVar = 0;
      useCut = 0.895-datMVA[useVar]->GetBinWidth(0)/2.0;
    }
    else if(mH == 500){
      useVar = 0;
      useCut = 0.915-datMVA[useVar]->GetBinWidth(0)/2.0;
    }
    else if(mH == 550){
      useVar = 0;
      useCut = 0.925-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 600){
      useVar = 0;
      useCut = 0.945-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
  }
  else if(nJetsType == 1){
    if     (mH == 115){
      useVar = 0;
      useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 120){
      useVar = 0;
      useCut = 0.665-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 130){
      useVar = 0;
      useCut = 0.885-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 140){
      useVar = 0;
      useCut = 0.725-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 150){
      useVar = 0;
      useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 160){
      useVar = 0;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 170){
      useVar = 0;
      useCut = 0.835-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 180){
      useVar = 0;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 190){
      useVar = 0;
      useCut = 0.745-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 200){
      useVar = 0;
      useCut = 0.735-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 250){
      useVar = 0;
      useCut = 0.655-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 300){
      useVar = 0;
      useCut = 0.765-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 350){
      useVar = 0;
      useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 400){
      useVar = 0;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 450){
      useVar = 0;
      useCut = 0.885-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 500){
      useVar = 0;
      useCut = 0.915-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 550){
      useVar = 0;
      useCut = 0.925-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 600){
      useVar = 0;
      useCut = 0.945-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
  }


  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          njets;
  UInt_t          run;
  UInt_t          lumi;
  UInt_t          event;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  jet1  = 0;
  LorentzVector*  jet2  = 0;
  LorentzVector*  jet3  = 0;
  Float_t         dPhi;
  Float_t         dR;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         met;
  Float_t         metPhi;
  Float_t         trackMet;
  Float_t         mt;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         dPhiLep1MET;
  Float_t         dPhiLep2MET;
  Float_t         dPhiDiLepMET;
  Float_t         dPhiDiLepJet1;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;
  Int_t           processId;
  Float_t         jetLowBtag;
  UInt_t          nSoftMuons;
  Float_t         jet1Btag;
  Float_t         jet2Btag;
  Int_t 	  lep1McId;
  Int_t 	  lep2McId;
  Int_t 	  lep1MotherMcId;
  Int_t 	  lep2MotherMcId;
  Float_t         nn = 0.0;
  Float_t         bdtg = 0.0;
  Float_t         higgsPt = -999;



  //**************************************************************************
  //Signal Loop
  //**************************************************************************
  
  signal->SetBranchAddress( "cuts"         , &cuts         );
  signal->SetBranchAddress( "dstype"       , &dstype       );
  signal->SetBranchAddress( "nvtx"         , &nvtx         );
  signal->SetBranchAddress( "njets"        , &njets        );
  signal->SetBranchAddress( "event"        , &event        );
  signal->SetBranchAddress( "scale1fb"     , &scale1fb     );
  signal->SetBranchAddress( "lep1"         , &lep1         );
  signal->SetBranchAddress( "lep2"         , &lep2         );
  signal->SetBranchAddress( "jet1"         , &jet1         );
  signal->SetBranchAddress( "jet2"         , &jet2         );
  signal->SetBranchAddress( "jet3"         , &jet3         );
  signal->SetBranchAddress( "dPhi"         , &dPhi         );
  signal->SetBranchAddress( "dR"           , &dR           );
  signal->SetBranchAddress( "dilep"        , &dilep        );
  signal->SetBranchAddress( "type"         , &type         );
  signal->SetBranchAddress( "pmet"         , &pmet         );
  signal->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  signal->SetBranchAddress( "met"          , &met          );
  signal->SetBranchAddress( "metPhi"       , &metPhi       );
  signal->SetBranchAddress( "trackMet"     , &trackMet     );
  signal->SetBranchAddress( "mt"           , &mt           );
  signal->SetBranchAddress( "mt1"          , &mt1          );
  signal->SetBranchAddress( "mt2"          , &mt2          );
  signal->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  signal->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  signal->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  signal->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  signal->SetBranchAddress( "lq1"          , &lq1          );
  signal->SetBranchAddress( "lq2"          , &lq2          );
  signal->SetBranchAddress( "lid1"         , &lid1         );
  signal->SetBranchAddress( "lid2"         , &lid2         );
  signal->SetBranchAddress( "lid3"         , &lid3         );
  signal->SetBranchAddress( "processId"    , &processId    );
  signal->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  signal->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  signal->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  signal->SetBranchAddress( "jet2Btag"     , &jet2Btag     );
  signal->SetBranchAddress(Form("nn_hzz%i"	  ,mH), &nn	      );
  signal->SetBranchAddress(Form("bdtg_hzz%i"   ,mH), &bdtg         );
  signal->SetBranchAddress( "higgsPt"      , &higgsPt      );

  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)signal->GetEntries());

   //VBF Jet Selection (Central Jet Veto)
    unsigned int Njet3 = njets;
    if(nJetsType == 2){
      if(jet3->pt() <= 30)					         Njet3 = 2;
      else if(jet3->pt() > 30 && (
        (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
        (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))   Njet3 = 0;
      else								 Njet3 = 2;
      if(njets < 2 || njets > 3) Njet3 = 0;
    }


    //************************************************************************************************
    // PreSelection
    //************************************************************************************************
    if( !(lep1->pt() > 30 && lep2->pt() > 20 )) continue;
    if( !((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection )) continue;
    if( !((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection )) continue;
    if( !(type == SmurfTree::mm || type == SmurfTree::ee )) continue;
    if( !(fabs(dilep->mass() - 91.1876) < 15 )) continue;
    if( !(lq1*lq2 < 0)) continue;
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( !(dilep->pt() > 40)) continue;
    if( jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 
    if( dilep->pt() <= 40.0         ) continue; // cut on low dilepton mass
    
    if( Njet3 != nJetsType          ) continue; // select n-jet type events
    if( !(TMath::Min(met,trackMet) > 50.0)) continue;

    
    //************************************************************************************************
    // weights and scale factors
    //************************************************************************************************
    Double_t weightFactor = 1.0;
    Double_t effCorrection = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1) * 
      leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2) *
      trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                              fabs(lep2->eta()), lep2->pt(), 
                                              TMath::Abs(lid1), TMath::Abs(lid2));
    weightFactor *= effCorrection;    
    
    if (processId == 10010) {
      weightFactor *= HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(higgsPt));
    }


    double myWeight = scaleFactorLum * scale1fb * weightFactor;

    int nSigBin = -1;
    // 10010 == ggH, 10001 == VBF, 121/122 == ttH, 24 == WH, 26 == ZH
    if(processId==10010)           nSigBin = 1;
    else if(processId==10001)      nSigBin = 2;
    else if(processId==24)         nSigBin = 3;
    else if(processId==26)         nSigBin = 3;
    else if(processId==121 ||
            processId==122)        nSigBin = 4;
    else  {return;}

    nSigAcc[0]  = nSigAcc[0]  + myWeight;
    nSigEAcc[0] = nSigEAcc[0] + myWeight*myWeight;
    nSigAcc[nSigBin]  = nSigAcc[nSigBin]  + myWeight;
    nSigEAcc[nSigBin] = nSigEAcc[nSigBin] + myWeight*myWeight;

    if(useVar == 0)
      histoSignal   ->Fill(TMath::Max(TMath::Min((double)nn  ,maxHis[0]-0.001),minHis[0]+0.001),   myWeight);
    else if(useVar == 1)
      histoSignal   ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);

    sigMVA[0][0]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),     myWeight);
    sigMVA[1][0]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);
    sigMVA[0][nSigBin]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),	  myWeight);
    sigMVA[1][nSigBin]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);


    //**************************************************************************
    //Define HZZ Cuts
    //**************************************************************************
    double pxzll = dilep->px() + met*cos(metPhi);
    double pyzll = dilep->py() + met*sin(metPhi);
    double mtHZZ = TMath::Power(sqrt(dilep->pt()*dilep->pt() + dilep->M()*dilep->M())+
                                sqrt(met*met + dilep->M()*dilep->M()),2) - (pow(pxzll,2) + pow(pyzll,2));
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;


    bool passCutBased = ( (0==0)
                          && ( TMath::Min(met,trackMet) > cutMETLow[njets][mHIndex])
                          && (mtHZZ > cutMTLow[njets][mHIndex] && mtHZZ < cutMTHigh[njets][mHIndex]) );

    if (nJetsType == 2) {
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passCutBased = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
                    jet1->eta()*jet2->eta() < 0 &&
		    (mH > 200 || dilep->mass() < 100.) &&
		    centrality == 1;
    }

    if(passCutBased == true) {
      nSigCut[0]  +=  myWeight;
      nSigECut[0] +=  myWeight*myWeight;
      nSigCut[nSigBin]  += myWeight;
      nSigECut[nSigBin] += myWeight*myWeight;
    }

    bool passMVACut = false;
    if     (useVar == 0 && nn         > useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passMVACut = true;}
    else if(useVar == 1 && bdtg       > useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passMVACut = true;}
    if     (useVar == 0 && nn         > useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passMVACut = true;}
    else if(useVar == 1 && bdtg       > useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passMVACut = true;}
  }
  for(int i=0; i<nSignalChannels; i++) nSigEAcc[i] = sqrt(nSigEAcc[i]);
  for(int i=0; i<nSignalChannels; i++) nSigECut[i] = sqrt(nSigECut[i]);
  for(int i=0; i<nSignalChannels; i++) nSigEMVA[i] = sqrt(nSigEMVA[i]);
  printf("--- Finished Signal loop\n");



  //**************************************************************************
  //Bkg Loop
  //**************************************************************************

  background->SetBranchAddress( "cuts"          , &cuts 	  );
  background->SetBranchAddress( "dstype"        , &dstype	  );
  background->SetBranchAddress( "nvtx"          , &nvtx 	  );
  background->SetBranchAddress( "njets"         , &njets	  );
  background->SetBranchAddress( "event"         , &event	  );
  background->SetBranchAddress( "scale1fb"      , &scale1fb	  );
  background->SetBranchAddress( "lep1"          , &lep1 	  );
  background->SetBranchAddress( "lep2"          , &lep2 	  );
  background->SetBranchAddress( "jet1"          , &jet1 	  );
  background->SetBranchAddress( "jet2"          , &jet2 	  );
  background->SetBranchAddress( "jet3"          , &jet3 	  );
  background->SetBranchAddress( "dPhi"          , &dPhi 	  );
  background->SetBranchAddress( "dR"            , &dR		  );
  background->SetBranchAddress( "dilep"         , &dilep	  );
  background->SetBranchAddress( "type"          , &type 	  );
  background->SetBranchAddress( "pmet"          , &pmet 	  );
  background->SetBranchAddress( "pTrackMet"     , &pTrackMet	  );
  background->SetBranchAddress( "met"           , &met  	  );
  background->SetBranchAddress( "metPhi"        , &metPhi         );
  background->SetBranchAddress( "trackMet"      , &trackMet  	  );
  background->SetBranchAddress( "mt"            , &mt		  );
  background->SetBranchAddress( "mt1"           , &mt1  	  );
  background->SetBranchAddress( "mt2"           , &mt2  	  );
  background->SetBranchAddress( "dPhiLep1MET"   , &dPhiLep1MET    );
  background->SetBranchAddress( "dPhiLep2MET"   , &dPhiLep2MET    );
  background->SetBranchAddress( "dPhiDiLepMET"  , &dPhiDiLepMET   );
  background->SetBranchAddress( "dPhiDiLepJet1" , &dPhiDiLepJet1  );
  background->SetBranchAddress( "lq1"           , &lq1  	  );
  background->SetBranchAddress( "lq2"           , &lq2  	  );
  background->SetBranchAddress( "lid1"          , &lid1 	  );
  background->SetBranchAddress( "lid2"          , &lid2 	  );
  background->SetBranchAddress( "lid3"          , &lid3 	  );
  background->SetBranchAddress( "processId"     , &processId	  );
  background->SetBranchAddress( "jetLowBtag"    , &jetLowBtag	  );
  background->SetBranchAddress( "nSoftMuons"    , &nSoftMuons	  );
  background->SetBranchAddress( "jet1Btag"      , &jet1Btag	  );
  background->SetBranchAddress( "jet2Btag"      , &jet2Btag	  );
  background->SetBranchAddress( "lep1McId"      , &lep1McId	  );
  background->SetBranchAddress( "lep2McId"      , &lep2McId	  );
  background->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId );
  background->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId );
  background->SetBranchAddress(Form("nn_hzz%i"     ,mH), &nn   );
  background->SetBranchAddress(Form("bdtg_hzz%i"   ,mH), &bdtg );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    //VBF Jet Selection (Central Jet Veto)
    unsigned int Njet3 = njets;
    if(nJetsType == 2){
      if(jet3->pt() <= 30)                                               Njet3 = 2;
      else if(jet3->pt() > 30 && (
        (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
        (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))   Njet3 = 0;
      else                                                               Njet3 = 2;
      if(njets < 2 || njets > 3) Njet3 = 0;
    }

    //************************************************************************************************
    // PreSelection
    //************************************************************************************************
    if( njets != nJetsType ) continue; // select n-jet type events
    if( !(lep1->pt() > 30 && lep2->pt() > 20 )) continue;
    if( !((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection )) continue;
    if( !((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection )) continue;
    if( !(type == SmurfTree::mm || type == SmurfTree::ee )) continue;
    if( !(fabs(dilep->mass() - 91.1876) < 15 )) continue;
    if( !(lq1*lq2 < 0)) continue;
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 
    if( !(TMath::Min(met,trackMet) > 50.0)) continue;

    //************************************************************************************************
    // Background Process
    //************************************************************************************************
    Int_t bkgIndex = -1;
    if (dstype == SmurfTree::zz) bkgIndex = 0;
    else if (dstype == SmurfTree::wz) bkgIndex = 1;
    else if (dstype == SmurfTree::qqww || dstype == SmurfTree::ggww ) bkgIndex = 2;
    else if (dstype == SmurfTree::ttbar || dstype == SmurfTree::tw ) bkgIndex = 3;
    else if (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm ) bkgIndex = 4;
    else if (dstype == SmurfTree::dytt ) bkgIndex = 5;
    else if (dstype == SmurfTree::wjets) bkgIndex = 6;
    else bkgIndex = 7;


    //************************************************************************************************
    // weights and scale factors
    //************************************************************************************************
    Double_t myWeight = 0.0;
    Double_t weightFactor = 1.0;
    
    if(dstype == SmurfTree::data || 
       ((cuts & SmurfTree::Lep1LooseMuV1) == SmurfTree::Lep1LooseMuV1) || ((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) ||
       ((cuts & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep2LooseMuV1) || ((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) || 
       ((cuts & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep3LooseMuV1) || ((cuts & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep2LooseEleV4)  ){
      weightFactor *= fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, ((cuts & SmurfTree::Lep1LooseMuV1) == SmurfTree::Lep1LooseMuV1), ((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4));
      weightFactor *= fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, ((cuts & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep2LooseMuV1), ((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4));
      
      if(dstype == SmurfTree::data) {
        bkgIndex = 6;
        myWeight = weightFactor;
      }
      else if(TMath::Abs(lep1McId)*TMath::Abs(lep2McId) > 0 || dstype == SmurfTree::wgamma){
    	weightFactor *= PileupReweightFactor(nvtx);
        bkgIndex               = 6;
        myWeight	       = -1.0 * scale1fb*scaleFactorLum*weightFactor;
      }
     else {
        myWeight = 0.0;
      }
    }
    else {
      //****************************************************************************************
      //For MC Bkg Estimates
      //****************************************************************************************

      weightFactor *= PileupReweightFactor(nvtx);

      if( bkgIndex == 0 || bkgIndex == 1 || bkgIndex == 5){
        if((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
          weightFactor  *= leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
        if((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
          weightFactor *= leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
    	double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    								 fabs(lep2->eta()), lep2->pt(), 
								 TMath::Abs(lid1), TMath::Abs(lid2));
        weightFactor *= trigEff;
      }

      //ZZ DileptonPt dependant KFactor
      if (bkgIndex == 0) {
        weightFactor *= ZZKFactor->GetBinContent( ZZKFactor->GetXaxis()->FindFixBin(dilep->pt()));
      }
      
      //Scale Factor from e-mu control region
      if (bkgIndex == 2 || bkgIndex == 3) {
        if (njets==0)
          weightFactor *= TopAndWWScaleFactor[nJetsType];
      }
      
      myWeight 	     = scale1fb*scaleFactorLum*weightFactor;
      printf("pileup: %i : %i : %4.5f %4.5f\n", event, nvtx,  PileupReweightFactor(nvtx), myWeight ) ;
    }

//     myWeight = 1.0;

    //************************************************************************************************
    // Yields and histograms
    //************************************************************************************************
    nBgdAcc  += myWeight;
    nBgdEAcc += myWeight*myWeight;
    nBgdAccDecays[bkgIndex]  += myWeight;
    nBgdEAccDecays[bkgIndex] += myWeight*myWeight;

    assert(bkgIndex >= 0 && bkgIndex < nBkgChannels);

    if(useVar == 0) histoBkg[bkgIndex]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),     myWeight);
    else if(useVar == 1) histoBkg[bkgIndex]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);

    bgdMVA[0]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),     myWeight);
    bgdMVA[1]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001), myWeight);
    bgdMVADecays[0][bkgIndex]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),     myWeight);
    bgdMVADecays[1][bkgIndex]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001), myWeight);

    //**************************************************************************
    //Define HZZ Cuts
    //**************************************************************************
    double pxzll = dilep->px() + met*cos(metPhi);
    double pyzll = dilep->py() + met*sin(metPhi);
    double mtHZZ = TMath::Power(sqrt(dilep->pt()*dilep->pt() + dilep->M()*dilep->M())+
                         sqrt(met*met + dilep->M()*dilep->M()),2) - (pow(pxzll,2) + pow(pyzll,2));
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;
    
    bool passCutBased = ( (0==0)
                          && ( TMath::Min(met,trackMet) > cutMETLow[njets][mHIndex])
                          && (mtHZZ > cutMTLow[njets][mHIndex] && mtHZZ < cutMTHigh[njets][mHIndex]) );

    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passCutBased = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
                    jet1->eta()*jet2->eta() < 0 &&
		    (mH > 200 || dilep->mass() < 100.) &&
		    centrality == 1;
    }

    if(passCutBased == true) {
      nBgdCut  += myWeight;
      nBgdECut += myWeight*myWeight;
      nBgdCutDecays[bkgIndex]  += myWeight;
      nBgdECutDecays[bkgIndex] += myWeight*myWeight;
    }

    bool passMVACut = false;
    if     (useVar == 0 && nn   > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[bkgIndex]  += myWeight; nBgdEMVADecays[bkgIndex]  += myWeight*myWeight; passMVACut = true;}
    else if(useVar == 1 && bdtg > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[bkgIndex]  += myWeight; nBgdEMVADecays[bkgIndex]  += myWeight*myWeight; passMVACut = true;}
  }
  nBgdEAcc = sqrt(nBgdEAcc);
  nBgdECut = sqrt(nBgdECut);
  nBgdEMVA = sqrt(nBgdEMVA);
  printf("--- Finished Bgdnal loop\n");



  //**************************************************************************
  //Data Loop
  //**************************************************************************

  data->SetBranchAddress( "cuts"         , &cuts         );
  data->SetBranchAddress( "nvtx"         , &nvtx         );
  data->SetBranchAddress( "njets"        , &njets        );
  data->SetBranchAddress( "run"          , &run          );
  data->SetBranchAddress( "lumi"         , &lumi         );
  data->SetBranchAddress( "event"        , &event        );
  data->SetBranchAddress( "scale1fb"     , &scale1fb     );
  data->SetBranchAddress( "lep1"         , &lep1         );
  data->SetBranchAddress( "lep2"         , &lep2         );
  data->SetBranchAddress( "jet1"         , &jet1         );
  data->SetBranchAddress( "jet2"         , &jet2         );
  data->SetBranchAddress( "jet3"         , &jet3         );
  data->SetBranchAddress( "dPhi"         , &dPhi         );
  data->SetBranchAddress( "dR"           , &dR           );
  data->SetBranchAddress( "dilep"        , &dilep        );
  data->SetBranchAddress( "type"         , &type         );
  data->SetBranchAddress( "pmet"         , &pmet         );
  data->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  data->SetBranchAddress( "met"          , &met          );
  data->SetBranchAddress( "metPhi"       , &metPhi       );
  data->SetBranchAddress( "trackMet"     , &trackMet     );
  data->SetBranchAddress( "mt"           , &mt           );
  data->SetBranchAddress( "mt1"          , &mt1          );
  data->SetBranchAddress( "mt2"          , &mt2          );
  data->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  data->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  data->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  data->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  data->SetBranchAddress( "lq1"          , &lq1          );
  data->SetBranchAddress( "lq2"          , &lq2          );
  data->SetBranchAddress( "lid1"         , &lid1         );
  data->SetBranchAddress( "lid2"         , &lid2         );
  data->SetBranchAddress( "lid3"         , &lid3         );
  data->SetBranchAddress( "processId"    , &processId    );
  data->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  data->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  data->SetBranchAddress( "jet1Btag"	 , &jet1Btag	 );
  data->SetBranchAddress( "jet2Btag"	 , &jet2Btag	 );
  data->SetBranchAddress(Form("nn_hzz%i",mH)  , &nn      );
  data->SetBranchAddress(Form("bdtg_hzz%i",mH), &bdtg    );

  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    //apply goodrun list
    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

    //Remove this for now because the EPS ntuple doesn't have jet3...
//    //VBF Jet Selection (Central Jet Veto)
//     unsigned int Njet3 = njets;
//     if(nJetsType == 2){
//       if(jet3->pt() <= 30)					         Njet3 = 2;
//       else if(jet3->pt() > 30 && (
//         (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
//         (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))   Njet3 = 0;
//       else								 Njet3 = 2;
//       if(njets < 2 || njets > 3) Njet3 = 0;
//     }


    //************************************************************************************************
    // PreSelection
    //************************************************************************************************

    if( njets != nJetsType ) continue; // select n-jet type events
    if( !(lep1->pt() > 30 && lep2->pt() > 20 )) continue;
    if( !((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection )) continue;
    if( !((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection )) continue;
    if( !(type == SmurfTree::mm || type == SmurfTree::ee )) continue;
    if( !(fabs(dilep->mass() - 91.1876) < 15 )) continue;
    if( !(lq1*lq2 < 0)) continue;
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 
    if( !(TMath::Min(met,trackMet) > 50.0)) continue;

    double myWeight = 1.0;

    //************************************************************************************************
    // Yields and histograms
    //************************************************************************************************
    nDatAcc = nDatAcc + myWeight;
         if(useVar == 0) histoData->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),   myWeight);
    else if(useVar == 1) histoData->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001), myWeight);

    datMVA[0]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[0]-0.001),minHis[0]+0.001),   myWeight);
    datMVA[1]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[1]-0.001),minHis[1]+0.001), myWeight);


    //**************************************************************************
    //Define HZZ Cuts
    //**************************************************************************
    double pxzll = dilep->px() + met*cos(metPhi);
    double pyzll = dilep->py() + met*sin(metPhi);
    double mtHZZ = TMath::Power(sqrt(dilep->pt()*dilep->pt() + dilep->M()*dilep->M())+
                         sqrt(met*met + dilep->M()*dilep->M()),2) - (pow(pxzll,2) + pow(pyzll,2));
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;
    
    bool passCutBased = ( (0==0)
                          && ( TMath::Min(met,trackMet) > cutMETLow[njets][mHIndex])
                          && (mtHZZ > cutMTLow[njets][mHIndex] && mtHZZ < cutMTHigh[njets][mHIndex]) );
    
    if(nJetsType == 2) {
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passCutBased = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
                    jet1->eta()*jet2->eta() < 0 &&
		    (mH > 200 || dilep->mass()< 100.) &&
		    centrality == 1;
    }


    if(passCutBased == true) {
      nDatCut += myWeight;
    }

    bool passMVACut = false;
         if(useVar == 0 && nn         > useCut) {nDatMVA += myWeight; passMVACut = true;}
    else if(useVar == 1 && bdtg       > useCut) {nDatMVA += myWeight; passMVACut = true;}
  }


  //**************************************************************************
  //Summary
  //**************************************************************************
  printf("--- Finished Data loop\n");
  printf("---\tacceptedSPresel  %8.3f +/- %8.3f events\n",nSigAcc[0],nSigEAcc[0]);
  printf("---\tacceptedSCuts    %8.3f +/- %8.3f events\n",nSigCut[0],nSigECut[0]);
  printf("---\tacceptedSANNCuts %8.3f +/- %8.3f events\n",nSigMVA[0],nSigEMVA[0]);
  printf("\n");
  printf("---\tacceptedBPresel  %8.3f +/- %8.3f events\n",nBgdAcc,nBgdEAcc);
  printf("---\tacceptedBCuts    %8.3f +/- %8.3f events\n",nBgdCut,nBgdECut);
  printf("---\tacceptedBANNCuts %8.3f +/- %8.3f events\n",nBgdMVA,nBgdEMVA);
  printf("\n");
  printf("---\tacceptedDPresel  %8.3f events\n",nDatAcc);
  printf("---\tacceptedDCuts    %8.3f events\n",nDatCut);
  printf("---\tacceptedDANNCuts %8.3f events\n",nDatMVA);
  printf("\n");

  for(int i=0; i<nBkgChannels; i++) {
    if(nBgdAccDecays[i] > 0.0) nBgdEAccDecays[i] = sqrt(nBgdEAccDecays[i])/nBgdAccDecays[i];
    if(nBgdCutDecays[i] > 0.0) nBgdECutDecays[i] = sqrt(nBgdECutDecays[i])/nBgdCutDecays[i];
    if(nBgdMVADecays[i] > 0.0) nBgdEMVADecays[i] = sqrt(nBgdEMVADecays[i])/nBgdMVADecays[i];
  }
  printf("            HZZ      ZZ     WZ     WW       top      DYee/mm    DYtt    wjets     other\n");
  printf("CLsAcc : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigAcc[0],nBgdAccDecays[0],nBgdAccDecays[1],nBgdAccDecays[2],nBgdAccDecays[3],nBgdAccDecays[4],nBgdAccDecays[5],nBgdAccDecays[6],nBgdAccDecays[7]);
  printf("CLsCut : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigCut[0],nBgdCutDecays[0],nBgdCutDecays[1],nBgdCutDecays[2],nBgdCutDecays[3],nBgdCutDecays[4],nBgdCutDecays[5],nBgdCutDecays[6],nBgdCutDecays[7]);
  printf("CLsMVA : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigMVA[0],nBgdMVADecays[0],nBgdMVADecays[1],nBgdMVADecays[2],nBgdMVADecays[3],nBgdMVADecays[4],nBgdMVADecays[5],nBgdMVADecays[6],nBgdMVADecays[7]);
  printf("CLsEAcc: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigEAcc[0]/nSigAcc[0],nBgdEAccDecays[0],nBgdEAccDecays[1],nBgdEAccDecays[2],nBgdEAccDecays[3],nBgdEAccDecays[4],nBgdEAccDecays[5],nBgdEAccDecays[6],nBgdEAccDecays[7]);
  printf("CLsECut: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigECut[0]/nSigCut[0],nBgdECutDecays[0],nBgdECutDecays[1],nBgdECutDecays[2],nBgdECutDecays[3],nBgdECutDecays[4],nBgdECutDecays[5],nBgdECutDecays[6],nBgdECutDecays[7]);
  printf("CLsEMVA: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigEMVA[0]/nSigMVA[0],nBgdEMVADecays[0],nBgdEMVADecays[1],nBgdEMVADecays[2],nBgdEMVADecays[3],nBgdEMVADecays[4],nBgdEMVADecays[5],nBgdEMVADecays[6],nBgdEMVADecays[7]);

  //----------------------------------------------------------------------------
  // Drawing part - c1
  //----------------------------------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","c1",100,100,700,800);
  c1->Divide(1,1);

  c1->cd(1); plotHistsInPad(sigMVA[0][0], bgdMVA[0]);
  c1->cd(2); plotHistsInPad(sigMVA[1][0], bgdMVA[1]);

  c1->SaveAs(c1Name);

  //----------------------------------------------------------------------------
  // Drawing part - c4
  //----------------------------------------------------------------------------
  TCanvas* c3 = new TCanvas("c3","c3",100,100,700,800);
  c3->SetGrid(1,1);
  TH2D* zone4 = new TH2D("zone4","discriminant vs. significance; Output; Significance",
			 1, -1, 1, 1, 1.e-3, 2);

  c3->SetLeftMargin(1.2*c3->GetLeftMargin());
  zone4->DrawCopy();

  TGraphErrors* g3MLP         = makeSignificanceCurve(sigMVA[0][0], bgdMVA[0], datMVA[0],"g3NN"	    );  
  TGraphErrors* g3BDTG        = makeSignificanceCurve(sigMVA[1][0], bgdMVA[1], datMVA[1],"g3BDTG"   );  

  TFile *f=new TFile("MVASig.root","UPDATE");
  f->WriteTObject(g3MLP, "NN", "WriteDelete");
  f->WriteTObject(g3BDTG, "BDT", "WriteDelete");
  f->Close();

  setGraph(g3BDTG      ,1, 25);
  c3->SaveAs(("plot_"+string(outTag.Data())+"_BDTG_Sig.gif").c_str());
  setGraph(g3MLP       ,0, 23);
  c3->SaveAs(("plot_"+string(outTag.Data())+"_NN_Sig.gif").c_str());

  TCanvas* c4 = new TCanvas("c4","c4",100,100,700,800);
  c4->cd();
  gPad->SetGrid(1,1);
  TGraphErrors* g4MLP         = makeGraphFromHists(sigMVA[0][0], bgdMVA[0], "g3NN"        );  
  TGraphErrors* g4BDTG        = makeGraphFromHists(sigMVA[1][0], bgdMVA[1], "g3BDTG"      );  
  setGraph(g4MLP       ,0, 23);
  setGraph(g4BDTG      ,1, 25);
  g4MLP       ->Draw("APXl");
  g4BDTG      ->Draw("PXl");
  TAxis *axisX=g4MLP->GetXaxis();
  TAxis *axisY=g4MLP->GetYaxis();
  axisX->SetTitle("Signal Efficiency");
  axisX->Draw();
  axisY->SetTitle("Background Efficiency");
  axisY->Draw();
  TLegend* leg4 = new TLegend(0.2, 0.7, 0.5, 0.9);
  leg4->SetFillColor(0);
  leg4->AddEntry(g4MLP       ,"NN"        ,"L");
  leg4->AddEntry(g4BDTG      ,"BTDG"      ,"L");
  leg4->Draw("same");
  c4->SaveAs(c2Name);

  char output[200];
  sprintf(output,"histo_tmva_%s_%dj_chan%d_mh%d.root",outTag.Data(),nJetsType,wwDecay,mH);     
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
  for(int i=0; i<nMVAVars; i++) {
    for(int j=0; j<nSignalChannels; j++){
      sigMVA[i][j]->Write();
    }
    bgdMVA[i]->Write();
    datMVA[i]->Write();
    for(int j=0; j<nBkgChannels; j++) {
      bgdMVADecays[i][j]->Write();
    }
    histoSignal->Write();
    for(UInt_t k=0; k<nBkgChannels; ++k) {
      histoBkg[k]->Write();
    } 
    histoData->Write();
  }
  outFilePlotsNote->Close();





  if(fillInfoNote == true){
    for (UInt_t k=1; k<nSignalChannels; ++k) {
      cout << sigMVA[useVar][k]->GetSumOfWeights() << " ";
    }
    for (UInt_t k=1; k<nBkgChannels; ++k) {
      cout << bgdMVADecays[useVar][k]->GetSumOfWeights() << " ";
    }
    cout << endl;


    TH1D *histo_ggH    = (TH1D*) sigMVA[useVar][1]->Clone("histo_ggH");
    TH1D *histo_qqH    = (TH1D*) sigMVA[useVar][2]->Clone("histo_qqH");
    TH1D *histo_VH     = (TH1D*) sigMVA[useVar][3]->Clone("histo_WH");
    TH1D *histo_ttH    = (TH1D*) sigMVA[useVar][5]->Clone("histo_ttH");

    TH1D *histo_Data   = (TH1D*) datMVA[useVar]->Clone("histo_Data");
    TH1D *histo_ZZ     = (TH1D*) bgdMVADecays[useVar][0]->Clone("histo_ZZ");
    TH1D *histo_WZ     = (TH1D*) bgdMVADecays[useVar][1]->Clone("histo_WZ");
    TH1D *histo_WW     = (TH1D*) bgdMVADecays[useVar][2]->Clone("histo_WW");
    TH1D *histo_Top    = (TH1D*) bgdMVADecays[useVar][3]->Clone("histo_Top");
    TH1D *histo_Zjets  = (TH1D*) bgdMVADecays[useVar][4]->Clone("histo_Zjets");
    TH1D *histo_DYtt   = (TH1D*) bgdMVADecays[useVar][5]->Clone("histo_DYtt");
    TH1D *histo_Wjets  = (TH1D*) bgdMVADecays[useVar][6]->Clone("histo_Wjets");
    histo_ggH    ->Rebin(10);
    histo_qqH    ->Rebin(10);
    histo_VH     ->Rebin(10);
    histo_ttH    ->Rebin(10);
    histo_Data   ->Rebin(10);
    histo_ZZ     ->Rebin(10);
    histo_WZ     ->Rebin(10);
    histo_WW	 ->Rebin(10);
    histo_Top	 ->Rebin(10);
    histo_Zjets  ->Rebin(10);
    histo_DYtt   ->Rebin(10);
    histo_Wjets  ->Rebin(10);
    char outputLimits[200];
    sprintf(outputLimits,"Smurf/HZZ/output/histo_limits_%s_%dj_chan%d_mh%d.root",outTag.Data(),nJetsType,wwDecay,mH);     
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
      histo_ggH    ->Write();
      histo_qqH    ->Write();
      histo_VH     ->Write();
      histo_ttH    ->Write();
      histo_Data   ->Write();
      histo_ZZ     ->Write();
      histo_WZ     ->Write();
      histo_WW	   ->Write();
      histo_Top	   ->Write();
      histo_Zjets  ->Write();
      histo_DYtt   ->Write();
      histo_Wjets  ->Write();
    outFileLimits->Close();

    double jeteff_E 	     = 1.02;
    double ZXS_E  	     = 1.30;
    double XS_QCDscale_ggH[3];
    double UEPS  = HiggsSignalPSUESystematics(mH,nJetsType);
    XS_QCDscale_ggH[0] = HiggsSignalQCDScaleKappa("QCDscale_ggH",mH,nJetsType);
    XS_QCDscale_ggH[1] = HiggsSignalQCDScaleKappa("QCDscale_ggH1in",mH,nJetsType);
    XS_QCDscale_ggH[2] = HiggsSignalQCDScaleKappa("QCDscale_ggH2in",mH,nJetsType);
    if     (nJetsType == 1) {
      jeteff_E  	= 1.05;
      ZXS_E     	= 1.30;
    }
    else if(nJetsType == 2) {
      jeteff_E  	= 1.10;
       ZXS_E     	= 1.30;
    }
    for(int i=0; i<nBkgChannels; i++) if(nBgdAccDecays[i] < 0) nBgdAccDecays[i] = 0.0;
    for(int i=0; i<nBkgChannels; i++) if(nBgdCutDecays[i] < 0) nBgdCutDecays[i] = 0.0;
    for(int i=0; i<nBkgChannels; i++) if(nBgdMVADecays[i] < 0) nBgdMVADecays[i] = 0.0;
    for(int i=0; i<nSignalChannels; i++) if(nSigAcc[i] <= 0) nSigAcc[i] = 0.001;
    for(int i=0; i<nSignalChannels; i++) if(nSigCut[i] <= 0) nSigCut[i] = 0.001;
    for(int i=0; i<nSignalChannels; i++) if(nSigMVA[i] <= 0) nSigMVA[i] = 0.001;
    double yieldE[12],yield[12];
    int nData;
    int nTotalBins = 1; // histo_VH->GetNbinsX();
    //********SHAPE*******************
    char outputLimitsShape[200];
    for(int i=1; i<=nTotalBins; i++){
      if(nTotalBins != 1){
        if(histo_ttH   ->GetBinContent(i) > 0) { yieldE[0] = histo_ttH   ->GetBinError(i)/histo_ttH   ->GetBinContent(i);} else {yieldE[0] = 0.0;}
        if(histo_VH    ->GetBinContent(i) > 0) { yieldE[1] = histo_VH    ->GetBinError(i)/histo_VH    ->GetBinContent(i);} else {yieldE[1] = 0.0;}
        if(histo_qqH   ->GetBinContent(i) > 0) { yieldE[2] = histo_qqH   ->GetBinError(i)/histo_qqH   ->GetBinContent(i);} else {yieldE[2] = 0.0;}
        if(histo_ggH   ->GetBinContent(i) > 0) { yieldE[3] = histo_ggH   ->GetBinError(i)/histo_ggH   ->GetBinContent(i);} else {yieldE[3] = 0.0;}
        if(histo_ZZ    ->GetBinContent(i) > 0) { yieldE[4] = histo_ZZ    ->GetBinError(i)/histo_ZZ    ->GetBinContent(i);} else {yieldE[4] = 0.0;}
        if(histo_WZ    ->GetBinContent(i) > 0) { yieldE[5] = histo_WZ    ->GetBinError(i)/histo_WZ    ->GetBinContent(i);} else {yieldE[5] = 0.0;}
        if(histo_WW    ->GetBinContent(i) > 0) { yieldE[6] = histo_WW    ->GetBinError(i)/histo_WW    ->GetBinContent(i);} else {yieldE[6] = 0.0;}
        if(histo_Top   ->GetBinContent(i) > 0) { yieldE[7] = histo_Top   ->GetBinError(i)/histo_Top   ->GetBinContent(i);} else {yieldE[7] = 0.0;}
        if(histo_Zjets ->GetBinContent(i) > 0) { yieldE[8] = histo_Zjets ->GetBinError(i)/histo_Zjets ->GetBinContent(i);} else {yieldE[8] = 0.0;}
        if(histo_DYtt  ->GetBinContent(i) > 0) { yieldE[9] = histo_DYtt  ->GetBinError(i)/histo_DYtt  ->GetBinContent(i);} else {yieldE[9] = 0.0;}
        if(histo_Wjets ->GetBinContent(i) > 0) { yieldE[10]= histo_Wjets ->GetBinError(i)/histo_Wjets ->GetBinContent(i);} else {yieldE[10]= 0.0;}

	yield[0] = histo_ttH   ->GetBinContent(i);
	yield[1] = histo_VH    ->GetBinContent(i);
	yield[2] = histo_qqH   ->GetBinContent(i);
	yield[3] = histo_ggH   ->GetBinContent(i);
	yield[4] = histo_ZZ    ->GetBinContent(i);
	yield[5] = histo_WZ    ->GetBinContent(i);
	yield[6] = histo_WW    ->GetBinContent(i);
	yield[7] = histo_Top   ->GetBinContent(i);
	yield[8] = histo_Zjets ->GetBinContent(i);
	yield[9] = histo_DYtt  ->GetBinContent(i);
	yield[10]= histo_Wjets ->GetBinContent(i);
	nData    = (int)histo_Data->GetBinContent(i);
      }
      else {
 	if(nSigAcc[4] > 0) {yieldE[0] = nSigEAcc[4]/nSigAcc[4];} else {yieldE[0] = 0.0;}
 	if(nSigAcc[3] > 0) {yieldE[1] = nSigEAcc[3]/nSigAcc[3];} else {yieldE[1] = 0.0;}
 	if(nSigAcc[2] > 0) {yieldE[2] = nSigEAcc[2]/nSigAcc[2];} else {yieldE[2] = 0.0;}
 	if(nSigAcc[1] > 0) {yieldE[3] = nSigEAcc[1]/nSigAcc[1];} else {yieldE[3] = 0.0;}
 	yieldE[4] = nBgdEAccDecays[0];
 	yieldE[5] = nBgdEAccDecays[1];
 	yieldE[6] = nBgdEAccDecays[2];
 	yieldE[7] = nBgdEAccDecays[3];
 	yieldE[8] = nBgdEAccDecays[4];
 	yieldE[9] = nBgdEAccDecays[5];
 	yieldE[10]= nBgdEAccDecays[6];
	yield[0]  = nSigAcc[4];
	yield[1]  = nSigAcc[3];
	yield[2]  = nSigAcc[2];
	yield[3]  = nSigAcc[1];
	yield[4]  = nBgdAccDecays[0];
	yield[5]  = nBgdAccDecays[1];
	yield[6]  = nBgdAccDecays[2];
	yield[7]  = nBgdAccDecays[3];
	yield[8]  = nBgdAccDecays[4];
	yield[9]  = nBgdAccDecays[5];;
	yield[10] = nBgdAccDecays[6];
	nData     = (int)nDatAcc;
      }
      if(nTotalBins != 1){
        sprintf(outputLimitsShape,"Smurf/HZZ/output/histo_limits_%s_%dj_chan%d_mh%d_shape_bin%d.txt",outTag.Data(),nJetsType,wwDecay,mH,i);     
      }
      else {
        sprintf(outputLimitsShape,"Smurf/HZZ/output/histo_limits_%s_%dj_chan%d_mh%d_shape.txt",outTag.Data(),nJetsType,wwDecay,mH);     
      }
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",nData);
      if(nTotalBins == 1){
        newcardShape << Form("shapes *	      * %s  histo_$PROCESS\n",outputLimits);
        newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
      }
      newcardShape << Form("bin 1 1 1 1 1 1 1 1 1 1 \n");
      newcardShape << Form("process VH qqH ggH ZZ WZ WW Top Zjets DYtt Wjets\n");
      newcardShape << Form("process -2 -1 0 1 2 3 4 5 6 7\n");
      newcardShape << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",yield[1],yield[2],yield[3],yield[4],yield[5],yield[6],yield[7],yield[8],yield[9],TMath::Max((double)yield[10],0.0));
      newcardShape << Form("lumi		           lnN 1.100 1.100 1.100 1.000 1.000 1.100 1.000 1.000 1.000 1.100\n"); 		      
      newcardShape << Form("CMS_eff_m                      lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 		      
      newcardShape << Form("CMS_trigger_m	           lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 		      
      newcardShape << Form("CMS_eff_e                      lnN 1.025 1.025 1.025 1.025 1.025 1.025 1.000 1.000 1.000 1.025\n"); 		      
      newcardShape << Form("CMS_trigger_e                  lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 		      
      newcardShape << Form("CMS_scale_m                    lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 		      
      newcardShape << Form("CMS_scale_e  	           lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n");
      newcardShape << Form("CMS_hww_met_resolution         lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n");
      newcardShape << Form("CMS_scale_j                    lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f 1.000 1.000 1.000 %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);        
      newcardShape << Form("CMS_hww_fake_em       	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.360 1.000\n");
      newcardShape << Form("UEPS 		           lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",UEPS);
      newcardShape << Form("pdf_gg		           lnN 1.000 1.000 1.100 1.000 1.040 1.000 1.000 1.000 1.000 1.000\n");
      newcardShape << Form("pdf_qqbar                      lnN 1.050 1.050 1.000 1.040 1.000 1.040 1.000 1.000 1.000 1.040\n");
      newcardShape << Form("QCDscale_ggH	           lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[0]);  
      newcardShape << Form("QCDscale_ggH1in	           lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[1]);  
      newcardShape << Form("QCDscale_ggH2in	           lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[2]);  
      newcardShape << Form("QCDscale_qqH	           lnN 1.000 1.010 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
      newcardShape << Form("QCDscale_VH 	           lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n"); 	      
      newcardShape << Form("QCDscale_VV 	           lnN 1.000 1.000 1.000 1.000 1.000 1.040 1.000 1.000 1.000 1.000\n");
      newcardShape << Form("QCDscale_ggH_EXTRAP            lnN 1.000 1.000 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
      newcardShape << Form("QCDscale_qqH_EXTRAP            lnN 1.000 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
      newcardShape << Form("QCDscale_VH_EXTRAP             lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
      newcardShape << Form("CMS_hww_%1dj_ttbar             lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,TopAndWWScaleFactorKappa[nJetsType]);	     
      newcardShape << Form("CMS_hww_%1dj_Z                 lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,ZXS_E); 		     
      newcardShape << Form("CMS_hww_%1dj_WW                lnN 1.000 1.000 1.000 %5.3f %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,TopAndWWScaleFactorKappa[nJetsType],TopAndWWScaleFactorKappa[nJetsType]); 		     
      newcardShape << Form("CMS_hww_stat_%1dj_VH_bin%d     lnN %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,i,yieldE[1]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_qqH_bin%d    lnN 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,i,yieldE[2]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_ggH_bin%d    lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,i,yieldE[3]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_ZZ_bin%d     lnN 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,i,yieldE[4]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_WZ_bin%d     lnN 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,i,yieldE[5]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_WW_bin%d     lnN 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000\n",nJetsType,i,yieldE[6]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_ttbar_bin%d  lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,i,yieldE[7]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_Z_bin%d      lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,i,yieldE[8]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_DYtt_bin%d   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000\n",nJetsType,i,yieldE[9]+1.0);
      newcardShape << Form("CMS_hww_stat_%1dj_Wjets_bin%d lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f\n",nJetsType,i,yieldE[10]+1.0);
      newcardShape.close();
    }

    //********CUT*******************
    char outputLimitsCut[200];
    sprintf(outputLimitsCut,"Smurf/HZZ/output/histo_limits_%s_%dj_chan%d_mh%d_cut.txt",outTag.Data(),nJetsType,wwDecay,mH);     
    ofstream newcardCut;
    newcardCut.open(outputLimitsCut);
    newcardCut << Form("imax 1 number of channels\n");
    newcardCut << Form("jmax * number of background\n");
    newcardCut << Form("kmax * number of nuisance parameters\n");

    //Actual Data
    newcardCut << Form("Observation %d\n",(int)nDatCut);
    //Expected from MC
    //newcardCut << Form("Observation %d\n",(int)(nBgdCutDecays[0]+nBgdCutDecays[1]+nBgdCutDecays[2]+nBgdCutDecays[3]+nBgdCutDecays[4]+nBgdCutDecays[5]+nBgdCutDecays[6]));

    newcardCut << Form("bin 1 1 1 1 1 1 1 1 1 1 \n");
    newcardCut << Form("process VH qqH ggH ZZ WZ WW Top Zjets DYtt Wjets\n");
    newcardCut << Form("process -2 -1 0 1 2 3 4 5 6 7\n");
    newcardCut << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigCut[3],nSigCut[2],nSigCut[1],nBgdCutDecays[0],nBgdCutDecays[1],nBgdCutDecays[2],nBgdCutDecays[3],nBgdCutDecays[4],nBgdCutDecays[5],TMath::Max((double)nBgdCutDecays[6],0.0));
    newcardCut << Form("lumi                  	   lnN 1.100 1.100 1.100 1.000 1.000 1.100 1.000 1.000 1.000 1.100\n"); 			
    newcardCut << Form("CMS_eff_m             	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardCut << Form("CMS_trigger_m         	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardCut << Form("CMS_eff_e             	   lnN 1.025 1.025 1.025 1.025 1.025 1.025 1.000 1.000 1.000 1.025\n"); 			
    newcardCut << Form("CMS_trigger_e         	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardCut << Form("CMS_scale_m           	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardCut << Form("CMS_scale_e           	   lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n"); 			
    newcardCut << Form("CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n"); 			
    newcardCut << Form("CMS_scale_j           	   lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f 1.000 1.000 1.000 %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E); 	 
    newcardCut << Form("CMS_hww_fake_em       	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.360 1.000\n");
    newcardCut << Form("UEPS 		           lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",UEPS);
    newcardCut << Form("pdf_gg                	   lnN 1.000 1.000 1.100 1.000 1.040 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("pdf_qqbar             	   lnN 1.050 1.050 1.000 1.040 1.000 1.040 1.000 1.000 1.000 1.040\n");
    newcardCut << Form("QCDscale_ggH          	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[0]);  
    newcardCut << Form("QCDscale_ggH1in       	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[1]);  
    newcardCut << Form("QCDscale_ggH2in       	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[2]);  
    newcardCut << Form("QCDscale_qqH          	   lnN 1.000 1.010 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_VH           	   lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n"); 		
    newcardCut << Form("QCDscale_VV           	   lnN 1.000 1.000 1.000 1.000 1.000 1.040 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_ggH_EXTRAP   	   lnN 1.000 1.000 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_qqH_EXTRAP   	   lnN 1.000 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_VH_EXTRAP    	   lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("CMS_hww_%1dj_ttbar	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,TopAndWWScaleFactorKappa[nJetsType]);    
    newcardCut << Form("CMS_hww_%1dj_Z  	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,ZXS_E); 	     
    newcardCut << Form("CMS_hww_%1dj_WW 	   lnN 1.000 1.000 1.000 %5.3f %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,TopAndWWScaleFactorKappa[nJetsType],TopAndWWScaleFactorKappa[nJetsType]); 		  
    newcardCut << Form("CMS_hww_stat_%1dj_VH	   lnN %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nSigECut[3]/TMath::Max((double)nSigCut[3],0.00001)+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_qqH	   lnN 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nSigECut[2]/TMath::Max((double)nSigCut[2],0.00001)+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_ggH	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nSigECut[1]/TMath::Max((double)nSigCut[1],0.00001)+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_ZZ	   lnN 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nBgdECutDecays[0]+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_WZ       lnN 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nBgdECutDecays[1]+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_WW	   lnN 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000\n",nJetsType,nBgdECutDecays[2]+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_ttbar    lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,nBgdECutDecays[3]+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_Z	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,nBgdECutDecays[4]+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_DYtt     lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000\n",nJetsType,nBgdECutDecays[5]+1.0);
    newcardCut << Form("CMS_hww_stat_%1dj_Wjets    lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f\n",nJetsType,nBgdECutDecays[6]+1.0);
    newcardCut.close();

    //***********MVA***************
    char outputLimitsMVA[200];
    sprintf(outputLimitsMVA,"Smurf/HZZ/output/histo_limits_%s_%dj_chan%d_mh%d_mva.txt",outTag.Data(),nJetsType,wwDecay,mH);     
    ofstream newcardMVA;
    newcardMVA.open(outputLimitsMVA);
    newcardMVA << Form("imax 1 number of channels\n");
    newcardMVA << Form("jmax * number of background\n");
    newcardMVA << Form("kmax * number of nuisance parameters\n");

    //Actual Data
    newcardMVA << Form("Observation %d\n",(int)nDatMVA);
    //Expected from MC
    //newcardMVA << Form("Observation %d\n",(int)(nBgdMVADecays[0]+nBgdMVADecays[1]+nBgdMVADecays[2]+nBgdMVADecays[3]+nBgdMVADecays[4]+nBgdMVADecays[5]+nBgdMVADecays[6]));

    newcardMVA << Form("bin 1 1 1 1 1 1 1 1 1 1 \n");
    newcardMVA << Form("process VH qqH ggH ZZ WZ WW Top Zjets DYtt Wjets \n");
    newcardMVA << Form("process -2 -1 0 1 2 3 4 5 6 7\n");
    newcardMVA << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigMVA[3],nSigMVA[2],nSigMVA[1],nBgdMVADecays[0],nBgdMVADecays[1],nBgdMVADecays[2],nBgdMVADecays[3],nBgdMVADecays[4],nBgdMVADecays[5],TMath::Max((double)nBgdMVADecays[6],0.0));
    newcardMVA << Form("lumi                  	   lnN 1.100 1.100 1.100 1.000 1.000 1.100 1.000 1.000 1.000 1.100\n");
    newcardMVA << Form("CMS_eff_m             	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardMVA << Form("CMS_trigger_m         	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardMVA << Form("CMS_eff_e             	   lnN 1.025 1.025 1.025 1.025 1.025 1.025 1.000 1.000 1.000 1.025\n"); 			
    newcardMVA << Form("CMS_trigger_e         	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardMVA << Form("CMS_scale_m           	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n"); 			
    newcardMVA << Form("CMS_scale_e           	   lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n"); 			
    newcardMVA << Form("CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n"); 			
    newcardMVA << Form("CMS_scale_j           	   lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f 1.000 1.000 1.000 %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E); 	  
    newcardMVA << Form("CMS_hww_fake_em       	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.360 1.000\n");
    newcardMVA << Form("UEPS 		           lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",UEPS);
    newcardMVA << Form("pdf_gg                	   lnN 1.000 1.000 1.100 1.000 1.040 1.000 1.000 1.000 1.000 1.000\n");
    newcardMVA << Form("pdf_qqbar             	   lnN 1.050 1.050 1.000 1.040 1.000 1.040 1.000 1.000 1.000 1.040\n");
    newcardMVA << Form("QCDscale_ggH          	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[0]);  
    newcardMVA << Form("QCDscale_ggH1in       	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[1]);  
    newcardMVA << Form("QCDscale_ggH2in       	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[2]);  
    newcardMVA << Form("QCDscale_qqH          	   lnN 1.000 1.010 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardMVA << Form("QCDscale_VH           	   lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n"); 		
    newcardMVA << Form("QCDscale_VV           	   lnN 1.000 1.000 1.000 1.000 1.000 1.040 1.000 1.000 1.000 1.000\n");
    newcardMVA << Form("QCDscale_ggH_EXTRAP   	   lnN 1.000 1.000 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardMVA << Form("QCDscale_qqH_EXTRAP   	   lnN 1.000 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardMVA << Form("QCDscale_VH_EXTRAP    	   lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardMVA << Form("CMS_hww_%1dj_ttbar	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,TopAndWWScaleFactorKappa[nJetsType]);    
    newcardMVA << Form("CMS_hww_%1dj_Z  	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,ZXS_E); 	     
    newcardMVA << Form("CMS_hww_%1dj_WW 	   lnN 1.000 1.000 1.000 %5.3f %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,TopAndWWScaleFactorKappa[nJetsType],TopAndWWScaleFactorKappa[nJetsType]); 		  
    newcardMVA << Form("CMS_hww_stat_%1dj_VH	   lnN %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nSigEMVA[3]/TMath::Max((double)nSigMVA[3],0.00001)+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_qqH	   lnN 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nSigEMVA[2]/TMath::Max((double)nSigMVA[2],0.00001)+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_ggH	   lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nSigEMVA[1]/TMath::Max((double)nSigMVA[1],0.00001)+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_ZZ	   lnN 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nBgdEMVADecays[0]+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_WZ       lnN 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,nBgdEMVADecays[1]+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_WW	   lnN 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000\n",nJetsType,nBgdEMVADecays[2]+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_ttbar    lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,nBgdEMVADecays[3]+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_Z	   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,nBgdEMVADecays[4]+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_DYtt     lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000\n",nJetsType,nBgdEMVADecays[5]+1.0);
    newcardMVA << Form("CMS_hww_stat_%1dj_Wjets    lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f\n",nJetsType,nBgdEMVADecays[6]+1.0);
    newcardMVA.close();
  }

}

//------------------------------------------------------------------------------
// makeSignificanceCurve
//------------------------------------------------------------------------------
TGraphErrors* makeSignificanceCurve(TH1D* sig, TH1D* bgd, TH1D* dat, const char* name)
{
  double xs [1000];
  double ys [1000];
  double dys[1000];
  double dxs[1000];

  int n = 0;
	
  double theValue[4] = {0., 0., 0., 0.};

  for (int bin=1; bin<=sig->GetNbinsX(); ++bin) {
		
    double s =      sig->Integral(bin, sig->GetNbinsX());
    double b =      bgd->Integral(bin, sig->GetNbinsX());
    int    d = (int)dat->Integral(bin, sig->GetNbinsX());

    if (b > 0.01) ys[n] = s/(sqrt(b+sigmaB*sigmaB*b*b));
    else          ys[n] = 0.0;

    xs[n] = sig->GetBinCenter(bin);

    if(verboseLevel == 1){
      if (b > 0.01 && s > 0.01) printf("%15s Bin-> Sig= %8.3f - S= %8.3f - B= %8.3f - S/B= %8.3f - Bin= %8.3f | D = %d\n",
        	        name,
        	        ys[n],s,
        	        b,s/b,xs[n],
			d);
    }

    dxs[n] = 1.e-6;
    dys[n] = 1.e-6;
	   
    if (ys[n] > theValue[0]) {
      theValue[0] = ys[n];
      theValue[1] = s;
      theValue[2] = b;
      theValue[3] = xs[n];	     
    }
    ++n;
  }

  if(theValue[2] <= 0) theValue[2] = 1000000.;
  printf("%15s -> Sig= %8.3f - S= %8.3f - B= %8.3f - S/B= %8.3f - Bin= %8.3f\n",
       name,
       theValue[0],theValue[1],
       theValue[2],theValue[1]/theValue[2],theValue[3]);
  TGraphErrors* g = new TGraphErrors(n, xs, ys, dxs, dys);
  g->SetName(name);
  return g;
}

//------------------------------------------------------------------------------
// setGraph
//------------------------------------------------------------------------------
void setGraph(TGraphErrors* g, int color, int marker)
{
  g->SetLineColor  ( color);
  g->SetLineWidth  (     2);
  g->SetMarkerColor( color);
  g->SetMarkerSize (   0.5);
  g->SetMarkerStyle(marker);

  g->Draw("L");
}

//------------------------------------------------------------------------------
// plotHistsInPad
//------------------------------------------------------------------------------
void plotHistsInPad(TH1D* h1, TH1D* h2)
{
  gPad->SetLogy();
  setHist(h1, 4, 20);
  setHist(h2, 2, 24);
  setPair(h1, h2);
  h1->DrawCopy("pe");
  h2->DrawCopy("pe same");

  TLegend* lg = new TLegend(0.65, 0.65, 0.93, 0.90);
  lg->SetFillColor(0);
  lg->AddEntry(h1," signal");
  lg->AddEntry(h2," background");
  lg->Draw("same");
}

//------------------------------------------------------------------------------
// setPair
//------------------------------------------------------------------------------
void setPair(TH1D* h1, TH1D* h2)
{
  double theMax = 10.*TMath::Max(h1->GetMaximum(),h2->GetMaximum());
  h1->SetMaximum(theMax);
  h2->SetMaximum(theMax);
}

//------------------------------------------------------------------------------
// setHist
//------------------------------------------------------------------------------
void setHist(TH1D* h, int color, int style)
{
  h->SetMarkerColor(color        );
  h->SetMarkerStyle(style        );
  h->SetLineColor  (color        );
  h->SetMarkerSize (0.5          );
  h->SetXTitle     (h->GetTitle());
}

//------------------------------------------------------------------------------
// makeGraphFromHists
//------------------------------------------------------------------------------
TGraphErrors* makeGraphFromHists(TH1D* hsig, TH1D* hbgd, const char* name)
{
  const int nbins = hsig->GetNbinsX();
  double xs[1000], ys[1000], dxs[1000], dys[1000];

  int i=0;
  for (int bin=1; bin <= nbins; ++bin) {
    double bgds = hbgd->Integral(bin, nbins);
    double sigs = hsig->Integral(bin, nbins);
    xs[i] = sigs/hsig->GetSumOfWeights();
    dxs[i] = sqrt(sigs)/hsig->GetSumOfWeights();
    ys[i] = bgds/hbgd->GetSumOfWeights();
    dys[i] = sqrt(bgds)/hbgd->GetSumOfWeights();
    ++i;
  }
  TGraphErrors* g = new TGraphErrors(i, xs, ys, dxs, dys);
  g->SetName(name);
  return g;
}
