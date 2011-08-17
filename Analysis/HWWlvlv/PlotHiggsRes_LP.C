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
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/SmurfTree.h"
#include "factors.h"
#include "HiggsQCDScaleSystematics.h"
#include "PSUESystematics.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/LeptonScaleLookup.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void plotHistsInPad(TH1D* h1, TH1D* h2);
void setHist(TH1D* h, int color, int style);
void setPair(TH1D* h1, TH1D* h2);
void setGraph(TGraphErrors* g, int color, int marker);
TGraphErrors* makeSignificanceCurve(TH1D* sig, TH1D* bgd, TH1D* dat, const char* name);
TGraphErrors* makeGraphFromHists   (TH1D* sig, TH1D* bgd, const char* name);
double DeltaPhi(double phi1, double phi2);

int    verboseLevel =   0;
const double sigmaB = 0.35;

//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void PlotHiggsRes_LP
(
 UInt_t  nJetsType   	 = 0,
 UInt_t  mH      	 = 300,
 TString outTag          = "default",
 TString sigInputFile    = "data/inputNtuple-train-data-standard-H150_WW_2l.root",
 TString bgdInputFile    = "data/inputNtuple-data-standard-HBCK_WW_2l-train.root",
 TString datInputFile    = "data/data-train.root",
 Int_t   wwDecay         = 0,
 bool fillInfoNote       = false,
 int  TheVerboseLevel    = 0
 )
{
  verboseLevel = TheVerboseLevel;
  bool useZjetsTemplates = false;
  bool useWWTemplates    = false;
  int rebinMVAHist = 10;

  TString sigFile1 = sigInputFile;
  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  //sigFile1.ReplaceAll("data/","output/" + outTag + "_");
  //bgdFile1.ReplaceAll("data/","output/" + outTag + "_");
  //datFile1.ReplaceAll("data/","output/" + outTag + "_");  
  //sigFile1 += ".root";
  //bgdFile1 += ".root";
  //datFile1 += ".root";

  char replace[200];
  TString c1Name = sigFile1;
  c1Name.ReplaceAll("data","pic");
  sprintf(replace,"_%d.eps",wwDecay);
  c1Name.ReplaceAll(".root",replace);

  TString c2Name = sigFile1;
  c2Name.ReplaceAll("data","pic");
  sprintf(replace,"_%d_effvsbkg.eps",wwDecay);
  c2Name.ReplaceAll(".root",replace);

  TString output = sigFile1;
  output.ReplaceAll("data/","histo_tmva_");
  sprintf(replace,"_chan%d.root",wwDecay);
  output.ReplaceAll(".root",replace);

  //char output[200];
  //sprintf(output,"histo_tmva_%s_%dj_chan%d_mh%d.root",outTag.Data(),nJetsType,wwDecay,mH);     
  unsigned int patternTopTag = SmurfTree::TopTag;

  int channel = -1;
  if     (mH == 115) channel = 0;
  else if(mH == 120) channel = 1;
  else if(mH == 130) channel = 2;
  else if(mH == 140) channel = 3;
  else if(mH == 150) channel = 4;
  else if(mH == 160) channel = 5;
  else if(mH == 170) channel = 6;
  else if(mH == 180) channel = 7;
  else if(mH == 190) channel = 8;
  else if(mH == 200) channel = 9;
  else if(mH == 210) channel = 10;
  else if(mH == 220) channel = 11;
  else if(mH == 230) channel = 12;
  else if(mH == 250) channel = 13;
  else if(mH == 300) channel = 14;
  else if(mH == 350) channel = 15;
  else if(mH == 400) channel = 16;
  else if(mH == 450) channel = 17;
  else if(mH == 500) channel = 18;
  else if(mH == 550) channel = 19;
  else if(mH == 600) channel = 20;

  if(channel == -1) return;

  float dilmass_cut = 10000;
   
  if     ( mH == 115 ) dilmass_cut =  70.0;
  else if( mH == 120 ) dilmass_cut =  70.0;
  else if( mH == 130 ) dilmass_cut =  80.0;
  else if( mH == 140 ) dilmass_cut =  90.0;
  else if( mH == 150 ) dilmass_cut = 100.0;
  else if( mH == 160 ) dilmass_cut = 100.0;
  else if( mH == 165 ) dilmass_cut = 100.0;
  else if( mH == 170 ) dilmass_cut = 100.0;
  else if( mH == 180 ) dilmass_cut = 110.0;
  else if( mH == 190 ) dilmass_cut = 120.0;
  else if( mH == 200 ) dilmass_cut = 130.0;
  else if( mH == 210 ) dilmass_cut = 140.0;
  else if( mH == 220 ) dilmass_cut = 150.0;
  else                 dilmass_cut = mH;

  double wwScaleFactor0jCut[21] = {1.07,1.07,1.07,1.05,0.98,0.99,0.98,0.98,0.99,0.98,
                                   1,1,1,1,1,1,1,1,1,1,1};
  double wwScaleFactor0jMVA[21] = {1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,
                                   1,1,1,1,1,1,1,1,1,1,1};
  double wwScaleFactor1jCut[21] = {1.11,1.11,1.12,1.11,1.09,1.09,1.10,1.12,1.12,1.14,
                                   1,1,1,1,1,1,1,1,1,1,1};
  double wwScaleFactor1jMVA[21] = {1.11,1.11,1.11,1.11,1.11,1.11,1.11,1.11,1.11,1.11,
                                   1,1,1,1,1,1,1,1,1,1,1};

  double zjScaleFactor[3][21] = {
  {1.675, 1.536, 1.013, 1.126, 0.520, 0.755, 1.265, 2.132, 0.795, 0.460, 1.000, 1.000, 1.000, 0.159, 0.788, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {0.381, 0.658, 0.438, 0.374, 0.847, 0.811, 1.498, 0.698, 1.416, 1.730, 1.000, 1.000, 1.000, 1.295, 1.278, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000}
  };
  double zjScaleFactorE[3][21] = {
  {2.893, 2.909, 2.307, 1.932, 5.153, 1.952, 2.691, 2.102, 3.029, 6.022, 1.000, 1.000, 7.083, 5.063, 0.609, 0.609, 0.609, 0.609, 0.609, 0.609, 0.609},
  {2.262, 2.243, 3.837, 2.743, 1.861, 1.259, 1.033, 1.388, 1.623, 1.451, 1.000, 1.000, 1.000, 0.319, 1.270, 0.498, 0.498, 0.498, 0.498, 0.498, 0.498},
  {0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378}
  };
  double zjScaleFactorWWE[3] = {
   0.609,
   0.498,
   0.378
  };
  double ZXS_E[3] = {0.0, 0.0, 0.0};
  double DYXS[3]  = {0.0, 0.0, 0.0};
  double VVXS[3]  = {0.0, 0.0, 0.0};

  cout << "Using dilepton mass < " << dilmass_cut << endl;

  TChain *chsignal = new TChain("tree");
  chsignal->Add(sigFile1);

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
  chdata->Add(datFile1);

  TTree *signal     = (TTree*) chsignal;
  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;

  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open("/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root");
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open("/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root");
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  LeptonScaleLookup trigLookup("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");

  TFile *fNvtxFile = TFile::Open("/data/smurf/data/LP2011/auxiliar/puReweighting.root");
  TH1D *fhDNvtx = (TH1D*)(fNvtxFile->Get("puWeights"));
  assert(fhDNvtx);
  fhDNvtx->SetDirectory(0);
  fNvtxFile->Close();
  delete fNvtxFile;

  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/data/EPS/auxiliar/ggHWW_KFactors_PowhegToHQT.root");
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

  TH1D *hDZjetsTemplate;
  if(useZjetsTemplates == true) printf("***********useZjetsTemplates = true***************\n");
  else                          printf("***********useZjetsTemplates = false***************\n");
  if(useZjetsTemplates == true){
    TFile *fZjetsTemplatesFile = TFile::Open("/data/smurf/data/EPS/auxiliar/histo_Zjets_Templates.root");
    char ZjetsHistName[100];
    sprintf(ZjetsHistName, "hDZjets%d_%d", mH,TMath::Min((int)nJetsType,1));
    hDZjetsTemplate = (TH1D*)(fZjetsTemplatesFile->Get(ZjetsHistName));
    if (hDZjetsTemplate) {
      hDZjetsTemplate->SetDirectory(0);
    }
    assert(hDZjetsTemplate);
    fZjetsTemplatesFile->Close();
    delete fZjetsTemplatesFile;
    if(hDZjetsTemplate->GetSumOfWeights() != 1.0){
      printf("hDZjetsTemplate norm = %f --> normalizing to 1\n",hDZjetsTemplate->GetSumOfWeights());
      hDZjetsTemplate->Scale(1./hDZjetsTemplate->GetSumOfWeights());
    }
  }

  double scaleFactorLum = 1.545;
  //double scaleFactorLum = 1.143;
  //double scaleFactorLum = 0.402;

  //----------------------------------------------------------------------------
  double cutMassHigh[21]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[21]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[21]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[21] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[21]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[21]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};

  //----------------------------------------------------------------------------
  const int nHist = 6;
  int    nBinHis[nHist] = { 200,  200,  200,  200,  200,  200};
  double minHis[nHist]  = {-1.0, -1.0, -1.0, -0.0, -1.0,  0.0};
  double maxHis[nHist]  = { 1.0,  1.0,  1.0,  1.0,  1.0,200.0};

  //----------------------------------------------------------------------------
  TH1D* sigMVA[nHist][6];
  TH1D* bgdMVA[nHist];
  TH1D* datMVA[nHist];
  const int nChan = 8;
  TH1D* bgdMVADecays[nHist][nChan];
  for(int i=0; i<nHist; i++) {
    for(int j=0; j<6; j++){
      sigMVA[i][j] = new TH1D(Form("sigMVA_%d_%d",i,j), Form("sigMVA_%d_%d",i,j), nBinHis[i], minHis[i], maxHis[i]);
      sigMVA[i][j]->Sumw2();
    }
    bgdMVA[i] = new TH1D(Form("bgdMVA_%d",i), Form("bgdMVA_%d",i), nBinHis[i], minHis[i], maxHis[i]);
    datMVA[i] = new TH1D(Form("datMVA_%d",i), Form("datMVA_%d",i), nBinHis[i], minHis[i], maxHis[i]);
    bgdMVA[i]->Sumw2();
    datMVA[i]->Sumw2();
    for(int j=0; j<nChan; j++) {
      bgdMVADecays[i][j] = new TH1D(Form("bgdMVADecays_%d_%d",i,j), Form("bgdMVADecays_%d_%d",i,j), nBinHis[i], minHis[i], maxHis[i]);
      bgdMVADecays[i][j]->Sumw2();
    }
  }

  TH1D* histos = new TH1D("histos", "histos", nBinHis[1], minHis[1], maxHis[1]);
  histos->Sumw2();
  TH1D* histo0 = (TH1D*) histos->Clone("histo0");
  TH1D* histo1 = (TH1D*) histos->Clone("histo1");
  TH1D* histo2 = (TH1D*) histos->Clone("histo2");
  TH1D* histo3 = (TH1D*) histos->Clone("histo3");
  TH1D* histo4 = (TH1D*) histos->Clone("histo4");
  TH1D* histo5 = (TH1D*) histos->Clone("histo5");
  histos->Scale(0.0);
  histo0->Scale(0.0);
  histo1->Scale(0.0);
  histo2->Scale(0.0);
  histo3->Scale(0.0);
  histo4->Scale(0.0);
  histo5->Scale(0.0);
  TH1D* histo_Zjets_MVAZBounding     = (TH1D*) histo1->Clone("histo_Zjets_MVAZBounding");
  TH1D* histo_Zjets_MVAZBoundingUp   = (TH1D*) histo1->Clone("histo_Zjets_MVAZBoundingUp");
  TH1D* histo_Zjets_MVAZBoundingDown = (TH1D*) histo1->Clone("histo_Zjets_MVAZBoundingDown");
  TH1D* histo_qqWW_MVAWWBoundingUp   = (TH1D*) histo1->Clone("histo_qqWW_MVAWWBoundingUp");
  TH1D* histo_qqWW_MVAWWBoundingDown = (TH1D*) histo1->Clone("histo_qqWW_MVAWWBoundingDown");
  TH1D* histo_ggWW_MVAWWBoundingUp   = (TH1D*) histo1->Clone("histo_ggWW_MVAWWBoundingUp");
  TH1D* histo_ggWW_MVAWWBoundingDown = (TH1D*) histo1->Clone("histo_ggWW_MVAWWBoundingDown");

  TH1D* histoS = new TH1D("histoS", "histoS", 50, 0, 5);
  TH1D* histoB = new TH1D("histoB", "histoB", 50, 0, 5);
  TH1D* histoD = new TH1D("histoD", "histoD", 50, 0, 5);
  histoS->Sumw2();
  histoB->Sumw2();
  histoD->Sumw2();

  //----------------------------------------------------------------------------

  // TMVA cuts (known a posteriori)
  int useVar = 0; // which MVA to be used
  double useCut = -999.0;
  if    (nJetsType == 0){
    if     (mH == 115){
      //useVar = 2;
      //useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.405-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 120){
      //useVar = 2;
      //useCut =  0.675-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut =  0.375-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 130){
      //useVar = 2;
      //useCut =  0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut =  0.375-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 140){
      //useVar = 2;
      //useCut = 0.735-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.485-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 150){
      //useVar = 2;
      //useCut = 0.795-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.585-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 160){
      //useVar = 2;
      //useCut = 0.945-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.875-datMVA[useVar]->GetBinWidth(0)/2.0;
    }
    else if(mH == 170){
      //useVar = 2;
      //useCut = 0.855-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.725-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 180){
      //useVar = 2;
      //useCut = 0.855-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.715-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 190){
      //useVar = 2;
      //useCut = 0.805-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.615-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 200){
      //useVar = 2;
      //useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.625-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 250){
      //useVar = 2;
      //useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.585-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 300){
      //useVar = 2;
      //useCut = 0.805-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.605-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 350){
      //useVar = 2;
      //useCut = 0.835-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.675-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 400){
      //useVar = 2;
      //useCut = 0.845-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.705-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 450){
      //useVar = 2;
      //useCut = 0.895-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
    }
    else if(mH == 500){
      //useVar = 2;
      //useCut = 0.915-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.845-datMVA[useVar]->GetBinWidth(0)/2.0;
    }
    else if(mH == 550){
      //useVar = 2;
      //useCut = 0.925-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.845-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 600){
      //useVar = 2;
      //useCut = 0.945-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.885-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
  }
  else if(nJetsType == 1){
    if     (mH == 115){
      useVar = 4;
      useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 120){
      useVar = 4;
      useCut = 0.665-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 130){
      //useVar = 2;
      //useCut = 0.885-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 140){
      useVar = 4;
      useCut = 0.725-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 150){
      useVar = 4;
      useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 160){
      useVar = 4;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 170){
      useVar = 4;
      useCut = 0.835-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 180){
      useVar = 4;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 190){
      useVar = 4;
      useCut = 0.745-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 200){
      useVar = 4;
      useCut = 0.735-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 250){
      useVar = 4;
      useCut = 0.655-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 300){
      useVar = 4;
      useCut = 0.765-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 350){
      useVar = 4;
      useCut = 0.785-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 400){
      useVar = 4;
      useCut = 0.815-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 450){
      useVar = 4;
      useCut = 0.885-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 500){
      useVar = 4;
      useCut = 0.915-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 550){
      useVar = 4;
      useCut = 0.925-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 600){
      useVar = 4;
      useCut = 0.945-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
  }
  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          njets;
  UInt_t          run;
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
  Float_t         bdt = 0.0;
  Float_t         bdtd = 0.0;
  Float_t         nn = 0.0;
  Float_t         knn = 0.0;
  Float_t         bdtg = 0.0;
  Float_t         higgsPt = -999;

  signal->SetBranchAddress( "cuts"         , &cuts         );
  signal->SetBranchAddress( "dstype"       , &dstype       );
  signal->SetBranchAddress( "nvtx"         , &nvtx         );
  signal->SetBranchAddress( "njets"        , &njets        );
  signal->SetBranchAddress( "run"          , &run          );
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
  signal->SetBranchAddress(Form("bdt_hww%i_%djet_ww"	  ,mH,nJetsType), &bdt        );
  signal->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"	  ,mH,nJetsType), &bdtd       );
  signal->SetBranchAddress(Form("nn_hww%i_%djet_ww"	  ,mH,nJetsType), &nn	      );
  signal->SetBranchAddress(Form("knn_hww%i_%djet_ww"	  ,mH,nJetsType), &knn        );
  signal->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"	  ,mH,nJetsType), &bdtg       );
  signal->SetBranchAddress( "higgsPt"      , &higgsPt      );

  float nSigAcc[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigCut[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigMVA[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigEAcc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigECut[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigEMVA[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)signal->GetEntries());

    unsigned int Njet3 = njets;
    //if(nJetsType == 2){
    //  if(jet3->pt() <= 30)					         Njet3 = 2;
    //  else if(jet3->pt() > 30 && (
    //    (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    //    (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))   Njet3 = 0;
    //  else								 Njet3 = 2;
    //  if(njets < 2 || njets > 3) Njet3 = 0;
    //}
    // WW Preselection
    //if( mt <= 80						       ) continue; // low mt cut
    //if( mt >= mH						       ) continue; // high mt cut
    //if( mt > 80                                                          ) continue;
    if( Njet3 != nJetsType          					 ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
    if( TMath::Min(pmet,pTrackMet) <= 20                                 ) continue; // cut on pmet for all lepton-pair flavors
    if( TMath::Min(pmet,pTrackMet) <= 40 && 
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on pmet for ee/mm lepton-pair
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    //if( lid3 != 0	                                                 ) continue; // cut on dileptons
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    //if( jetLowBtag >= 2.1	    					 ) continue; // cut on anti b-tagging
    //if( nSoftMuons != 0                                                ) continue; // cut on soft muons veto
    //if( jet1Btag >= 2.1             					 ) continue; // cut on jet1Btag
    //if( jet2Btag >= 2.1             					 ) continue; // cut on jet2Btag
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1

    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);
    if(wwDecayCut == false) continue;

    double add = 1.0;
    add = add*nVtxScaleFactor(fhDNvtx,nvtx);

    add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
    add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);

    add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
        					      fabs(lep2->eta()), lep2->pt(), 
        					      TMath::Abs( lid1), TMath::Abs(lid2));

    // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
    // add = 1.0;
    // add = add * enhancementFactor(mH,true);
    if (processId == 10010) {
      float theMax = 0.00;
      add = add * HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(TMath::Max(higgsPt, theMax)));
      // add = add * enhancementFactor(mH,false);
    }
    double myWeight = scaleFactorLum * scale1fb * add;

    int nSigBin = -1;
    // GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
    if     (processId==121 ||
    	    processId==122)   nSigBin = 1;
    else if(processId==24)    nSigBin = 2;
    else if(processId==26)    nSigBin = 3;
    else if(processId==10001) nSigBin = 4;
    else if(processId==10010) nSigBin = 5;
    else  {return;}

    bool passAllCuts = dilep->mass()         < cutMassHigh[channel] &&
        	       mt		     > cutMTLow[channel] &&
        	       mt		     < cutMTHigh[channel] &&
        	       lep1->pt()	     > cutPtMaxLow[channel] &&
        	       lep2->pt()	     > cutPtMinLow[channel] &&
        	       dPhi*180.0/TMath::Pi()< cutDeltaphilHigh[channel];
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    (mH > 200 || dilep->mass() < 100.) &&
		    centrality == 1;
    }
    if(passAllCuts == true) {
      nSigCut[0]  = nSigCut[0]  + myWeight;
      nSigECut[0] = nSigECut[0] + myWeight*myWeight;
      nSigCut[nSigBin]  = nSigCut[nSigBin]  + myWeight;
      nSigECut[nSigBin] = nSigECut[nSigBin] + myWeight*myWeight;
      sigMVA[5][0]      ->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), myWeight);
      sigMVA[5][nSigBin]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), myWeight);
    }

    bool passMVAPreselCuts = mt > 80 && mt < mH;
    if(passMVAPreselCuts == true){
      nSigAcc[0]  = nSigAcc[0]  + myWeight;
      nSigEAcc[0] = nSigEAcc[0] + myWeight*myWeight;
      nSigAcc[nSigBin]  = nSigAcc[nSigBin]  + myWeight;
      nSigEAcc[nSigBin] = nSigEAcc[nSigBin] + myWeight*myWeight;

      if     (useVar == 0)
        histos   ->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),	myWeight);
      else if(useVar == 1)
        histos   ->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),	myWeight);
      else if(useVar == 2)
        histos   ->Fill(TMath::Max(TMath::Min((double)nn  ,maxHis[2]-0.001),minHis[2]+0.001),	myWeight);
      else if(useVar == 4)
        histos   ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),	myWeight);

      sigMVA[0][0]->Fill(TMath::Max(TMath::Min((double)bdt,maxHis[0]-0.001),minHis[0]+0.001),	 myWeight);
      sigMVA[1][0]->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);
      sigMVA[2][0]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),	 myWeight);
      sigMVA[3][0]->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),	 myWeight);
      sigMVA[4][0]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeight);
      sigMVA[0][nSigBin]->Fill(TMath::Max(TMath::Min((double)bdt,maxHis[0]-0.001),minHis[0]+0.001),    myWeight);
      sigMVA[1][nSigBin]->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);
      sigMVA[2][nSigBin]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),	    myWeight);
      sigMVA[3][nSigBin]->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),    myWeight);
      sigMVA[4][nSigBin]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeight);

      bool passFinalCut = false;
      if     (useVar == 0 && bdt	> useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 1 && bdtd	> useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 2 && nn 	> useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 3 && knn	> useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 4 && bdtg	> useCut) {nSigMVA[0] += myWeight; nSigEMVA[0] += myWeight*myWeight; passFinalCut = true;}
      if     (useVar == 0 && bdt	> useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 1 && bdtd	> useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 2 && nn 	> useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 3 && knn	> useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 4 && bdtg	> useCut) {nSigMVA[nSigBin] += myWeight; nSigEMVA[nSigBin] += myWeight*myWeight; passFinalCut = true;}
      if(passFinalCut == true){
        double myVar = dR;
        //if(jet1->pt()>15&&jet2->pt()>15) myVar = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi();
        histoS->Fill(myVar,myWeight);
      }
    } // passMVAPreselCuts
  }
  for(int i=0; i<6; i++) nSigEAcc[i] = sqrt(nSigEAcc[i]);
  for(int i=0; i<6; i++) nSigECut[i] = sqrt(nSigECut[i]);
  for(int i=0; i<6; i++) nSigEMVA[i] = sqrt(nSigEMVA[i]);
  printf("--- Finished Signal loop\n");

  background->SetBranchAddress( "cuts"          , &cuts 	  );
  background->SetBranchAddress( "dstype"        , &dstype	  );
  background->SetBranchAddress( "nvtx"          , &nvtx 	  );
  background->SetBranchAddress( "njets"         , &njets	  );
  background->SetBranchAddress( "run"           , &run	          );
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
  background->SetBranchAddress(Form("bdt_hww%i_%djet_ww"      ,mH,nJetsType), &bdt	  );
  background->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"     ,mH,nJetsType), &bdtd	  );
  background->SetBranchAddress(Form("nn_hww%i_%djet_ww"       ,mH,nJetsType), &nn	  );
  background->SetBranchAddress(Form("knn_hww%i_%djet_ww"      ,mH,nJetsType), &knn	  );
  background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"     ,mH,nJetsType), &bdtg	  );

  float nBgdAcc = 0.0;
  float nBgdCut = 0.0;
  float nBgdMVA = 0.0;
  float nBgdEAcc = 0.0;
  float nBgdECut = 0.0;
  float nBgdEMVA = 0.0;
  float nBgdAccDecays[nChan]  = {0.,0.,0.,0.,0.,0.,0.,0.};
  float nBgdCutDecays[nChan]  = {0.,0.,0.,0.,0.,0.,0.,0.};
  float nBgdMVADecays[nChan]  = {0.,0.,0.,0.,0.,0.,0.,0.};
  float nBgdEAccDecays[nChan] = {0.,0.,0.,0.,0.,0.,0.,0.};
  float nBgdECutDecays[nChan] = {0.,0.,0.,0.,0.,0.,0.,0.};
  float nBgdEMVADecays[nChan] = {0.,0.,0.,0.,0.,0.,0.,0.};
  for (UInt_t i=0; i<background->GetEntries(); i++) {

    background->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    if(dstype == SmurfTree::data &&
      (cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    //if(dstype == SmurfTree::data && run >= 170826) continue;
    //if(dstype == SmurfTree::data && run <  170826) continue;

    unsigned int Njet3 = njets;
    //if(nJetsType == 2){
    //  if(jet3->pt() <= 30)                                               Njet3 = 2;
    //  else if(jet3->pt() > 30 && (
    //    (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    //    (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))   Njet3 = 0;
    // else                                                               Njet3 = 2;
    // if(njets < 2 || njets > 3) Njet3 = 0;
    //}
    // WW Preselection
    bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
                        ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
                         dstype != SmurfTree::data;

    if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
    //if( mt <= 80						       ) continue; // low mt cut
    //if( mt >= mH						       ) continue; // high mt cut
    //if( mt > 80                                                          ) continue;
    if( Njet3 != nJetsType          					 ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
    if( TMath::Min(pmet,pTrackMet) <= 20                                 ) continue; // cut on pmet for all lepton-pair flavors
    if( TMath::Min(pmet,pTrackMet) <= 40 && 
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on pmet for ee/mm lepton-pair
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    //if( lid3 != 0	                                                 ) continue; // cut on dileptons
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    //if( jetLowBtag >= 2.1	    					 ) continue; // cut on anti b-tagging
    //if( nSoftMuons != 0		    			         ) continue; // cut on soft muons veto
    //if( jet1Btag >= 2.1             					 ) continue; // cut on jet1Btag
    //if( jet2Btag >= 2.1             					 ) continue; // cut on jet2Btag
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1

    int fDecay = 0;
    if     (dstype == SmurfTree::wjets ) fDecay = 5;
    else if(dstype == SmurfTree::ttbar ) fDecay = 3;
    else if(dstype == SmurfTree::dyee  ) fDecay = 4;
    else if(dstype == SmurfTree::dymm  ) fDecay = 4;
    else if(dstype == SmurfTree::dytt  ) fDecay = 7;
    else if(dstype == SmurfTree::tw    ) fDecay = 3;
    else if(dstype == SmurfTree::qqww  ) fDecay = 0;
    else if(dstype == SmurfTree::wz    ) fDecay = 2;
    else if(dstype == SmurfTree::zz    ) fDecay = 2;
    else if(dstype == SmurfTree::ggww  ) fDecay = 1;
    else if(dstype == SmurfTree::wgamma) fDecay = 6;
    else if(dstype == SmurfTree::data  ) fDecay = 5;
    else                                 {printf("bad dstype: %d\n",dstype); assert(0);}
    if(dstype == SmurfTree::wz || dstype == SmurfTree::zz) {
      if(lep1MotherMcId == 23 && lep2MotherMcId == 23) {
        fDecay = 4;
      }
    }

    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);
    if(wwDecayCut == false) continue;

    double myWeight = 1.0;
    double add      = 1.0;
    int nFake = 0;
    if(dstype == SmurfTree::data ){
      if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    } else {
      if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    }
    if(nFake > 1){
      myWeight = 0.0;
    }
    else if(nFake == 1){
      if(dstype == SmurfTree::data){
        add = add*fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        							      (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        							      (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        fDecay  	      = 5;
        myWeight	      = add;
      }
      else if(TMath::Abs(lep1McId)*TMath::Abs(lep2McId) > 0 || dstype == SmurfTree::wgamma){
    	add = 1.0;
    	add = add*fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    										      (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    										      (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

    	add = add*nVtxScaleFactor(fhDNvtx,nvtx);

    	add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
    	add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);

    	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
        fDecay  	       = 5;
        myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
      }
      else {
        myWeight = 0.0;
      }
    }
    else if(dstype == SmurfTree::data) myWeight = 0.0;
    else if(dstype != SmurfTree::data){
      add = 1.0;
      add = add*nVtxScaleFactor(fhDNvtx,nvtx);

      add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
      add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
      add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							fabs(lep2->eta()), lep2->pt(), 
    							TMath::Abs( lid1), TMath::Abs(lid2));
      if(fDecay == 4  && (type   == SmurfTree::mm   || type   == SmurfTree::ee)
                      && (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm)) {
    	if(njets == 0) add=add*3.02; 
    	if(njets == 1) add=add*2.81; 
    	if(njets >= 2) add=add*4.84; 
      }
      if(fDecay == 3) {
    	if(njets == 0) add=add*1.51;
    	if(njets == 1) add=add*1.16; 
    	if(njets >= 2) add=add*1.00; 
      }
      if(fDecay == 5) {
    	add=add*1.95; 
      }
      if((fDecay == 0 || fDecay == 1)){     
        if(njets == 0) add=add*wwScaleFactor0jMVA[channel]; 
        else	       add=add*wwScaleFactor1jMVA[channel]; 
      }
      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;
      myWeight = scale1fb*scaleFactorLum*add;
    }

    if(myWeight == 0) continue;

    bool passAllCuts = dilep->mass()         < cutMassHigh[channel] &&
        	       mt		     > cutMTLow[channel] &&
        	       mt		     < cutMTHigh[channel] &&
        	       lep1->pt()	     > cutPtMaxLow[channel] &&
        	       lep2->pt()	     > cutPtMinLow[channel] &&
        	       dPhi*180.0/TMath::Pi()< cutDeltaphilHigh[channel];
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    (mH > 200 || dilep->mass() < 100.) &&
		    centrality == 1;
    }
    if(passAllCuts == true) {
      double newWeight = myWeight;
      if((fDecay == 0 || fDecay == 1)){ // only for WW
	if(njets == 0) newWeight=newWeight*wwScaleFactor0jCut[channel]/wwScaleFactor0jMVA[channel];
	else           newWeight=newWeight*wwScaleFactor1jCut[channel]/wwScaleFactor1jMVA[channel];      
      }
      if((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
         (type   == SmurfTree::mm   || type   == SmurfTree::ee)){     
        newWeight=newWeight*zjScaleFactor[TMath::Min((int)nJetsType,2)][channel];
	DYXS[1] += newWeight;
      }
      else if(fDecay == 4){
	VVXS[1] += newWeight;
      }
      nBgdCut  = nBgdCut  + newWeight;
      nBgdECut = nBgdECut + newWeight*newWeight;
      nBgdCutDecays[fDecay]  = nBgdCutDecays[fDecay]  + newWeight;
      nBgdECutDecays[fDecay] = nBgdECutDecays[fDecay] + newWeight*newWeight;
      bgdMVA[5]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), newWeight);
      bgdMVADecays[5][fDecay]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), newWeight);
    }

    bool passMVAPreselCuts = mt > 80 && mt < mH;
    if(passMVAPreselCuts == true){
      nBgdAcc  = nBgdAcc  + myWeight;
      nBgdEAcc = nBgdEAcc + myWeight*myWeight;
      nBgdAccDecays[fDecay]  = nBgdAccDecays[fDecay]  + myWeight;
      nBgdEAccDecays[fDecay] = nBgdEAccDecays[fDecay] + myWeight*myWeight;
      double myWeightMVA = myWeight;
      if((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
	 (type   == SmurfTree::mm   || type   == SmurfTree::ee)){
	DYXS[0] += myWeight;
	if(useZjetsTemplates == true) {histo_Zjets_MVAZBoundingUp->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);myWeightMVA = 0.0;}
      }
      else if(fDecay == 4){
	VVXS[0] += myWeight;
      }

      if      (fDecay == 4 || fDecay == 7){ // Z+jets
	if     (useVar == 0) histo1->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),   myWeightMVA);
	else if(useVar == 1) histo1->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeightMVA);
	else if(useVar == 2) histo1->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeightMVA);
	else if(useVar == 4) histo1->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeightMVA);
      }
      else if(fDecay == 3){ // ttbar
	if     (useVar == 0) histo2->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),   myWeightMVA);
	else if(useVar == 1) histo2->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeightMVA);
	else if(useVar == 2) histo2->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeightMVA);
	else if(useVar == 4) histo2->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeightMVA);
      }
      else if(fDecay == 2){ // WZ/ZZ
	if     (useVar == 0) histo3->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),   myWeightMVA);
	else if(useVar == 1) histo3->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeightMVA);
	else if(useVar == 2) histo3->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeightMVA);
	else if(useVar == 4) histo3->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeightMVA);
      }
      else if(fDecay == 5 || fDecay == 6){ // W+X
	if     (useVar == 0) histo4->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),   myWeightMVA);
	else if(useVar == 1) histo4->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeightMVA);
	else if(useVar == 2) histo4->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeightMVA);
	else if(useVar == 4) histo4->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeightMVA);
      }
      else if(fDecay == 0 || fDecay == 1){ // WW
	if     (useVar == 0) histo0->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),  myWeightMVA);
	else if(useVar == 1) histo0->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),  myWeightMVA);
	else if(useVar == 2) histo0->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),    myWeightMVA);
	else if(useVar == 4) histo0->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
      }

      bgdMVA[0]->Fill(TMath::Max(TMath::Min((double)bdt,maxHis[0]-0.001),minHis[0]+0.001),    myWeightMVA);
      bgdMVA[1]->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeightMVA);
      bgdMVA[2]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeightMVA);
      bgdMVA[3]->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),    myWeightMVA);
      bgdMVA[4]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeightMVA);
      bgdMVADecays[0][fDecay]->Fill(TMath::Max(TMath::Min((double)bdt,maxHis[0]-0.001),minHis[0]+0.001),    myWeightMVA);
      bgdMVADecays[1][fDecay]->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeightMVA);
      bgdMVADecays[2][fDecay]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeightMVA);
      bgdMVADecays[3][fDecay]->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),    myWeightMVA);
      bgdMVADecays[4][fDecay]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeightMVA);

      bool passFinalCut = false;
      if     (useVar == 0 && bdt        > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 1 && bdtd       > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 2 && nn         > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 3 && knn        > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 4 && bdtg       > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      if(passFinalCut == true){
	if((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
           (type   == SmurfTree::mm   || type   == SmurfTree::ee)){     
	  DYXS[2] += myWeight;
	}
	else if(fDecay == 4){
	  VVXS[2] += myWeight;
	}
	double myVar = dR;
	//if(jet1->pt()>15&&jet2->pt()>15) myVar = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi();
	histoB->Fill(myVar,myWeight);
      }
    } // passMVAPreselCuts
  }
  nBgdEAcc = sqrt(nBgdEAcc);
  nBgdECut = sqrt(nBgdECut);
  nBgdEMVA = sqrt(nBgdEMVA);
  printf("--- Finished Bgdnal loop\n");

  data->SetBranchAddress( "cuts"         , &cuts         );
  data->SetBranchAddress( "dstype"       , &dstype	 );
  data->SetBranchAddress( "nvtx"         , &nvtx         );
  data->SetBranchAddress( "njets"        , &njets        );
  data->SetBranchAddress( "run"          , &run          );
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
  data->SetBranchAddress(Form("bdt_hww%i_%djet_ww"	,mH,nJetsType), &bdt	    );
  data->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"	,mH,nJetsType), &bdtd	    );
  data->SetBranchAddress(Form("nn_hww%i_%djet_ww"	,mH,nJetsType), &nn	    );
  data->SetBranchAddress(Form("knn_hww%i_%djet_ww"	,mH,nJetsType), &knn	    );
  data->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"	,mH,nJetsType), &bdtg	    );

  float nDatAcc = 0.0;
  float nDatCut = 0.0;
  float nDatMVA = 0.0;
  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    if((cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    //if(dstype == SmurfTree::data && run >= 170826) continue;
    //if(dstype == SmurfTree::data && run <  170826) continue;
    //if(event==582505290||event==325565599) continue;

    unsigned int Njet3 = njets;
    //if(nJetsType == 2){
    //  if(jet3->pt() <= 30)					         Njet3 = 2;
    //  else if(jet3->pt() > 30 && (
    //    (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    //    (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))   Njet3 = 0;
    //  else								 Njet3 = 2;
    //  if(njets < 2 || njets > 3) Njet3 = 0;
    //}
    // WW Preselection
    //if( mt <= 80						       ) continue; // low mt cut
    //if( mt >= mH						       ) continue; // high mt cut
    //if( mt > 80                                                          ) continue;
    if( Njet3 != nJetsType          					 ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
    if( TMath::Min(pmet,pTrackMet) <= 20                                 ) continue; // cut on pmet for all lepton-pair flavors
    if( TMath::Min(pmet,pTrackMet) <= 40 && 
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on pmet for ee/mm lepton-pair
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    //if( lid3 != 0	                                                 ) continue; // cut on dileptons
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    //if( jetLowBtag >= 2.1	    					 ) continue; // cut on anti b-tagging
    //if( nSoftMuons != 0                                                ) continue; // cut on soft muons veto
    //if( jet1Btag >= 2.1             					 ) continue; // cut on jet1Btag
    //if( jet2Btag >= 2.1             					 ) continue; // cut on jet2Btag
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1

    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);
    if(wwDecayCut == false) continue;

    double myWeight = 1.0;

    bool passAllCuts = dilep->mass()         < cutMassHigh[channel] &&
        	       mt		     > cutMTLow[channel] &&
        	       mt		     < cutMTHigh[channel] &&
        	       lep1->pt()	     > cutPtMaxLow[channel] &&
        	       lep2->pt()	     > cutPtMinLow[channel] &&
        	       dPhi*180.0/TMath::Pi()< cutDeltaphilHigh[channel];
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    (mH > 200 || dilep->mass()< 100.) &&
		    centrality == 1;
    }
    if(passAllCuts == true) {
      nDatCut = nDatCut + myWeight;
      datMVA[5]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), myWeight);
    }

    bool passMVAPreselCuts = mt > 80 && mt < mH;
    if(passMVAPreselCuts == true){
      nDatAcc = nDatAcc + myWeight;
      if     (useVar == 0) histo5->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),   myWeight);
      else if(useVar == 1) histo5->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);
      else if(useVar == 2) histo5->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeight);
      else if(useVar == 4) histo5->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeight);

      datMVA[0]->Fill(TMath::Max(TMath::Min((double)bdt,maxHis[0]-0.001),minHis[0]+0.001),    myWeight);
      datMVA[1]->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);
      datMVA[2]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeight);
      datMVA[3]->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),    myWeight);
      datMVA[4]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      bool passFinalCut = false;
      if     (useVar == 0 && bdt        > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 1 && bdtd       > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 2 && nn         > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 3 && knn        > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 4 && bdtg       > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      if(passFinalCut == true){
	double myVar = dR;
	//if(jet1->pt()>15&&jet2->pt()>15) myVar = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi();
	histoD->Fill(myVar,myWeight);
      }
    } // passMVAPreselCuts
  }

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
  for(int i=0; i<nChan; i++) {
    if(nBgdAccDecays[i] > 0.0) nBgdEAccDecays[i] = sqrt(nBgdEAccDecays[i])/nBgdAccDecays[i];
    if(nBgdCutDecays[i] > 0.0) nBgdECutDecays[i] = sqrt(nBgdECutDecays[i])/nBgdCutDecays[i];
    if(nBgdMVADecays[i] > 0.0) nBgdEMVADecays[i] = sqrt(nBgdEMVADecays[i])/nBgdMVADecays[i];
  }
  printf("            HWW      qqww     ggww     VV       top      dyll     wjets    vg       Ztautau\n");
  printf("CLsAcc : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigAcc[0],nBgdAccDecays[0],nBgdAccDecays[1],nBgdAccDecays[2],nBgdAccDecays[3],nBgdAccDecays[4],nBgdAccDecays[5],nBgdAccDecays[6],nBgdAccDecays[7]);
  printf("CLsCut : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigCut[0],nBgdCutDecays[0],nBgdCutDecays[1],nBgdCutDecays[2],nBgdCutDecays[3],nBgdCutDecays[4],nBgdCutDecays[5],nBgdCutDecays[6],nBgdCutDecays[7]);
  printf("CLsMVA : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigMVA[0],nBgdMVADecays[0],nBgdMVADecays[1],nBgdMVADecays[2],nBgdMVADecays[3],nBgdMVADecays[4],nBgdMVADecays[5],nBgdMVADecays[6],nBgdMVADecays[7]);
  printf("CLsEAcc: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigEAcc[0]/nSigAcc[0],nBgdEAccDecays[0],nBgdEAccDecays[1],nBgdEAccDecays[2],nBgdEAccDecays[3],nBgdEAccDecays[4],nBgdEAccDecays[5],nBgdEAccDecays[6],nBgdEAccDecays[7]);
  printf("CLsECut: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigECut[0]/nSigCut[0],nBgdECutDecays[0],nBgdECutDecays[1],nBgdECutDecays[2],nBgdECutDecays[3],nBgdECutDecays[4],nBgdECutDecays[5],nBgdECutDecays[6],nBgdECutDecays[7]);
  printf("CLsEMVA: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSigEMVA[0]/nSigMVA[0],nBgdEMVADecays[0],nBgdEMVADecays[1],nBgdEMVADecays[2],nBgdEMVADecays[3],nBgdEMVADecays[4],nBgdEMVADecays[5],nBgdEMVADecays[6],nBgdEMVADecays[7]);

  // uncertainties for DY
  // VERY CAREFUL HERE!!!, if we don't use bdtd or change binning!,it has to be changed
  if(useZjetsTemplates == true){
    hDZjetsTemplate->Scale(nBgdAccDecays[4]);
    histo_Zjets_MVAZBounding->Add(hDZjetsTemplate);
    histo_Zjets_MVAZBounding->Rebin(rebinMVAHist);
    histo_Zjets_MVAZBoundingUp->Rebin(rebinMVAHist);
    histo_Zjets_MVAZBoundingDown->Rebin(rebinMVAHist);
    for(int i=1; i<=histo_Zjets_MVAZBounding->GetNbinsX(); i++){
      double mean = histo_Zjets_MVAZBounding  ->GetBinContent(i);
      double up   = histo_Zjets_MVAZBoundingUp->GetBinContent(i);
      double d    = TMath::Abs(mean-up);
      if     (mean-up >0) histo_Zjets_MVAZBoundingDown->SetBinContent(i,TMath::Max(mean+d,0.0));
      else                histo_Zjets_MVAZBoundingDown->SetBinContent(i,TMath::Max(mean-d,0.0));
    }
    if(histo_Zjets_MVAZBoundingDown->GetSumOfWeights() > 0)
      histo_Zjets_MVAZBoundingDown->Scale(histo_Zjets_MVAZBoundingUp->GetSumOfWeights()/histo_Zjets_MVAZBoundingDown->GetSumOfWeights());
    TH1D* histoVV = (TH1D*) histo1->Clone("histoVV");
    histoVV->Rebin(rebinMVAHist);
    histo_Zjets_MVAZBoundingUp->Add(histoVV);
    histo_Zjets_MVAZBoundingDown->Add(histoVV);

    bgdMVA[4]->Add(hDZjetsTemplate);
    bgdMVADecays[4][4]->Add(hDZjetsTemplate);
    histo1->Add(hDZjetsTemplate);
  }
  if(useWWTemplates == true){
    printf("WW Bins: %d\n",histo_qqWW_MVAWWBoundingUp->GetNbinsX());
    for(int i=1; i<=histo_qqWW_MVAWWBoundingUp->GetNbinsX(); i++){
      histo_qqWW_MVAWWBoundingUp  ->SetBinContent(i,bgdMVADecays[useVar][0]->GetBinContent(i)*(1.0+(i-(histo_qqWW_MVAWWBoundingUp  ->GetNbinsX()/2)+1)*0.2/30.));
      histo_ggWW_MVAWWBoundingUp  ->SetBinContent(i,bgdMVADecays[useVar][1]->GetBinContent(i)*(1.0+(i-(histo_ggWW_MVAWWBoundingUp  ->GetNbinsX()/2)+1)*0.2/30.));
    }
    for(int i=1; i<=histo_ggWW_MVAWWBoundingDown->GetNbinsX(); i++){
      histo_qqWW_MVAWWBoundingDown->SetBinContent(i,bgdMVADecays[useVar][0]->GetBinContent(i)*(1.0+((histo_qqWW_MVAWWBoundingDown->GetNbinsX()/2)+1-i)*0.2/30.));
      histo_ggWW_MVAWWBoundingDown->SetBinContent(i,bgdMVADecays[useVar][1]->GetBinContent(i)*(1.0+((histo_ggWW_MVAWWBoundingDown->GetNbinsX()/2)+1-i)*0.2/30.));
    }
    histo_qqWW_MVAWWBoundingUp  ->Scale(bgdMVADecays[useVar][0]->GetSumOfWeights()/histo_qqWW_MVAWWBoundingUp  ->GetSumOfWeights());
    histo_qqWW_MVAWWBoundingDown->Scale(bgdMVADecays[useVar][0]->GetSumOfWeights()/histo_qqWW_MVAWWBoundingDown->GetSumOfWeights());
    histo_ggWW_MVAWWBoundingUp  ->Scale(bgdMVADecays[useVar][1]->GetSumOfWeights()/histo_ggWW_MVAWWBoundingUp  ->GetSumOfWeights());
    histo_ggWW_MVAWWBoundingDown->Scale(bgdMVADecays[useVar][1]->GetSumOfWeights()/histo_ggWW_MVAWWBoundingDown->GetSumOfWeights());
  }

  if(TMath::Abs(DYXS[0]+VVXS[0]-nBgdAccDecays[4]) > 0.01) {printf("Problem: %f %f %f\n",DYXS[0],VVXS[0],nBgdAccDecays[4]); assert(0);}
  if(TMath::Abs(DYXS[1]+VVXS[1]-nBgdCutDecays[4]) > 0.01) {printf("Problem: %f %f %f\n",DYXS[1],VVXS[1],nBgdCutDecays[4]); assert(0);}
  if(TMath::Abs(DYXS[2]+VVXS[2]-nBgdMVADecays[4]) > 0.01) {printf("Problem: %f %f %f\n",DYXS[2],VVXS[2],nBgdMVADecays[4]); assert(0);}
  if(nBgdAccDecays[4] > 0.0) ZXS_E[0] = sqrt(DYXS[0]*DYXS[0]*zjScaleFactorWWE[TMath::Min((int)nJetsType,2)]*zjScaleFactorWWE[TMath::Min((int)nJetsType,2)]+
                                             VVXS[0]*VVXS[0]*0.10*0.10)/nBgdAccDecays[4]; else ZXS_E[0] = 0;
  if(nBgdCutDecays[4] > 0.0) ZXS_E[1] = sqrt(DYXS[1]*DYXS[1]*zjScaleFactorE[TMath::Min((int)nJetsType,2)][channel]*zjScaleFactorE[TMath::Min((int)nJetsType,2)][channel]+
                                             VVXS[1]*VVXS[1]*0.10*0.10)/nBgdCutDecays[4]; else ZXS_E[1] = 0;
  if(nBgdMVADecays[4] > 0.0) ZXS_E[2] = sqrt(DYXS[2]*DYXS[2]*zjScaleFactorE[TMath::Min((int)nJetsType,2)][channel]*zjScaleFactorE[TMath::Min((int)nJetsType,2)][channel]+
                                             VVXS[2]*VVXS[2]*0.10*0.10)/nBgdMVADecays[4]; else ZXS_E[2] = 0;
  //----------------------------------------------------------------------------
  // Drawing part - c1
  //----------------------------------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","c1",100,100,700,800);
  c1->Divide(2,3);

  c1->cd(1); plotHistsInPad(sigMVA[0][0], bgdMVA[0]);
  c1->cd(2); plotHistsInPad(sigMVA[1][0], bgdMVA[1]);
  c1->cd(3); plotHistsInPad(sigMVA[2][0], bgdMVA[2]);
  c1->cd(4); plotHistsInPad(sigMVA[3][0], bgdMVA[3]);
  c1->cd(5); plotHistsInPad(sigMVA[4][0], bgdMVA[4]);
  c1->cd(6); plotHistsInPad(sigMVA[5][0], bgdMVA[5]);

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

  TGraphErrors* g3BDT         = makeSignificanceCurve(sigMVA[0][0], bgdMVA[0], datMVA[0],"g3BDT"    );
  TGraphErrors* g3BDTD        = makeSignificanceCurve(sigMVA[1][0], bgdMVA[1], datMVA[1],"g3BDTD"   );    
  TGraphErrors* g3MLP         = makeSignificanceCurve(sigMVA[2][0], bgdMVA[2], datMVA[2],"g3NN"	    );  
  TGraphErrors* g3KNN         = makeSignificanceCurve(sigMVA[3][0], bgdMVA[3], datMVA[3],"g3KNN"    );
  TGraphErrors* g3BDTG        = makeSignificanceCurve(sigMVA[4][0], bgdMVA[4], datMVA[4],"g3BDTG"   );  
  TGraphErrors* g3MT          = makeSignificanceCurve(sigMVA[5][0], bgdMVA[5], datMVA[5],"g3MT"	    );  
  setGraph(g3BDT       ,1, 20);
  setGraph(g3BDTD      ,2, 21);
  setGraph(g3MLP       ,3, 23);
  setGraph(g3KNN       ,4, 24);
  setGraph(g3BDTG      ,5, 25);
  setGraph(g3MT        ,6, 26);

  TCanvas* c4 = new TCanvas("c4","c4",100,100,700,800);
  c4->cd();
  gPad->SetGrid(1,1);
  TGraphErrors* g4BDT         = makeGraphFromHists(sigMVA[0][0], bgdMVA[0], "g3BDT"       );
  TGraphErrors* g4BDTD        = makeGraphFromHists(sigMVA[1][0], bgdMVA[1], "g3BDTD"      );	   
  TGraphErrors* g4MLP         = makeGraphFromHists(sigMVA[2][0], bgdMVA[2], "g3NN"        );  
  TGraphErrors* g4KNN         = makeGraphFromHists(sigMVA[3][0], bgdMVA[3], "g3KNN"       );
  TGraphErrors* g4BDTG        = makeGraphFromHists(sigMVA[4][0], bgdMVA[4], "g3BDTG"      );  
  TGraphErrors* g4MT          = makeGraphFromHists(sigMVA[5][0], bgdMVA[5], "g3MT"        );  
  setGraph(g4BDT       ,1, 20);
  setGraph(g4BDTD      ,2, 21);
  setGraph(g4MLP       ,3, 23);
  setGraph(g4KNN       ,4, 24);
  setGraph(g4BDTG      ,5, 25);
  setGraph(g4MT        ,6, 26);
  g4BDT       ->Draw("APXl");
  g4BDTD      ->Draw("PXl");
  g4MLP       ->Draw("PXl");
  g4KNN       ->Draw("PXl");
  g4BDTG      ->Draw("PXl");
  g4MT        ->Draw("PXl");
  TAxis *axisX=g4BDT->GetXaxis();
  TAxis *axisY=g4BDT->GetYaxis();
  axisX->SetTitle("Signal Efficiency");
  axisX->Draw();
  axisY->SetTitle("Background Efficiency");
  axisY->Draw();
  TLegend* leg4 = new TLegend(0.2, 0.7, 0.5, 0.9);
  leg4->SetFillColor(0);
  leg4->AddEntry(g4BDT       ,"BDT"       ,"L");
  leg4->AddEntry(g4BDTD      ,"BDTD"      ,"L");
  leg4->AddEntry(g4MLP       ,"NN"        ,"L");
  leg4->AddEntry(g4KNN       ,"KNN"       ,"L");
  leg4->AddEntry(g4BDTG      ,"BDTG"      ,"L");
  leg4->AddEntry(g4MT	     ,"MT"        ,"L");
  leg4->Draw("same");
  c4->SaveAs(c2Name);

  char finalStateName[10];
  sprintf(finalStateName,"ll");
  if	 (wwDecay == 0) sprintf(finalStateName,"mm");
  else if(wwDecay == 1) sprintf(finalStateName,"me");
  else if(wwDecay == 2) sprintf(finalStateName,"em");
  else if(wwDecay == 3) sprintf(finalStateName,"ee");
  else if(wwDecay == 5) sprintf(finalStateName,"sf");
  else if(wwDecay == 6) sprintf(finalStateName,"of");

  //char output[200];
  //sprintf(output,"histo_tmva_%s_%dj_chan%d_mh%d.root",outTag.Data(),nJetsType,wwDecay,mH);     
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
  for(int i=0; i<nHist; i++) {
    histoS->Write();
    histoB->Write();
    histoD->Write();    
    for(int j=0; j<4; j++){
      sigMVA[i][j]->Write();
    }
    bgdMVA[i]->Write();
    datMVA[i]->Write();
    for(int j=0; j<5; j++) {
      bgdMVADecays[i][j]->Write();
    }
    histos->Write();
    histo0->Write();
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();
  }
  outFilePlotsNote->Close();

  if(fillInfoNote == true){
    cout << sigMVA[useVar][1]->GetSumOfWeights() << " ";
    cout << sigMVA[useVar][2]->GetSumOfWeights() << " ";
    cout << sigMVA[useVar][3]->GetSumOfWeights() << " ";
    cout << sigMVA[useVar][4]->GetSumOfWeights() << " ";
    cout << sigMVA[useVar][5]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][0]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][1]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][2]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][3]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][4]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][5]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][6]->GetSumOfWeights() << " ";
    cout << bgdMVADecays[useVar][7]->GetSumOfWeights() << endl;
    TH1D *histo_ttH    = (TH1D*) sigMVA[useVar][1]->Clone("histo_ttH");
    TH1D *histo_ZH     = (TH1D*) sigMVA[useVar][2]->Clone("histo_ZH");
    TH1D *histo_WH     = (TH1D*) sigMVA[useVar][3]->Clone("histo_WH");
    TH1D *histo_qqH    = (TH1D*) sigMVA[useVar][4]->Clone("histo_qqH");
    TH1D *histo_ggH    = (TH1D*) sigMVA[useVar][5]->Clone("histo_ggH");
    TH1D *histo_Data   = (TH1D*) datMVA[useVar]->Clone("histo_Data");
    TH1D *histo_qqWW   = (TH1D*) bgdMVADecays[useVar][0]->Clone("histo_qqWW");
    TH1D *histo_ggWW   = (TH1D*) bgdMVADecays[useVar][1]->Clone("histo_ggWW");
    TH1D *histo_VV     = (TH1D*) bgdMVADecays[useVar][2]->Clone("histo_VV");
    TH1D *histo_Top    = (TH1D*) bgdMVADecays[useVar][3]->Clone("histo_Top");
    TH1D *histo_Zjets  = (TH1D*) bgdMVADecays[useVar][4]->Clone("histo_Zjets");
    TH1D *histo_Wjets  = (TH1D*) bgdMVADecays[useVar][5]->Clone("histo_Wjets");
    TH1D *histo_Wgamma = (TH1D*) bgdMVADecays[useVar][6]->Clone("histo_Wgamma");
    TH1D *histo_Ztt    = (TH1D*) bgdMVADecays[useVar][7]->Clone("histo_Ztt");
    histo_ggH    ->Rebin(rebinMVAHist);
    histo_qqH    ->Rebin(rebinMVAHist);
    histo_WH     ->Rebin(rebinMVAHist);
    histo_ZH     ->Rebin(rebinMVAHist);
    histo_ttH    ->Rebin(rebinMVAHist);
    histo_Data   ->Rebin(rebinMVAHist);
    histo_qqWW   ->Rebin(rebinMVAHist);
    histo_ggWW   ->Rebin(rebinMVAHist);
    histo_VV	 ->Rebin(rebinMVAHist);
    histo_Top	 ->Rebin(rebinMVAHist);
    histo_Zjets  ->Rebin(rebinMVAHist);
    histo_Wjets  ->Rebin(rebinMVAHist);
    histo_Wgamma ->Rebin(rebinMVAHist);
    histo_Ztt    ->Rebin(rebinMVAHist);
    char outputLimits[200];
    //sprintf(outputLimits,"output/histo_limits_%s_%dj_chan%d_mh%d.root",outTag.Data(),nJetsType,wwDecay,mH);     
    sprintf(outputLimits,"hww%s_%dj.input.root",finalStateName,nJetsType);     
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
      histo_ggH    ->Write();
      histo_qqH    ->Write();
      histo_WH     ->Write();
      histo_ZH     ->Write();
      histo_ttH    ->Write();
      histo_Data   ->Write();
      histo_qqWW   ->Write();
      histo_ggWW   ->Write();
      histo_VV	   ->Write();
      histo_Top	   ->Write();
      histo_Zjets  ->Write();
      histo_Wjets  ->Write();
      histo_Wgamma ->Write();
      histo_Ztt    ->Write();
      if(wwDecay == 5 && useZjetsTemplates == true){
        histo_Zjets_MVAZBoundingUp->Write();
        histo_Zjets_MVAZBoundingDown->Write();
      }
      if(useWWTemplates == true){
        histo_qqWW_MVAWWBoundingUp  ->Rebin(rebinMVAHist);
        histo_qqWW_MVAWWBoundingDown->Rebin(rebinMVAHist);
        histo_qqWW_MVAWWBoundingUp  ->Write();
        histo_qqWW_MVAWWBoundingDown->Write();
        histo_ggWW_MVAWWBoundingUp  ->Rebin(rebinMVAHist);
        histo_ggWW_MVAWWBoundingDown->Rebin(rebinMVAHist);
        histo_ggWW_MVAWWBoundingUp  ->Write();
        histo_ggWW_MVAWWBoundingDown->Write();
      }
    outFileLimits->Close();

    double jeteff_E 	     = 1.02;
    double topXS_E  	     = 1.24;
    double wwXS_E            = 1.15;
    double wwXS_E_jet_extrap = 0.954;
    double XS_QCDscale_ggH[3];
    double UEPS  = HiggsSignalPSUESystematics(mH,nJetsType);
    XS_QCDscale_ggH[0] = HiggsSignalQCDScaleKappa("QCDscale_ggH",mH,nJetsType);
    XS_QCDscale_ggH[1] = HiggsSignalQCDScaleKappa("QCDscale_ggH1in",mH,nJetsType);
    XS_QCDscale_ggH[2] = HiggsSignalQCDScaleKappa("QCDscale_ggH2in",mH,nJetsType);
    if     (nJetsType == 1) {
      jeteff_E  	= 1.05;
      topXS_E   	= 1.09;
      wwXS_E    	= 1.30;
      wwXS_E_jet_extrap = 1.206;
    }
    else if(nJetsType == 2) {
      jeteff_E  	= 1.10;
      topXS_E   	= 1.50;
      wwXS_E    	= 1.30;
      wwXS_E_jet_extrap = 1.00;
    }
    for(int i=0; i<8; i++) if(nBgdAccDecays[i] < 0) nBgdAccDecays[i] = 0.0;
    for(int i=0; i<8; i++) if(nBgdCutDecays[i] < 0) nBgdCutDecays[i] = 0.0;
    for(int i=0; i<8; i++) if(nBgdMVADecays[i] < 0) nBgdMVADecays[i] = 0.0;
    for(int i=0; i<6; i++) if(nSigAcc[i] <= 0) nSigAcc[i] = 0.000;
    for(int i=0; i<6; i++) if(nSigCut[i] <= 0) nSigCut[i] = 0.000;
    for(int i=0; i<6; i++) if(nSigMVA[i] <= 0) nSigMVA[i] = 0.000;
    double yieldE[13],yield[13];
    int nData;
    int nTotalBins = 1; // histo_VH->GetNbinsX();
    //********SHAPE*******************
    char outputLimitsShape[200];
    for(int i=1; i<=nTotalBins; i++){
      if(nTotalBins != 1){
        if(histo_ttH   ->GetBinContent(i) > 0) { yieldE[0] = histo_ttH   ->GetBinError(i)/histo_ttH   ->GetBinContent(i);} else {yieldE[0] = 0.0;}
        if(histo_ZH    ->GetBinContent(i) > 0) { yieldE[1] = histo_ZH    ->GetBinError(i)/histo_ZH    ->GetBinContent(i);} else {yieldE[1] = 0.0;}
        if(histo_WH    ->GetBinContent(i) > 0) { yieldE[2] = histo_WH    ->GetBinError(i)/histo_WH    ->GetBinContent(i);} else {yieldE[2] = 0.0;}
        if(histo_qqH   ->GetBinContent(i) > 0) { yieldE[3] = histo_qqH   ->GetBinError(i)/histo_qqH   ->GetBinContent(i);} else {yieldE[3] = 0.0;}
        if(histo_ggH   ->GetBinContent(i) > 0) { yieldE[4] = histo_ggH   ->GetBinError(i)/histo_ggH   ->GetBinContent(i);} else {yieldE[4] = 0.0;}
        if(histo_qqWW  ->GetBinContent(i) > 0) { yieldE[5] = histo_qqWW  ->GetBinError(i)/histo_qqWW  ->GetBinContent(i);} else {yieldE[5] = 0.0;}
        if(histo_ggWW  ->GetBinContent(i) > 0) { yieldE[6] = histo_ggWW  ->GetBinError(i)/histo_ggWW  ->GetBinContent(i);} else {yieldE[6] = 0.0;}
        if(histo_VV    ->GetBinContent(i) > 0) { yieldE[7] = histo_VV    ->GetBinError(i)/histo_VV    ->GetBinContent(i);} else {yieldE[7] = 0.0;}
        if(histo_Top   ->GetBinContent(i) > 0) { yieldE[8] = histo_Top   ->GetBinError(i)/histo_Top   ->GetBinContent(i);} else {yieldE[8] = 0.0;}
        if(histo_Zjets ->GetBinContent(i) > 0) { yieldE[9] = histo_Zjets ->GetBinError(i)/histo_Zjets ->GetBinContent(i);} else {yieldE[9] = 0.0;}
        if(histo_Wjets ->GetBinContent(i) > 0) { yieldE[10]= histo_Wjets ->GetBinError(i)/histo_Wjets ->GetBinContent(i);} else {yieldE[10]= 0.0;}
        if(histo_Wgamma->GetBinContent(i) > 0) { yieldE[11]= histo_Wgamma->GetBinError(i)/histo_Wgamma->GetBinContent(i);} else {yieldE[11]= 0.0;}
        if(histo_Ztt   ->GetBinContent(i) > 0) { yieldE[12]= histo_Ztt   ->GetBinError(i)/histo_Ztt   ->GetBinContent(i);} else {yieldE[12]= 0.0;}
	yield[0] = histo_ttH   ->GetBinContent(i);
	yield[1] = histo_ZH    ->GetBinContent(i);
	yield[2] = histo_WH    ->GetBinContent(i);
	yield[3] = histo_qqH   ->GetBinContent(i);
	yield[4] = histo_ggH   ->GetBinContent(i);
	yield[5] = histo_qqWW  ->GetBinContent(i);
	yield[6] = histo_ggWW  ->GetBinContent(i);
	yield[7] = histo_VV    ->GetBinContent(i);
	yield[8] = histo_Top   ->GetBinContent(i);
	yield[9] = histo_Zjets ->GetBinContent(i);
	yield[10]= histo_Wjets ->GetBinContent(i);
	yield[11]= histo_Wgamma->GetBinContent(i);
	yield[12]= histo_Ztt   ->GetBinContent(i);
	nData    = (int)histo_Data->GetBinContent(i);
      }
      else {
 	if(nSigAcc[1] > 0) {yieldE[0] = nSigEAcc[1]/nSigAcc[1];} else {yieldE[0] = 0.0;}
 	if(nSigAcc[2] > 0) {yieldE[1] = nSigEAcc[2]/nSigAcc[2];} else {yieldE[1] = 0.0;}
 	if(nSigAcc[3] > 0) {yieldE[2] = nSigEAcc[3]/nSigAcc[3];} else {yieldE[2] = 0.0;}
 	if(nSigAcc[4] > 0) {yieldE[3] = nSigEAcc[4]/nSigAcc[4];} else {yieldE[3] = 0.0;}
 	if(nSigAcc[5] > 0) {yieldE[4] = nSigEAcc[5]/nSigAcc[5];} else {yieldE[4] = 0.0;}
 	yieldE[5] = nBgdEAccDecays[0];
 	yieldE[6] = nBgdEAccDecays[1];
 	yieldE[7] = nBgdEAccDecays[2];
 	yieldE[8] = nBgdEAccDecays[3];
 	yieldE[9] = nBgdEAccDecays[4];
 	yieldE[10]= nBgdEAccDecays[5];
 	yieldE[11]= nBgdEAccDecays[6];
 	yieldE[12]= nBgdEAccDecays[7];
	yield[0]  = nSigAcc[1];
	yield[1]  = nSigAcc[2];
	yield[2]  = nSigAcc[3];
	yield[3]  = nSigAcc[4];
	yield[4]  = nSigAcc[5];
	yield[5]  = nBgdAccDecays[0];
	yield[6]  = nBgdAccDecays[1];
	yield[7]  = nBgdAccDecays[2];
	yield[8]  = nBgdAccDecays[3];
	yield[9]  = nBgdAccDecays[4];
	yield[10]  = nBgdAccDecays[5];;
	yield[11]  = nBgdAccDecays[6];
	yield[12]  = nBgdAccDecays[7];
	nData     = (int)nDatAcc;
      }
      if(nTotalBins != 1){
        sprintf(outputLimitsShape,"output/histo_limits_%s_%dj_chan%d_mh%d_shape_bin%d.txt",outTag.Data(),nJetsType,wwDecay,mH,i);     
      }
      else {
        sprintf(outputLimitsShape,"output/histo_limits_%s_%dj_chan%d_mh%d_shape.txt",outTag.Data(),nJetsType,wwDecay,mH);     
      }
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",nData);
      if(nTotalBins == 1){
        if((wwDecay == 5 && useZjetsTemplates == true) || useWWTemplates == true)
          newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
        else
          newcardShape << Form("shapes *   *   %s  histo_$PROCESS\n",outputLimits);
        newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
      }
      newcardShape << Form("bin 1 1 1 1 1 1 1 1 1 1 1 1\n");
      newcardShape << Form("process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt\n");
      newcardShape << Form("process -3 -2 -1 0 1 2 3 4 5 6 7 8\n");
      newcardShape << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",yield[1],yield[2],yield[3],yield[4],yield[5],yield[6],yield[7],yield[8],yield[9],TMath::Max((double)yield[10],0.0),yield[11],yield[12]);
      newcardShape << Form("lumi                             lnN 1.045 1.045 1.045 1.045   -     -   1.045   -     -     -   1.045 1.045\n");				 
      newcardShape << Form("CMS_eff_m                        lnN 1.030 1.030 1.030 1.030 1.030 1.030 1.030   -     -     -   1.030 1.030\n");				 
      newcardShape << Form("CMS_eff_e                        lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040   -     -     -   1.040 1.040\n");				 
      newcardShape << Form("CMS_scale_m                      lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.015   -     -     -   1.015 1.015\n");				 
      newcardShape << Form("CMS_scale_e                      lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");
      newcardShape << Form("CMS_hww_met_resolution           lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");
      newcardShape << Form("CMS_scale_j                      lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -     -   %5.3f %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);	     
      newcardShape << Form("CMS_hww_fake_em       	     lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  \n");
      newcardShape << Form("UEPS 		             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",UEPS);
      newcardShape << Form("pdf_gg		             lnN   -     -     -   1.100   -   1.040   -     -     -     -     -     -  \n");
      newcardShape << Form("pdf_qqbar                        lnN 1.050 1.050 1.050   -   1.040   -   1.040   -     -     -   1.040 1.040\n");
      newcardShape << Form("QCDscale_ggH	             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[0]);  
      newcardShape << Form("QCDscale_ggH1in	             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[1]);  
      newcardShape << Form("QCDscale_ggH2in	             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[2]);  
      newcardShape << Form("QCDscale_qqH	             lnN   -     -   1.010   -     -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_VH 	             lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");			 
      newcardShape << Form("QCDscale_VV           	     lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_V           	     lnN   -     -     -     -     -     -     -     -     -     -     -   1.100\n");
      newcardShape << Form("QCDscale_ggVV	             lnN   -     -     -     -     -   1.500   -     -     -     -     -     -  \n");
      newcardShape << Form("CMS_QCDscale_WW_EXTRAP           lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",wwXS_E_jet_extrap);
      newcardShape << Form("QCDscale_ggH_ACEPT               lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_qqH_ACEPT               lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_VH_ACEPT                lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("CMS_hww_%1dj_ttbar               lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",nJetsType,topXS_E); 	
      newcardShape << Form("CMS_hww%s_%1dj_Z                 lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,ZXS_E[0]+1.0);			
      newcardShape << Form("CMS_hww_%1dj_WW                  lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",nJetsType,wwXS_E,wwXS_E);				
      if(wwDecay == 5 && useZjetsTemplates == true){
      newcardShape << Form("MVAZBounding                    shape   -     -     -     -    -     -     -     -    1.0    -     -     -   \n");			   
      }
      if(useWWTemplates == true){
      newcardShape << Form("MVAWWBounding                   shape   -     -     -     -   1.0   1.0    -     -     -     -     -     -   \n");			   
      }
      newcardShape << Form("CMS_hww%s_stat_%1dj_ZH_bin%d     lnN %5.3f   -     -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[1]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_WH_bin%d     lnN   -   %5.3f   -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[2]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_qqH_bin%d    lnN   -     -   %5.3f   -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[3]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_ggH_bin%d    lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[4]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_WW_bin%d     lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[5]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_ggWW_bin%d   lnN   -     -     -     -     -   %5.3f   -     -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[6]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_VV_bin%d     lnN   -     -     -     -     -     -   %5.3f   -     -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[7]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_ttbar_bin%d  lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",finalStateName,nJetsType,i,yieldE[8]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_Z_bin%d      lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,i,yieldE[9]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_Wjets_bin%d  lnN   -     -     -     -     -     -     -     -     -   %5.3f   -     -  \n",finalStateName,nJetsType,i,yieldE[10]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_Wgamma_bin%d lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",finalStateName,nJetsType,i,yieldE[11]+1.0);
      newcardShape << Form("CMS_hww%s_stat_%1dj_Ztt_bin%d    lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",finalStateName,nJetsType,i,yieldE[12]+1.0);
      newcardShape.close();
    }
    //********CUT*******************
    char outputLimitsCut[200];
    sprintf(outputLimitsCut,"output/histo_limits_%s_%dj_chan%d_mh%d_cut.txt",outTag.Data(),nJetsType,wwDecay,mH);     
    ofstream newcardCut;
    newcardCut.open(outputLimitsCut);
    newcardCut << Form("imax 1 number of channels\n");
    newcardCut << Form("jmax * number of background\n");
    newcardCut << Form("kmax * number of nuisance parameters\n");
    newcardCut << Form("Observation %d\n",(int)nDatCut);
    //newcardCut << Form("Observation %d\n",(int)(nBgdCutDecays[0]+nBgdCutDecays[1]+nBgdCutDecays[2]+nBgdCutDecays[3]+nBgdCutDecays[4]+nBgdCutDecays[5]+nBgdCutDecays[6]));
    newcardCut << Form("bin 1 1 1 1 1 1 1 1 1 1 1 1\n");
    newcardCut << Form("process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt\n");
    newcardCut << Form("process -3 -2 -1 0 1 2 3 4 5 6 7 8\n");
    newcardCut << Form("rate  %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigCut[2],nSigCut[3],nSigCut[4],nSigCut[5],nBgdCutDecays[0],nBgdCutDecays[1],nBgdCutDecays[2],nBgdCutDecays[3],nBgdCutDecays[4],TMath::Max((double)nBgdCutDecays[5],0.0),nBgdCutDecays[6],nBgdCutDecays[7]);
    newcardCut << Form("lumi                  	   lnN 1.045 1.045 1.045 1.045   -     -   1.045   -     -     -   1.045 1.045\n");			   
    newcardCut << Form("CMS_eff_m             	   lnN 1.030 1.030 1.030 1.030 1.030 1.030 1.030   -     -     -   1.030 1.030\n");			   
    newcardCut << Form("CMS_eff_e             	   lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040   -     -     -   1.040 1.040\n");			   
    newcardCut << Form("CMS_scale_m           	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.015   -     -     -   1.015 1.015\n");			   
    newcardCut << Form("CMS_scale_e           	   lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");			   
    newcardCut << Form("CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");			   
    newcardCut << Form("CMS_scale_j           	   lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -     -   %5.3f %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);		    
    newcardCut << Form("CMS_hww_fake_em       	   lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  \n");
    newcardCut << Form("UEPS 		           lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",UEPS);
    newcardCut << Form("pdf_gg                	   lnN   -     -     -   1.100   -   1.040   -     -     -     -     -     -  \n");
    newcardCut << Form("pdf_qqbar             	   lnN 1.050 1.050 1.050   -   1.040   -   1.040   -     -     -   1.040 1.040\n");
    newcardCut << Form("QCDscale_ggH          	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[0]);  
    newcardCut << Form("QCDscale_ggH1in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[1]);  
    newcardCut << Form("QCDscale_ggH2in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[2]);  
    newcardCut << Form("QCDscale_qqH          	   lnN   -     -   1.010   -     -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_VH           	   lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");		    
    newcardCut << Form("QCDscale_VV           	   lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_V           	   lnN   -     -     -     -     -     -     -     -     -     -     -   1.100\n");
    newcardCut << Form("QCDscale_ggVV         	   lnN   -     -     -     -     -   1.500   -     -     -     -     -     -  \n");
    newcardCut << Form("CMS_QCDscale_WW_EXTRAP     lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",wwXS_E_jet_extrap);
    newcardCut << Form("QCDscale_ggH_ACEPT    	   lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_qqH_ACEPT    	   lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_VH_ACEPT     	   lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("CMS_hww_%1dj_ttbar	   lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",nJetsType,topXS_E);	 
    newcardCut << Form("CMS_hww%s_%1dj_Z           lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,ZXS_E[1]+1.0);		 
    newcardCut << Form("CMS_hww_%1dj_WW            lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",nJetsType,wwXS_E,wwXS_E);		      
    newcardCut << Form("CMS_hww%s_stat_%1dj_ZH	   lnN %5.3f   -     -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[2]/TMath::Max((double)nSigCut[2],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_WH	   lnN   -   %5.3f   -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[3]/TMath::Max((double)nSigCut[3],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_qqH	   lnN   -     -   %5.3f   -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[4]/TMath::Max((double)nSigCut[4],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_ggH	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[5]/TMath::Max((double)nSigCut[5],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_WW	   lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[0]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_ggWW   lnN   -     -     -     -     -   %5.3f   -     -     -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[1]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_VV	   lnN   -     -     -     -     -     -   %5.3f   -     -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[2]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_ttbar  lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[3]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Z      lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[4]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Wjets  lnN   -     -     -     -     -     -     -     -     -   %5.3f   -     -  \n",finalStateName,nJetsType,nBgdECutDecays[5]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Wgamma lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",finalStateName,nJetsType,nBgdECutDecays[6]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Ztt    lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",finalStateName,nJetsType,nBgdECutDecays[7]+1.0);
    newcardCut.close();

    //***********MVA***************
    char outputLimitsMVA[200];
    sprintf(outputLimitsMVA,"output/histo_limits_%s_%dj_chan%d_mh%d_mva.txt",outTag.Data(),nJetsType,wwDecay,mH);     
    ofstream newcardMVA;
    newcardMVA.open(outputLimitsMVA);
    newcardMVA << Form("imax 1 number of channels\n");
    newcardMVA << Form("jmax * number of background\n");
    newcardMVA << Form("kmax * number of nuisance parameters\n");
    newcardMVA << Form("Observation %d\n",(int)nDatMVA);
    //newcardMVA << Form("Observation %d\n",(int)(nBgdMVADecays[0]+nBgdMVADecays[1]+nBgdMVADecays[2]+nBgdMVADecays[3]+nBgdMVADecays[4]+nBgdMVADecays[5]+nBgdMVADecays[6]));
    newcardMVA << Form("bin 1 1 1 1 1 1 1 1 1 1 1 1\n");
    newcardMVA << Form("process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt\n");
    newcardMVA << Form("process -3 -2 -1 0 1 2 3 4 5 6 7 8\n");
    newcardMVA << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigMVA[2],nSigMVA[3],nSigMVA[4],nSigMVA[5],nBgdMVADecays[0],nBgdMVADecays[1],nBgdMVADecays[2],nBgdMVADecays[3],nBgdMVADecays[4],TMath::Max((double)nBgdMVADecays[5],0.0),nBgdMVADecays[6],nBgdMVADecays[7]);
    newcardMVA << Form("lumi                  	   lnN 1.045 1.045 1.045 1.045   -     -   1.045   -     -     -   1.045 1.045\n");
    newcardMVA << Form("CMS_eff_m             	   lnN 1.030 1.030 1.030 1.030 1.030 1.030 1.030   -     -     -   1.030 1.030\n");			    
    newcardMVA << Form("CMS_eff_e             	   lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040   -     -     -   1.040 1.040\n");			    
    newcardMVA << Form("CMS_scale_m           	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.015   -     -     -   1.015 1.015\n");			    
    newcardMVA << Form("CMS_scale_e           	   lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");			    
    newcardMVA << Form("CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");			    
    newcardMVA << Form("CMS_scale_j           	   lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -     -   %5.3f %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);	      
    newcardMVA << Form("CMS_hww_fake_em       	   lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  \n");
    newcardMVA << Form("UEPS 		           lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",UEPS);
    newcardMVA << Form("pdf_gg                	   lnN   -     -     -   1.100   -   1.040   -     -     -     -     -     -  \n");
    newcardMVA << Form("pdf_qqbar             	   lnN 1.050 1.050 1.050   -   1.040   -   1.040   -     -     -   1.040 1.040\n");
    newcardMVA << Form("QCDscale_ggH          	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[0]);  
    newcardMVA << Form("QCDscale_ggH1in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[1]);  
    newcardMVA << Form("QCDscale_ggH2in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[2]);  
    newcardMVA << Form("QCDscale_qqH          	   lnN   -     -   1.010   -     -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_VH           	   lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");		    
    newcardMVA << Form("QCDscale_VV           	   lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_V           	   lnN   -     -     -     -     -     -   1.040   -     -     -     -   1.100\n");
    newcardMVA << Form("QCDscale_ggVV         	   lnN   -     -     -     -     -   1.500   -     -     -     -     -     -  \n");
    newcardMVA << Form("CMS_QCDscale_WW_EXTRAP     lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",wwXS_E_jet_extrap);
    newcardMVA << Form("QCDscale_ggH_ACEPT    	   lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_qqH_ACEPT    	   lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_VH_ACEPT     	   lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("CMS_hww_%1dj_ttbar	   lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",nJetsType,topXS_E);	 
    newcardMVA << Form("CMS_hww%s_%1dj_Z           lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,ZXS_E[2]+1.0);		 
    newcardMVA << Form("CMS_hww_%1dj_WW            lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",nJetsType,wwXS_E,wwXS_E);		      
    newcardMVA << Form("CMS_hww%s_stat_%1dj_ZH	   lnN %5.3f   -     -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigEMVA[2]/TMath::Max((double)nSigMVA[2],0.00001)+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_WH	   lnN   -   %5.3f   -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigEMVA[3]/TMath::Max((double)nSigMVA[3],0.00001)+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_qqH	   lnN   -     -   %5.3f   -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigEMVA[4]/TMath::Max((double)nSigMVA[4],0.00001)+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_ggH	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigEMVA[5]/TMath::Max((double)nSigMVA[5],0.00001)+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_WW	   lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nBgdEMVADecays[0]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_ggWW   lnN   -     -     -     -     -   %5.3f   -     -     -     -     -     -  \n",finalStateName,nJetsType,nBgdEMVADecays[1]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_VV	   lnN   -     -     -     -     -     -   %5.3f   -     -     -     -     -  \n",finalStateName,nJetsType,nBgdEMVADecays[2]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_ttbar  lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",finalStateName,nJetsType,nBgdEMVADecays[3]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_Z      lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,nBgdEMVADecays[4]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_Wjets  lnN   -     -     -     -     -     -     -     -     -   %5.3f   -     -  \n",finalStateName,nJetsType,nBgdEMVADecays[5]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_Wgamma lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",finalStateName,nJetsType,nBgdEMVADecays[6]+1.0);
    newcardMVA << Form("CMS_hww%s_stat_%1dj_Ztt    lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",finalStateName,nJetsType,nBgdEMVADecays[7]+1.0);
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

//------------------------------------------------------------------------------
// DeltaPhi
//------------------------------------------------------------------------------
double DeltaPhi(double phi1, double phi2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].
  double dphi = TMath::Abs(phi1-phi2);
  while (dphi>TMath::Pi())
    dphi = TMath::Abs(dphi - TMath::TwoPi());
  return(dphi);
}
