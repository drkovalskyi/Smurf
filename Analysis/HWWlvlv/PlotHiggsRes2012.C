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
#include "/home/ceballos/releases/CMSSW_5_2_3_patch3/src/Smurf/Core/SmurfTree.h"
#include "factors.h"
#include "HiggsQCDScaleSystematics_8TeV.h"
#include "PSUESystematics_8TeV.h"
#include "PDFgHHSystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_2_3_patch3/src/Smurf/Core/LeptonScaleLookup.h"
#include "DYBkgScaleFactors_8TeV.h"
#include "TopBkgScaleFactors_8TeV.h"
#include "TopVBFBkgScaleFactors_8TeV.h"
#include "WWBkgScaleFactors_8TeV.h"
#include "OtherBkgScaleFactors_8TeV.h"
#include "HWWCuts.h"
#include "HiggsSM4Systematics_8TeV.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void plotHistsInPad(TH1D* h1, TH1D* h2);
void setHist(TH1D* h, int color, int style);
void setPair(TH1D* h1, TH1D* h2);
void setGraph(TGraphErrors* g, int color, int marker);
TGraphErrors* makeSignificanceCurve(TH1D* sig, TH1D* bgd, TH1D* dat, const char* name);
TGraphErrors* makeGraphFromHists   (TH1D* sig, TH1D* bgd, const char* name);

int    verboseLevel =   0;
const double sigmaB = 0.35;

void PlotHiggsRes2012
(
 UInt_t  nJetsType   	 = 0,
 UInt_t  mH      	 = 300,
 TString outTag          = "default",
 TString sigInputFile    = "data/inputNtuple-train-data-standard-H150_WW_2l.root",
 TString bgdInputFile    = "data/inputNtuple-data-standard-HBCK_WW_2l-train.root",
 TString datInputFile    = "data/data-train.root",
 Int_t   wwDecay         = 0,
 bool fillInfoNote       = false,
 TString systInputFile   = "",
 int period              = 2,
 int  TheVerboseLevel    = 0
 )
{
  bool wwPresel = false; double ptlMin = 10.0;
  if(mH == 0) {wwPresel = true; mH = 160; ptlMin = 20.0;}
  Bool_t useDYMVA = kFALSE;
  if(wwPresel == false && mH <= 140) useDYMVA = kTRUE;

  int category = 0;
  if(period > 9) {period = period - 10; category = 1;}

  bool signalInjection = false;

  bool isFermioPhobic = false;
  bool isSM4          = false;

  verboseLevel = TheVerboseLevel;
  bool useZjetsTemplates   = true;
  bool useWWTemplates      = true;
  bool useStatTemplates    = true;
  bool useExpTemplates     = true;
  bool useJESTemplates     = true;
  bool useWJetsTemplates   = true;
  bool useWJetsMCTemplates = false; // default is false in 2012 analysis
  bool useTopTemplates     = true;
  bool useggHTemplates     = false; // default is false in 2012 analysis
  bool useWgammaTemplates  = false; // this is intentional
  if(nJetsType != 0 && nJetsType != 1){
    useZjetsTemplates	= false;
    useWWTemplates	= false;
    useStatTemplates	= false;
    useExpTemplates	= false;
    useJESTemplates	= false;
    useWJetsTemplates	= false;
    useWJetsMCTemplates = false;
    useTopTemplates	= false;
    useggHTemplates	= false;
    useWgammaTemplates  = false;
  }
  int rebinMVAHist        = 10;
  unsigned int rebinSmurf = 4;
  TH1D *histSmurfRebin;

  int binVarA = 10;
  int binVarB = 7;
  bool is2DAna = false;
  if(binVarA > 0 && binVarB > 0) is2DAna = true;

  bool makeZjetsTemplates = false;
  if(makeZjetsTemplates == true) useZjetsTemplates = false;

  bool useAlternativeStatTemplates = false;

  if(wwDecay != 0 && wwDecay != 3 && wwDecay != 5) useZjetsTemplates = false;

  TString sigFile1 = sigInputFile;
  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  char replace[200];
  TString c1Name = sigFile1;
  c1Name.ReplaceAll("data2012","pic");
  sprintf(replace,"_%d_8TeV.eps",wwDecay);
  c1Name.ReplaceAll(".root",replace);

  TString c2Name = sigFile1;
  c2Name.ReplaceAll("data2012","pic");
  sprintf(replace,"_%d_effvsbkg_8TeV.eps",wwDecay);
  c2Name.ReplaceAll(".root",replace);

  TString output = sigFile1;
  output.ReplaceAll("data2012/","histo_tmva_");
  sprintf(replace,"_chan%d_8TeV.root",wwDecay+10*category);
  output.ReplaceAll(".root",replace);

  unsigned int patternTopTag = SmurfTree::TopTag;

  int channel = HiggsMassIndex(mH)-1;

  if(channel == -1) return;

  float dilmass_cut = DileptonMassPreselectionCut(mH);
  if(wwPresel == true) dilmass_cut = 99999.;
  float mtUpperCut = mH;
  float mtLowerCut = 80;
  if     (is2DAna == true && mH <= 250) {dilmass_cut = 200; mtLowerCut = 70; mtUpperCut = 250; useZjetsTemplates = false; useExpTemplates = false;}
  else if(is2DAna == true && mH >  250) {dilmass_cut = mH;  mtLowerCut = 70; mtUpperCut = mH;  useZjetsTemplates = false; useExpTemplates = false;}

  char finalStateName[10];
  sprintf(finalStateName,"ll");
  if	 (wwDecay == 0) sprintf(finalStateName,"mm");
  else if(wwDecay == 1) sprintf(finalStateName,"me");
  else if(wwDecay == 2) sprintf(finalStateName,"em");
  else if(wwDecay == 3) sprintf(finalStateName,"ee");
  else if(wwDecay == 5) sprintf(finalStateName,"sf");
  else if(wwDecay == 6) sprintf(finalStateName,"of");
  if(category == 1) sprintf(finalStateName,"%slt",finalStateName);

  //----------------------------------------------------------------------------
  // These are used to compute the DY Bkg systematics uncertainties
  // DYXS, VVXS give the normalization for the DY Bkg and the WW,WZ Bkg's
  // ZXS_E is the systematic uncertainty in the normalization
  // The indices parameterize the different versions of the analysis:
  // [0] : MVA Shape Analysis
  // [1] : Cut-Based Analysis
  // [2] : MVA Cut Analysis
  //----------------------------------------------------------------------------
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

  TChain *chsystInputFile = new TChain("tree");
  chsystInputFile->Add(systInputFile);
  TTree *treeSyst = (TTree*) chsystInputFile;

  // k-factor ggH systematics
  // 0: cc
  // 1: --
  // 2: -c
  // 3: c-
  // 4: c+
  // 5: +c
  // 6: ++
  // 7: -+
  // 8: +-

  TString effPath      = "";
  TString fakePath     = "";
  TString puPath       = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  double scaleFactorLum = 2.121;
  if	 (period == 0){ // Full2012-Summer12-V9-3500ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    scaleFactorLum     = 3.553;minRun =      0;maxRun = 999999;
  }
  else if(period == -1){ // Full2012-Summer12-V7
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-02_V3/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-02_V3/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X//auxiliar/puWeights_Summer12.root";
    scaleFactorLum     = 0.818;minRun =      0;maxRun = 999999;
  }
  else if(period == 1){ // Full2012-Summer12-V9-5000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-06_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_5000ipb_71mb.root";
    scaleFactorLum     = 5.098;minRun =      0;maxRun = 999999;
  }
  else if(period == 2){ // Full2012-Summer12-V9-5000ipb-NewTraining
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-06_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_5000ipb_71mb.root";
    scaleFactorLum     = 5.098;minRun =      0;maxRun = 999999;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  //----------------------------------------------------------------------------
  // Lepton Efficiency Scale Factors, Fake Rates, 
  // Trigger Efficiencies
  // Pileup Weights
  //----------------------------------------------------------------------------
  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open(fakePath.Data());
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open(fakePath.Data());
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  LeptonScaleLookup trigLookup(effPath.Data());

  TFile *fPU2012File = TFile::Open(puPath.Data());
  TH1D *fhDPU2012 = (TH1D*)(fPU2012File->Get("puWeights"));
  assert(fhDPU2012);
  fhDPU2012->SetDirectory(0);
  delete fPU2012File;

  //Fake rate systematics
  TFile *fLeptonFRFileSyst = TFile::Open(fakePath.Data());
  TH2D *fhDFRMuSyst = (TH2D*)(fLeptonFRFileSyst->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
  TH2D *fhDFRElSyst = (TH2D*)(fLeptonFRFileSyst->Get("ElectronFakeRate_V4_ptThreshold50_PtEta"));
  assert(fhDFRMuSyst);
  assert(fhDFRElSyst);
  fhDFRMuSyst->SetDirectory(0);
  fhDFRElSyst->SetDirectory(0);
  fLeptonFRFileSyst->Close();
  delete fLeptonFRFileSyst;

  //----------------------------------------------------------------------------
  // ggH pT spectrum Weights
  //----------------------------------------------------------------------------
  int newMH = mH;
  if(newMH == 110) newMH = 115; // there is no correction for mh=110!
  if(newMH >  600) newMH = 600; // there is no correction for mh>600!
  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/data/Winter11_4700ipb/auxiliar/ggHWW_KFactors_PowhegToHQT_WithAdditionalMassPoints.root");
  TH1D *HiggsPtKFactor,*HiggsPtKFactorSyst[8];
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", newMH);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  for(int i=0; i<8; i++) HiggsPtKFactorSyst[i] = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%d_QCDscaleSys%1d",newMH,i+1)));
  if (HiggsPtKFactor) {
    HiggsPtKFactor->SetDirectory(0);
    for(int i=0; i<8; i++) if(HiggsPtKFactorSyst[i]) HiggsPtKFactorSyst[i]->SetDirectory(0);
  }
  assert(HiggsPtKFactor);
  //for(int i=0; i<8; i++) assert(HiggsPtKFactorSyst[i]); // only 1 and 6 are available currently (0 and 5 in c++ version)
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;

  //----------------------------------------------------------------------------
  // Load MVA distribution histograms for the DY bkg systematics from file
  // These are produced separately and saved into a file, and loaded when 
  // it is needed.
  // The nominal MVA shape is taken from the MET sideband 
  // (20 < Met < norminal met cut)
  //----------------------------------------------------------------------------
  TH1D *hDZjetsMCTemplate,*hDZjetsDATemplate;
  if(useZjetsTemplates == true) printf("***********useZjetsTemplates = true***************\n");
  else                          printf("***********useZjetsTemplates = false***************\n");
  if(useZjetsTemplates == true){
    TFile *fZjetsTemplatesFile;
    if     (period == 0 || period == -1) fZjetsTemplatesFile = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/histo_Zjets_Templates_3500ipb.root");
    else if(period == 1                ) fZjetsTemplatesFile = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/histo_Zjets_Templates_5000ipb.root");
    else if(period == 2                ) fZjetsTemplatesFile = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/histo_Zjets_Templates_5000ipb_PostICHEP.root");
    else                            assert(0);
    char ZjetsHistName[100];
    sprintf(ZjetsHistName, "hDMCZjets%d_%d", mH,TMath::Min((int)nJetsType,1));
    hDZjetsMCTemplate = (TH1D*)(fZjetsTemplatesFile->Get(ZjetsHistName));
    if (hDZjetsMCTemplate) {
      hDZjetsMCTemplate->SetDirectory(0);
    }
    assert(hDZjetsMCTemplate);
    sprintf(ZjetsHistName, "hDDAZjets%d_%d", mH,TMath::Min((int)nJetsType,1));
    hDZjetsDATemplate = (TH1D*)(fZjetsTemplatesFile->Get(ZjetsHistName));
    if (hDZjetsDATemplate) {
      hDZjetsDATemplate->SetDirectory(0);
    }
    assert(hDZjetsDATemplate);
    fZjetsTemplatesFile->Close();
    delete fZjetsTemplatesFile;
    if(hDZjetsMCTemplate->GetSumOfWeights() != 1.0){
      printf("hDZjetsMCTemplate(%d) norm = %f --> normalizing to 1\n",mH,hDZjetsMCTemplate->GetSumOfWeights());
      hDZjetsMCTemplate->Scale(1./hDZjetsMCTemplate->GetSumOfWeights());
    }
    if(hDZjetsDATemplate->GetSumOfWeights() != 1.0){
      printf("hDZjetsDATemplate(%d) norm = %f --> normalizing to 1\n",mH,hDZjetsDATemplate->GetSumOfWeights());
      hDZjetsDATemplate->Scale(1./hDZjetsDATemplate->GetSumOfWeights());
    }
  }

  //----------------------------------------------------------------------------
  // Define m_Higgs dependant cuts
  //----------------------------------------------------------------------------
  double theCutMassHigh	     = cutMassHigh (mH);
  double theCutPtMaxLow	     = cutPtMaxLow (mH);
  double theCutDeltaphilHigh = cutDeltaphiHigh (mH);
  double theCutMTLow         = cutMTLow (mH);
  double theCutMTHigh        = cutMTHigh (mH);

  //----------------------------------------------------------------------------
  // Define Histogram ranges
  //----------------------------------------------------------------------------
  const int nHist = 6;
  int    nBinHis        = 200;
  double minHis[nHist]  = {-1.0, -1.0, -1.0, -0.0, -1.0,  0.0};
  double maxHis[nHist]  = { 1.0,  1.0,  1.0,  1.0,  1.0,200.0};

  if(is2DAna == true) {nBinHis = binVarA*binVarB; rebinMVAHist = 1; rebinSmurf = 1;}

  //----------------------------------------------------------------------------
  // Define MVA output histograms
  //----------------------------------------------------------------------------
  TH1D* sigMVA[nHist][6];
  TH1D* bgdMVA[nHist];
  TH1D* datMVA[nHist];
  const int nChan = 8;
  TH1D* bgdMVADecays[nHist][nChan];
  for(int i=0; i<nHist; i++) {
    for(int j=0; j<6; j++){
      sigMVA[i][j] = new TH1D(Form("sigMVA_%d_%d",i,j), Form("sigMVA_%d_%d",i,j), nBinHis, minHis[i], maxHis[i]);
      sigMVA[i][j]->Sumw2();
    }
    bgdMVA[i] = new TH1D(Form("bgdMVA_%d",i), Form("bgdMVA_%d",i), nBinHis, minHis[i], maxHis[i]);
    datMVA[i] = new TH1D(Form("datMVA_%d",i), Form("datMVA_%d",i), nBinHis, minHis[i], maxHis[i]);
    bgdMVA[i]->Sumw2();
    datMVA[i]->Sumw2();
    for(int j=0; j<nChan; j++) {
      bgdMVADecays[i][j] = new TH1D(Form("bgdMVADecays_%d_%d",i,j), Form("bgdMVADecays_%d_%d",i,j), nBinHis, minHis[i], maxHis[i]);
      bgdMVADecays[i][j]->Sumw2();
    }
  }

  TH1D* histos = new TH1D("histos", "histos", nBinHis, minHis[1], maxHis[1]); histos->Sumw2(); histos->Scale(0.0);
  TH1D* histo0 = new TH1D("histo0", "histo0", nBinHis, minHis[1], maxHis[1]); histo0->Sumw2(); histo0->Scale(0.0);
  TH1D* histo1 = new TH1D("histo1", "histo1", nBinHis, minHis[1], maxHis[1]); histo1->Sumw2(); histo1->Scale(0.0);
  TH1D* histo2 = new TH1D("histo2", "histo2", nBinHis, minHis[1], maxHis[1]); histo2->Sumw2(); histo2->Scale(0.0);
  TH1D* histo3 = new TH1D("histo3", "histo3", nBinHis, minHis[1], maxHis[1]); histo3->Sumw2(); histo3->Scale(0.0);
  TH1D* histo4 = new TH1D("histo4", "histo4", nBinHis, minHis[1], maxHis[1]); histo4->Sumw2(); histo4->Scale(0.0);
  TH1D* histo5 = new TH1D("histo5", "histo5", nBinHis, minHis[1], maxHis[1]); histo5->Sumw2(); histo5->Scale(0.0);

  //----------------------------------------------------------------------------
  // Define MVA output systematics histograms
  //----------------------------------------------------------------------------
  TH1D* histo_Zjets_CMS_MVAZBounding = new TH1D( Form("histo_Zjets_CMS_hww%s_%1dj_MVAZBounding",finalStateName,nJetsType), Form("histo_Zjets_CMS_hww%s_%1dj_MVAZBounding",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Zjets_CMS_MVAZBounding->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZBoundingUp = new TH1D( Form("histo_Zjets_CMS_hww%s_%1dj_MVAZBoundingUp",finalStateName,nJetsType), Form("histo_Zjets_CMS_hww%s_%1dj_MVAZBoundingUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Zjets_CMS_MVAZBoundingUp->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZBoundingDown = new TH1D( Form("histo_Zjets_CMS_hww%s_%1dj_MVAZBoundingDown",finalStateName,nJetsType), Form("histo_Zjets_CMS_hww%s_%1dj_MVAZBoundingDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Zjets_CMS_MVAZBoundingDown->Sumw2();
  TH1D* histoVV = new TH1D("histoVV", "histoVV", nBinHis, minHis[1], maxHis[1]);  histoVV->Sumw2();

  TH1D* histo_Wjets_CMS_MVAWBoundingUp = new TH1D( Form("histo_Wjets_CMS_hww_MVAWBoundingUp"), Form("histo_Wjets_CMS_hww_MVAWBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Wjets_CMS_MVAWBoundingUp->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWBoundingDown = new TH1D( Form("histo_Wjets_CMS_hww_MVAWBoundingDown"), Form("histo_Wjets_CMS_hww_MVAWBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Wjets_CMS_MVAWBoundingDown->Sumw2();

  TH1D* histo_Wjets_CMS_hww_MVAWMCBoundingUp = new TH1D( Form("histo_Wjets_CMS_hww_MVAWMCBoundingUp"), Form("histo_Wjets_CMS_hww_MVAWMCBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Wjets_CMS_hww_MVAWMCBoundingUp->Sumw2();
  TH1D* histo_Wjets_CMS_hww_MVAWMCBoundingDown = new TH1D( Form("histo_Wjets_CMS_hww_MVAWMCBoundingDown"), Form("histo_Wjets_CMS_hww_MVAWMCBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Wjets_CMS_hww_MVAWMCBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAWWBoundingUp = new TH1D( Form("histo_qqWW_CMS_hww_MVAWWBoundingUp"), Form("histo_qqWW_CMS_hww_MVAWWBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_MVAWWBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWWBoundingDown = new TH1D( Form("histo_qqWW_CMS_hww_MVAWWBoundingDown"), Form("histo_qqWW_CMS_hww_MVAWWBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_MVAWWBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAWWNLOBoundingUp = new TH1D( Form("histo_qqWW_CMS_hww_MVAWWNLOBoundingUp"), Form("histo_qqWW_CMS_hww_MVAWWNLOBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_MVAWWNLOBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWWNLOBoundingDown = new TH1D( Form("histo_qqWW_CMS_hww_MVAWWNLOBoundingDown"), Form("histo_qqWW_CMS_hww_MVAWWNLOBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_MVAWWNLOBoundingDown->Sumw2();

  TH1D* histo_Top_CMS_MVATopBoundingUp = new TH1D( Form("histo_Top_CMS_hww_MVATopBoundingUp"), Form("histo_Top_CMS_hww_MVATopBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_MVATopBoundingUp->Sumw2();
  TH1D* histo_Top_CMS_MVATopBoundingDown = new TH1D( Form("histo_Top_CMS_hww_MVATopBoundingDown"), Form("histo_Top_CMS_hww_MVATopBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_MVATopBoundingDown->Sumw2();

  TH1D* histo_Wgamma = new TH1D( Form("histo_Wgamma"), Form("histo_Wgamma"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma->Sumw2();

  TH1D* histo_ttH_CMS_MVAttHStatBounding_8TeVUp = new TH1D( Form("histo_ttH_CMS_hww%s_%1dj_MVAttHStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_ttH_CMS_hww%s_%1dj_MVAttHStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_MVAttHStatBounding_8TeVUp->Sumw2();
  TH1D* histo_ttH_CMS_MVAttHStatBounding_8TeVDown = new TH1D( Form("histo_ttH_CMS_hww%s_%1dj_MVAttHStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_ttH_CMS_hww%s_%1dj_MVAttHStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_MVAttHStatBounding_8TeVDown->Sumw2();
  TH1D* histo_ZH_CMS_MVAZHStatBounding_8TeVUp = new TH1D( Form("histo_ZH_CMS_hww%s_%1dj_MVAZHStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_ZH_CMS_hww%s_%1dj_MVAZHStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_MVAZHStatBounding_8TeVUp->Sumw2();
  TH1D* histo_ZH_CMS_MVAZHStatBounding_8TeVDown = new TH1D( Form("histo_ZH_CMS_hww%s_%1dj_MVAZHStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_ZH_CMS_hww%s_%1dj_MVAZHStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_MVAZHStatBounding_8TeVDown->Sumw2();
  TH1D* histo_WH_CMS_MVAWHStatBounding_8TeVUp = new TH1D( Form("histo_WH_CMS_hww%s_%1dj_MVAWHStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_WH_CMS_hww%s_%1dj_MVAWHStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_MVAWHStatBounding_8TeVUp->Sumw2();
  TH1D* histo_WH_CMS_MVAWHStatBounding_8TeVDown = new TH1D( Form("histo_WH_CMS_hww%s_%1dj_MVAWHStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_WH_CMS_hww%s_%1dj_MVAWHStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_MVAWHStatBounding_8TeVDown->Sumw2();
  TH1D* histo_qqH_CMS_MVAqqHStatBounding_8TeVUp = new TH1D( Form("histo_qqH_CMS_hww%s_%1dj_MVAqqHStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_qqH_CMS_hww%s_%1dj_MVAqqHStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_MVAqqHStatBounding_8TeVUp->Sumw2();
  TH1D* histo_qqH_CMS_MVAqqHStatBounding_8TeVDown = new TH1D( Form("histo_qqH_CMS_hww%s_%1dj_MVAqqHStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_qqH_CMS_hww%s_%1dj_MVAqqHStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_MVAqqHStatBounding_8TeVDown->Sumw2();
  TH1D* histo_ggH_CMS_MVAggHStatBounding_8TeVUp = new TH1D( Form("histo_ggH_CMS_hww%s_%1dj_MVAggHStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_ggH_CMS_hww%s_%1dj_MVAggHStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_MVAggHStatBounding_8TeVUp->Sumw2();
  TH1D* histo_ggH_CMS_MVAggHStatBounding_8TeVDown = new TH1D( Form("histo_ggH_CMS_hww%s_%1dj_MVAggHStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_ggH_CMS_hww%s_%1dj_MVAggHStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_MVAggHStatBounding_8TeVDown->Sumw2();
  TH1D* histo_qqWW_CMS_MVAqqWWStatBounding_8TeVUp = new TH1D( Form("histo_qqWW_CMS_hww%s_%1dj_MVAqqWWStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_qqWW_CMS_hww%s_%1dj_MVAqqWWStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_MVAqqWWStatBounding_8TeVUp->Sumw2();
  TH1D* histo_qqWW_CMS_MVAqqWWStatBounding_8TeVDown = new TH1D( Form("histo_qqWW_CMS_hww%s_%1dj_MVAqqWWStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_qqWW_CMS_hww%s_%1dj_MVAqqWWStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_MVAqqWWStatBounding_8TeVDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBounding_8TeVUp = new TH1D( Form("histo_ggWW_CMS_hww%s_%1dj_MVAggWWStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_ggWW_CMS_hww%s_%1dj_MVAggWWStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_MVAggWWStatBounding_8TeVUp->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBounding_8TeVDown = new TH1D( Form("histo_ggWW_CMS_hww%s_%1dj_MVAggWWStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_ggWW_CMS_hww%s_%1dj_MVAggWWStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_MVAggWWStatBounding_8TeVDown->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBounding_8TeVUp = new TH1D( Form("histo_VV_CMS_hww%s_%1dj_MVAVVStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_VV_CMS_hww%s_%1dj_MVAVVStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_MVAVVStatBounding_8TeVUp->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBounding_8TeVDown = new TH1D( Form("histo_VV_CMS_hww%s_%1dj_MVAVVStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_VV_CMS_hww%s_%1dj_MVAVVStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_MVAVVStatBounding_8TeVDown->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBounding_8TeVUp = new TH1D( Form("histo_Top_CMS_hww%s_%1dj_MVATopStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_Top_CMS_hww%s_%1dj_MVATopStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_MVATopStatBounding_8TeVUp->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBounding_8TeVDown = new TH1D( Form("histo_Top_CMS_hww%s_%1dj_MVATopStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_Top_CMS_hww%s_%1dj_MVATopStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_MVATopStatBounding_8TeVDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBounding_8TeVUp = new TH1D( Form("histo_Zjets_CMS_hww%s_%1dj_MVAZjetsStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_Zjets_CMS_hww%s_%1dj_MVAZjetsStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Zjets_CMS_MVAZjetsStatBounding_8TeVUp->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBounding_8TeVDown = new TH1D( Form("histo_Zjets_CMS_hww%s_%1dj_MVAZjetsStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_Zjets_CMS_hww%s_%1dj_MVAZjetsStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Zjets_CMS_MVAZjetsStatBounding_8TeVDown->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp = new TH1D( Form("histo_Wjets_CMS_hww%s_%1dj_MVAWjetsStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_Wjets_CMS_hww%s_%1dj_MVAWjetsStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown = new TH1D( Form("histo_Wjets_CMS_hww%s_%1dj_MVAWjetsStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_Wjets_CMS_hww%s_%1dj_MVAWjetsStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp = new TH1D( Form("histo_Wgamma_CMS_hww%s_%1dj_MVAWgammaStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_Wgamma_CMS_hww%s_%1dj_MVAWgammaStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown = new TH1D( Form("histo_Wgamma_CMS_hww%s_%1dj_MVAWgammaStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_Wgamma_CMS_hww%s_%1dj_MVAWgammaStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown->Sumw2();
  TH1D* histo_Ztt_CMS_MVAZttStatBounding_8TeVUp = new TH1D( Form("histo_Ztt_CMS_hww%s_%1dj_MVAZttStatBounding_8TeVUp",finalStateName,nJetsType), Form("histo_Ztt_CMS_hww%s_%1dj_MVAZttStatBounding_8TeVUp",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_MVAZttStatBounding_8TeVUp->Sumw2();
  TH1D* histo_Ztt_CMS_MVAZttStatBounding_8TeVDown = new TH1D( Form("histo_Ztt_CMS_hww%s_%1dj_MVAZttStatBounding_8TeVDown",finalStateName,nJetsType), Form("histo_Ztt_CMS_hww%s_%1dj_MVAZttStatBounding_8TeVDown",finalStateName,nJetsType), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_MVAZttStatBounding_8TeVDown->Sumw2();

  TH1D* histo_ttH_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_ttH_CMS_hww_MVALepEffBoundingUp"), Form("histo_ttH_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_ttH_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_ttH_CMS_hww_MVALepEffBoundingDown"), Form("histo_ttH_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_ZH_CMS_hww_MVALepEffBoundingUp"), Form("histo_ZH_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_ZH_CMS_hww_MVALepEffBoundingDown"), Form("histo_ZH_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_WH_CMS_hww_MVALepEffBoundingUp"), Form("histo_WH_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_WH_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_WH_CMS_hww_MVALepEffBoundingDown"), Form("histo_WH_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_qqH_CMS_hww_MVALepEffBoundingUp"), Form("histo_qqH_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_qqH_CMS_hww_MVALepEffBoundingDown"), Form("histo_qqH_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_ggH_CMS_hww_MVALepEffBoundingUp"), Form("histo_ggH_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_ggH_CMS_hww_MVALepEffBoundingDown"), Form("histo_ggH_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_qqWW_CMS_hww_MVALepEffBoundingUp"), Form("histo_qqWW_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_qqWW_CMS_hww_MVALepEffBoundingDown"), Form("histo_qqWW_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_ggWW_CMS_hww_MVALepEffBoundingUp"), Form("histo_ggWW_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_ggWW_CMS_hww_MVALepEffBoundingDown"), Form("histo_ggWW_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_VV_CMS_hww_MVALepEffBoundingUp"), Form("histo_VV_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_VV_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_VV_CMS_hww_MVALepEffBoundingDown"), Form("histo_VV_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_Wgamma_CMS_hww_MVALepEffBoundingUp"), Form("histo_Wgamma_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_Wgamma_CMS_hww_MVALepEffBoundingDown"), Form("histo_Wgamma_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVALepEffBoundingUp = new TH1D( Form("histo_Ztt_CMS_hww_MVALepEffBoundingUp"), Form("histo_Ztt_CMS_hww_MVALepEffBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVALepEffBoundingUp->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVALepEffBoundingDown = new TH1D( Form("histo_Ztt_CMS_hww_MVALepEffBoundingDown"), Form("histo_Ztt_CMS_hww_MVALepEffBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVALepEffBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_ttH_CMS_hww_MVALepResBoundingUp"), Form("histo_ttH_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_ttH_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_ttH_CMS_hww_MVALepResBoundingDown"), Form("histo_ttH_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_ZH_CMS_hww_MVALepResBoundingUp"), Form("histo_ZH_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_ZH_CMS_hww_MVALepResBoundingDown"), Form("histo_ZH_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_WH_CMS_hww_MVALepResBoundingUp"), Form("histo_WH_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_WH_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_WH_CMS_hww_MVALepResBoundingDown"), Form("histo_WH_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_qqH_CMS_hww_MVALepResBoundingUp"), Form("histo_qqH_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_qqH_CMS_hww_MVALepResBoundingDown"), Form("histo_qqH_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_ggH_CMS_hww_MVALepResBoundingUp"), Form("histo_ggH_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_ggH_CMS_hww_MVALepResBoundingDown"), Form("histo_ggH_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_qqWW_CMS_hww_MVALepResBoundingUp"), Form("histo_qqWW_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_qqWW_CMS_hww_MVALepResBoundingDown"), Form("histo_qqWW_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_ggWW_CMS_hww_MVALepResBoundingUp"), Form("histo_ggWW_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_ggWW_CMS_hww_MVALepResBoundingDown"), Form("histo_ggWW_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_VV_CMS_hww_MVALepResBoundingUp"), Form("histo_VV_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_VV_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_VV_CMS_hww_MVALepResBoundingDown"), Form("histo_VV_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_Top_CMS_hww_MVALepResBoundingUp"), Form("histo_Top_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_Top_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_Top_CMS_hww_MVALepResBoundingDown"), Form("histo_Top_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_Wgamma_CMS_hww_MVALepResBoundingUp"), Form("histo_Wgamma_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_Wgamma_CMS_hww_MVALepResBoundingDown"), Form("histo_Wgamma_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVALepResBoundingUp = new TH1D( Form("histo_Ztt_CMS_hww_MVALepResBoundingUp"), Form("histo_Ztt_CMS_hww_MVALepResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVALepResBoundingUp->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVALepResBoundingDown = new TH1D( Form("histo_Ztt_CMS_hww_MVALepResBoundingDown"), Form("histo_Ztt_CMS_hww_MVALepResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVALepResBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_ttH_CMS_hww_MVAMETResBoundingUp"), Form("histo_ttH_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_ttH_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_ttH_CMS_hww_MVAMETResBoundingDown"), Form("histo_ttH_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_ZH_CMS_hww_MVAMETResBoundingUp"), Form("histo_ZH_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_ZH_CMS_hww_MVAMETResBoundingDown"), Form("histo_ZH_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_WH_CMS_hww_MVAMETResBoundingUp"), Form("histo_WH_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_WH_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_WH_CMS_hww_MVAMETResBoundingDown"), Form("histo_WH_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_qqH_CMS_hww_MVAMETResBoundingUp"), Form("histo_qqH_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_qqH_CMS_hww_MVAMETResBoundingDown"), Form("histo_qqH_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_ggH_CMS_hww_MVAMETResBoundingUp"), Form("histo_ggH_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_ggH_CMS_hww_MVAMETResBoundingDown"), Form("histo_ggH_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_qqWW_CMS_hww_MVAMETResBoundingUp"), Form("histo_qqWW_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_qqWW_CMS_hww_MVAMETResBoundingDown"), Form("histo_qqWW_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_ggWW_CMS_hww_MVAMETResBoundingUp"), Form("histo_ggWW_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_ggWW_CMS_hww_MVAMETResBoundingDown"), Form("histo_ggWW_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_VV_CMS_hww_MVAMETResBoundingUp"), Form("histo_VV_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_VV_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_VV_CMS_hww_MVAMETResBoundingDown"), Form("histo_VV_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_Top_CMS_hww_MVAMETResBoundingUp"), Form("histo_Top_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_Top_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_Top_CMS_hww_MVAMETResBoundingDown"), Form("histo_Top_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_Wgamma_CMS_hww_MVAMETResBoundingUp"), Form("histo_Wgamma_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_Wgamma_CMS_hww_MVAMETResBoundingDown"), Form("histo_Wgamma_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVAMETResBoundingUp = new TH1D( Form("histo_Ztt_CMS_hww_MVAMETResBoundingUp"), Form("histo_Ztt_CMS_hww_MVAMETResBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVAMETResBoundingUp->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVAMETResBoundingDown = new TH1D( Form("histo_Ztt_CMS_hww_MVAMETResBoundingDown"), Form("histo_Ztt_CMS_hww_MVAMETResBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVAMETResBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_ttH_CMS_hww_MVAJESBoundingUp"), Form("histo_ttH_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_ttH_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_ttH_CMS_hww_MVAJESBoundingDown"), Form("histo_ttH_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ttH_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_ZH_CMS_hww_MVAJESBoundingUp"), Form("histo_ZH_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_ZH_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_ZH_CMS_hww_MVAJESBoundingDown"), Form("histo_ZH_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ZH_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_WH_CMS_hww_MVAJESBoundingUp"), Form("histo_WH_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_WH_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_WH_CMS_hww_MVAJESBoundingDown"), Form("histo_WH_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_WH_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_qqH_CMS_hww_MVAJESBoundingUp"), Form("histo_qqH_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_qqH_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_qqH_CMS_hww_MVAJESBoundingDown"), Form("histo_qqH_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqH_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_ggH_CMS_hww_MVAJESBoundingUp"), Form("histo_ggH_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_ggH_CMS_hww_MVAJESBoundingDown"), Form("histo_ggH_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_qqWW_CMS_hww_MVAJESBoundingUp"), Form("histo_qqWW_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_qqWW_CMS_hww_MVAJESBoundingDown"), Form("histo_qqWW_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_qqWW_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_ggWW_CMS_hww_MVAJESBoundingUp"), Form("histo_ggWW_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_ggWW_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_ggWW_CMS_hww_MVAJESBoundingDown"), Form("histo_ggWW_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggWW_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_VV_CMS_hww_MVAJESBoundingUp"), Form("histo_VV_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_VV_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_VV_CMS_hww_MVAJESBoundingDown"), Form("histo_VV_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_VV_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_Top_CMS_hww_MVAJESBoundingUp"), Form("histo_Top_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_Top_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_Top_CMS_hww_MVAJESBoundingDown"), Form("histo_Top_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Top_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_Wgamma_CMS_hww_MVAJESBoundingUp"), Form("histo_Wgamma_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_Wgamma_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_Wgamma_CMS_hww_MVAJESBoundingDown"), Form("histo_Wgamma_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Wgamma_CMS_hww_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVAJESBoundingUp = new TH1D( Form("histo_Ztt_CMS_hww_MVAJESBoundingUp"), Form("histo_Ztt_CMS_hww_MVAJESBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVAJESBoundingUp->Sumw2();
  TH1D* histo_Ztt_CMS_hww_MVAJESBoundingDown = new TH1D( Form("histo_Ztt_CMS_hww_MVAJESBoundingDown"), Form("histo_Ztt_CMS_hww_MVAJESBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_Ztt_CMS_hww_MVAJESBoundingDown->Sumw2();

  TH1D* histo_ggH_CMS_hww_MVAggHBoundingUp = new TH1D( Form("histo_ggH_CMS_hww_MVAggHBoundingUp"), Form("histo_ggH_CMS_hww_MVAggHBoundingUp"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVAggHBoundingUp->Sumw2();
  TH1D* histo_ggH_CMS_hww_MVAggHBoundingDown = new TH1D( Form("histo_ggH_CMS_hww_MVAggHBoundingDown"), Form("histo_ggH_CMS_hww_MVAggHBoundingDown"), nBinHis, minHis[1], maxHis[1]); histo_ggH_CMS_hww_MVAggHBoundingDown->Sumw2();

  TH1D* histoS = new TH1D("histoS", "histoS", 180, 0, 180);
  TH1D* histoB = new TH1D("histoB", "histoB", 180, 0, 180);
  TH1D* histoD = new TH1D("histoD", "histoD", 180, 0, 180);
  histoS->Sumw2();
  histoB->Sumw2();
  histoD->Sumw2();

  //----------------------------------------------------------------------------
  // Define MVA cut values (known a posteriori)
  //----------------------------------------------------------------------------
  int useVar = 0; // which MVA to be used
  double useCut = -999.0;
  if    (nJetsType == 0){
    if     (mH == 110){
      //useVar = 2;
      //useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.405-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 115 || mH == 118 ||  mH == 122 || mH == 124 || mH == 125 || mH == 126 || mH == 128 || mH == 135 || mH == 145){
      //useVar = 2;
      //useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
      useVar = 4;
      useCut = 0.405-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 120 || mH == 155 || mH == 700 || mH == 800 || mH == 900 || mH == 1000){
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
    if     (mH == 110){
      useVar = 4;
      useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 115 || mH == 118 ||  mH == 122 || mH == 124 || mH == 125 || mH == 126 || mH == 128 || mH == 135 || mH == 145){
      useVar = 4;
      useCut = 0.685-datMVA[useVar]->GetBinWidth(0)/2.0;
    }	
    else if(mH == 120 || mH == 155 || mH == 700 || mH == 800 || mH == 900 || mH == 1000){
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
  if(nJetsType !=2 && useVar != 4) assert(0);
  //----------------------------------------------------------------------------

  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          npu;
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
  Float_t         jet3Btag;
  Int_t 	  lep1McId;
  Int_t 	  lep2McId;
  Int_t 	  lep1MotherMcId;
  Int_t 	  lep2MotherMcId;
  Float_t         bdt = 0.0;
  Float_t         bdtd = 0.0;
  Float_t         nn = 0.0;
  Float_t         knn = 0.0;
  Float_t         bdtg = 0.0;
  Float_t         bdtg_aux0 = 0.0;
  Float_t         bdtg_aux1 = 0.0;
  Float_t         bdtg_aux2 = 0.0;
  Float_t         higgsPt = -999;
  Float_t         bdtg_wjets = 0.0;
  //Float_t         knn_wjets = 0.0;
  Float_t         sfWeightPU;
  Float_t         sfWeightEff;
  Float_t         sfWeightTrig;
  Float_t         sfWeightHPt;
  Float_t         dymva = -100.0;

  //****************************************************************************
  //
  // Loop Over Signal Sample
  //
  //****************************************************************************

  signal->SetBranchAddress( "cuts"         , &cuts         );
  signal->SetBranchAddress( "dstype"       , &dstype       );
  signal->SetBranchAddress( "nvtx"         , &nvtx         );
  signal->SetBranchAddress( "npu"          , &npu          );
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
  signal->SetBranchAddress( "jet3Btag"     , &jet3Btag     );
  signal->SetBranchAddress(Form("bdt_hww%i_%djet_ww"	  ,mH,nJetsType), &bdt        );
  signal->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"	  ,mH,nJetsType), &bdtd       );
  signal->SetBranchAddress(Form("nn_hww%i_%djet_ww"	  ,mH,nJetsType), &nn	      );
  signal->SetBranchAddress(Form("knn_hww%i_%djet_ww"	  ,mH,nJetsType), &knn        );
  signal->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"	  ,mH,nJetsType), &bdtg       );
  signal->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux0",mH,nJetsType), &bdtg_aux0  );
  signal->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux1",mH,nJetsType), &bdtg_aux1  );
  signal->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux2",mH,nJetsType), &bdtg_aux2  );
  signal->SetBranchAddress( "higgsPt"      , &higgsPt      );
  //signal->SetBranchAddress(Form("bdtg_hww%i_%djet_wjets"  ,mH,nJetsType), &bdtg_wjets );
  //signal->SetBranchAddress(Form("knn_hww%i_%djet_wjets"   ,mH,nJetsType), &knn_wjets  );
  signal->SetBranchAddress( "sfWeightPU"     , &sfWeightPU     );
  signal->SetBranchAddress( "sfWeightEff"    , &sfWeightEff    );
  signal->SetBranchAddress( "sfWeightTrig"   , &sfWeightTrig   );
  signal->SetBranchAddress( "sfWeightHPt"    , &sfWeightHPt    );
  signal->SetBranchAddress( "dymva"          , &dymva          );

  float nSigAcc[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigCut[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigMVA[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigEAcc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigECut[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigEMVA[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)signal->GetEntries());

    bool lId = (cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
    if(category == 1) lId = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2LooseEleV1)    == SmurfTree::Lep2LooseEleV1   ) ||
    			    ((cuts & SmurfTree::Lep1LooseEleV1)    == SmurfTree::Lep1LooseEleV1    && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection);
    if(!lId) continue;

    //----------------------------------------------------------------------------
    //for data require that the event fired one of the designated signal triggers
    //----------------------------------------------------------------------------
    if(dstype == SmurfTree::data &&
      (cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;

    //----------------------------------------------------------------------------
    //define jet bin. 
    //For 2-Jet Bin:
    //  - apply "central jet veto"
    //    ie. no third jet above 30 GeV that lies between the 
    //    first two jets in pseudorapidity
    //----------------------------------------------------------------------------
    unsigned int Njet3 = njets;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(jet3->pt() <= 30)					           Njet3 = 2;
      else if(jet3->pt() > 30 && (
    	(jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    	(jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 0;
      else							           Njet3 = 2;
      if(njets < 2 || njets > 3)                                           Njet3 = 0;
    }

    //----------------------------------------------------------------------------
    //For Jet Energy scale systematics
    //----------------------------------------------------------------------------
    bool passJetCut[3] = {Njet3 == nJetsType, false, false};
    if(nJetsType == 0 && 			 jet1->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 4.7  			       ) passJetCut[1] = true;
    if(nJetsType == 0 && 			 jet1->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 4.7  			       ) passJetCut[2] = true;
    if(nJetsType == 1 && jet1->pt()*1.05 > 30 && jet2->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 4.7 && TMath::Abs(jet2->eta()) < 4.7) passJetCut[1] = true;
    if(nJetsType == 1 && jet1->pt()*0.95 > 30 && jet2->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 4.7 && TMath::Abs(jet2->eta()) < 4.7) passJetCut[2] = true;

    double minmet = TMath::Min(pmet,pTrackMet);
    bool passMET = minmet > 20.;
    if(useDYMVA == kFALSE){
      if     (njets == 0) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0   || type == SmurfTree::em || type == SmurfTree::me);
    } else {
      if     (njets == 0) passMET = passMET && (dymva >  0.60 || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (dymva >  0.30 || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0    || type == SmurfTree::em || type == SmurfTree::me);
    }

    bool passNewCuts = true;
    if(dilep->pt() <= 45) passNewCuts = false;

    //----------------------------------------------------------------------------
    //To create the DY bkg systematics MVA shapes (makeZjetsTemplates == true)
    //we loosen the MET cut to "minmet > 20 GeV". 
    //----------------------------------------------------------------------------
    if     ( makeZjetsTemplates == true && passMET == true ) passMET = false;
    else if( makeZjetsTemplates == true && passMET == false) passMET = minmet > 20.;

    //----------------------------------------------------------------------------
    // WW Preselection
    //----------------------------------------------------------------------------
    if( passJetCut[0]==false&&passJetCut[1]==false&&passJetCut[2]==false ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= ptlMin	    					 ) continue; // cut on trailing lepton pt
    if( passNewCuts == false                                             ) continue; // cut on new pt cuts
    if( passMET == false                                                 ) continue; // cut on pmet
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    //if( lid3 != 0	                                                 ) continue; // cut on dileptons
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJetCut = kTRUE;
    if(useDYMVA == kFALSE){
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    }
    if(njets >= 2) dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

    //----------------------------------------------------------------------------
    // Define Final States:
    // 5: Same Flavor
    // 6: Opposite Flavor
    //----------------------------------------------------------------------------
    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);
    if(wwDecayCut == false) continue;

    if(is2DAna == true) {
      double varA = (dilep->mass()-12.0)/(dilmass_cut-12.0);
      double varB = (mt-mtLowerCut)/(mtUpperCut-mtLowerCut);
      bdtg = Unroll2VarTo1Var(varA,varB,binVarA,binVarB,false);
    }
    if(bdtg < -1.0) bdtg = -0.999;
    if(bdtg > +1.0) bdtg =  0.999;

    //----------------------------------------------------------------------------
    // Define event weights    
    // Apply lepton efficiency scale factors, trigger efficiencies
    //----------------------------------------------------------------------------
    double add = 1.0;
    add = add*nPUScaleFactor2012(fhDPU2012,npu);

    double addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
                       leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
    add = add*addLepEff;

    add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
        					      fabs(lep2->eta()), lep2->pt(), 
        					      TMath::Abs( lid1), TMath::Abs(lid2));
    double addggH = 1.0;
    double addggHSyst[8] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
    // add = 1.0;
    if(isFermioPhobic == true) {
      add = add * enhancementFactor(mH,2);
      if(processId == 121 || processId == 122) add = 0.0;
    }

    if(processId != 10010 && isSM4 == true) add = 0.0;

    if(isSM4 == true) add = add * enhancementFactor(mH,1); // BR(H->WW) enhancement factor

    if (processId == 10010) {
      //float theMax = 0.00;
      //addggH = HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(TMath::Max(higgsPt, theMax)));
      //for(int ns=0; ns<8; ns++) if(HiggsPtKFactorSyst[ns]) addggHSyst[ns] = HiggsPtKFactorSyst[ns]->GetBinContent( HiggsPtKFactorSyst[ns]->GetXaxis()->FindFixBin(TMath::Max(higgsPt, theMax)));

      if(isFermioPhobic == true) addggH = 0.0;
      if(isSM4 == true) {
        addggH = addggH * enhancementFactor(mH,0); // ggH enhancement factor
	for(int ns=0; ns<8; ns++) addggHSyst[ns] = addggHSyst[ns] * enhancementFactor(mH,0); // ggH enhancement factor
      }
    }

    add = add*addggH;
    double myWeight = scaleFactorLum * scale1fb * add;

    if(myWeight == 0) continue;

    //----------------------------------------------------------------------------
    // Classify Signal Events by production mechanism
    //----------------------------------------------------------------------------
    int nSigBin = -1;
    // GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
    if     (processId==121 ||
    	    processId==122)   nSigBin = 1;
    else if(processId==24)    nSigBin = 2;
    else if(processId==26)    nSigBin = 3;
    else if(processId==10001) nSigBin = 4;
    else if(processId==10010) nSigBin = 5;
    else  {return;}

    //----------------------------------------------------------------------------
    //
    // Higgs Signal Selection Cuts
    //
    //----------------------------------------------------------------------------
    double theCutPtMinLow = cutPtMinLow (mH, type);
    bool passAllCuts = dilep->mass()         < theCutMassHigh &&
        	       mt		     > theCutMTLow &&
        	       mt		     < theCutMTHigh &&
        	       lep1->pt()	     > theCutPtMaxLow &&
        	       lep2->pt()	     > theCutPtMinLow &&
        	       dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh &&
		       passJetCut[0]==true;

    //----------------------------------------------------------------------------
    // VBF selection cuts for 2-Jet bin
    //----------------------------------------------------------------------------
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    mt > 30 && mt < mtUpperCut &&
		    centrality == 1 &&
		    passJetCut[0]==true;
      if(is2DAna == true) {
        passAllCuts = passAllCuts &&
                      dilep->mass()         < theCutMassHigh &&
        	      mt 		    < theCutMTHigh &&
        	      lep1->pt() 	    > theCutPtMaxLow &&
        	      lep2->pt() 	    > theCutPtMinLow &&
        	      dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh;
      }
    }

    //----------------------------------------------------------------------------
    // Add signal yields for Cut-Based analysis
    //----------------------------------------------------------------------------
    if(passAllCuts == true) {
      nSigCut[0]  = nSigCut[0]  + myWeight;
      nSigECut[0] = nSigECut[0] + myWeight*myWeight;
      nSigCut[nSigBin]  = nSigCut[nSigBin]  + myWeight;
      nSigECut[nSigBin] = nSigECut[nSigBin] + myWeight*myWeight;
      sigMVA[5][0]      ->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), myWeight);
      sigMVA[5][nSigBin]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), myWeight);
    }

    //----------------------------------------------------------------------------
    // Add signal yields and fill MVA output distributions for MVA Shape analysis
    // Apply mT_Higgs cut for MVA Shape analysis
    //----------------------------------------------------------------------------
    bool passMVAPreselCuts = mt > mtLowerCut && mt < mtUpperCut; if(wwPresel == true) passMVAPreselCuts = true;
    if(is2DAna == true) passMVAPreselCuts = passMVAPreselCuts && lep1->pt() > theCutPtMaxLow;
    if(passMVAPreselCuts == true && passJetCut[0] == true){
      nSigAcc[0]  = nSigAcc[0]  + myWeight;
      nSigEAcc[0] = nSigEAcc[0] + myWeight*myWeight;
      nSigAcc[nSigBin]  = nSigAcc[nSigBin]  + myWeight;
      nSigEAcc[nSigBin] = nSigEAcc[nSigBin] + myWeight*myWeight;

      if(bdtg>-1.0&&bdtg<1){
        double myVar = dPhi*180.0/TMath::Pi(); myVar = CalcGammaMRstar(*lep1,*lep2);
        histoS->Fill(myVar,myWeight);
      }

      if     (useVar == 0) histos   ->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),       myWeight);
      else if(useVar == 1) histos   ->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),       myWeight);
      else if(useVar == 2) histos   ->Fill(TMath::Max(TMath::Min((double)nn  ,maxHis[2]-0.001),minHis[2]+0.001),       myWeight);
      else if(useVar == 3) histos   ->Fill(TMath::Max(TMath::Min((double)knn ,maxHis[3]-0.001),minHis[3]+0.001),       myWeight);
      else if(useVar == 4) histos   ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),       myWeight);
      else if(useVar == 5) histos   ->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001), myWeight);

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
      
      // WARNING, THIS IS ONLY GOOD FOR BDTG!
      if(useExpTemplates == true){
        //----------------------------------------------------------------------------
        // MVA Shape systematics for:
        // 1: Lepton efficiencies
        // 2: Lepton Resolution (bdtg_aux0,bdtg_aux1 branch)
        // 2: Met Resolution    (bdtg_aux2 branch)
        //----------------------------------------------------------------------------
        double addLepEffUp = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1,1)*
                             leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2,1);
        double addLepEffDown = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1,-1)*
                               leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2,-1);
        if     (nSigBin == 1){
	  histo_ttH_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_ttH_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_ttH_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ttH_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ttH_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(nSigBin == 2){
	  histo_ZH_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_ZH_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_ZH_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ZH_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ZH_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(nSigBin == 3){
	  histo_WH_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_WH_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_WH_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_WH_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_WH_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(nSigBin == 4){
	  histo_qqH_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_qqH_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_qqH_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_qqH_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_qqH_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(nSigBin == 5){
	  histo_ggH_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_ggH_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_ggH_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ggH_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ggH_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
      }
      // WARNING, THIS IS ONLY GOOD FOR BDTG!
      //----------------------------------------------------------------------------
      // MVA Shape systematics for:
      // Gluon Fusion Higgs Missing Higher Order Corrections (weights [0] & [5])
      //----------------------------------------------------------------------------
      if(useggHTemplates == true && nSigBin == 5 && addggH != 0){
	histo_ggH_CMS_hww_MVAggHBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addggHSyst[0]/addggH);
	histo_ggH_CMS_hww_MVAggHBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addggHSyst[5]/addggH);
      }

      //----------------------------------------------------------------------------
      // Add signal yields for MVA Cut analysis
      //----------------------------------------------------------------------------
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
    } // passMVAPreselCuts

    //----------------------------------------------------------------------------
    // MVA Shape systematics for:
    // Jet Energy Scale
    //----------------------------------------------------------------------------
    if(passMVAPreselCuts == true && useJESTemplates == true && (passJetCut[1] == true || passJetCut[2] == true)){ // Begin JES
      if     (nSigBin == 1){
        if(passJetCut[1] == true) histo_ttH_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_ttH_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(nSigBin == 2){
        if(passJetCut[1] == true) histo_ZH_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_ZH_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(nSigBin == 3){
        if(passJetCut[1] == true) histo_WH_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_WH_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(nSigBin == 4){
        if(passJetCut[1] == true) histo_qqH_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_qqH_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(nSigBin == 5){
        if(passJetCut[1] == true) histo_ggH_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_ggH_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
    } // End JES
  } // end loop over signal events

  for(int i=0; i<6; i++) nSigEAcc[i] = sqrt(nSigEAcc[i]);
  for(int i=0; i<6; i++) nSigECut[i] = sqrt(nSigECut[i]);
  for(int i=0; i<6; i++) nSigEMVA[i] = sqrt(nSigEMVA[i]);
  printf("--- Finished Signal loop\n");

  //****************************************************************************
  //
  // Loop Over Background Sample
  //
  //****************************************************************************
  background->SetBranchAddress( "cuts"          , &cuts 	  );
  background->SetBranchAddress( "dstype"        , &dstype	  );
  background->SetBranchAddress( "nvtx"          , &nvtx 	  );
  background->SetBranchAddress( "npu"           , &npu 	          );
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
  background->SetBranchAddress( "jet3Btag"      , &jet3Btag	  );
  background->SetBranchAddress( "lep1McId"      , &lep1McId	  );
  background->SetBranchAddress( "lep2McId"      , &lep2McId	  );
  background->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId );
  background->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId );
  background->SetBranchAddress(Form("bdt_hww%i_%djet_ww"      ,mH,nJetsType), &bdt	  );
  background->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"     ,mH,nJetsType), &bdtd	  );
  background->SetBranchAddress(Form("nn_hww%i_%djet_ww"       ,mH,nJetsType), &nn	  );
  background->SetBranchAddress(Form("knn_hww%i_%djet_ww"      ,mH,nJetsType), &knn	  );
  background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"     ,mH,nJetsType), &bdtg	  );
  background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux0",mH,nJetsType), &bdtg_aux0  );
  background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux1",mH,nJetsType), &bdtg_aux1  );
  background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux2",mH,nJetsType), &bdtg_aux2  );
  //background->SetBranchAddress(Form("bdtg_hww%i_%djet_wjets"  ,mH,nJetsType), &bdtg_wjets );
  //background->SetBranchAddress(Form("knn_hww%i_%djet_wjets"   ,mH,nJetsType), &knn_wjets  );
  background->SetBranchAddress( "dymva"          , &dymva          );

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
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    //----------------------------------------------------------------------------
    //for data require that the event fired one of the designated signal triggers
    //----------------------------------------------------------------------------
    if(dstype == SmurfTree::data &&
      (cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dstype == SmurfTree::data && run <  minRun) continue;
    if(dstype == SmurfTree::data && run >  maxRun) continue;

    //----------------------------------------------------------------------------
    //define jet bin. 
    //For 2-Jet Bin:
    //  - apply "central jet veto"
    //    ie. no third jet above 30 GeV that lies between the 
    //    first two jets in pseudorapidity
    //----------------------------------------------------------------------------
    unsigned int Njet3 = njets;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(jet3->pt() <= 30)					           Njet3 = 2;
      else if(jet3->pt() > 30 && (
    	(jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    	(jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 0;
      else							           Njet3 = 2;
      if(njets < 2 || njets > 3)                                           Njet3 = 0;
    }

    //----------------------------------------------------------------------------
    //For Jet Energy scale systematics
    //----------------------------------------------------------------------------
    bool passJetCut[3] = {Njet3 == nJetsType, false, false};
    if(nJetsType == 0 && 			 jet1->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 4.7  			       ) passJetCut[1] = true;
    if(nJetsType == 0 && 			 jet1->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 4.7  			       ) passJetCut[2] = true;
    if(nJetsType == 1 && jet1->pt()*1.05 > 30 && jet2->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 4.7 && TMath::Abs(jet2->eta()) < 4.7) passJetCut[1] = true;
    if(nJetsType == 1 && jet1->pt()*0.95 > 30 && jet2->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 4.7 && TMath::Abs(jet2->eta()) < 4.7) passJetCut[2] = true;

    double minmet = TMath::Min(pmet,pTrackMet);
    bool passMET = minmet > 20.;
    if(useDYMVA == kFALSE){
      if     (njets == 0) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0   || type == SmurfTree::em || type == SmurfTree::me);
    } else {
      if     (njets == 0) passMET = passMET && (dymva >  0.60 || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (dymva >  0.30 || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0    || type == SmurfTree::em || type == SmurfTree::me);
    }

    bool passNewCuts = true;
    if(dilep->pt() <= 45) passNewCuts = false;

    //----------------------------------------------------------------------------
    //To create the DY bkg systematics MVA shapes (makeZjetsTemplates == true)
    //we loosen the MET cut to "minmet > 20 GeV". 
    //----------------------------------------------------------------------------
    if     ( makeZjetsTemplates == true && passMET == true ) passMET = false;
    else if( makeZjetsTemplates == true && passMET == false) passMET = minmet > 20.;

    //----------------------------------------------------------------------------
    // WW Preselection
    //----------------------------------------------------------------------------
    bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
                        ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
                         dstype != SmurfTree::data;

    if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
    if( passJetCut[0]==false&&passJetCut[1]==false&&passJetCut[2]==false ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= ptlMin	    					 ) continue; // cut on trailing lepton pt
    if( passNewCuts == false                                             ) continue; // cut on new pt cuts
    if( passMET == false                                                 ) continue; // cut on pmet
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    //if( lid3 != 0	                                                 ) continue; // cut on dileptons
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJetCut = kTRUE;
    if(useDYMVA == kFALSE){
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    }
    if(njets >= 2) dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

    //----------------------------------------------------------------------------
    // Define Background Type
    // 0 : qqWW
    // 1 : ggWW
    // 2 : WZ,ZZ
    // 3 : Top (ttbar, tW)
    // 4 : DY
    // 5 : W+Jets (fakes)
    // 6 : W+gamma , W+gamma*
    // 7 : DY->tautau
    //----------------------------------------------------------------------------
    int fDecay = 0;
    if     (dstype == SmurfTree::wjets  	 ) fDecay = 5;
    else if(dstype == SmurfTree::ttbar  	 ) fDecay = 3;
    else if(dstype == SmurfTree::dyee   	 ) fDecay = 4;
    else if(dstype == SmurfTree::dymm   	 ) fDecay = 4;
    else if(dstype == SmurfTree::dytt   	 ) fDecay = 7;
    else if(dstype == SmurfTree::tw     	 ) fDecay = 3;
    else if(dstype == SmurfTree::qqww   	 ) fDecay = 0;
    else if(dstype == SmurfTree::wz     	 ) fDecay = 2;
    else if(dstype == SmurfTree::zz     	 ) fDecay = 2;
    else if(dstype == SmurfTree::www     	 ) fDecay = 2;
    else if(dstype == SmurfTree::ggww   	 ) fDecay = 1;
    else if(dstype == SmurfTree::wgamma 	 ) fDecay = 6;
    else if(dstype == SmurfTree::wgstar 	 ) fDecay = 6;
    else if(dstype == SmurfTree::data   	 ) fDecay = 5;
    else if(dstype == SmurfTree::dyttDataDriven  ) fDecay = 7;
    else if(dstype == SmurfTree::qcd             ) fDecay = 7;
    else                                 {printf("bad dstype: %d\n",dstype); assert(0);}
    if(dstype == SmurfTree::wz || dstype == SmurfTree::zz) {
      if(lep1MotherMcId == 23 && lep2MotherMcId == 23) {
        fDecay = 4;
      }
    }

    //----------------------------------------------------------------------------
    // Define Final States:
    // 5: Same Flavor
    // 6: Opposite Flavor
    //----------------------------------------------------------------------------
    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);
    if(wwDecayCut == false) continue;

    if(is2DAna == true) {
      double varA = (dilep->mass()-12.0)/(dilmass_cut-12.0);
      double varB = (mt-mtLowerCut)/(mtUpperCut-mtLowerCut);
      bdtg = Unroll2VarTo1Var(varA,varB,binVarA,binVarB,false);
    }
    if(bdtg < -1.0) bdtg = -0.999;
    if(bdtg > +1.0) bdtg =  0.999;

    //----------------------------------------------------------------------------
    // Define event weights    
    //----------------------------------------------------------------------------
    double myWeight = 1.0;
    double add      = 1.0;

    //----------------------------------------------------------------------------
    // First classify event into tight+tight, tight+fail, fail+fail 
    //----------------------------------------------------------------------------
    int nFake = 0;
    if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
    if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
    if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    if(category == 1) {
      if(((cuts & SmurfTree::Lep1LooseEleV1)  == SmurfTree::Lep1LooseEleV1)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake--;
      if(((cuts & SmurfTree::Lep2LooseEleV1)  == SmurfTree::Lep2LooseEleV1)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake--;
    }
    if(nFake < 0) assert(0);

    bool isRealLepton = false;
    if((TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) &&
       (TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13)) isRealLepton = true;
    double addLepEff = 1.0;
    double addFR     = 1.0;

    //----------------------------------------------------------------------------
    // Explicitly neglect fail+fail events
    //----------------------------------------------------------------------------
    if(nFake > 1){
      myWeight = 0.0;
    }
    //----------------------------------------------------------------------------
    // Tight+Fail events
    //----------------------------------------------------------------------------
    else if(nFake == 1){
      //----------------------------------------------------------------------------
      // For data, apply fake rate to generate bkg prediction
      //----------------------------------------------------------------------------
      if(dstype == SmurfTree::data){
        addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        							          (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        							          (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        add = addFR;
	if(category == 1) add = add*1.10; // HACK!!!
	fDecay  	      = 5;
        myWeight	      = add;
      }
      //----------------------------------------------------------------------------
      // For real lepton or w+gamma:
      // apply fake rates, and give negative weight to subtract non jet-fake
      // contamination in the tight+fail sample
      //----------------------------------------------------------------------------
      else if(isRealLepton == true || dstype == SmurfTree::wgamma){
    	addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    		        						  (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    									  (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

        add = addFR;
    	add = add*nPUScaleFactor2012(fhDPU2012,npu);

        addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
                    leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
        add = add*addLepEff;

    	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
	if(category == 1) add = add*1.10; // HACK!!!
        fDecay  	       = 5;
        myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
      }
      //----------------------------------------------------------------------------
      // Neglect fake lepton events in MC
      //----------------------------------------------------------------------------
      else {
        myWeight = 0.0;
      }
    }

    //----------------------------------------------------------------------------
    // Neglect any tight+tight events from data
    //----------------------------------------------------------------------------
    else if(dstype == SmurfTree::data) myWeight = 0.0;

    //----------------------------------------------------------------------------
    // DY->TauTau Embedded sample events: apply normalization weight.
    //----------------------------------------------------------------------------
    else if(dstype== SmurfTree::dyttDataDriven || dstype == SmurfTree::qcd) {
      myWeight = ZttScaleFactor(nvtx,period,scale1fb)*scaleFactorLum;
    }

    //----------------------------------------------------------------------------
    // Regular Tight+Tight events from Monte Carlo
    //----------------------------------------------------------------------------
    else if(dstype != SmurfTree::data){
      //----------------------------------------------------------------------------      
      // Apply lepton efficiency scale factors, trigger efficiencies,
      // Pileup weights
      //----------------------------------------------------------------------------
      add = 1.0;
      add = add*nPUScaleFactor2012(fhDPU2012,npu);

      addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
        	  leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
      add = add*addLepEff;
      add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							fabs(lep2->eta()), lep2->pt(), 
    							TMath::Abs( lid1), TMath::Abs(lid2));

      //----------------------------------------------------------------------------      
      // Apply DY Bkg Scale Factors
      //----------------------------------------------------------------------------
      if(fDecay == 4  && (type   == SmurfTree::mm   || type   == SmurfTree::ee)
                      && (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm)) {
    	if(njets == 0) add=add*DYBkgScaleFactor(0,0); 
    	if(njets == 1) add=add*DYBkgScaleFactor(0,1); 
    	if(njets >= 2) add=add*DYBkgScaleFactor(0,2);
      }

      //----------------------------------------------------------------------------      
      // Apply Top Bkg Scale Factors
      //----------------------------------------------------------------------------
      if(fDecay == 3) {
    	if(njets == 0) add=add*TopBkgScaleFactor(0);
    	if(njets == 1) add=add*TopBkgScaleFactor(1); 
    	if(njets >= 2) add=add*TopBkgScaleFactor(2); 
      }

      //----------------------------------------------------------------------------      
      // Apply W+Jets Bkg Scale Factor for MC (not nominally used)
      //----------------------------------------------------------------------------
      if(fDecay == 5) add=add*WJetsMCScaleFactor(); 

      //----------------------------------------------------------------------------      
      // Apply W+gamma* normalization scale factor
      //----------------------------------------------------------------------------
      if(dstype == SmurfTree::wgstar) add=add*WGstarScaleFactor();

      //----------------------------------------------------------------------------      
      // Apply WW Bkg Scale Factors 
      // Don't do this for the WW selection (wwPresel == true)
      //----------------------------------------------------------------------------
      if((fDecay == 0 || fDecay == 1) && wwPresel == false){     
        if(njets == 0) add=add*WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0); 
        else	       add=add*WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1); 
      }
      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;
      myWeight = scale1fb*scaleFactorLum*add;
    }

    if(category == 1){
      bool passCuts = (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection ||
                      (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection;
      if(passCuts == false) myWeight = 0;
    }
 
    //----------------------------------------------------------------------------      
    // Explicitly neglect events with 0 weight
    //----------------------------------------------------------------------------

    if(myWeight == 0) continue;

    //----------------------------------------------------------------------------
    //
    // Higgs Signal Selection Cuts
    //
    //----------------------------------------------------------------------------
    double theCutPtMinLow = cutPtMinLow (mH, type);
    bool passAllCuts = dilep->mass()         < theCutMassHigh &&
        	       mt		     > theCutMTLow &&
        	       mt		     < theCutMTHigh &&
        	       lep1->pt()	     > theCutPtMaxLow &&
        	       lep2->pt()	     > theCutPtMinLow &&
        	       dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh &&
		       passJetCut[0]==true;

    //----------------------------------------------------------------------------
    // VBF selection cuts for 2-Jet bin
    //----------------------------------------------------------------------------
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    mt > 30 && mt < mtUpperCut &&
		    centrality == 1 &&
		    passJetCut[0]==true;
      if(is2DAna == true) {
        passAllCuts = passAllCuts &&
                      dilep->mass()         < theCutMassHigh &&
        	      mt 		    < theCutMTHigh &&
        	      lep1->pt() 	    > theCutPtMaxLow &&
        	      lep2->pt() 	    > theCutPtMinLow &&
        	      dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh;
      }
    }

    //----------------------------------------------------------------------------
    // Add bkg yields for Cut-Based analysis
    //----------------------------------------------------------------------------
    if(passAllCuts == true) {
      double newWeight = myWeight;
      //----------------------------------------------------------------------------
      // For Cut-Based Analysis, use Cut-based analysis specfic 
      // scale factors for WW Bkg. We originally applied WW scale factors for MVA
      // selection, and correct it by the ratio here
      //----------------------------------------------------------------------------      
      if((fDecay == 0 || fDecay == 1) && wwPresel == false){ // only for WW
	if(njets == 0) newWeight=newWeight*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),0)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0);
	else           newWeight=newWeight*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),1)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1);	   
      }

      //----------------------------------------------------------------------------
      // For Cut-Based Analysis, use mass dependant  
      // scale factors for DY Bkg. We originally applied DY scale factors for the
      // WW pre-selection, and correct it by the ratio here.
      //----------------------------------------------------------------------------      
      if(fDecay == 4&& (dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
                       (type   == SmurfTree::mm   || type   == SmurfTree::ee)){
	if(((nJetsType != 2 && mH <= 300) || nJetsType == 2) && (wwDecay == 4 || wwDecay == 5)){
          //newWeight=newWeight*DYBkgScaleFactor(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2))/DYBkgScaleFactor(0,TMath::Min((int)nJetsType,2));
          newWeight=0.0;
	}
	DYXS[1] += newWeight;
      }
      else if(fDecay == 4){
	VVXS[1] += newWeight;
      }

      if(fDecay == 3 && njets >= 2) newWeight=newWeight*TopVBFBkgScaleFactor(mH)/TopBkgScaleFactor(2);

      nBgdCut  = nBgdCut  + newWeight;
      nBgdECut = nBgdECut + newWeight*newWeight;
      nBgdCutDecays[fDecay]  = nBgdCutDecays[fDecay]  + newWeight;
      nBgdECutDecays[fDecay] = nBgdECutDecays[fDecay] + newWeight*newWeight;
      bgdMVA[5]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), newWeight);
      bgdMVADecays[5][fDecay]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), newWeight);
    }

    //----------------------------------------------------------------------------
    // Add bkg yields and fill MVA output distributions for MVA Shape analysis
    // Apply mT_Higgs cut for MVA Shape analysis
    //----------------------------------------------------------------------------
    bool passMVAPreselCuts = mt > mtLowerCut && mt < mtUpperCut; if(wwPresel == true) passMVAPreselCuts = true;
    if(is2DAna == true) passMVAPreselCuts = passMVAPreselCuts && lep1->pt() > theCutPtMaxLow;
    if(passMVAPreselCuts == true && passJetCut[0] == true){
      double myWeightMVA = myWeight;

      if(fDecay == 4 && (dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
                        (type   == SmurfTree::mm   || type   == SmurfTree::ee) && makeZjetsTemplates == false){
	if(nJetsType != 2 && mH <= 300 && (wwDecay == 4 || wwDecay == 5)){
          myWeightMVA = 0.0;
	}
      }
      //----------------------------------------------------------------------------
      // The systematics shape for the DY Bkg is derived from the MC with the
      // nominal MET cut.
      //----------------------------------------------------------------------------
      if(fDecay == 4 && (dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
	                (type   == SmurfTree::mm   || type   == SmurfTree::ee) && makeZjetsTemplates == false){
	DYXS[0] = DYXS[0] + myWeightMVA;
	myWeightMVA = 0.0;
      }
      else if(fDecay == 4 && makeZjetsTemplates == false){
	VVXS[0] = VVXS[0] + myWeightMVA;
	histoVV->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeightMVA);
      }
      else if(fDecay == 4 && (dstype == SmurfTree::wz || dstype == SmurfTree::zz)){
	histoVV->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeightMVA);
      }
      nBgdAcc  = nBgdAcc  + myWeightMVA;
      nBgdEAcc = nBgdEAcc + myWeightMVA*myWeightMVA;
      nBgdAccDecays[fDecay]  = nBgdAccDecays[fDecay]  + myWeightMVA;
      nBgdEAccDecays[fDecay] = nBgdEAccDecays[fDecay] + myWeightMVA*myWeightMVA;

      //----------------------------------------------------------------------------
      // For DY Bkg
      // When makeZjetsTemplates == true, we use only the DY->ee/mm samples
      // Otherwise, include the DY->tautau bkg into "histo1" as well
      //----------------------------------------------------------------------------
      if      (((fDecay == 4 || fDecay == 7) && makeZjetsTemplates == false) ||
               ((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) && fDecay == 4 &&
	        (type   == SmurfTree::mm   || type   == SmurfTree::ee) && makeZjetsTemplates == true)){ // Z+jets
	if     (useVar == 0) histo1->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),  myWeightMVA);
	else if(useVar == 1) histo1->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),  myWeightMVA);
	else if(useVar == 2) histo1->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),    myWeightMVA);
	else if(useVar == 3) histo1->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),  myWeightMVA);
	else if(useVar == 4) histo1->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
	else if(useVar == 5) histo1->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
      }
      //----------------------------------------------------------------------------
      // For Top Bkg
      //----------------------------------------------------------------------------
      else if(fDecay == 3){
	if     (useVar == 0) histo2->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),  myWeightMVA);
	else if(useVar == 1) histo2->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),  myWeightMVA);
	else if(useVar == 2) histo2->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),    myWeightMVA);
	else if(useVar == 3) histo2->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),  myWeightMVA);
	else if(useVar == 4) histo2->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
	else if(useVar == 5) histo2->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
      }
      //----------------------------------------------------------------------------
      // For WZ,ZZ Bkg
      //----------------------------------------------------------------------------
      else if(fDecay == 2){ 
	if     (useVar == 0) histo3->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),  myWeightMVA);
	else if(useVar == 1) histo3->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),  myWeightMVA);
	else if(useVar == 2) histo3->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),    myWeightMVA);
	else if(useVar == 3) histo3->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),  myWeightMVA);
	else if(useVar == 4) histo3->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
	else if(useVar == 5) histo3->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
      }
      //----------------------------------------------------------------------------
      // For W+Jets, W+gamma, W+gamma* Bkg
      //----------------------------------------------------------------------------
      else if(fDecay == 5 || fDecay == 6){ 
	if     (useVar == 0) histo4->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),  myWeightMVA);
	else if(useVar == 1) histo4->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),  myWeightMVA);
	else if(useVar == 2) histo4->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),    myWeightMVA);
	else if(useVar == 3) histo4->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),  myWeightMVA);
	else if(useVar == 4) histo4->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
	else if(useVar == 5) histo4->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
      }
      //----------------------------------------------------------------------------
      // For WW Bkg
      //----------------------------------------------------------------------------
      else if(fDecay == 0 || fDecay == 1){ // WW
	if     (useVar == 0) histo0->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),  myWeightMVA);
	else if(useVar == 1) histo0->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),  myWeightMVA);
	else if(useVar == 2) histo0->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),    myWeightMVA);
	else if(useVar == 3) histo0->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),  myWeightMVA);
	else if(useVar == 4) histo0->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
	else if(useVar == 5) histo0->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001),  myWeightMVA);
      }

      if(bdtg>-1.0&&bdtg<1){
        double myVar = dPhi*180.0/TMath::Pi(); myVar = CalcGammaMRstar(*lep1,*lep2);
        histoB->Fill(myVar,myWeightMVA);
      }

      //----------------------------------------------------------------------------
      // For all backgrounds together     
      //----------------------------------------------------------------------------
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

      //----------------------------------------------------------------------------
      // WARNING, THIS IS ONLY GOOD FOR BDTG!
      //----------------------------------------------------------------------------

      //----------------------------------------------------------------------------
      // W+Jet Bkg systematics:
      // Jet pT spectrum systematics shape 
      //----------------------------------------------------------------------------
      if(useWJetsTemplates == true && fDecay == 5){
        double addFRS=fakeRate(lep1->pt(), lep1->eta(), fhDFRMuSyst, fhDFRElSyst, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                								  (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	addFRS=addFRS*fakeRate(lep2->pt(), lep2->eta(), fhDFRMuSyst, fhDFRElSyst, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                								  (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        histo_Wjets_CMS_MVAWBoundingUp       ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addFRS/addFR);
      }

      if(useExpTemplates == true){
        //----------------------------------------------------------------------------
        // MVA Shape systematics for:
        // 1: Lepton efficiencies
        // 2: Lepton Resolution (bdtg_aux0,bdtg_aux1 branch)
        // 2: Met Resolution    (bdtg_aux2 branch)
        //----------------------------------------------------------------------------
        double addLepEffUp = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1,1)*
                             leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2,1);
        double addLepEffDown = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1,-1)*
                               leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2,-1);
        if     (fDecay == 0){
	  histo_qqWW_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_qqWW_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_qqWW_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_qqWW_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_qqWW_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(fDecay == 1){
	  histo_ggWW_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_ggWW_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_ggWW_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ggWW_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_ggWW_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(fDecay == 2){
	  histo_VV_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_VV_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_VV_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_VV_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_VV_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(fDecay == 3){
	  histo_Top_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Top_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Top_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(fDecay == 4){
	}
        else if(fDecay == 5){
	}
        else if(fDecay == 6 && useWgammaTemplates == false){
	  histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_Wgamma_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Wgamma_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
        else if(fDecay == 7){
	  histo_Ztt_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_Ztt_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_Ztt_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Ztt_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Ztt_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	}
      }

      //----------------------------------------------------------------------------
      // Add signal yields for MVA Cut analysis
      //----------------------------------------------------------------------------
      bool passFinalCut = false;
      if     (useVar == 0 && bdt        > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 1 && bdtd       > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 2 && nn         > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 3 && knn        > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      else if(useVar == 4 && bdtg       > useCut) {nBgdMVA += myWeight; nBgdEMVA += myWeight*myWeight; nBgdMVADecays[fDecay]  += myWeight; nBgdEMVADecays[fDecay]  += myWeight*myWeight; passFinalCut = true;}
      if(passFinalCut == true){
	if(fDecay == 4 && (dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
                          (type   == SmurfTree::mm   || type   == SmurfTree::ee)){     
	  DYXS[2] += myWeight;
	}
	else if(fDecay == 4){
	  VVXS[2] += myWeight;
	}
      }
    } // passMVAPreselCuts

    //----------------------------------------------------------------------------
    // MVA Shape systematics for:
    // Jet Energy Scale
    //----------------------------------------------------------------------------
    if(passMVAPreselCuts == true && useJESTemplates == true && (passJetCut[1] == true || passJetCut[2] == true)){
      if     (fDecay == 0){
        if(passJetCut[1] == true) histo_qqWW_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_qqWW_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(fDecay == 1){
        if(passJetCut[1] == true) histo_ggWW_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_ggWW_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(fDecay == 2){
        if(passJetCut[1] == true) histo_VV_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_VV_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(fDecay == 3){
       if(passJetCut[1] == true) histo_Top_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
       if(passJetCut[2] == true) histo_Top_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(fDecay == 4){
      }
      else if(fDecay == 5){
      }
      else if(fDecay == 6 && useWgammaTemplates == false){
        if(passJetCut[1] == true) histo_Wgamma_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_Wgamma_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
      else if(fDecay == 7){
        if(passJetCut[1] == true) histo_Ztt_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_Ztt_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
    } // End JES
  } // end loop over bkg events
  nBgdEAcc = sqrt(nBgdEAcc);
  nBgdECut = sqrt(nBgdECut);
  nBgdEMVA = sqrt(nBgdEMVA);
  printf("--- Finished Bgdnal loop\n");

  //****************************************************************************
  //
  // Loop Over Background Systematics Sample
  //
  //****************************************************************************

  treeSyst->SetBranchAddress( "cuts"          , &cuts 	  	);
  treeSyst->SetBranchAddress( "dstype"        , &dstype	  	);
  treeSyst->SetBranchAddress( "nvtx"          , &nvtx 	  	);
  treeSyst->SetBranchAddress( "npu"           , &npu 	        );
  treeSyst->SetBranchAddress( "njets"         , &njets	        );
  treeSyst->SetBranchAddress( "run"           , &run	        );
  treeSyst->SetBranchAddress( "event"         , &event	        );
  treeSyst->SetBranchAddress( "scale1fb"      , &scale1fb	);
  treeSyst->SetBranchAddress( "lep1"	      , &lep1		);
  treeSyst->SetBranchAddress( "lep2"	      , &lep2		);
  treeSyst->SetBranchAddress( "jet1"	      , &jet1		);
  treeSyst->SetBranchAddress( "jet2"	      , &jet2		);
  treeSyst->SetBranchAddress( "jet3"	      , &jet3		);
  treeSyst->SetBranchAddress( "dPhi"	      , &dPhi		);
  treeSyst->SetBranchAddress( "dR"	      , &dR		);
  treeSyst->SetBranchAddress( "dilep"	      , &dilep  	);
  treeSyst->SetBranchAddress( "type"	      , &type		);
  treeSyst->SetBranchAddress( "pmet"	      , &pmet		);
  treeSyst->SetBranchAddress( "pTrackMet"     , &pTrackMet	);
  treeSyst->SetBranchAddress( "met"	      , &met		);
  treeSyst->SetBranchAddress( "mt"	      , &mt		);
  treeSyst->SetBranchAddress( "mt1"	      , &mt1		);
  treeSyst->SetBranchAddress( "mt2"	      , &mt2		);
  treeSyst->SetBranchAddress( "dPhiLep1MET"   , &dPhiLep1MET	);
  treeSyst->SetBranchAddress( "dPhiLep2MET"   , &dPhiLep2MET	);
  treeSyst->SetBranchAddress( "dPhiDiLepMET"  , &dPhiDiLepMET	);
  treeSyst->SetBranchAddress( "dPhiDiLepJet1" , &dPhiDiLepJet1  );
  treeSyst->SetBranchAddress( "lq1"	      , &lq1		);
  treeSyst->SetBranchAddress( "lq2"	      , &lq2		);
  treeSyst->SetBranchAddress( "lid1"	      , &lid1		);
  treeSyst->SetBranchAddress( "lid2"	      , &lid2		);
  treeSyst->SetBranchAddress( "lid3"	      , &lid3		);
  treeSyst->SetBranchAddress( "processId"     , &processId	);
  treeSyst->SetBranchAddress( "jetLowBtag"    , &jetLowBtag	);
  treeSyst->SetBranchAddress( "nSoftMuons"    , &nSoftMuons	);
  treeSyst->SetBranchAddress( "jet1Btag"      , &jet1Btag	);
  treeSyst->SetBranchAddress( "jet2Btag"      , &jet2Btag	);
  treeSyst->SetBranchAddress( "jet3Btag"      , &jet3Btag	);
  treeSyst->SetBranchAddress( "lep1McId"      , &lep1McId	);
  treeSyst->SetBranchAddress( "lep2McId"      , &lep2McId	);
  treeSyst->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId );
  treeSyst->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId );
  treeSyst->SetBranchAddress(Form("bdt_hww%i_%djet_ww"      ,mH,nJetsType), &bdt	);
  treeSyst->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"     ,mH,nJetsType), &bdtd	);
  treeSyst->SetBranchAddress(Form("nn_hww%i_%djet_ww"	    ,mH,nJetsType), &nn 	);
  treeSyst->SetBranchAddress(Form("knn_hww%i_%djet_ww"      ,mH,nJetsType), &knn	);
  treeSyst->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"     ,mH,nJetsType), &bdtg	);
  treeSyst->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux0",mH,nJetsType), &bdtg_aux0  );
  treeSyst->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux1",mH,nJetsType), &bdtg_aux1  );
  treeSyst->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux2",mH,nJetsType), &bdtg_aux2  );
  //treeSyst->SetBranchAddress(Form("bdtg_hww%i_%djet_wjets"  ,mH,nJetsType), &bdtg_wjets );
  //treeSyst->SetBranchAddress(Form("knn_hww%i_%djet_wjets"   ,mH,nJetsType), &knn_wjets  );
  treeSyst->SetBranchAddress( "dymva"          , &dymva          );
  printf("--- Finished treeSyst loop\n");

  float nSystAcc  = 0.0;
  float nSystEAcc = 0.0;
  float nSystAccDecays[nChan+1]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float nSystEAccDecays[nChan+1] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for (UInt_t i=0; i<treeSyst->GetEntries(); i++) {

    treeSyst->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)treeSyst->GetEntries());

    //----------------------------------------------------------------------------
    //for data require that the event fired one of the designated signal triggers
    //----------------------------------------------------------------------------
    if(dstype == SmurfTree::data &&
      (cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dstype == SmurfTree::data && run <  minRun) continue;
    if(dstype == SmurfTree::data && run >  maxRun) continue;

    //----------------------------------------------------------------------------
    //define jet bin. 
    //For 2-Jet Bin:
    //  - apply "central jet veto"
    //    ie. no third jet above 30 GeV that lies between the 
    //    first two jets in pseudorapidity
    //----------------------------------------------------------------------------
    unsigned int Njet3 = njets;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(jet3->pt() <= 30)					           Njet3 = 2;
      else if(jet3->pt() > 30 && (
    	(jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    	(jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 0;
      else							           Njet3 = 2;
      if(njets < 2 || njets > 3)                                           Njet3 = 0;
    }

    //----------------------------------------------------------------------------
    //For Jet Energy scale systematics
    //----------------------------------------------------------------------------
    bool passJetCut[3] = {Njet3 == nJetsType, false, false};
    if(nJetsType == 0 && 			 jet1->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 4.7  			       ) passJetCut[1] = true;
    if(nJetsType == 0 && 			 jet1->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 4.7  			       ) passJetCut[2] = true;
    if(nJetsType == 1 && jet1->pt()*1.05 > 30 && jet2->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 4.7 && TMath::Abs(jet2->eta()) < 4.7) passJetCut[1] = true;
    if(nJetsType == 1 && jet1->pt()*0.95 > 30 && jet2->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 4.7 && TMath::Abs(jet2->eta()) < 4.7) passJetCut[2] = true;

    double minmet = TMath::Min(pmet,pTrackMet);
    bool passMET = minmet > 20.;
    if(useDYMVA == kFALSE){
      if     (njets == 0) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0   || type == SmurfTree::em || type == SmurfTree::me);
    } else {
      if     (njets == 0) passMET = passMET && (dymva >  0.60 || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (dymva >  0.30 || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0    || type == SmurfTree::em || type == SmurfTree::me);
    }

    bool passNewCuts = true;
    if(dilep->pt() <= 45) passNewCuts = false;

    //----------------------------------------------------------------------------
    //To create the DY bkg systematics MVA shapes (makeZjetsTemplates == true)
    //we loosen the MET cut to "minmet > 20 GeV". 
    //----------------------------------------------------------------------------
    if     ( makeZjetsTemplates == true && passMET == true ) passMET = false;
    else if( makeZjetsTemplates == true && passMET == false) passMET = minmet > 20.;

    //----------------------------------------------------------------------------
    // WW Preselection
    //----------------------------------------------------------------------------
    bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
                        ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
                         dstype != SmurfTree::data;

    if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
    if( passJetCut[0]==false&&passJetCut[1]==false&&passJetCut[2]==false ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= ptlMin	    					 ) continue; // cut on trailing lepton pt
    if( passNewCuts == false                                             ) continue; // cut on new pt cuts
    if( passMET == false                                                 ) continue; // cut on pmet
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJetCut = kTRUE;
    if(useDYMVA == kFALSE){
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    }
    if(njets >= 2) dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

    //----------------------------------------------------------------------------
    // Define Background Type
    // We have used a trick to resuse the dstype variable for systematics samples
    // 0 : MC@NLO qqWW nominal
    // 1 : MC@NLO qqWW Up
    // 8 : MC@NLO qqWW Down
    // 2 : WZ,ZZ
    // 3 : Top (ttbar, tW)
    // 4 : DY
    // 5 : W+Jets (fakes)
    // 6 : W+gamma , W+gamma*
    // 7 : DY->tautau
    // 
    //----------------------------------------------------------------------------
    int fDecay = 0;
    if     (dstype == SmurfTree::wjets  	 ) fDecay = 5;
    else if(dstype == SmurfTree::ttbar  	 ) fDecay = 3;
    else if(dstype == SmurfTree::dyee   	 ) fDecay = 4;
    else if(dstype == SmurfTree::dymm   	 ) fDecay = 4;
    else if(dstype == SmurfTree::dytt   	 ) fDecay = 7;
    else if(dstype == SmurfTree::ggzz   	 ) fDecay = 8;
    else if(dstype == SmurfTree::tw     	 ) fDecay = 3;
    else if(dstype == SmurfTree::qqww   	 ) fDecay = 0;
    else if(dstype == SmurfTree::wz     	 ) fDecay = 2;
    else if(dstype == SmurfTree::zz     	 ) fDecay = 2;
    else if(dstype == SmurfTree::www     	 ) fDecay = 2;
    else if(dstype == SmurfTree::ggww   	 ) fDecay = 1;
    else if(dstype == SmurfTree::wgamma 	 ) fDecay = 6;
    else if(dstype == SmurfTree::wgstar 	 ) fDecay = 6;
    else if(dstype == SmurfTree::data   	 ) fDecay = 5;
    else if(dstype == SmurfTree::dyttDataDriven  ) fDecay = 7;
    else if(dstype == SmurfTree::qcd             ) fDecay = 7;
    else                                 {printf("bad dstype: %d\n",dstype); assert(0);}
    if(dstype == SmurfTree::wz || dstype == SmurfTree::zz) {
      if(lep1MotherMcId == 23 && lep2MotherMcId == 23) {
        fDecay = 4;
      }
    }

    //----------------------------------------------------------------------------
    // Define Final States:
    // 5: Same Flavor
    // 6: Opposite Flavor
    //----------------------------------------------------------------------------
    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);
    if(wwDecayCut == false) continue;

    if(is2DAna == true) {
      double varA = (dilep->mass()-12.0)/(dilmass_cut-12.0);
      double varB = (mt-mtLowerCut)/(mtUpperCut-mtLowerCut);
      bdtg = Unroll2VarTo1Var(varA,varB,binVarA,binVarB,false);
    }
    if(bdtg < -1.0) bdtg = -0.999;
    if(bdtg > +1.0) bdtg =  0.999;

    //----------------------------------------------------------------------------
    // Define event weights    
    //----------------------------------------------------------------------------
    double myWeight = 1.0;
    double add      = 1.0;

    //----------------------------------------------------------------------------
    // First classify event into tight+tight, tight+fail, fail+fail 
    //----------------------------------------------------------------------------
    int nFake = 0;
    if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
    if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
    if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    if(category == 1) {
      if(((cuts & SmurfTree::Lep1LooseEleV1)  == SmurfTree::Lep1LooseEleV1)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake--;
      if(((cuts & SmurfTree::Lep2LooseEleV1)  == SmurfTree::Lep2LooseEleV1)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake--;
    }
    if(nFake < 0) assert(0);

    double addLepEff = 1.0;
    double addFR     = 1.0;
    bool isRealLepton = false;
    if((TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) &&
       (TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13)) isRealLepton = true;

    //----------------------------------------------------------------------------
    // Explicitly neglect fail+fail events
    //----------------------------------------------------------------------------
    if(nFake > 1){
      myWeight = 0.0;
    }

    //----------------------------------------------------------------------------
    // Tight+Fail events
    //----------------------------------------------------------------------------
    else if(nFake == 1){
      //----------------------------------------------------------------------------
      // For data, apply fake rate to generate bkg prediction
      //----------------------------------------------------------------------------
      if(dstype == SmurfTree::data){
        addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        							          (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        							          (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        add = addFR;
	if(category == 1) add = add*1.10; // HACK!!!
	fDecay  	      = 5;
        myWeight	      = add;
	// WE SHOULD NOT HAVE DATA EVENTS HERE
	assert(0);
      }

      //----------------------------------------------------------------------------
      // THIS IS SPECIAL FOR WGAMMA SYSTEMATICS OR WJETS SYSTEMATICS!!!
      // Apply nominal fake rates from data to W+jets MC
      //
      // This is not used at the moment. Eventually we want to do photon->electron
      // fake rates, but this is not implemented yet.
      //----------------------------------------------------------------------------
      else if((dstype == SmurfTree::wgamma && useWgammaTemplates == true ) ||
              (dstype == SmurfTree::wjets  && useWJetsMCTemplates == true)){
    	addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    		        						  (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    									  (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

        add = addFR;
    	add = add*nPUScaleFactor2012(fhDPU2012,npu);

        addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
                    leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
        add = add*addLepEff;

    	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
        if     (dstype == SmurfTree::wgamma) fDecay = 6;
	else if(dstype == SmurfTree::wjets)  fDecay = 5;
	else				     assert(0); 
	if(category == 1) add = add*1.10; // HACK!!!
        myWeight	       = 1.0 * scale1fb*scaleFactorLum*add;
      }
 
      //----------------------------------------------------------------------------
      // For real lepton or w+gamma:
      // apply fake rates, and give negative weight to subtract non jet-fake
      // contamination in the tight+fail sample
      //----------------------------------------------------------------------------
      else if(isRealLepton == true || (dstype == SmurfTree::wgamma && useWgammaTemplates == false)){
    	addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    		        						  (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    									  (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

        add = addFR;
    	add = add*nPUScaleFactor2012(fhDPU2012,npu);

        addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
                    leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
        add = add*addLepEff;

    	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
	if(category == 1) add = add*1.10; // HACK!!!
        fDecay  	       = 7; // SENT TO A GARBAGGE CATEGORY!!!
        myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
      }

      //----------------------------------------------------------------------------
      // Neglect fake lepton events in MC
      //----------------------------------------------------------------------------
      else {
        myWeight = 0.0;
      }
    }

    //----------------------------------------------------------------------------
    // Neglect any tight+tight events from data
    //----------------------------------------------------------------------------
    else if(dstype == SmurfTree::data) myWeight = 0.0;

    //----------------------------------------------------------------------------
    // DY->TauTau Embedded sample events: apply normalization weight.
    //----------------------------------------------------------------------------
    else if(dstype== SmurfTree::dyttDataDriven || dstype == SmurfTree::qcd) {
      myWeight = ZttScaleFactor(nvtx,period,scale1fb)*scaleFactorLum;
      assert(0);
    }

    //----------------------------------------------------------------------------
    // Regular Tight+Tight events from Monte Carlo
    //----------------------------------------------------------------------------
    else if(dstype != SmurfTree::data){
      //----------------------------------------------------------------------------      
      // Apply lepton efficiency scale factors, trigger efficiencies,
      // Pileup weights
      //----------------------------------------------------------------------------
      add = nPUScaleFactor2012(fhDPU2012,npu);

      addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
        	  leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
      add = add*addLepEff;
      add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							fabs(lep2->eta()), lep2->pt(), 
    							TMath::Abs( lid1), TMath::Abs(lid2));

      //----------------------------------------------------------------------------      
      // Apply DY Bkg Scale Factors
      //----------------------------------------------------------------------------
      if(fDecay == 4  && (type   == SmurfTree::mm   || type   == SmurfTree::ee)
                      && (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm)) {
    	if(njets == 0) add=add*DYBkgScaleFactor(0,0); 
    	if(njets == 1) add=add*DYBkgScaleFactor(0,1); 
    	if(njets >= 2) add=add*DYBkgScaleFactor(0,2); 
      }

      //----------------------------------------------------------------------------      
      // Apply Top Bkg Scale Factors
      //----------------------------------------------------------------------------
      if(fDecay == 3) {
    	if(njets == 0) add=add*TopBkgScaleFactor(0);
    	if(njets == 1) add=add*TopBkgScaleFactor(1); 
    	if(njets >= 2) add=add*TopBkgScaleFactor(2); 
      }

      //----------------------------------------------------------------------------      
      // Apply 0.0 scale factor to W+Jets
      //----------------------------------------------------------------------------
      if(fDecay == 5) add=add*0.0; 

      //----------------------------------------------------------------------------      
      // Apply W+gamma* normalization scale factor
      //----------------------------------------------------------------------------
      if(dstype == SmurfTree::wgstar) add=add*WGstarScaleFactor();

      //----------------------------------------------------------------------------      
      // Apply WW Bkg Scale Factors 
      // Don't do this for the WW selection (wwPresel == true)
      //----------------------------------------------------------------------------
      if((fDecay == 0 || fDecay == 1) && wwPresel == false){     
        if(njets == 0) add=add*WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0); 
        else	       add=add*WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1); 
      }
      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;
      myWeight = scale1fb*scaleFactorLum*add;
    }

    if(category == 1){
      bool passCuts = (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection ||
                      (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection;
      if(passCuts == false) myWeight = 0;
    }
 
    //----------------------------------------------------------------------------      
    // Explicitly neglect events with 0 weight
    //----------------------------------------------------------------------------

    if(myWeight == 0) continue;

    //----------------------------------------------------------------------------
    // Add bkg yields and fill MVA output distributions for MVA Shape analysis
    // Apply mT_Higgs cut for MVA Shape analysis
    //----------------------------------------------------------------------------
    bool passMVAPreselCuts = mt > mtLowerCut && mt < mtUpperCut; if(wwPresel == true) passMVAPreselCuts = true;
    if(is2DAna == true) passMVAPreselCuts = passMVAPreselCuts && lep1->pt() > theCutPtMaxLow;
    if(passMVAPreselCuts == true && passJetCut[0] == true){
      nSystAcc  = nSystAcc  + myWeight;
      nSystEAcc = nSystEAcc + myWeight*myWeight;
      nSystAccDecays[fDecay]  = nSystAccDecays[fDecay]  + myWeight;
      nSystEAccDecays[fDecay] = nSystEAccDecays[fDecay] + myWeight*myWeight;

      //----------------------------------------------------------------------------
      //systematics shapes for WW bkg (MC@NLO nominal, up, and down)
      //----------------------------------------------------------------------------
      if(useWWTemplates == true){
        if     (fDecay == 1) histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);
	else if(fDecay == 8) histo_qqWW_CMS_MVAWWNLOBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);
	else if(fDecay == 0) histo_qqWW_CMS_MVAWWBoundingUp     ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);
      }

      //----------------------------------------------------------------------------
      //systematics shapes for Top bkg (Madgraph ttbar+jets)
      //----------------------------------------------------------------------------
      if(useTopTemplates == true && fDecay == 3){
        histo_Top_CMS_MVATopBoundingUp->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);
      }

      //----------------------------------------------------------------------------
      //systematics shapes for WJets bkg (Madgraph W+jets MC)
      //----------------------------------------------------------------------------
      if(useWJetsMCTemplates == true && fDecay == 5){
        histo_Wjets_CMS_hww_MVAWMCBoundingUp->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);
      }

      //----------------------------------------------------------------------------
      //systematics shapes for Wgamma bkg
      //----------------------------------------------------------------------------
      if(useWgammaTemplates == true && fDecay == 6){
        histo_Wgamma->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),myWeight);
	if(useExpTemplates == true){
          double addLepEffUp = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1,1)*
          		       leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2,1);
          double addLepEffDown = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1,-1)*
          			 leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2,-1);
	  histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffUp  /addLepEff);
	  histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight*addLepEffDown/addLepEff);
	  histo_Wgamma_CMS_hww_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Wgamma_CMS_hww_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
	  histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        }
      }
    } // passMVAPreselCuts

    //----------------------------------------------------------------------------
    //Jet energy scale systematics shapes for Wgamma bkg
    //----------------------------------------------------------------------------
    if(passMVAPreselCuts == true && useJESTemplates == true && (passJetCut[1] == true || passJetCut[2] == true)){
      if(fDecay == 6 && useWgammaTemplates == true){
        if(passJetCut[1] == true) histo_Wgamma_CMS_hww_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
        if(passJetCut[2] == true) histo_Wgamma_CMS_hww_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      }
    }
  } // end loop over systematics bkg events

  nSystEAcc = sqrt(nSystEAcc);
  for(int i=0; i<nChan+1; i++) nSystEAccDecays[i] = sqrt(nSystEAccDecays[i]);
  printf("---\tacceptedSystPresel  %8.3f +/- %8.3f events\n",nSystAcc,nSystEAcc);
  printf("              qqww     ggww     VV       top      dyll     wjets    vg       Ztautau    ggZZ\n");
  printf("CLsSystAcc : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSystAccDecays[0],nSystAccDecays[1],nSystAccDecays[2],nSystAccDecays[3],nSystAccDecays[4],nSystAccDecays[5],nSystAccDecays[6],nSystAccDecays[7],nSystAccDecays[8]);
  printf("CLsSystEAcc: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",nSystEAccDecays[0],nSystEAccDecays[1],nSystEAccDecays[2],nSystEAccDecays[3],nSystEAccDecays[4],nSystEAccDecays[5],nSystEAccDecays[6],nSystEAccDecays[7],nSystEAccDecays[8]);

  //****************************************************************************
  //
  // Loop Over Data Sample
  //
  //****************************************************************************
  data->SetBranchAddress( "cuts"         , &cuts         );
  data->SetBranchAddress( "dstype"       , &dstype	 );
  data->SetBranchAddress( "nvtx"         , &nvtx         );
  data->SetBranchAddress( "npu"          , &npu          );
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
  data->SetBranchAddress( "jet3Btag"	 , &jet3Btag	 );
  data->SetBranchAddress(Form("bdt_hww%i_%djet_ww"	,mH,nJetsType), &bdt	    );
  data->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"	,mH,nJetsType), &bdtd	    );
  data->SetBranchAddress(Form("nn_hww%i_%djet_ww"	,mH,nJetsType), &nn	    );
  data->SetBranchAddress(Form("knn_hww%i_%djet_ww"	,mH,nJetsType), &knn	    );
  data->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"	,mH,nJetsType), &bdtg	    );
  //data->SetBranchAddress(Form("bdtg_hww%i_%djet_wjets"  ,mH,nJetsType), &bdtg_wjets );
  //data->SetBranchAddress(Form("knn_hww%i_%djet_wjets"   ,mH,nJetsType), &knn_wjets  );
  data->SetBranchAddress( "sfWeightPU"     , &sfWeightPU     );
  data->SetBranchAddress( "sfWeightEff"    , &sfWeightEff    );
  data->SetBranchAddress( "sfWeightTrig"   , &sfWeightTrig   );
  data->SetBranchAddress( "sfWeightHPt"    , &sfWeightHPt    );
  data->SetBranchAddress( "dymva"          , &dymva          );

  float nDatAcc = 0.0;
  float nDatCut = 0.0;
  float nDatMVA = 0.0;
  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    bool lId = (cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
    if(category == 1) lId = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2LooseEleV1)    == SmurfTree::Lep2LooseEleV1   ) ||
    			    ((cuts & SmurfTree::Lep1LooseEleV1)    == SmurfTree::Lep1LooseEleV1    && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection);
    if(!lId) continue;

    //----------------------------------------------------------------------------
    //for data require that the event fired one of the designated signal triggers
    //----------------------------------------------------------------------------
    if(dstype == SmurfTree::data &&
      (cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dstype == SmurfTree::data && run <  minRun) continue;
    if(dstype == SmurfTree::data && run >  maxRun) continue;

    //----------------------------------------------------------------------------
    //for data require that the event fired one of the designated signal triggers
    //----------------------------------------------------------------------------
    unsigned int Njet3 = njets;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(jet3->pt() <= 30)					           Njet3 = 2;
      else if(jet3->pt() > 30 && (
    	(jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    	(jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 0;
      else							           Njet3 = 2;
      if(njets < 2 || njets > 3)                                           Njet3 = 0;
    }
    bool passJetCut[1] = {Njet3 == nJetsType};

    double minmet = TMath::Min(pmet,pTrackMet);
    bool passMET = minmet > 20.;
    if(useDYMVA == kFALSE){
      if     (njets == 0) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (minmet > 45. || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0   || type == SmurfTree::em || type == SmurfTree::me);
    } else {
      if     (njets == 0) passMET = passMET && (dymva >  0.60 || type == SmurfTree::em || type == SmurfTree::me);
      else if(njets == 1) passMET = passMET && (dymva >  0.30 || type == SmurfTree::em || type == SmurfTree::me);
      else                passMET = passMET && (met > 45.0    || type == SmurfTree::em || type == SmurfTree::me);
    }

    bool passNewCuts = true;
    if(dilep->pt() <= 45) passNewCuts = false;

    //----------------------------------------------------------------------------
    //To create the DY bkg systematics MVA shapes (makeZjetsTemplates == true)
    //we loosen the MET cut to "minmet > 20 GeV". 
    //----------------------------------------------------------------------------
    if     ( makeZjetsTemplates == true && passMET == true ) passMET = false;
    else if( makeZjetsTemplates == true && passMET == false) passMET = minmet > 20.;

    //----------------------------------------------------------------------------
    // WW Preselection
    //----------------------------------------------------------------------------
    if( passJetCut[0]==false                                             ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= ptlMin	    					 ) continue; // cut on trailing lepton pt
    if( passNewCuts == false                                             ) continue; // cut on new pt cuts
    if( passMET == false                                                 ) continue; // cut on pmet
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    //if( lid3 != 0	                                                 ) continue; // cut on dileptons
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJetCut = kTRUE;
    if(useDYMVA == kFALSE){
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    }
    if(njets >= 2) dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
    if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

    //----------------------------------------------------------------------------
    // Define Final States:
    // 5: Same Flavor
    // 6: Opposite Flavor
    //----------------------------------------------------------------------------
    bool wwDecayCut = true;
    if     (wwDecay == 0) wwDecayCut = (type == SmurfTree::mm);
    else if(wwDecay == 1) wwDecayCut = (type == SmurfTree::ee);
    else if(wwDecay == 2) wwDecayCut = (type == SmurfTree::me);
    else if(wwDecay == 3) wwDecayCut = (type == SmurfTree::em);
    else if(wwDecay == 5) wwDecayCut = (type == SmurfTree::mm || type == SmurfTree::ee);
    else if(wwDecay == 6) wwDecayCut = (type == SmurfTree::me || type == SmurfTree::em);

    if(wwDecay == 5 && makeZjetsTemplates == true) wwDecayCut = true; // we use em events too

    if(wwDecayCut == false) continue;

    if(is2DAna == true) {
      double varA = (dilep->mass()-12.0)/(dilmass_cut-12.0);
      double varB = (mt-mtLowerCut)/(mtUpperCut-mtLowerCut);
      bdtg = Unroll2VarTo1Var(varA,varB,binVarA,binVarB,false);
    }
    if(bdtg < -1.0) bdtg = -0.999;
    if(bdtg > +1.0) bdtg =  0.999;

    //----------------------------------------------------------------------------
    // Weights for signal injection study
    //----------------------------------------------------------------------------
    double myWeight = 1.0;
    if(signalInjection == true) {
      myWeight = sfWeightPU*sfWeightEff*sfWeightTrig*sfWeightHPt*scaleFactorLum*scale1fb;
    }
    if(wwDecay == 5 && makeZjetsTemplates == true && (type == SmurfTree::me || type == SmurfTree::em)) myWeight = -1.0;

    if(myWeight == 0) continue;

    //----------------------------------------------------------------------------
    //
    // Higgs Signal Selection Cuts
    //
    //----------------------------------------------------------------------------
    double theCutPtMinLow = cutPtMinLow (mH, type);
    bool passAllCuts = dilep->mass()         < theCutMassHigh &&
        	       mt		     > theCutMTLow &&
        	       mt		     < theCutMTHigh &&
        	       lep1->pt()	     > theCutPtMaxLow &&
        	       lep2->pt()	     > theCutPtMinLow &&
        	       dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh;

    //----------------------------------------------------------------------------
    // VBF selection cuts for 2-Jet bin
    //----------------------------------------------------------------------------
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    mt > 30 && mt < mtUpperCut &&
		    centrality == 1;
      if(is2DAna == true) {
        passAllCuts = passAllCuts &&
                      dilep->mass()         < theCutMassHigh &&
        	      mt 		    < theCutMTHigh &&
        	      lep1->pt() 	    > theCutPtMaxLow &&
        	      lep2->pt() 	    > theCutPtMinLow &&
        	      dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh;
      }
    }

    //----------------------------------------------------------------------------
    // Data yields for Cut-Based analysis
    //----------------------------------------------------------------------------
    if(passAllCuts == true) {
      nDatCut = nDatCut + myWeight;
      datMVA[5]->Fill(TMath::Max(TMath::Min((double)mt,maxHis[5]-0.001),minHis[5]+0.001), myWeight);
    }

    //----------------------------------------------------------------------------
    // Data yields and fill MVA output distributions for MVA Shape analysis
    // Apply mT_Higgs cut for MVA Shape analysis
    //----------------------------------------------------------------------------
    bool passMVAPreselCuts = mt > mtLowerCut && mt < mtUpperCut; if(wwPresel == true) passMVAPreselCuts = true;
    if(is2DAna == true) passMVAPreselCuts = passMVAPreselCuts && lep1->pt() > theCutPtMaxLow;
    if(passMVAPreselCuts == true){
      nDatAcc = nDatAcc + myWeight;
      if(bdtg>-1.0&&bdtg<1){
        double myVar = dPhi*180.0/TMath::Pi(); myVar = CalcGammaMRstar(*lep1,*lep2);
        histoD->Fill(myVar,myWeight);
      }
      if     (useVar == 0) histo5->Fill(TMath::Max(TMath::Min((double)bdt ,maxHis[0]-0.001),minHis[0]+0.001),       myWeight);
      else if(useVar == 1) histo5->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),       myWeight);
      else if(useVar == 2) histo5->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),         myWeight);
      else if(useVar == 3) histo5->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),        myWeight);
      else if(useVar == 4) histo5->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),       myWeight);
      else if(useVar == 5) histo5->Fill(TMath::Max(TMath::Min((double)bdtg_wjets,maxHis[4]-0.001),minHis[4]+0.001), myWeight);
      //printf("SSS %6d %15d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %1d %6.2f %6.2f\n",run,event,bdtg,lep1->pt(),lep2->pt(),dPhi,dR,dilep->mass(),mt,type,dPhiDiLepMET,dPhiDiLepJet1);
      datMVA[0]->Fill(TMath::Max(TMath::Min((double)bdt,maxHis[0]-0.001),minHis[0]+0.001),    myWeight);
      datMVA[1]->Fill(TMath::Max(TMath::Min((double)bdtd,maxHis[1]-0.001),minHis[1]+0.001),   myWeight);
      datMVA[2]->Fill(TMath::Max(TMath::Min((double)nn,maxHis[2]-0.001),minHis[2]+0.001),     myWeight);
      datMVA[3]->Fill(TMath::Max(TMath::Min((double)knn,maxHis[3]-0.001),minHis[3]+0.001),    myWeight);
      datMVA[4]->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis[4]-0.001),minHis[4]+0.001),   myWeight);

      //----------------------------------------------------------------------------
      // Data yields for MVA Cut analysis
      //----------------------------------------------------------------------------
      bool passFinalCut = false;
      if     (useVar == 0 && bdt        > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 1 && bdtd       > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 2 && nn         > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 3 && knn        > useCut) {nDatMVA += myWeight; passFinalCut = true;}
      else if(useVar == 4 && bdtg       > useCut) {nDatMVA += myWeight; passFinalCut = true;}
    } // passMVAPreselCuts
  } // end loop over data events

  // remove VV component when making Zjets templates
  if(wwDecay == 5 && makeZjetsTemplates == true) {
    printf("Z data template data-VV: %f - %f = ",histo5->GetSumOfWeights(),histoVV->GetSumOfWeights());
    histo5->Add(histoVV,-1.0);
    printf("%f\n",histo5->GetSumOfWeights());
  }

  // cut-based 0/1 jet bin analyses
  if(nJetsType != 2 && mH <= 300 && (wwDecay == 4 || wwDecay == 5)){
    DYXS[1] = DYBkgScaleFactor(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2));
    nBgdCutDecays[4] += DYXS[1];
    nBgdCut += DYXS[1];
  }
  if(nJetsType == 2 && (wwDecay == 4 || wwDecay == 5)){
    DYXS[1] = DYBkgScaleFactor(TMath::Min(TMath::Max((int)mH,115),300),TMath::Min((int)nJetsType,2));
    nBgdCutDecays[4] += DYXS[1];
    nBgdCut += DYXS[1];
  }

  // shape-based 0/1 jet bin analyses
  if(nJetsType != 2 && mH <= 300 && (wwDecay == 4 || wwDecay == 5) && makeZjetsTemplates == false){
    DYXS[0] = DYBkgScaleFactorBDT(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2));
  }
  nBgdAccDecays[4] += DYXS[0];
  nBgdAcc += DYXS[0];
  //****************************************************************************
  //
  // Print Summary Information
  //
  //****************************************************************************

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

  //----------------------------------------------------------------------------
  // Systematic Uncertainties for DY
  // VERY CAREFUL HERE!!!, if we don't use bdtg or change binning!,it has to be changed
  // 
  //----------------------------------------------------------------------------
  if(useZjetsTemplates == true){
    //----------------------------------------------------------------------------
    // hDZjetsMCTemplate is the nominal MVA shape for the DY bkg
    // normalize the histogram to the predicted DY background yield
    // hDZjetsDATemplate is the alternative MVA shape for the DY bkg
    // normalize the histogram to the predicted DY background yield
    // 
    // histo_Zjets_CMS_MVAZBounding   : nominal MVA shape
    //----------------------------------------------------------------------------
    hDZjetsMCTemplate->Scale(DYXS[0]);
    hDZjetsDATemplate->Scale(DYXS[0]);
    histo_Zjets_CMS_MVAZBounding  ->Add(hDZjetsMCTemplate);
    histo_Zjets_CMS_MVAZBoundingUp->Add(hDZjetsDATemplate);

    //rebin the MVA shape histograms
    histo_Zjets_CMS_MVAZBounding->Rebin(rebinMVAHist);
    histo_Zjets_CMS_MVAZBoundingUp->Rebin(rebinMVAHist);
    histo_Zjets_CMS_MVAZBoundingDown->Rebin(rebinMVAHist);

    //----------------------------------------------------------------------------
    // Construct Opposite Bounding Shape
    // histo_Zjets_CMS_MVAZBoundingDown: mirror the difference between nominal 
    //                                   shape and the systematics shape
    //----------------------------------------------------------------------------
    for(int i=1; i<=histo_Zjets_CMS_MVAZBounding->GetNbinsX(); i++){
      double mean = histo_Zjets_CMS_MVAZBounding  ->GetBinContent(i);
      double up   = histo_Zjets_CMS_MVAZBoundingUp->GetBinContent(i);
      double diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_Zjets_CMS_MVAZBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else                histo_Zjets_CMS_MVAZBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    }
    if(histo_Zjets_CMS_MVAZBoundingDown->GetSumOfWeights() > 0)
       histo_Zjets_CMS_MVAZBoundingDown->Scale(histo_Zjets_CMS_MVAZBoundingUp->GetSumOfWeights()/histo_Zjets_CMS_MVAZBoundingDown->GetSumOfWeights());

    //----------------------------------------------------------------------------
    // Add in the WZ,ZZ to the "Zjets" MVA shape, 
    // but leave them as the nominal shape
    //----------------------------------------------------------------------------
    printf("VV norm %f\n",histoVV->GetSumOfWeights());
    histoVV->Rebin(rebinMVAHist);
    histo_Zjets_CMS_MVAZBoundingUp  ->Add(histoVV);
    histo_Zjets_CMS_MVAZBoundingDown->Add(histoVV);

    bgdMVA[4]->Add(hDZjetsMCTemplate);
    bgdMVADecays[4][4]->Add(hDZjetsMCTemplate);
    histo1->Add(hDZjetsMCTemplate);
    printf("Zjets norm up/default/down: %f/%f/%f\n",histo_Zjets_CMS_MVAZBoundingUp->GetSumOfWeights(),bgdMVADecays[4][4]->GetSumOfWeights(),histo_Zjets_CMS_MVAZBoundingDown->GetSumOfWeights());
  }

  //----------------------------------------------------------------------------
  // Compute DY Bkg Uncertainties on the yield
  //----------------------------------------------------------------------------
  if(nBgdAccDecays[4] > 0 && makeZjetsTemplates == false && TMath::Abs(DYXS[0]+VVXS[0]-nBgdAccDecays[4])/nBgdAccDecays[4] > 0.001) {printf("Problem: %f %f %f\n",DYXS[0],VVXS[0],nBgdAccDecays[4]); assert(0);}
  if(nBgdCutDecays[4] > 0 && makeZjetsTemplates == false && TMath::Abs(DYXS[1]+VVXS[1]-nBgdCutDecays[4])/nBgdCutDecays[4] > 0.001) {printf("Problem: %f %f %f\n",DYXS[1],VVXS[1],nBgdCutDecays[4]); assert(0);}
  if(nBgdMVADecays[4] > 0 && makeZjetsTemplates == false && TMath::Abs(DYXS[2]+VVXS[2]-nBgdMVADecays[4])/nBgdMVADecays[4] > 0.001) {printf("Problem: %f %f %f\n",DYXS[2],VVXS[2],nBgdMVADecays[4]); assert(0);}
  int MhForError = mH;
  if(mH >= 110 && mH < 115) MhForError = 115;
  if(mH > 300) MhForError = 0;
  if(nJetsType == 2 && mH > 300) MhForError = 300;
  if(nBgdAccDecays[4] > 0.0) ZXS_E[0] = sqrt(DYXS[0]*DYXS[0]*(DYBkgScaleFactorBDTKappa(MhForError,TMath::Min((int)nJetsType,2))-1.0)*(DYBkgScaleFactorBDTKappa(MhForError,TMath::Min((int)nJetsType,2))-1.0)+
                                             VVXS[0]*VVXS[0]*0.10*0.10)/nBgdAccDecays[4]; else ZXS_E[0] = 0;
  if(nBgdCutDecays[4] > 0.0) ZXS_E[1] = sqrt(DYXS[1]*DYXS[1]*(DYBkgScaleFactorKappa   (MhForError,TMath::Min((int)nJetsType,2))-1.0)*(DYBkgScaleFactorKappa   (MhForError,TMath::Min((int)nJetsType,2))-1.0)+
                                             VVXS[1]*VVXS[1]*0.10*0.10)/nBgdCutDecays[4]; else ZXS_E[1] = 0;
  if(nBgdMVADecays[4] > 0.0) ZXS_E[2] = sqrt(DYXS[2]*DYXS[2]*(DYBkgScaleFactorKappa   (MhForError,TMath::Min((int)nJetsType,2))-1.0)*(DYBkgScaleFactorKappa   (MhForError,TMath::Min((int)nJetsType,2))-1.0)+
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
    TH1D* histo_ttH    = new TH1D("histo_ttH",	"histo_ttH",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_ZH     = new TH1D("histo_ZH",	"histo_ZH",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_WH     = new TH1D("histo_WH",	"histo_WH",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_qqH    = new TH1D("histo_qqH",	"histo_qqH",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_ggH    = new TH1D("histo_ggH",	"histo_ggH",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_Data   = new TH1D("histo_Data",	"histo_Data",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_qqWW   = new TH1D("histo_qqWW",	"histo_qqWW",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_ggWW   = new TH1D("histo_ggWW",	"histo_ggWW",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_VV     = new TH1D("histo_VV",	"histo_VV",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_Top    = new TH1D("histo_Top",	"histo_Top",	nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_Zjets  = new TH1D("histo_Zjets","histo_Zjets",  nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_Wjets  = new TH1D("histo_Wjets","histo_Wjets",  nBinHis, minHis[1], maxHis[1]);
    TH1D* histo_Ztt    = new TH1D("histo_Ztt",	"histo_Ztt",	nBinHis, minHis[1], maxHis[1]);

    histo_ttH   ->Add(sigMVA[useVar][1]);
    histo_ZH	->Add(sigMVA[useVar][2]);
    histo_WH	->Add(sigMVA[useVar][3]);
    histo_qqH	->Add(sigMVA[useVar][4]);
    histo_ggH	->Add(sigMVA[useVar][5]);
    histo_Data  ->Add(datMVA[useVar]);
    histo_qqWW  ->Add(bgdMVADecays[useVar][0]);
    histo_ggWW  ->Add(bgdMVADecays[useVar][1]);
    histo_VV	->Add(bgdMVADecays[useVar][2]);
    histo_Top	->Add(bgdMVADecays[useVar][3]);
    histo_Zjets ->Add(bgdMVADecays[useVar][4]);
    histo_Wjets ->Add(bgdMVADecays[useVar][5]);
    histo_Ztt	->Add(bgdMVADecays[useVar][7]);

    if(useWgammaTemplates == true){
      double scaleWg = histo_Wgamma->GetSumOfWeights();
      if(scaleWg > 0){
      	histo_Wgamma->Scale(bgdMVADecays[useVar][6]->GetSumOfWeights()/scaleWg);
      	if(useExpTemplates == true){ // Allow for different normalization
      	  histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ->Scale(bgdMVADecays[useVar][6]->GetSumOfWeights()/scaleWg);
      	  histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Scale(bgdMVADecays[useVar][6]->GetSumOfWeights()/scaleWg);
      	  histo_Wgamma_CMS_hww_MVALepResBoundingUp  ->Scale(bgdMVADecays[useVar][6]->GetSumOfWeights()/scaleWg);
      	  histo_Wgamma_CMS_hww_MVALepResBoundingDown->Scale(bgdMVADecays[useVar][6]->GetSumOfWeights()/scaleWg);
      	  histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ->Scale(bgdMVADecays[useVar][6]->GetSumOfWeights()/scaleWg);
      	}
      }
      else {
      	histo_Wgamma->Add(bgdMVADecays[useVar][6]);
      	if(useExpTemplates == true){ // Allow for different normalization
      	  histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ->Add(bgdMVADecays[useVar][6]);
      	  histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Add(bgdMVADecays[useVar][6]);
      	  histo_Wgamma_CMS_hww_MVALepResBoundingUp  ->Add(bgdMVADecays[useVar][6]);
      	  histo_Wgamma_CMS_hww_MVALepResBoundingDown->Add(bgdMVADecays[useVar][6]);
      	  histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ->Add(bgdMVADecays[useVar][6]);
      	}
      }
    }
    else {
      histo_Wgamma->Add(bgdMVADecays[useVar][6]);
    }
    
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

    if(rebinSmurf > 1){
       histSmurfRebin = SmurfRebin(histo_Wgamma ,rebinSmurf);histo_Wgamma->Scale(0.0);histo_Wgamma->Add(histSmurfRebin);
       if(wwDecay == 6){
         histSmurfRebin = SmurfRebin(histo_Zjets,rebinSmurf);histo_Zjets ->Scale(0.0);histo_Zjets ->Add(histSmurfRebin);   
       }
    }

    for(int i=1; i<=histo_Wjets->GetNbinsX(); i++){
      if(histo_Wjets->GetBinContent(i) < 0) histo_Wjets->SetBinContent(i,0.000001);
    }
    // We need to renormalize
    if(bgdMVADecays[useVar][5]->GetSumOfWeights() > 0) histo_Wjets->Scale(bgdMVADecays[useVar][5]->GetSumOfWeights()/histo_Wjets ->GetSumOfWeights());

    if(signalInjection == true){
      for(int i=1; i<=histo_Data->GetNbinsX(); i++){
        double SplusB = histo_Data->GetBinContent(i)+
	                histo_qqWW->GetBinContent(i)+histo_ggWW->GetBinContent(i)+
	                histo_VV->GetBinContent(i)+histo_Top->GetBinContent(i)+
	                histo_Zjets->GetBinContent(i)+histo_Wjets->GetBinContent(i)+
	                histo_Wgamma->GetBinContent(i)+histo_Ztt->GetBinContent(i);
        histo_Data->SetBinContent(i,(int)(SplusB+0.5));
        histo_Data->SetBinError(i,sqrt((int)(SplusB+0.5)));
      }
    }

    //----------------------------------------------------------------------------
    //
    // Fill Data Cards and generate MVA shape histograms
    // 
    //----------------------------------------------------------------------------
    char outputLimits[200];
    sprintf(outputLimits,"hww%s_%dj.input_8TeV.root",finalStateName,nJetsType);
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();

    //----------------------------------------------------------------------------
    // Nominal Shapes
    //----------------------------------------------------------------------------
    histo_ggH	 ->Write();
    histo_qqH	 ->Write();
    histo_WH	 ->Write();
    histo_ZH	 ->Write();
    histo_ttH	 ->Write();
    histo_Data   ->Write();
    histo_qqWW   ->Write();
    histo_ggWW   ->Write();
    histo_VV	 ->Write();
    histo_Top	 ->Write();
    histo_Zjets  ->Write();
    histo_Wjets  ->Write();
    histo_Wgamma ->Write();
    histo_Ztt	 ->Write();
    cout << histo_ttH->GetSumOfWeights() << " ";
    cout << histo_ZH ->GetSumOfWeights() << " ";
    cout << histo_WH ->GetSumOfWeights() << " ";
    cout << histo_qqH->GetSumOfWeights() << " ";
    cout << histo_ggH->GetSumOfWeights() << " ";
    cout << histo_qqWW  ->GetSumOfWeights() << " ";
    cout << histo_ggWW  ->GetSumOfWeights() << " ";
    cout << histo_VV    ->GetSumOfWeights() << " ";
    cout << histo_Top   ->GetSumOfWeights() << " ";
    cout << histo_Zjets ->GetSumOfWeights() << " ";
    cout << histo_Wjets ->GetSumOfWeights() << " ";
    cout << histo_Wgamma->GetSumOfWeights() << " ";
    cout << histo_Ztt   ->GetSumOfWeights() << endl;
    //----------------------------------------------------------------------------
    // DY Bkg systematics shapes
    //----------------------------------------------------------------------------
    if(useZjetsTemplates == true){
      histo_Zjets_CMS_MVAZBoundingUp->Write();
      histo_Zjets_CMS_MVAZBoundingDown->Write();
    }

    //----------------------------------------------------------------------------
    // WW Bkg systematics shapes
    // Rebin, normalize histograms to nominal WW bkg yield
    // MVAWWBounding : mirror the difference between nominal shape (Madgraph) and 
    //                 systematics shape (MC@NLO)
    //----------------------------------------------------------------------------
    if(useWWTemplates == true){
      histo_qqWW_CMS_MVAWWBoundingUp  ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_MVAWWBoundingDown->Rebin(rebinMVAHist);
      histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_MVAWWNLOBoundingDown->Rebin(rebinMVAHist);
      histo_qqWW_CMS_MVAWWBoundingUp     ->Scale(histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWWBoundingUp     ->GetSumOfWeights());
      histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Scale(histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetSumOfWeights());
      histo_qqWW_CMS_MVAWWNLOBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWWNLOBoundingDown->GetSumOfWeights());
      for(int i=1; i<=histo_qqWW->GetNbinsX(); i++){
        double mean = histo_qqWW  ->GetBinContent(i);
        double up   = histo_qqWW_CMS_MVAWWBoundingUp->GetBinContent(i);
        double diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_qqWW_CMS_MVAWWBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_qqWW_CMS_MVAWWBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        double meanNLO = histo_qqWW_CMS_MVAWWBoundingUp->GetBinContent(i);
	if(meanNLO > 0) histo_qqWW_CMS_MVAWWNLOBoundingUp  ->SetBinContent(i,histo_qqWW->GetBinContent(i)*histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetBinContent(i)/meanNLO);
	if(meanNLO > 0) histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinContent(i,histo_qqWW->GetBinContent(i)*histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinContent(i)/meanNLO);
      }
      histo_qqWW_CMS_MVAWWBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWWBoundingDown->GetSumOfWeights());
      histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Scale(histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetSumOfWeights());
      histo_qqWW_CMS_MVAWWNLOBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWWNLOBoundingDown->GetSumOfWeights());

      histo_qqWW_CMS_MVAWWBoundingUp  ->Write();
      histo_qqWW_CMS_MVAWWBoundingDown->Write();
      histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Write();
      histo_qqWW_CMS_MVAWWNLOBoundingDown->Write();
    }

    //----------------------------------------------------------------------------
    // Top Bkg systematics shapes
    // Rebin, normalize histograms to nominal WW bkg yield
    // MVATopBounding : mirror the difference between nominal shape (POWHEG) and 
    //                 systematics shape (Madgraph)
    //----------------------------------------------------------------------------
    if(useTopTemplates == true){
      histo_Top_CMS_MVATopBoundingUp  ->Scale(histo_Top->GetSumOfWeights()/histo_Top_CMS_MVATopBoundingUp  ->GetSumOfWeights());
      histo_Top_CMS_MVATopBoundingUp  ->Rebin(rebinMVAHist);
      histo_Top_CMS_MVATopBoundingDown->Rebin(rebinMVAHist);
      for(int i=1; i<=histo_Top->GetNbinsX(); i++){
        double mean = histo_Top  ->GetBinContent(i);
        double up   = histo_Top_CMS_MVATopBoundingUp->GetBinContent(i);
        double diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Top_CMS_MVATopBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Top_CMS_MVATopBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
      }
      histo_Top_CMS_MVATopBoundingDown->Scale(histo_Top->GetSumOfWeights()/histo_Top_CMS_MVATopBoundingDown->GetSumOfWeights());
      histo_Top_CMS_MVATopBoundingUp  ->Write();
      histo_Top_CMS_MVATopBoundingDown->Write();
    }

    //----------------------------------------------------------------------------
    // Systematics Shapes due to MC statistics
    // Left 1/3 increases by 1sigma, right 1/3 decreases by 1sigma
    //----------------------------------------------------------------------------
    if(useStatTemplates == true){
      histo_ttH_CMS_MVAttHStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_ttH_CMS_MVAttHStatBounding_8TeVDown         ->Rebin(rebinMVAHist);
      histo_ZH_CMS_MVAZHStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_ZH_CMS_MVAZHStatBounding_8TeVDown	   ->Rebin(rebinMVAHist);
      histo_WH_CMS_MVAWHStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_WH_CMS_MVAWHStatBounding_8TeVDown	   ->Rebin(rebinMVAHist);
      histo_qqH_CMS_MVAqqHStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_qqH_CMS_MVAqqHStatBounding_8TeVDown         ->Rebin(rebinMVAHist);
      histo_ggH_CMS_MVAggHStatBounding_8TeVUp      ->Rebin(rebinMVAHist);
      histo_ggH_CMS_MVAggHStatBounding_8TeVDown         ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_MVAqqWWStatBounding_8TeVUp    ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_MVAqqWWStatBounding_8TeVDown       ->Rebin(rebinMVAHist);
      histo_ggWW_CMS_MVAggWWStatBounding_8TeVUp    ->Rebin(rebinMVAHist);
      histo_ggWW_CMS_MVAggWWStatBounding_8TeVDown       ->Rebin(rebinMVAHist);
      histo_VV_CMS_MVAVVStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_VV_CMS_MVAVVStatBounding_8TeVDown	   ->Rebin(rebinMVAHist);
      histo_Top_CMS_MVATopStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_Top_CMS_MVATopStatBounding_8TeVDown         ->Rebin(rebinMVAHist);
      histo_Zjets_CMS_MVAZjetsStatBounding_8TeVUp  ->Rebin(rebinMVAHist);
      histo_Zjets_CMS_MVAZjetsStatBounding_8TeVDown     ->Rebin(rebinMVAHist);
      histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp  ->Rebin(rebinMVAHist);
      histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown     ->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown   ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_MVAZttStatBounding_8TeVUp	   ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_MVAZttStatBounding_8TeVDown         ->Rebin(rebinMVAHist);
      if(histo_ttH_CMS_MVAttHStatBounding_8TeVUp->GetNbinsX() != histo_ttH->GetNbinsX()) {printf("PROBLEMMMMM\n");return;}
      for(int i=1; i<=histo_ttH->GetNbinsX(); i++){
        double factorUp = +1.0; double factorDown = -1.0;
	if(useAlternativeStatTemplates == true){
	  if     (i<1.0*histo_ttH->GetNbinsX()/3.0) {factorUp = -1.0; factorDown = +1.0;}
	  else if(i<2.0*histo_ttH->GetNbinsX()/3.0) {factorUp = +0.0; factorDown = +0.0;}
	  else                                      {factorUp = +1.0; factorDown = -1.0;}
	}
    	histo_ttH_CMS_MVAttHStatBounding_8TeVUp      ->SetBinContent(i,TMath::Max(histo_ttH   ->GetBinContent(i)+factorUp  *histo_ttH   ->GetBinError(i),0.000001));
    	histo_ttH_CMS_MVAttHStatBounding_8TeVDown         ->SetBinContent(i,TMath::Max(histo_ttH   ->GetBinContent(i)+factorDown*histo_ttH   ->GetBinError(i),0.000001));
    	histo_ZH_CMS_MVAZHStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_ZH    ->GetBinContent(i)+factorUp  *histo_ZH    ->GetBinError(i),0.000001));
    	histo_ZH_CMS_MVAZHStatBounding_8TeVDown           ->SetBinContent(i,TMath::Max(histo_ZH    ->GetBinContent(i)+factorDown*histo_ZH    ->GetBinError(i),0.000001));
    	histo_WH_CMS_MVAWHStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_WH    ->GetBinContent(i)+factorUp  *histo_WH    ->GetBinError(i),0.000001));
    	histo_WH_CMS_MVAWHStatBounding_8TeVDown           ->SetBinContent(i,TMath::Max(histo_WH    ->GetBinContent(i)+factorDown*histo_WH    ->GetBinError(i),0.000001));
    	histo_qqH_CMS_MVAqqHStatBounding_8TeVUp      ->SetBinContent(i,TMath::Max(histo_qqH   ->GetBinContent(i)+factorUp  *histo_qqH   ->GetBinError(i),0.000001));
    	histo_qqH_CMS_MVAqqHStatBounding_8TeVDown         ->SetBinContent(i,TMath::Max(histo_qqH   ->GetBinContent(i)+factorDown*histo_qqH   ->GetBinError(i),0.000001));
    	histo_ggH_CMS_MVAggHStatBounding_8TeVUp      ->SetBinContent(i,TMath::Max(histo_ggH   ->GetBinContent(i)+factorUp  *histo_ggH   ->GetBinError(i),0.000001));
    	histo_ggH_CMS_MVAggHStatBounding_8TeVDown         ->SetBinContent(i,TMath::Max(histo_ggH   ->GetBinContent(i)+factorDown*histo_ggH   ->GetBinError(i),0.000001));
    	histo_qqWW_CMS_MVAqqWWStatBounding_8TeVUp    ->SetBinContent(i,TMath::Max(histo_qqWW  ->GetBinContent(i)+factorUp  *histo_qqWW  ->GetBinError(i),0.000001));
    	histo_qqWW_CMS_MVAqqWWStatBounding_8TeVDown       ->SetBinContent(i,TMath::Max(histo_qqWW  ->GetBinContent(i)+factorDown*histo_qqWW  ->GetBinError(i),0.000001));
    	histo_ggWW_CMS_MVAggWWStatBounding_8TeVUp    ->SetBinContent(i,TMath::Max(histo_ggWW  ->GetBinContent(i)+factorUp  *histo_ggWW  ->GetBinError(i),0.000001));
    	histo_ggWW_CMS_MVAggWWStatBounding_8TeVDown       ->SetBinContent(i,TMath::Max(histo_ggWW  ->GetBinContent(i)+factorDown*histo_ggWW  ->GetBinError(i),0.000001));
    	histo_VV_CMS_MVAVVStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_VV    ->GetBinContent(i)+factorUp  *histo_VV    ->GetBinError(i),0.000001));
    	histo_VV_CMS_MVAVVStatBounding_8TeVDown           ->SetBinContent(i,TMath::Max(histo_VV    ->GetBinContent(i)+factorDown*histo_VV    ->GetBinError(i),0.000001));
    	histo_Top_CMS_MVATopStatBounding_8TeVUp      ->SetBinContent(i,TMath::Max(histo_Top   ->GetBinContent(i)+factorUp  *histo_Top   ->GetBinError(i),0.000001));
    	histo_Top_CMS_MVATopStatBounding_8TeVDown         ->SetBinContent(i,TMath::Max(histo_Top   ->GetBinContent(i)+factorDown*histo_Top   ->GetBinError(i),0.000001));
    	histo_Zjets_CMS_MVAZjetsStatBounding_8TeVUp  ->SetBinContent(i,TMath::Max(histo_Zjets ->GetBinContent(i)+factorUp  *histo_Zjets ->GetBinError(i),0.000001));
    	histo_Zjets_CMS_MVAZjetsStatBounding_8TeVDown     ->SetBinContent(i,TMath::Max(histo_Zjets ->GetBinContent(i)+factorDown*histo_Zjets ->GetBinError(i),0.000001));
    	histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp  ->SetBinContent(i,TMath::Max(histo_Wjets ->GetBinContent(i)+factorUp  *histo_Wjets ->GetBinError(i),0.000001));
    	histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown     ->SetBinContent(i,TMath::Max(histo_Wjets ->GetBinContent(i)+factorDown*histo_Wjets ->GetBinError(i),0.000001));
    	histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp->SetBinContent(i,TMath::Max(histo_Wgamma->GetBinContent(i)+factorUp  *histo_Wgamma->GetBinError(i),0.000001));
    	histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown   ->SetBinContent(i,TMath::Max(histo_Wgamma->GetBinContent(i)+factorDown*histo_Wgamma->GetBinError(i),0.000001));
    	histo_Ztt_CMS_MVAZttStatBounding_8TeVUp      ->SetBinContent(i,TMath::Max(histo_Ztt   ->GetBinContent(i)+factorUp  *histo_Ztt   ->GetBinError(i),0.000001));
    	histo_Ztt_CMS_MVAZttStatBounding_8TeVDown         ->SetBinContent(i,TMath::Max(histo_Ztt   ->GetBinContent(i)+factorDown*histo_Ztt   ->GetBinError(i),0.000001));
      }
      histo_ttH_CMS_MVAttHStatBounding_8TeVUp	    ->Write();
      histo_ttH_CMS_MVAttHStatBounding_8TeVDown          ->Write();
      histo_ZH_CMS_MVAZHStatBounding_8TeVUp	    ->Write();
      histo_ZH_CMS_MVAZHStatBounding_8TeVDown	    ->Write();
      histo_WH_CMS_MVAWHStatBounding_8TeVUp	    ->Write();
      histo_WH_CMS_MVAWHStatBounding_8TeVDown	    ->Write();
      histo_qqH_CMS_MVAqqHStatBounding_8TeVUp	    ->Write();
      histo_qqH_CMS_MVAqqHStatBounding_8TeVDown          ->Write();
      histo_ggH_CMS_MVAggHStatBounding_8TeVUp	    ->Write();
      histo_ggH_CMS_MVAggHStatBounding_8TeVDown          ->Write();
      histo_qqWW_CMS_MVAqqWWStatBounding_8TeVUp     ->Write();
      histo_qqWW_CMS_MVAqqWWStatBounding_8TeVDown        ->Write();
      histo_ggWW_CMS_MVAggWWStatBounding_8TeVUp     ->Write();
      histo_ggWW_CMS_MVAggWWStatBounding_8TeVDown        ->Write();
      histo_VV_CMS_MVAVVStatBounding_8TeVUp	    ->Write();
      histo_VV_CMS_MVAVVStatBounding_8TeVDown	    ->Write();
      histo_Top_CMS_MVATopStatBounding_8TeVUp	    ->Write();
      histo_Top_CMS_MVATopStatBounding_8TeVDown          ->Write();
      histo_Zjets_CMS_MVAZjetsStatBounding_8TeVUp   ->Write();
      histo_Zjets_CMS_MVAZjetsStatBounding_8TeVDown      ->Write();
      histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp   ->Write();
      histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown      ->Write();
      histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp ->Write();
      histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown    ->Write();
      histo_Ztt_CMS_MVAZttStatBounding_8TeVUp	    ->Write();
      histo_Ztt_CMS_MVAZttStatBounding_8TeVDown          ->Write();
    }

    //----------------------------------------------------------------------------
    // W+Jets Bkg Systematics Shapes 
    // MVAWBounding   : jet pt spectrum systematics, mirror the difference between
    //                  the nominal and the "up" systematic fake rates
    // MVAWMCBounding : MC closure test, mirror difference between the nominal and
    //                  the shape from W+Jets MC with data fake rate applied
    //----------------------------------------------------------------------------
    if(useWJetsTemplates == true){
      histo_Wjets_CMS_MVAWBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_Wjets_CMS_MVAWBoundingDown	 ->Rebin(rebinMVAHist);
      double mean,up,diff;
      if(histo_Wjets_CMS_MVAWBoundingUp->GetNbinsX() != histo_Wjets->GetNbinsX()) {printf("Different binning in W!\n"); return;}
      histo_Wjets_CMS_MVAWBoundingUp  ->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingUp  ->GetSumOfWeights());
      for(int i=1; i<=histo_Wjets_CMS_MVAWBoundingUp->GetNbinsX(); i++){
        mean = histo_Wjets                   ->GetBinContent(i);
        up   = histo_Wjets_CMS_MVAWBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
      }
      histo_Wjets_CMS_MVAWBoundingDown->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingDown->GetSumOfWeights());
      histo_Wjets_CMS_MVAWBoundingUp  ->Write();
      histo_Wjets_CMS_MVAWBoundingDown->Write();
    }
    if(useWJetsMCTemplates == true){
      histo_Wjets_CMS_hww_MVAWMCBoundingUp 	->Rebin(rebinMVAHist);
      histo_Wjets_CMS_hww_MVAWMCBoundingDown->Rebin(rebinMVAHist);
      double mean,up,diff;
      if(histo_Wjets_CMS_hww_MVAWMCBoundingUp->GetNbinsX() != histo_Wjets->GetNbinsX()) {printf("Different binning in W!\n"); return;}
      if(histo_Wjets_CMS_hww_MVAWMCBoundingUp  ->GetSumOfWeights() == 0) assert(0);
      histo_Wjets_CMS_hww_MVAWMCBoundingUp  ->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_hww_MVAWMCBoundingUp  ->GetSumOfWeights());
      for(int i=1; i<=histo_Wjets_CMS_hww_MVAWMCBoundingUp->GetNbinsX(); i++){
        mean = histo_Wjets                   ->GetBinContent(i);
        up   = histo_Wjets_CMS_hww_MVAWMCBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Wjets_CMS_hww_MVAWMCBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Wjets_CMS_hww_MVAWMCBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
      }
      histo_Wjets_CMS_hww_MVAWMCBoundingDown->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_hww_MVAWMCBoundingDown->GetSumOfWeights());
      histo_Wjets_CMS_hww_MVAWMCBoundingUp  ->Write();
      histo_Wjets_CMS_hww_MVAWMCBoundingDown->Write();
    }

    //----------------------------------------------------------------------------
    // ggH Systematics Shapes : missing higher order corrections
    //----------------------------------------------------------------------------
    if(useggHTemplates == true){
      histo_ggH_CMS_hww_MVAggHBoundingUp  ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVAggHBoundingDown->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVAggHBoundingUp  ->Scale(histo_ggH->GetSumOfWeights()/histo_ggH_CMS_hww_MVAggHBoundingUp  ->GetSumOfWeights());
      histo_ggH_CMS_hww_MVAggHBoundingDown->Scale(histo_ggH->GetSumOfWeights()/histo_ggH_CMS_hww_MVAggHBoundingDown->GetSumOfWeights());
      histo_ggH_CMS_hww_MVAggHBoundingUp  ->Write();
      histo_ggH_CMS_hww_MVAggHBoundingDown->Write();
    }
 
    //----------------------------------------------------------------------------
    // Systematics Shapes from 
    // - Lepton Efficiencies
    // - Lepton Resolution
    // - MET Resolution
    //----------------------------------------------------------------------------
    if(useExpTemplates == true){
      histo_ttH_CMS_hww_MVALepEffBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ttH_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVALepEffBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVALepEffBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVALepEffBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVALepEffBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_hww_MVALepEffBoundingUp	 ->Rebin(rebinMVAHist);  if(mH <= 200) histo_qqWW_CMS_hww_MVALepEffBoundingUp  ->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVALepEffBoundingUp  ->GetSumOfWeights());
      histo_qqWW_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);  if(mH <= 200) histo_qqWW_CMS_hww_MVALepEffBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVALepEffBoundingDown->GetSumOfWeights());
      histo_ggWW_CMS_hww_MVALepEffBoundingUp	 ->Rebin(rebinMVAHist);  if(mH <= 200) histo_ggWW_CMS_hww_MVALepEffBoundingUp  ->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVALepEffBoundingUp  ->GetSumOfWeights());
      histo_ggWW_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);  if(mH <= 200) histo_ggWW_CMS_hww_MVALepEffBoundingDown->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVALepEffBoundingDown->GetSumOfWeights());
      histo_VV_CMS_hww_MVALepEffBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_VV_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_hww_MVALepEffBoundingUp	 ->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVALepEffBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVALepEffBoundingDown	 ->Rebin(rebinMVAHist);
      if(rebinSmurf > 1){
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ,rebinSmurf);histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ->Scale(0.0);histo_Wgamma_CMS_hww_MVALepEffBoundingUp  ->Add(histSmurfRebin); 
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVALepEffBoundingDown,rebinSmurf);histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Scale(0.0);histo_Wgamma_CMS_hww_MVALepEffBoundingDown->Add(histSmurfRebin);
      }
      histo_ttH_CMS_hww_MVALepEffBoundingUp 	 ->Write();
      histo_ttH_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_ZH_CMS_hww_MVALepEffBoundingUp  	 ->Write();
      histo_ZH_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_WH_CMS_hww_MVALepEffBoundingUp  	 ->Write();
      histo_WH_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_qqH_CMS_hww_MVALepEffBoundingUp 	 ->Write();
      histo_qqH_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_ggH_CMS_hww_MVALepEffBoundingUp 	 ->Write();
      histo_ggH_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_qqWW_CMS_hww_MVALepEffBoundingUp	 ->Write();
      histo_qqWW_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_ggWW_CMS_hww_MVALepEffBoundingUp	 ->Write();
      histo_ggWW_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_VV_CMS_hww_MVALepEffBoundingUp  	 ->Write();
      histo_VV_CMS_hww_MVALepEffBoundingDown	 ->Write();
      histo_Wgamma_CMS_hww_MVALepEffBoundingUp	 ->Write();
      histo_Wgamma_CMS_hww_MVALepEffBoundingDown ->Write();
      histo_Ztt_CMS_hww_MVALepEffBoundingUp 	 ->Write();
      histo_Ztt_CMS_hww_MVALepEffBoundingDown	 ->Write();

      histo_ttH_CMS_hww_MVALepResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ttH_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVALepResBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVALepResBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVALepResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVALepResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_hww_MVALepResBoundingUp	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_qqWW_CMS_hww_MVALepResBoundingUp  ->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVALepResBoundingUp  ->GetSumOfWeights());
      histo_qqWW_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_qqWW_CMS_hww_MVALepResBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVALepResBoundingDown->GetSumOfWeights());
      histo_ggWW_CMS_hww_MVALepResBoundingUp	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_ggWW_CMS_hww_MVALepResBoundingUp  ->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVALepResBoundingUp  ->GetSumOfWeights());
      histo_ggWW_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_ggWW_CMS_hww_MVALepResBoundingDown->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVALepResBoundingDown->GetSumOfWeights());
      histo_VV_CMS_hww_MVALepResBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_VV_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Top_CMS_hww_MVALepResBoundingUp 	 ->Rebin(rebinMVAHist);  histo_Top_CMS_hww_MVALepResBoundingUp    ->Scale(histo_Top->GetSumOfWeights()/ histo_Top_CMS_hww_MVALepResBoundingUp  ->GetSumOfWeights());
      histo_Top_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);  histo_Top_CMS_hww_MVALepResBoundingDown  ->Scale(histo_Top->GetSumOfWeights()/ histo_Top_CMS_hww_MVALepResBoundingDown->GetSumOfWeights());
      histo_Wgamma_CMS_hww_MVALepResBoundingUp	 ->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVALepResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVALepResBoundingDown	 ->Rebin(rebinMVAHist);
      if(rebinSmurf > 1){
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVALepResBoundingUp  ,rebinSmurf);histo_Wgamma_CMS_hww_MVALepResBoundingUp  ->Scale(0.0);histo_Wgamma_CMS_hww_MVALepResBoundingUp  ->Add(histSmurfRebin);
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVALepResBoundingDown,rebinSmurf);histo_Wgamma_CMS_hww_MVALepResBoundingDown->Scale(0.0);histo_Wgamma_CMS_hww_MVALepResBoundingDown->Add(histSmurfRebin);
      }
      histo_ttH_CMS_hww_MVALepResBoundingUp 	 ->Write();
      histo_ttH_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_ZH_CMS_hww_MVALepResBoundingUp  	 ->Write();
      histo_ZH_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_WH_CMS_hww_MVALepResBoundingUp  	 ->Write();
      histo_WH_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_qqH_CMS_hww_MVALepResBoundingUp 	 ->Write();
      histo_qqH_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_ggH_CMS_hww_MVALepResBoundingUp 	 ->Write();
      histo_ggH_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_qqWW_CMS_hww_MVALepResBoundingUp	 ->Write();
      histo_qqWW_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_ggWW_CMS_hww_MVALepResBoundingUp	 ->Write();
      histo_ggWW_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_VV_CMS_hww_MVALepResBoundingUp  	 ->Write();
      histo_VV_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_Top_CMS_hww_MVALepResBoundingUp 	 ->Write();
      histo_Top_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_Wgamma_CMS_hww_MVALepResBoundingUp	 ->Write();
      histo_Wgamma_CMS_hww_MVALepResBoundingDown	 ->Write();
      histo_Ztt_CMS_hww_MVALepResBoundingUp 	 ->Write();
      histo_Ztt_CMS_hww_MVALepResBoundingDown	 ->Write();

      histo_ttH_CMS_hww_MVAMETResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ttH_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVAMETResBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVAMETResBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVAMETResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVAMETResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_hww_MVAMETResBoundingUp	 ->Rebin(rebinMVAHist); 
      histo_qqWW_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist); 
      histo_ggWW_CMS_hww_MVAMETResBoundingUp	 ->Rebin(rebinMVAHist); 
      histo_ggWW_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist); 
      histo_VV_CMS_hww_MVAMETResBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_VV_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Top_CMS_hww_MVAMETResBoundingUp 	 ->Rebin(rebinMVAHist); 
      histo_Top_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist); 
      histo_Wgamma_CMS_hww_MVAMETResBoundingUp	 ->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVAMETResBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVAMETResBoundingDown	 ->Rebin(rebinMVAHist);
      if(histo_qqWW_CMS_hww_MVAMETResBoundingUp->GetNbinsX() != histo_qqWW->GetNbinsX()) {printf("Different binning in METRes!\n"); return;}
      double mean,up,diff;
      for(int i=1; i<=histo_ttH_CMS_hww_MVAMETResBoundingUp->GetNbinsX(); i++){
      
        mean = histo_ttH                        ->GetBinContent(i);
        up   = histo_ttH_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_ttH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_ttH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_ZH                        ->GetBinContent(i);
        up   = histo_ZH_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_ZH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_ZH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_WH                        ->GetBinContent(i);
        up   = histo_WH_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_qqH                        ->GetBinContent(i);
        up   = histo_qqH_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_qqH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_qqH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_ggH                        ->GetBinContent(i);
        up   = histo_ggH_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_ggH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_ggH_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_qqWW                        ->GetBinContent(i);
        up   = histo_qqWW_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_qqWW_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_qqWW_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_ggWW                        ->GetBinContent(i);
        up   = histo_ggWW_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_ggWW_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_ggWW_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_VV                        ->GetBinContent(i);
        up   = histo_VV_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_VV_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_VV_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_Top                        ->GetBinContent(i);
        up   = histo_Top_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Top_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Top_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_Wgamma                        ->GetBinContent(i);
        up   = histo_Wgamma_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Wgamma_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Wgamma_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_Ztt                        ->GetBinContent(i);
        up   = histo_Ztt_CMS_hww_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Ztt_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Ztt_CMS_hww_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
      }

      if(mH <= 200) histo_qqWW_CMS_hww_MVAMETResBoundingUp  ->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVAMETResBoundingUp  ->GetSumOfWeights());
      if(mH <= 200) histo_qqWW_CMS_hww_MVAMETResBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVAMETResBoundingDown->GetSumOfWeights());
      if(mH <= 200) histo_ggWW_CMS_hww_MVAMETResBoundingUp  ->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVAMETResBoundingUp  ->GetSumOfWeights());
      if(mH <= 200) histo_ggWW_CMS_hww_MVAMETResBoundingDown->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVAMETResBoundingDown->GetSumOfWeights());
      histo_Top_CMS_hww_MVAMETResBoundingUp    ->Scale(histo_Top->GetSumOfWeights()/ histo_Top_CMS_hww_MVAMETResBoundingUp  ->GetSumOfWeights());
      histo_Top_CMS_hww_MVAMETResBoundingDown  ->Scale(histo_Top->GetSumOfWeights()/ histo_Top_CMS_hww_MVAMETResBoundingDown->GetSumOfWeights());

      if(rebinSmurf > 1){
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ,rebinSmurf);histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ->Scale(0.0);histo_Wgamma_CMS_hww_MVAMETResBoundingUp  ->Add(histSmurfRebin);
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVAMETResBoundingDown,rebinSmurf);histo_Wgamma_CMS_hww_MVAMETResBoundingDown->Scale(0.0);histo_Wgamma_CMS_hww_MVAMETResBoundingDown->Add(histSmurfRebin);
      }
      histo_ttH_CMS_hww_MVAMETResBoundingUp 	 ->Write();
      histo_ttH_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_ZH_CMS_hww_MVAMETResBoundingUp  	 ->Write();
      histo_ZH_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_WH_CMS_hww_MVAMETResBoundingUp  	 ->Write();
      histo_WH_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_qqH_CMS_hww_MVAMETResBoundingUp 	 ->Write();
      histo_qqH_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_ggH_CMS_hww_MVAMETResBoundingUp 	 ->Write();
      histo_ggH_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_qqWW_CMS_hww_MVAMETResBoundingUp	 ->Write();
      histo_qqWW_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_ggWW_CMS_hww_MVAMETResBoundingUp	 ->Write();
      histo_ggWW_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_VV_CMS_hww_MVAMETResBoundingUp  	 ->Write();
      histo_VV_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_Top_CMS_hww_MVAMETResBoundingUp 	 ->Write();
      histo_Top_CMS_hww_MVAMETResBoundingDown	 ->Write();
      histo_Wgamma_CMS_hww_MVAMETResBoundingUp	 ->Write();
      histo_Wgamma_CMS_hww_MVAMETResBoundingDown ->Write();
      histo_Ztt_CMS_hww_MVAMETResBoundingUp 	 ->Write();
      histo_Ztt_CMS_hww_MVAMETResBoundingDown	 ->Write();
    }

    //----------------------------------------------------------------------------
    // Systematics Shapes from 
    // - Jet Energy Scale 
    //----------------------------------------------------------------------------
    if(useJESTemplates == true){
      histo_ttH_CMS_hww_MVAJESBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ttH_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVAJESBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_ZH_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVAJESBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_WH_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVAJESBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_qqH_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVAJESBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_ggH_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      histo_qqWW_CMS_hww_MVAJESBoundingUp	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_qqWW_CMS_hww_MVAJESBoundingUp  ->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights());
      histo_qqWW_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_qqWW_CMS_hww_MVAJESBoundingDown->Scale(histo_qqWW->GetSumOfWeights()/ histo_qqWW_CMS_hww_MVAJESBoundingDown->GetSumOfWeights());
      histo_ggWW_CMS_hww_MVAJESBoundingUp	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_ggWW_CMS_hww_MVAJESBoundingUp  ->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights());
      histo_ggWW_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist); if(mH <= 200) histo_ggWW_CMS_hww_MVAJESBoundingDown->Scale(histo_ggWW->GetSumOfWeights()/ histo_ggWW_CMS_hww_MVAJESBoundingDown->GetSumOfWeights());
      histo_VV_CMS_hww_MVAJESBoundingUp  	 ->Rebin(rebinMVAHist);
      histo_VV_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      histo_Top_CMS_hww_MVAJESBoundingUp 	 ->Rebin(rebinMVAHist); histo_Top_CMS_hww_MVAJESBoundingUp    ->Scale(histo_Top->GetSumOfWeights()/ histo_Top_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights());
      histo_Top_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist); histo_Top_CMS_hww_MVAJESBoundingDown  ->Scale(histo_Top->GetSumOfWeights()/ histo_Top_CMS_hww_MVAJESBoundingDown->GetSumOfWeights());
      histo_Wgamma_CMS_hww_MVAJESBoundingUp	 ->Rebin(rebinMVAHist);
      histo_Wgamma_CMS_hww_MVAJESBoundingDown->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVAJESBoundingUp 	 ->Rebin(rebinMVAHist);
      histo_Ztt_CMS_hww_MVAJESBoundingDown	 ->Rebin(rebinMVAHist);
      if(histo_ttH_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights() == 0) {histo_ttH_CMS_hww_MVAJESBoundingUp  ->SetBinContent(1,0.000001);}
      if(histo_ttH_CMS_hww_MVAJESBoundingDown->GetSumOfWeights() == 0) {histo_ttH_CMS_hww_MVAJESBoundingDown->SetBinContent(1,0.000001);}
      if(rebinSmurf > 1){
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVAJESBoundingUp  ,rebinSmurf);histo_Wgamma_CMS_hww_MVAJESBoundingUp  ->Scale(0.0);histo_Wgamma_CMS_hww_MVAJESBoundingUp  ->Add(histSmurfRebin);
        histSmurfRebin = SmurfRebin(histo_Wgamma_CMS_hww_MVAJESBoundingDown,rebinSmurf);histo_Wgamma_CMS_hww_MVAJESBoundingDown->Scale(0.0);histo_Wgamma_CMS_hww_MVAJESBoundingDown->Add(histSmurfRebin);
      }
      histo_ttH_CMS_hww_MVAJESBoundingUp 	 ->Write();
      histo_ttH_CMS_hww_MVAJESBoundingDown	 ->Write();
      if(histo_ZH_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights() == 0) {histo_ZH_CMS_hww_MVAJESBoundingUp  ->SetBinContent(1,0.000001);}
      if(histo_ZH_CMS_hww_MVAJESBoundingDown->GetSumOfWeights() == 0) {histo_ZH_CMS_hww_MVAJESBoundingDown->SetBinContent(1,0.000001);}
      histo_ZH_CMS_hww_MVAJESBoundingUp  	 ->Write();
      histo_ZH_CMS_hww_MVAJESBoundingDown	 ->Write();
      if(histo_WH_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights() == 0) {histo_WH_CMS_hww_MVAJESBoundingUp  ->SetBinContent(1,0.000001);}
      if(histo_WH_CMS_hww_MVAJESBoundingDown->GetSumOfWeights() == 0) {histo_WH_CMS_hww_MVAJESBoundingDown->SetBinContent(1,0.000001);}
      histo_WH_CMS_hww_MVAJESBoundingUp  	 ->Write();
      histo_WH_CMS_hww_MVAJESBoundingDown	 ->Write();
      histo_qqH_CMS_hww_MVAJESBoundingUp 	 ->Write();
      histo_qqH_CMS_hww_MVAJESBoundingDown	 ->Write();
      histo_ggH_CMS_hww_MVAJESBoundingUp 	 ->Write();
      histo_ggH_CMS_hww_MVAJESBoundingDown	 ->Write();
      histo_qqWW_CMS_hww_MVAJESBoundingUp	 ->Write();
      histo_qqWW_CMS_hww_MVAJESBoundingDown	 ->Write();
      histo_ggWW_CMS_hww_MVAJESBoundingUp	 ->Write();
      histo_ggWW_CMS_hww_MVAJESBoundingDown	 ->Write();
      histo_VV_CMS_hww_MVAJESBoundingUp  	 ->Write();
      histo_VV_CMS_hww_MVAJESBoundingDown	 ->Write();
      histo_Top_CMS_hww_MVAJESBoundingUp 	 ->Write();
      histo_Top_CMS_hww_MVAJESBoundingDown	 ->Write();
      if(histo_Wgamma_CMS_hww_MVAJESBoundingUp  ->GetSumOfWeights() == 0) {histo_Wgamma_CMS_hww_MVAJESBoundingUp  ->SetBinContent(1,0.000001);}
      if(histo_Wgamma_CMS_hww_MVAJESBoundingDown->GetSumOfWeights() == 0) {histo_Wgamma_CMS_hww_MVAJESBoundingDown->SetBinContent(1,0.000001);}
      histo_Wgamma_CMS_hww_MVAJESBoundingUp	 ->Write();
      histo_Wgamma_CMS_hww_MVAJESBoundingDown->Write();
      histo_Ztt_CMS_hww_MVAJESBoundingUp 	 ->Write();
      histo_Ztt_CMS_hww_MVAJESBoundingDown	 ->Write();
    }
    outFileLimits->Close();

    //----------------------------------------------------------------------------
    // 
    // Systematics For Yields
    //
    //----------------------------------------------------------------------------
    double theoryUncXS_HighMH = 1.0;
    if(mH > 200) theoryUncXS_HighMH = 1.0+1.5*(mH/1000.0)*(mH/1000.0)*(mH/1000.0);
    double wwXS_E_jet_extrap = 1.060;
    double jeteff_E 	     = 1.02;
    double topXS_E  	     = TopBkgScaleFactorKappa(nJetsType); if(nJetsType==2) topXS_E = TopVBFBkgScaleFactorKappa(mH);
    double wwXS_E_MVA        = WWBkgScaleFactorKappaMVA     (TMath::Max((int)mH,115),TMath::Min((int)nJetsType,1)); if(mH > 200) wwXS_E_MVA = 1.000;
    double wwXS_E_Cut        = WWBkgScaleFactorKappaCutBased(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,1)); if(mH > 200) wwXS_E_Cut = 1.000;
    char theWWThString[20]; sprintf(theWWThString,"CMS_hww_%1dj_WW_8TeV",nJetsType); if(mH > 200) sprintf(theWWThString,"CMS_hww_WW");

    double XS_QCDscale_WW[3] = {1.0, 1.0, 1.0};
    if(mH > 200) {XS_QCDscale_WW[0] = 1.042; XS_QCDscale_WW[1] = 0.978; XS_QCDscale_WW[2] = 1.000;}

    double pdf_ggH = PDFgHHSystematics(mH);

    double XS_QCDscale_ggH[3];
    double UEPS  = HiggsSignalPSUESystematics(mH,nJetsType);
    XS_QCDscale_ggH[0] = HiggsSignalQCDScaleKappa("QCDscale_ggH",mH,nJetsType);
    XS_QCDscale_ggH[1] = HiggsSignalQCDScaleKappa("QCDscale_ggH1in",mH,nJetsType);
    XS_QCDscale_ggH[2] = HiggsSignalQCDScaleKappa("QCDscale_ggH2in",mH,nJetsType);
    
    double XS_PDF_VH = 1.05; 
    double XS_QCDscale_qqH = 1.01; double XS_QCDscale_VH = 1.02; 
    if(isFermioPhobic == true) {XS_QCDscale_qqH += 0.05; XS_QCDscale_VH += 0.05; XS_PDF_VH += 0.00;}
    
    double gamma_Hff = 1.0; double gamma_HVV = 1.0; double gamma_Hgluglu = 1.0;
    if(isSM4 == true) {
      gamma_Hff     = HiggsSM4Systematics_HiggsBRErr_Hff    (mH);
      gamma_HVV     = HiggsSM4Systematics_HiggsBRErr_HVV    (mH);
      gamma_Hgluglu = HiggsSM4Systematics_HiggsBRErr_Hgluglu(mH);
    }

    //                                eff_m,eff_e,scale_m,scale_e,hww_met_resolution
    double theExpUncertainties[5] = {1.030, 1.040, 1.015, 1.020, 1.020};
    if(mH <= 200) for(int i=0; i<5; i++) theExpUncertainties[i] = 1.0;

    if     (nJetsType == 1) {
      jeteff_E  	= 1.05;
      if(mH > 200) {XS_QCDscale_WW[0] = 1.000; XS_QCDscale_WW[1] = 1.076; XS_QCDscale_WW[2] = 0.914;}
    }
    else if(nJetsType == 2) {
      jeteff_E  	= 1.10;
      if(mH > 200) {XS_QCDscale_WW[0] = 1.000; XS_QCDscale_WW[1] = 1.000; XS_QCDscale_WW[2] = 1.420;}
    }
    double lumiErr = 1.000; if(mH > 200) lumiErr = 1.050;

    for(int i=0; i<8; i++) if(nBgdAccDecays[i] < 0) nBgdAccDecays[i] = 0.0;
    for(int i=0; i<8; i++) if(nBgdCutDecays[i] < 0) nBgdCutDecays[i] = 0.0;
    for(int i=0; i<8; i++) if(nBgdMVADecays[i] < 0) nBgdMVADecays[i] = 0.0;
    for(int i=0; i<6; i++) if(nSigAcc[i] <= 0) nSigAcc[i] = 0.000;
    for(int i=0; i<6; i++) if(nSigCut[i] <= 0) nSigCut[i] = 0.000;
    for(int i=0; i<6; i++) if(nSigMVA[i] <= 0) nSigMVA[i] = 0.000;
    double yieldE[13],yield[13];
    int nData;
    int nTotalBins = 1; // histo_VH->GetNbinsX();

    //----------------------------------------------------------------------------
    // Yields for MVA Shape Analysis
    //----------------------------------------------------------------------------
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
	if(signalInjection == true){
	  nData = (int)histo_Data->GetSumOfWeights();
	}
	else {
	  if(histo_Data->GetSumOfWeights() != nDatAcc) {
	    printf("histo_Data != nDatAcc\n");
	    assert(0);
	  }
	}
      }
      if(nTotalBins != 1){
        sprintf(outputLimitsShape,"output/histo_limits_%s_%dj_chan%d_mh%d_shape_bin%d_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH,i);
	if(category == 1) sprintf(outputLimitsShape,"output/histo_limits_%s_%dj_chan%d_mh%d_shape_bin%d_lt_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH,i);
      }
      else {
        sprintf(outputLimitsShape,"output/histo_limits_%s_%dj_chan%d_mh%d_shape_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH);
	if(category == 1) sprintf(outputLimitsShape,"output/histo_limits_%s_%dj_chan%d_mh%d_shape_lt_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH);
      }
      char theZttString[20];
      if(histo_Ztt->GetSumOfWeights() > 0) sprintf(theZttString,"1.000");
      else                                 sprintf(theZttString,"  -  ");
      char theZHString[20];
      if(histo_ZH->GetSumOfWeights() > 0) sprintf(theZHString,"1.000");
      else                                sprintf(theZHString,"  -  ");
      char theWHString[20];
      if(histo_WH->GetSumOfWeights() > 0) sprintf(theWHString,"1.000");
      else                                sprintf(theWHString,"  -  ");
      char theqqHString[20];
      if(histo_qqH->GetSumOfWeights() > 0) sprintf(theqqHString,"1.000");
      else                                 sprintf(theqqHString,"  -  ");
      char theggHString[20];
      if(histo_ggH->GetSumOfWeights() > 0) sprintf(theggHString,"1.000");
      else                                 sprintf(theggHString,"  -  ");
      char theWgammaString[20];
      if(histo_Wgamma->GetSumOfWeights() > 0) sprintf(theWgammaString,"1.000");
      else                                    sprintf(theWgammaString,"  -  ");

      //----------------------------------------------------------------------------
      // Produce output cards for MVA Shape analysis
      //----------------------------------------------------------------------------
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",nData);
      if(nTotalBins == 1){
        if(useZjetsTemplates == true || useWWTemplates      == true || useStatTemplates  == true || useExpTemplates   == true ||
	   useWJetsTemplates == true || useWJetsMCTemplates == true || useggHTemplates   == true || useTopTemplates   == true)
          newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
        else
          newcardShape << Form("shapes *   *   %s  histo_$PROCESS\n",outputLimits);
        newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
      }
      newcardShape << Form("bin j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s\n",nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName);
      newcardShape << Form("process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt\n");
      newcardShape << Form("process -3 -2 -1 0 1 2 3 4 5 6 7 8\n");
      newcardShape << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",yield[1],yield[2],yield[3],yield[4],yield[5],yield[6],yield[7],yield[8],yield[9],TMath::Max((double)yield[10],0.0),yield[11],yield[12]);
      newcardShape << Form("lumi_8TeV                        lnN 1.050 1.050 1.050 1.050 %5.3f %5.3f 1.050   -     -     -   1.050 1.050\n",lumiErr,lumiErr);				 
      if(useExpTemplates == true){
      newcardShape << Form("CMS_hww_MVALepEffBounding          shape   %s   %s   %s   %s   1.000 1.000 1.000   -     -     -   %s %s\n",theZHString,theWHString,theqqHString,theggHString,theWgammaString,theZttString);			   
      newcardShape << Form("CMS_hww_MVALepResBounding          shape   %s   %s   %s   %s   1.000 1.000 1.000 1.000   -     -   %s %s\n",theZHString,theWHString,theqqHString,theggHString,theWgammaString,theZttString);			   
      newcardShape << Form("CMS_hww_MVAMETResBounding          shape   %s   %s   %s   %s   1.000 1.000 1.000 1.000   -     -   %s %s\n",theZHString,theWHString,theqqHString,theggHString,theWgammaString,theZttString);			   
      }
      else {
      newcardShape << Form("CMS_eff_m                        lnN 1.030 1.030 1.030 1.030 %5.3f %5.3f 1.030   -     -	 -   1.030 1.030\n",theExpUncertainties[0],theExpUncertainties[0]);
      newcardShape << Form("CMS_eff_e                        lnN 1.040 1.040 1.040 1.040 %5.3f %5.3f 1.040   -     -	 -   1.040 1.040\n",theExpUncertainties[1],theExpUncertainties[1]);			      
      newcardShape << Form("CMS_scale_m                      lnN 1.015 1.015 1.015 1.015 %5.3f %5.3f 1.015   -     -	 -   1.015 1.015\n",theExpUncertainties[2],theExpUncertainties[2]);			      
      newcardShape << Form("CMS_scale_e                      lnN 1.020 1.020 1.020 1.020 %5.3f %5.3f 1.020   -     -	 -   1.020 1.020\n",theExpUncertainties[3],theExpUncertainties[3]);
      newcardShape << Form("CMS_hww_met_resolution           lnN 1.020 1.020 1.020 1.020 %5.3f %5.3f 1.020   -     -	 -   1.020 1.020\n",theExpUncertainties[4],theExpUncertainties[4]);
      }
      if(useJESTemplates == true){
      newcardShape << Form("CMS_hww_MVAJESBounding             shape   %s   %s   %s   %s   1.000 1.000 1.000 1.000   -     -   %s %s\n",theZHString,theWHString,theqqHString,theggHString,theWgammaString,theZttString);			   
      }
      else {
      newcardShape << Form("CMS_scale_j                      lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -     -   %5.3f %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);	     
      }
      newcardShape << Form("FakeRate              	     lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  \n");
      if(useWJetsTemplates == true){
        newcardShape << Form("CMS_hww_MVAWBounding          shape  -     -     -     -     -     -     -     -     -   1.000   -     -  \n");
      }
      if(useWJetsMCTemplates == true){
        newcardShape << Form("CMS_hww_MVAWMCBounding        shape  -     -     -     -     -     -     -     -     -   1.000   -     -  \n");
      }
      if(useggHTemplates == true && histo_ggH->GetSumOfWeights() > 0.0){
        newcardShape << Form("CMS_hww_MVAggHBounding            shape  -     -     -   1.000   -     -     -     -     -     -    -     -  \n");
      }
      newcardShape << Form("UEPS 		             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",UEPS);
      newcardShape << Form("theoryUncXS_HighMH               lnN %5.3f %5.3f %5.3f %5.3f   -     -     -     -	   -     -     -     -  \n",theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH);
      newcardShape << Form("pdf_gg		             lnN   -     -     -   %5.3f   -   1.040   -     -     -     -     -     -  \n",pdf_ggH);
      newcardShape << Form("pdf_qqbar                        lnN %5.3f %5.3f %5.3f   -   1.040   -   1.040   -     -     -   1.040 1.040\n",XS_PDF_VH,XS_PDF_VH,XS_PDF_VH);
      newcardShape << Form("QCDscale_ggH	             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[0]);  
      newcardShape << Form("QCDscale_ggH1in	             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[1]);  
      newcardShape << Form("QCDscale_ggH2in	             lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[2]);  
      newcardShape << Form("QCDscale_qqH	             lnN   -     -   %5.3f   -     -     -     -     -     -     -     -     -  \n",XS_QCDscale_qqH);
      newcardShape << Form("QCDscale_VH 	             lnN %5.3f %5.3f   -     -     -     -     -     -     -     -     -     -  \n",XS_QCDscale_VH,XS_QCDscale_VH);			 
      newcardShape << Form("QCDscale_WW	                     lnN   -     -     -     -   %5.3f   -     -     -     -	 -     -     -  \n",XS_QCDscale_WW[0]);  
      newcardShape << Form("QCDscale_WW1in		     lnN   -	 -     -     -	 %5.3f	 -     -     -	   -	 -     -     -  \n",XS_QCDscale_WW[1]);  
      newcardShape << Form("QCDscale_WW2in		     lnN   -	 -     -     -	 %5.3f	 -     -     -	   -	 -     -     -  \n",XS_QCDscale_WW[2]);  
      newcardShape << Form("QCDscale_VV           	     lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_Vgamma                  lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",1.30);
      newcardShape << Form("QCDscale_ggVV	             lnN   -     -     -     -     -   1.300   -     -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_WW_EXTRAP               lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",wwXS_E_jet_extrap);
      newcardShape << Form("QCDscale_ggH_ACCEPT              lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_qqH_ACCEPT              lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("QCDscale_VH_ACCEPT               lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");
      newcardShape << Form("CMS_hww_%1dj_ttbar_8TeV          lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",nJetsType,topXS_E); 	
      newcardShape << Form("CMS_hww%s_%1dj_Z_8TeV            lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,ZXS_E[0]+1.0);			
      newcardShape << Form("%s                               lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",theWWThString,wwXS_E_MVA,wwXS_E_MVA);				
      newcardShape << Form("CMS_hww_Ztt           	     lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",ZttScaleFactorKappa());
      if(useZjetsTemplates == true){
        newcardShape << Form("CMS_hww%s_%1dj_MVAZBounding           shape   -     -      -    -      -     -	 -     -    2.0    -	 -     -   \n",finalStateName,nJetsType);			  
      }
      if(useTopTemplates == true){
        newcardShape << Form("CMS_hww_MVATopBounding                shape   -     -      -    -      -    -	 -    1.0    -     -	 -     -   \n");			  
      }
      if(useWWTemplates == true){
        newcardShape << Form("CMS_hww_MVAWWBounding                 shape   -     -      -    -     1.0    -	 -     -     -     -	 -     -   \n");			 
        newcardShape << Form("CMS_hww_MVAWWNLOBounding              shape   -     -      -    -     1.0    -	 -     -     -     -	 -     -   \n");			 
      }
      if(useStatTemplates == true){
	if(histo_ZH->GetSumOfWeights() > 0)
      	newcardShape << Form("CMS_hww%s_%1dj_MVAZHStatBounding_8TeV      shape  %s     -      -    -      -	   -	 -     -     -	   -	 -     -  \n",finalStateName,nJetsType,theZHString);
	if(histo_WH->GetSumOfWeights() > 0)
      	newcardShape << Form("CMS_hww%s_%1dj_MVAWHStatBounding_8TeV      shape   -     %s     -    -      -	   -	 -     -     -	   -	 -     -  \n",finalStateName,nJetsType,theWHString);
	if(histo_qqH->GetSumOfWeights() > 0)
      	newcardShape << Form("CMS_hww%s_%1dj_MVAqqHStatBounding_8TeV     shape   -     -    1.0    -      -     -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType);
	if(histo_ggH->GetSumOfWeights() > 0)
      	newcardShape << Form("CMS_hww%s_%1dj_MVAggHStatBounding_8TeV     shape   -     -	 -    1.0    -     -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType);
      	newcardShape << Form("CMS_hww%s_%1dj_MVAqqWWStatBounding_8TeV    shape   -     -	 -     -    1.0    -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType);
	if(histo_ggWW->GetSumOfWeights() > 0)
      	newcardShape << Form("CMS_hww%s_%1dj_MVAggWWStatBounding_8TeV    shape   -     -	 -     -     -    1.0    -     -     -     -	 -     -  \n",finalStateName,nJetsType);
      	newcardShape << Form("CMS_hww%s_%1dj_MVAVVStatBounding_8TeV      shape   -     -	 -     -     -     -    1.0    -     -     -	 -     -  \n",finalStateName,nJetsType);
      	newcardShape << Form("CMS_hww%s_%1dj_MVATopStatBounding_8TeV     shape   -     -	 -     -     -     -	 -    1.0    -     -	 -     -  \n",finalStateName,nJetsType);
      	newcardShape << Form("CMS_hww%s_%1dj_MVAZjetsStatBounding_8TeV   shape   -     -	 -     -     -     -	 -     -    1.0    -	 -     -  \n",finalStateName,nJetsType);
      	newcardShape << Form("CMS_hww%s_%1dj_MVAWjetsStatBounding_8TeV   shape   -     -	 -     -     -     -	 -     -     -    1.0    -     -  \n",finalStateName,nJetsType);
	if(histo_Wgamma->GetSumOfWeights() > 0)
      	newcardShape << Form("CMS_hww%s_%1dj_MVAWgammaStatBounding_8TeV  shape   -     -	 -     -     -     -	 -     -     -     -    %s    -  \n",finalStateName,nJetsType,theWgammaString);
	if(histo_Ztt->GetSumOfWeights() > 0)
	newcardShape << Form("CMS_hww%s_%1dj_MVAZttStatBounding_8TeV     shape   -     -	 -     -     -     -	 -     -     -     -	 -    %s  \n",finalStateName,nJetsType,theZttString);
      }
      else {
      	newcardShape << Form("CMS_hww%s_stat_%1dj_ZH_bin%d_8TeV     lnN %5.3f   -	 -     -     -     -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[1]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_WH_bin%d_8TeV     lnN   -   %5.3f   -     -     -     -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[2]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_qqH_bin%d_8TeV    lnN   -     -   %5.3f   -     -     -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[3]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_ggH_bin%d_8TeV    lnN   -     -	 -   %5.3f   -     -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[4]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_WW_bin%d_8TeV     lnN   -     -	 -     -   %5.3f   -	 -     -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[5]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_ggWW_bin%d_8TeV   lnN   -     -	 -     -     -   %5.3f   -     -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[6]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_VV_bin%d_8TeV     lnN   -     -	 -     -     -     -   %5.3f   -     -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[7]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_ttbar_bin%d_8TeV  lnN   -     -	 -     -     -     -	 -   %5.3f   -     -	 -     -  \n",finalStateName,nJetsType,i,yieldE[8]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_Z_bin%d_8TeV      lnN   -     -	 -     -     -     -	 -     -   %5.3f   -	 -     -  \n",finalStateName,nJetsType,i,yieldE[9]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_Wjets_bin%d_8TeV  lnN   -     -	 -     -     -     -	 -     -     -   %5.3f   -     -  \n",finalStateName,nJetsType,i,yieldE[10]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_Wgamma_bin%d_8TeV lnN   -     -	 -     -     -     -	 -     -     -     -   %5.3f   -  \n",finalStateName,nJetsType,i,yieldE[11]+1.0);
      	newcardShape << Form("CMS_hww%s_stat_%1dj_Ztt_bin%d_8TeV    lnN   -     -	 -     -     -     -	 -     -     -     -	 -   %5.3f\n",finalStateName,nJetsType,i,yieldE[12]+1.0);
      }
      if(isSM4 == true){
        newcardShape << Form("gamma_Hff                           lnN   -     -	    -   %5.3f   -     -	    -     -     -     -	    -     -\n",gamma_Hff);
        newcardShape << Form("gamma_HVV                           lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -\n",gamma_HVV);
        newcardShape << Form("gamma_Hgluglu                       lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -\n",gamma_Hgluglu);
      }
      newcardShape.close();
    }
    if(signalInjection == true){
      nDatCut = (nDatCut+
                 nBgdCutDecays[0]+nBgdCutDecays[1]+nBgdCutDecays[2]+nBgdCutDecays[3]+
		 nBgdCutDecays[4]+nBgdCutDecays[5]+nBgdCutDecays[6]+nBgdCutDecays[7])+0.5;
    }

    //----------------------------------------------------------------------------
    // Produce output cards for cut-based analysis
    //----------------------------------------------------------------------------
    char outputLimitsCut[200];
    sprintf(outputLimitsCut,"output/histo_limits_%s_%dj_chan%d_mh%d_cut_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH); 
    if(category == 1) sprintf(outputLimitsCut,"output/histo_limits_%s_%dj_chan%d_mh%d_cut_lt_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH);  
    ofstream newcardCut;
    newcardCut.open(outputLimitsCut);
    newcardCut << Form("imax 1 number of channels\n");
    newcardCut << Form("jmax * number of background\n");
    newcardCut << Form("kmax * number of nuisance parameters\n");
    newcardCut << Form("Observation %d\n",(int)nDatCut);
    newcardCut << Form("bin j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s\n",nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName);
    newcardCut << Form("process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt\n");
    newcardCut << Form("process -3 -2 -1 0 1 2 3 4 5 6 7 8\n");
    newcardCut << Form("rate  %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigCut[2],nSigCut[3],nSigCut[4],nSigCut[5],nBgdCutDecays[0],nBgdCutDecays[1],nBgdCutDecays[2],nBgdCutDecays[3],nBgdCutDecays[4],TMath::Max((double)nBgdCutDecays[5],0.0),nBgdCutDecays[6],nBgdCutDecays[7]);
    newcardCut << Form("lumi_8TeV                  lnN 1.050 1.050 1.050 1.050 %5.3f %5.3f 1.050   -     -     -   1.050 1.050\n",lumiErr,lumiErr);				 
    newcardCut << Form("CMS_eff_m		   lnN 1.030 1.030 1.030 1.030 %5.3f %5.3f 1.030   -	 -     -   1.030 1.030\n",theExpUncertainties[0],theExpUncertainties[0]);
    newcardCut << Form("CMS_eff_e		   lnN 1.040 1.040 1.040 1.040 %5.3f %5.3f 1.040   -	 -     -   1.040 1.040\n",theExpUncertainties[1],theExpUncertainties[1]);			    
    newcardCut << Form("CMS_scale_m		   lnN 1.015 1.015 1.015 1.015 %5.3f %5.3f 1.015   -	 -     -   1.015 1.015\n",theExpUncertainties[2],theExpUncertainties[2]);			    
    newcardCut << Form("CMS_scale_e		   lnN 1.020 1.020 1.020 1.020 %5.3f %5.3f 1.020   -	 -     -   1.020 1.020\n",theExpUncertainties[3],theExpUncertainties[3]);
    newcardCut << Form("CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 %5.3f %5.3f 1.020   -	 -     -   1.020 1.020\n",theExpUncertainties[4],theExpUncertainties[4]);
    newcardCut << Form("CMS_scale_j           	   lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -     -   %5.3f %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);		    
    newcardCut << Form("FakeRate              	   lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  \n");
    newcardCut << Form("UEPS 		           lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",UEPS);
    newcardCut << Form("theoryUncXS_HighMH         lnN %5.3f %5.3f %5.3f %5.3f   -     -     -     -	 -     -     -     -  \n",theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH);
    newcardCut << Form("pdf_gg                	   lnN   -     -     -   %5.3f   -   1.040   -     -     -     -     -     -  \n",pdf_ggH);
    newcardCut << Form("pdf_qqbar             	   lnN %5.3f %5.3f %5.3f   -   1.040   -   1.040   -     -     -   1.040 1.040\n",XS_PDF_VH,XS_PDF_VH,XS_PDF_VH);
    newcardCut << Form("QCDscale_ggH          	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[0]);  
    newcardCut << Form("QCDscale_ggH1in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[1]);  
    newcardCut << Form("QCDscale_ggH2in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[2]);  
    newcardCut << Form("QCDscale_qqH          	   lnN   -     -   %5.3f   -	 -     -     -     -	 -     -     -     -  \n",XS_QCDscale_qqH);
    newcardCut << Form("QCDscale_VH           	   lnN %5.3f %5.3f   -     -	 -     -     -     -	 -     -     -     -  \n",XS_QCDscale_VH,XS_QCDscale_VH);		   
    newcardCut << Form("QCDscale_WW		   lnN   -     -     -     -   %5.3f   -     -     -	 -     -     -     -  \n",XS_QCDscale_WW[0]);  
    newcardCut << Form("QCDscale_WW1in  	   lnN   -     -     -     -   %5.3f   -     -     -	 -     -     -     -  \n",XS_QCDscale_WW[1]);  
    newcardCut << Form("QCDscale_WW2in  	   lnN   -     -     -     -   %5.3f   -     -     -	 -     -     -     -  \n",XS_QCDscale_WW[2]);  
    newcardCut << Form("QCDscale_VV           	   lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_Vgamma            lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",1.30);
    newcardCut << Form("QCDscale_ggVV         	   lnN   -     -     -     -     -   1.300   -     -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_WW_EXTRAP         lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",wwXS_E_jet_extrap);
    newcardCut << Form("QCDscale_ggH_ACCEPT   	   lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_qqH_ACCEPT   	   lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("QCDscale_VH_ACCEPT    	   lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");
    newcardCut << Form("CMS_hww_%1dj_ttbar_8TeV	   lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",nJetsType,topXS_E);	 
    newcardCut << Form("CMS_hww%s_%1dj_Z_8TeV      lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,ZXS_E[1]+1.0);		 
    newcardCut << Form("%s                         lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",theWWThString,wwXS_E_Cut,wwXS_E_Cut);		      
    newcardCut << Form("CMS_hww_Ztt                lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",ZttScaleFactorKappa());
    newcardCut << Form("CMS_hww%s_stat_%1dj_ZH_8TeV	lnN %5.3f   -     -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[2]/TMath::Max((double)nSigCut[2],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_WH_8TeV	lnN   -   %5.3f   -     -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[3]/TMath::Max((double)nSigCut[3],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_qqH_8TeV    lnN   -     -   %5.3f   -     -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[4]/TMath::Max((double)nSigCut[4],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_ggH_8TeV	lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nSigECut[5]/TMath::Max((double)nSigCut[5],0.00001)+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_WW_8TeV	lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[0]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_ggWW_8TeV   lnN   -     -     -     -     -   %5.3f   -     -     -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[1]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_VV_8TeV	lnN   -     -     -     -     -     -   %5.3f   -     -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[2]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_ttbar_8TeV  lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[3]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Z_8TeV      lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,nBgdECutDecays[4]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Wjets_8TeV  lnN   -     -     -     -     -     -     -     -     -   %5.3f   -     -  \n",finalStateName,nJetsType,nBgdECutDecays[5]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Wgamma_8TeV lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",finalStateName,nJetsType,nBgdECutDecays[6]+1.0);
    newcardCut << Form("CMS_hww%s_stat_%1dj_Ztt_8TeV    lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",finalStateName,nJetsType,nBgdECutDecays[7]+1.0);
    if(isSM4 == true){
      newcardCut << Form("gamma_Hff                  lnN   -	 -     -   %5.3f   -	 -     -     -     -	 -     -     -  \n",gamma_Hff);
      newcardCut << Form("gamma_HVV                  lnN   -	 -     -   %5.3f   -	 -     -     -     -	 -     -     -  \n",gamma_HVV);
      newcardCut << Form("gamma_Hgluglu              lnN   -	 -     -   %5.3f   -	 -     -     -     -	 -     -     -  \n",gamma_Hgluglu);
    }
    newcardCut.close();

    if(signalInjection == true){
      nDatMVA = (nDatMVA+
                 nBgdMVADecays[0]+nBgdMVADecays[1]+nBgdMVADecays[2]+nBgdMVADecays[3]+
		 nBgdMVADecays[4]+nBgdMVADecays[5]+nBgdMVADecays[6]+nBgdMVADecays[7])+0.5;
    }

    //----------------------------------------------------------------------------
    // Produce output cards for MVA cut analysis
    //----------------------------------------------------------------------------
    char outputLimitsMVA[200];
    sprintf(outputLimitsMVA,"output/histo_limits_%s_%dj_chan%d_mh%d_mva_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH);
    if(category == 1)    sprintf(outputLimitsMVA,"output/histo_limits_%s_%dj_chan%d_mh%d_mva_lt_8TeV.txt",outTag.Data(),nJetsType,wwDecay,mH);
    ofstream newcardMVA;
    newcardMVA.open(outputLimitsMVA);
    newcardMVA << Form("imax 1 number of channels\n");
    newcardMVA << Form("jmax * number of background\n");
    newcardMVA << Form("kmax * number of nuisance parameters\n");
    newcardMVA << Form("Observation %d\n",(int)nDatMVA);
    newcardMVA << Form("bin j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s j%1d%s\n",nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName,nJetsType,finalStateName);
    newcardMVA << Form("process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt\n");
    newcardMVA << Form("process -3 -2 -1 0 1 2 3 4 5 6 7 8\n");
    newcardMVA << Form("rate  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigMVA[2],nSigMVA[3],nSigMVA[4],nSigMVA[5],nBgdMVADecays[0],nBgdMVADecays[1],nBgdMVADecays[2],nBgdMVADecays[3],nBgdMVADecays[4],TMath::Max((double)nBgdMVADecays[5],0.0),nBgdMVADecays[6],nBgdMVADecays[7]);
    newcardMVA << Form("lumi_8TeV             	   lnN 1.050 1.050 1.050 1.050   -     -   1.050   -     -     -   1.050 1.050\n");
    newcardMVA << Form("CMS_eff_m             	   lnN 1.030 1.030 1.030 1.030 1.030 1.030 1.030   -     -     -   1.030 1.030\n");			    
    newcardMVA << Form("CMS_eff_e             	   lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040   -     -     -   1.040 1.040\n");			    
    newcardMVA << Form("CMS_scale_m           	   lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.015   -     -     -   1.015 1.015\n");			    
    newcardMVA << Form("CMS_scale_e           	   lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");			    
    newcardMVA << Form("CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020\n");			    
    newcardMVA << Form("CMS_scale_j           	   lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -     -   %5.3f %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);	      
    newcardMVA << Form("FakeRate              	   lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  \n");
    newcardMVA << Form("UEPS 		           lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",UEPS);
    newcardMVA << Form("theoryUncXS_HighMH         lnN %5.3f %5.3f %5.3f %5.3f   -     -     -     -	 -     -     -     -  \n",theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH);
    newcardMVA << Form("pdf_gg                	   lnN   -     -     -   %5.3f   -   1.040   -     -     -     -     -     -  \n",pdf_ggH);
    newcardMVA << Form("pdf_qqbar             	   lnN %5.3f %5.3f %5.3f   -   1.040   -   1.040   -     -     -   1.040 1.040\n",XS_PDF_VH,XS_PDF_VH,XS_PDF_VH);
    newcardMVA << Form("QCDscale_ggH          	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[0]);  
    newcardMVA << Form("QCDscale_ggH1in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[1]);  
    newcardMVA << Form("QCDscale_ggH2in       	   lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",XS_QCDscale_ggH[2]);  
    newcardMVA << Form("QCDscale_qqH          	   lnN   -     -   %5.3f   -	 -     -     -     -	 -     -     -     -  \n",XS_QCDscale_qqH);
    newcardMVA << Form("QCDscale_VH           	   lnN %5.3f %5.3f   -     -	 -     -     -     -	 -     -     -     -  \n",XS_QCDscale_VH,XS_QCDscale_VH);		    
    newcardMVA << Form("QCDscale_WW		   lnN   -     -     -     -   %5.3f   -     -     -	 -     -     -     -  \n",XS_QCDscale_WW[0]);  
    newcardMVA << Form("QCDscale_WW1in  	   lnN   -     -     -     -   %5.3f   -     -     -	 -     -     -     -  \n",XS_QCDscale_WW[1]);  
    newcardMVA << Form("QCDscale_WW2in  	   lnN   -     -     -     -   %5.3f   -     -     -	 -     -     -     -  \n",XS_QCDscale_WW[2]);  
    newcardMVA << Form("QCDscale_VV           	   lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_Vgamma            lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",1.30);
    newcardMVA << Form("QCDscale_ggVV         	   lnN   -     -     -     -     -   1.300   -     -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_WW_EXTRAP         lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",wwXS_E_jet_extrap);
    newcardMVA << Form("QCDscale_ggH_ACCEPT   	   lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_qqH_ACCEPT   	   lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("QCDscale_VH_ACCEPT    	   lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  \n");
    newcardMVA << Form("CMS_hww_%1dj_ttbar	   lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",nJetsType,topXS_E);	 
    newcardMVA << Form("CMS_hww%s_%1dj_Z           lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",finalStateName,nJetsType,ZXS_E[2]+1.0);		 
    newcardMVA << Form("%s                         lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",theWWThString,wwXS_E_Cut,wwXS_E_Cut);		      
    newcardMVA << Form("CMS_hww_Ztt                lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",ZttScaleFactorKappa());
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
    if(isSM4 == true){
      newcardMVA << Form("gamma_Hff			     lnN   -	 -     -   %5.3f   -	 -     -     -     -	 -     -     -  \n",gamma_Hff);
      newcardMVA << Form("gamma_HVV			     lnN   -	 -     -   %5.3f   -	 -     -     -     -	 -     -     -  \n",gamma_HVV);
      newcardMVA << Form("gamma_Hgluglu                      lnN   -	 -     -   %5.3f   -	 -     -     -     -	 -     -     -  \n",gamma_Hgluglu);
    }
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
