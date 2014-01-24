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
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/HWWCuts.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double DeltaPhi(double phi1, double phi2);

int    verboseLevel =   0;
const double sigmaB = 0.35;

//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void SummarizeFakeLeptonBkgPrediction
(
//  UInt_t  nJetsType   	 = 0,
 UInt_t  mH      	 = 0,
 TString wjetsMCInputFile = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/wjets.root",
 TString bgdInputFile     = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/backgroundA_skim6.root",
 int period               = 3
 )
{

  bool wwPresel = false;
  if(mH == 0) {wwPresel = true; }

  TString wjetsMCFile1 = wjetsMCInputFile;
  TString bgdFile1 = bgdInputFile;

  unsigned int patternTopTag = SmurfTree::TopTag;

  float dilmass_cut = DileptonMassPreselectionCut(mH);
  if(wwPresel == true) dilmass_cut = 99999.;


  cout << "Using dilepton mass < " << dilmass_cut << endl;

  TChain *chwjetsMC = new TChain("tree");
  chwjetsMC->Add(wjetsMCFile1);

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TTree *wjetsMC = (TTree*) chwjetsMC;
  TTree *background = (TTree*) chbackground;

  double scaleFactorLum = 1.0;
  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  if	 (period == 0){ // Full2011-Fall11-V9-3500ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    scaleFactorLum     =3.553;minRun =      0;maxRun = 999999;
  }
  else if(period == -1){ // Full2011-Fall11-V9 NoJetId
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    scaleFactorLum     =3.553;minRun =      0;maxRun = 999999;
  }
  else if(period == 1){ // Full2011-Fall11-V9-5000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-06_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_5000ipb_71mb.root";
    scaleFactorLum     =5.098;minRun =      0;maxRun = 999999;
  }
  else if(period == 2){ //  Full2012-Summer12-V9-12000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-07_HCP_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_12p1ifb.root";
    scaleFactorLum     = 12.1;minRun =      0;maxRun = 999999;
  }
  else if(period == 3){ //  Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-09/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    scaleFactorLum     = 19.365;minRun =      0;maxRun = 999999;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open(fakePath.Data());
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
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

  TFile *fPUS4File = TFile::Open(puPath.Data());
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

  TFile *fLeptonFRFileSyst = TFile::Open(fakePath.Data());
  TH2D *fhDFRMuSystUp = (TH2D*)(fLeptonFRFileSyst->Get("MuonFakeRate_M2_ptThreshold50_PtEta"));
  TH2D *fhDFRMuSystDown = (TH2D*)(fLeptonFRFileSyst->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  TH2D *fhDFRElSystUp = (TH2D*)(fLeptonFRFileSyst->Get("ElectronFakeRate_V4_ptThreshold50_PtEta"));
  TH2D *fhDFRElSystDown = (TH2D*)(fLeptonFRFileSyst->Get("ElectronFakeRate_V4_ptThreshold20_PtEta"));
  assert(fhDFRMuSystUp);
  assert(fhDFRMuSystDown);
  assert(fhDFRElSystUp);
  assert(fhDFRElSystDown);
  fhDFRMuSystUp->SetDirectory(0);
  fhDFRMuSystDown->SetDirectory(0);
  fhDFRElSystUp->SetDirectory(0);
  fhDFRElSystDown->SetDirectory(0);
  fLeptonFRFileSyst->Close();
  delete fLeptonFRFileSyst;


  //----------------------------------------------------------------------------
  double theCutMassHigh	     = cutMassHigh (mH);
  double theCutPtMaxLow	     = cutPtMaxLow (mH);
  double theCutDeltaphilHigh = cutDeltaphiHigh (mH);
  double theCutMTLow         = cutMTLow (mH);
  double theCutMTHigh        = cutMTHigh (mH);

  cout << "theCutMassHigh: " << theCutMassHigh << endl;
  cout << "theCutPtMaxLow: " << theCutPtMaxLow << endl;
  cout << "theCutDeltaphilHigh: " << theCutDeltaphilHigh << endl;
  cout << "theCutMTLow: " << theCutMTLow << endl;
  cout << "theCutMTHigh: " << theCutMTHigh << endl;

  //----------------------------------------------------------------------------
  const int nHist = 6;
  int    nBinHis[nHist] = { 200,  200,  200,  200,  200,  200};
  double minHis[nHist]  = {-1.0, -1.0, -1.0, -0.0, -1.0,  0.0};
  double maxHis[nHist]  = { 1.0,  1.0,  1.0,  1.0,  1.0,200.0};


  //****************************************************************************
  // Yields and Histograms
  //****************************************************************************
  enum { kFakeElectron, kFakeMuon, kFakeLepton };
  enum { kMuMu, kEleEle, kEleMu, kMuEle, kSameFlavor, kDifferentFlavor, kAllFinalStates };
  enum { kZeroJetBin, kOneJetBin, kVBFBin };
  //pt eta bins: { [all],[pt<20,0-1],[pt<20,1-1.479][pt<20,1.479-2.5],[pt>20,0-1],[pt>20,1-1.479][pt>20,1.479-2.5]}

  //[Fake Electron/Muon][Final State][Jet Bin][Pt/Eta Bin]
  vector<vector<vector<vector<double> > > > FakeLeptonBkgYields;
  vector<vector<vector<vector<double> > > > FakeLeptonBkgYieldsErrSqr;
  vector<vector<vector<vector<double> > > > FakeLeptonBkgYieldsSystUp;
  vector<vector<vector<vector<double> > > > FakeLeptonBkgYieldsSystUpErrSqr;
  vector<vector<vector<vector<double> > > > FakeLeptonBkgYieldsSystDown;
  vector<vector<vector<vector<double> > > > FakeLeptonBkgYieldsSystDownErrSqr;

  //[Fake Electron/Muon][Final State][Jet Bin][Variable]
  vector<vector<vector<vector<TH1F*> > > > FakeLeptonBkgHists_Data;
  vector<vector<vector<vector<TH1F*> > > > FakeLeptonBkgHists_MC;

  for(UInt_t i = 0; i < 3; i++) {
    vector<vector<vector<double> > > FakeLeptonBkgYields_tmp1;
    vector<vector<vector<double> > > FakeLeptonBkgYieldsErrSqr_tmp1;
    vector<vector<vector<double> > > FakeLeptonBkgYieldsSystUp_tmp1;
    vector<vector<vector<double> > > FakeLeptonBkgYieldsSystUpErrSqr_tmp1;
    vector<vector<vector<double> > > FakeLeptonBkgYieldsSystDown_tmp1;
    vector<vector<vector<double> > > FakeLeptonBkgYieldsSystDownErrSqr_tmp1;
    vector<vector<vector<TH1F*> > > FakeLeptonBkgHists_Data_tmp1;    
    vector<vector<vector<TH1F*> > > FakeLeptonBkgHists_MC_tmp1;    
    for(UInt_t j = 0; j < 7; j++) {
    vector<vector<double> > FakeLeptonBkgYields_tmp2;
    vector<vector<double> > FakeLeptonBkgYieldsErrSqr_tmp2;
    vector<vector<double> > FakeLeptonBkgYieldsSystUp_tmp2;
    vector<vector<double> > FakeLeptonBkgYieldsSystUpErrSqr_tmp2;
    vector<vector<double> > FakeLeptonBkgYieldsSystDown_tmp2;
    vector<vector<double> > FakeLeptonBkgYieldsSystDownErrSqr_tmp2;
    vector<vector<TH1F*> > FakeLeptonBkgHists_Data_tmp2;    
    vector<vector<TH1F*> > FakeLeptonBkgHists_MC_tmp2;    
      for(UInt_t k = 0; k < 3; k++) {
        vector<double> FakeLeptonBkgYields_tmp3;
        vector<double> FakeLeptonBkgYieldsErrSqr_tmp3;
        vector<double> FakeLeptonBkgYieldsSystUp_tmp3;
        vector<double> FakeLeptonBkgYieldsSystUpErrSqr_tmp3;
        vector<double> FakeLeptonBkgYieldsSystDown_tmp3;
        vector<double> FakeLeptonBkgYieldsSystDownErrSqr_tmp3;
        vector<TH1F*> FakeLeptonBkgHists_Data_tmp3;    

        for(UInt_t l = 0; l < 7; l++) {
          FakeLeptonBkgYields_tmp3.push_back(0.0);
          FakeLeptonBkgYieldsErrSqr_tmp3.push_back(0.0);
          FakeLeptonBkgYieldsSystUp_tmp3.push_back(0.0);
          FakeLeptonBkgYieldsSystUpErrSqr_tmp3.push_back(0.0);
          FakeLeptonBkgYieldsSystDown_tmp3.push_back(0.0);
          FakeLeptonBkgYieldsSystDownErrSqr_tmp3.push_back(0.0);
        }
        
        string fakeLeptonLabel; 
        string finalStateLabel;
        string jetbinLabel;
        if (i==0) fakeLeptonLabel = "FakeElectron"; 
        else if (i==1) fakeLeptonLabel = "FakeMuon"; 
        else if (i==2) fakeLeptonLabel = "FakeLepton"; 
        else { cout << "fake lepton type index not recognized.\n"; assert(0); }
        if (j==0) finalStateLabel = "mm"; 
        else if (j==1) finalStateLabel = "ee"; 
        else if (j==2) finalStateLabel = "em"; 
        else if (j==3) finalStateLabel = "me"; 
        else if (j==4) finalStateLabel = "SF"; 
        else if (j==5) finalStateLabel = "DF"; 
        else if (j==6) finalStateLabel = "AllFinalStates"; 
        else { cout << "final state type index not recognized.\n"; assert(0); }
        if (k==0) jetbinLabel = "ZeroJetBin"; 
        else if (k==1) jetbinLabel = "OneJetBin"; 
        else if (k==2) jetbinLabel = "VBFBin"; 
        else { cout << "jet bin index not recognized.\n"; assert(0); }

        char hname[50];
        sprintf(hname,"hMVA_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());           FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 20, -1, 1));  FakeLeptonBkgHists_Data_tmp3[0]->Sumw2();
        sprintf(hname,"hFakeLeptonPt_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());  FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 10, 0, 100)); FakeLeptonBkgHists_Data_tmp3[1]->Sumw2();
        sprintf(hname,"hFakeLeptonEta_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str()); FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 10, 0, 2.5)); FakeLeptonBkgHists_Data_tmp3[2]->Sumw2();
        sprintf(hname,"hPtMax_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());         FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 10, 0, 100)); FakeLeptonBkgHists_Data_tmp3[3]->Sumw2();
        sprintf(hname,"hPtMin_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());         FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 10, 0, 100)); FakeLeptonBkgHists_Data_tmp3[4]->Sumw2();
        sprintf(hname,"hMet_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());           FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 15, 0, 150)); FakeLeptonBkgHists_Data_tmp3[5]->Sumw2();
        sprintf(hname,"hDeltaPhi_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());      FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 15, 0, 180)); FakeLeptonBkgHists_Data_tmp3[6]->Sumw2();
        sprintf(hname,"hDileptonMass_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());  FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 20, 0, 200)); FakeLeptonBkgHists_Data_tmp3[7]->Sumw2();
        sprintf(hname,"hMtHiggs_FakeLeptonBkg_Data_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());       FakeLeptonBkgHists_Data_tmp3.push_back(new TH1F(hname,"", 20, 0, 200)); FakeLeptonBkgHists_Data_tmp3[8]->Sumw2();

        vector<TH1F*> FakeLeptonBkgHists_MC_tmp3;    

        sprintf(hname,"hMVA_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());      
        FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 20, -1, 1));  FakeLeptonBkgHists_MC_tmp3[0]->Sumw2();
        sprintf(hname,"hFakeLeptonPt_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());  FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 10, 0, 100)); FakeLeptonBkgHists_MC_tmp3[1]->Sumw2();
        sprintf(hname,"hFakeLeptonEta_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str()); FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 10, 0, 2.5)); FakeLeptonBkgHists_MC_tmp3[2]->Sumw2();
        sprintf(hname,"hPtMax_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());         FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 10, 0, 100)); FakeLeptonBkgHists_MC_tmp3[3]->Sumw2();
        sprintf(hname,"hPtMin_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());         FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 10, 0, 100)); FakeLeptonBkgHists_MC_tmp3[4]->Sumw2();
        sprintf(hname,"hMet_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());           FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 15, 0, 150)); FakeLeptonBkgHists_MC_tmp3[5]->Sumw2();
        sprintf(hname,"hDeltaPhi_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());      FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 15, 0, 180)); FakeLeptonBkgHists_MC_tmp3[6]->Sumw2();
        sprintf(hname,"hDileptonMass_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());  FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 20, 0, 200)); FakeLeptonBkgHists_MC_tmp3[7]->Sumw2();
        sprintf(hname,"hMtHiggs_FakeLeptonBkg_MC_%s_%s_%s",fakeLeptonLabel.c_str(),finalStateLabel.c_str(),jetbinLabel.c_str());       FakeLeptonBkgHists_MC_tmp3.push_back(new TH1F(hname,"", 20, 0, 200)); FakeLeptonBkgHists_MC_tmp3[8]->Sumw2();
              
        FakeLeptonBkgYields_tmp2.push_back(FakeLeptonBkgYields_tmp3);
        FakeLeptonBkgYieldsErrSqr_tmp2.push_back(FakeLeptonBkgYieldsErrSqr_tmp3);
        FakeLeptonBkgYieldsSystUp_tmp2.push_back(FakeLeptonBkgYieldsSystUp_tmp3);
        FakeLeptonBkgYieldsSystUpErrSqr_tmp2.push_back(FakeLeptonBkgYieldsSystUpErrSqr_tmp3);
        FakeLeptonBkgYieldsSystDown_tmp2.push_back(FakeLeptonBkgYieldsSystDown_tmp3);
        FakeLeptonBkgYieldsSystDownErrSqr_tmp2.push_back(FakeLeptonBkgYieldsSystDownErrSqr_tmp3);
        FakeLeptonBkgHists_Data_tmp2.push_back(FakeLeptonBkgHists_Data_tmp3);
        FakeLeptonBkgHists_MC_tmp2.push_back(FakeLeptonBkgHists_MC_tmp3);
      }
      FakeLeptonBkgYields_tmp1.push_back(FakeLeptonBkgYields_tmp2);
      FakeLeptonBkgYieldsErrSqr_tmp1.push_back(FakeLeptonBkgYieldsErrSqr_tmp2);
      FakeLeptonBkgYieldsSystUp_tmp1.push_back(FakeLeptonBkgYieldsSystUp_tmp2);
      FakeLeptonBkgYieldsSystUpErrSqr_tmp1.push_back(FakeLeptonBkgYieldsSystUpErrSqr_tmp2);
      FakeLeptonBkgYieldsSystDown_tmp1.push_back(FakeLeptonBkgYieldsSystDown_tmp2);
      FakeLeptonBkgYieldsSystDownErrSqr_tmp1.push_back(FakeLeptonBkgYieldsSystDownErrSqr_tmp2);
      FakeLeptonBkgHists_Data_tmp1.push_back(FakeLeptonBkgHists_Data_tmp2);
      FakeLeptonBkgHists_MC_tmp1.push_back(FakeLeptonBkgHists_MC_tmp2);
      
    }
    FakeLeptonBkgYields.push_back(FakeLeptonBkgYields_tmp1);
    FakeLeptonBkgYieldsErrSqr.push_back(FakeLeptonBkgYieldsErrSqr_tmp1);
    FakeLeptonBkgYieldsSystUp.push_back(FakeLeptonBkgYieldsSystUp_tmp1);
    FakeLeptonBkgYieldsSystUpErrSqr.push_back(FakeLeptonBkgYieldsSystUpErrSqr_tmp1);
    FakeLeptonBkgYieldsSystDown.push_back(FakeLeptonBkgYieldsSystDown_tmp1);
    FakeLeptonBkgYieldsSystDownErrSqr.push_back(FakeLeptonBkgYieldsSystDownErrSqr_tmp1);
    FakeLeptonBkgHists_Data.push_back(FakeLeptonBkgHists_Data_tmp1);
    FakeLeptonBkgHists_MC.push_back(FakeLeptonBkgHists_MC_tmp1);
  }

  cout << "Start\n";

  for (UInt_t nJetsType=0; nJetsType < 3; nJetsType++) {
    
    //----------------------------------------------------------------------------
    UInt_t          cuts;
    Int_t           dstype;
    UInt_t          nvtx;
    Float_t         npu;
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
    Int_t           type;
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
    Float_t         bdtg_aux0 = 0.0;
    Float_t         bdtg_aux1 = 0.0;
    Float_t         bdtg_aux2 = 0.0;
    Float_t         higgsPt = -999;
    Float_t         bdtg_wjets = 0.0;
    //Float_t         knn_wjets = 0.0;
    Float_t         sfWeightPU;
    Float_t         sfWeightEff;
    Float_t         sfWeightTrig;

    background->SetBranchAddress( "cuts"          , &cuts 	  );
    background->SetBranchAddress( "dstype"        , &dstype	  );
    background->SetBranchAddress( "nvtx"          , &nvtx 	  );
    background->SetBranchAddress( "npu"           , &npu          );
    background->SetBranchAddress( "njets"         , &njets	  );
    background->SetBranchAddress( "run"           , &run          );
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
    background->SetBranchAddress(Form("bdt_hww%i_%djet_ww"      ,130,nJetsType), &bdt	  );
    background->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"     ,130,nJetsType), &bdtd	  );
    background->SetBranchAddress(Form("nn_hww%i_%djet_ww"       ,130,nJetsType), &nn	  );
    background->SetBranchAddress(Form("knn_hww%i_%djet_ww"      ,130,nJetsType), &knn	  );
    background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"     ,130,nJetsType), &bdtg	  );
    background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux0",130,nJetsType), &bdtg_aux0  );
    background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux1",130,nJetsType), &bdtg_aux1  );
    background->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux2",130,nJetsType), &bdtg_aux2  );


    for (UInt_t i=0; i<background->GetEntries(); i++) {

      background->GetEntry(i);
      if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

      if(dstype == SmurfTree::data &&
         (cuts & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
      if(dstype == SmurfTree::data && run <  minRun) continue;
      if(dstype == SmurfTree::data && run >  maxRun) continue;

      unsigned int Njet3 = njets;
      if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
        if(jet3->pt() <= 30)					           Njet3 = 2;
        else if(jet3->pt() > 30 && (
                  (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
                  (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 0;
        else							           Njet3 = 2;
        if(njets < 2 || njets > 3)                                           Njet3 = 0;
        if(TMath::Abs(jet1->eta()) >= 4.5 || TMath::Abs(jet2->eta()) >= 4.5) Njet3 = 0;
      }
      bool passJetCut[3] = {Njet3 == nJetsType, false, false};
      if(nJetsType == 0 && 			 jet1->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 5.0  			       ) passJetCut[1] = true;
      if(nJetsType == 0 && 			 jet1->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 5.0  			       ) passJetCut[2] = true;
      if(nJetsType == 1 && jet1->pt()*1.05 > 30 && jet2->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 5.0 && TMath::Abs(jet2->eta()) < 5.0) passJetCut[1] = true;
      if(nJetsType == 1 && jet1->pt()*0.95 > 30 && jet2->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 5.0 && TMath::Abs(jet2->eta()) < 5.0) passJetCut[2] = true;

      double minmet = TMath::Min(pmet,pTrackMet);
      bool passMET = minmet > 20. &&
        (minmet > 37.+nvtx/2.0 || type == SmurfTree::em || type == SmurfTree::me);

      bool passNewCuts = true;
      if(lep2->pt() <= 15 && (type == SmurfTree::mm||type == SmurfTree::ee)) passNewCuts = false;
      if(dilep->pt() <= 45) passNewCuts = false;

      // WW Preselection
      bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
        ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
        dstype != SmurfTree::data;


      if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
      if( passJetCut[0]==false&&passJetCut[1]==false&&passJetCut[2]==false ) continue; // select n-jet type events
      if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
       if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
       if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
      if( dilep->mass() <= 20.0  &&
          (type == SmurfTree::mm || 
           type == SmurfTree::ee)      					 ) continue; // cut on low dilepton mass for ee/mm
       if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
       if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
       if( passNewCuts == false                                             ) continue; // cut on new pt cuts
       if( passMET == false                                                 ) continue; // cut on pmet
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

      bool dPhiDiLepJetCut = true;
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
      if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

      int fDecay = 0;
      if     (dstype == SmurfTree::wjets  	 ) fDecay = 5;
      else if(dstype == SmurfTree::ttbar  	 ) fDecay = 3;
      else if(dstype == SmurfTree::dyee   	 ) fDecay = 4;
      else if(dstype == SmurfTree::dymm   	 ) fDecay = 4;
      else if(dstype == SmurfTree::dytt   	 ) fDecay = 7;
      else if(dstype == SmurfTree::tw     	 ) fDecay = 3;
      else if(dstype == SmurfTree::qqww   	 ) fDecay = 0;
      else if(dstype == SmurfTree::qqww2j        ) fDecay = 0;
      else if(dstype == SmurfTree::wz     	 ) fDecay = 2;
      else if(dstype == SmurfTree::zz     	 ) fDecay = 2;
      else if(dstype == SmurfTree::www           ) fDecay = 2;
      else if(dstype == SmurfTree::ggww   	 ) fDecay = 1;
      else if(dstype == SmurfTree::wgamma 	 ) fDecay = 6;
      else if(dstype == SmurfTree::wgstar 	 ) fDecay = 6;
      else if(dstype == SmurfTree::data   	 ) fDecay = 5;
      else if(dstype == SmurfTree::dyttDataDriven) fDecay = 7;
      else if(dstype == SmurfTree::qcd           ) fDecay = 7;
      else                                 {printf("bad dstype: %d\n",dstype); assert(0);}
      if(dstype == SmurfTree::wz || dstype == SmurfTree::zz) {
        if(lep1MotherMcId == 23 && lep2MotherMcId == 23) {
          fDecay = 4;
        }
      }


      //****************************************************************************************
      //Find the right bin
      //****************************************************************************************
      Int_t JetBinIndex = -1;
      if(njets == 0) JetBinIndex = kZeroJetBin;
      else if(njets == 1) JetBinIndex = kOneJetBin;
      else JetBinIndex = kVBFBin;
    
      Int_t FakeLeptonType = -1;
      Int_t FakeLeptonIndex = -1;
      Double_t FakeLeptonPt = 0;
      Double_t FakeLeptonEta = 0;
      Int_t FakeLeptonPtEtaBin = -1;

      double myWeight = 1.0;
      double myWeightSystUp = 1.0;
      double myWeightSystDown = 1.0;
      double add      = 1.0;
      double addSystUp      = 1.0;
      double addSystDown      = 1.0;
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
      bool isRealLepton = false;
      if((TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) &&
         (TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13)) isRealLepton = true;
      double addLepEff = 1.0;
      double addFR     = 1.0;
      double addFRSystUp     = 1.0;
      double addFRSystDown     = 1.0;
      if(nFake > 1){
        myWeight = 0.0;
        myWeightSystUp = 0.0;
        myWeightSystDown = 0.0;
      }
      else if(nFake == 1){
        if(dstype == SmurfTree::data){
          addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = addFR;


          addFRSystUp =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMuSystUp, fhDFRElSystUp, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFRSystUp = addFRSystUp*fakeRate(lep2->pt(), lep2->eta(), fhDFRMuSystUp, fhDFRElSystUp, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          addSystUp = addFRSystUp;
          addFRSystDown =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMuSystDown, fhDFRElSystDown, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFRSystDown = addFRSystDown*fakeRate(lep2->pt(), lep2->eta(), fhDFRMuSystDown, fhDFRElSystDown, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          addSystDown = addFRSystDown;

          fDecay  	      = 5;
          myWeight	      = add;
          myWeightSystUp      = addSystUp;
          myWeightSystDown    = addSystDown;

        }
        else if(isRealLepton == true || dstype == SmurfTree::wgamma){
          addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);


          addFRSystUp =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMuSystUp, fhDFRElSystUp, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFRSystUp = addFRSystUp*fakeRate(lep2->pt(), lep2->eta(), fhDFRMuSystUp, fhDFRElSystUp, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          addFRSystDown =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMuSystDown, fhDFRElSystDown, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFRSystDown = addFRSystDown*fakeRate(lep2->pt(), lep2->eta(), fhDFRMuSystDown, fhDFRElSystDown, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);


          add = addFR;
          addSystUp = addFRSystUp;
          addSystDown = addFRSystDown;
          add = add*nPUScaleFactor2012(fhDPUS4,npu);
          addSystUp = addSystUp*nPUScaleFactor2012(fhDPUS4,npu);
          addSystDown = addSystDown*nPUScaleFactor2012(fhDPUS4,npu);

          addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
            leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
          add = add*addLepEff;
          addSystUp   = addSystUp*addLepEff;
          addSystDown = addSystDown*addLepEff;

          add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                                            fabs(lep2->eta()), lep2->pt(), 
                                                            TMath::Abs( lid1), TMath::Abs(lid2));
          addSystUp = addSystUp*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                                            fabs(lep2->eta()), lep2->pt(), 
                                                            TMath::Abs( lid1), TMath::Abs(lid2));
          addSystDown = addSystDown*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                                            fabs(lep2->eta()), lep2->pt(), 
                                                            TMath::Abs( lid1), TMath::Abs(lid2));
          fDecay  	       = 5;
          myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
          myWeightSystUp       = -1.0 * scale1fb*scaleFactorLum*addSystUp;
          myWeightSystDown     = -1.0 * scale1fb*scaleFactorLum*addSystDown;
        }
        else {
          myWeight = 0.0;
          myWeightSystUp = 0.0;
          myWeightSystDown = 0.0;
        }

        if (myWeight > 0) {
//           cout << dstype << " : " << myWeight << " : " << scale1fb*scaleFactorLum << " " << add << " : " 
//                << fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) << " " 
//                << fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection, (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) << " "
//                << endl;
//           cout << "CUTS: " 
//                << ( MinPreselCut == false                                            )  << " "
//                << ( passJetCut[0]==false&&passJetCut[1]==false&&passJetCut[2]==false )  << " "
//                << ( dilep->mass() > dilmass_cut 					 )  << " "
//                << ( lq1*lq2 > 0                 					 )  << " "
//                << ( dilep->mass() <= 12.0       					 )  << " "
//                << ( dilep->mass() <= 20.0  &&
//                     (type == SmurfTree::mm || 
//                      type == SmurfTree::ee)      					 )  << " "
//                << ( lep1->pt() <= 20	    					 )  << " "
//                << ( lep2->pt() <= 10	    					 )  << " "
//                << ( passNewCuts == false                                             )  << " "
//                << ( passMET == false                                                 )  << " "
//                << (fabs(dilep->mass()-91.1876) <= 15 &&
//                    (type == SmurfTree::mm || 
//                     type == SmurfTree::ee)                                            )  << " "
//                << ( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)  << " "
//                << ( (cuts & patternTopTag) == patternTopTag                          )  << " "
//                << ( dPhiDiLepJetCut == false                                         )  << " "
//                << endl;

//           cout << pmet << " " << pTrackMet << " : " << minmet << " -- " << type << " : " << passMET << endl;
         }

        //*********************************
        //Find Fake Lepton Type
        //*********************************
        if ( (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) {
          FakeLeptonType = 13;
          FakeLeptonIndex = 0;         
        } else if ((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) {
          FakeLeptonType = 11;
          FakeLeptonIndex = 0;                 
        } else if ((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) {
          FakeLeptonType = 13;
          FakeLeptonIndex = 1;                 
        } else if ((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) {
          FakeLeptonType = 11;
          FakeLeptonIndex = 1;
        } else {
          cout << "No Fake type found \n";
        }  

        if (FakeLeptonIndex == 0) {
          FakeLeptonPt = lep1->pt();
          FakeLeptonEta = lep1->eta();
        }
        else if (FakeLeptonIndex == 1) {
          FakeLeptonPt = lep2->pt();
          FakeLeptonEta = lep2->eta();
        }

        if (FakeLeptonType == 11) {
            if (FakeLeptonPt <= 20) {
            if (fabs(FakeLeptonEta) < 1.0) {
              FakeLeptonPtEtaBin = 1;
            } else if (fabs(FakeLeptonEta) < 1.479) {
              FakeLeptonPtEtaBin = 2;
            } else if (fabs(FakeLeptonEta) <= 2.5) {
              FakeLeptonPtEtaBin = 3;
            } else assert(0);
          } else {
            if (fabs(FakeLeptonEta) < 1.0) {
              FakeLeptonPtEtaBin = 4;
            } else if (fabs(FakeLeptonEta) < 1.479) {
              FakeLeptonPtEtaBin = 5;
            } else if (fabs(FakeLeptonEta) <= 2.5) {
              FakeLeptonPtEtaBin = 6;
            } else assert(0);
          }
        } else if (FakeLeptonType == 13) {
          if (FakeLeptonPt <= 14.5) {
            if (fabs(FakeLeptonEta) < 1.479) {
              FakeLeptonPtEtaBin = 1;
            } else {
              FakeLeptonPtEtaBin = 2;
            } 
          } else if (FakeLeptonPt <= 20.0) {
            if (fabs(FakeLeptonEta) < 1.479) {
              FakeLeptonPtEtaBin = 3;
            } else {
              FakeLeptonPtEtaBin = 4;
            }           
          } else {
            if (fabs(FakeLeptonEta) < 1.479) {
              FakeLeptonPtEtaBin = 5;
            } else {
              FakeLeptonPtEtaBin = 6;
            }          
          }
        }
    
      }
      else if(dstype == SmurfTree::data) {
        myWeight = 0.0;
        myWeightSystUp   = 0.0;
        myWeightSystDown = 0.0;
      }
      else if(dstype== SmurfTree::dyttDataDriven || dstype == SmurfTree::qcd) {
        myWeight = 0.0;
        myWeightSystUp   = 0.0;
        myWeightSystDown = 0.0;
      }
      else if(dstype != SmurfTree::data){
        myWeight = 0.0;
        myWeightSystUp   = 0.0;
        myWeightSystDown = 0.0;
     }

      if(myWeight == 0) continue;
  

      double theCutPtMinLow = cutPtMinLow (mH, type);
      bool passAllCuts = dilep->mass()         < theCutMassHigh &&
        mt		     > theCutMTLow &&
        mt		     < theCutMTHigh &&
        lep1->pt()	     > theCutPtMaxLow &&
        lep2->pt()	     > theCutPtMinLow &&
        dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh &&
        passJetCut[0]==true;
//       if(nJetsType == 2){
//         int centrality = 0;
//         if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
//             (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
//            ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
//             (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
//         passAllCuts = (*jet1+*jet2).M() > 500. &&
//           TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
//           (mH > 200 || dilep->mass() < 100.) &&
//           centrality == 1 &&
//           passJetCut[0]==true;
//       }    

      if(passAllCuts == true) {
        double newWeight = myWeight;
        double newWeightSystUp = myWeightSystUp;
        double newWeightSystDown = myWeightSystDown;

        if((fDecay == 0 || fDecay == 1) && wwPresel == false){ // only for WW
          if(njets == 0) {
            newWeight=newWeight*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),0)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0);
            newWeightSystUp=newWeightSystUp*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),0)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0);
            newWeightSystDown=newWeightSystDown*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),0)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0);
          }
          else {
            newWeight=newWeight*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),1)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1);	   
            newWeightSystUp=newWeightSystUp*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),1)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1);	   
            newWeightSystDown=newWeightSystDown*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),1)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1);	   
          }
        }
        if((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
           (type   == SmurfTree::mm   || type   == SmurfTree::ee)){
          if(nJetsType != 2){
            newWeight=newWeight*DYBkgScaleFactor(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2))/DYBkgScaleFactor(0,TMath::Min((int)nJetsType,2));
            newWeightSystUp=newWeightSystUp*DYBkgScaleFactor(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2))/DYBkgScaleFactor(0,TMath::Min((int)nJetsType,2));
            newWeightSystDown=newWeightSystDown*DYBkgScaleFactor(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2))/DYBkgScaleFactor(0,TMath::Min((int)nJetsType,2));
          }
        }
        else if(fDecay == 4){
        }

        if (type == SmurfTree::mm) {
          FakeLeptonBkgYields[kFakeLepton][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYields[kFakeLepton][kMuMu][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuMu][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

          FakeLeptonBkgYieldsSystUp[kFakeLepton][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kMuMu][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kMuMu][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kMuMu][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kMuMu][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;

        }
        if (type == SmurfTree::ee) {
          FakeLeptonBkgYields[kFakeLepton][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYields[kFakeLepton][kEleEle][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleEle][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

          FakeLeptonBkgYieldsSystUp[kFakeLepton][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kEleEle][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kEleEle][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kEleEle][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kEleEle][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
        }
        if (type == SmurfTree::em) {
          FakeLeptonBkgYields[kFakeLepton][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYields[kFakeLepton][kEleMu][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kDifferentFlavor][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleMu][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kDifferentFlavor][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

          FakeLeptonBkgYieldsSystUp[kFakeLepton][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kEleMu][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kEleMu][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kEleMu][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kEleMu][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
        }
        if (type == SmurfTree::me) {

          FakeLeptonBkgYields[kFakeLepton][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
          FakeLeptonBkgYields[kFakeLepton][kMuEle][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kDifferentFlavor][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuEle][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kDifferentFlavor][JetBinIndex][0] += newWeight*newWeight;
          FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

          FakeLeptonBkgYieldsSystUp[kFakeLepton][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kMuEle][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUp[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kMuEle][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystUpErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kMuEle][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDown[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kMuEle][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          FakeLeptonBkgYieldsSystDownErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
        }

        if (FakeLeptonType == 11) {
          if (type == SmurfTree::mm) {
            FakeLeptonBkgYields[kFakeElectron][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeElectron][kMuMu][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

            FakeLeptonBkgYieldsSystUp[kFakeElectron][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kMuMu][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kMuMu][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kMuMu][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kMuMu][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgYields[kFakeElectron][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeElectron][kEleEle][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;
 
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kEleEle][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kEleEle][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kEleEle][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kEleEle][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
         }
          if (type == SmurfTree::em) {
            FakeLeptonBkgYields[kFakeElectron][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeElectron][kEleMu][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kDifferentFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kDifferentFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;
 
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kEleMu][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kEleMu][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kEleMu][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kEleMu][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
         }
          if (type == SmurfTree::me) {
            FakeLeptonBkgYields[kFakeElectron][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeElectron][kMuEle][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kDifferentFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kDifferentFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

            FakeLeptonBkgYieldsSystUp[kFakeElectron][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kMuEle][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kMuEle][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kMuEle][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kMuEle][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          }
        }
        if (FakeLeptonType == 13) {
          if (type == SmurfTree::mm) {
            FakeLeptonBkgYields[kFakeMuon][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeMuon][kMuMu][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

            FakeLeptonBkgYieldsSystUp[kFakeMuon][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kMuMu][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kMuMu][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kMuMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kMuMu][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kMuMu][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgYields[kFakeMuon][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeMuon][kEleEle][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

            FakeLeptonBkgYieldsSystUp[kFakeMuon][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kEleEle][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kEleEle][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kEleEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kEleEle][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kEleEle][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgYields[kFakeMuon][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeMuon][kEleMu][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kDifferentFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kDifferentFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

            FakeLeptonBkgYieldsSystUp[kFakeMuon][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kEleMu][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kEleMu][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kEleMu][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kEleMu][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kEleMu][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgYields[kFakeMuon][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kDifferentFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeight*newWeight;
            FakeLeptonBkgYields[kFakeMuon][kMuEle][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kDifferentFlavor][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kDifferentFlavor][JetBinIndex][0] += newWeight*newWeight;
            FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeight*newWeight;

            FakeLeptonBkgYieldsSystUp[kFakeMuon][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kMuEle][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kMuEle][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystUp*newWeightSystUp;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kMuEle][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][FakeLeptonPtEtaBin] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kMuEle][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kMuEle][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
            FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][0] += newWeightSystDown*newWeightSystDown;
          }
        }

        for (UInt_t l=1; l < 9; ++l) {
          Double_t varValue = -9999;
          if (l==1) {
            varValue = FakeLeptonPt;
          }
          if (l==2) {
            varValue = fabs(FakeLeptonEta);
          }
          if (l==3) {
            varValue = lep1->pt();
          }
          if (l==4) {
            varValue = lep2->pt();
          }
          if (l==5) {
            varValue = minmet;
          }
          if (l==6) {
            varValue = dPhi*180.0/TMath::Pi();
          }
          if (l==7) {
            varValue = dilep->mass();
          }
          if (l==8) {
            varValue = mt;
          }
      
          if (type == SmurfTree::mm) {
            FakeLeptonBkgHists_Data[kFakeLepton][kMuMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgHists_Data[kFakeLepton][kEleEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgHists_Data[kFakeLepton][kEleMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgHists_Data[kFakeLepton][kMuEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }

          if (FakeLeptonType == 11) {
            if (type == SmurfTree::mm) {
              FakeLeptonBkgHists_Data[kFakeElectron][kMuMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::ee) {
              FakeLeptonBkgHists_Data[kFakeElectron][kEleEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::em) {
              FakeLeptonBkgHists_Data[kFakeElectron][kEleMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::me) {
              FakeLeptonBkgHists_Data[kFakeElectron][kMuEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
          }
          if (FakeLeptonType == 13) {
            if (type == SmurfTree::mm) {
              FakeLeptonBkgHists_Data[kFakeMuon][kMuMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::ee) {
              FakeLeptonBkgHists_Data[kFakeMuon][kEleEle][JetBinIndex][l]->Fill(varValue, newWeight);
              FakeLeptonBkgHists_Data[kFakeMuon][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);            
            }        
            if (type == SmurfTree::em) {
              FakeLeptonBkgHists_Data[kFakeMuon][kEleMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::me) {
              FakeLeptonBkgHists_Data[kFakeMuon][kMuEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }        
          }     
        }
      }  

      bool passMVAPreselCuts = mt > 80 && mt < mH; if(wwPresel == true) passMVAPreselCuts = true;
      if(passMVAPreselCuts == true && passJetCut[0] == true){

        if (type == SmurfTree::mm) {
          FakeLeptonBkgHists_Data[kFakeLepton][kMuMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }
        if (type == SmurfTree::ee) {
          FakeLeptonBkgHists_Data[kFakeLepton][kEleEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }
        if (type == SmurfTree::em) {
          FakeLeptonBkgHists_Data[kFakeLepton][kEleMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }
        if (type == SmurfTree::me) {
          FakeLeptonBkgHists_Data[kFakeLepton][kMuEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_Data[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }

        if (FakeLeptonType == 11) {
          if (type == SmurfTree::mm) {
            FakeLeptonBkgHists_Data[kFakeElectron][kMuMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgHists_Data[kFakeElectron][kEleEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgHists_Data[kFakeElectron][kEleMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgHists_Data[kFakeElectron][kMuEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
        }
        if (FakeLeptonType == 13) {
          if (type == SmurfTree::mm) {
            FakeLeptonBkgHists_Data[kFakeMuon][kMuMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgHists_Data[kFakeMuon][kEleEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgHists_Data[kFakeMuon][kEleMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgHists_Data[kFakeMuon][kMuEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_Data[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
        }
      } // passMVAPreselCuts
    }

    printf("--- Finished Bgdnal loop\n");
    wjetsMC->SetBranchAddress( "cuts"          , &cuts 	  );
    wjetsMC->SetBranchAddress( "dstype"        , &dstype	  );
    wjetsMC->SetBranchAddress( "nvtx"          , &nvtx 	  );
    wjetsMC->SetBranchAddress( "npu"           , &npu 	          );
    wjetsMC->SetBranchAddress( "njets"         , &njets	  );
    wjetsMC->SetBranchAddress( "run"           , &run	          );
    wjetsMC->SetBranchAddress( "event"         , &event	  );
    wjetsMC->SetBranchAddress( "scale1fb"      , &scale1fb	  );
    wjetsMC->SetBranchAddress( "lep1"          , &lep1 	  );
    wjetsMC->SetBranchAddress( "lep2"          , &lep2 	  );
    wjetsMC->SetBranchAddress( "jet1"          , &jet1 	  );
    wjetsMC->SetBranchAddress( "jet2"          , &jet2 	  );
    wjetsMC->SetBranchAddress( "jet3"          , &jet3 	  );
    wjetsMC->SetBranchAddress( "dPhi"          , &dPhi 	  );
    wjetsMC->SetBranchAddress( "dR"            , &dR		  );
    wjetsMC->SetBranchAddress( "dilep"         , &dilep	  );
    wjetsMC->SetBranchAddress( "type"          , &type 	  );
    wjetsMC->SetBranchAddress( "pmet"          , &pmet 	  );
    wjetsMC->SetBranchAddress( "pTrackMet"     , &pTrackMet	  );
    wjetsMC->SetBranchAddress( "met"           , &met  	  );
    wjetsMC->SetBranchAddress( "mt"            , &mt		  );
    wjetsMC->SetBranchAddress( "mt1"           , &mt1  	  );
    wjetsMC->SetBranchAddress( "mt2"           , &mt2  	  );
    wjetsMC->SetBranchAddress( "dPhiLep1MET"   , &dPhiLep1MET    );
    wjetsMC->SetBranchAddress( "dPhiLep2MET"   , &dPhiLep2MET    );
    wjetsMC->SetBranchAddress( "dPhiDiLepMET"  , &dPhiDiLepMET   );
    wjetsMC->SetBranchAddress( "dPhiDiLepJet1" , &dPhiDiLepJet1  );
    wjetsMC->SetBranchAddress( "lq1"           , &lq1  	  );
    wjetsMC->SetBranchAddress( "lq2"           , &lq2  	  );
    wjetsMC->SetBranchAddress( "lid1"          , &lid1 	  );
    wjetsMC->SetBranchAddress( "lid2"          , &lid2 	  );
    wjetsMC->SetBranchAddress( "lid3"          , &lid3 	  );
    wjetsMC->SetBranchAddress( "processId"     , &processId	  );
    wjetsMC->SetBranchAddress( "jetLowBtag"    , &jetLowBtag	  );
    wjetsMC->SetBranchAddress( "nSoftMuons"    , &nSoftMuons	  );
    wjetsMC->SetBranchAddress( "jet1Btag"      , &jet1Btag	  );
    wjetsMC->SetBranchAddress( "jet2Btag"      , &jet2Btag	  );
    wjetsMC->SetBranchAddress( "lep1McId"      , &lep1McId	  );
    wjetsMC->SetBranchAddress( "lep2McId"      , &lep2McId	  );
    wjetsMC->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId );
    wjetsMC->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId );
    wjetsMC->SetBranchAddress(Form("bdt_hww%i_%djet_ww"      ,130,nJetsType), &bdt	  );
    wjetsMC->SetBranchAddress(Form("bdtd_hww%i_%djet_ww"     ,130,nJetsType), &bdtd	  );
    wjetsMC->SetBranchAddress(Form("nn_hww%i_%djet_ww"       ,130,nJetsType), &nn	  );
    wjetsMC->SetBranchAddress(Form("knn_hww%i_%djet_ww"      ,130,nJetsType), &knn	  );
    wjetsMC->SetBranchAddress(Form("bdtg_hww%i_%djet_ww"     ,130,nJetsType), &bdtg	  );
    wjetsMC->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux0",130,nJetsType), &bdtg_aux0  );
    wjetsMC->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux1",130,nJetsType), &bdtg_aux1  );
    wjetsMC->SetBranchAddress(Form("bdtg_hww%i_%djet_ww_aux2",130,nJetsType), &bdtg_aux2  );

    for (UInt_t i=0; i<wjetsMC->GetEntries(); i++) {

      wjetsMC->GetEntry(i);
      if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)wjetsMC->GetEntries());

      if(dstype != SmurfTree::wjets ) continue;
      if(dstype == SmurfTree::data && run <  minRun) continue;
      if(dstype == SmurfTree::data && run >  maxRun) continue;

      unsigned int Njet3 = njets;
      if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
        if(jet3->pt() <= 30)					           Njet3 = 2;
        else if(jet3->pt() > 30 && (
                  (jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
                  (jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 0;
        else							           Njet3 = 2;
        if(njets < 2 || njets > 3)                                           Njet3 = 0;
        if(TMath::Abs(jet1->eta()) >= 4.5 || TMath::Abs(jet2->eta()) >= 4.5) Njet3 = 0;
      }
      bool passJetCut[3] = {Njet3 == nJetsType, false, false};
      if(nJetsType == 0 && 			 jet1->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 5.0  			       ) passJetCut[1] = true;
      if(nJetsType == 0 && 			 jet1->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 5.0  			       ) passJetCut[2] = true;
      if(nJetsType == 1 && jet1->pt()*1.05 > 30 && jet2->pt()*1.05 < 30 && TMath::Abs(jet1->eta()) < 5.0 && TMath::Abs(jet2->eta()) < 5.0) passJetCut[1] = true;
      if(nJetsType == 1 && jet1->pt()*0.95 > 30 && jet2->pt()*0.95 < 30 && TMath::Abs(jet1->eta()) < 5.0 && TMath::Abs(jet2->eta()) < 5.0) passJetCut[2] = true;

      double minmet = TMath::Min(pmet,pTrackMet);
      bool passMET = minmet > 20. &&
        (minmet > 37.+nvtx/2.0 || type == SmurfTree::em || type == SmurfTree::me);

      bool passNewCuts = true;
      if(lep2->pt() <= 15 && (type == SmurfTree::mm||type == SmurfTree::ee)) passNewCuts = false;
      if(dilep->pt() <= 45) passNewCuts = false;

      // WW Preselection
      bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
        ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
        dstype != SmurfTree::data;

      if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
      if( passJetCut[0]==false&&passJetCut[1]==false&&passJetCut[2]==false ) continue; // select n-jet type events
      if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
      if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
      if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
      if( dilep->mass() <= 20.0  &&
          (type == SmurfTree::mm || 
           type == SmurfTree::ee)      					 ) continue; // cut on low dilepton mass for ee/mm
      if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
      if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
      if( passNewCuts == false                                             ) continue; // cut on new pt cuts
      if( passMET == false                                                 ) continue; // cut on pmet
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

      bool dPhiDiLepJetCut = true;
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
      if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

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
      else if(dstype == SmurfTree::www           ) fDecay = 2;
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

      //****************************************************************************************
      //Find the right bin
      //****************************************************************************************
      Int_t JetBinIndex = -1;
      if(njets == 0) JetBinIndex = kZeroJetBin;
      else if(njets == 1) JetBinIndex = kOneJetBin;
      else JetBinIndex = kVBFBin;
    
      Int_t FakeLeptonType = -1;
      Int_t FakeLeptonIndex = -1;
      Int_t FakeLeptonPtEtaBin = -1;
      Double_t FakeLeptonPt = 0;
      Double_t FakeLeptonEta = 0;
 
      if(!(TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) && (TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13)) {
        if (type == SmurfTree::ee || type == SmurfTree::em) FakeLeptonType = 11;
        else if (type == SmurfTree::mm || type == SmurfTree::me) FakeLeptonType = 13;
        FakeLeptonIndex = 0;
        FakeLeptonPt = lep1->pt();
        FakeLeptonEta = lep1->eta();            
      } else if(!(TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13) && (TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) ) {
        if (type == SmurfTree::ee || type == SmurfTree::me) FakeLeptonType = 11;
        else if (type == SmurfTree::mm || type == SmurfTree::em) FakeLeptonType = 13;
        FakeLeptonIndex = 1;
        FakeLeptonPt = lep2->pt();
        FakeLeptonEta = lep2->eta();            
      }

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

      bool isRealLepton = false;
      if((TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) &&
         (TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13)) isRealLepton = true;
      double addLepEff = 1.0;
      double addFR     = 1.0;


      //Both leptons pass selection
      if( (cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) {
      
        add = 1.0;
        add = add*nPUScaleFactor2012(fhDPUS4,npu);

        addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
          leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
        add = add*addLepEff;
        add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                                          fabs(lep2->eta()), lep2->pt(), 
                                                          TMath::Abs( lid1), TMath::Abs(lid2));
        if(fDecay == 4  && (type   == SmurfTree::mm   || type   == SmurfTree::ee)
           && (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm)) {
          if(njets == 0) add=add*DYBkgScaleFactor(0,0); 
          if(njets == 1) add=add*DYBkgScaleFactor(0,1); 
          if(njets >= 2) add=add*DYBkgScaleFactor(0,2); 
        }
        if(fDecay == 3) {
          if(njets == 0) add=add*TopBkgScaleFactor(0);
          if(njets == 1) add=add*TopBkgScaleFactor(1); 
          if(njets >= 2) add=add*TopBkgScaleFactor(2); 
        }

        if(dstype == SmurfTree::wgstar) add=add*WGstarScaleFactor(type,met);

        if((fDecay == 0 || fDecay == 1) && wwPresel == false){     
          if(njets == 0) add=add*WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0); 
          else	       add=add*WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1); 
        }
        // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
        // add = 1.0;
        myWeight = scale1fb*scaleFactorLum*add;
      }

      if(myWeight == 0) continue;

      double theCutPtMinLow = cutPtMinLow (mH, type);
      bool passAllCuts = dilep->mass()         < theCutMassHigh &&
        mt		     > theCutMTLow &&
        mt		     < theCutMTHigh &&
        lep1->pt()	     > theCutPtMaxLow &&
        lep2->pt()	     > theCutPtMinLow &&
        dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh &&
        passJetCut[0]==true;
      if(nJetsType == 2){
        int centrality = 0;
        if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
            (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
           ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
            (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
        passAllCuts = (*jet1+*jet2).M() > 500. &&
          TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
          (mH > 200 || dilep->mass() < 100.) &&
          centrality == 1 &&
          passJetCut[0]==true;
      }    

      if(passAllCuts == true) {
        double newWeight = myWeight;
        if((fDecay == 0 || fDecay == 1) && wwPresel == false){ // only for WW
          if(njets == 0) newWeight=newWeight*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),0)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),0);
          else           newWeight=newWeight*WWBkgScaleFactorCutBased(TMath::Max((int)mH,115),1)/WWBkgScaleFactorMVA(TMath::Max((int)mH,115),1);	   
        }
        if((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
           (type   == SmurfTree::mm   || type   == SmurfTree::ee)){
          if(nJetsType != 2){
            newWeight=newWeight*DYBkgScaleFactor(TMath::Max((int)mH,115),TMath::Min((int)nJetsType,2))/DYBkgScaleFactor(0,TMath::Min((int)nJetsType,2));
          }
        }
        else if(fDecay == 4){
        }



        for (UInt_t l=1; l < 9; ++l) {
          Double_t varValue = -9999;
          if (l==1) {
            varValue = FakeLeptonPt;
          }
          if (l==2) {
            varValue = fabs(FakeLeptonEta);
          }
          if (l==3) {
            varValue = lep1->pt();
          }
          if (l==4) {
            varValue = lep2->pt();
          }
          if (l==5) {
            varValue = minmet;
          }
          if (l==6) {
            varValue = dPhi*180.0/TMath::Pi();
          }
          if (l==7) {
            varValue = dilep->mass();
          }
          if (l==8) {
            varValue = mt;
          }

          if (type == SmurfTree::mm) {
            FakeLeptonBkgHists_MC[kFakeLepton][kMuMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgHists_MC[kFakeLepton][kEleEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgHists_MC[kFakeLepton][kEleMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgHists_MC[kFakeLepton][kMuEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
            FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
          }

          if (FakeLeptonType == 11) {
            if (type == SmurfTree::mm) {
              FakeLeptonBkgHists_MC[kFakeElectron][kMuMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::ee) {
              FakeLeptonBkgHists_MC[kFakeElectron][kEleEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::em) {
              FakeLeptonBkgHists_MC[kFakeElectron][kEleMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::me) {
              FakeLeptonBkgHists_MC[kFakeElectron][kMuEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
          }
          if (FakeLeptonType == 13) {
            if (type == SmurfTree::mm) {
              FakeLeptonBkgHists_MC[kFakeMuon][kMuMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::ee) {
              FakeLeptonBkgHists_MC[kFakeMuon][kEleEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kSameFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::em) {
              FakeLeptonBkgHists_MC[kFakeMuon][kEleMu][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
            if (type == SmurfTree::me) {
              FakeLeptonBkgHists_MC[kFakeMuon][kMuEle][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kDifferentFlavor][JetBinIndex][l]->Fill(varValue,  newWeight);  
              FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][l]->Fill(varValue,  newWeight);  
            }
          }     
        }
      }


      bool passMVAPreselCuts = mt > 80 && mt < mH; if(wwPresel == true) passMVAPreselCuts = true;
      if(passMVAPreselCuts == true && passJetCut[0] == true){


        if (type == SmurfTree::mm) {
          FakeLeptonBkgHists_MC[kFakeLepton][kMuMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }
        if (type == SmurfTree::ee) {
          FakeLeptonBkgHists_MC[kFakeLepton][kEleEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }
        if (type == SmurfTree::em) {
          FakeLeptonBkgHists_MC[kFakeLepton][kEleMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }
        if (type == SmurfTree::me) {
          FakeLeptonBkgHists_MC[kFakeLepton][kMuEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          FakeLeptonBkgHists_MC[kFakeLepton][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
        }

        if (FakeLeptonType == 11) {
          if (type == SmurfTree::mm) {
            FakeLeptonBkgHists_MC[kFakeElectron][kMuMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgHists_MC[kFakeElectron][kEleEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgHists_MC[kFakeElectron][kEleMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgHists_MC[kFakeElectron][kMuEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeElectron][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
        }
        if (FakeLeptonType == 13) {
          if (type == SmurfTree::mm) {
            FakeLeptonBkgHists_MC[kFakeMuon][kMuMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::ee) {
            FakeLeptonBkgHists_MC[kFakeMuon][kEleEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kSameFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::em) {
            FakeLeptonBkgHists_MC[kFakeMuon][kEleMu][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
          if (type == SmurfTree::me) {
            FakeLeptonBkgHists_MC[kFakeMuon][kMuEle][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kDifferentFlavor][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
            FakeLeptonBkgHists_MC[kFakeMuon][kAllFinalStates][JetBinIndex][0]->Fill(TMath::Max(TMath::Min((double)bdtg,1.0-0.001),-1.0+0.001),  myWeight);
          }
        }
      } // passMVAPreselCuts
    }
    printf("--- Finished Bgdnal loop\n");
  } // loop over NJetBins


  //**********************************************************************
  //Save Fake Lepton Bkg Plots
  //**********************************************************************
  TFile *fileOutput = new TFile("FakeLeptonBkgPlots.root", "RECREATE");
  for(UInt_t i = 0; i < 3; ++i) {
    for(UInt_t j = 0; j < 7; ++j) {
      for(UInt_t k = 0; k < 3; ++k) {
        for(UInt_t l = 0; l < 9; ++l) {
          fileOutput->WriteTObject(FakeLeptonBkgHists_Data[i][j][k][l], FakeLeptonBkgHists_Data[i][j][k][l]->GetName(), "WriteDelete");
          fileOutput->WriteTObject(FakeLeptonBkgHists_MC[i][j][k][l], FakeLeptonBkgHists_MC[i][j][k][l]->GetName(), "WriteDelete");
        }
      }
    }
  }
  fileOutput->Close();

  //**********************************************************************
  //Print Tables
  //**********************************************************************
  ofstream fResultTexTable("FakeLeptonBkgEstimate.tex");
  char buffer[200];

  fResultTexTable << setw(30) << left << "Analysis/Final State"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$\\mu\\mu$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$\\mu$e"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "e$\\mu$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "ee"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Total"
                  << " \\\\"
                  << endl;

  fResultTexTable << "\\hline" << endl;

  fResultTexTable << setw(30) << left << "0-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kMuMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kMuEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kEleMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kEleEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;

  fResultTexTable << " \\\\"
                  << endl;



  fResultTexTable << setw(30) << left << "1-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kMuMu][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuMu][kOneJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kMuEle][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuEle][kOneJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kEleMu][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleMu][kOneJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kEleEle][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleEle][kOneJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][0]));
  fResultTexTable << setw(20) << left << buffer;

  fResultTexTable << " \\\\"
                  << endl;



  fResultTexTable << setw(30) << left << "2-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kMuMu][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuMu][kVBFBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kMuEle][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kMuEle][kVBFBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kEleMu][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleMu][kVBFBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kEleEle][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kEleEle][kVBFBin][0]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeLepton][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeLepton][kAllFinalStates][kVBFBin][0]));
  fResultTexTable << setw(20) << left << buffer;

  fResultTexTable << " \\\\"
                  << endl;

  fResultTexTable << "\\hline" << endl;

    

  //**********************************************************************
  //Print More Detailed Tables
  //**********************************************************************
  ofstream fResultTexTableDetailed("FakeLeptonBkgEstimate.MoreBins.tex");

  fResultTexTableDetailed << "**************\n" 
                  << "Fake Electrons\n"
                  << "**************\n";


  fResultTexTableDetailed << setw(30) << left << "Analysis/Final State"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$\\mu\\mu$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$\\mu$e"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "e$\\mu$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "ee"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Total"
                  << " \\\\"
                  << endl;

  fResultTexTableDetailed << "\\hline" << endl;

  fResultTexTableDetailed << setw(30) << left << "0-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;

  fResultTexTableDetailed << " \\\\"
                  << endl;



  fResultTexTableDetailed << setw(30) << left << "1-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;

  fResultTexTableDetailed << " \\\\"
                  << endl;



  fResultTexTableDetailed << setw(30) << left << "2-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;

  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << "\\hline" << endl;

  //*****************************************************************************************************
  //*****************************************************************************************************
  //*****************************************************************************************************
  fResultTexTableDetailed << "**************\n" 
                          << "Fake Muons\n"
                          << "**************\n";


  fResultTexTableDetailed << setw(30) << left << "Analysis/Final State"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$\\mu\\mu$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$\\mu$e"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "e$\\mu$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "ee"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Total"
                  << " \\\\"
                  << endl;

  fResultTexTableDetailed << "\\hline" << endl;

  fResultTexTableDetailed << setw(30) << left << "0-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;

  fResultTexTableDetailed << " \\\\"
                  << endl;



  fResultTexTableDetailed << setw(30) << left << "1-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kOneJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;

  fResultTexTableDetailed << " \\\\"
                  << endl;



  fResultTexTableDetailed << setw(30) << left << "2-Jet Bin"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kVBFBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;

  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << "\\hline" << endl;






  //*****************************************************************************************************
  //*****************************************************************************************************
  //*****************************************************************************************************
  fResultTexTableDetailed << "**************\n" 
                          << "Fake Electrons \n"
                          << "**************\n";




  fResultTexTableDetailed << setw(30) << left << "FinalState/Kinematic Bin"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} \\le 20 , |\\eta| < 1.0$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} \\le  20 , 1.0 \\le |\\eta| < 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} \\le  20 , 1.479 \\le |\\eta| < 2.5$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} > 20 , |\\eta| < 1.0$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} > 20 , 1.0 \\le |\\eta| < 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} > 20 , 1.479 \\le |\\eta| < 2.5$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Total"
                  << " \\\\"
                  << endl;

  fResultTexTableDetailed << "\\hline" << endl;

  fResultTexTableDetailed << setw(30) << left << "$\\mu\\mu$"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "ee"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "e$\\mu$"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kEleMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kEleMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "$\\mu$e"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kMuEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kMuEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "All Final States"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;






  //*****************************************************************************************************
  //*****************************************************************************************************
  //*****************************************************************************************************
  fResultTexTableDetailed << "**************\n" 
                          << "Fake Muons \n"
                          << "**************\n";




  fResultTexTableDetailed << setw(30) << left << "FinalState/Kinematic Bin"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} \\le 14.5 , |\\eta| < 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} \\le 14.5 , |\\eta| \\ge 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ 14.5 > p_{T} \\le 20 , |\\eta| < 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ 14.5 > p_{T} \\le 20 , |\\eta| \\ge 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} > 20 , |\\eta| < 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "$ p_{T} > 20 , |\\eta| \\ge 1.479$"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Total"
                  << " \\\\"
                  << endl;

  fResultTexTableDetailed << "\\hline" << endl;

  fResultTexTableDetailed << setw(30) << left << "$\\mu\\mu$"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "ee"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "e$\\mu$"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kEleMu][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kEleMu][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "$\\mu$e"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kMuEle][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kMuEle][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  fResultTexTableDetailed << setw(30) << left << "All Final States"
                  << setw(3) << left << "&";
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][1],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][1]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][2],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][2]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][3],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][3]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][4],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][4]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][5],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][5]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][6],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][6]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << setw(3) << left << "&" ;
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][0]));
  fResultTexTableDetailed << setw(20) << left << buffer;
  fResultTexTableDetailed << " \\\\"
                  << endl;

  ofstream fResultSystTexTable("FakeLeptonBkgEstimateJetPtSpectrumSystematics.tex");

  fResultSystTexTable << "**************\n" 
                      << "Fake Muons \n"
                      << "**************\n";
  fResultSystTexTable << setw(45) << left << "JetBin"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "0-Jet"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "1-Jet"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "2-Jet"
                  << " \\\\"
                  << endl;

  fResultSystTexTable << "\\hline" << endl;

  fResultSystTexTable << setw(45) << left << "Nominal Fake Rates"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kOneJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeMuon][kAllFinalStates][kVBFBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << setw(45) << left << "Jet pT Spectrum Systematics Up"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][kOneJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystUpErrSqr[kFakeMuon][kAllFinalStates][kVBFBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << setw(45) << left << "Jet pT Spectrum Systematics Down"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][kZeroJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][kOneJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystDownErrSqr[kFakeMuon][kAllFinalStates][kVBFBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << setw(45) << left << "Systematic Uncertainty"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"+ %.0f - %.0f", 
          100*fabs(FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][kZeroJetBin][0] - 
               FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0])
          /FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0],
          100*fabs(FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][kZeroJetBin][0] - 
               FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0])
          /FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kZeroJetBin][0]
    );
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"+%.0f - %.0f", 
          100*fabs(FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][kOneJetBin][0] - 
               FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kOneJetBin][0])
          /FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kOneJetBin][0],
          100*fabs(FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][kOneJetBin][0] - 
               FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kOneJetBin][0])
          /FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kOneJetBin][0]
    );
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"+%.0f - %.0f", 
          100*fabs(FakeLeptonBkgYieldsSystUp[kFakeMuon][kAllFinalStates][kVBFBin][0] - 
               FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kVBFBin][0])
          /FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kVBFBin][0],
          100*fabs(FakeLeptonBkgYieldsSystDown[kFakeMuon][kAllFinalStates][kVBFBin][0] - 
               FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kVBFBin][0])
          /FakeLeptonBkgYields[kFakeMuon][kAllFinalStates][kVBFBin][0]
    );
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << "**************\n" 
                      << "Fake Electrons \n"
                      << "**************\n";
  fResultSystTexTable << setw(45) << left << "JetBin"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "0-Jet"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "1-Jet"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "2-Jet"
                  << " \\\\"
                  << endl;

  fResultSystTexTable << "\\hline" << endl;

  fResultSystTexTable << setw(45) << left << "Nominal Fake Rates"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kOneJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsErrSqr[kFakeElectron][kAllFinalStates][kVBFBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << setw(45) << left << "Jet pT Spectrum Systematics Up"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][kOneJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystUpErrSqr[kFakeElectron][kAllFinalStates][kVBFBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << setw(45) << left << "Jet pT Spectrum Systematics Down"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][kZeroJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][kZeroJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][kOneJetBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][kOneJetBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"%.1f +/- %.1f",FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][kVBFBin][0],TMath::Sqrt(FakeLeptonBkgYieldsSystDownErrSqr[kFakeElectron][kAllFinalStates][kVBFBin][0]));
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;

  fResultSystTexTable << setw(45) << left << "Systematic Uncertainty"
                  << setw(3) << left << "&";
  
  sprintf(buffer,"+%.0f - %.0f", 
          100*fabs(FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][kZeroJetBin][0] - 
               FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0])
          /FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0],
          100*fabs(FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][kZeroJetBin][0] - 
               FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0])
          /FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kZeroJetBin][0]
    );
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"+%.0f - %.0f", 
          100*fabs(FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][kOneJetBin][0] - 
               FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kOneJetBin][0])
          /FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kOneJetBin][0],
          100*fabs(FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][kOneJetBin][0] - 
               FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kOneJetBin][0])
          /FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kOneJetBin][0]
    );
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << setw(3) << left << "&" ;

  sprintf(buffer,"+%.0f - %.0f", 
          100*fabs(FakeLeptonBkgYieldsSystUp[kFakeElectron][kAllFinalStates][kVBFBin][0] - 
               FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kVBFBin][0])
          /FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kVBFBin][0],
          100*fabs(FakeLeptonBkgYieldsSystDown[kFakeElectron][kAllFinalStates][kVBFBin][0] - 
               FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kVBFBin][0])
          /FakeLeptonBkgYields[kFakeElectron][kAllFinalStates][kVBFBin][0]
    );
  fResultSystTexTable << setw(20) << left << buffer;
  fResultSystTexTable << " \\\\"
                  << endl;
}
