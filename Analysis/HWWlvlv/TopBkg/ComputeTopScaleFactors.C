#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/HWWCuts.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"

const int verboseLevel =   1;

//------------------------------------------------------------------------------
// dataEstimations
//------------------------------------------------------------------------------
void ComputeTopScaleFactors
(
 Int_t period = 0,
 TString bgdInputFile    = "",
 TString dataInputFile   = ""
 )
{

  //*******************************************************************************
  //Settings 
  //*******************************************************************************
  bool WWXSSel = false;
  double ptLepMin = 10.0;
  if(WWXSSel == true) ptLepMin = 20.;
  Bool_t useDYMVA = false;

  const Int_t nmass = 32;
  const Double_t mH[nmass] = {  0,  0,110,115,118,120,122,124,125,126,
                              128,130,135,140,145,150,155,160,170,180,
			      190,200,250,300,350,400,450,500,550,600,
			      700,800};  
  int nVBFLoop = nmass;
  if(WWXSSel == true) {nVBFLoop = 1;}

  double lumi = 1;
  enum { kOther, kTTBAR, kTW, kData };

  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  if	 (period == 0){ // Full2011-Fall11-V9-3500ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    lumi     = 3.553;minRun =      0;maxRun = 999999;
    bgdInputFile  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_3500ipb/backgroundA_skim2.root";
    dataInputFile = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_3500ipb/data_skim2.root";
  }
  else if(period == -1){ // Full2011-Fall11-V9 NoJetId
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    lumi     = 3.553;minRun =      0;maxRun = 999999;
    bgdInputFile  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_NoJetId/backgroundA_skim2.root";
    dataInputFile = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_NoJetId/data_skim2.root";
  }
  else if(period == 1){ // Full2011-Fall11-V9-5000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-06_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_5000ipb_71mb.root";
    lumi     = 5.098;minRun =      0;maxRun = 999999;
    bgdInputFile  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/backgroundA_skim2.root";
    dataInputFile = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/data_skim2.root";
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  //*******************************************************************************
  //Load Auxiliary Inputs
  //*******************************************************************************
  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;
  LeptonScaleLookup trigLookup(effPath.Data());

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
 
  TFile *fPUS4File = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

  //*******************************************************************************
  //Yields and Histograms
  //*******************************************************************************
  vector<vector<double> > btag_central_2j_den,btag_central_2j_num,btag_central_2j_den_error,btag_central_2j_num_error;
  vector<double>          btag_central_All_2j_den,btag_central_All_2j_num,btag_central_All_2j_den_error,btag_central_All_2j_num_error;
  vector<vector<double> > btag_lowpt_1j_den;
  vector<vector<double> > btag_lowpt_1j_num;
  vector<vector<double> > btag_lowpt_0j_den;
  vector<vector<double> > btag_lowpt_0j_num;
  vector<vector<double> > btag_highestpt_2j_den,btag_highestpt_2j_num,btag_highestpt_1j_den,btag_highestpt_1j_num;
  vector<vector<double> > btag_lowpt_1j_den_error,btag_lowpt_1j_num_error,btag_lowpt_0j_den_error,btag_lowpt_0j_num_error;
  vector<vector<double> > btag_highestpt_2j_den_error,btag_highestpt_2j_num_error,btag_highestpt_1j_den_error,btag_highestpt_1j_num_error;
  vector<vector<double> > btag_vbf_2j_den0,btag_vbf_2j_num0,btag_vbf_2j_den0_error,btag_vbf_2j_num0_error;
  vector<vector<double> > btag_vbf_2j_den1,btag_vbf_2j_num1,btag_vbf_2j_den1_error,btag_vbf_2j_num1_error;
  vector<vector<double> > btag_vbf_2j_den2,btag_vbf_2j_num2,btag_vbf_2j_den2_error,btag_vbf_2j_num2_error;
  vector<vector<double> > btag_vbf_2j_den3,btag_vbf_2j_num3,btag_vbf_2j_den3_error,btag_vbf_2j_num3_error;
  vector<vector<double> > btag_vbfSF_2j_den0,btag_vbfSF_2j_num0,btag_vbfSF_2j_den0_error,btag_vbfSF_2j_num0_error;
  vector<vector<double> > btag_vbfSF_2j_den1,btag_vbfSF_2j_num1,btag_vbfSF_2j_den1_error,btag_vbfSF_2j_num1_error;
  vector<vector<double> > btag_vbfSF_2j_den2,btag_vbfSF_2j_num2,btag_vbfSF_2j_den2_error,btag_vbfSF_2j_num2_error;
  vector<vector<double> > btag_vbfSF_2j_den3,btag_vbfSF_2j_num3,btag_vbfSF_2j_den3_error,btag_vbfSF_2j_num3_error;
  vector<vector<double> > btag_vbfOF_2j_den0,btag_vbfOF_2j_num0,btag_vbfOF_2j_den0_error,btag_vbfOF_2j_num0_error;
  vector<vector<double> > btag_vbfOF_2j_den1,btag_vbfOF_2j_num1,btag_vbfOF_2j_den1_error,btag_vbfOF_2j_num1_error;
  vector<vector<double> > btag_vbfOF_2j_den2,btag_vbfOF_2j_num2,btag_vbfOF_2j_den2_error,btag_vbfOF_2j_num2_error;
  vector<vector<double> > btag_vbfOF_2j_den3,btag_vbfOF_2j_num3,btag_vbfOF_2j_den3_error,btag_vbfOF_2j_num3_error;

  for(int j=0; j<5; j++){
    btag_central_All_2j_den.push_back(0),btag_central_All_2j_num.push_back(0),btag_central_All_2j_den_error.push_back(0),btag_central_All_2j_num_error.push_back(0);
  }

  for(int i=0; i<nmass; i++){
    vector<double> tmpbtag_vbf_2j_den0,tmpbtag_vbf_2j_num0,tmpbtag_vbf_2j_den0_error,tmpbtag_vbf_2j_num0_error;
    vector<double> tmpbtag_vbf_2j_den1,tmpbtag_vbf_2j_num1,tmpbtag_vbf_2j_den1_error,tmpbtag_vbf_2j_num1_error;
    vector<double> tmpbtag_vbf_2j_den2,tmpbtag_vbf_2j_num2,tmpbtag_vbf_2j_den2_error,tmpbtag_vbf_2j_num2_error;
    vector<double> tmpbtag_vbf_2j_den3,tmpbtag_vbf_2j_num3,tmpbtag_vbf_2j_den3_error,tmpbtag_vbf_2j_num3_error;
    vector<double> tmpbtag_vbfSF_2j_den0,tmpbtag_vbfSF_2j_num0,tmpbtag_vbfSF_2j_den0_error,tmpbtag_vbfSF_2j_num0_error;
    vector<double> tmpbtag_vbfSF_2j_den1,tmpbtag_vbfSF_2j_num1,tmpbtag_vbfSF_2j_den1_error,tmpbtag_vbfSF_2j_num1_error;
    vector<double> tmpbtag_vbfSF_2j_den2,tmpbtag_vbfSF_2j_num2,tmpbtag_vbfSF_2j_den2_error,tmpbtag_vbfSF_2j_num2_error;
    vector<double> tmpbtag_vbfSF_2j_den3,tmpbtag_vbfSF_2j_num3,tmpbtag_vbfSF_2j_den3_error,tmpbtag_vbfSF_2j_num3_error;
    vector<double> tmpbtag_vbfOF_2j_den0,tmpbtag_vbfOF_2j_num0,tmpbtag_vbfOF_2j_den0_error,tmpbtag_vbfOF_2j_num0_error;
    vector<double> tmpbtag_vbfOF_2j_den1,tmpbtag_vbfOF_2j_num1,tmpbtag_vbfOF_2j_den1_error,tmpbtag_vbfOF_2j_num1_error;
    vector<double> tmpbtag_vbfOF_2j_den2,tmpbtag_vbfOF_2j_num2,tmpbtag_vbfOF_2j_den2_error,tmpbtag_vbfOF_2j_num2_error;
    vector<double> tmpbtag_vbfOF_2j_den3,tmpbtag_vbfOF_2j_num3,tmpbtag_vbfOF_2j_den3_error,tmpbtag_vbfOF_2j_num3_error;
    for(int j=0; j<5; j++){
      tmpbtag_vbf_2j_den0.push_back(0),tmpbtag_vbf_2j_num0.push_back(0),tmpbtag_vbf_2j_den0_error.push_back(0),tmpbtag_vbf_2j_num0_error.push_back(0);
      tmpbtag_vbf_2j_den1.push_back(0),tmpbtag_vbf_2j_num1.push_back(0),tmpbtag_vbf_2j_den1_error.push_back(0),tmpbtag_vbf_2j_num1_error.push_back(0);
      tmpbtag_vbf_2j_den2.push_back(0),tmpbtag_vbf_2j_num2.push_back(0),tmpbtag_vbf_2j_den2_error.push_back(0),tmpbtag_vbf_2j_num2_error.push_back(0);
      tmpbtag_vbf_2j_den3.push_back(0),tmpbtag_vbf_2j_num3.push_back(0),tmpbtag_vbf_2j_den3_error.push_back(0),tmpbtag_vbf_2j_num3_error.push_back(0);
      tmpbtag_vbfSF_2j_den0.push_back(0),tmpbtag_vbfSF_2j_num0.push_back(0),tmpbtag_vbfSF_2j_den0_error.push_back(0),tmpbtag_vbfSF_2j_num0_error.push_back(0);
      tmpbtag_vbfSF_2j_den1.push_back(0),tmpbtag_vbfSF_2j_num1.push_back(0),tmpbtag_vbfSF_2j_den1_error.push_back(0),tmpbtag_vbfSF_2j_num1_error.push_back(0);
      tmpbtag_vbfSF_2j_den2.push_back(0),tmpbtag_vbfSF_2j_num2.push_back(0),tmpbtag_vbfSF_2j_den2_error.push_back(0),tmpbtag_vbfSF_2j_num2_error.push_back(0);
      tmpbtag_vbfSF_2j_den3.push_back(0),tmpbtag_vbfSF_2j_num3.push_back(0),tmpbtag_vbfSF_2j_den3_error.push_back(0),tmpbtag_vbfSF_2j_num3_error.push_back(0);
      tmpbtag_vbfOF_2j_den0.push_back(0),tmpbtag_vbfOF_2j_num0.push_back(0),tmpbtag_vbfOF_2j_den0_error.push_back(0),tmpbtag_vbfOF_2j_num0_error.push_back(0);
      tmpbtag_vbfOF_2j_den1.push_back(0),tmpbtag_vbfOF_2j_num1.push_back(0),tmpbtag_vbfOF_2j_den1_error.push_back(0),tmpbtag_vbfOF_2j_num1_error.push_back(0);
      tmpbtag_vbfOF_2j_den2.push_back(0),tmpbtag_vbfOF_2j_num2.push_back(0),tmpbtag_vbfOF_2j_den2_error.push_back(0),tmpbtag_vbfOF_2j_num2_error.push_back(0);
      tmpbtag_vbfOF_2j_den3.push_back(0),tmpbtag_vbfOF_2j_num3.push_back(0),tmpbtag_vbfOF_2j_den3_error.push_back(0),tmpbtag_vbfOF_2j_num3_error.push_back(0);
    }

    btag_vbf_2j_den0.push_back(tmpbtag_vbf_2j_den0),
      btag_vbf_2j_num0.push_back(tmpbtag_vbf_2j_num0);
    btag_vbf_2j_den0_error.push_back(tmpbtag_vbf_2j_den0_error),
      btag_vbf_2j_num0_error.push_back(tmpbtag_vbf_2j_num0_error);
    btag_vbf_2j_den1.push_back(tmpbtag_vbf_2j_den1),
      btag_vbf_2j_num1.push_back(tmpbtag_vbf_2j_num1);
    btag_vbf_2j_den1_error.push_back(tmpbtag_vbf_2j_den1_error),
      btag_vbf_2j_num1_error.push_back(tmpbtag_vbf_2j_num1_error);
    btag_vbf_2j_den2.push_back(tmpbtag_vbf_2j_den2),
      btag_vbf_2j_num2.push_back(tmpbtag_vbf_2j_num2);
    btag_vbf_2j_den2_error.push_back(tmpbtag_vbf_2j_den2_error),
      btag_vbf_2j_num2_error.push_back(tmpbtag_vbf_2j_num2_error);
    btag_vbf_2j_den3.push_back(tmpbtag_vbf_2j_den3),
      btag_vbf_2j_num3.push_back(tmpbtag_vbf_2j_num3);
    btag_vbf_2j_den3_error.push_back(tmpbtag_vbf_2j_den3_error),
      btag_vbf_2j_num3_error.push_back(tmpbtag_vbf_2j_num3_error);
    btag_vbfSF_2j_den0.push_back(tmpbtag_vbfSF_2j_den0),
      btag_vbfSF_2j_num0.push_back(tmpbtag_vbfSF_2j_num0);
    btag_vbfSF_2j_den0_error.push_back(tmpbtag_vbfSF_2j_den0_error),
      btag_vbfSF_2j_num0_error.push_back(tmpbtag_vbfSF_2j_num0_error);
    btag_vbfSF_2j_den1.push_back(tmpbtag_vbfSF_2j_den1),
      btag_vbfSF_2j_num1.push_back(tmpbtag_vbfSF_2j_num1);
    btag_vbfSF_2j_den1_error.push_back(tmpbtag_vbfSF_2j_den1_error),
      btag_vbfSF_2j_num1_error.push_back(tmpbtag_vbfSF_2j_num1_error);
    btag_vbfSF_2j_den2.push_back(tmpbtag_vbfSF_2j_den2),
      btag_vbfSF_2j_num2.push_back(tmpbtag_vbfSF_2j_num2);
    btag_vbfSF_2j_den2_error.push_back(tmpbtag_vbfSF_2j_den2_error),
      btag_vbfSF_2j_num2_error.push_back(tmpbtag_vbfSF_2j_num2_error);
    btag_vbfSF_2j_den3.push_back(tmpbtag_vbfSF_2j_den3),
      btag_vbfSF_2j_num3.push_back(tmpbtag_vbfSF_2j_num3);
    btag_vbfSF_2j_den3_error.push_back(tmpbtag_vbfSF_2j_den3_error),
      btag_vbfSF_2j_num3_error.push_back(tmpbtag_vbfSF_2j_num3_error);
    btag_vbfOF_2j_den0.push_back(tmpbtag_vbfOF_2j_den0),
      btag_vbfOF_2j_num0.push_back(tmpbtag_vbfOF_2j_num0);
    btag_vbfOF_2j_den0_error.push_back(tmpbtag_vbfOF_2j_den0_error),
      btag_vbfOF_2j_num0_error.push_back(tmpbtag_vbfOF_2j_num0_error);
    btag_vbfOF_2j_den1.push_back(tmpbtag_vbfOF_2j_den1),
      btag_vbfOF_2j_num1.push_back(tmpbtag_vbfOF_2j_num1);
    btag_vbfOF_2j_den1_error.push_back(tmpbtag_vbfOF_2j_den1_error),
      btag_vbfOF_2j_num1_error.push_back(tmpbtag_vbfOF_2j_num1_error);
    btag_vbfOF_2j_den2.push_back(tmpbtag_vbfOF_2j_den2),
      btag_vbfOF_2j_num2.push_back(tmpbtag_vbfOF_2j_num2);
    btag_vbfOF_2j_den2_error.push_back(tmpbtag_vbfOF_2j_den2_error),
      btag_vbfOF_2j_num2_error.push_back(tmpbtag_vbfOF_2j_num2_error);
    btag_vbfOF_2j_den3.push_back(tmpbtag_vbfOF_2j_den3),
      btag_vbfOF_2j_num3.push_back(tmpbtag_vbfOF_2j_num3);
    btag_vbfOF_2j_den3_error.push_back(tmpbtag_vbfOF_2j_den3_error),
      btag_vbfOF_2j_num3_error.push_back(tmpbtag_vbfOF_2j_num3_error);
  }

  for(int i=0; i<4; i++){
    vector<double> tmpbtag_central_2j_den,tmpbtag_central_2j_num,tmpbtag_central_2j_den_error,tmpbtag_central_2j_num_error;
    vector<double> tmpbtag_central_All_2j_den,tmpbtag_central_All_2j_num,tmpbtag_central_All_2j_den_error,tmpbtag_central_All_2j_num_error;
    vector<double> tmpbtag_lowpt_1j_den,tmpbtag_lowpt_1j_num,tmpbtag_lowpt_0j_den,tmpbtag_lowpt_0j_num;
    vector<double> tmpbtag_highestpt_2j_den,tmpbtag_highestpt_2j_num,tmpbtag_highestpt_1j_den,tmpbtag_highestpt_1j_num;
    vector<double> tmpbtag_lowpt_1j_den_error,tmpbtag_lowpt_1j_num_error,tmpbtag_lowpt_0j_den_error,tmpbtag_lowpt_0j_num_error;
    vector<double> tmpbtag_highestpt_2j_den_error,tmpbtag_highestpt_2j_num_error,tmpbtag_highestpt_1j_den_error,tmpbtag_highestpt_1j_num_error;

    for(int j=0; j<5; j++){
      tmpbtag_central_2j_den.push_back(0),tmpbtag_central_2j_num.push_back(0),tmpbtag_central_2j_den_error.push_back(0),tmpbtag_central_2j_num_error.push_back(0);
      tmpbtag_lowpt_1j_den.push_back(0),tmpbtag_lowpt_1j_num.push_back(0),tmpbtag_lowpt_0j_den.push_back(0),tmpbtag_lowpt_0j_num.push_back(0);
      tmpbtag_highestpt_2j_den.push_back(0),tmpbtag_highestpt_2j_num.push_back(0),tmpbtag_highestpt_1j_den.push_back(0),tmpbtag_highestpt_1j_num.push_back(0);
      tmpbtag_lowpt_1j_den_error.push_back(0),tmpbtag_lowpt_1j_num_error.push_back(0),tmpbtag_lowpt_0j_den_error.push_back(0),tmpbtag_lowpt_0j_num_error.push_back(0);
      tmpbtag_highestpt_2j_den_error.push_back(0),tmpbtag_highestpt_2j_num_error.push_back(0),tmpbtag_highestpt_1j_den_error.push_back(0),tmpbtag_highestpt_1j_num_error.push_back(0);
    }
    btag_central_2j_den.push_back(tmpbtag_central_2j_den),
      btag_central_2j_num.push_back(tmpbtag_central_2j_num);
    btag_central_2j_den_error.push_back(tmpbtag_central_2j_den_error),
      btag_central_2j_num_error.push_back(tmpbtag_central_2j_num_error);

    btag_lowpt_1j_den.push_back(tmpbtag_lowpt_1j_den),
      btag_lowpt_1j_num.push_back(tmpbtag_lowpt_1j_num),
      btag_lowpt_0j_den.push_back(tmpbtag_lowpt_0j_den),
      btag_lowpt_0j_num.push_back(tmpbtag_lowpt_0j_num);
    
    btag_highestpt_2j_den.push_back(tmpbtag_highestpt_2j_den),
      btag_highestpt_2j_num.push_back(tmpbtag_highestpt_2j_num),
      btag_highestpt_1j_den.push_back(tmpbtag_highestpt_1j_den),
      btag_highestpt_1j_num.push_back(tmpbtag_highestpt_1j_num);
    
    btag_lowpt_1j_den_error.push_back(tmpbtag_lowpt_1j_den_error),
      btag_lowpt_1j_num_error.push_back(tmpbtag_lowpt_1j_num_error),
      btag_lowpt_0j_den_error.push_back(tmpbtag_lowpt_0j_den_error),
      btag_lowpt_0j_num_error.push_back(tmpbtag_lowpt_0j_num_error);
    
    btag_highestpt_2j_den_error.push_back(tmpbtag_highestpt_2j_den_error),
      btag_highestpt_2j_num_error.push_back(tmpbtag_highestpt_2j_num_error),
      btag_highestpt_1j_den_error.push_back(tmpbtag_highestpt_1j_den_error),
      btag_highestpt_1j_num_error.push_back(tmpbtag_highestpt_1j_num_error);       
  }

  unsigned int patternTopTagNotInJets = SmurfTree::TopTagNotInJets;

  //*******************************************************************************
  //Background Events
  //*******************************************************************************
  Float_t dymva = -100.0;
  bgdEvent.tree_->SetBranchAddress(Form("dymva"), &dymva);
  int nBgd=bgdEvent.tree_->GetEntries();
  for (int i=0; i<nBgd; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.tree_->GetEntry(i);
    if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue;
    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::wjets           ) fDecay = 3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar           ) fDecay = 5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee            ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm            ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt            ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::qcd             ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw              ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::qqww            ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz              ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz              ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww            ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma          ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::data            ) fDecay =  1;
    else {fDecay = 0;cout << bgdEvent.dstype_ << endl;assert(0);}

    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);

    double minmet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool passMET = minmet > 20.;
    if(useDYMVA == false){
      if     (bgdEvent.njets_ == 0) passMET = passMET && (minmet > 45. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      else if(bgdEvent.njets_ == 1) passMET = passMET && (minmet > 45. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (bgdEvent.met_ > 45.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    } else {
      if     (bgdEvent.njets_ == 0) passMET = passMET && (dymva >  0.60 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      else if(bgdEvent.njets_ == 1) passMET = passMET && (dymva >  0.30 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (bgdEvent.met_ > 45.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }

    bool passNewCuts = bgdEvent.dilep_.Pt() > 45;

    // begin computing weights
    double theWeight = 0.0;
    double add       = 1.0;
    int nFake = 0;
    if(bgdEvent.dstype_ == SmurfTree::data ){
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    } else {
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    }
    bool isRealLepton = false;
    if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
       (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;
    if(nFake > 1){
      theWeight = 0.0;
    }
    else if(nFake == 1){
      if(bgdEvent.dstype_ == SmurfTree::data){
    	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        										(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        										(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	theWeight	       = add;
      }
      else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){
        add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	add = add*nPUScaleFactor2012(fhDPUS4,bgdEvent.npu_);
        add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	        					         fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
							         TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        add = add*trigEff;
	fDecay  	       = 3;
	theWeight	       = -1.0 * bgdEvent.scale1fb_*lumi*add;
      }
    }
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven || bgdEvent.dstype_ == SmurfTree::qcd) {
      theWeight = ZttScaleFactor(bgdEvent.nvtx_,period,bgdEvent.scale1fb_)*lumi;
    }
    else if(bgdEvent.dstype_ != SmurfTree::data){
      double add1 = nPUScaleFactor2012(fhDPUS4,bgdEvent.npu_);

      double add2 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
      add2   = add2*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
      double add3 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
        					            fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        					            TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      add = add1*add2*add3;

      if(fDecay == 9 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) {
        //For DY Background
        if(bgdEvent.njets_ == 0) add=add*DYBkgScaleFactor(0,0);
        if(bgdEvent.njets_ == 1) add=add*DYBkgScaleFactor(0,1);
        if(bgdEvent.njets_ >= 2) add=add*DYBkgScaleFactor(0,2);
      }
      if(fDecay == 3)  add=add*WJetsMCScaleFactor();
      if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor();
      theWeight = bgdEvent.scale1fb_*lumi*add;
    }

    bool dPhiDiLepJetCut = true;
    if(useDYMVA == false){
      if(bgdEvent.njets_ <= 1) dPhiDiLepJetCut = bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
	                                         bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
      else                     dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
	                                                   bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
    }
    if(bgdEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
	                                                 bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;

    if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
       bgdEvent.dstype_ != SmurfTree::data) &&
      bgdEvent.dilep_.M()   > 12 &&
      (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
      charge == 0 &&
      bgdEvent.lep1_.Pt() > 20. &&
      bgdEvent.lep2_.Pt() > ptLepMin &&
      passMET == true &&
      passNewCuts == true &&
      (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
      dPhiDiLepJetCut == true &&
      1 == 1
      ){
      int classType = kOther;
      if(fDecay == 5 ) classType = kTTBAR;
      if(fDecay == 13) classType = kTW;


      if(bgdEvent.jet1Btag_ >= 2.10 && bgdEvent.njets_ == 1){
        btag_lowpt_1j_den[classType][4]		 += theWeight;
        btag_lowpt_1j_den[classType][bgdEvent.type_] += theWeight;
        btag_lowpt_1j_den_error[classType][4]	       += theWeight*theWeight;
        btag_lowpt_1j_den_error[classType][bgdEvent.type_] += theWeight*theWeight;
        if((bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
          btag_lowpt_1j_num[classType][4]		   += theWeight;
          btag_lowpt_1j_num[classType][bgdEvent.type_] += theWeight;
          btag_lowpt_1j_num_error[classType][4]		 += theWeight*theWeight;
          btag_lowpt_1j_num_error[classType][bgdEvent.type_] += theWeight*theWeight;
        }
      }

     if(bgdEvent.njets_ == 0){
        btag_lowpt_0j_den[classType][4]		 += theWeight;
        btag_lowpt_0j_den[classType][bgdEvent.type_] += theWeight;
        btag_lowpt_0j_den_error[classType][4]	       += theWeight*theWeight;
        btag_lowpt_0j_den_error[classType][bgdEvent.type_] += theWeight*theWeight;
        if((bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
          btag_lowpt_0j_num[classType][4]		   += theWeight;
          btag_lowpt_0j_num[classType][bgdEvent.type_] += theWeight;
          btag_lowpt_0j_num_error[classType][4]		 += theWeight*theWeight;
          btag_lowpt_0j_num_error[classType][bgdEvent.type_] += theWeight*theWeight;
        }
      }

      if((bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && bgdEvent.jet2Btag_ >= 2.10 && bgdEvent.njets_ == 2){
        btag_highestpt_2j_den[classType][4] 	     += theWeight;
        btag_highestpt_2j_den[classType][bgdEvent.type_] += theWeight;
        btag_highestpt_2j_den_error[classType][4] 	           += theWeight*theWeight;
        btag_highestpt_2j_den_error[classType][bgdEvent.type_] += theWeight*theWeight;
        if(bgdEvent.jet1Btag_ >= 2.10){
          btag_highestpt_2j_num[classType][4]	       += theWeight;
          btag_highestpt_2j_num[classType][bgdEvent.type_] += theWeight;
          btag_highestpt_2j_num_error[classType][4]	             += theWeight*theWeight;
          btag_highestpt_2j_num_error[classType][bgdEvent.type_] += theWeight*theWeight;
        }
      }

      if((bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && bgdEvent.njets_ == 1){
        btag_highestpt_1j_den[classType][4] 	     += theWeight;
        btag_highestpt_1j_den[classType][bgdEvent.type_] += theWeight;
        btag_highestpt_1j_den_error[classType][4] 	           += theWeight*theWeight;
        btag_highestpt_1j_den_error[classType][bgdEvent.type_] += theWeight*theWeight;
        if(bgdEvent.jet1Btag_ >= 2.10){
          btag_highestpt_1j_num[classType][4]	       += theWeight;
          btag_highestpt_1j_num[classType][bgdEvent.type_] += theWeight;
          btag_highestpt_1j_num_error[classType][4]	             += theWeight*theWeight;
          btag_highestpt_1j_num_error[classType][bgdEvent.type_] += theWeight*theWeight;
        }
      }

      Double_t etaMin = TMath::Abs(bgdEvent.jet1_.Eta()); Double_t bTagMax[2] = {bgdEvent.jet1Btag_, bgdEvent.jet2Btag_};
      if(etaMin > TMath::Abs(bgdEvent.jet2_.Eta())) {
        etaMin = TMath::Abs(bgdEvent.jet2_.Eta());
	bTagMax[0] = bgdEvent.jet2Btag_; bTagMax[1] = bgdEvent.jet1Btag_;
      }

      unsigned int Njet3 = bgdEvent.njets_;
      if(bgdEvent.jet3_.pt() <= 30)									     Njet3 = 2;
      else if(bgdEvent.jet3_.pt() > 30 && (
        (bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() < 0) ||
        (bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() < 0)))   Njet3 = 0;
      else												     Njet3 = 2;
      if(bgdEvent.njets_ < 2 || bgdEvent.njets_ > 3)							     Njet3 = 0;

      if((bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && bgdEvent.nSoftMuons_ == 0 && bTagMax[1] < 2.10 && 
          bgdEvent.jet3Btag_ < 2.1 && Njet3 == 2 && (bgdEvent.jet1_+bgdEvent.jet2_).M() > 0){
        int nEta = TMath::Min(etaMin,2.499)/2.5*5;
        btag_central_2j_den[classType][nEta] 	     += theWeight;
        btag_central_2j_den_error[classType][nEta]   += theWeight*theWeight;
        if(bTagMax[0] >= 2.10){
          btag_central_2j_num[classType][nEta]	     += theWeight;
          btag_central_2j_num_error[classType][nEta] += theWeight*theWeight;
        }
        btag_central_All_2j_den[classType]	     += theWeight;
        btag_central_All_2j_den_error[classType]     += theWeight*theWeight;
        if(bTagMax[0] >= 2.10){
          btag_central_All_2j_num[classType]         += theWeight;
          btag_central_All_2j_num_error[classType]   += theWeight*theWeight;
        }
      }

      int centrality = 0;
      if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
          (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
         ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
          (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
      bool passWBFSel = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 450. && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 && 
                         centrality == 1;
      if(WWXSSel == true) {passWBFSel = true;}

      if((bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && bgdEvent.nSoftMuons_ == 0 && bTagMax[1] < 2.10 && bgdEvent.jet3Btag_ < 2.1 && Njet3 == 2){
        int nEta = TMath::Min(etaMin,2.499)/2.5*5;
        for(int ivbf=0; ivbf<nVBFLoop; ivbf++){
	  bool passAllCuts = true;
	  if     (ivbf==1){
	    passAllCuts = passWBFSel == true;
	  }
	  else if(ivbf>=2){
	    passAllCuts = bgdEvent.dilep_.M() <  DileptonMassPreselectionCut(mH[ivbf]) && 
	                  bgdEvent.mt_ > 30.0 && bgdEvent.mt_ < mH[ivbf] && 
			  passWBFSel == true;
	  }
	  if(passAllCuts == true){
            if(bTagMax[0] >= 2.10){
	      if       (classType == kOther) {
        	btag_vbf_2j_num0[ivbf][nEta]	   += theWeight;
        	btag_vbf_2j_num0_error[ivbf][nEta] += theWeight*theWeight;
		if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_num0[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_num0_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_num0[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_num0_error[ivbf][nEta] += theWeight*theWeight;		  
		}
	      } else if(classType == kTTBAR) {
        	btag_vbf_2j_num1[ivbf][nEta]	   += theWeight;
        	btag_vbf_2j_num1_error[ivbf][nEta] += theWeight*theWeight;
		if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_num1[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_num1_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_num1[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_num1_error[ivbf][nEta] += theWeight*theWeight;		  
		}
	      } else if(classType == kTW) {
        	btag_vbf_2j_num2[ivbf][nEta]	   += theWeight;
        	btag_vbf_2j_num2_error[ivbf][nEta] += theWeight*theWeight;
		if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_num2[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_num2_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_num2[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_num2_error[ivbf][nEta] += theWeight*theWeight;		  
		}
	      }
            } else {
	      if       (classType == kOther) {
        	btag_vbf_2j_den0[ivbf][nEta]	   += theWeight;
        	btag_vbf_2j_den0_error[ivbf][nEta] += theWeight*theWeight;
		if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_den0[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_den0_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_den0[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_den0_error[ivbf][nEta] += theWeight*theWeight;		  
		}
	      } else if(classType == kTTBAR) {
        	btag_vbf_2j_den1[ivbf][nEta]	   += theWeight;
        	btag_vbf_2j_den1_error[ivbf][nEta] += theWeight*theWeight;
		if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_den1[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_den1_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_den1[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_den1_error[ivbf][nEta] += theWeight*theWeight;		  
		}
	      } else if(classType == kTW) {
        	btag_vbf_2j_den2[ivbf][nEta]	   += theWeight;
        	btag_vbf_2j_den2_error[ivbf][nEta] += theWeight*theWeight;
		if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_den2[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_den2_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_den2[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_den2_error[ivbf][nEta] += theWeight*theWeight;		  
		}
	      }
	    }
	  } // passAllCuts
	} // Loop over masses
      } // minimum preselection
    }
  } // end background loop

  //*******************************************************************************
  //Data Events
  //*******************************************************************************
  dymva = -100.0;
  dataEvent.tree_->SetBranchAddress(Form("dymva"), &dymva);
  int nData=dataEvent.tree_->GetEntries();
  for (int i=0; i<nData; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nData);
    dataEvent.tree_->GetEntry(i);

    if(!((dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)) continue;

    if((dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue;
    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    double minmet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
    bool passMET = minmet > 20.;
    if(useDYMVA == false){
      if     (dataEvent.njets_ == 0) passMET = passMET && (minmet > 45. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else if(dataEvent.njets_ == 1) passMET = passMET && (minmet > 45. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (dataEvent.met_ > 45.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    } else {
      if     (dataEvent.njets_ == 0) passMET = passMET && (dymva >  0.60 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else if(dataEvent.njets_ == 1) passMET = passMET && (dymva >  0.30 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (dataEvent.met_ > 45.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }

    bool passNewCuts = dataEvent.dilep_.Pt() > 45;

    bool dPhiDiLepJetCut = true;
    if(useDYMVA == false){
      if(dataEvent.njets_ <= 1) dPhiDiLepJetCut = dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
	                                          dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
      else                     dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
	                                                   dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    }
    if(dataEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
	                                                  dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;

    double theWeight = 1.0;
    if(
      dataEvent.dilep_.M()   > 12 &&
      (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
      charge == 0 &&
      dataEvent.lep1_.Pt() > 20. &&
      dataEvent.lep2_.Pt() > ptLepMin &&
      passMET == true &&
      passNewCuts == true &&
      (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
      dPhiDiLepJetCut == true &&
      1 == 1
      ) {
      int classType = kData;
      if(dataEvent.jet1Btag_ >= 2.10 && dataEvent.njets_ == 1){
        btag_lowpt_1j_den[classType][4]		  += theWeight;
        btag_lowpt_1j_den[classType][dataEvent.type_] += theWeight;
        btag_lowpt_1j_den_error[classType][4]	        += theWeight*theWeight;
        btag_lowpt_1j_den_error[classType][dataEvent.type_] += theWeight*theWeight;
        if((dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
          btag_lowpt_1j_num[classType][4]		    += theWeight;
          btag_lowpt_1j_num[classType][dataEvent.type_] += theWeight;
          btag_lowpt_1j_num_error[classType][4]		  += theWeight*theWeight;
          btag_lowpt_1j_num_error[classType][dataEvent.type_] += theWeight*theWeight;
        }
      }

      if(dataEvent.njets_ == 0){
        btag_lowpt_0j_den[classType][4]		  += theWeight;
        btag_lowpt_0j_den[classType][dataEvent.type_] += theWeight;
        btag_lowpt_0j_den_error[classType][4]	        += theWeight*theWeight;
        btag_lowpt_0j_den_error[classType][dataEvent.type_] += theWeight*theWeight;
        if((dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
          btag_lowpt_0j_num[classType][4]		    += theWeight;
          btag_lowpt_0j_num[classType][dataEvent.type_] += theWeight;
          btag_lowpt_0j_num_error[classType][4]		  += theWeight*theWeight;
          btag_lowpt_0j_num_error[classType][dataEvent.type_] += theWeight*theWeight;
        }
      }

      if((dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && dataEvent.jet2Btag_ >= 2.10 && dataEvent.njets_ == 2){
        btag_highestpt_2j_den[classType][4] 	      += theWeight;
        btag_highestpt_2j_den[classType][dataEvent.type_] += theWeight;
        btag_highestpt_2j_den_error[classType][4] 	            += theWeight*theWeight;
        btag_highestpt_2j_den_error[classType][dataEvent.type_] += theWeight*theWeight;
        if(dataEvent.jet1Btag_ >= 2.10){
          btag_highestpt_2j_num[classType][4]	        += theWeight;
          btag_highestpt_2j_num[classType][dataEvent.type_] += theWeight;
          btag_highestpt_2j_num_error[classType][4]	              += theWeight*theWeight;
          btag_highestpt_2j_num_error[classType][dataEvent.type_] += theWeight*theWeight;
        }
      }

      if((dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && dataEvent.njets_ == 1){
        btag_highestpt_1j_den[classType][4] 	      += theWeight;
        btag_highestpt_1j_den[classType][dataEvent.type_] += theWeight;
        btag_highestpt_1j_den_error[classType][4] 	            += theWeight*theWeight;
        btag_highestpt_1j_den_error[classType][dataEvent.type_] += theWeight*theWeight;
        if(dataEvent.jet1Btag_ >= 2.10){
          btag_highestpt_1j_num[classType][4]	        += theWeight;
          btag_highestpt_1j_num[classType][dataEvent.type_] += theWeight;
          btag_highestpt_1j_num_error[classType][4]	              += theWeight*theWeight;
          btag_highestpt_1j_num_error[classType][dataEvent.type_] += theWeight*theWeight;
        }
      }

      Double_t etaMin = TMath::Abs(dataEvent.jet1_.Eta()); Double_t bTagMax[2] = {dataEvent.jet1Btag_, dataEvent.jet2Btag_};
      if(etaMin > TMath::Abs(dataEvent.jet2_.Eta())) {
        etaMin = TMath::Abs(dataEvent.jet2_.Eta());
	bTagMax[0] = dataEvent.jet2Btag_; bTagMax[1] = dataEvent.jet1Btag_;
      }

      unsigned int Njet3 = dataEvent.njets_;
      if(dataEvent.jet3_.pt() <= 30)									         Njet3 = 2;
      else if(dataEvent.jet3_.pt() > 30 && (
        (dataEvent.jet1_.eta()-dataEvent.jet3_.eta() > 0 && dataEvent.jet2_.eta()-dataEvent.jet3_.eta() < 0) ||
        (dataEvent.jet2_.eta()-dataEvent.jet3_.eta() > 0 && dataEvent.jet1_.eta()-dataEvent.jet3_.eta() < 0)))   Njet3 = 0;
      else												         Njet3 = 2;
      if(dataEvent.njets_ < 2 || dataEvent.njets_ > 3)							         Njet3 = 0;

      if((dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && dataEvent.nSoftMuons_ == 0 && bTagMax[1] < 2.10 && 
          dataEvent.jet3Btag_ < 2.1 && Njet3 == 2 && (dataEvent.jet1_+dataEvent.jet2_).M() > 0){
        int nEta = TMath::Min(etaMin,2.499)/2.5*5;
        btag_central_2j_den[classType][nEta] 	     += theWeight;
        btag_central_2j_den_error[classType][nEta]   += theWeight*theWeight;
        if(bTagMax[0] >= 2.10){
          btag_central_2j_num[classType][nEta]	     += theWeight;
          btag_central_2j_num_error[classType][nEta] += theWeight*theWeight;
        }
        btag_central_All_2j_den[classType]	     += theWeight;
        btag_central_All_2j_den_error[classType]     += theWeight*theWeight;
        if(bTagMax[0] >= 2.10){
          btag_central_All_2j_num[classType]         += theWeight;
          btag_central_All_2j_num_error[classType]   += theWeight*theWeight;
        }
      }

      int centrality = 0;
      if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
          (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
         ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
          (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1;
      bool passWBFSel = (dataEvent.jet1_+dataEvent.jet2_).M() > 450. && TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 3.5 && 
                         centrality == 1;
      if(WWXSSel == true) {passWBFSel = true;}

      if((dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && dataEvent.nSoftMuons_ == 0 && bTagMax[1] < 2.10 && dataEvent.jet3Btag_ < 2.1 && Njet3 == 2){
        int nEta = TMath::Min(etaMin,2.499)/2.5*5;
        for(int ivbf=0; ivbf<nVBFLoop; ivbf++){
	  bool passAllCuts = true;
	  if     (ivbf==1){
	    passAllCuts = passWBFSel == true;
	  }
	  else if(ivbf>=2){
	    passAllCuts = dataEvent.dilep_.M() <  DileptonMassPreselectionCut(mH[ivbf]) && 
	                  dataEvent.mt_ > 30.0 && dataEvent.mt_ < mH[ivbf] && 
			  passWBFSel == true;
	  }
	  if(passAllCuts == true){
            if(bTagMax[0] >= 2.10){
              btag_vbf_2j_num3[ivbf][nEta]	 += theWeight;
              btag_vbf_2j_num3_error[ivbf][nEta] += theWeight*theWeight;
		if(dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee){
        	  btag_vbfSF_2j_num3[ivbf][nEta]       += theWeight;
        	  btag_vbfSF_2j_num3_error[ivbf][nEta] += theWeight*theWeight;		  
		} else {
        	  btag_vbfOF_2j_num3[ivbf][nEta]       += theWeight;
        	  btag_vbfOF_2j_num3_error[ivbf][nEta] += theWeight*theWeight;		  
		}
            } else {
              btag_vbf_2j_den3[ivbf][nEta]       += theWeight;
              btag_vbf_2j_den3_error[ivbf][nEta] += theWeight*theWeight;
	      if(dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee){
                btag_vbfSF_2j_den3[ivbf][nEta]       += theWeight;
                btag_vbfSF_2j_den3_error[ivbf][nEta] += theWeight*theWeight;		
	      } else {
                btag_vbfOF_2j_den3[ivbf][nEta]       += theWeight;
                btag_vbfOF_2j_den3_error[ivbf][nEta] += theWeight*theWeight;		
	      }
	    }
	  } // passAllCuts
	} // Loop over masses
      } // minimum preselection
    }

  } // End loop data

  //*******************************************************************************
  //Print Summary 
  //*******************************************************************************
  char *classLabel[5] = {"mm ", "me ", "em " , "ee ", "all"};

  double btagSF = 1.0;

  //*******************************************************************************
  //2-Jet Bin : BTag Efficiency for central jet
  //*******************************************************************************
  printf("**********eff central jet 2-j**********\n");
  double effttMC_btag_central_2j[5],effttMC_btag_central_2j_error[5],effttMC_btag_central_tt_2j[5],effttMC_btag_central_tt_2j_error[5];
  double effttDA_btag_central_2j[5],effttDA_btag_central_2j_error[5],effttMC_btag_central_tw_2j[5],effttMC_btag_central_tw_2j_error[5];
  double TopBkgScaleFactor_2Jet_central,TopBkgScaleFactorUncertainty_2Jet_central;
  printf("channel               (data/background/top)-num             (data/background/top)-den\n");
  for(int i=0; i<5; i++) {

    //MC efficiencies
    effttMC_btag_central_tt_2j[i] = (btag_central_2j_num[1][i])/(btag_central_2j_den[1][i]);
    effttMC_btag_central_tw_2j[i] = (btag_central_2j_num[2][i])/(btag_central_2j_den[2][i]);
    effttMC_btag_central_2j[i]    = (btag_central_2j_num[1][i] + btag_central_2j_num[2][i]) / (btag_central_2j_den[1][i] + btag_central_2j_den[2][i]);
    
    effttMC_btag_central_tt_2j_error[i] = sqrt((1.0-effttMC_btag_central_tt_2j[i])*effttMC_btag_central_tt_2j[i]/(btag_central_2j_den[1][i])*
                                                 (btag_central_2j_den_error[1][i])/(btag_central_2j_den[1][i]));    
    effttMC_btag_central_tw_2j_error[i] = sqrt((1.0-effttMC_btag_central_tw_2j[i])*effttMC_btag_central_tw_2j[i]/(btag_central_2j_den[2][i])*
                                                 (btag_central_2j_den_error[2][i])/(btag_central_2j_den[2][i]));
    effttMC_btag_central_2j_error[i]    = sqrt((1.0-effttMC_btag_central_2j[i])*effttMC_btag_central_2j[i]/(btag_central_2j_den[1][i]+btag_central_2j_den[2][i])*
                                                 (btag_central_2j_den_error[1][i]+btag_central_2j_den_error[2][i])/(btag_central_2j_den[1][i]+btag_central_2j_den[2][i]));

    //Data efficiencies
    effttDA_btag_central_2j[i] = (btag_central_2j_num[3][i]-btag_central_2j_num[0][i]-btag_central_2j_num[2][i])/
                                 (btag_central_2j_den[3][i]-btag_central_2j_den[0][i]-btag_central_2j_den[2][i]);    
    //effttDA_btag_central_2j[i] = (btag_central_2j_num[3][i]-btag_central_2j_num[0][i])/
    //                             (btag_central_2j_den[3][i]-btag_central_2j_den[0][i]);
    effttDA_btag_central_2j[i] = TMath::Min(effttDA_btag_central_2j[i],0.999);
    effttDA_btag_central_2j_error[i] = sqrt((1-effttDA_btag_central_2j[i])*effttDA_btag_central_2j[i]/btag_central_2j_den[3][i]);

    printf("details_central(%d) --> %5.0f/%7.2f/%7.2f  - %5.0f/%7.2f/%7.2f\n",i,
           btag_central_2j_num[3][i],btag_central_2j_num[0][i],btag_central_2j_num[1][i],
	   btag_central_2j_den[3][i],btag_central_2j_den[0][i],btag_central_2j_den[1][i]);
  }
  printf("channel                      eff_data              eff_tt            ScaleFactor\n");
  for(int i=0; i<5; i++) {
    TopBkgScaleFactor_2Jet_central = effttDA_btag_central_2j[i]/effttMC_btag_central_tt_2j[i];
    TopBkgScaleFactorUncertainty_2Jet_central = effttDA_btag_central_2j_error[i]/effttMC_btag_central_tt_2j[i];
    printf("scaleFactor_central(%d) --> %6.3f +/- %6.3f | %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",i,
                           effttDA_btag_central_2j[i],effttDA_btag_central_2j_error[i],effttMC_btag_central_tt_2j[i],effttMC_btag_central_tt_2j_error[i],
			   TopBkgScaleFactor_2Jet_central,TopBkgScaleFactorUncertainty_2Jet_central);
  }

  // Overall
  double effttMC_btag_central_All_2j,effttMC_btag_central_All_2j_error,effttMC_btag_central_All_tt_2j,effttMC_btag_central_All_tt_2j_error;
  double effttDA_btag_central_All_2j,effttDA_btag_central_All_2j_error,effttMC_btag_central_All_tw_2j,effttMC_btag_central_All_tw_2j_error;
  effttMC_btag_central_All_tt_2j = (btag_central_All_2j_num[1])/(btag_central_All_2j_den[1]);
  effttMC_btag_central_All_tw_2j = (btag_central_All_2j_num[2])/(btag_central_All_2j_den[2]);
  effttMC_btag_central_All_2j    = (btag_central_All_2j_num[1] + btag_central_All_2j_num[2]) / (btag_central_All_2j_den[1] + btag_central_All_2j_den[2]);
  
  effttMC_btag_central_All_tt_2j_error = sqrt((1.0-effttMC_btag_central_All_tt_2j)*effttMC_btag_central_All_tt_2j/(btag_central_All_2j_den[1])*
  					       (btag_central_All_2j_den_error[1])/(btag_central_All_2j_den[1]));    
  effttMC_btag_central_All_tw_2j_error = sqrt((1.0-effttMC_btag_central_All_tw_2j)*effttMC_btag_central_All_tw_2j/(btag_central_All_2j_den[2])*
  					       (btag_central_All_2j_den_error[2])/(btag_central_All_2j_den[2]));
  effttMC_btag_central_All_2j_error    = sqrt((1.0-effttMC_btag_central_All_2j)*effttMC_btag_central_All_2j/(btag_central_All_2j_den[1]+btag_central_All_2j_den[2])*
  					       (btag_central_All_2j_den_error[1]+btag_central_All_2j_den_error[2])/(btag_central_All_2j_den[1]+btag_central_All_2j_den[2]));
  effttDA_btag_central_All_2j = (btag_central_All_2j_num[3]-btag_central_All_2j_num[0]-btag_central_All_2j_num[2])/
  			       (btag_central_All_2j_den[3]-btag_central_All_2j_den[0]-btag_central_All_2j_den[2]);    
  //effttDA_btag_central_All_2j = (btag_central_All_2j_num[3]-btag_central_All_2j_num[0])/
  //				 (btag_central_All_2j_den[3]-btag_central_All_2j_den[0]);    
  effttDA_btag_central_All_2j_error = sqrt((1-effttDA_btag_central_All_2j)*effttDA_btag_central_All_2j/btag_central_All_2j_den[3]);

  printf("details_central(all)         --> %5.0f/%7.2f/%7.2f  - %5.0f/%7.2f/%7.2f\n",
  	 btag_central_All_2j_num[3],btag_central_All_2j_num[0],btag_central_All_2j_num[1],
         btag_central_All_2j_den[3],btag_central_All_2j_den[0],btag_central_All_2j_den[1]);
  double TopBkgScaleFactor_2Jet_central_All = effttDA_btag_central_All_2j/effttMC_btag_central_All_tt_2j;
  double TopBkgScaleFactorUncertainty_2Jet_central_All = effttDA_btag_central_All_2j_error/effttMC_btag_central_All_tt_2j;
  printf("scaleFactor_central_All(all) --> %6.3f +/- %6.3f | %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",
  			 effttDA_btag_central_All_2j,effttDA_btag_central_All_2j_error,effttMC_btag_central_All_tt_2j,effttMC_btag_central_All_tt_2j_error,
        		 TopBkgScaleFactor_2Jet_central_All,TopBkgScaleFactorUncertainty_2Jet_central_All);

  //*******************************************************************************
  //2-Jet Bin : BTag Efficiency for vbf jet
  //*******************************************************************************
  printf("**********eff vbf jet 2-j**********\n");
  double evtMC_vbf_2j[5],evtMC_vbf_2j_error[5],evtDA_vbf_2j[5],evtDA_vbf_2j_error[5];
  double TopBkgScaleFactor_2Jet_vbf[nmass],TopBkgScaleFactorUncertainty_2Jet_vbf[nmass];
  for(int imass = 0; imass<nVBFLoop; imass++){
    for(int i=0; i<5; i++) {

      evtMC_vbf_2j[i] = (btag_vbf_2j_num1[imass][i]+btag_vbf_2j_num2[imass][i])*(1.0-effttMC_btag_central_tt_2j[i])/effttMC_btag_central_tt_2j[i];
      evtDA_vbf_2j[i] = TMath::Max(btag_vbf_2j_num3[imass][i]-btag_vbf_2j_num0[imass][i],0.0)*(1.0-effttDA_btag_central_2j[i])/effttDA_btag_central_2j[i];

      evtMC_vbf_2j_error[i] = (btag_vbf_2j_num1[imass][i]+btag_vbf_2j_num2[imass][i])*effttMC_btag_central_tt_2j_error[i]/effttMC_btag_central_tt_2j[i]/effttMC_btag_central_tt_2j[i];
      double errDa[2] = {sqrt(TMath::Max(btag_vbf_2j_num3[imass][i],0.0))*(1.0-effttDA_btag_central_2j[i])/effttDA_btag_central_2j[i],
                	 TMath::Max(btag_vbf_2j_num3[imass][i]-btag_vbf_2j_num0[imass][i],0.0)*effttDA_btag_central_2j_error[i]/effttDA_btag_central_2j[i]/effttDA_btag_central_2j[i]};
      evtDA_vbf_2j_error[i] = sqrt(errDa[0]*errDa[0]+errDa[1]*errDa[1]);
    }
    for(int i=0; i<5; i++) printf("eventsTaggedVBF(%d) --> top: %6.3f  data: %5d  bck:%6.3f\n",i,btag_vbf_2j_num1[imass][i]+btag_vbf_2j_num2[imass][i],
                           (int)btag_vbf_2j_num3[imass][i],btag_vbf_2j_num0[imass][i]);
    for(int i=0; i<5; i++) printf("eventsNonTaggedVBF(%d) --> %6.3f  --> %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",i,btag_vbf_2j_den1[imass][i]+btag_vbf_2j_den2[imass][i],
               evtMC_vbf_2j[i],evtMC_vbf_2j_error[i],evtDA_vbf_2j[i],evtDA_vbf_2j_error[i]);

    double NonTaggedTopMCPred = 0; for(int i=0; i<5; i++) NonTaggedTopMCPred = NonTaggedTopMCPred + evtMC_vbf_2j[i];
    double NonTaggedTopMC     = 0; for(int i=0; i<5; i++) NonTaggedTopMC     = NonTaggedTopMC     + btag_vbf_2j_den1[imass][i]+btag_vbf_2j_den2[imass][i];
    double NonTaggedTopDA     = 0; for(int i=0; i<5; i++) NonTaggedTopDA     = NonTaggedTopDA     + evtDA_vbf_2j[i];
    double NonTaggedTopDA_error = 0; for(int i=0; i<5; i++) NonTaggedTopDA_error = NonTaggedTopDA_error + evtDA_vbf_2j_error[i]*evtDA_vbf_2j_error[i];
    TopBkgScaleFactor_2Jet_vbf[imass] = NonTaggedTopDA/NonTaggedTopMC;
    TopBkgScaleFactorUncertainty_2Jet_vbf[imass] = sqrt(NonTaggedTopDA_error)/NonTaggedTopMC;
    printf("data/MC/MCPred(%6.3f/%6.3f/%6.3f) -> scaleFactorVBF(%3d): %6.3f +/- %6.3f\n",NonTaggedTopDA,NonTaggedTopMC,NonTaggedTopMCPred,(int)mH[imass],TopBkgScaleFactor_2Jet_vbf[imass],TopBkgScaleFactorUncertainty_2Jet_vbf[imass]);
  }
  printf("**********eff vbf jet 2-j SF**********\n");
  double evtMC_vbfSF_2j[5],evtMC_vbfSF_2j_error[5],evtDA_vbfSF_2j[5],evtDA_vbfSF_2j_error[5];
  double TopBkgScaleFactor_2Jet_vbfSF[nmass],TopBkgScaleFactorUncertainty_2Jet_vbfSF[nmass];
  for(int imass = 0; imass<nVBFLoop; imass++){
    for(int i=0; i<5; i++) {

      evtMC_vbfSF_2j[i] = (btag_vbfSF_2j_num1[imass][i]+btag_vbfSF_2j_num2[imass][i])*(1.0-effttMC_btag_central_tt_2j[i])/effttMC_btag_central_tt_2j[i];
      evtDA_vbfSF_2j[i] = TMath::Max(btag_vbfSF_2j_num3[imass][i]-btag_vbfSF_2j_num0[imass][i],0.0)*(1.0-effttDA_btag_central_2j[i])/effttDA_btag_central_2j[i];

      evtMC_vbfSF_2j_error[i] = (btag_vbfSF_2j_num1[imass][i]+btag_vbfSF_2j_num2[imass][i])*effttMC_btag_central_tt_2j_error[i]/effttMC_btag_central_tt_2j[i]/effttMC_btag_central_tt_2j[i];
      double errDa[2] = {sqrt(TMath::Max(btag_vbfSF_2j_num3[imass][i],0.0))*(1.0-effttDA_btag_central_2j[i])/effttDA_btag_central_2j[i],
                	 TMath::Max(btag_vbfSF_2j_num3[imass][i]-btag_vbfSF_2j_num0[imass][i],0.0)*effttDA_btag_central_2j_error[i]/effttDA_btag_central_2j[i]/effttDA_btag_central_2j[i]};
      evtDA_vbfSF_2j_error[i] = sqrt(errDa[0]*errDa[0]+errDa[1]*errDa[1]);
    }
    for(int i=0; i<5; i++) printf("eventsTaggedVBFSF(%d) --> top: %6.3f  data: %5d  bck:%6.3f\n",i,btag_vbfSF_2j_num1[imass][i]+btag_vbfSF_2j_num2[imass][i],
                           (int)btag_vbfSF_2j_num3[imass][i],btag_vbfSF_2j_num0[imass][i]);
    for(int i=0; i<5; i++) printf("eventsNonTaggedVBFSF(%d) --> %6.3f  --> %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",i,btag_vbfSF_2j_den1[imass][i]+btag_vbfSF_2j_den2[imass][i],
               evtMC_vbfSF_2j[i],evtMC_vbfSF_2j_error[i],evtDA_vbfSF_2j[i],evtDA_vbfSF_2j_error[i]);

    double NonTaggedTopMCPred = 0; for(int i=0; i<5; i++) NonTaggedTopMCPred = NonTaggedTopMCPred + evtMC_vbfSF_2j[i];
    double NonTaggedTopMC     = 0; for(int i=0; i<5; i++) NonTaggedTopMC     = NonTaggedTopMC + btag_vbfSF_2j_den1[imass][i]+btag_vbfSF_2j_den2[imass][i];
    double NonTaggedTopDA     = 0; for(int i=0; i<5; i++) NonTaggedTopDA     = NonTaggedTopDA + evtDA_vbfSF_2j[i];
    double NonTaggedTopDA_error = 0; for(int i=0; i<5; i++) NonTaggedTopDA_error = NonTaggedTopDA_error + evtDA_vbfSF_2j_error[i]*evtDA_vbfSF_2j_error[i];
    TopBkgScaleFactor_2Jet_vbfSF[imass] = NonTaggedTopDA/NonTaggedTopMC;
    TopBkgScaleFactorUncertainty_2Jet_vbfSF[imass] = sqrt(NonTaggedTopDA_error)/NonTaggedTopMC;
    printf("data/MC/MCPred(%6.3f/%6.3f/%6.3f) -> scaleFactorVBFSF(%3d): %6.3f +/- %6.3f\n",NonTaggedTopDA,NonTaggedTopMC,NonTaggedTopMCPred,(int)mH[imass],TopBkgScaleFactor_2Jet_vbfSF[imass],TopBkgScaleFactorUncertainty_2Jet_vbfSF[imass]);
  }
  printf("**********eff vbf jet 2-j OF**********\n");
  double evtMC_vbfOF_2j[5],evtMC_vbfOF_2j_error[5],evtDA_vbfOF_2j[5],evtDA_vbfOF_2j_error[5];
  double TopBkgScaleFactor_2Jet_vbfOF[nmass],TopBkgScaleFactorUncertainty_2Jet_vbfOF[nmass];
  for(int imass = 0; imass<nVBFLoop; imass++){
    for(int i=0; i<5; i++) {

      evtMC_vbfOF_2j[i] = (btag_vbfOF_2j_num1[imass][i]+btag_vbfOF_2j_num2[imass][i])*(1.0-effttMC_btag_central_tt_2j[i])/effttMC_btag_central_tt_2j[i];
      evtDA_vbfOF_2j[i] = TMath::Max(btag_vbfOF_2j_num3[imass][i]-btag_vbfOF_2j_num0[imass][i],0.0)*(1.0-effttDA_btag_central_2j[i])/effttDA_btag_central_2j[i];

      evtMC_vbfOF_2j_error[i] = (btag_vbfOF_2j_num1[imass][i]+btag_vbfOF_2j_num2[imass][i])*effttMC_btag_central_tt_2j_error[i]/effttMC_btag_central_tt_2j[i]/effttMC_btag_central_tt_2j[i];
      double errDa[2] = {sqrt(TMath::Max(btag_vbfOF_2j_num3[imass][i],0.0))*(1.0-effttDA_btag_central_2j[i])/effttDA_btag_central_2j[i],
                	 TMath::Max(btag_vbfOF_2j_num3[imass][i]-btag_vbfOF_2j_num0[imass][i],0.0)*effttDA_btag_central_2j_error[i]/effttDA_btag_central_2j[i]/effttDA_btag_central_2j[i]};
      evtDA_vbfOF_2j_error[i] = sqrt(errDa[0]*errDa[0]+errDa[1]*errDa[1]);
    }
    for(int i=0; i<5; i++) printf("eventsTaggedVBFOF(%d) --> top: %6.3f  data: %5d  bck:%6.3f\n",i,btag_vbfOF_2j_num1[imass][i]+btag_vbfOF_2j_num2[imass][i],
                           (int)btag_vbfOF_2j_num3[imass][i],btag_vbfOF_2j_num0[imass][i]);
    for(int i=0; i<5; i++) printf("eventsNonTaggedVBFOF(%d) --> %6.3f  --> %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",i,btag_vbfOF_2j_den1[imass][i]+btag_vbfOF_2j_den2[imass][i],
               evtMC_vbfOF_2j[i],evtMC_vbfOF_2j_error[i],evtDA_vbfOF_2j[i],evtDA_vbfOF_2j_error[i]);

    double NonTaggedTopMCPred = 0; for(int i=0; i<5; i++) NonTaggedTopMCPred = NonTaggedTopMCPred + evtMC_vbfOF_2j[i];
    double NonTaggedTopMC     = 0; for(int i=0; i<5; i++) NonTaggedTopMC     = NonTaggedTopMC + btag_vbfOF_2j_den1[imass][i]+btag_vbfOF_2j_den2[imass][i];
    double NonTaggedTopDA     = 0; for(int i=0; i<5; i++) NonTaggedTopDA     = NonTaggedTopDA + evtDA_vbfOF_2j[i];
    double NonTaggedTopDA_error = 0; for(int i=0; i<5; i++) NonTaggedTopDA_error = NonTaggedTopDA_error + evtDA_vbfOF_2j_error[i]*evtDA_vbfOF_2j_error[i];
    TopBkgScaleFactor_2Jet_vbfOF[imass] = NonTaggedTopDA/NonTaggedTopMC;
    TopBkgScaleFactorUncertainty_2Jet_vbfOF[imass] = sqrt(NonTaggedTopDA_error)/NonTaggedTopMC;
    printf("data/MC/MCPred(%6.3f/%6.3f/%6.3f) -> scaleFactorVBFOF(%3d): %6.3f +/- %6.3f\n",NonTaggedTopDA,NonTaggedTopMC,NonTaggedTopMCPred,(int)mH[imass],TopBkgScaleFactor_2Jet_vbfOF[imass],TopBkgScaleFactorUncertainty_2Jet_vbfOF[imass]);
  }
  //*******************************************************************************
  //2-Jet Bin : BTag Efficiency for highest pt jet
  //*******************************************************************************
  printf("**********eff highest pt jet 2-j**********\n");
  double effttMC_btag_highestpt_2j[5],effttMC_btag_highestpt_2j_error[5],effttMC_btag_highestpt_tt_2j[5],effttMC_btag_highestpt_tt_2j_error[5];
  double effttDA_btag_highestpt_2j[5],effttDA_btag_highestpt_2j_error[5],effttMC_btag_highestpt_tw_2j[5],effttMC_btag_highestpt_tw_2j_error[5];

  for(int i=0; i<5; i++) {

    //MC efficiencies
    effttMC_btag_highestpt_tt_2j[i] = (btag_highestpt_2j_num[1][i])/(btag_highestpt_2j_den[1][i]);
    effttMC_btag_highestpt_tw_2j[i] = (btag_highestpt_2j_num[2][i])/(btag_highestpt_2j_den[2][i]);
    effttMC_btag_highestpt_2j[i]    = (btag_highestpt_2j_num[1][i] + btag_highestpt_2j_num[2][i]) / (btag_highestpt_2j_den[1][i] + btag_highestpt_2j_den[2][i]);
    
    effttMC_btag_highestpt_tt_2j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tt_2j[i])*effttMC_btag_highestpt_tt_2j[i]/(btag_highestpt_2j_den[1][i])*
                                                 (btag_highestpt_2j_den_error[1][i])/(btag_highestpt_2j_den[1][i]));    
    effttMC_btag_highestpt_tw_2j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tw_2j[i])*effttMC_btag_highestpt_tw_2j[i]/(btag_highestpt_2j_den[2][i])*
                                                 (btag_highestpt_2j_den_error[2][i])/(btag_highestpt_2j_den[2][i]));
    effttMC_btag_highestpt_2j_error[i]    = sqrt((1.0-effttMC_btag_highestpt_2j[i])*effttMC_btag_highestpt_2j[i]/(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i])*
                                                 (btag_highestpt_2j_den_error[1][i]+btag_highestpt_2j_den_error[2][i])/(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i]));

    //Data efficiencies
    effttDA_btag_highestpt_2j[i] = (btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i]-btag_highestpt_2j_num[2][i]*btagSF)/
      (btag_highestpt_2j_den[3][i]-btag_highestpt_2j_den[0][i]-btag_highestpt_2j_den[2][i]*btagSF);    
    effttDA_btag_highestpt_2j_error[i] = sqrt((1-effttDA_btag_highestpt_2j[i])*effttDA_btag_highestpt_2j[i]/btag_highestpt_2j_den[3][i]);
  }

  for(int i=0; i<5; i++) printf("scaleFactor2j(%d) --> %6.3f +/- %6.3f\n",i,
             (btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),
             sqrt(btag_highestpt_2j_num[3][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]));
  //double TopBkgScaleFactor_2Jet = (btag_highestpt_2j_num[3][4]-btag_highestpt_2j_num[0][4])/(btag_highestpt_2j_num[1][4]+btag_highestpt_2j_num[2][4]);
  //double TopBkgScaleFactorUncertainty_2Jet = sqrt(btag_highestpt_2j_num[3][4])/(btag_highestpt_2j_num[1][4]+btag_highestpt_2j_num[2][4]);

  for(int i=0; i<5; i++) {
    printf("numerator  (%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_highestpt_2j_num[3][i],btag_highestpt_2j_num[0][i],(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),btag_highestpt_2j_num[1][i],btag_highestpt_2j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_highestpt_2j_den[3][i],btag_highestpt_2j_den[0][i],(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i]),btag_highestpt_2j_den[1][i],btag_highestpt_2j_den[2][i]);
  }

  printf("channel         eff_tttw             eff_tt                eff_tw               eff_data                         ScaleFactor\n");
  for(int i=0; i<5; i++) {
    printf("eff (%s): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f      : scaleFactor2j(%s) --> %6.3f +/- %6.3f\n",classLabel[i],
           effttMC_btag_highestpt_2j[i]   ,effttMC_btag_highestpt_2j_error[i],effttMC_btag_highestpt_tt_2j[i],effttMC_btag_highestpt_tt_2j_error[i],
           effttMC_btag_highestpt_tw_2j[i],effttMC_btag_highestpt_tw_2j_error[i],effttDA_btag_highestpt_2j[i],effttDA_btag_highestpt_2j_error[i],
           classLabel[i],(btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),
           sqrt(btag_highestpt_2j_num[3][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]));
//     printf("scaleFactor2j(%s) --> %6.3f +/- %6.3f\n",classLabel[i],(btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),
//            sqrt(btag_highestpt_2j_num[3][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]));    
  }

  printf("****************************************************************************************************************************************\n");

  //*******************************************************************************
  //1-Jet Bin : BTag Efficiency for highest pt jet
  //*******************************************************************************
  printf("**********eff highest pt jet 1-j**********\n");
  double effttMC_btag_highestpt_1j[5],effttMC_btag_highestpt_1j_error[5],effttMC_btag_highestpt_tt_1j[5],effttMC_btag_highestpt_tt_1j_error[5];
  double effttDA_btag_highestpt_1j[5],effttDA_btag_highestpt_1j_error[5],effttMC_btag_highestpt_tw_1j[5],effttMC_btag_highestpt_tw_1j_error[5];


  for(int i=0; i<5; i++) {
    //MC btag efficiencies
    effttMC_btag_highestpt_tt_1j[i] = btag_highestpt_1j_num[1][i] / btag_highestpt_1j_den[1][i];
    effttMC_btag_highestpt_tw_1j[i] = btag_highestpt_1j_num[2][i] / btag_highestpt_1j_den[2][i];
    effttMC_btag_highestpt_1j[i]    = (btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i])/(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i]);
    
    effttMC_btag_highestpt_1j_error[i]    = sqrt((1.0-effttMC_btag_highestpt_1j[i])*effttMC_btag_highestpt_1j[i]/(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])*
                                                 (btag_highestpt_1j_den_error[1][i]+btag_highestpt_1j_den_error[2][i])/(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i]));
    effttMC_btag_highestpt_tt_1j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tt_1j[i])*effttMC_btag_highestpt_tt_1j[i]/(btag_highestpt_1j_den[1][i])*
                                                 (btag_highestpt_1j_den_error[1][i])/(btag_highestpt_1j_den[1][i]));
    effttMC_btag_highestpt_tw_1j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tw_1j[i])*effttMC_btag_highestpt_tw_1j[i]/(btag_highestpt_1j_den[2][i])*
                                                 (btag_highestpt_1j_den_error[2][i])/(btag_highestpt_1j_den[2][i]));
    
    //Data btag efficiencies
    effttDA_btag_highestpt_1j[i] = (btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i]-btag_highestpt_1j_num[2][i]*btagSF)/(btag_highestpt_1j_den[3][i]-btag_highestpt_1j_den[0][i]-btag_highestpt_1j_den[2][i]*btagSF);
    effttDA_btag_highestpt_1j_error[i] = sqrt((1-effttDA_btag_highestpt_1j[i])*effttDA_btag_highestpt_1j[i]/btag_highestpt_1j_den[3][i]);    
  }

  for(int i=0; i<5; i++) {
    printf("numerator  (%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_highestpt_1j_num[3][i],btag_highestpt_1j_num[0][i],(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]),btag_highestpt_1j_num[1][i],btag_highestpt_1j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_highestpt_1j_den[3][i],btag_highestpt_1j_den[0][i],(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i]),btag_highestpt_1j_den[1][i],btag_highestpt_1j_den[2][i]);
  }

  printf("channel       eff_tttw           eff_tt              eff_tw               eff_data              \n");
  for(int i=0; i<5; i++) {
    printf("eff (%d): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f  \n",i,
           effttMC_btag_highestpt_1j[i]   ,effttMC_btag_highestpt_1j_error[i],effttMC_btag_highestpt_tt_1j[i],effttMC_btag_highestpt_tt_1j_error[i],
           effttMC_btag_highestpt_tw_1j[i],effttMC_btag_highestpt_tw_1j_error[i],effttDA_btag_highestpt_1j[i],effttDA_btag_highestpt_1j_error[i]);
  }
  
  // we use the combined efficiency obtained in the 2-j bin, instead of the obtained final state by final state
  double estimationMC_btag_highestpt_1j[5];     
  double estimationMC_btag_highestpt_1j_err[5]; 
  double estimationDA_btag_highestpt_1j[5]; 
  double estimationDA_btag_highestpt_1j_error[5]; 

  for(int i=0; i<5; i++) {
    estimationMC_btag_highestpt_1j[i] = (1-effttMC_btag_highestpt_2j[4])/effttMC_btag_highestpt_2j[4]*(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]);
    estimationMC_btag_highestpt_1j_err[i] = effttMC_btag_highestpt_tt_1j_error[i]/effttMC_btag_highestpt_2j[4]/effttMC_btag_highestpt_2j[4]*(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]);
    estimationDA_btag_highestpt_1j[i] = (1-effttDA_btag_highestpt_2j[4])/effttDA_btag_highestpt_2j[4]*(btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i]);

    estimationDA_btag_highestpt_1j_error[i] = sqrt(((btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i])*effttDA_btag_highestpt_2j_error[4]/effttDA_btag_highestpt_2j[4]/effttDA_btag_highestpt_2j[4])*
                                                   ((btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i])*effttDA_btag_highestpt_2j_error[4]/effttDA_btag_highestpt_2j[4]/effttDA_btag_highestpt_2j[4])+
                                                   (1-effttDA_btag_highestpt_2j[4])/effttDA_btag_highestpt_2j[4]*
                                                   (1-effttDA_btag_highestpt_2j[4])/effttDA_btag_highestpt_2j[4]*btag_highestpt_1j_num[3][i]);
    
  }

  printf("Predicted ttbar+tW background for 1jet analysis (fails btag):  top background scale factor\n");
  printf("               MC(tt + tW)     Predicted from MC     | Prediction from Data   |    Scale Factor \n");
  for(int i=0; i<5; i++) {
    printf("top 1-jet(%s):    %7.3f        %7.3f +/- %6.3f       | %7.3f +/- %6.3f        |    %5.3f +/- %5.3f\n",classLabel[i],
           (btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]),
           estimationMC_btag_highestpt_1j[i],estimationMC_btag_highestpt_1j_err[i],
           estimationDA_btag_highestpt_1j[i],estimationDA_btag_highestpt_1j_error[i],
           estimationDA_btag_highestpt_1j[i]      /((btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i])),
           estimationDA_btag_highestpt_1j_error[i]/((btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i])));
  }

  printf("**********************************************************\n");

  //*******************************************************************************
  //0-Jet Bin
  //*******************************************************************************

  printf("**********eff low pt jet 1-j**********\n");

  //Single Top Scale Factor for 0-Jet bin events which failed b-tagging
  double btagSFTW = estimationDA_btag_highestpt_1j[4] /((btag_highestpt_1j_den[1][4]+btag_highestpt_1j_den[2][4])-(btag_highestpt_1j_num[1][4]+btag_highestpt_1j_num[2][4]));
  printf("btagSFTW = %f\n",btagSFTW);
  double TopBkgScaleFactor_1Jet = btagSFTW;
  double TopBkgScaleFactorUncertainty_1Jet = estimationDA_btag_highestpt_1j_error[4]/((btag_highestpt_1j_den[1][4]+btag_highestpt_1j_den[2][4])-(btag_highestpt_1j_num[1][4]+btag_highestpt_1j_num[2][4]));

  double effttMC_btag_lowpt_1j[5],effttMC_btag_lowpt_1j_error[5],effttMC_btag_lowpt_tt_1j[5],effttMC_btag_lowpt_tt_1j_error[5];
  double effttDA_btag_lowpt_1j[5],effttDA_btag_lowpt_1j_error[5],effttMC_btag_lowpt_tw_1j[5],effttMC_btag_lowpt_tw_1j_error[5];

  for(int i=0; i<5; i++) {
    //MC btag efficiencies
    effttMC_btag_lowpt_1j[i]    = (btag_lowpt_1j_num[1][i]+btag_lowpt_1j_num[2][i])/(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]);
    effttMC_btag_lowpt_tt_1j[i] = (btag_lowpt_1j_num[1][i])/(btag_lowpt_1j_den[1][i]                          );
    effttMC_btag_lowpt_tw_1j[i] = (btag_lowpt_1j_num[2][i])/(btag_lowpt_1j_den[2][i]);
    effttMC_btag_lowpt_1j_error[i]    = sqrt((1.0-effttMC_btag_lowpt_1j[i])*effttMC_btag_lowpt_1j[i]/(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i])*
                                             (btag_lowpt_1j_den_error[1][i]+btag_lowpt_1j_den_error[2][i])/(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]));
    effttMC_btag_lowpt_tt_1j_error[i] = sqrt((1.0-effttMC_btag_lowpt_tt_1j[i])*effttMC_btag_lowpt_tt_1j[i]/(btag_lowpt_1j_den[1][i])*
                                             (btag_lowpt_1j_den_error[1][i])/(btag_lowpt_1j_den[1][i]));
    effttMC_btag_lowpt_tw_1j_error[i] = sqrt((1.0-effttMC_btag_lowpt_tw_1j[i])*effttMC_btag_lowpt_tw_1j[i]/(btag_lowpt_1j_den[2][i])*
                                             (btag_lowpt_1j_den_error[2][i])/(btag_lowpt_1j_den[2][i]));

    //Data btag efficiencies
    effttDA_btag_lowpt_1j[i]       = (btag_lowpt_1j_num[3][i]-btag_lowpt_1j_num[0][i]-btag_lowpt_1j_num[2][i]*btagSF*btagSFTW)/ (btag_lowpt_1j_den[3][i]-btag_lowpt_1j_den[0][i]-btag_lowpt_1j_den[2][i]*btagSF*btagSFTW);
    effttDA_btag_lowpt_1j_error[i] = sqrt((1-effttDA_btag_lowpt_1j[i])*effttDA_btag_lowpt_1j[i]/btag_lowpt_1j_den[3][i]);
  }

  for(int i=0; i<5; i++) {
    printf("numerator(%s)   --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_lowpt_1j_num[3][i],btag_lowpt_1j_num[0][i],(btag_lowpt_1j_num[1][i]+btag_lowpt_1j_num[2][i]),btag_lowpt_1j_num[1][i],btag_lowpt_1j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_lowpt_1j_den[3][i],btag_lowpt_1j_den[0][i],(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]),btag_lowpt_1j_den[1][i],btag_lowpt_1j_den[2][i]);
  }

  printf("channel       eff_tttw         eff_tt            eff_tw           eff_data              \n");
  for(int i=0; i<5; i++) {
    printf("eff (%d): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f\n",i,
           effttMC_btag_lowpt_1j[i]   ,effttMC_btag_lowpt_1j_error[i],effttMC_btag_lowpt_tt_1j[i],effttMC_btag_lowpt_tt_1j_error[i],
           effttMC_btag_lowpt_tw_1j[i],effttMC_btag_lowpt_tw_1j_error[i],effttDA_btag_lowpt_1j[i],effttDA_btag_lowpt_1j_error[i]);
  }

  printf("**********************************************************\n");

  //*******************************************************************************
  //Closure Test 0-Jet Bin
  //*******************************************************************************

  printf("**********eff low pt jet 0-j**********\n");
  double ftw_b[5]; 
  for(int i=0; i<5; i++) {
    ftw_b[i] = effttMC_btag_lowpt_tw_1j[i];
    printf("ftw_b(%d) = %5.3f \n",i,ftw_b[i]);
  }

  double N_top_expected_0j[5]; 
  double fttbar[5]; 

  for(int i=0; i<5; i++) {
    N_top_expected_0j[i] = (btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i])-(btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i]);
    fttbar[i] = (btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]*ftw_b[i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]);
// fttbar[i] = btag_lowpt_0j_den[1][i]/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]);
  }

  double effMC_btag_lowpt_tt_0j_expected[5];  
  double effMC_btag_lowpt_tw_0j_expected[5];  
  double effMC_btag_lowpt_tt_0j[5];           
  double effMC_btag_lowpt_tw_0j[5];           

  for(int i=0; i<5; i++) {
    effMC_btag_lowpt_tt_0j_expected[i] = btag_lowpt_0j_num[1][i]/btag_lowpt_0j_den[1][i];
    effMC_btag_lowpt_tw_0j_expected[i] = btag_lowpt_0j_num[2][i]/btag_lowpt_0j_den[2][i];
    effMC_btag_lowpt_tt_0j[i]          = 1-(1-effttMC_btag_lowpt_tt_1j[i])*(1-effttMC_btag_lowpt_tt_1j[i]);
    effMC_btag_lowpt_tw_0j[i]          = effttMC_btag_lowpt_tt_1j[i];    
  }

  printf("channel        ttbar MC ( 0Jet / Extrapolated from 1Jet )            tW MC ( 0Jet / Extrapolated from 1Jet) \n");
  for(int i=0; i<5; i++) { 
    printf("(%s),                  %5.3f/%5.3f                                              %5.3f/%5.3f\n",
           classLabel[i],effMC_btag_lowpt_tt_0j_expected[i],effMC_btag_lowpt_tt_0j[i],effMC_btag_lowpt_tw_0j_expected[i],effMC_btag_lowpt_tw_0j[i]);
  }

  // begin get closure test closing!!!!!!!!!!!!!!!
  // double ftw_2b[5];
  // for(int i=0; i<5; i++) {
  //   ftw_2b[i] = (effMC_btag_lowpt_tw_0j_expected[i]-effttMC_btag_lowpt_tt_1j[i])/(effMC_btag_lowpt_tt_0j[i]-effttMC_btag_lowpt_tt_1j[i]);
  //   printf("ftw_2b(%d) = %5.3f - ",i,ftw_2b[i]);
  //   fttbar[i] = (btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]*ftw_2b[i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]);
  //   printf("fttbar(%d) = %5.3f | ",i,fttbar[i]);
  // }
  // printf("\n");
  // end get closure test closing!!!!!!!!!!!!!!!


  double sigma_ftop[2]={0.00,0.17};
  double effMC_btag_lowpt_0j[5]; 
  double effDA_btag_lowpt_0j[5]; 
  double effMC_btag_lowpt_0j_error[5];
  double effDA_btag_lowpt_0j_error[5];

  for(int i=0; i<5; i++) {
    effMC_btag_lowpt_0j[i] = fttbar[i]*effMC_btag_lowpt_tt_0j[i]+(1-fttbar[i])*effMC_btag_lowpt_tw_0j[i];
    effDA_btag_lowpt_0j[i] = fttbar[i]*(1-(1-effttDA_btag_lowpt_1j[i])*(1-effttDA_btag_lowpt_1j[i]))+(1-fttbar[i])*effttDA_btag_lowpt_1j[i];

    effMC_btag_lowpt_0j_error[i] = sqrt(sigma_ftop[0]*sigma_ftop[0]*(effMC_btag_lowpt_0j[i]*(1-effMC_btag_lowpt_0j[i]))*(effMC_btag_lowpt_0j[i]*(1-effMC_btag_lowpt_0j[i]))+
    					effttMC_btag_lowpt_tt_1j_error[i]*effttMC_btag_lowpt_tt_1j_error[i]*(fttbar[i]*(1-2*effMC_btag_lowpt_0j[i])+1)*(fttbar[i]*(1-2*effMC_btag_lowpt_0j[i])+1));
    effDA_btag_lowpt_0j_error[i] = sqrt(sigma_ftop[1]*sigma_ftop[1]*(effDA_btag_lowpt_0j[i]*(1-effDA_btag_lowpt_0j[i]))*(effDA_btag_lowpt_0j[i]*(1-effDA_btag_lowpt_0j[i]))+
    					effttDA_btag_lowpt_1j_error[i]*effttDA_btag_lowpt_1j_error[i]*(fttbar[i]*(1-2*effDA_btag_lowpt_0j[i])+1)*(fttbar[i]*(1-2*effDA_btag_lowpt_0j[i])+1));
    
  }

  printf("top tagging efficiency\n");
  printf("Channel    fttbar        Eff toptag(MC)    Eff toptag(MC extrapolated)       Eff toptab Data \n");
  for(int i=0; i<5; i++) {
    printf("(%s)       %5.3f,        : %6.3f                 %6.3f +/- %6.3f             %6.3f +/- %6.3f\n",
           classLabel[i],fttbar[i],
           (btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]),
           effMC_btag_lowpt_0j[i],effMC_btag_lowpt_0j_error[i],
           effDA_btag_lowpt_0j[i],effDA_btag_lowpt_0j_error[i]);
  }

  double sigma_0f_bck = 0.20;
  double estimationMC_btag_lowpt_0j[5]; 
  double estimationDA_btag_lowpt_0j[5]; 
  double estimationMC_btag_lowpt_0j_error[5]; 
  double estimationDA_btag_lowpt_0j_error[5]; 
  // we use the combined efficiency obtained in the 1-j bin, instead of the obtained final state by final state
  for(int i=0; i<5; i++) {
    estimationMC_btag_lowpt_0j[i] = (1-effMC_btag_lowpt_0j[4])/effMC_btag_lowpt_0j[4]*(btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i]);
    estimationDA_btag_lowpt_0j[i] = (1-effDA_btag_lowpt_0j[4])/effDA_btag_lowpt_0j[4]*(btag_lowpt_0j_num[3][i]-btag_lowpt_0j_num[0][i]);
    estimationMC_btag_lowpt_0j_error[i] = sqrt(((btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i])*effMC_btag_lowpt_0j_error[4]/effMC_btag_lowpt_0j[4]/effMC_btag_lowpt_0j[4])*
                                               ((btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i])*effMC_btag_lowpt_0j_error[4]/effMC_btag_lowpt_0j[4]/effMC_btag_lowpt_0j[4]));
    
    estimationDA_btag_lowpt_0j_error[i] = sqrt(((btag_lowpt_0j_num[3][i]-btag_lowpt_0j_num[0][i])*effDA_btag_lowpt_0j_error[4]/effDA_btag_lowpt_0j[4]/effDA_btag_lowpt_0j[4])*
                                               ((btag_lowpt_0j_num[3][i]-btag_lowpt_0j_num[0][i])*effDA_btag_lowpt_0j_error[4]/effDA_btag_lowpt_0j[4]/effDA_btag_lowpt_0j[4])+
                                               (1-estimationDA_btag_lowpt_0j[i])/estimationDA_btag_lowpt_0j[i]*
                                               (1-estimationDA_btag_lowpt_0j[i])/estimationDA_btag_lowpt_0j[i]*btag_lowpt_0j_num[3][i]+
                                               TMath::Power(sigma_0f_bck*btag_lowpt_0j_num[0][i]*(1-effMC_btag_lowpt_0j[4])/effMC_btag_lowpt_0j[4],2));
  }

  printf("0-Jet Top Background\n");
  printf("Channel    top-tagged region    bkgs top-tagged          MC tt+tW              MC top            Data top   \n");
  printf("              Event Count       region (non-top)    non-top-tagged count     estimation         estimation  \n");
  for(int i=0; i<5; i++) {
    printf("(%s)              %3d               %6.3f                %6.3f             %6.3f +/- %6.3f    %6.3f +/- %6.3f\n",
           classLabel[i],
           (int)btag_lowpt_0j_num[3][i],btag_lowpt_0j_num[0][i],N_top_expected_0j[i],
           estimationMC_btag_lowpt_0j[i],estimationMC_btag_lowpt_0j_error[i],
           estimationDA_btag_lowpt_0j[i],estimationDA_btag_lowpt_0j_error[i]);
  }

  // Writing results in latex format
  double NonTaggedTopMC = 0; for(int i=0; i<5; i++) NonTaggedTopMC = NonTaggedTopMC + btag_vbf_2j_den1[1][i]+btag_vbf_2j_den2[1][i];
  double NonTaggedTopDA = 0; for(int i=0; i<5; i++) NonTaggedTopDA = NonTaggedTopDA + evtDA_vbf_2j[i];
  double NonTaggedTopDA_error = 0; for(int i=0; i<5; i++) NonTaggedTopDA_error = NonTaggedTopDA_error + evtDA_vbf_2j_error[i]*evtDA_vbf_2j_error[i];
  int Chan = 4;
  printf("*****************************************\n");
  printf("\\begin{table}\n");
  printf("\\begin{center}\n");
  printf("{\\tiny\n");
  printf("\\hline\n");
  printf(" Sample & 0-jet & 1-jet & 2-jet \\\\\n");
  printf("\\hline\n");
  printf("estimated top events in simulation  & %5.1f $\\pm$ %5.1f &  %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f \\\\\n",N_top_expected_0j[Chan],estimationMC_btag_lowpt_0j_error[Chan],
  (btag_highestpt_1j_den[1][Chan]+btag_highestpt_1j_den[2][Chan])-(btag_highestpt_1j_num[1][Chan]+btag_highestpt_1j_num[2][Chan]),estimationMC_btag_highestpt_1j_err[Chan],
  NonTaggedTopMC,NonTaggedTopMC*0.04);
  printf("tagging efficiency     (\\%)         & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & - \\\\ \n",effDA_btag_lowpt_0j[Chan]*100,effDA_btag_lowpt_0j_error[Chan]*100,effttDA_btag_highestpt_1j[Chan]*100,effttDA_btag_highestpt_1j_error[Chan]*100);
  printf("data events in control region       & %4d & %4d & -\n",(int)btag_lowpt_0j_num[3][Chan],(int)btag_highestpt_1j_num[3][Chan]);
  printf("background events in control region & %5.1f $\\pm$ %5.1f &  %5.1f $\\pm$ %5.1f & - \\\\ \n",btag_lowpt_0j_num[0][Chan],btag_lowpt_0j_num[0][Chan]*0.15,
  btag_highestpt_1j_num[0][Chan],btag_highestpt_1j_num[0][Chan]*0.20);
  printf("top estimation in data              &  %5.1f $\\pm$ %5.1f &  %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f \\\\\n",estimationDA_btag_lowpt_0j[Chan],estimationDA_btag_lowpt_0j_error[Chan],
  estimationDA_btag_highestpt_1j[Chan],estimationDA_btag_highestpt_1j_error[Chan],NonTaggedTopDA,sqrt(NonTaggedTopDA_error));
  printf("data/simulation scale factor        &  %5.2f $\\pm$ %5.2f &  %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\\n",
  estimationDA_btag_lowpt_0j[Chan] / N_top_expected_0j[Chan],  estimationDA_btag_lowpt_0j_error[Chan] / N_top_expected_0j[Chan],
  TopBkgScaleFactor_1Jet,TopBkgScaleFactorUncertainty_1Jet,TopBkgScaleFactor_2Jet_vbf[1],TopBkgScaleFactorUncertainty_2Jet_vbf[1]);
  printf("\\hline\n");
  printf("\\end{tabular}\n");
  printf("}\n");
  printf("\\end{center}\n");
  printf("\\end{table}\n");
  printf("*****************************************\n");

  // Writing scale factors in *.h files
  double TopBkgScaleFactor_0Jet = estimationDA_btag_lowpt_0j[4] / N_top_expected_0j[4];
  double TopBkgScaleFactorUncertainty_0Jet = estimationDA_btag_lowpt_0j_error[4] / N_top_expected_0j[4];

  char TopVBFBkgScaleFactorsName[200];
  sprintf(TopVBFBkgScaleFactorsName,"TopVBFBkgScaleFactors.h");
  ofstream outVBF(TopVBFBkgScaleFactorsName);

  outVBF << "static Double_t TopVBFBkgScaleFactor(Int_t mH) {" << endl;
  outVBF << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outVBF << mH[i+1];
    if (i < nmass-1-1) outVBF << ",";    
  }
  outVBF << "};" << endl;
  outVBF << "  Double_t TopVBFBkgSF[" << nmass-1 << "] = { " << endl;
  for (UInt_t i = 0; i < nmass-1; ++i) {
    outVBF << TopBkgScaleFactor_2Jet_vbf[i+1];
    if (i < nmass-1-1) outVBF << ",";
  }
  outVBF << "};" << endl;
  outVBF << "  Int_t massIndex = -1;" << endl;
  outVBF << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outVBF << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outVBF << "  }" << endl;
  outVBF << "  if (massIndex >= 0) {" << endl;  
  outVBF << "    return TopVBFBkgSF[massIndex];" << endl;
  outVBF << "  } else { assert(0); }" << endl;
  outVBF << "  return 1.0;" << endl;
  outVBF << "}" << endl;

  outVBF << "static Double_t TopVBFBkgScaleFactorKappa(Int_t mH) {" << endl;
  outVBF << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outVBF << mH[i+1];
    if (i < nmass-1-1) outVBF << ",";    
  }
  outVBF << "};" << endl;
  outVBF << "  Double_t TopVBFBkgKappa[" << nmass-1 << "] = { " << endl;
  for (UInt_t i = 0; i < nmass-1; ++i) {
    if(TopBkgScaleFactor_2Jet_vbf[i+1] == 0) TopBkgScaleFactor_2Jet_vbf[i+1] = 1.0;
    outVBF << (1.0 + TopBkgScaleFactorUncertainty_2Jet_vbf[i+1]/TopBkgScaleFactor_2Jet_vbf[i+1]);    
    if (i < nmass-1-1) outVBF << ",";
  }
  outVBF << "};" << endl;
  outVBF << "  Int_t massIndex = -1;" << endl;
  outVBF << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outVBF << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outVBF << "  }" << endl;
  outVBF << "  if (massIndex >= 0) {" << endl;  
  outVBF << "    return TopVBFBkgKappa[massIndex];" << endl;
  outVBF << "  } else { assert(0); }" << endl;
  outVBF << "  return 1.0;" << endl;
  outVBF << "}" << endl;

  outVBF << endl;


  outVBF << "static Double_t TopVBFOFYield(Int_t mH) {" << endl;
  outVBF << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outVBF << mH[i+1];
    if (i < nmass-1-1) outVBF << ",";    
  }
  outVBF << "};" << endl;
  outVBF << "  Double_t TopVBFBkgOF[" << nmass-1 << "] = { " << endl;
  for (UInt_t i = 0; i < nmass-1; ++i) {
    outVBF << TopBkgScaleFactor_2Jet_vbfOF[i+1];
    if (i < nmass-1-1) outVBF << ",";
  }
  outVBF << "};" << endl;
  outVBF << "  Int_t massIndex = -1;" << endl;
  outVBF << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outVBF << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outVBF << "  }" << endl;
  outVBF << "  if (massIndex >= 0) {" << endl;  
  outVBF << "    return TopVBFBkgOF[massIndex];" << endl;
  outVBF << "  } else { assert(0); }" << endl;
  outVBF << "  return 1.0;" << endl;
  outVBF << "}" << endl;

  outVBF << "static Double_t TopVBFOFYieldKappa(Int_t mH) {" << endl;
  outVBF << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outVBF << mH[i+1];
    if (i < nmass-1-1) outVBF << ",";    
  }
  outVBF << "};" << endl;
  outVBF << "  Double_t TopVBFOFBkgKappa[" << nmass-1 << "] = { " << endl;
  for (UInt_t i = 0; i < nmass-1; ++i) {
    if(TopBkgScaleFactor_2Jet_vbfOF[i+1] == 0) TopBkgScaleFactor_2Jet_vbfOF[i+1] = 1.0;
    outVBF << (1.0 + TopBkgScaleFactorUncertainty_2Jet_vbfOF[i+1]/TopBkgScaleFactor_2Jet_vbfOF[i+1]);    
    if (i < nmass-1-1) outVBF << ",";
  }
  outVBF << "};" << endl;
  outVBF << "  Int_t massIndex = -1;" << endl;
  outVBF << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outVBF << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outVBF << "  }" << endl;
  outVBF << "  if (massIndex >= 0) {" << endl;  
  outVBF << "    return TopVBFOFBkgKappa[massIndex];" << endl;
  outVBF << "  } else { assert(0); }" << endl;
  outVBF << "  return 1.0;" << endl;
  outVBF << "}" << endl;

  outVBF << endl;


  outVBF << "static Double_t TopVBFSFYield(Int_t mH) {" << endl;
  outVBF << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outVBF << mH[i+1];
    if (i < nmass-1-1) outVBF << ",";    
  }
  outVBF << "};" << endl;
  outVBF << "  Double_t TopVBFBkgSF[" << nmass-1 << "] = { " << endl;
  for (UInt_t i = 0; i < nmass-1; ++i) {
    outVBF << TopBkgScaleFactor_2Jet_vbfSF[i+1];
    if (i < nmass-1-1) outVBF << ",";
  }
  outVBF << "};" << endl;
  outVBF << "  Int_t massIndex = -1;" << endl;
  outVBF << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outVBF << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outVBF << "  }" << endl;
  outVBF << "  if (massIndex >= 0) {" << endl;  
  outVBF << "    return TopVBFBkgSF[massIndex];" << endl;
  outVBF << "  } else { assert(0); }" << endl;
  outVBF << "  return 1.0;" << endl;
  outVBF << "}" << endl;

  outVBF << "static Double_t TopVBFSFYieldKappa(Int_t mH) {" << endl;
  outVBF << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outVBF << mH[i+1];
    if (i < nmass-1-1) outVBF << ",";    
  }
  outVBF << "};" << endl;
  outVBF << "  Double_t TopVBFSFBkgKappa[" << nmass-1 << "] = { " << endl;
  for (UInt_t i = 0; i < nmass-1; ++i) {
    if(TopBkgScaleFactor_2Jet_vbfSF[i+1] == 0) TopBkgScaleFactor_2Jet_vbfSF[i+1] = 1.0;
    outVBF << (1.0 + TopBkgScaleFactorUncertainty_2Jet_vbfSF[i+1]/TopBkgScaleFactor_2Jet_vbfSF[i+1]);    
    if (i < nmass-1-1) outVBF << ",";
  }
  outVBF << "};" << endl;
  outVBF << "  Int_t massIndex = -1;" << endl;
  outVBF << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outVBF << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outVBF << "  }" << endl;
  outVBF << "  if (massIndex >= 0) {" << endl;  
  outVBF << "    return TopVBFSFBkgKappa[massIndex];" << endl;
  outVBF << "  } else { assert(0); }" << endl;
  outVBF << "  return 1.0;" << endl;
  outVBF << "}" << endl;

  outVBF << endl;


  ofstream outf("TopBkgScaleFactors.h");

  outf << "static Double_t TopBkgScaleFactor(Int_t jetBin) {" << endl;
  outf << "  assert(jetBin >=0 && jetBin <= 2);" << endl;
  
  outf << "  Double_t TopBkgScaleFactor[3] = { " 
       << TopBkgScaleFactor_0Jet << ", "
       << TopBkgScaleFactor_1Jet << ", "
       << TopBkgScaleFactor_2Jet_vbf[0] << "  "
       << " };" << endl;
  outf << "  return TopBkgScaleFactor[jetBin];" << endl;
  outf << "}" << endl;
  outf << endl;


  outf << "static Double_t TopBkgScaleFactorKappa(Int_t jetBin) {" << endl;
  
  outf << "  assert(jetBin >=0 && jetBin <= 2);" << endl;
  outf << "  Double_t TopBkgScaleFactorKappa[3] = { " 
       << (1.0 + TopBkgScaleFactorUncertainty_0Jet/TopBkgScaleFactor_0Jet) << ", "
       << (1.0 + TopBkgScaleFactorUncertainty_1Jet/TopBkgScaleFactor_1Jet) << ", "
       << (1.0 + TopBkgScaleFactorUncertainty_2Jet_vbf[0]/TopBkgScaleFactor_2Jet_vbf[0]) << "  "
       << " };" << endl;

  outf << "  return TopBkgScaleFactorKappa[jetBin];" << endl;
  outf << "}" << endl;
  outf << endl;

  outf.close();

  return;
}
