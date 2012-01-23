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
#include "THStack.h"
#include "TPaveText.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h"
#include "Smurf/Analysis/HWWlvlv/TopBkgScaleFactors.h"
#include "Smurf/Analysis/HWWlvlv/WWBkgScaleFactors.h"
#include "Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors.h"
#include "Smurf/Analysis/HWWlvlv/HWWCuts.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double DeltaPhi(double phi1, double phi2);

int    verboseLevel =   0;
const double sigmaB = 0.35;




//*************************************************************************************************
//Draws the signal and background histograms together and makes gif file
//*************************************************************************************************
void DrawDataBkgHistogram(TH1F* data, THStack* bkg, TLegend *legend , string histname, 
                             Double_t minY = -999 , Double_t maxY = -999,
                             Bool_t useLogY = kFALSE) {

  string filename = histname;
  if (useLogY) filename += "_logY";
//  filename += ".eps";

  TCanvas *cv = new TCanvas(histname.c_str(), histname.c_str(), 0,0,800,600);
  if (useLogY) cv->SetLogy();

  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);

  THStack *tmpBkg = (THStack*)bkg->Clone("bkgStack");

  //set best Y-range
  double MaxY = data->GetMaximum();

  tmpBkg->Draw("hist");
  if (tmpBkg->GetMaximum() > MaxY)
    MaxY = tmpBkg->GetMaximum();

  tmpBkg->GetYaxis()->SetTitleOffset(1.1);
  data->GetYaxis()->SetTitleOffset(1.1);
  tmpBkg->GetXaxis()->SetTitleOffset(1.0);
  data->GetXaxis()->SetTitleOffset(1.0);

  tmpBkg->SetMinimum(0.01);
  tmpBkg->SetMaximum(MaxY*1.2);  
  data->SetMinimum(0.01);
  data->SetMaximum(MaxY*1.2);
  
  if (minY != -999) {
    tmpBkg->SetMinimum(minY);
    data->SetMinimum(minY);    
  }
  if (maxY != -999) {
    tmpBkg->SetMaximum(maxY);
    data->SetMaximum(maxY);
  }

  //CMS Preliminary label
  TPaveText *prelimLabel = new TPaveText(0.21,0.85,0.41,0.90,"NDC");
  prelimLabel->SetTextColor(kBlack);
  prelimLabel->SetFillColor(kWhite);
  prelimLabel->SetBorderSize(0);
  prelimLabel->SetTextAlign(12);
  prelimLabel->SetTextSize(0.03);
  prelimLabel->AddText("CMS Preliminary 2011 #sqrt{s} = 7 TeV");
  prelimLabel->Draw();

  //Luminosity label
  TPaveText *tb = new TPaveText(0.21,0.77,0.41,0.82,"NDC");
  tb->SetTextColor(kBlack);
  tb->SetFillColor(kWhite);
  tb->SetBorderSize(0);
  tb->SetTextAlign(12);
  tb->AddText((string("#int#font[12]{L}dt = 4.6 pb^{ -1}")).c_str());
  tb->Draw();


  tmpBkg->Draw("samehist");
  data->Draw("sameE1");

  legend->Draw();
  cv->SaveAs((filename + ".gif").c_str());
  cv->SaveAs((filename + ".eps").c_str());


}








//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void MakeSameSignControlSampleTables
(
//  UInt_t  nJetsType   	 = 0,
 UInt_t  mH      	 = 0,
 TString dataInputFile = "/build/sixie/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/data_2l_skim2.root",
 TString bgdInputFile    = "/build/sixie/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/backgroundC_skim2.root",
 int period              = 2
 )
{

  bool wwPresel = false;
  if(mH == 0) {wwPresel = true; }

  TString dataFile1 = dataInputFile;
  TString bgdFile1 = bgdInputFile;

  unsigned int patternTopTag = SmurfTree::TopTag;

  float dilmass_cut = DileptonMassPreselectionCut(mH);
  if(wwPresel == true) dilmass_cut = 99999.;


  cout << "Using dilepton mass < " << dilmass_cut << endl;

  TChain *chdata = new TChain("tree");
  chdata->Add(dataFile1);

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TTree *data = (TTree*) chdata;
  TTree *background = (TTree*) chbackground;

  TString effPath  = "/data/smurf/data/LP2011/auxiliar/efficiency_results_v6_42x.root";
  TString fakePath = "/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.LP2011.root";
  TString puPath   = "/data/smurf/data/LP2011/auxiliar/puWeights_PU4_68mb.root";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  double scaleFactorLum = 2.121;
  if	 (period == 0){ // Run2011A
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Run2011A.root";
    //scaleFactorLum     = 2.1;minRun =      0;maxRun = 173692;
    scaleFactorLum     = 1.1;minRun =      0;maxRun = 167913;
  }
  else if(period == 1){ // Run2011B
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    //puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Run2011B.root";
    //scaleFactorLum     = 1.9;minRun = 173693;maxRun = 999999;
    puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Full2011.root";
    scaleFactorLum     = 3.6;minRun = 167914;maxRun = 999999;
  }
  else if(period == 2){ // Full2011
    effPath  = "/build/sixie/Thesis/auxiliar/Winter11_4700ipb/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/build/sixie/Thesis/auxiliar/Winter11_4700ipb/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/build/sixie/Thesis/auxiliar/Winter11_4700ipb/PileupReweighting.Summer11DYmm_To_Full2011.root";
    scaleFactorLum     = 4.7;minRun =      0;maxRun = 999999;
  }
  else if(period == 12){ // Full2011 with MVAIDIsoCombinedSameSigWP Lepton Selection
    effPath  = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/efficiency_results_MVAIDIsoCombinedSameSigWP_Full2011.root";
    fakePath = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/FakeRates_MVAIDIsoCombinedSameSigWP.root";
    puPath   = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Run2011B.root";
    scaleFactorLum     = 4.6;minRun =      0;maxRun = 999999;
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

  TFile *fPUS4File = TFile::Open(puPath.Data());
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;


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

  cout << "here01\n";

  //****************************************************************************
  // Yields and Histograms
  //****************************************************************************
  enum { kFakeElectron, kFakeMuon, kFakeLepton };
  enum { kMuMu, kEleEle, kEleMu, kMuEle, kSameFlavor, kDifferentFlavor, kAllFinalStates };
  enum { kZeroJetBin, kOneJetBin, kVBFBin };
  enum { kData, kFakeLeptonBkg, kWgamma, kWZ, kOther };
  //pt eta bins: { [all],[pt<20,0-1],[pt<20,1-1.479][pt<20,1.479-2.5],[pt>20,0-1],[pt>20,1-1.479][pt>20,1.479-2.5]}

  //[Fake Electron/Muon][Final State][Jet Bin][BkgType]

  vector<vector<vector<vector<double> > > > SameSignYields;
  vector<vector<vector<vector<double> > > > SameSignYieldsErrSqr;
  vector<vector<vector<vector<TH1F*> > > > SameSignPtMax;
  vector<vector<vector<vector<TH1F*> > > > SameSignPtMin;
  vector<vector<vector<vector<TH1F*> > > > SameSignDileptonMass;
  vector<vector<vector<vector<TH1F*> > > > SameSignDileptonPt;

  for(UInt_t i = 0; i < 3; ++i) {
    vector<vector<vector<double> > > SameSignYields_tmp1;
    vector<vector<vector<double> > > SameSignYieldsErrSqr_tmp1;
    vector<vector<vector<TH1F*> > > SameSignPtMax_tmp1;
    vector<vector<vector<TH1F*> > > SameSignPtMin_tmp1;
    vector<vector<vector<TH1F*> > > SameSignDileptonMass_tmp1;
    vector<vector<vector<TH1F*> > > SameSignDileptonPt_tmp1;
 
    for(UInt_t j = 0; j < 7; ++j) {
    vector<vector<double> > SameSignYields_tmp2;
    vector<vector<double> > SameSignYieldsErrSqr_tmp2;
    vector<vector<TH1F*> > SameSignPtMax_tmp2;
    vector<vector<TH1F*> > SameSignPtMin_tmp2;
    vector<vector<TH1F*> > SameSignDileptonMass_tmp2;
    vector<vector<TH1F*> > SameSignDileptonPt_tmp2;
 
      for(UInt_t k = 0; k < 3; ++k) {
        vector<double> SameSignYields_tmp3;
        vector<double> SameSignYieldsErrSqr_tmp3;
        vector<TH1F*>  SameSignPtMax_tmp3;
        vector<TH1F*>  SameSignPtMin_tmp3;
        vector<TH1F*>  SameSignDileptonMass_tmp3;
        vector<TH1F*>  SameSignDileptonPt_tmp3;

        for(UInt_t l = 0; l < 5; ++l) {
          char tempbuffer[200];
          string FakeBinLabel;
          if (i==0) FakeBinLabel = "FakeElectron";
          if (i==1) FakeBinLabel = "FakeMuon";
          if (i==2) FakeBinLabel = "FakeLepton";
          string FinalStateBinLabel;
          if (j==0) FinalStateBinLabel = "MuMu";
          if (j==0) FinalStateBinLabel = "EleEle";
          if (j==0) FinalStateBinLabel = "EleMu";
          if (j==0) FinalStateBinLabel = "MuEle";
          if (j==0) FinalStateBinLabel = "SameFlavor";
          if (j==0) FinalStateBinLabel = "DifferentFlavor";          
          if (j==0) FinalStateBinLabel = "AllFinalStates";
          string JetBinLabel;
          if (k == 0) JetBinLabel = "0JetBin";
          if (k == 1) JetBinLabel = "1JetBin";
          if (k == 2) JetBinLabel = "VBFBin";
          string TypeLabel;
          if (l==0) TypeLabel = "Data";
          if (l==1) TypeLabel = "FakeLepton";
          if (l==2) TypeLabel = "Wgamma";
          if (l==3) TypeLabel = "WZ";
          if (l==4) TypeLabel = "Other";
          sprintf(tempbuffer,"%s_%s_%s_%s", FakeBinLabel.c_str(), FinalStateBinLabel.c_str(), JetBinLabel.c_str(), TypeLabel.c_str());

          SameSignYields_tmp3.push_back(0.0);
          SameSignYieldsErrSqr_tmp3.push_back(0.0);

          SameSignPtMax_tmp3.push_back(new TH1F(("SameSignPtMax" + string(tempbuffer)).c_str(), ";PtMax [GeV/c]; Number of Events", 80, 0, 400));
          SameSignPtMin_tmp3.push_back(new TH1F(("SameSignPtMin" + string(tempbuffer)).c_str(), ";PtMin [GeV/c]; Number of Events", 80, 0, 400));
          SameSignDileptonMass_tmp3.push_back(new TH1F(("SameSignDileptonMass" + string(tempbuffer)).c_str(), ";DileptonMass [GeV/c^{2}]; Number of Events", 80, 0, 400));
          SameSignDileptonPt_tmp3.push_back(new TH1F(("SameSignPtMax" + string(tempbuffer)).c_str(), ";DileptonPt [GeV/c]; Number of Events", 80, 0, 400));


        }
        
        SameSignYields_tmp2.push_back(SameSignYields_tmp3);
        SameSignYieldsErrSqr_tmp2.push_back(SameSignYieldsErrSqr_tmp3);
        SameSignPtMax_tmp2.push_back(SameSignPtMax_tmp3);
        SameSignPtMin_tmp2.push_back(SameSignPtMin_tmp3);
        SameSignDileptonMass_tmp2.push_back(SameSignDileptonMass_tmp3);
        SameSignDileptonPt_tmp2.push_back(SameSignDileptonPt_tmp3);
      }
      SameSignYields_tmp1.push_back(SameSignYields_tmp2);
      SameSignYieldsErrSqr_tmp1.push_back(SameSignYieldsErrSqr_tmp2);
      SameSignPtMax_tmp1.push_back(SameSignPtMax_tmp2);
      SameSignPtMin_tmp1.push_back(SameSignPtMin_tmp2);
      SameSignDileptonMass_tmp1.push_back(SameSignDileptonMass_tmp2);
      SameSignDileptonPt_tmp1.push_back(SameSignDileptonPt_tmp2);
      
    }
    SameSignYields.push_back(SameSignYields_tmp1);
    SameSignYieldsErrSqr.push_back(SameSignYieldsErrSqr_tmp1);
    SameSignPtMax.push_back(SameSignPtMax_tmp1);
    SameSignPtMin.push_back(SameSignPtMin_tmp1);
    SameSignDileptonMass.push_back(SameSignDileptonMass_tmp1);
    SameSignDileptonPt.push_back(SameSignDileptonPt_tmp1);
  }


  cout << "Start\n";


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



//   Int_t nJetsType = 0;
//   {
//    for (Int_t nJetsType=0; nJetsType < 1; nJetsType++) {
  for (Int_t nJetsType=0; nJetsType < 3; nJetsType++) {
    
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
      if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

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
      if( lq1*lq2 < 0                 					 ) continue; // cut on opposite-sign leptons

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

      //adhoc cut to "fix" problem with wgamma-star MC sample (there's a weird bump in lep1.pt that needs to be studied in more detail)
      if( lep1->pt() >= 100	    					 ) continue; // cut on leading lepton pt

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
      Double_t FakeLeptonPt = 0;
      Double_t FakeLeptonEta = 0;
      Int_t FakeLeptonPtEtaBin = -1;

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
      if(nFake > 1){
        myWeight = 0.0;
      }
      else if(nFake == 1){
        if(dstype == SmurfTree::data){
          addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = addFR;
          fDecay  	      = 5;
          myWeight	      = add;

        }
        else if(isRealLepton == true || dstype == SmurfTree::wgamma){
          addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
                                 (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
                                 (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

          add = addFR;
          add = add*nPUScaleFactor(fhDPUS4,npu);

          addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
            leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
          add = add*addLepEff;

          add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                                            fabs(lep2->eta()), lep2->pt(), 
                                                            TMath::Abs( lid1), TMath::Abs(lid2));
          fDecay  	       = 5;
          myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;



        }
        else {
          myWeight = 0.0;
        }



 //        if (myWeight > 0) {
//           cout << dstype << " : " << myWeight << " : " << scale1fb*scaleFactorLum << " " << add << " : " 
//                << fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) << " " 
//                << fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection, (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) << " "
//                << endl;
//         }


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

        if (FakeLeptonPt <= 20) {
          if (fabs(FakeLeptonEta) < 1.0) {
            FakeLeptonPtEtaBin = 1;
          } else if (fabs(FakeLeptonEta) < 1.479) {
            FakeLeptonPtEtaBin = 2;
          } else if (fabs(FakeLeptonEta) < 2.5) {
            FakeLeptonPtEtaBin = 3;
          }
        } else {
          if (fabs(FakeLeptonEta) < 1.0) {
            FakeLeptonPtEtaBin = 4;
          } else if (fabs(FakeLeptonEta) < 1.479) {
            FakeLeptonPtEtaBin = 5;
          } else if (fabs(FakeLeptonEta) < 2.5) {
            FakeLeptonPtEtaBin = 6;
          }
        }
    
      }
      else if(dstype == SmurfTree::data) myWeight = 0.0;
      else if(dstype== SmurfTree::dyttDataDriven || dstype == SmurfTree::qcd) {
        myWeight = ZttScaleFactor(nvtx,period,scale1fb)*scaleFactorLum;
      }
      else if(dstype != SmurfTree::data){
        add = 1.0;
        add = add*nPUScaleFactor(fhDPUS4,npu);

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
        if(fDecay == 5) add=add*WJetsMCScaleFactor(); 

        if(dstype == SmurfTree::wgstar) add=add*WGstarScaleFactor();

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
//       if(nJetsType == 2){
//         int centrality = 0;
//         if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
//             (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
//            ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
//             (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
//         passAllCuts = (*jet1+*jet2).M() > 450. &&
//           TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
//           (mH > 200 || dilep->mass() < 100.) &&
//           centrality == 1 &&
//           passJetCut[0]==true;
//       }    

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


        Int_t BkgTypeBin = -1;
        if (fDecay == 5) BkgTypeBin = kFakeLeptonBkg;
        else if (fDecay == 6) BkgTypeBin = kWgamma;
        else if (fDecay == 2) BkgTypeBin = kWZ;
        else BkgTypeBin = kOther;        

        if (type == SmurfTree::mm) {
          SameSignYields[kFakeLepton][kMuMu][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kMuMu][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          
          SameSignPtMax[kFakeLepton][kMuMu][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kMuMu][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kMuMu][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kMuMu][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 

        }
        if (type == SmurfTree::ee) {
          SameSignYields[kFakeLepton][kEleEle][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kEleEle][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

          SameSignPtMax[kFakeLepton][kEleEle][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kEleEle][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kEleEle][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kEleEle][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
        }
        if (type == SmurfTree::em) {
          SameSignYields[kFakeLepton][kEleMu][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kEleMu][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

          SameSignPtMax[kFakeLepton][kEleMu][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kEleMu][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kEleMu][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kEleMu][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
        }
        if (type == SmurfTree::me) {
          SameSignYields[kFakeLepton][kMuEle][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kMuEle][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
          SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

          SameSignPtMax[kFakeLepton][kMuEle][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kMuEle][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kMuEle][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kMuEle][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
        }



        if (FakeLeptonType == 11) {
          if (type == SmurfTree::mm) {
            SameSignYields[kFakeElectron][kMuMu][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kMuMu][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeElectron][kMuMu][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kMuMu][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kMuMu][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kMuMu][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          }
          if (type == SmurfTree::ee) {
            SameSignYields[kFakeElectron][kEleEle][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kEleEle][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeElectron][kEleEle][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kEleEle][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kEleEle][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kEleEle][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          }
          if (type == SmurfTree::em) {
            SameSignYields[kFakeElectron][kEleMu][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kEleMu][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeElectron][kEleMu][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kEleMu][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kEleMu][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kEleMu][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          }
          if (type == SmurfTree::me) {
            SameSignYields[kFakeElectron][kMuEle][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kMuEle][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeElectron][kMuEle][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kMuEle][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kMuEle][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kMuEle][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeElectron][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
           }
        }
        if (FakeLeptonType == 13) {
          if (type == SmurfTree::mm) {
            SameSignYields[kFakeMuon][kMuMu][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kMuMu][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeMuon][kMuMu][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kMuMu][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kMuMu][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kMuMu][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
           }
          if (type == SmurfTree::ee) {
            SameSignYields[kFakeMuon][kEleEle][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kEleEle][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeMuon][kEleEle][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kEleEle][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kEleEle][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kEleEle][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          }
          if (type == SmurfTree::em) {
            SameSignYields[kFakeMuon][kEleMu][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kEleMu][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeMuon][kEleMu][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kEleMu][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kEleMu][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kEleMu][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          }
          if (type == SmurfTree::me) {
            SameSignYields[kFakeMuon][kMuEle][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYields[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kMuEle][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kDifferentFlavor][JetBinIndex][BkgTypeBin] += newWeight*newWeight;
            SameSignYieldsErrSqr[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin] += newWeight*newWeight;

            SameSignPtMax[kFakeMuon][kMuEle][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMax[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep1->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kMuEle][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignPtMin[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(lep2->pt(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kMuEle][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonMass[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->mass(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kMuEle][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kSameFlavor][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
            SameSignDileptonPt[kFakeMuon][kAllFinalStates][JetBinIndex][BkgTypeBin]->Fill(dilep->pt(),newWeight); 
          }
        }
     
      }
    }  


    printf("--- Finished Bgdnal loop\n");




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

    float nDatAcc = 0.0;
    float nDatCut = 0.0;
    float nDatMVA = 0.0;
    for (UInt_t i=0; i<data->GetEntries(); i++) {
    
      data->GetEntry(i);
      if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

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
      bool passJetCut[1] = {Njet3 == nJetsType};

      double minmet = TMath::Min(pmet,pTrackMet);
      bool passMET = minmet > 20. &&
        (minmet > 37.+nvtx/2.0 || type == SmurfTree::em || type == SmurfTree::me);

      bool passNewCuts = true;
      if(lep2->pt() <= 15 && (type == SmurfTree::mm||type == SmurfTree::ee)) passNewCuts = false;
      if(dilep->pt() <= 45) passNewCuts = false;



      // WW Preselection
      if( passJetCut[0]==false                                             ) continue; // select n-jet type events
      if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
      if( lq1*lq2 < 0                 					 ) continue; // cut on opposite-sign leptons
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
      //if( nSoftMuons != 0                                                ) continue; // cut on soft muons veto
      //if( jet1Btag >= 2.1             					 ) continue; // cut on jet1Btag
      //if( jet2Btag >= 2.1             					 ) continue; // cut on jet2Btag
      if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

      bool dPhiDiLepJetCut = true;
      if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
      else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
      if( dPhiDiLepJetCut == false                                         ) continue; // cut on dPhiDiLepJetCut

      //adhoc cut to "fix" problem with wgamma-star MC sample (there's a weird bump in lep1.pt that needs to be studied in more detail)
      if( lep1->pt() >= 100	    					 ) continue; // cut on leading lepton pt


      //****************************************************************************************
      //Find the right bin
      //****************************************************************************************
      Int_t JetBinIndex = -1;
      if(njets == 0) JetBinIndex = kZeroJetBin;
      else if(njets == 1) JetBinIndex = kOneJetBin;
      else JetBinIndex = kVBFBin;



      double myWeight = 1.0;

      double theCutPtMinLow = cutPtMinLow (mH, type);
      bool passAllCuts = dilep->mass()         < theCutMassHigh &&
        mt		     > theCutMTLow &&
        mt		     < theCutMTHigh &&
        lep1->pt()	     > theCutPtMaxLow &&
        lep2->pt()	     > theCutPtMinLow &&
        dPhi*180.0/TMath::Pi()< theCutDeltaphilHigh;
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
        if (type == SmurfTree::mm) {
          SameSignYields[kFakeLepton][kMuMu][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kSameFlavor][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][kData]++;

          SameSignPtMax[kFakeLepton][kMuMu][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMin[kFakeLepton][kMuMu][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignDileptonMass[kFakeLepton][kMuMu][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonPt[kFakeLepton][kMuMu][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->pt()); 
        }
        if (type == SmurfTree::ee) {
          SameSignYields[kFakeLepton][kEleEle][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kSameFlavor][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][kData]++;

          SameSignPtMax[kFakeLepton][kEleEle][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMin[kFakeLepton][kEleEle][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignDileptonMass[kFakeLepton][kEleEle][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonPt[kFakeLepton][kEleEle][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->pt()); 
        }
        if (type == SmurfTree::em) {
          SameSignYields[kFakeLepton][kEleMu][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kSameFlavor][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][kData]++;

          SameSignPtMax[kFakeLepton][kEleMu][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMin[kFakeLepton][kEleMu][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignDileptonMass[kFakeLepton][kEleMu][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonPt[kFakeLepton][kEleMu][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->pt()); 
        }
        if (type == SmurfTree::me) {
          SameSignYields[kFakeLepton][kMuEle][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kSameFlavor][JetBinIndex][kData]++;
          SameSignYields[kFakeLepton][kAllFinalStates][JetBinIndex][kData]++;

          SameSignPtMax[kFakeLepton][kMuEle][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMax[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep1->pt()); 
          SameSignPtMin[kFakeLepton][kMuEle][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignPtMin[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(lep2->pt()); 
          SameSignDileptonMass[kFakeLepton][kMuEle][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonMass[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->mass()); 
          SameSignDileptonPt[kFakeLepton][kMuEle][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kSameFlavor][JetBinIndex][kData]->Fill(dilep->pt()); 
          SameSignDileptonPt[kFakeLepton][kAllFinalStates][JetBinIndex][kData]->Fill(dilep->pt()); 
        }
      }
    }



  } // loop over NJetBins





  //**********************************************************************
  //Print Tables
  //**********************************************************************
  ofstream fResultTexTable("SameSignControlSampleYields.tex");
  char buffer[200];

  fResultTexTable << "\\hline" << endl;
  fResultTexTable << setw(30) << left << "\\multicolumn{2}{|c|}{0-Jet Bin} \\\\" << endl;
  fResultTexTable << "\\hline" << endl;
  fResultTexTable << setw(30) << left << "Background Process"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Yield"
                  << " \\\\"
                  << endl;
  fResultTexTable << "\\hline" << endl;

  
  fResultTexTable << setw(30) << left << "Data"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f",SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kData]);
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;


  fResultTexTable << setw(30) << left << "Fake Lepton Background"
                  << setw(3) << left << "&"  ; 
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg],2)));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;

  fResultTexTable << setw(30) << left << "W+$\\gamma$, W+$\\gamma^{\\ast}$"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma],2)));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;

  fResultTexTable << setw(30) << left << "WZ, ZZ"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ],2)) );
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;

  fResultTexTable << setw(30) << left << "Other"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;

  fResultTexTable << "\\hline" << endl;

  fResultTexTable << setw(30) << left << "Total"
                  << setw(3) << left << "&";     
  sprintf(buffer,"%.1f +/- %.1f",
          SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg] + SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma] + SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ] + SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther],
          TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg],2) + SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma],2) + SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ],2) + SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;
  fResultTexTable << "\\hline" << endl;






  fResultTexTable << "\\hline" << endl;
  fResultTexTable << setw(30) << left << "\\multicolumn{2}{|c|}{1-Jet Bin} \\\\" << endl;
  fResultTexTable << "\\hline" << endl;
  fResultTexTable << setw(30) << left << "Background Process"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Yield"
                  << " \\\\"
                  << endl;
  fResultTexTable << "\\hline" << endl;

  
  fResultTexTable << setw(30) << left << "Data"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f",SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kData]);
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;


  fResultTexTable << setw(30) << left << "Fake Lepton Background"
                  << setw(3) << left << "&"  ; 
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kFakeLeptonBkg], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kFakeLeptonBkg]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;

  fResultTexTable << setw(30) << left << "W+$\\gamma$, W+$\\gamma^{\\ast}$"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kWgamma], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kWgamma]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;

  fResultTexTable << setw(30) << left << "WZ, ZZ"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kWZ],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kWZ]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;


  fResultTexTable << setw(30) << left << "Other"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kOther],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kOther]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;
  fResultTexTable << "\\hline" << endl;

  fResultTexTable << setw(30) << left << "Total"
                  << setw(3) << left << "&";     
  sprintf(buffer,"%.1f +/- %.1f",
          SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kFakeLeptonBkg] + SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kWgamma] + SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kWZ] + SameSignYields[kFakeLepton][kAllFinalStates][kOneJetBin][kOther],
          TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kFakeLeptonBkg] + SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kWgamma] + SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kWZ] + SameSignYieldsErrSqr[kFakeLepton][kAllFinalStates][kOneJetBin][kOther]));
  fResultTexTable << setw(20) << left << buffer;
  fResultTexTable << " \\\\" 
                  << endl;
  fResultTexTable << "\\hline" << endl;



  //**********************************************************************
  //Print Tables
  //**********************************************************************
  ofstream fMoreResultTexTable("SameSignControlSampleYieldsMoreBins.tex");

  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "\\multicolumn{2}{|c|}{0-Jet Bin, ee Final State} \\\\" << endl;
  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "Background Process"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Yield"
                  << " \\\\"
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;

  
  fMoreResultTexTable << setw(30) << left << "Data"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f",SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kData]);
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;


  fMoreResultTexTable << setw(30) << left << "Fake Lepton Background"
                  << setw(3) << left << "&"  ; 
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kFakeLeptonBkg], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kFakeLeptonBkg],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "W+$\\gamma$, W+$\\gamma^{\\ast}$"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWgamma], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWgamma],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "WZ, ZZ"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWZ],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWZ],2)) );
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "Other"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kOther],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << "\\hline" << endl;

  fMoreResultTexTable << setw(30) << left << "Total"
                  << setw(3) << left << "&";     
  sprintf(buffer,"%.1f +/- %.1f",
          SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kFakeLeptonBkg] + SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWgamma] + SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWZ] + SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kOther],
          TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kFakeLeptonBkg],2) + SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWgamma],2) + SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kEleEle][kZeroJetBin][kWZ],2) + SameSignYieldsErrSqr[kFakeLepton][kEleEle][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;




  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "\\multicolumn{2}{|c|}{0-Jet Bin, e$\\mu$ Final State} \\\\" << endl;
  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "Background Process"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Yield"
                  << " \\\\"
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;

  
  fMoreResultTexTable << setw(30) << left << "Data"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f",SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kData]);
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;


  fMoreResultTexTable << setw(30) << left << "Fake Lepton Background"
                  << setw(3) << left << "&"  ; 
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kFakeLeptonBkg], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kFakeLeptonBkg],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "W+$\\gamma$, W+$\\gamma^{\\ast}$"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWgamma], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWgamma],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "WZ, ZZ"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWZ],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWZ],2)) );
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "Other"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kOther],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << "\\hline" << endl;

  fMoreResultTexTable << setw(30) << left << "Total"
                  << setw(3) << left << "&";     
  sprintf(buffer,"%.1f +/- %.1f",
          SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kFakeLeptonBkg] + SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWgamma] + SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWZ] + SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kOther],
          TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kFakeLeptonBkg],2) + SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWgamma],2) + SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kEleMu][kZeroJetBin][kWZ],2) + SameSignYieldsErrSqr[kFakeLepton][kEleMu][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;






  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "\\multicolumn{2}{|c|}{0-Jet Bin, $\\mu$e Final State} \\\\" << endl;
  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "Background Process"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Yield"
                  << " \\\\"
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;

  
  fMoreResultTexTable << setw(30) << left << "Data"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f",SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kData]);
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;


  fMoreResultTexTable << setw(30) << left << "Fake Lepton Background"
                  << setw(3) << left << "&"  ; 
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kFakeLeptonBkg], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kFakeLeptonBkg],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "W+$\\gamma$, W+$\\gamma^{\\ast}$"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWgamma], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWgamma],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "WZ, ZZ"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWZ],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWZ],2)) );
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "Other"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kOther],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << "\\hline" << endl;

  fMoreResultTexTable << setw(30) << left << "Total"
                  << setw(3) << left << "&";     
  sprintf(buffer,"%.1f +/- %.1f",
          SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kFakeLeptonBkg] + SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWgamma] + SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWZ] + SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kOther],
          TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kFakeLeptonBkg],2) + SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWgamma],2) + SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kMuEle][kZeroJetBin][kWZ],2) + SameSignYieldsErrSqr[kFakeLepton][kMuEle][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;



  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "\\multicolumn{2}{|c|}{0-Jet Bin, $\\mu\\mu$ Final State} \\\\" << endl;
  fMoreResultTexTable << "\\hline" << endl;
  fMoreResultTexTable << setw(30) << left << "Background Process"
                  << setw(3) << left << "&"   
                  << setw(20) << left << "Yield"
                  << " \\\\"
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;

  
  fMoreResultTexTable << setw(30) << left << "Data"
                  << setw(3) << left << "&"   ;
  sprintf(buffer,"%.1f",SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kData]);
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;


  fMoreResultTexTable << setw(30) << left << "Fake Lepton Background"
                  << setw(3) << left << "&"  ; 
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kFakeLeptonBkg], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kFakeLeptonBkg],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "W+$\\gamma$, W+$\\gamma^{\\ast}$"
                  << setw(3) << left << "&"   ;

  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWgamma], TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWgamma],2)));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "WZ, ZZ"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWZ],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWZ],2)) );
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << setw(30) << left << "Other"
                  << setw(3) << left << "&";   
  sprintf(buffer,"%.1f +/- %.1f",SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kOther],TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;

  fMoreResultTexTable << "\\hline" << endl;

  fMoreResultTexTable << setw(30) << left << "Total"
                  << setw(3) << left << "&";     
  sprintf(buffer,"%.1f +/- %.1f",
          SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kFakeLeptonBkg] + SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWgamma] + SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWZ] + SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kOther],
          TMath::Sqrt(SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kFakeLeptonBkg] + pow(0.36*SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kFakeLeptonBkg],2) + SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kWgamma] + pow(0.30*SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWgamma],2) + SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kWZ] + pow(0.10*SameSignYields[kFakeLepton][kMuMu][kZeroJetBin][kWZ],2) + SameSignYieldsErrSqr[kFakeLepton][kMuMu][kZeroJetBin][kOther]));
  fMoreResultTexTable << setw(20) << left << buffer;
  fMoreResultTexTable << " \\\\" 
                  << endl;
  fMoreResultTexTable << "\\hline" << endl;



  //**********************************************************************
  //Make Plots
  //**********************************************************************
  vector<string> legends;
  legends.push_back("FakeLepton"); 
  legends.push_back("Wgamma"); 
  legends.push_back("WZ");
  legends.push_back("Other");

  TLegend *tmpLegend = 0;


  //**********************************************************************
  //PtMax
  //**********************************************************************
  tmpLegend = new TLegend(0.73,0.55,0.93,0.90);
  tmpLegend->AddEntry(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg], "FakeLepton", "F");
  tmpLegend->AddEntry(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma], "Wgamma", "F");
  tmpLegend->AddEntry(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ], "WZ", "F");
  tmpLegend->AddEntry(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther], "Other", "F");

  
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillColor(kGray+1);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillColor(kMagenta);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillColor(kAzure-2);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillColor(kAzure-9);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillStyle(1001);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillStyle(1001);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillStyle(1001);
  SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillStyle(1001);

  THStack *SameSignPtMaxStack = new THStack("SameSignPtMax","");
  SameSignPtMaxStack->Add(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]);
  SameSignPtMaxStack->Add(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]);
  SameSignPtMaxStack->Add(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]);
  SameSignPtMaxStack->Add(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]);

  DrawDataBkgHistogram(SameSignPtMax[kFakeLepton][kAllFinalStates][kZeroJetBin][kData], 
                       SameSignPtMaxStack,tmpLegend,
                       "SameSignPtMax",-999, -999,kFALSE); 

  //**********************************************************************
  //PtMin
  //**********************************************************************
  tmpLegend = new TLegend(0.73,0.55,0.93,0.90);
  tmpLegend->AddEntry(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg], "FakeLepton", "F");
  tmpLegend->AddEntry(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma], "Wgamma", "F");
  tmpLegend->AddEntry(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ], "WZ", "F");
  tmpLegend->AddEntry(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther], "Other", "F");

  
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillColor(kGray+1);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillColor(kMagenta);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillColor(kAzure-2);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillColor(kAzure-9);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillStyle(1001);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillStyle(1001);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillStyle(1001);
  SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillStyle(1001);

  THStack *SameSignPtMinStack = new THStack("SameSignPtMin","");
  SameSignPtMinStack->Add(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]);
  SameSignPtMinStack->Add(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]);
  SameSignPtMinStack->Add(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]);
  SameSignPtMinStack->Add(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]);

  DrawDataBkgHistogram(SameSignPtMin[kFakeLepton][kAllFinalStates][kZeroJetBin][kData], 
                       SameSignPtMinStack,tmpLegend,
                       "SameSignPtMin",-999, -999,kFALSE); 

  //**********************************************************************
  //DileptonMass
  //**********************************************************************
  tmpLegend = new TLegend(0.73,0.55,0.93,0.90);
  tmpLegend->AddEntry(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg], "FakeLepton", "F");
  tmpLegend->AddEntry(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma], "Wgamma", "F");
  tmpLegend->AddEntry(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ], "WZ", "F");
  tmpLegend->AddEntry(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther], "Other", "F");

  
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillColor(kGray+1);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillColor(kMagenta);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillColor(kAzure-2);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillColor(kAzure-9);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillStyle(1001);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillStyle(1001);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillStyle(1001);
  SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillStyle(1001);

  THStack *SameSignDileptonMassStack = new THStack("SameSignDileptonMass","");
  SameSignDileptonMassStack->Add(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]);
  SameSignDileptonMassStack->Add(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]);
  SameSignDileptonMassStack->Add(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]);
  SameSignDileptonMassStack->Add(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]);

  DrawDataBkgHistogram(SameSignDileptonMass[kFakeLepton][kAllFinalStates][kZeroJetBin][kData], 
                       SameSignDileptonMassStack,tmpLegend,
                       "SameSignDileptonMass",-999, -999,kFALSE); 

  //**********************************************************************
  //DileptonPt
  //**********************************************************************
  tmpLegend = new TLegend(0.73,0.55,0.93,0.90);
  tmpLegend->AddEntry(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg], "FakeLepton", "F");
  tmpLegend->AddEntry(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma], "Wgamma", "F");
  tmpLegend->AddEntry(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ], "WZ", "F");
  tmpLegend->AddEntry(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther], "Other", "F");

  
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillColor(kGray+1);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillColor(kMagenta);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillColor(kAzure-2);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillColor(kAzure-9);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]->SetFillStyle(1001);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]->SetFillStyle(1001);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]->SetFillStyle(1001);
  SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]->SetFillStyle(1001);

  THStack *SameSignDileptonPtStack = new THStack("SameSignDileptonPt","");
  SameSignDileptonPtStack->Add(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kFakeLeptonBkg]);
  SameSignDileptonPtStack->Add(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWgamma]);
  SameSignDileptonPtStack->Add(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kWZ]);
  SameSignDileptonPtStack->Add(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kOther]);

  DrawDataBkgHistogram(SameSignDileptonPt[kFakeLepton][kAllFinalStates][kZeroJetBin][kData], 
                       SameSignDileptonPtStack,tmpLegend,
                       "SameSignDileptonPt",-999, -999,kFALSE); 



    
}

