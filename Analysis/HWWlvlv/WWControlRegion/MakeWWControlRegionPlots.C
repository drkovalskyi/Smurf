//root -l -q -b WWControlRegion/MakeWWControlRegionPlots.C+\(0,\"default\",\"data/backgroundC_skim2.root\"\,\"data/data_2l_skim2.root\",4,13\) ;

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
#include "/home/ceballos/releases/CMSSW_4_2_8_patch4/src/Smurf/Core/SmurfTree.h"
#include "factors.h"
#include "HiggsQCDScaleSystematics.h"
#include "PSUESystematics.h"
#include "PDFgHHSystematics.h"
#include "/home/ceballos/releases/CMSSW_4_2_8_patch4/src/Smurf/Core/LeptonScaleLookup.h"
#include "DYBkgScaleFactors.h"
#include "TopBkgScaleFactors.h"
#include "WWBkgScaleFactors.h"
#include "OtherBkgScaleFactors.h"
#include "HWWCuts.h"
#include "StandardPlot.C"

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
// PlotHiggsRes
//------------------------------------------------------------------------------
void FillWWControlRegionPlots
(
 UInt_t  nJetsType   	 = 0,
 TString outputLabel     = "default",
 TString bgdInputFile    = "data/inputNtuple-data-standard-HBCK_WW_2l-train.root",
 TString datInputFile    = "data/data-train.root",
 Int_t   wwDecay         = 0,
 int period              = 13,
 int  TheVerboseLevel    = 0
 )
{

  verboseLevel = TheVerboseLevel;

  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  unsigned int patternTopTag = SmurfTree::TopTag;
  float dilmass_cut = DileptonMassPreselectionCut(0);

  cout << "dilcut: " <<  dilmass_cut << endl;

  char finalStateName[10];
  sprintf(finalStateName,"ll");
  if	 (wwDecay == 0) sprintf(finalStateName,"mm");
  else if(wwDecay == 1) sprintf(finalStateName,"me");
  else if(wwDecay == 2) sprintf(finalStateName,"em");
  else if(wwDecay == 3) sprintf(finalStateName,"ee");
  else if(wwDecay == 5) sprintf(finalStateName,"sf");
  else if(wwDecay == 6) sprintf(finalStateName,"of");

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

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
  chdata->Add(datFile1);

  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;


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
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Full2011.root";
    scaleFactorLum     = 4.7;minRun =      0;maxRun = 999999;
  }
  else if(period == 3){ // Full2011
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root";
    scaleFactorLum     = 4.7;minRun =      0;maxRun = 999999;
  }
  else if(period == 13){ // Full2011
    effPath  = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/efficiency_results_MVAIDIsoCombinedDetIsoSameSigWP_Full2011.root";
    fakePath = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/FakeRates_MVAIDIsoCombinedDetIsoSameSigWP.root";
    puPath   = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root";
    scaleFactorLum     = 4.7;minRun =      0;maxRun = 999999;
  }
  else if(period == 14){ // Full2011
    effPath  = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/efficiency_results_MVAIDIsoCombinedSameSigWP_Full2011.root";
    fakePath = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/FakeRates_MVAIDIsoCombinedSameSigWP.root";
    puPath   = "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root";
    scaleFactorLum     = 4.7;minRun =      0;maxRun = 999999;
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

  TFile *fPUS4File = TFile::Open(puPath.Data());
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

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



  //***********************************************************************************************
  //Define Histograms & Yields
  //***********************************************************************************************

  //----------------------------------------------------------------------------
  const int nHist = 8;
  int    nBinHis[nHist]  = {    200,     200,            200,        180,  200,          200, 200, 200};
  double minHis[nHist]   = {      0,       0,              0,          0,    0,            0,   0,   0};
  double maxHis[nHist]   = {    100,     100,            200,        180,  200,          100, 7.0, 900};
  string HistName[nHist] = {"PtMax", "PtMin", "DileptonMass", "DeltaPhi", "Mt", "DileptonPt", 
                            "DeltaEtaVBFJets", "DijetMassVBFJets"};
  string HistAxisLabel[nHist] = {"Leading Lepton p_{T} [GeV/c]", "Trailing Lepton p_{T} [GeV/c]",
                                 "M_{ll} [GeV/c^{2}]", "#Delta#phi(l,l)",
                                 "m_{T}^{ll E_{T}^{miss}}",
                                 "p_{T}^{ll} [GeV/c]",
                                 "#Delta#eta(j,j)",
                                 "M_{jj} [GeV/c^{2}]" };

  const int nBkgProcesses = 5;
  string bkgProcessLabel[nBkgProcesses] = {"WW","WZZZ","Top","DYeemm","WJets"};
  string JetBinLabel[3] = {"0Jet","1Jet","2Jet"};

  //----------------------------------------------------------------------------
  TH1D* dataDistributions[nHist];
  TH1D* bkgDistributions[nHist][nBkgProcesses];
  for(int i=0; i<nHist; i++) {
    dataDistributions[i] = new TH1D( ("dataDistributions_" + HistName[i] + "_" + JetBinLabel[nJetsType] + "_" + string(finalStateName)).c_str(), 
                                       (";"+HistAxisLabel[i]+";Number of Events").c_str(), 
                                       nBinHis[i], minHis[i], maxHis[i]);
    dataDistributions[i]->Sumw2();

    for(int j=0; j<nBkgProcesses; j++){
      bkgDistributions[i][j] = new TH1D( ("bkgDistributions_" + HistName[i]+"_"+bkgProcessLabel[j]+ "_" + JetBinLabel[nJetsType] + "_" + string(finalStateName)).c_str(), 
                                      (";"+HistAxisLabel[i]+";Number of Events").c_str(), 
                                      nBinHis[i], minHis[i], maxHis[i]);
      bkgDistributions[i][j]->Sumw2();
    }
  }




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
  Float_t         higgsPt = -999;
  Float_t         sfWeightPU;
  Float_t         sfWeightEff;
  Float_t         sfWeightTrig;
  Float_t         sfWeightHPt;




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
  background->SetBranchAddress( "lep1McId"      , &lep1McId	  );
  background->SetBranchAddress( "lep2McId"      , &lep2McId	  );
  background->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId );
  background->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId );

  float nSigAcc[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigCut[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigMVA[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigEAcc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigECut[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float nSigEMVA[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  const int nChan = 8;
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
    //  - VBF jets must have |eta| < 4.5 (instead of 5.0)
    //----------------------------------------------------------------------------
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
 
    //----------------------------------------------------------------------------
    //For Jet Energy scale systematics
    //----------------------------------------------------------------------------
    bool passJetCut[3] = {Njet3 == nJetsType, false, false};

    double minmet = TMath::Min(pmet,pTrackMet);
    bool passMET = minmet > 20. &&
                  (minmet > 37.+nvtx/2.0 || type == SmurfTree::em || type == SmurfTree::me);

    bool passNewCuts = true;
    if(lep2->pt() <= 15 && (type == SmurfTree::mm||type == SmurfTree::ee)) passNewCuts = false;
    if(dilep->pt() <= 45) passNewCuts = false;

    //----------------------------------------------------------------------------
    // WW Preselection
    //----------------------------------------------------------------------------
    bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
                        ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
                         dstype != SmurfTree::data;





    if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
    if( passJetCut[0]==false                                             ) continue; // select n-jet type events
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
    else if(dstype == SmurfTree::ggww   	 ) fDecay = 1;
    else if(dstype == SmurfTree::wgamma 	 ) fDecay = 6;
    else if(dstype == SmurfTree::wgstar 	 ) fDecay = 6;
    else if(dstype == SmurfTree::data   	 ) fDecay = 5;
    else if(dstype == SmurfTree::dyttDataDriven  ) fDecay = 7;
    else if(dstype == SmurfTree::qcd             ) fDecay = 7;
    else                                 {printf("bad dstype: %d\n",dstype); assert(0);}

    //----------------------------------------------------------------------------
    // Categorize by Background type for histograms
    //----------------------------------------------------------------------------
    int BkgType = 0;
    if(dstype == SmurfTree::qqww  ) BkgType = 0;
    else if(dstype == SmurfTree::ggww  ) BkgType = 0;
    else if(dstype == SmurfTree::wz    ) BkgType = 1;
    else if(dstype == SmurfTree::zz    ) BkgType = 1;
    else if(dstype == SmurfTree::tw    ) BkgType = 2;
    else if(dstype == SmurfTree::ttbar ) BkgType = 2;
    else if(dstype == SmurfTree::dyee  ) BkgType = 3;
    else if(dstype == SmurfTree::dymm  ) BkgType = 3;
    else if(dstype == SmurfTree::dytt  ) BkgType = 3;
    else if(dstype == SmurfTree::wjets ) BkgType = 4;
    else if(dstype == SmurfTree::data  ) BkgType = 4;
    else if(dstype == SmurfTree::wgamma) BkgType = 4;


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

    //if(nJetsType == 0) bdtg = (bdtg+1.2)*TMath::Max(TMath::Min((double)knn_wjets,0.99999),0.00001)-1.0;
    //else               bdtg = (bdtg+1.0)*TMath::Max(TMath::Min((double)bdtg_wjets+1,1.99999),0.00001)/2.0-1.0;

    //bdtg = TMath::Max(TMath::Min((bdtg+1.0)/2.0,1.0),0.0)*
    //       TMath::Max(TMath::Min((bdtg+1.0)/2.0,1.0),0.0)+
    //	   TMath::Max(TMath::Min((bdtg_wjets+1.0)/2.0,1.0),0.0)*
    //       TMath::Max(TMath::Min((bdtg_wjets+1.0)/2.0,1.0),0.0)-1.0;
    // bdtg = (bdtg+bdtg_wjets)/2.0;
    //if(nJetsType == 0 && bdtg_wjets <= 0.8) continue;
    //if(nJetsType == 1 && bdtg_wjets <= 0.9) continue;
    //bdtg = TMath::Min(knn-0.5,0.999)/0.5;
    //bdtg = (dilep->mass()-12.0)/(dilmass_cut-12.0);
    //if(bdtg<=0) bdtg = 0.001; if(bdtg>=1) bdtg = 0.999; bdtg = (bdtg-0.5)*2.0;

    //----------------------------------------------------------------------------
    // Define event weights    
    //----------------------------------------------------------------------------
    double myWeight = 1.0;
    double add      = 1.0;

    //----------------------------------------------------------------------------
    // First classify event into tight+tight, tight+fail, fail+fail 
    //----------------------------------------------------------------------------
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

      //don't do pileup reweighting for DYtautau embedded sample
      if (!(dstype== SmurfTree::dyttDataDriven || dstype == SmurfTree::qcd)) {
        add = add*nPUScaleFactor(fhDPUS4,npu);
      }

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
      //----------------------------------------------------------------------------
      if((fDecay == 0 || fDecay == 1) ){     
        if(njets == 0) add=add*WWBkgScaleFactorCutBased(115,0); 
        else	       add=add*WWBkgScaleFactorCutBased(115,1); 
      }
      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;
      myWeight = scale1fb*scaleFactorLum*add;
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
    double theCutPtMinLow = cutPtMinLow (0, type);
    bool passAllCuts = (lep2->pt() > theCutPtMinLow);

//     //----------------------------------------------------------------------------
//     // VBF selection cuts for 2-Jet bin
//     //----------------------------------------------------------------------------
//     if(nJetsType == 2){
//       int centrality = 0;
//       if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
//           (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
//          ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
//           (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
//       passAllCuts = (*jet1+*jet2).M() > 450. &&
//                     TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
// 		    (mH > 200 || dilep->mass() < 100.) &&
// 		    centrality == 1 &&
// 		    passJetCut[0]==true;
//     }


    //----------------------------------------------------------------------------
    // Add bkg yields for Cut-Based analysis
    //----------------------------------------------------------------------------
    if(passAllCuts == true) {
      double newWeight = myWeight;

      //----------------------------------------------------------------------------
      // For Cut-Based Analysis, use mass dependant  
      // scale factors for DY Bkg. We originally applied DY scale factors for the
      // WW pre-selection, and correct it by the ratio here.
      //----------------------------------------------------------------------------      
      if((dstype == SmurfTree::dymm || dstype == SmurfTree::dyee) &&
         (type   == SmurfTree::mm   || type   == SmurfTree::ee)){
	DYXS[1] += newWeight; 
      }
      else if(fDecay == 4){
	VVXS[1] += newWeight; 
      }
      nBgdCut  = nBgdCut  + newWeight;
      nBgdECut = nBgdECut + newWeight*newWeight;
      nBgdCutDecays[fDecay]  = nBgdCutDecays[fDecay]  + newWeight;
      nBgdECutDecays[fDecay] = nBgdECutDecays[fDecay] + newWeight*newWeight;

      //----------------------------------------------------------------------------      
      // Make WW Control Region Plots
      //----------------------------------------------------------------------------
      bkgDistributions[0][BkgType]->Fill(lep1->pt(), myWeight);
      bkgDistributions[1][BkgType]->Fill(lep2->pt(), myWeight);
      bkgDistributions[2][BkgType]->Fill(dilep->mass(), myWeight);
      bkgDistributions[3][BkgType]->Fill(dPhi*180.0/TMath::Pi(), myWeight);
      bkgDistributions[4][BkgType]->Fill(mt, myWeight);
      bkgDistributions[5][BkgType]->Fill(dilep->pt(), myWeight);


    }



  } //end loop over bkg events

  nBgdEAcc = sqrt(nBgdEAcc);
  nBgdECut = sqrt(nBgdECut);
  nBgdEMVA = sqrt(nBgdEMVA);
  printf("--- Finished Bgdnal loop\n");



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
  data->SetBranchAddress( "sfWeightPU"     , &sfWeightPU     );
  data->SetBranchAddress( "sfWeightEff"    , &sfWeightEff    );
  data->SetBranchAddress( "sfWeightTrig"   , &sfWeightTrig   );
  data->SetBranchAddress( "sfWeightHPt"    , &sfWeightHPt   );

  float nDatAcc = 0.0;
  float nDatCut = 0.0;
  float nDatMVA = 0.0;
  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

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
    //  - VBF jets must have |eta| < 4.5 (instead of 5.0)
    //----------------------------------------------------------------------------
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


    //----------------------------------------------------------------------------
    // WW Preselection
    //----------------------------------------------------------------------------
    if( passJetCut[0]==false                                             ) continue; // select n-jet type events
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
    //if( nSoftMuons != 0                                                ) continue; // cut on soft muons veto
    //if( jet1Btag >= 2.1             					 ) continue; // cut on jet1Btag
    //if( jet2Btag >= 2.1             					 ) continue; // cut on jet2Btag
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJetCut = true;
    if(njets <= 1) dPhiDiLepJetCut = jet1->pt() <= 15. || dPhiDiLepJet1*180.0/TMath::Pi() < 165.         || type == SmurfTree::em || type == SmurfTree::me;
    else           dPhiDiLepJetCut = DeltaPhi((*jet1+*jet2).Phi(),dilep->phi())*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me;
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

    //if(nJetsType == 0) bdtg = (bdtg+1.2)*TMath::Max(TMath::Min((double)knn_wjets,0.99999),0.00001)-1.0;
    //else               bdtg = (bdtg+1.0)*TMath::Max(TMath::Min((double)bdtg_wjets+1,1.99999),0.00001)/2.0-1.0;

    //bdtg = TMath::Max(TMath::Min((bdtg+1.0)/2.0,1.0),0.0)*
    //       TMath::Max(TMath::Min((bdtg+1.0)/2.0,1.0),0.0)+
    //	   TMath::Max(TMath::Min((bdtg_wjets+1.0)/2.0,1.0),0.0)*
    //       TMath::Max(TMath::Min((bdtg_wjets+1.0)/2.0,1.0),0.0)-1.0;
    // bdtg = (bdtg+bdtg_wjets)/2.0;
    //if(nJetsType == 0 && bdtg_wjets <= 0.8) continue;
    //if(nJetsType == 1 && bdtg_wjets <= 0.9) continue;
    //bdtg = TMath::Min(knn-0.5,0.999)/0.5;
    //bdtg = (dilep->mass()-12.0)/(dilmass_cut-12.0);
    //if(bdtg<=0) bdtg = 0.001; if(bdtg>=1) bdtg = 0.999; bdtg = (bdtg-0.5)*2.0;


    double myWeight = 1.0;

    //----------------------------------------------------------------------------
    //
    // Higgs Signal Selection Cuts
    //
    //----------------------------------------------------------------------------
    double theCutPtMinLow = cutPtMinLow (0, type);
    bool passAllCuts = lep2->pt()	     > theCutPtMinLow ;

//     //----------------------------------------------------------------------------
//     // VBF selection cuts for 2-Jet bin
//     //----------------------------------------------------------------------------
//     if(nJetsType == 2){
//       int centrality = 0;
//       if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
//           (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
//          ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
//           (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
//       passAllCuts = (*jet1+*jet2).M() > 450. &&
//         TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
//         (mH > 200 || dilep->mass()< 100.) &&
//         centrality == 1;
//     }

    //----------------------------------------------------------------------------
    // Data yields for Cut-Based analysis
    //----------------------------------------------------------------------------
    if(passAllCuts == true) {
      nDatCut = nDatCut + myWeight;

      dataDistributions[0]->Fill(lep1->pt(), myWeight);
      dataDistributions[1]->Fill(lep2->pt(), myWeight);
      dataDistributions[2]->Fill(dilep->mass(), myWeight);
      dataDistributions[3]->Fill(dPhi*180.0/TMath::Pi(), myWeight);
      dataDistributions[4]->Fill(mt, myWeight);
      dataDistributions[5]->Fill(dilep->pt(), myWeight);
    }

  } //end loop over data events



  //****************************************************************************
  //
  // Print Summary Information
  //
  //****************************************************************************

  printf("--- Finished Data loop\n");
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

  //****************************************************************************
  // Compute Total Systematic Uncertainty of the Background
  //****************************************************************************
  float BkgSystematicUncertainty[nChan]  = {0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=0; i<nChan; i++) {
    if (i==0 || i==1) BkgSystematicUncertainty[i] = WWBkgScaleFactorKappaCutBased(115,nJetsType) - 1.0;
    if (i==2) BkgSystematicUncertainty[i] = 0.05;
    if (i==3) BkgSystematicUncertainty[i] = TopBkgScaleFactorKappa(nJetsType) - 1.0;


    if (i==4) BkgSystematicUncertainty[i] = DYBkgScaleFactorKappa(0,nJetsType) - 1.0;
    if (i==5) BkgSystematicUncertainty[i] = 0.36;
    if (i==6) BkgSystematicUncertainty[i] = 0.30;
    if (i==7) BkgSystematicUncertainty[i] = 0.20;
  }

  Double_t nBgdCutSystematicUncertaintySqr = 0;
  for(int i=0; i<nChan; i++) {
    Double_t tempSystematicUncertainty = 0;        
    nBgdCutSystematicUncertaintySqr += pow(nBgdCutDecays[i]*BkgSystematicUncertainty[i],2);
  }
  Double_t nBgdCutSystematicUncertainty = sqrt(nBgdCutSystematicUncertaintySqr);
  Double_t nBgdTotalUncertainty = sqrt(pow(nBgdCutSystematicUncertainty,2) + pow(nBgdECut,2));

  //****************************************************************************
  //
  // Save Distribution Plots
  //
  //****************************************************************************
  TFile *fileOutput = new TFile("WWControlRegionPlots.root", "UPDATE");
  for(int i=0; i<nHist; i++) {
    fileOutput->WriteTObject(dataDistributions[i], dataDistributions[i]->GetName(), "WriteDelete");
    for(int j=0; j<nBkgProcesses; j++){
      fileOutput->WriteTObject(bkgDistributions[i][j], bkgDistributions[i][j]->GetName(), "WriteDelete");
    }
  }
  fileOutput->Close();
 

  //****************************************************************************
  //
  // Print Summary Tex Table
  //
  //****************************************************************************
  ofstream fResultTexTable("WWControlRegionTable_" + outputLabel + ".tex");
  char buffer[200];

  fResultTexTable << setw(10) << left << " "
                  << setw(3) << left << "&"
                  << setw(25) << left << "data"
                  << setw(3) << left << "&"
                  << setw(25) << left << "All Bkg"
                  << setw(3) << left << "&"
                  << setw(25) << left << "$qq \\to \\WW$"
                  << setw(3) << left << "&"
                  << setw(25) << left << "$gg \\to \\WW$"
                  << setw(3) << left << "&"
                  << setw(25) << left << "$\\ttbar+tW$"
                  << setw(3) << left << "&"
                  << setw(25) << left << "$\\Wjets$"
                  << " \\\\"
                  << endl;
                  
  if (nJetsType == 0) {
    fResultTexTable << setw(10) << left << "0-jet" ;
  } else if (nJetsType == 1) {
    fResultTexTable << setw(10) << left << "1-jet" ;
  } else if (nJetsType == 2) {
    fResultTexTable << setw(10) << left << "2-jet" ;
  } else {
    cout << "jet bin not recognized\n";
    return;
  }

  //Data Yields
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << nDatCut ;

  //Total Bkg Yields
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCut, nBgdTotalUncertainty);
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;

  //qqWW Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[0], sqrt(pow(nBgdCutDecays[0]*BkgSystematicUncertainty[0],2) + pow(nBgdECutDecays[0],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;
  //ggWW Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[1], sqrt(pow(nBgdCutDecays[1]*BkgSystematicUncertainty[1],2) + pow(nBgdECutDecays[1],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;
  //top Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[3], sqrt(pow(nBgdCutDecays[3]*BkgSystematicUncertainty[3],2) + pow(nBgdECutDecays[3],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;

  //W+jets Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[5], sqrt(pow(nBgdCutDecays[5]*BkgSystematicUncertainty[5],2) + pow(nBgdECutDecays[5],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;

  fResultTexTable << "\\\\" << endl;



  fResultTexTable << setw(10) << left << " "
                  << setw(3) << left << "&"
                  << setw(25) << left << "$WZ$/$ZZ$"
                  << setw(3) << left << "&"
                  << setw(25) << left << "$\\dyll$"
                  << setw(3) << left << "&"
                  << setw(25) << left << "$W+\\gamma$"
                  << setw(3) << left << "&"
                  << setw(25) << left << "\\dytt"              
                  << " \\\\"
                  << endl;
  if (nJetsType == 0) {
    fResultTexTable << setw(10) << left << "0-jet" ;
  } else if (nJetsType == 1) {
    fResultTexTable << setw(10) << left << "1-jet" ;
  } else if (nJetsType == 2) {
    fResultTexTable << setw(10) << left << "2-jet" ;
  } else {
    cout << "jet bin not recognized\n";
    return;
  }

  //VV Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[2], sqrt(pow(nBgdCutDecays[2]*BkgSystematicUncertainty[2],2) + pow(nBgdECutDecays[2],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;

  //DY Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[4], sqrt(pow(nBgdCutDecays[4]*BkgSystematicUncertainty[4],2) + pow(nBgdECutDecays[4],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;


  //W+Gamma Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[6], sqrt(pow(nBgdCutDecays[6]*BkgSystematicUncertainty[6],2) + pow(nBgdECutDecays[6],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;

  //Ztautau Bkg
  sprintf(buffer,"$%.2f \\pm %.2f$", nBgdCutDecays[7], sqrt(pow(nBgdCutDecays[7]*BkgSystematicUncertainty[7],2) + pow(nBgdECutDecays[7],2)));
  fResultTexTable << setw(3) << left << "&"
                  << setw(25) << left << buffer;
  fResultTexTable << "\\\\" << endl;

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



void MakePlot(int nsel = 0, int ReBin = 1, 
              char XTitle[300] = "N_{jets}", 
              char units[300] = "", 
              string histname = "", 
              string jetbinLabel = "",
              string finalstateLabel = "",
              char outputName[300] = "njets",
              bool isLogY = false, int MassH = 160, double lumi = 1.55, Double_t lowCutValue = -999, Double_t highCutValue = -999) {

    gROOT->LoadMacro("/home/sixie/CMSSW_analysis/src/EWKAna/Thesis/Selection/StandardPlot.C");

    TFile* file = new TFile("WWControlRegionPlots.root", "read");

    StandardPlot myPlot;
    myPlot.setLumi(lumi);
    myPlot.setLabel(XTitle);
    myPlot.addLabel("");
    myPlot.setUnits(units);

    TH1F* hWW    = (TH1F*)file->Get(("bkgDistributions_" + histname + "_WW" + "_" + jetbinLabel + "_" + finalstateLabel).c_str());
    TH1F* hVV    = (TH1F*)file->Get(("bkgDistributions_" + histname + "_WZZZ" + "_" + jetbinLabel+ "_" + finalstateLabel).c_str());
    TH1F* hTop   = (TH1F*)file->Get(("bkgDistributions_" + histname + "_Top" + "_" + jetbinLabel+ "_" + finalstateLabel).c_str());
    TH1F* hZJets = (TH1F*)file->Get(("bkgDistributions_" + histname + "_DYeemm" + "_" + jetbinLabel+ "_" + finalstateLabel).c_str());
    TH1F* hWJets = (TH1F*)file->Get(("bkgDistributions_" + histname + "_WJets" + "_" + jetbinLabel+ "_" + finalstateLabel).c_str());
    assert(hWW);
    assert(hVV);
    assert(hTop);
    assert(hZJets);
    assert(hWJets);
    
    if(hWW->GetSumOfWeights(   ) > 0) myPlot.setMCHist(iWW,    (TH1F*)hWW   ->Clone("hWW"));
    if(hVV->GetSumOfWeights()    > 0) myPlot.setMCHist(iVV,    (TH1F*)hVV   ->Clone("hVV")); 
    if(hTop->GetSumOfWeights()   > 0) myPlot.setMCHist(iTop,   (TH1F*)hTop  ->Clone("hTop"));
    if(hZJets->GetSumOfWeights() > 0) myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
    if(hWJets->GetSumOfWeights() > 0) myPlot.setMCHist(iWJets, (TH1F*)hWJets->Clone("hWJets"));
    
    TH1F* hData  = (TH1F*)file->Get(("dataDistributions_" + histname + "_" + jetbinLabel + "_" + finalstateLabel).c_str());
    assert(hData);
    myPlot.setDataHist(hData);

    if (lowCutValue != -999) myPlot.setLowCutValue(lowCutValue);
    if (highCutValue != -999) myPlot.setHighCutValue(highCutValue);
    myPlot.Draw(outputName, isLogY, ReBin );  // Can pass a rebin 

    return;

}





void MakeWWControlRegionPlots() {

  FillWWControlRegionPlots(0,"0Jet_All","data/ntuples_115train_0jets_backgroundC_skim2.root","data/ntuples_115train_0jets_data_2l_skim2.root",4,2);
  FillWWControlRegionPlots(1,"1Jet_All","data/ntuples_115train_0jets_backgroundC_skim2.root","data/ntuples_115train_0jets_data_2l_skim2.root",4,2);
  FillWWControlRegionPlots(0,"0Jet_SF","data/ntuples_115train_0jets_backgroundC_skim2.root","data/ntuples_115train_0jets_data_2l_skim2.root",5,2);
  FillWWControlRegionPlots(1,"1Jet_SF","data/ntuples_115train_0jets_backgroundC_skim2.root","data/ntuples_115train_0jets_data_2l_skim2.root",5,2);
  FillWWControlRegionPlots(0,"0Jet_DF","data/ntuples_115train_0jets_backgroundC_skim2.root","data/ntuples_115train_0jets_data_2l_skim2.root",6,2);
  FillWWControlRegionPlots(1,"1Jet_DF","data/ntuples_115train_0jets_backgroundC_skim2.root","data/ntuples_115train_0jets_data_2l_skim2.root",6,2);

//     MakePlot(0, 4, "p_{T}^{Leading Lepton}","GeV/c","PtMax", "Selection_HWW160Selection_PtMax",0,160,1.55, 30, -999 );
//     MakePlot(0, 4, "p_{T}^{Trailing Lepton}","GeV/c","PtMin", "Selection_HWW160Selection_PtMin",0,160,1.55, 25, -999 );
//     MakePlot(0, 4, "#Delta#phi(l,l)","degrees","DeltaPhi", "Selection_HWW160Selection_DeltaPhi",0,160,1.55 , -999, 60);
//     MakePlot(0, 4, "m_{T}^{ll E_{T}^{miss}}","GeV/c^{2}","Mt", "Selection_HWW160Selection_Mt",0,160,1.55 , 90, 160);
//     MakePlot(0, 4, "p_{T}^{ll}","GeV/c","DileptonPt", "Selection_HWW160Selection_DileptonPt",0,160,1.55);


  MakePlot(0, 4, "M_{ll}","GeV/c^{2}","DileptonMass","0Jet","sf","DileptonMass_0Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "M_{ll}","GeV/c^{2}","DileptonMass","0Jet","of","DileptonMass_0Jet_OF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "M_{ll}","GeV/c^{2}","DileptonMass","1Jet","sf","DileptonMass_1Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "M_{ll}","GeV/c^{2}","DileptonMass","1Jet","of","DileptonMass_1Jet_OF",0,160,4.7, -999, -999 );

  MakePlot(0, 4, "p_{T}^{Trailing Lepton}","GeV/c","PtMin","0Jet","sf","PtMin_0Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "p_{T}^{Trailing Lepton}","GeV/c","PtMin","0Jet","of","PtMin_0Jet_OF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "p_{T}^{Trailing Lepton}","GeV/c","PtMin","1Jet","sf","PtMin_1Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "p_{T}^{Trailing Lepton}","GeV/c","PtMin","1Jet","of","PtMin_1Jet_OF",0,160,4.7, -999, -999 );

  MakePlot(0, 4, "p_{T}^{Leading Lepton}","GeV/c","PtMax","0Jet","sf","PtMax_0Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "p_{T}^{Leading Lepton}","GeV/c","PtMax","0Jet","of","PtMax_0Jet_OF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "p_{T}^{Leading Lepton}","GeV/c","PtMax","1Jet","sf","PtMax_1Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "p_{T}^{Leading Lepton}","GeV/c","PtMax","1Jet","of","PtMax_1Jet_OF",0,160,4.7, -999, -999 );

  MakePlot(0, 4, "#Delta#phi(l,l)","degrees","DeltaPhi","0Jet","sf","DeltaPhi_0Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "#Delta#phi(l,l)","degrees","DeltaPhi","0Jet","of","DeltaPhi_0Jet_OF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "#Delta#phi(l,l)","degrees","DeltaPhi","1Jet","sf","DeltaPhi_1Jet_SF",0,160,4.7, -999, -999 );
  MakePlot(0, 4, "#Delta#phi(l,l)","degrees","DeltaPhi","1Jet","of","DeltaPhi_1Jet_OF",0,160,4.7, -999, -999 );


}
