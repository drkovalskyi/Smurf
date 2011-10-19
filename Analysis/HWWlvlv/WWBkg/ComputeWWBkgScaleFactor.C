#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
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
#include "Smurf/Core/HWWCuts.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h"
#include "Smurf/Analysis/HWWlvlv/TopBkgScaleFactors.h"

//------------------------------------------------------------------------------
// optimalCuts_42x
//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void ComputeWWBkgScaleFactor (
  char* bgdInputFile    = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/backgroundA_skim1.root",
  char* dataInputFile   = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/data_2l_skim1.root"
)
{

// ***********************************************************************************************
// EPS Data : Cross-checked with Giuseppe's macro
// ***********************************************************************************************
//   TChain *chdata = new TChain("tree");
//   chdata->Add("/data/smurf/data/EPS/tas/data-met20-1092ipb.root");
//   TTree *data     = (TTree*) chdata;

//   TChain *chbackground = new TChain("tree");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/qqww.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/ggww.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/ttbar.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/tw.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/dymm.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/dyee.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/dytt.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/zz.root");
//   chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/wz.root");
//   chbackground->Add("/data/smurf/data/EPS/tas/data-met20-1092ipb.root");
//   TTree *background = (TTree*) chbackground;
//   double scaleFactorLum = 1.092;
  

// ***********************************************************************************************
// Full Run2011A Data 
// ***********************************************************************************************
  TChain *chdata = new TChain("tree");
  chdata->Add(dataInputFile);
  TTree *data     = (TTree*) chdata;

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdInputFile);
  TTree *background = (TTree*) chbackground;
  
  double scaleFactorLum = 2.121;

  bool isOldAna = true;
  TString effPath  = "/data/smurf/data/LP2011/auxiliar/efficiency_results_v6_42x.root";
  TString fakePath = "/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.LP2011.root";
  if(isOldAna == false){
    effPath  = "/data/smurf/data/Winter11/auxiliar/efficiency_results_SmurfBDTGWithIPInfoElectrons.root";
    fakePath = "/data/smurf/data/Winter11/auxiliar/FakeRates_SmurfMVAWithIPInfoElectron.root";
  }
  // ***********************************************************************************************
  // Load Auxiliary Files
  // ***********************************************************************************************

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

  TFile *fPUS4File = TFile::Open("/data/smurf/data/LP2011/auxiliar/puWeights_PU4_68mb.root");
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;




  //***********************************************************************************************
  //Define Histograms & Yields
  //***********************************************************************************************

   const Int_t nmass = 9;
  const Double_t mH[nmass] = {115,120,130,140,150,160,170,180,190};  
  //scale factors for 4 classes of selections
  enum { kCutBasedZeroJet, kCutBasedOneJet, kMVAZeroJet, kMVAOneJet };
    

  vector<vector<Double_t> > Yield_WWControlRegion_Data, Yield_WWControlRegion_Bkg, Yield_WWControlRegion_WWMC; 
  vector<vector<Double_t> > YieldUncertainty_WWControlRegion_Data, YieldUncertainty_WWControlRegion_Bkg, YieldUncertainty_WWControlRegion_WWMC; 
  vector<vector<Double_t> > Yield_WWControlRegion_Top, Yield_WWControlRegion_DY, Yield_WWControlRegion_WJets, Yield_WWControlRegion_Other; 
  vector<vector<Double_t> > YieldUncertainty_WWControlRegion_Top, YieldUncertainty_WWControlRegion_DY, YieldUncertainty_WWControlRegion_WJets, YieldUncertainty_WWControlRegion_Other; 
  vector<vector<Double_t> > WWScaleFactor, WWScaleFactorUncertainty; 


  //4 classes of selections
  for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {

    vector<Double_t>  tmpYield_WWControlRegion_Data, tmpYield_WWControlRegion_Bkg, tmpYield_WWControlRegion_WWMC; 
    vector<Double_t>  tmpYieldUncertainty_WWControlRegion_Data, tmpYieldUncertainty_WWControlRegion_Bkg, tmpYieldUncertainty_WWControlRegion_WWMC; 
    vector<Double_t>  tmpYield_WWControlRegion_Top, tmpYield_WWControlRegion_DY, tmpYield_WWControlRegion_WJets, tmpYield_WWControlRegion_Other; 
    vector<Double_t>  tmpYieldUncertainty_WWControlRegion_Top, tmpYieldUncertainty_WWControlRegion_DY, tmpYieldUncertainty_WWControlRegion_WJets, tmpYieldUncertainty_WWControlRegion_Other; 
    vector<Double_t>  tmpWWScaleFactor, tmpWWScaleFactorUncertainty; 

    for(Int_t imass=0; imass<nmass; imass++) {
      tmpYield_WWControlRegion_Data.push_back(0), tmpYield_WWControlRegion_Bkg.push_back(0), tmpYield_WWControlRegion_WWMC.push_back(0);
      tmpYieldUncertainty_WWControlRegion_Data.push_back(0), tmpYieldUncertainty_WWControlRegion_Bkg.push_back(0), tmpYieldUncertainty_WWControlRegion_WWMC.push_back(0); 
      tmpYield_WWControlRegion_Top.push_back(0), tmpYield_WWControlRegion_DY.push_back(0), tmpYield_WWControlRegion_WJets.push_back(0), tmpYield_WWControlRegion_Other.push_back(0);  
      tmpYieldUncertainty_WWControlRegion_Top.push_back(0), tmpYieldUncertainty_WWControlRegion_DY.push_back(0), tmpYieldUncertainty_WWControlRegion_WJets.push_back(0), tmpYieldUncertainty_WWControlRegion_Other.push_back(0); 
      tmpWWScaleFactor.push_back(0),tmpWWScaleFactorUncertainty.push_back(0) ;
    }

    Yield_WWControlRegion_Data.push_back(tmpYield_WWControlRegion_Data);
    Yield_WWControlRegion_Bkg.push_back(tmpYield_WWControlRegion_Bkg);
    Yield_WWControlRegion_WWMC.push_back(tmpYield_WWControlRegion_WWMC);
    YieldUncertainty_WWControlRegion_Data.push_back(tmpYieldUncertainty_WWControlRegion_Data);
    YieldUncertainty_WWControlRegion_Bkg.push_back(tmpYieldUncertainty_WWControlRegion_Bkg);
    YieldUncertainty_WWControlRegion_WWMC.push_back(tmpYieldUncertainty_WWControlRegion_WWMC);
    Yield_WWControlRegion_Top.push_back(tmpYield_WWControlRegion_Top);
    Yield_WWControlRegion_DY.push_back(tmpYield_WWControlRegion_DY);
    Yield_WWControlRegion_WJets.push_back(tmpYield_WWControlRegion_WJets);
    Yield_WWControlRegion_Other.push_back(tmpYield_WWControlRegion_Other);
    YieldUncertainty_WWControlRegion_Top.push_back(tmpYieldUncertainty_WWControlRegion_Top);
    YieldUncertainty_WWControlRegion_DY.push_back(tmpYieldUncertainty_WWControlRegion_DY);
    YieldUncertainty_WWControlRegion_WJets.push_back(tmpYieldUncertainty_WWControlRegion_WJets);
    YieldUncertainty_WWControlRegion_Other.push_back(tmpYieldUncertainty_WWControlRegion_Other);
    WWScaleFactor.push_back(tmpWWScaleFactor);
    WWScaleFactorUncertainty.push_back(tmpWWScaleFactorUncertainty);
  }



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
  Float_t         met;
  Float_t         trackMet;
  Float_t         pmet;
  Float_t         pTrackMet;
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
  Float_t         higgsPt = -999;

  //***********************************************************************************************
  //Data
  //***********************************************************************************************

  data->SetBranchAddress( "cuts"         , &cuts         );
  data->SetBranchAddress( "dstype"       , &dstype       );
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
  data->SetBranchAddress( "met"          , &met         );
  data->SetBranchAddress( "trackMet"     , &trackMet    );
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
  data->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  data->SetBranchAddress( "jet2Btag"     , &jet2Btag     );
  data->SetBranchAddress( "higgsPt"      , &higgsPt      );

  for (UInt_t i=0; i<data->GetEntries(); i++) {

    if (i%100000 == 0 )
      printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());
    data->GetEntry(i);

    if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
          && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
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
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1



    for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {
      for(Int_t imass=0; imass<nmass; imass++) {

        Bool_t passMassCut = kFALSE;        
        if (classIndex == kMVAZeroJet || classIndex == kMVAOneJet) {
          passMassCut = ( dilep->mass() > DileptonMassPreselectionCut(mH[imass]));
        } else {
          passMassCut = ( dilep->mass() > 100.0 );
        }
        Bool_t passPtMaxCut = (lep1->pt() > cutPtMaxLow(mH[imass]));
        Bool_t passPtMinCut = (lep1->pt() > cutPtMinLow(mH[imass]));
        Bool_t passJetBinCut = kFALSE;
        if (classIndex == kCutBasedZeroJet || classIndex == kMVAZeroJet) {
          passJetBinCut = (njets == 0);
        } else {
          passJetBinCut = (njets == 1);
        }

        if (!(passMassCut && passPtMaxCut && passPtMinCut && passJetBinCut)) continue;

        //WW Sideband Selection
        if ( !( (cuts & SmurfTree::TopTag) == SmurfTree::TopTag )
          ) {
          Yield_WWControlRegion_Data[classIndex][imass]++;
          YieldUncertainty_WWControlRegion_Data[classIndex][imass]++;
        }

      }
    }

  }

  //***********************************************************************************************
  //Background
  //***********************************************************************************************
  background->SetBranchAddress( "cuts"         , &cuts         );
  background->SetBranchAddress( "dstype"       , &dstype       );
  background->SetBranchAddress( "nvtx"         , &nvtx         );
  background->SetBranchAddress( "npu"          , &npu          );
  background->SetBranchAddress( "njets"        , &njets        );
  background->SetBranchAddress( "run"          , &run          );
  background->SetBranchAddress( "event"        , &event        );
  background->SetBranchAddress( "scale1fb"     , &scale1fb     );
  background->SetBranchAddress( "lep1"         , &lep1         );
  background->SetBranchAddress( "lep2"         , &lep2         );
  background->SetBranchAddress( "jet1"         , &jet1         );
  background->SetBranchAddress( "jet2"         , &jet2         );
  background->SetBranchAddress( "jet3"         , &jet3         );
  background->SetBranchAddress( "dPhi"         , &dPhi         );
  background->SetBranchAddress( "dR"           , &dR           );
  background->SetBranchAddress( "dilep"        , &dilep        );
  background->SetBranchAddress( "type"         , &type         );
  background->SetBranchAddress( "met"          , &met         );
  background->SetBranchAddress( "trackMet"     , &trackMet    );
  background->SetBranchAddress( "pmet"         , &pmet         );
  background->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  background->SetBranchAddress( "met"          , &met          );
  background->SetBranchAddress( "mt"           , &mt           );
  background->SetBranchAddress( "mt1"          , &mt1          );
  background->SetBranchAddress( "mt2"          , &mt2          );
  background->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  background->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  background->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  background->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  background->SetBranchAddress( "lq1"          , &lq1          );
  background->SetBranchAddress( "lq2"          , &lq2          );
  background->SetBranchAddress( "lid1"         , &lid1         );
  background->SetBranchAddress( "lid2"         , &lid2         );
  background->SetBranchAddress( "lid3"         , &lid3         );
  background->SetBranchAddress( "processId"    , &processId    );
  background->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  background->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  background->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  background->SetBranchAddress( "jet2Btag"     , &jet2Btag     );
  background->SetBranchAddress( "higgsPt"      , &higgsPt      );
  background->SetBranchAddress( "lep1McId"     , &lep1McId     );
  background->SetBranchAddress( "lep2McId"     , &lep2McId     );

  for (UInt_t i=0; i<background->GetEntries(); i++) {

    if (i%100000 == 0 )
      printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());
    background->GetEntry(i);


    bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
      ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
      dstype != SmurfTree::data;
    if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
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
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    if( (cuts & SmurfTree::TopTag) == SmurfTree::TopTag                  ) continue; // cut on btagging
    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1


    int BkgType = 0;
    if(dstype == SmurfTree::qqww                 ) BkgType = 0;
    else if(dstype == SmurfTree::ggww            ) BkgType = 0;
    else if(dstype == SmurfTree::wz              ) BkgType = 1;
    else if(dstype == SmurfTree::zz              ) BkgType = 1;
    else if(dstype == SmurfTree::tw              ) BkgType = 2;
    else if(dstype == SmurfTree::ttbar           ) BkgType = 2;
    else if(dstype == SmurfTree::dyee            ) BkgType = 3;
    else if(dstype == SmurfTree::dymm            ) BkgType = 3;
    else if(dstype == SmurfTree::wjets           ) BkgType = 4;
    else if(dstype == SmurfTree::data            ) BkgType = 4;
    else if(dstype == SmurfTree::wgamma          ) BkgType = 5;
    else if(dstype == SmurfTree::dytt            ) BkgType = 5;
    else if(dstype == SmurfTree::dyttDataDriven  ) BkgType = 5;
    else if(dstype == SmurfTree::qcd             ) BkgType = 5;


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
        BkgType = 4;
        myWeight = add;
      }
      else if(TMath::Abs(lep1McId)*TMath::Abs(lep2McId) > 0 
              || dstype == SmurfTree::wgamma
              || (0==0)
        ){
    	add = 1.0;
    	add = add*fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    										      (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    										      (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

        add = add*nPUScaleFactor(fhDPUS4,npu);
    	add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
    	add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);

    	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
        BkgType = 4;
        myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
      }
      else {
        myWeight = 0.0;
      }
    }
    else if(dstype == SmurfTree::data) myWeight = 0.0;
    else if(dstype== SmurfTree::dyttDataDriven) {
      //do not apply any MC->data scale factors for data driven dy->tautau 
      myWeight = scale1fb*scaleFactorLum;
    }
    else if(dstype == SmurfTree::qcd) { //Z->tautau embedded sample
      myWeight = 0.0189712*scaleFactorLum;
    }
    else if(dstype != SmurfTree::data){

      //normal MC
      //require both leptons pass full selection
      if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
            && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
 

      add = 1.0;
      add = add*nPUScaleFactor(fhDPUS4,npu);
      add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
      add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
      add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							fabs(lep2->eta()), lep2->pt(), 
    							TMath::Abs( lid1), TMath::Abs(lid2));

      //DY Scale Factors
      if(BkgType == 3  && (type   == SmurfTree::mm   || type   == SmurfTree::ee)
                      && (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm)) {
    	if(njets == 0) add=add*DYBkgScaleFactor(0,0); 
    	if(njets == 1) add=add*DYBkgScaleFactor(0,1); 
    	if(njets >= 2) add=add*DYBkgScaleFactor(0,2); 
      }
      
      //Top Scale Factors
      if(BkgType == 2) {
    	if(njets == 0) add=add*TopBkgScaleFactor(0);
    	if(njets == 1) add=add*TopBkgScaleFactor(1); 
    	if(njets >= 2) add=add*TopBkgScaleFactor(2); 
      }

      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;
      myWeight = scale1fb*scaleFactorLum*add;      
    }

    if(myWeight == 0) continue;


    for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {
      for(Int_t imass=0; imass<nmass; imass++) {

        Bool_t passMassCut = kFALSE;        
        if (classIndex == kMVAZeroJet || classIndex == kMVAOneJet) {
          passMassCut = ( dilep->mass() > DileptonMassPreselectionCut(mH[imass]));
        } else {
          passMassCut = ( dilep->mass() > 100.0 );
        }
        Bool_t passPtMaxCut = (lep1->pt() > cutPtMaxLow(mH[imass]));
        Bool_t passPtMinCut = (lep1->pt() > cutPtMinLow(mH[imass]));
        Bool_t passJetBinCut = kFALSE;
        if (classIndex == kCutBasedZeroJet || classIndex == kMVAZeroJet) {
          passJetBinCut = (njets == 0);
        } else {
          passJetBinCut = (njets == 1);
        }

        if (!(passMassCut && passPtMaxCut && passPtMinCut && passJetBinCut)) continue;

  
        //WW Monte Carlo Yields
        if (BkgType == 0 && myWeight >= 0) {          
          Yield_WWControlRegion_WWMC[classIndex][imass] += myWeight;
          if (dstype == SmurfTree::qqww && classIndex == kCutBasedZeroJet && imass == 1) {
          }
          YieldUncertainty_WWControlRegion_WWMC[classIndex][imass] += myWeight*myWeight;
        } else {
          
          //Backgrounds
          
          if (BkgType == 2) {
            Yield_WWControlRegion_Top[classIndex][imass] += myWeight;
            YieldUncertainty_WWControlRegion_Top[classIndex][imass] += myWeight*myWeight;
          }
          else if (BkgType == 3) {
            Yield_WWControlRegion_DY[classIndex][imass] += myWeight;
            YieldUncertainty_WWControlRegion_DY[classIndex][imass] += myWeight*myWeight;
          }
          else if (BkgType == 4) { 
            Yield_WWControlRegion_WJets[classIndex][imass] += myWeight;
            YieldUncertainty_WWControlRegion_WJets[classIndex][imass] += myWeight*myWeight;
          }
          else if (BkgType != 0) {
            Yield_WWControlRegion_Other[classIndex][imass] += myWeight;
            YieldUncertainty_WWControlRegion_Other[classIndex][imass] += myWeight*myWeight;
          }

        }

      }
    }



  }

  //********************************************************
  // Correct for Systematic Uncertainties
  //********************************************************
  for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {
    for(Int_t imass=0; imass<nmass; imass++) {
      YieldUncertainty_WWControlRegion_WJets[classIndex][imass] += pow( (0.36 * Yield_WWControlRegion_WJets[classIndex][imass]) , 2);
      YieldUncertainty_WWControlRegion_DY[classIndex][imass] += pow( (1.00 * Yield_WWControlRegion_DY[classIndex][imass]) , 2);
      YieldUncertainty_WWControlRegion_Top[classIndex][imass] += pow( (0.20 * Yield_WWControlRegion_Top[classIndex][imass]) , 2);
      YieldUncertainty_WWControlRegion_Other[classIndex][imass] += pow( (0.20 * Yield_WWControlRegion_Other[classIndex][imass]) , 2);

    }
  }

  //Do Square Roots
  for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {
    for(Int_t imass=0; imass<nmass; imass++) {

      YieldUncertainty_WWControlRegion_Data[classIndex][imass] = TMath::Sqrt(YieldUncertainty_WWControlRegion_Data[classIndex][imass]);
      YieldUncertainty_WWControlRegion_WJets[classIndex][imass] = TMath::Sqrt(YieldUncertainty_WWControlRegion_WJets[classIndex][imass]);
      YieldUncertainty_WWControlRegion_DY[classIndex][imass] = TMath::Sqrt(YieldUncertainty_WWControlRegion_DY[classIndex][imass]);
      YieldUncertainty_WWControlRegion_Top[classIndex][imass] = TMath::Sqrt(YieldUncertainty_WWControlRegion_Top[classIndex][imass]);
      YieldUncertainty_WWControlRegion_Other[classIndex][imass] = TMath::Sqrt(YieldUncertainty_WWControlRegion_Other[classIndex][imass]);

    }
  }

  //********************************************************
  // Add All Backgrounds Together
  //********************************************************
  //Do Square Roots
  for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {
    for(Int_t imass=0; imass<nmass; imass++) {
      Yield_WWControlRegion_Bkg[classIndex][imass] = Yield_WWControlRegion_DY[classIndex][imass] +
        Yield_WWControlRegion_WJets[classIndex][imass] + 
        Yield_WWControlRegion_Top[classIndex][imass] + 
        Yield_WWControlRegion_Other[classIndex][imass];

      YieldUncertainty_WWControlRegion_Bkg[classIndex][imass] = TMath::Sqrt( pow(YieldUncertainty_WWControlRegion_DY[classIndex][imass], 2) +
                                                                             pow(YieldUncertainty_WWControlRegion_WJets[classIndex][imass], 2) +
                                                                             pow(YieldUncertainty_WWControlRegion_Top[classIndex][imass], 2) +
                                                                             pow(YieldUncertainty_WWControlRegion_Other[classIndex][imass], 2)
        );
    }
  }


  //********************************************************
  // Print Summary
  //********************************************************
  for(UInt_t classIndex = 0; classIndex < 4; ++classIndex) {
    for(Int_t imass=0; imass<nmass; imass++) {

      double SF = (Yield_WWControlRegion_Data[classIndex][imass] - Yield_WWControlRegion_Bkg[classIndex][imass]) / Yield_WWControlRegion_WWMC[classIndex][imass];
      double SFUncertainty = SF*TMath::Sqrt( 
        ( Yield_WWControlRegion_Data[classIndex][imass] + pow(YieldUncertainty_WWControlRegion_Bkg[classIndex][imass],2)) / 
        pow(Yield_WWControlRegion_Data[classIndex][imass] - Yield_WWControlRegion_Bkg[classIndex][imass] ,2) 
        +
        pow( YieldUncertainty_WWControlRegion_Bkg[classIndex][imass] / Yield_WWControlRegion_WWMC[classIndex][imass] , 2)
        );
      WWScaleFactor[classIndex][imass] = SF;
      WWScaleFactorUncertainty[classIndex][imass] = SFUncertainty;

      //if (!(classIndex == kCutBasedZeroJet && imass == 1)) continue;

      cout << "mH = 120 Analysis : WW Control Region" << endl;
      cout << "Data Yield : " << Yield_WWControlRegion_Data[classIndex][imass] 
           << " +/- " << YieldUncertainty_WWControlRegion_Data[classIndex][imass] << endl;
      cout << "DY Bkg : " 
           << Yield_WWControlRegion_DY[classIndex][imass] << " +/- " 
           << YieldUncertainty_WWControlRegion_DY[classIndex][imass] << endl;
      cout << "WJets Bkg : " 
           << Yield_WWControlRegion_WJets[classIndex][imass] << " +/- " 
           << YieldUncertainty_WWControlRegion_WJets[classIndex][imass] << endl;
      cout << "Top Bkg : " 
           << Yield_WWControlRegion_Top[classIndex][imass] << " +/- " 
           << YieldUncertainty_WWControlRegion_Top[classIndex][imass] << endl;
      cout << "ZZ,WZ,DYtautau Bkg : " 
           << Yield_WWControlRegion_Other[classIndex][imass] << " +/- " 
           << YieldUncertainty_WWControlRegion_Other[classIndex][imass] << endl;
      cout << "Total Bkg : " 
           << Yield_WWControlRegion_Bkg[classIndex][imass] << " +/- " 
           << YieldUncertainty_WWControlRegion_Bkg[classIndex][imass] << endl;
      cout << "******************************************************************\n";
      cout << "WW Yield from MC : " 
           << Yield_WWControlRegion_WWMC[classIndex][imass] << " +/- " 
           << YieldUncertainty_WWControlRegion_WWMC[classIndex][imass] << endl;
      cout << "ScaleFactor : " 
           << WWScaleFactor[classIndex][imass] << " +/- " 
           << WWScaleFactorUncertainty[classIndex][imass] << endl;
    }
  }
  

  //********************************************************
  // Output WW Scale Factor
  //********************************************************
  ofstream outf("WWBkgScaleFactors.h");
  outf << "Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {" << endl;
  outf << "assert(jetBin >= 0 && jetBin <= 1);" << endl;

  outf << "  Int_t mHiggs[" << nmass << "] = {";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << mH[i];
    if (i < nmass-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t WWBkgScaleFactorHiggsSelection[2][" << nmass << "] = { " << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << WWScaleFactor[kCutBasedZeroJet][i];    
    if (i < nmass-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << WWScaleFactor[kCutBasedOneJet][i];    
    if (i < nmass-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return 1.0;" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;



  outf << "Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {" << endl;
  outf << "assert(jetBin >= 0 && jetBin <= 1);" << endl;

  outf << "  Int_t mHiggs[" << nmass << "] = {";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << mH[i];
    if (i < nmass-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t WWBkgScaleFactorHiggsSelection[2][" << nmass << "] = { " << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << WWScaleFactor[kMVAZeroJet][i];    
    if (i < nmass-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << WWScaleFactor[kMVAOneJet][i];    
    if (i < nmass-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return 1.0;" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;


  outf << "Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {" << endl;
  outf << "assert(jetBin >= 0 && jetBin <= 1);" << endl;

  outf << "  Int_t mHiggs[" << nmass << "] = {";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << mH[i];
    if (i < nmass-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t WWBkgScaleFactorKappaHiggsSelection[2][" << nmass << "] = { " << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << ( 1.0 + WWScaleFactorUncertainty[kCutBasedZeroJet][i] / WWScaleFactor[kCutBasedZeroJet][i]);    
    if (i < nmass) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << ( 1.0 + WWScaleFactorUncertainty[kCutBasedOneJet][i] / WWScaleFactor[kCutBasedOneJet][i]);    
    if (i < nmass-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return 1.0;" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;



  outf << "Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {" << endl;
  outf << "assert(jetBin >= 0 && jetBin <= 1);" << endl;

  outf << "  Int_t mHiggs[" << nmass << "] = {";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << mH[i];
    if (i < nmass-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t WWBkgScaleFactorKappaHiggsSelection[2][" << nmass << "] = { " << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << ( 1.0 + WWScaleFactorUncertainty[kMVAZeroJet][i] / WWScaleFactor[kMVAZeroJet][i]);    
    if (i < nmass-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (Int_t i = 0; i < nmass ; ++i) {
    outf << ( 1.0 + WWScaleFactorUncertainty[kMVAOneJet][i] / WWScaleFactor[kMVAOneJet][i]);    
    if (i < nmass-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return 1.0;" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;

  outf.close();



  return;

}


