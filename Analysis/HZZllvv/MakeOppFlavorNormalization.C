//root -l Smurf/Analysis/HZZllvv/MakeOppFlavorNormalization.C+'(0,300,"","/data/smurf/sixie/data/Run2011_Summer11_EPSHZZV0/mitf-alljets/background.root","/data/smurf/sixie/data/Run2011_Summer11_EPSHZZV0/mitf-alljets/data_2l.goodlumi1092ipb.root")'
//root -l Smurf/Analysis/HZZllvv/MakeOppFlavorNormalization.C+'(0,300,"","/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/background.root","/data/smurf/data/EPS/tas/data-met20-1092ipb.root")'
//root -l Smurf/Analysis/HZZllvv/MakeOppFlavorNormalization.C+'(0,300,"","/data/smurf/sixie/data/Run2011_Summer11_SmurfV6/mitf-alljets/background.root","/data/smurf/sixie/data/Run2011_Summer11_SmurfV6/mitf-alljets/data_2l.goodlumi1092ipb.root")'


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
#include "Smurf/Analysis/HZZllvv/factors.h"
#include "Smurf/Analysis/HZZllvv/HiggsQCDScaleSystematics.h"
#include "Smurf/Analysis/HZZllvv/PSUESystematics.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HZZllvv/PileupReweighting.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void MakeOppFlavorNormalization
(
 UInt_t  nJetsType   	 = 0,
 UInt_t  mH      	 = 250,
 TString outTag          = "default",
 TString bgdInputFile    = "data/ntuples_250train_hzzjets_zz.root",
 TString datInputFile    = "data/ntuples_250train_hzzjets_data_2l.root"
 )
{

  //************************************************************************************************
  // Setup
  //************************************************************************************************
 
  TString bgdFile1    = bgdInputFile;
  TString datFile1    = datInputFile;

  //************************************************************************************************
  // Load Trees
  //************************************************************************************************

   TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
  chdata->Add(datFile1);

  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;



  //************************************************************************************************
  // PU Reweighting
  //************************************************************************************************
  TFile *fPUWeightsFile = TFile::Open("Smurf/Analysis/HZZllvv/PUWeights.root");
  TH1D *fPUWeightsS3 = (TH1D*)(fPUWeightsFile->Get("PUWeights_BX0_S3"));
  TH1D *fPUWeightsS4 = (TH1D*)(fPUWeightsFile->Get("PUWeights_BX0_S4"));
  fPUWeightsS3->SetDirectory(0);
  fPUWeightsS4->SetDirectory(0);
  fPUWeightsFile->Close();
  delete fPUWeightsFile;
  TH1D *myPUWeights = 0;

  //************************************************************************************************
  // Lepton Efficiencies
  //************************************************************************************************

  LeptonScaleLookup trigLookup("Smurf/Analysis/HZZllvv/efficiency_42X.root");

  TFile *fLeptonFRFileM = TFile::Open("/data/smurf/data/EPS/auxiliar/FakeRates_SmurfV6.root");
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open("/data/smurf/data/EPS/auxiliar/FakeRates_SmurfV6.root");
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

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
  // Luminosity
  //************************************************************************************************

  double scaleFactorLum = 1.092;


  //************************************************************************************************
  // higgs masses
  //************************************************************************************************
  Int_t mHIndex = -1;
  if (mH == 200) mHIndex = 0;
  if (mH == 250) mHIndex = 1;
  if (mH == 275) mHIndex = 2;
  if (mH == 300) mHIndex = 3;
  if (mH == 325) mHIndex = 4;
  if (mH == 350) mHIndex = 5;
  if (mH == 375) mHIndex = 6;
  if (mH == 400) mHIndex = 7;
  if (mH == 425) mHIndex = 8;
  if (mH == 450) mHIndex = 9;
  if (mH == 475) mHIndex = 10;
  if (mH == 500) mHIndex = 11;
  if (mH == 525) mHIndex = 12;
  if (mH == 550) mHIndex = 13;
  if (mH == 575) mHIndex = 14;
  if (mH == 600) mHIndex = 15;
  if (mHIndex == -1 ) return;


  //************************************************************************************************
  // mH dependant Cuts
  // [NJETS][mH]
  // double NJETS[3] = { 0, 1, 2 };
  // double mH[16] = { 200, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600 };
  //************************************************************************************************
  double cutMETLow[3][16]        = { { 50, 60, 60, 70 , 70, 70, 70, 70, 50, 60, 80, 90,110,120,130,150},
                                     { 60, 60, 60, 100,100,100,100,100, 50, 60, 80, 90,110,120,130,150},
                                     { 80, 80, 80, 100,100,100,100,100, 50, 60, 80, 90,110,120,130,150} };
  double cutMTLow[3][16]         = { { 180,220,220,260,320,320,320,320,110,120,120,120,120,120,120,120},
                                     { 180,180,180,260,300,300,300,300,110,120,120,120,120,120,120,120},
                                     { 180,180,180,240,300,300,300,300,110,120,120,120,120,120,120,120} };
  double cutMTHigh[3][16]        = { { 220,260,260,320,450,450,450,450,450,180,190,200,210,220,230,250}, 
                                     { 200,260,260,320,450,450,450,450,450,180,190,200,210,220,230,250}, 
                                     { 200,260,260,320,450,450,450,450,450,180,190,200,210,220,230,250} };
  double cutDPhiMetJet[3][16]    = { { 0.62,0.62,0.62,0.28,0.28,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14},
                                     { 0.62,0.62,0.62,0.28,0.28,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14},
                                     { 0.62,0.62,0.62,0.28,0.28,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14} };


  //----------------------------------------------------------------------------

  //************************************************************************************************
  //Yields
  //************************************************************************************************
  const int nBkgChannels = 8; //number of different background channels
  TH1D* bkgMet[nBkgChannels];
  for(int j=0; j<nBkgChannels; j++) {
    bkgMet[j] = new TH1D(Form("bkgMet_%d",j), Form("bkgMet_%d",j), 100, 0, 200);
    bkgMet[j]->Sumw2();
  }
  TH1D* dataMet = new TH1D( "dataMet", "bkgMet", 100, 0, 200);


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
  UInt_t          npu;


  //**************************************************************************
  //Background Loop
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
  background->SetBranchAddress( "npu"           , &npu );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    //************************************************************************************************
    // PreSelection
    //************************************************************************************************
    bool passTightPlusDenominatorSelection = ( dstype == SmurfTree::data && 
                                               (((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
                                                 (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) 
                                                ||
                                                ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && 
                                                 (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
                                                 ));
    bool passTightPlusTightSelection = ( ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) &&
                                         ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) );
    
    
    if( njets != nJetsType ) continue; // select n-jet type events
    if( !(lep1->pt() > 20 && lep2->pt() > 20 )) continue;
    if( !( (dstype != SmurfTree::data && passTightPlusTightSelection) || passTightPlusDenominatorSelection)) continue;
    if( !(type == SmurfTree::em || type == SmurfTree::me )) continue;
    if( !(dilep->mass() > 12 )) continue;
    if( !(lq1*lq2 < 0)) continue;
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( njets < 2 && jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 

    //************************************************************************************************
    // Background Process
    //************************************************************************************************
    Int_t bkgIndex = -1;
    if (dstype == SmurfTree::zz || dstype == SmurfTree::ggzz) bkgIndex = 0;
    else if (dstype == SmurfTree::wz) bkgIndex = 1;
    else if (dstype == SmurfTree::qqww || dstype == SmurfTree::ggww ) bkgIndex = 2;
    else if (dstype == SmurfTree::ttbar || dstype == SmurfTree::tw ) bkgIndex = 3;
    else if (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm ) bkgIndex = 4;
    else if (dstype == SmurfTree::dytt ) bkgIndex = 5;
    else if (dstype == SmurfTree::wjets || dstype == SmurfTree::data) bkgIndex = 6;
    else bkgIndex = 7;

    if (bkgIndex == 7) { cout << "other bkg: " << dstype << " " << endl; }

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
        myWeight = weightFactor;

        cout << "data: " << dstype << " " << weightFactor 
             << " : " << ((cuts & SmurfTree::Lep1LooseMuV1) == SmurfTree::Lep1LooseMuV1) << " " << ((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) << " , " << ((cuts & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep2LooseMuV1) << " " << ((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) << " "
             << " : " << fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, ((cuts & SmurfTree::Lep1LooseMuV1) == SmurfTree::Lep1LooseMuV1), ((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4)) << " " << fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, ((cuts & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep2LooseMuV1), ((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4)) 
             << endl;


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

      //Some MC's may be S4, some are S3
      if (dstype == SmurfTree::zz || dstype == SmurfTree::ggzz ||
          dstype == SmurfTree::wz || dstype == SmurfTree::qqww || 
          dstype == SmurfTree::ggww ||dstype == SmurfTree::ttbar 
          || dstype == SmurfTree::tw || dstype == SmurfTree::dyee 
          || dstype == SmurfTree::dymm || dstype == SmurfTree::dytt) {
        myPUWeights = fPUWeightsS4;
      } else {
        myPUWeights = fPUWeightsS3;
      }
      weightFactor *= myPUWeights->GetBinContent(myPUWeights->GetXaxis()->FindFixBin(npu));
//       weightFactor *= PileupReweightFactor(nvtx);
      
      
      if((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
        weightFactor  *= trigLookup.GetExpectedLeptonSF( lep1->eta(), lep1->pt(), lid1);
      if((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
        weightFactor  *= trigLookup.GetExpectedLeptonSF( lep2->eta(), lep2->pt(), lid2);
      double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
                                                               fabs(lep2->eta()), lep2->pt(), 
                                                               TMath::Abs(lid1), TMath::Abs(lid2));
      weightFactor *= trigEff;
      
      //ZZ DileptonPt dependant KFactor
      if (dstype == SmurfTree::zz) {
        weightFactor *= ZZKFactor->GetBinContent( ZZKFactor->GetXaxis()->FindFixBin(dilep->pt()));
      }
      
      myWeight 	     = scale1fb*scaleFactorLum*weightFactor;
    }
    
    //************************************************************************************************
    // Yields and histograms
    //************************************************************************************************
    assert(bkgIndex >= 0 && bkgIndex < nBkgChannels);
    bkgMet[bkgIndex]->Fill(TMath::Min(met,trackMet), myWeight);


  }
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

    //************************************************************************************************
    // PreSelection
    //************************************************************************************************

    if( njets != nJetsType ) continue; // select n-jet type events
    if( !(lep1->pt() > 20 && lep2->pt() > 20 )) continue;
    if( !((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection )) continue;
    if( !((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection )) continue;
    if( !(type == SmurfTree::em || type == SmurfTree::me )) continue;
    if( !(dilep->mass() > 12 )) continue;
    if( !(lq1*lq2 < 0)) continue;
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( njets < 2 && jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 

    double myWeight = 1.0;

    //************************************************************************************************
    // Yields and histograms
    //************************************************************************************************
    dataMet->Fill(TMath::Min(met,trackMet));


  }

  //**************************************************************************
  //Compute Yields for different Met cuts
  //**************************************************************************
  const Int_t nMetBins = 5;
  Double_t metThresholds[nMetBins] = { 20.0, 30.0, 40.0, 50.0, 60.0 };
  Double_t bkgYields[nBkgChannels+1][nMetBins];
  Double_t bkgYieldsErrorSqr[nBkgChannels+1][nMetBins];
  for (UInt_t i=0; i<nBkgChannels+1; ++i) {
    for (UInt_t j=0; j<nMetBins; ++j) {
      bkgYields[i][j] = 0.0;
      bkgYieldsErrorSqr[i][j] = 0.0;
    }
  }
  Double_t dataYields[nMetBins] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  //loop over different met cuts
  for (UInt_t j=0; j<nMetBins; ++j) {

    //bkg yields
    for (UInt_t bkgIndex=0; bkgIndex<nBkgChannels; ++bkgIndex) {
      Double_t tmpYield = 0;
      Double_t tmpYieldErrorSqr = 0;
      for (UInt_t bin=0; bin < bkgMet[bkgIndex]->GetXaxis()->GetNbins()+2; ++bin) {
        if (bkgMet[bkgIndex]->GetXaxis()->GetBinUpEdge(bin) >  metThresholds[j]) {
          tmpYield += bkgMet[bkgIndex]->GetBinContent(bin);
          tmpYieldErrorSqr += pow(bkgMet[bkgIndex]->GetBinError(bin),2);
        }
      }
      bkgYields[bkgIndex][j] = tmpYield;
      bkgYieldsErrorSqr[bkgIndex][j] = tmpYieldErrorSqr;
    }
    
    //add up all bkgs
    double tmpTotal = 0.0;
    double tmpTotalErrorSqr = 0.0;
    for (UInt_t bkgIndex=0; bkgIndex<nBkgChannels; ++bkgIndex) {
      tmpTotal += bkgYields[bkgIndex][j];
      tmpTotalErrorSqr += bkgYieldsErrorSqr[bkgIndex][j];
    }
    bkgYields[nBkgChannels][j] = tmpTotal;
    bkgYieldsErrorSqr[nBkgChannels][j] = tmpTotalErrorSqr;

    //data yield
    double tmpDataTotal = 0.0;
    for (UInt_t bin=0; bin < dataMet->GetXaxis()->GetNbins()+2; ++bin) {
      if (dataMet->GetXaxis()->GetBinUpEdge(bin) >  metThresholds[j]) {
        tmpDataTotal += dataMet->GetBinContent(bin);
      }
    }
    dataYields[j] = tmpDataTotal;
  }
  

  //**************************************************************************
  //Summary
  //**************************************************************************

  printf("| *%d jet* | *ZZ* | *WZ* | *WW* | *top* | *Ztt* | *wjets* | *TotalBkg* | *Data* | *Data/MC Ratio* |\n", nJetsType);

  for(UInt_t i=0; i<nMetBins; ++i) {
    Double_t numerator = dataYields[i] - ( bkgYields[0][i] + bkgYields[1][i] + bkgYields[6][i] );
    Double_t denomiantor = bkgYields[2][i] + bkgYields[3][i] + bkgYields[5][i];
    double nrelerr = sqrt(dataYields[i] + bkgYieldsErrorSqr[0][i] + bkgYieldsErrorSqr[1][i] + bkgYieldsErrorSqr[6][i])/numerator;
    double drelerr = sqrt(bkgYieldsErrorSqr[2][i] + bkgYieldsErrorSqr[3][i] + bkgYieldsErrorSqr[5][i])/denomiantor;
    double ratio = numerator/denomiantor;
    double ratioerr = ratio*sqrt(nrelerr*nrelerr + drelerr*drelerr);


    printf("| *Met > %8.3f* | %8.3f +/- %8.3f | %8.3f +/- %8.3f | %8.3f +/- %8.3f | %8.3f +/- %8.3f | %8.3f +/- %8.3f | %8.3f +/- %8.3f | %8.3f  +/- %8.3f| %8.3f +/- %8.3f | %8.3f +/- %8.3f |\n", metThresholds[i], bkgYields[0][i], TMath::Sqrt(bkgYieldsErrorSqr[0][i]) , bkgYields[1][i], TMath::Sqrt(bkgYieldsErrorSqr[1][i]), bkgYields[2][i], TMath::Sqrt(bkgYieldsErrorSqr[2][i]), bkgYields[3][i], TMath::Sqrt(bkgYieldsErrorSqr[3][i]), bkgYields[5][i], TMath::Sqrt(bkgYieldsErrorSqr[5][i]), bkgYields[6][i], TMath::Sqrt(bkgYieldsErrorSqr[6][i]), bkgYields[nBkgChannels][i], TMath::Sqrt(bkgYieldsErrorSqr[nBkgChannels][i]), dataYields[i], TMath::Sqrt(dataYields[i]) , ratio , ratioerr );
  }



}

