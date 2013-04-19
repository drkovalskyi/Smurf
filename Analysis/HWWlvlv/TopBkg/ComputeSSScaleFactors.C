#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors.h"
#include "Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h"
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
void ComputeSSScaleFactors
(
 Int_t period = 2,
 TString bgdInputFile    = "/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/backgroundA_skim1.root",
 TString dataInputFile   = "/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/data_2l_skim1.root"
 )
{
  //*******************************************************************************
  //Settings 
  //*******************************************************************************
  int category = 0;
  if(period > 9) {period = period - 10; category = 1;}
  double lumi = 1;
  enum { kOther, kTTBAR, kTW, kData };
  unsigned int patternTopVeto = SmurfTree::TopVeto;

  TString effPath  = "/data/smurf/data/LP2011/auxiliar/efficiency_results_v6_42x.root";
  TString fakePath = "/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.LP2011.root";
  TString puPath   = "/data/smurf/data/LP2011/auxiliar/puWeights_PU4_68mb.root";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  if	 (period == 0){ // Run2011A
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Run2011A.root";
    //lumi     = 2.1;minRun =      0;maxRun = 173692;
    lumi     = 1.1;minRun =      0;maxRun = 167913;
  }
  else if(period == 1){ // Run2011B
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    //puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Run2011B.root";
    //lumi     = 1.9;minRun = 173693;maxRun = 999999;
    puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Full2011.root";
    lumi     = 3.6;minRun = 167914;maxRun = 999999;
  }
  else if(period == 2){ // Full2011-Fall11-V7
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Full2011.root";
    lumi     = 4.924;minRun =      0;maxRun = 999999;
  }
  else if(period == 3){ // Full2011-Fall11-V7
    effPath  = "/data/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root";
    lumi     = 4.924;minRun =      0;maxRun = 999999;
  }
  else if(period == 4){ // Full2011-Fall11-V8
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV8_42X/auxiliar/efficiency_results_MVAIDIsoCombinedDetIsoSameSigWP_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV8_42X/auxiliar/FakeRates_MVAIDIsoCombinedDetIsoSameSigWP.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV8_42X/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root";
    lumi     = 4.924;minRun =      0;maxRun = 999999;
  }
  else if(period == 5){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_MVAIDIsoCombinedDetIsoSameSigWP_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_MVAIDIsoCombinedDetIsoSameSigWP.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root";
    lumi     = 4.924;minRun =      0;maxRun = 999999;
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
  vector<double> btag_highestpt_2j_den,btag_highestpt_2j_num;
  vector<double> btag_highestpt_2j_den_error,btag_highestpt_2j_num_error;
  vector<double> qcharge_zsel_den,qcharge_zsel_num;
  vector<double> qcharge_zsel_den_error,qcharge_zsel_num_error;
  vector<double> vhselNjet2_den,vhselNjet2_num;
  vector<double> vhselNjet2_den_error,vhselNjet2_num_error;
  vector<double> vhselNjet1_den,vhselNjet1_num;
  vector<double> vhselNjet1_den_error,vhselNjet1_num_error;
  vector<double> evt_topNjet2_B,evt_topNjet2_E;
  vector<double> evt_topNjet1_B,evt_topNjet1_E;

  for(int i=0; i<4; i++){
    btag_highestpt_2j_den.push_back(0),btag_highestpt_2j_num.push_back(0),btag_highestpt_2j_den_error.push_back(0),btag_highestpt_2j_num_error.push_back(0);
    qcharge_zsel_den.push_back(0),qcharge_zsel_num.push_back(0),qcharge_zsel_den_error.push_back(0),qcharge_zsel_num_error.push_back(0);
    vhselNjet2_den.push_back(0),vhselNjet2_num.push_back(0),vhselNjet2_den_error.push_back(0),vhselNjet2_num_error.push_back(0);
    vhselNjet1_den.push_back(0),vhselNjet1_num.push_back(0),vhselNjet1_den_error.push_back(0),vhselNjet1_num_error.push_back(0);
    evt_topNjet2_B.push_back(0),evt_topNjet2_E.push_back(0);
    evt_topNjet1_B.push_back(0),evt_topNjet1_E.push_back(0);
  }

  unsigned int patternTopTagNotInJets = SmurfTree::TopTagNotInJets;

  //*******************************************************************************
  //Background Events
  //*******************************************************************************
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
    if(bgdEvent.lep1_.Pt() <= 20.0) continue;
    if(bgdEvent.lep2_.Pt() <= 10.0) continue;
    if(bgdEvent.dilep_.M() <= 12) continue;

    bool passIdCuts = true;
    if(category == 1){
      passIdCuts = (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection ||
                   (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection;
    }
    if(passIdCuts == false) continue;

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
    else if(bgdEvent.dstype_ == SmurfTree::qqwwPWG         ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqww2j          ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz              ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz              ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww            ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma          ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::data            ) fDecay =  1;
    else {fDecay = 0;cout << bgdEvent.dstype_ << endl;assert(0);}

    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);

    double minmet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool passMETA = minmet > 20. &&
         	   (minmet > 37.+bgdEvent.nvtx_/2.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    bool passMETB = minmet > 20. &&
         	   (minmet > 22.+bgdEvent.nvtx_/2.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);

    // begin computing weights
    double theWeight = 0.0;
    double add       = 1.0;
    int nFake = 0;
    if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
    if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
    if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
    if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
    if(category == 1) {
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV1)  == SmurfTree::Lep1LooseEleV1)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake--;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV1)  == SmurfTree::Lep2LooseEleV1)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake--;
    }
    if(nFake < 0) assert(0);

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
	if(category == 1) add = add*1.10; // HACK!!!!
	theWeight	       = add;
      }
      else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){
        add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	add = add*nPUScaleFactor(fhDPUS4,bgdEvent.npu_);
        add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	        					         fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
							         TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        add = add*trigEff;
	if(category == 1) add = add*1.10; // HACK!!!!
	fDecay  	       = 3;
	theWeight	       = -1.0 * bgdEvent.scale1fb_*lumi*add;
      }
    }
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven || bgdEvent.dstype_ == SmurfTree::qcd) {
      theWeight = ZttScaleFactor(bgdEvent.nvtx_,period,bgdEvent.scale1fb_)*lumi;
    }
    else if(bgdEvent.dstype_ != SmurfTree::data){
      double add1 = nPUScaleFactor(fhDPUS4,bgdEvent.npu_);

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

    if(theWeight == 0.0) continue;

    int classType = kOther;
    if(fDecay == 5 ) classType = kTTBAR;
    if(fDecay == 13) classType = kTW;

    bool isTopSel = kFALSE; bool isZSel = kFALSE;
    bool isVHSSSelNjets1 = kFALSE; bool isVHSSSelNjets2 = kFALSE;
    if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
       bgdEvent.dstype_ != SmurfTree::data) &&
      charge == 0 &&
      passMETA == true &&
      (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
      bgdEvent.njets_ >= 2 && bgdEvent.njets_ <= 3 &&
      TMath::Abs(bgdEvent.jet1_.Eta()) < 2.4 &&
      TMath::Abs(bgdEvent.jet2_.Eta()) < 2.4 &&
      (bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets &&
      1 == 1
      ) isTopSel = kTRUE;

    if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
       bgdEvent.dstype_ != SmurfTree::data) &&
      fabs(bgdEvent.dilep_.M()-91.1876) < 10. && bgdEvent.type_ == SmurfTree::ee && 
      1 == 1
      ) isZSel = kTRUE;

    double deltaPhiQQL[3] = {DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep1_.Phi()),DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep2_.Phi()),0.0};
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
       bgdEvent.dstype_ != SmurfTree::data) &&
      passMETB == true &&
      bgdEvent.type_ != SmurfTree::mm &&
     (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ != SmurfTree::ee) && 
      bgdEvent.njets_ >= 2 && bgdEvent.njets_ <= 3 &&
      TMath::Abs(bgdEvent.jet1_.Eta()) < 2.4 &&
      TMath::Abs(bgdEvent.jet2_.Eta()) < 2.4 &&
      (bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets &&
      bgdEvent.jet2Btag_ < 2.1 && bgdEvent.jet3Btag_ < 2.1 &&
      bgdEvent.dilep_.M() < 200 &&
     (bgdEvent.jet1_+bgdEvent.jet2_).M() > 60. && (bgdEvent.jet1_+bgdEvent.jet2_).M() < 110. &&
      bgdEvent.mt_ > 70. && bgdEvent.mt_ < 200. && 
      deltaPhiQQL[2]*180.0/TMath::Pi() < 110.0 &&
      1 == 1
      ) isVHSSSelNjets2 = kTRUE;

    deltaPhiQQL[0] = DeltaPhi(bgdEvent.jet1_.Phi(),bgdEvent.lep1_.Phi()); deltaPhiQQL[1] = DeltaPhi(bgdEvent.jet1_.Phi(),bgdEvent.lep2_.Phi());
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    if(
      minmet > 37.+bgdEvent.nvtx_/2.0 &&
      bgdEvent.type_ != SmurfTree::mm &&
     (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ != SmurfTree::ee) && 
      bgdEvent.njets_ == 1 &&
      TMath::Abs(bgdEvent.jet1_.Eta()) < 2.4 &&
      (bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets &&
      bgdEvent.dilep_.M() < 200 &&
      bgdEvent.mt_ > 70. && bgdEvent.mt_ < 200. && 
      deltaPhiQQL[2]*180.0/TMath::Pi() < 85.0 &&
      1 == 1
      ) isVHSSSelNjets1 = kTRUE;

     if(isRealLepton == false &&
     	(bgdEvent.dstype_ == SmurfTree::ttbar || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
     	 bgdEvent.dstype_ == SmurfTree::qqww  || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
     	 bgdEvent.dstype_ == SmurfTree::wgstar)) {isZSel = kFALSE; isVHSSSelNjets1 = kFALSE; isVHSSSelNjets2 = kFALSE; isTopSel = kFALSE;}

    if(isTopSel == kTRUE && bgdEvent.jet2Btag_ >= 2.1 && bgdEvent.jet3Btag_ < 2.1){
      btag_highestpt_2j_den[classType]	       += theWeight;
      btag_highestpt_2j_den_error[classType]   += theWeight*theWeight;
      if(bgdEvent.jet1Btag_ >= 2.1){
    	btag_highestpt_2j_num[classType]       += theWeight;
    	btag_highestpt_2j_num_error[classType] += theWeight*theWeight;
      }
    }

    if(isZSel == kTRUE){
      if     (fabs(bgdEvent.lep1_.Eta()) < 1.479 && fabs(bgdEvent.lep2_.Eta()) < 1.479){
        qcharge_zsel_den[0] 	   += 1.0;
        qcharge_zsel_den_error[0]  += 1.0*1.0;
      }
      else if(fabs(bgdEvent.lep1_.Eta()) > 1.479 && fabs(bgdEvent.lep2_.Eta()) > 1.479){
        qcharge_zsel_den[1] 	   += 1.0;
        qcharge_zsel_den_error[1]  += 1.0*1.0;
      }
      if(charge != 0){
	if     (fabs(bgdEvent.lep1_.Eta()) < 1.479 && fabs(bgdEvent.lep2_.Eta()) < 1.479){
          qcharge_zsel_num[0] 	     += 1.0;
          qcharge_zsel_num_error[0]  += 1.0*1.0;
	}
	else if(fabs(bgdEvent.lep1_.Eta()) > 1.479 && fabs(bgdEvent.lep2_.Eta()) > 1.479){
          qcharge_zsel_num[1] 	     += 1.0;
          qcharge_zsel_num_error[1]  += 1.0*1.0;
	}
      }
    }

    if(isVHSSSelNjets2 == kTRUE){
      if     (charge == 0 && bgdEvent.jet1Btag_ >= 2.1) {
        vhselNjet2_den[classType]	   += theWeight;
        vhselNjet2_den_error[classType] += theWeight*theWeight;
	if      (bgdEvent.type_ == SmurfTree::ee){
	  if(fabs(bgdEvent.lep2_.Eta()) < 1.479) evt_topNjet2_B[classType]	   += theWeight;
	  else                                   evt_topNjet2_E[classType]	   += theWeight;
	}
	else if(abs(bgdEvent.lid1_) == 11){
	  if(fabs(bgdEvent.lep1_.Eta()) < 1.479) evt_topNjet2_B[classType]	   += theWeight;
	  else                                   evt_topNjet2_E[classType]	   += theWeight;
	}
	else if(abs(bgdEvent.lid2_) == 11){
	  if(fabs(bgdEvent.lep2_.Eta()) < 1.479) evt_topNjet2_B[classType]	   += theWeight;
	  else                                   evt_topNjet2_E[classType]	   += theWeight;
	}	  
      } 
      else if(charge != 0 && bgdEvent.jet1Btag_  < 2.1 && (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto){
    	vhselNjet2_num[classType]       += theWeight;
    	vhselNjet2_num_error[classType] += theWeight*theWeight;
      }
    }

    if(isVHSSSelNjets1 == kTRUE){
      if     (charge == 0 && bgdEvent.jet1Btag_ >= 2.1) {
        vhselNjet1_den[classType]	   += theWeight;
        vhselNjet1_den_error[classType] += theWeight*theWeight;
	if      (bgdEvent.type_ == SmurfTree::ee){
	  if(fabs(bgdEvent.lep2_.Eta()) < 1.479) evt_topNjet1_B[classType]	   += theWeight;
	  else                                   evt_topNjet1_E[classType]	   += theWeight;
	}
	else if(abs(bgdEvent.lid1_) == 11){
	  if(fabs(bgdEvent.lep1_.Eta()) < 1.479) evt_topNjet1_B[classType]	   += theWeight;
	  else                                   evt_topNjet1_E[classType]	   += theWeight;
	}
	else if(abs(bgdEvent.lid2_) == 11){
	  if(fabs(bgdEvent.lep2_.Eta()) < 1.479) evt_topNjet1_B[classType]	   += theWeight;
	  else                                   evt_topNjet1_E[classType]	   += theWeight;
	}	  
      } 
      else if(charge != 0 && bgdEvent.jet1Btag_  < 2.1 && (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto){
    	vhselNjet1_num[classType]       += theWeight;
    	vhselNjet1_num_error[classType] += theWeight*theWeight;
      }
    }

  } // end background loop

  //*******************************************************************************
  //Data Events
  //*******************************************************************************
  int nData=dataEvent.tree_->GetEntries();
  for (int i=0; i<nData; ++i) {

    if (i%10000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nData);
    dataEvent.tree_->GetEntry(i);

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
    if(category == 1) lId = ((dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2LooseEleV1)    == SmurfTree::Lep2LooseEleV1   ) ||
    			    ((dataEvent.cuts_ & SmurfTree::Lep1LooseEleV1)    == SmurfTree::Lep1LooseEleV1    && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection);
    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue;
    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;
    if(dataEvent.lep1_.Pt() <= 20.0) continue;
    if(dataEvent.lep2_.Pt() <= 10.0) continue;
    if(dataEvent.dilep_.M() <= 12) continue;

    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    double minmet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
    bool passMETA = minmet > 20. &&
         	   (minmet > 37.+dataEvent.nvtx_/2.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    bool passMETB = minmet > 20. &&
         	   (minmet > 22.+dataEvent.nvtx_/2.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);

    double theWeight = 1.0;
    int classType = kData;
    bool isTopSel = kFALSE; bool isZSel = kFALSE;
    bool isVHSSSelNjets1 = kFALSE; bool isVHSSSelNjets2 = kFALSE;
    if(
      (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
      (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
      charge == 0 &&
      passMETA == true &&
      (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
      dataEvent.njets_ >= 2 && dataEvent.njets_ <= 3 &&
      TMath::Abs(dataEvent.jet1_.Eta()) < 2.4 &&
      TMath::Abs(dataEvent.jet2_.Eta()) < 2.4 &&
      (dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets &&
      1 == 1
      ) isTopSel = kTRUE;

    if(
      (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
      (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
      fabs(dataEvent.dilep_.M()-91.1876) < 10. && dataEvent.type_ == SmurfTree::ee && 
      1 == 1
      ) isZSel = kTRUE;

    double deltaPhiQQL[3] = {DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep1_.Phi()),DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep2_.Phi()),0.0};
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    if(
      passMETB == true &&
      dataEvent.type_ != SmurfTree::mm &&
     (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
      dataEvent.njets_ >= 2 && dataEvent.njets_ <= 3 &&
      TMath::Abs(dataEvent.jet1_.Eta()) < 2.4 &&
      TMath::Abs(dataEvent.jet2_.Eta()) < 2.4 &&
      (dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets &&
      dataEvent.jet2Btag_ < 2.1 && dataEvent.jet3Btag_ < 2.1 &&
      dataEvent.dilep_.M() < 200 &&
     (dataEvent.jet1_+dataEvent.jet2_).M() > 60. && (dataEvent.jet1_+dataEvent.jet2_).M() < 110. &&
      dataEvent.mt_ > 70. && dataEvent.mt_ < 200. && 
      deltaPhiQQL[2]*180.0/TMath::Pi() < 110.0 &&
      1 == 1
      ) isVHSSSelNjets2 = kTRUE;

    deltaPhiQQL[0] = DeltaPhi(dataEvent.jet1_.Phi(),dataEvent.lep1_.Phi()); deltaPhiQQL[1] = DeltaPhi(dataEvent.jet1_.Phi(),dataEvent.lep2_.Phi());
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    if(
      minmet > 37.+dataEvent.nvtx_/2.0 &&
      dataEvent.type_ != SmurfTree::mm &&
     (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
      dataEvent.njets_ == 1 &&
      TMath::Abs(dataEvent.jet1_.Eta()) < 2.4 &&
      (dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets &&
      dataEvent.dilep_.M() < 200 &&
      dataEvent.mt_ > 70. && dataEvent.mt_ < 200. && 
      deltaPhiQQL[2]*180.0/TMath::Pi() < 85.0 &&
      1 == 1
      ) isVHSSSelNjets1 = kTRUE;

    if(isTopSel == kTRUE && dataEvent.jet2Btag_ >= 2.1 && dataEvent.jet3Btag_ < 2.1){
      btag_highestpt_2j_den[classType]	       += theWeight;
      btag_highestpt_2j_den_error[classType]   += theWeight*theWeight;
      if(dataEvent.jet1Btag_ >= 2.1){
    	btag_highestpt_2j_num[classType]       += theWeight;
    	btag_highestpt_2j_num_error[classType] += theWeight*theWeight;
      }
    }

    if(isZSel == kTRUE){
      if     (fabs(dataEvent.lep1_.Eta()) < 1.479 && fabs(dataEvent.lep2_.Eta()) < 1.479){
        qcharge_zsel_den[2] 	   += 1.0;
        qcharge_zsel_den_error[2]  += 1.0*1.0;
      }
      else if(fabs(dataEvent.lep1_.Eta()) > 1.479 && fabs(dataEvent.lep2_.Eta()) > 1.479){
        qcharge_zsel_den[3] 	   += 1.0;
        qcharge_zsel_den_error[3]  += 1.0*1.0;
      }
      if(charge != 0){
	if     (fabs(dataEvent.lep1_.Eta()) < 1.479 && fabs(dataEvent.lep2_.Eta()) < 1.479){
          qcharge_zsel_num[2] 	     += 1.0;
          qcharge_zsel_num_error[2]  += 1.0*1.0;
	}
	else if(fabs(dataEvent.lep1_.Eta()) > 1.479 && fabs(dataEvent.lep2_.Eta()) > 1.479){
          qcharge_zsel_num[3] 	     += 1.0;
          qcharge_zsel_num_error[3]  += 1.0*1.0;
	}
      }
    }

    if(isVHSSSelNjets2 == kTRUE){
      if     (charge == 0 && dataEvent.jet1Btag_ >= 2.1) {
        vhselNjet2_den[classType]	   += theWeight;
        vhselNjet2_den_error[classType] += theWeight*theWeight;
	if      (dataEvent.type_ == SmurfTree::ee){
	  if(fabs(dataEvent.lep2_.Eta()) < 1.479) evt_topNjet2_B[classType]	   += theWeight;
	  else                                    evt_topNjet2_E[classType]	   += theWeight;
	}
	else if(abs(dataEvent.lid1_) == 11){
	  if(fabs(dataEvent.lep1_.Eta()) < 1.479) evt_topNjet2_B[classType]	   += theWeight;
	  else                                    evt_topNjet2_E[classType]	   += theWeight;
	}
	else if(abs(dataEvent.lid2_) == 11){
	  if(fabs(dataEvent.lep2_.Eta()) < 1.479) evt_topNjet2_B[classType]	   += theWeight;
	  else                                    evt_topNjet2_E[classType]	   += theWeight;
	}	  
      } 
      else if(charge != 0 && dataEvent.jet1Btag_  < 2.1 && (dataEvent.cuts_ & patternTopVeto) == patternTopVeto){
    	vhselNjet2_num[classType]       += theWeight;
    	vhselNjet2_num_error[classType] += theWeight*theWeight;
      }
    }

    if(isVHSSSelNjets1 == kTRUE){
      if     (charge == 0 && dataEvent.jet1Btag_ >= 2.1) {
        vhselNjet1_den[classType]	+= theWeight;
        vhselNjet1_den_error[classType] += theWeight*theWeight;
	if      (dataEvent.type_ == SmurfTree::ee){
	  if(fabs(dataEvent.lep2_.Eta()) < 1.479) evt_topNjet1_B[classType]	   += theWeight;
	  else                                    evt_topNjet1_E[classType]	   += theWeight;
	}
	else if(abs(dataEvent.lid1_) == 11){
	  if(fabs(dataEvent.lep1_.Eta()) < 1.479) evt_topNjet1_B[classType]	   += theWeight;
	  else                                    evt_topNjet1_E[classType]	   += theWeight;
	}
	else if(abs(dataEvent.lid2_) == 11){
	  if(fabs(dataEvent.lep2_.Eta()) < 1.479) evt_topNjet1_B[classType]	   += theWeight;
	  else                                    evt_topNjet1_E[classType]	   += theWeight;
	}	  
      } 
      else if(charge != 0 && dataEvent.jet1Btag_  < 2.1 && (dataEvent.cuts_ & patternTopVeto) == patternTopVeto){
    	vhselNjet1_num[classType]       += theWeight;
    	vhselNjet1_num_error[classType] += theWeight*theWeight;
      }
    }

  } // End loop data

  //*******************************************************************************
  //Print Summary 
  //*******************************************************************************
  //*******************************************************************************
  //2-Jet Bin : BTag Efficiency for highest pt jet
  //*******************************************************************************
  printf("**********eff highest pt jet 2-j**********\n");
  double effttMC_btag_highestpt_2j,effttMC_btag_highestpt_2j_error,effttMC_btag_highestpt_2j_tt,effttMC_btag_highestpt_2j_tt_error;
  double effttDA_btag_highestpt_2j,effttDA_btag_highestpt_2j_error,effttMC_btag_highestpt_2j_tw,effttMC_btag_highestpt_2j_tw_error;

  //MC efficiencies
  effttMC_btag_highestpt_2j_tt = (btag_highestpt_2j_num[1])/(btag_highestpt_2j_den[1]);
  effttMC_btag_highestpt_2j_tw = (btag_highestpt_2j_num[2])/(btag_highestpt_2j_den[2]);
  effttMC_btag_highestpt_2j    = (btag_highestpt_2j_num[1] + btag_highestpt_2j_num[2]) / (btag_highestpt_2j_den[1] + btag_highestpt_2j_den[2]);
  
  effttMC_btag_highestpt_2j_tt_error = sqrt((1.0-effttMC_btag_highestpt_2j_tt)*effttMC_btag_highestpt_2j_tt/(btag_highestpt_2j_den[1])*
  					       (btag_highestpt_2j_den_error[1])/(btag_highestpt_2j_den[1]));	
  effttMC_btag_highestpt_2j_tw_error = sqrt((1.0-effttMC_btag_highestpt_2j_tw)*effttMC_btag_highestpt_2j_tw/(btag_highestpt_2j_den[2])*
  					       (btag_highestpt_2j_den_error[2])/(btag_highestpt_2j_den[2]));
  effttMC_btag_highestpt_2j_error    = sqrt((1.0-effttMC_btag_highestpt_2j)*effttMC_btag_highestpt_2j/(btag_highestpt_2j_den[1]+btag_highestpt_2j_den[2])*
  					       (btag_highestpt_2j_den_error[1]+btag_highestpt_2j_den_error[2])/(btag_highestpt_2j_den[1]+btag_highestpt_2j_den[2]));

  //Data efficiencies
  effttDA_btag_highestpt_2j = (btag_highestpt_2j_num[3]-btag_highestpt_2j_num[0]-btag_highestpt_2j_num[2])/
    (btag_highestpt_2j_den[3]-btag_highestpt_2j_den[0]-btag_highestpt_2j_den[2]);    
  effttDA_btag_highestpt_2j_error = sqrt((1-effttDA_btag_highestpt_2j)*effttDA_btag_highestpt_2j/btag_highestpt_2j_den[3]);

  printf("                 tt          tw           data           bck\n");
  printf("num/den --> %7.2f/%7.2f=%6.3f - %7.2f/%7.2f=%6.3f - %7.2f/%7.2f=%6.3f - %7.2f/%7.2f\n",
         btag_highestpt_2j_num[1],btag_highestpt_2j_den[1],btag_highestpt_2j_num[1]/btag_highestpt_2j_den[1],btag_highestpt_2j_num[2],btag_highestpt_2j_den[2],btag_highestpt_2j_num[2]/btag_highestpt_2j_den[2],
	 btag_highestpt_2j_num[3],btag_highestpt_2j_den[3], btag_highestpt_2j_num[3]/btag_highestpt_2j_den[3],btag_highestpt_2j_num[0],btag_highestpt_2j_den[0]);

  printf("scaleFactor2j --> %6.3f +/- %6.3f - %6.3f +/- %6.3f ==> %6.3f +/- %6.3f\n",
         effttMC_btag_highestpt_2j,effttMC_btag_highestpt_2j_error,
	 effttDA_btag_highestpt_2j,effttDA_btag_highestpt_2j_error,
	 effttDA_btag_highestpt_2j/effttMC_btag_highestpt_2j,
	 effttDA_btag_highestpt_2j/effttMC_btag_highestpt_2j*
	 sqrt(effttDA_btag_highestpt_2j_error/effttDA_btag_highestpt_2j*effttDA_btag_highestpt_2j_error/effttDA_btag_highestpt_2j+
	      effttMC_btag_highestpt_2j_error/effttMC_btag_highestpt_2j*effttMC_btag_highestpt_2j_error/effttMC_btag_highestpt_2j));
 
  //*******************************************************************************
  //Wrong charge scale factor
  //*******************************************************************************
  printf("**********eff wrong charge**********\n");
  double effZMC_qcharge_zsel[2],effZMC_qcharge_zsel_error[2];
  double effZDA_qcharge_zsel[2],effZDA_qcharge_zsel_error[2];

  //MC efficiencies
  effZMC_qcharge_zsel[0]    = (qcharge_zsel_num[0]) / (qcharge_zsel_den[0]);
  effZMC_qcharge_zsel[1]    = (qcharge_zsel_num[1]) / (qcharge_zsel_den[1]);
  
  effZMC_qcharge_zsel_error[0]    = sqrt((1.0-effZMC_qcharge_zsel[0])*effZMC_qcharge_zsel[0]/(qcharge_zsel_den[0])*
  					       (qcharge_zsel_den_error[0])/(qcharge_zsel_den[0]));
  effZMC_qcharge_zsel_error[1]    = sqrt((1.0-effZMC_qcharge_zsel[1])*effZMC_qcharge_zsel[1]/(qcharge_zsel_den[1])*
  					       (qcharge_zsel_den_error[1])/(qcharge_zsel_den[1]));

  //Data efficiencies
  effZDA_qcharge_zsel[0]    = (qcharge_zsel_num[2]) / (qcharge_zsel_den[2]);
  effZDA_qcharge_zsel[1]    = (qcharge_zsel_num[3]) / (qcharge_zsel_den[3]);
  
  effZDA_qcharge_zsel_error[0]    = sqrt((1.0-effZDA_qcharge_zsel[0])*effZDA_qcharge_zsel[0]/(qcharge_zsel_den[2])*
  					       (qcharge_zsel_den_error[0])/(qcharge_zsel_den[2]));
  effZDA_qcharge_zsel_error[1]    = sqrt((1.0-effZDA_qcharge_zsel[1])*effZDA_qcharge_zsel[1]/(qcharge_zsel_den[3])*
  					       (qcharge_zsel_den_error[1])/(qcharge_zsel_den[3]));

  printf("WrongChargeB --> %8.5f +/- %8.5f  - %8.5f +/- %5.5f --> %6.3f +/- %6.3f\n",
             effZDA_qcharge_zsel[0],effZDA_qcharge_zsel_error[0],effZMC_qcharge_zsel[0],effZMC_qcharge_zsel_error[0],
	     effZDA_qcharge_zsel[0]/effZMC_qcharge_zsel[0],
	     effZDA_qcharge_zsel[0]/effZMC_qcharge_zsel[0]*
	     sqrt(effZDA_qcharge_zsel_error[0]/effZDA_qcharge_zsel[0]*effZDA_qcharge_zsel_error[0]/effZDA_qcharge_zsel[0]+
	          effZMC_qcharge_zsel_error[0]/effZMC_qcharge_zsel[0]*effZMC_qcharge_zsel_error[0]/effZMC_qcharge_zsel[0]));
  printf("WrongChargeE --> %8.5f +/- %8.5f  - %8.5f +/- %5.5f --> %6.3f +/- %6.3f\n",
             effZDA_qcharge_zsel[1],effZDA_qcharge_zsel_error[1],effZMC_qcharge_zsel[1],effZMC_qcharge_zsel_error[1],
	     effZDA_qcharge_zsel[1]/effZMC_qcharge_zsel[1],
	     effZDA_qcharge_zsel[1]/effZMC_qcharge_zsel[1]*
	     sqrt(effZDA_qcharge_zsel_error[1]/effZDA_qcharge_zsel[1]*effZDA_qcharge_zsel_error[1]/effZDA_qcharge_zsel[1]+
	          effZMC_qcharge_zsel_error[1]/effZMC_qcharge_zsel[1]*effZMC_qcharge_zsel_error[1]/effZMC_qcharge_zsel[1]));
  //*******************************************************************************
  //VH background
  //*******************************************************************************
  bool useDefaultWrongCharge = kFALSE;
  if(useDefaultWrongCharge == kTRUE){
    effZDA_qcharge_zsel[0] = 0.00303; effZDA_qcharge_zsel_error[0] = 0.00014;
    effZMC_qcharge_zsel[0] = 0.00359; effZMC_qcharge_zsel_error[0] = 0.00004;
    effZDA_qcharge_zsel[1] = 0.02426; effZDA_qcharge_zsel_error[1] = 0.00121;
    effZMC_qcharge_zsel[1] = 0.02100; effZMC_qcharge_zsel_error[1] = 0.00024;
  }

  bool doNjets2 = kTRUE;
  if(doNjets2 == kTRUE){
  printf("**********VH background Njet2**********\n");
  printf("                                  tt                 tw                   data                bck\n");
  printf("evt_top B/E --> %6.3f/%6.3f/%6.3f - %6.3f/%6.3f/%6.3f - %6.3f/%6.3f/%6.3f - %6.3f/%6.3f/%6.3f\n",evt_topNjet2_B[1],evt_topNjet2_E[1],evt_topNjet2_B[1]/(evt_topNjet2_E[1]+evt_topNjet2_B[1]),
        evt_topNjet2_B[2],evt_topNjet2_E[2],evt_topNjet2_B[2]/(evt_topNjet2_E[2]+evt_topNjet2_B[2]), evt_topNjet2_B[3],evt_topNjet2_E[3],evt_topNjet2_B[3]/(evt_topNjet2_E[3]+evt_topNjet2_B[3]), evt_topNjet2_B[0],evt_topNjet2_E[0],evt_topNjet2_B[0]/(evt_topNjet2_E[0]+evt_topNjet2_B[0]));

  double evt_topNjet2_B_avg = (evt_topNjet2_B[1]+evt_topNjet2_B[2])/(evt_topNjet2_B[1]+evt_topNjet2_B[2]+evt_topNjet2_E[1]+evt_topNjet2_E[2]);
  double effZMC_qcharge_zsel_avg = evt_topNjet2_B_avg*effZMC_qcharge_zsel[0]+(1.0-evt_topNjet2_B_avg)*effZMC_qcharge_zsel[1];
  double effZDA_qcharge_zsel_avg = evt_topNjet2_B_avg*effZDA_qcharge_zsel[0]+(1.0-evt_topNjet2_B_avg)*effZDA_qcharge_zsel[1];
  double effZMC_qcharge_zsel_error_avg = sqrt(TMath::Power(evt_topNjet2_B_avg*effZMC_qcharge_zsel_error[0],2)+TMath::Power((1-evt_topNjet2_B_avg)*effZMC_qcharge_zsel_error[1],2));
  double effZDA_qcharge_zsel_error_avg = sqrt(TMath::Power(evt_topNjet2_B_avg*effZDA_qcharge_zsel_error[0],2)+TMath::Power((1-evt_topNjet2_B_avg)*effZDA_qcharge_zsel_error[1],2));
  printf("Wrong charge average  --> (DA) %7.5f +/- %7.5f  - (MC) %7.5f +/- %7.5f --> %6.4f +/- %6.4f\n",effZDA_qcharge_zsel_avg,effZDA_qcharge_zsel_error_avg,
  effZMC_qcharge_zsel_avg,effZMC_qcharge_zsel_error_avg,effZDA_qcharge_zsel_avg/effZMC_qcharge_zsel_avg,effZDA_qcharge_zsel_avg/effZMC_qcharge_zsel_avg*
  sqrt(TMath::Power(effZDA_qcharge_zsel_error_avg/effZDA_qcharge_zsel_avg,2)+TMath::Power(effZMC_qcharge_zsel_error_avg/effZMC_qcharge_zsel_avg,2)));

  //vhselNjet2_num[1]+=0.130;
  printf("                          tt                tw                data                 bck\n");
  printf("top-like btag,q1*q2<0  --> %8.3f +/- %6.3f  - %8.3f +/- %6.3f - %8.3f +/ %6.3f - %8.3f +/- %6.3f\n",
         vhselNjet2_den[1],sqrt(vhselNjet2_den_error[1]),vhselNjet2_den[2],sqrt(vhselNjet2_den_error[2]),
	 vhselNjet2_den[3],sqrt(vhselNjet2_den_error[3]),vhselNjet2_den[0],sqrt(vhselNjet2_den_error[0]));
  printf("vh-like nobtag,q1*q2>0 --> %8.3f +/- %6.3f  - %8.3f +/- %6.3f - %8.3f +/ %6.3f - %8.3f +/- %6.3f\n",
         vhselNjet2_num[1],sqrt(vhselNjet2_num_error[1]),vhselNjet2_num[2],sqrt(vhselNjet2_num_error[2]),
	 vhselNjet2_num[3],sqrt(vhselNjet2_num_error[3]),vhselNjet2_num[0],sqrt(vhselNjet2_num_error[0]));
  double extrapolateVH_MC = (vhselNjet2_den[1]+vhselNjet2_den[2])*(1.0-effttMC_btag_highestpt_2j)/effttMC_btag_highestpt_2j*effZMC_qcharge_zsel_avg/(1.0-effZMC_qcharge_zsel_avg);

  double extrapolateVH_MC_error = (vhselNjet2_den_error[1]+vhselNjet2_den_error[2])/(vhselNjet2_den[1]+vhselNjet2_den[2])*extrapolateVH_MC+
  (vhselNjet2_den[1]+vhselNjet2_den[2])*effZMC_qcharge_zsel_avg/(1.0-effZMC_qcharge_zsel_avg)*effttMC_btag_highestpt_2j_error/effttMC_btag_highestpt_2j/effttMC_btag_highestpt_2j+
  (vhselNjet2_den[1]+vhselNjet2_den[2])*(1.0-effttMC_btag_highestpt_2j)/effttMC_btag_highestpt_2j*effZMC_qcharge_zsel_error_avg/(1.0-effZMC_qcharge_zsel_avg)/(1.0-effZMC_qcharge_zsel_avg);

  double extrapolateVH_DA = (vhselNjet2_den[3]-vhselNjet2_den[0])*(1.0-effttDA_btag_highestpt_2j)/effttDA_btag_highestpt_2j*effZDA_qcharge_zsel_avg/(1.0-effZDA_qcharge_zsel_avg);

  double extrapolateVH_DA_error = (vhselNjet2_den_error[1]+vhselNjet2_den_error[2])/(vhselNjet2_den[1]+vhselNjet2_den[2])*extrapolateVH_DA+
  (vhselNjet2_den[3]-vhselNjet2_den[0])*effZDA_qcharge_zsel_avg/(1.0-effZDA_qcharge_zsel_avg)*effttDA_btag_highestpt_2j_error/effttDA_btag_highestpt_2j/effttDA_btag_highestpt_2j+
  (vhselNjet2_den[3]-vhselNjet2_den[0])*(1.0-effttDA_btag_highestpt_2j)/effttDA_btag_highestpt_2j*effZDA_qcharge_zsel_error_avg/(1.0-effZDA_qcharge_zsel_avg)/(1.0-effZDA_qcharge_zsel_avg);

  printf("\nvh-like extrapolation  --> (MCest) %8.3f +/- %6.3f  - (DA) %8.3f +/- %6.3f --> %8.3f +/- %6.3f\n",
         extrapolateVH_MC,extrapolateVH_MC_error,extrapolateVH_DA,extrapolateVH_DA_error,
	 extrapolateVH_DA/extrapolateVH_MC,extrapolateVH_DA/extrapolateVH_MC*extrapolateVH_DA_error/extrapolateVH_DA);
  printf("  vh-like MC closure     --> (MC) %8.3f +/- %6.3f  vs. (MCest) %8.3f +/- %6.3f --> %8.3f +/- %6.3f\n",
         vhselNjet2_num[1]+vhselNjet2_num[2],sqrt(vhselNjet2_num_error[1]+vhselNjet2_num_error[2]),extrapolateVH_MC,extrapolateVH_MC_error,(vhselNjet2_num[1]+vhselNjet2_num[2])/extrapolateVH_MC,
	(vhselNjet2_num[1]+vhselNjet2_num[2])/extrapolateVH_MC*sqrt(TMath::Power(sqrt(vhselNjet2_num_error[1]+vhselNjet2_num_error[2])/(vhselNjet2_num[1]+vhselNjet2_num[2]),2)+
	                                                  TMath::Power(extrapolateVH_MC_error/extrapolateVH_MC,2)));
  printf("  final  extrapolation   --> (DA) %8.3f +/- %6.3f  vs. (MC) %8.3f +/- %6.3f --> %8.3f +/- %6.3f\n",
         extrapolateVH_DA,extrapolateVH_DA_error,vhselNjet2_num[1]+vhselNjet2_num[2],sqrt(vhselNjet2_num_error[1]+vhselNjet2_num_error[2]),extrapolateVH_DA/(vhselNjet2_num[1]+vhselNjet2_num[2]),
	 extrapolateVH_DA/(vhselNjet2_num[1]+vhselNjet2_num[2])*sqrt(TMath::Power(sqrt(vhselNjet2_num_error[1]+vhselNjet2_num_error[2])/(vhselNjet2_num[1]+vhselNjet2_num[2]),2)+
	                                                  TMath::Power(extrapolateVH_DA_error/extrapolateVH_DA,2)));
  }  
  bool doNjets1 = kTRUE;
  if(doNjets1 == kTRUE){
  printf("**********VH background Njet1**********\n");
  printf("                                  tt                 tw                   data                bck\n");
  printf("evt_top B/E --> %6.3f/%6.3f/%6.3f - %6.3f/%6.3f/%6.3f - %6.3f/%6.3f/%6.3f - %6.3f/%6.3f/%6.3f\n",evt_topNjet1_B[1],evt_topNjet1_E[1],evt_topNjet1_B[1]/(evt_topNjet1_E[1]+evt_topNjet1_B[1]),
        evt_topNjet1_B[2],evt_topNjet1_E[2],evt_topNjet1_B[2]/(evt_topNjet1_E[2]+evt_topNjet1_B[2]), evt_topNjet1_B[3],evt_topNjet1_E[3],evt_topNjet1_B[3]/(evt_topNjet1_E[3]+evt_topNjet1_B[3]), evt_topNjet1_B[0],evt_topNjet1_E[0],evt_topNjet1_B[0]/(evt_topNjet1_E[0]+evt_topNjet1_B[0]));

  double evt_topNjet1_B_avg = (evt_topNjet1_B[1]+evt_topNjet1_B[2])/(evt_topNjet1_B[1]+evt_topNjet1_B[2]+evt_topNjet1_E[1]+evt_topNjet1_E[2]);
  double effZMC_qcharge_zsel_avg = evt_topNjet1_B_avg*effZMC_qcharge_zsel[0]+(1.0-evt_topNjet1_B_avg)*effZMC_qcharge_zsel[1];
  double effZDA_qcharge_zsel_avg = evt_topNjet1_B_avg*effZDA_qcharge_zsel[0]+(1.0-evt_topNjet1_B_avg)*effZDA_qcharge_zsel[1];
  double effZMC_qcharge_zsel_error_avg = sqrt(TMath::Power(evt_topNjet1_B_avg*effZMC_qcharge_zsel_error[0],2)+TMath::Power((1-evt_topNjet1_B_avg)*effZMC_qcharge_zsel_error[1],2));
  double effZDA_qcharge_zsel_error_avg = sqrt(TMath::Power(evt_topNjet1_B_avg*effZDA_qcharge_zsel_error[0],2)+TMath::Power((1-evt_topNjet1_B_avg)*effZDA_qcharge_zsel_error[1],2));
  printf("Wrong charge average  --> (DA) %7.5f +/- %7.5f  - (MC) %7.5f +/- %7.5f --> %6.4f +/- %6.4f\n",effZDA_qcharge_zsel_avg,effZDA_qcharge_zsel_error_avg,
  effZMC_qcharge_zsel_avg,effZMC_qcharge_zsel_error_avg,effZDA_qcharge_zsel_avg/effZMC_qcharge_zsel_avg,effZDA_qcharge_zsel_avg/effZMC_qcharge_zsel_avg*
  sqrt(TMath::Power(effZDA_qcharge_zsel_error_avg/effZDA_qcharge_zsel_avg,2)+TMath::Power(effZMC_qcharge_zsel_error_avg/effZMC_qcharge_zsel_avg,2)));

  //vhselNjet1_num[1]+=0.488;
  printf("                          tt                tw                data                 bck\n");
  printf("top-like btag,q1*q2<0  --> %8.3f +/- %6.3f  - %8.3f +/- %6.3f - %8.3f +/ %6.3f - %8.3f +/- %6.3f\n",
         vhselNjet1_den[1],sqrt(vhselNjet1_den_error[1]),vhselNjet1_den[2],sqrt(vhselNjet1_den_error[2]),
	 vhselNjet1_den[3],sqrt(vhselNjet1_den_error[3]),vhselNjet1_den[0],sqrt(vhselNjet1_den_error[0]));
  printf("vh-like nobtag,q1*q2>0 --> %8.3f +/- %6.3f  - %8.3f +/- %6.3f - %8.3f +/ %6.3f - %8.3f +/- %6.3f\n",
         vhselNjet1_num[1],sqrt(vhselNjet1_num_error[1]),vhselNjet1_num[2],sqrt(vhselNjet1_num_error[2]),
	 vhselNjet1_num[3],sqrt(vhselNjet1_num_error[3]),vhselNjet1_num[0],sqrt(vhselNjet1_num_error[0]));
  double extrapolateVH_MC = (vhselNjet1_den[1]+vhselNjet1_den[2])*(1.0-effttMC_btag_highestpt_2j)/effttMC_btag_highestpt_2j*effZMC_qcharge_zsel_avg/(1.0-effZMC_qcharge_zsel_avg);

  double extrapolateVH_MC_error = (vhselNjet1_den_error[1]+vhselNjet1_den_error[2])/(vhselNjet1_den[1]+vhselNjet1_den[2])*extrapolateVH_MC+
  (vhselNjet1_den[1]+vhselNjet1_den[2])*effZMC_qcharge_zsel_avg/(1.0-effZMC_qcharge_zsel_avg)*effttMC_btag_highestpt_2j_error/effttMC_btag_highestpt_2j/effttMC_btag_highestpt_2j+
  (vhselNjet1_den[1]+vhselNjet1_den[2])*(1.0-effttMC_btag_highestpt_2j)/effttMC_btag_highestpt_2j*effZMC_qcharge_zsel_error_avg/(1.0-effZMC_qcharge_zsel_avg)/(1.0-effZMC_qcharge_zsel_avg);

  double extrapolateVH_DA = (vhselNjet1_den[3]-vhselNjet1_den[0])*(1.0-effttDA_btag_highestpt_2j)/effttDA_btag_highestpt_2j*effZDA_qcharge_zsel_avg/(1.0-effZDA_qcharge_zsel_avg);

  double extrapolateVH_DA_error = (vhselNjet1_den_error[1]+vhselNjet1_den_error[2])/(vhselNjet1_den[1]+vhselNjet1_den[2])*extrapolateVH_DA+
  (vhselNjet1_den[3]-vhselNjet1_den[0])*effZDA_qcharge_zsel_avg/(1.0-effZDA_qcharge_zsel_avg)*effttDA_btag_highestpt_2j_error/effttDA_btag_highestpt_2j/effttDA_btag_highestpt_2j+
  (vhselNjet1_den[3]-vhselNjet1_den[0])*(1.0-effttDA_btag_highestpt_2j)/effttDA_btag_highestpt_2j*effZDA_qcharge_zsel_error_avg/(1.0-effZDA_qcharge_zsel_avg)/(1.0-effZDA_qcharge_zsel_avg);

  printf("\nvh-like extrapolation  --> (MCest) %8.3f +/- %6.3f  - (DA) %8.3f +/- %6.3f --> %8.3f +/- %6.3f\n",
         extrapolateVH_MC,extrapolateVH_MC_error,extrapolateVH_DA,extrapolateVH_DA_error,
	 extrapolateVH_DA/extrapolateVH_MC,extrapolateVH_DA/extrapolateVH_MC*extrapolateVH_DA_error/extrapolateVH_DA);
  printf("  vh-like MC closure     --> (MC) %8.3f +/- %6.3f  vs. (MCest) %8.3f +/- %6.3f --> %8.3f +/- %6.3f\n",
         vhselNjet1_num[1]+vhselNjet1_num[2],sqrt(vhselNjet1_num_error[1]+vhselNjet1_num_error[2]),extrapolateVH_MC,extrapolateVH_MC_error,(vhselNjet1_num[1]+vhselNjet1_num[2])/extrapolateVH_MC,
	(vhselNjet1_num[1]+vhselNjet1_num[2])/extrapolateVH_MC*sqrt(TMath::Power(sqrt(vhselNjet1_num_error[1]+vhselNjet1_num_error[2])/(vhselNjet1_num[1]+vhselNjet1_num[2]),2)+
	                                                  TMath::Power(extrapolateVH_MC_error/extrapolateVH_MC,2)));
  printf("  final  extrapolation   --> (DA) %8.3f +/- %6.3f  vs. (MC) %8.3f +/- %6.3f --> %8.3f +/- %6.3f\n",
         extrapolateVH_DA,extrapolateVH_DA_error,vhselNjet1_num[1]+vhselNjet1_num[2],sqrt(vhselNjet1_num_error[1]+vhselNjet1_num_error[2]),extrapolateVH_DA/(vhselNjet1_num[1]+vhselNjet1_num[2]),
	 extrapolateVH_DA/(vhselNjet1_num[1]+vhselNjet1_num[2])*sqrt(TMath::Power(sqrt(vhselNjet1_num_error[1]+vhselNjet1_num_error[2])/(vhselNjet1_num[1]+vhselNjet1_num[2]),2)+
	                                                  TMath::Power(extrapolateVH_DA_error/extrapolateVH_DA,2)));
  }  
  return;
}
