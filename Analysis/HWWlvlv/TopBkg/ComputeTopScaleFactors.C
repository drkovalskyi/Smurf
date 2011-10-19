#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"
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
void ComputeTopScaleFactors
(
 TString bgdInputFile    = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/backgroundA_skim1.root",
 TString dataInputFile   = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/data_2l_skim1.root"
 )
{

  //*******************************************************************************
  //Settings 
  //*******************************************************************************
  enum { kOther, kTTBAR, kTW, kData };
  double lumi = 2.121;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);


  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  //*******************************************************************************
  //Load Auxiliary Inputs
  //*******************************************************************************
  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/LP2011/auxiliar/efficiency_results_v6_42x.root");
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

  LeptonScaleLookup trigLookup("/data/smurf/data/LP2011/auxiliar/efficiency_results_v6_42x.root");
  
  TFile *fNvtxFile = TFile::Open("/data/smurf/data/LP2011/auxiliar/puReweighting.root");
  TH1D *fhDNvtx = (TH1D*)(fNvtxFile->Get("puWeights"));
  assert(fhDNvtx);
  fhDNvtx->SetDirectory(0);
  fNvtxFile->Close();
  delete fNvtxFile;

  TFile *fPUS3File = TFile::Open("/data/smurf/data/LP2011/auxiliar/puWeights_PU3_68mb.root");
  TH1D *fhDPUS3 = (TH1D*)(fPUS3File->Get("puWeights"));
  assert(fhDPUS3);
  fhDPUS3->SetDirectory(0);
  delete fPUS3File;

  TFile *fPUS4File = TFile::Open("/data/smurf/data/LP2011/auxiliar/puWeights_PU4_68mb.root");
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;


  //*******************************************************************************
  //Yields and Histograms
  //*******************************************************************************
  vector<vector<double> >  btag_lowpt_1j_den;
  vector<vector<double> > btag_lowpt_1j_num;
  vector<vector<double> > btag_lowpt_0j_den;
  vector<vector<double> > btag_lowpt_0j_num;
  vector<vector<double> >  btag_highestpt_2j_den,btag_highestpt_2j_num,btag_highestpt_1j_den,btag_highestpt_1j_num;
  vector<vector<double> >  btag_lowpt_1j_den_error,btag_lowpt_1j_num_error,btag_lowpt_0j_den_error,btag_lowpt_0j_num_error;
  vector<vector<double> >  btag_highestpt_2j_den_error,btag_highestpt_2j_num_error,btag_highestpt_1j_den_error,btag_highestpt_1j_num_error;

  for(int i=0; i<4; i++){
    vector<double> tmpbtag_lowpt_1j_den,tmpbtag_lowpt_1j_num,tmpbtag_lowpt_0j_den,tmpbtag_lowpt_0j_num;
    vector<double> tmpbtag_highestpt_2j_den,tmpbtag_highestpt_2j_num,tmpbtag_highestpt_1j_den,tmpbtag_highestpt_1j_num;
    vector<double> tmpbtag_lowpt_1j_den_error,tmpbtag_lowpt_1j_num_error,tmpbtag_lowpt_0j_den_error,tmpbtag_lowpt_0j_num_error;
    vector<double> tmpbtag_highestpt_2j_den_error,tmpbtag_highestpt_2j_num_error,tmpbtag_highestpt_1j_den_error,tmpbtag_highestpt_1j_num_error;
    for(int j=0; j<5; j++){
      tmpbtag_lowpt_1j_den.push_back(0),tmpbtag_lowpt_1j_num.push_back(0),tmpbtag_lowpt_0j_den.push_back(0),tmpbtag_lowpt_0j_num.push_back(0);
      tmpbtag_highestpt_2j_den.push_back(0),tmpbtag_highestpt_2j_num.push_back(0),tmpbtag_highestpt_1j_den.push_back(0),tmpbtag_highestpt_1j_num.push_back(0);
      tmpbtag_lowpt_1j_den_error.push_back(0),tmpbtag_lowpt_1j_num_error.push_back(0),tmpbtag_lowpt_0j_den_error.push_back(0),tmpbtag_lowpt_0j_num_error.push_back(0);
      tmpbtag_highestpt_2j_den_error.push_back(0),tmpbtag_highestpt_2j_num_error.push_back(0),tmpbtag_highestpt_1j_den_error.push_back(0),tmpbtag_highestpt_1j_num_error.push_back(0);
    }

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


  unsigned int patternTopVeto         = SmurfTree::TopVeto;
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

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::wjets           ) fDecay = 3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar           ) fDecay = 5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee            ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm            ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt            ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::qcd             ) fDecay = 11;
    else if(bgdEvent.dstype_ == SmurfTree::tw              ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::qqww            ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz              ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz              ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww            ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma          ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::data            ) fDecay =  1;
    else {fDecay = 0;cout << bgdEvent.dstype_ << endl;assert(0);}


    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);

    double minmet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool passMET = minmet > 20. &&
      (minmet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);


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
	add = add*nPUScaleFactor(fhDPUS4,bgdEvent.npu_);
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
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven) {
      //do not apply any MC->data scale factors for data driven dy->tautau 
      theWeight = bgdEvent.scale1fb_*lumi;
    }
    else if(bgdEvent.dstype_ == SmurfTree::qcd) { //Z->tautau embedded sample
      theWeight = 0.0189712*lumi;
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
      if(fDecay == 3) {
        //for wjets MC
        add=add*1.95; 
      }
      theWeight = bgdEvent.scale1fb_*lumi*add;
    }

    if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
       bgdEvent.dstype_ != SmurfTree::data) &&
      bgdEvent.dilep_.M()   > 12 &&
      (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
      charge == 0 &&
      bgdEvent.lep1_.Pt() > 20. &&
      bgdEvent.lep2_.Pt() > 10. &&
      passMET == true &&
      (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
      (bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
      1 == 1
      ){
      int classType = kOther;
      if(fDecay == 5 ) classType = kTTBAR;
      if(fDecay == 13) classType = kTW;


      if(bgdEvent.jet1Btag_ >= 2.1 && bgdEvent.njets_ == 1){
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

      if((bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && bgdEvent.jet2Btag_ >= 2.1 && bgdEvent.njets_ == 2){
        btag_highestpt_2j_den[classType][4] 	     += theWeight;
        btag_highestpt_2j_den[classType][bgdEvent.type_] += theWeight;
        btag_highestpt_2j_den_error[classType][4] 	           += theWeight*theWeight;
        btag_highestpt_2j_den_error[classType][bgdEvent.type_] += theWeight*theWeight;
        if(bgdEvent.jet1Btag_ >= 2.1){
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
        if(bgdEvent.jet1Btag_ >= 2.1){
          btag_highestpt_1j_num[classType][4]	       += theWeight;
          btag_highestpt_1j_num[classType][bgdEvent.type_] += theWeight;
          btag_highestpt_1j_num_error[classType][4]	             += theWeight*theWeight;
          btag_highestpt_1j_num_error[classType][bgdEvent.type_] += theWeight*theWeight;
        }
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

    if((dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue;

    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    int Njet3 = 0;
    if(dataEvent.jet3_.Pt() <= 30)									     Njet3 = 0;
    else if(dataEvent.jet3_.Pt() > 30 && (
      (dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta() < 0) ||
      (dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta() < 0))) Njet3 = 2;
    else												     Njet3 = 1;
    int centrality = 0;
    if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
       ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 


    double minmet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
    bool passMET = minmet > 20. &&
      (minmet > 40. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);


    double theWeight = 1.0;
    if(
      (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
      (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
      dataEvent.dilep_.M()   > 12 &&
      (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
      charge == 0 &&
      dataEvent.lep1_.Pt() > 20. &&
      dataEvent.lep2_.Pt() > 10. &&
      passMET == true &&
      (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
      (dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
      1 == 1
      ) {
      int classType = kData;
      if(dataEvent.jet1Btag_ >= 2.1 && dataEvent.njets_ == 1){
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
      if((dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && dataEvent.jet2Btag_ >= 2.1 && dataEvent.njets_ == 2){
        btag_highestpt_2j_den[classType][4] 	      += theWeight;
        btag_highestpt_2j_den[classType][dataEvent.type_] += theWeight;
        btag_highestpt_2j_den_error[classType][4] 	            += theWeight*theWeight;
        btag_highestpt_2j_den_error[classType][dataEvent.type_] += theWeight*theWeight;
        if(dataEvent.jet1Btag_ >= 2.1){
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
        if(dataEvent.jet1Btag_ >= 2.1){
          btag_highestpt_1j_num[classType][4]	        += theWeight;
          btag_highestpt_1j_num[classType][dataEvent.type_] += theWeight;
          btag_highestpt_1j_num_error[classType][4]	              += theWeight*theWeight;
          btag_highestpt_1j_num_error[classType][dataEvent.type_] += theWeight*theWeight;
        }
      }
    }


  } // End loop data


  //*******************************************************************************
  //Print Summary 
  //*******************************************************************************
  char *classLabel[5] = {"mm ", "me ", "em " , "ee ", "all"};

  double btagSF = 1.0;

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


  for(int i=0; i<5; i++) {
    printf("numerator  (%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_highestpt_2j_num[3][i],btag_highestpt_2j_num[0][i],(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),btag_highestpt_2j_num[1][i],btag_highestpt_2j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_highestpt_2j_den[3][i],btag_highestpt_2j_den[0][i],(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i]),btag_highestpt_2j_den[1][i],btag_highestpt_2j_den[2][i]);
  }
  printf("\n\n");
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
  printf("\n\n");
  printf("****************************************************************************************************************************************\n");
  printf("\n\n");
  printf("\n\n");
    

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
  printf("\n\n");
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


  printf("\n\n");
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

  printf("\n\n");
  printf("**********************************************************\n");
  printf("\n\n");



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

  
  printf("\n\n");
  for(int i=0; i<5; i++) {
    printf("numerator(%s)   --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_lowpt_1j_num[3][i],btag_lowpt_1j_num[0][i],(btag_lowpt_1j_num[1][i]+btag_lowpt_1j_num[2][i]),btag_lowpt_1j_num[1][i],btag_lowpt_1j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i],btag_lowpt_1j_den[3][i],btag_lowpt_1j_den[0][i],(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]),btag_lowpt_1j_den[1][i],btag_lowpt_1j_den[2][i]);
  }
  printf("\n\n");
  printf("channel       eff_tttw         eff_tt            eff_tw           eff_data              \n");
  for(int i=0; i<5; i++) {
    printf("eff (%d): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f\n",i,
           effttMC_btag_lowpt_1j[i]   ,effttMC_btag_lowpt_1j_error[i],effttMC_btag_lowpt_tt_1j[i],effttMC_btag_lowpt_tt_1j_error[i],
           effttMC_btag_lowpt_tw_1j[i],effttMC_btag_lowpt_tw_1j_error[i],effttDA_btag_lowpt_1j[i],effttDA_btag_lowpt_1j_error[i]);
  }

  printf("\n\n");
  printf("**********************************************************\n");
  printf("\n\n");



  //*******************************************************************************
  //Closure Test 0-Jet Bin
  //*******************************************************************************

  printf("**********eff low pt jet 0-j**********\n");
  double ftw_b[5]; 
  for(int i=0; i<5; i++) {
    ftw_b[i] = effttMC_btag_lowpt_tw_1j[i];
    printf("ftw_b(%d) = %5.3f \n",i,ftw_b[i]);
  }
  printf("\n");
  printf("\n");


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
  printf("\n");
  printf("\n");

    
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
  printf("\n");
  printf("\n");

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
    printf("(%s)              %2d               %6.3f                %6.3f             %6.3f +/- %6.3f    %6.3f +/- %6.3f\n",
           classLabel[i],
           (int)btag_lowpt_0j_num[3][i],btag_lowpt_0j_num[0][i],N_top_expected_0j[i],
           estimationMC_btag_lowpt_0j[i],estimationMC_btag_lowpt_0j_error[i],
           estimationDA_btag_lowpt_0j[i],estimationDA_btag_lowpt_0j_error[i]);
  }

  double TopBkgScaleFactor_0Jet = estimationDA_btag_lowpt_0j[4] / N_top_expected_0j[4];
  double TopBkgScaleFactorUncertainty_0Jet = estimationDA_btag_lowpt_0j_error[4] / N_top_expected_0j[4];

  ofstream outf("TopBkgScaleFactors.h");

  outf << "Double_t TopBkgScaleFactor(Int_t jetBin) {" << endl;
  outf << "  assert(jetBin >=0 && jetBin <= 2);" << endl;
  
  outf << "  Double_t TopBkgScaleFactor[3] = { " 
       << TopBkgScaleFactor_0Jet << ", "
       << TopBkgScaleFactor_1Jet << ", "
       << 1.188 << " "
       << " };" << endl;
  outf << "  return TopBkgScaleFactor[jetBin];" << endl;
  outf << "}" << endl;
  outf << endl;


  outf << "Double_t TopBkgScaleFactorKappa(Int_t jetBin) {" << endl;
  
  outf << "  assert(jetBin >=0 && jetBin <= 2);" << endl;
  outf << "  Double_t TopBkgScaleFactorKappa[3] = { " 
       << (1.0 + TopBkgScaleFactorUncertainty_0Jet/TopBkgScaleFactor_0Jet) << ", "
       << (1.0 + TopBkgScaleFactorUncertainty_1Jet/TopBkgScaleFactor_1Jet) << ", "
       << (1.0 + 0.50) << " "
       << " };" << endl;

  outf << "  return TopBkgScaleFactorKappa[jetBin];" << endl;
  outf << "}" << endl;
  outf << endl;

  outf.close();

  return;
}
