/*
   This script produces the boost histograms for the Matrix Element calculation based on the smurf ntuples
   boostproducer(smurfFDir, outputFDir);
   0 - location of the smurf ntuples
   1 - output directory
   ==
   run by root -l boostproducer.C+
   it creates a file at outputFDir directory called "Boost.root"
   ==

*/

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRint.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCut.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

/// DON'T CHANGE ORDER
enum Selection {
  BaseLine           = 1UL<<0,  // pt(reco)>20/10, acceptance,!STA muon, mll>12
  ChargeMatch        = 1UL<<1,  // q1*q2<0
  Lep1FullSelection  = 1UL<<2,  // full id, isolation, d0, dz etc
  Lep1LooseEleV1     = 1UL<<3,  // electron fakeable object selection is passed V1
  Lep1LooseEleV2     = 1UL<<4,  // electron fakeable object selection is passed V2
  Lep1LooseEleV3     = 1UL<<5,  // electron fakeable object selection is passed V3
  Lep1LooseEleV4     = 1UL<<6,  // electron fakeable object selection is passed V4
  Lep1LooseMuV1      = 1UL<<7,  // muon fakeable object selection (relIso<1.0)
  Lep1LooseMuV2      = 1UL<<8,  // muon fakeable object selection (relIso<0.4)
  Lep2FullSelection  = 1UL<<9,  // full id, isolation, d0, dz etc
  Lep2LooseEleV1     = 1UL<<10, // electron fakeable object selection is passed V1
  Lep2LooseEleV2     = 1UL<<11, // electron fakeable object selection is passed V2
  Lep2LooseEleV3     = 1UL<<12, // electron fakeable object selection is passed V3
  Lep2LooseEleV4     = 1UL<<13, // electron fakeable object selection is passed V4
  Lep2LooseMuV1      = 1UL<<14, // muon fakeable object selection (relIso<1.0)
  Lep2LooseMuV2      = 1UL<<15, // muon fakeable object selection (relIso<0.4)
  FullMET            = 1UL<<16, // full met selection
  ZVeto              = 1UL<<17, // event is not in the Z-mass peak for ee/mm final states
  TopTag             = 1UL<<18, // soft muon and b-jet tagging for the whole event regardless of n-jets (non-zero means tagged)
  TopVeto            = 1UL<<19, // soft muon and b-jet tagging for the whole event regardless of n-jets (zero means tagged)
  OneBJet            = 1UL<<20, // 1-jet events, where the jet is b-tagged (top control sample with one b-quark missing)
  TopTagNotInJets    = 1UL<<21, // soft muon and b-jet tagging for areas outside primary jets (non-zero means tagged)
  ExtraLeptonVeto    = 1UL<<22, // extra lepton veto, DR(muon-electron)>=0.3
  Lep3FullSelection  = 1UL<<23,  // full id, isolation, d0, dz etc
  Lep3LooseEleV1     = 1UL<<24, // electron fakeable object selection is passed V1
  Lep3LooseEleV2     = 1UL<<25, // electron fakeable object selection is passed V2
  Lep3LooseEleV3     = 1UL<<26, // electron fakeable object selection is passed V3
  Lep3LooseEleV4     = 1UL<<27, // electron fakeable object selection is passed V4
  Lep3LooseMuV1      = 1UL<<28, // muon fakeable object selection (relIso<1.0)
  Lep3LooseMuV2      = 1UL<<29, // muon fakeable object selection (relIso<0.4)
  Trigger            = 1UL<<30  // passed a set of triggers
};

UInt_t ww = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ZVeto|TopVeto|ExtraLeptonVeto;
void singleboost(TString inputFileName, TFile* outputFile, TString histName);


//###################
//# main function
//###################


void boostproducer(TString smurfFDir = "/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/", TString outputFDir = "./") 
{
  
  TString outputFileName = Form("%s/Boost.root", outputFDir.Data());
  TFile *outputFile = new TFile(outputFileName, "RECREATE");

  // write boost histograms for the given boost 
  singleboost(Form("%s/qqww.root", smurfFDir.Data()), outputFile, "ww");
  singleboost(Form("%s/hww110.root", smurfFDir.Data()), outputFile, "hww110");
  singleboost(Form("%s/hww115.root", smurfFDir.Data()), outputFile, "hww115");
  singleboost(Form("%s/hww120.root", smurfFDir.Data()), outputFile, "hww120");
  singleboost(Form("%s/hww125.root", smurfFDir.Data()), outputFile, "hww125");
  singleboost(Form("%s/hww130.root", smurfFDir.Data()), outputFile, "hww130");
  singleboost(Form("%s/hww135.root", smurfFDir.Data()), outputFile, "hww135");
  singleboost(Form("%s/hww140.root", smurfFDir.Data()), outputFile, "hww140");
  singleboost(Form("%s/hww145.root", smurfFDir.Data()), outputFile, "hww145");
  singleboost(Form("%s/hww150.root", smurfFDir.Data()), outputFile, "hww150");
  singleboost(Form("%s/hww155.root", smurfFDir.Data()), outputFile, "hww155");
  singleboost(Form("%s/hww160.root", smurfFDir.Data()), outputFile, "hww160");
  singleboost(Form("%s/hww170.root", smurfFDir.Data()), outputFile, "hww170");
  singleboost(Form("%s/hww180.root", smurfFDir.Data()), outputFile, "hww180");
  singleboost(Form("%s/hww190.root", smurfFDir.Data()), outputFile, "hww190");
  singleboost(Form("%s/hww200.root", smurfFDir.Data()), outputFile, "hww200");
  singleboost(Form("%s/hww250.root", smurfFDir.Data()), outputFile, "hww250");
  singleboost(Form("%s/hww300.root", smurfFDir.Data()), outputFile, "hww300");
  singleboost(Form("%s/hww350.root", smurfFDir.Data()), outputFile, "hww350");
  singleboost(Form("%s/hww400.root", smurfFDir.Data()), outputFile, "hww400");
  singleboost(Form("%s/hww450.root", smurfFDir.Data()), outputFile, "hww450");
  singleboost(Form("%s/hww500.root", smurfFDir.Data()), outputFile, "hww500");
  singleboost(Form("%s/hww550.root", smurfFDir.Data()), outputFile, "hww550");
  singleboost(Form("%s/hww600.root", smurfFDir.Data()), outputFile, "hww600");
  singleboost(Form("%s/hww700.root", smurfFDir.Data()), outputFile, "hww700");
  singleboost(Form("%s/hww800.root", smurfFDir.Data()), outputFile, "hww800");
  singleboost(Form("%s/hww900.root", smurfFDir.Data()), outputFile, "hww900");
  singleboost(Form("%s/hww1000.root", smurfFDir.Data()), outputFile, "hww1000");
  
  outputFile->Close();
}

void singleboost(TString inputFileName, TFile *outputFile, TString histName) {

  TFile* fin = new TFile(inputFileName, "READ");
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  gROOT->cd();
  
  TH1F *hkx_0j = new TH1F(Form("%s_kx", histName.Data()), Form("%s_kx", histName.Data()), 25, -50, 50);
  TH1F *hky_0j = new TH1F(Form("%s_ky", histName.Data()), Form("%s_ky", histName.Data()), 25, -50, 50);
  TH1F *hkt_0j = new TH1F(Form("%s_kt", histName.Data()), Form("%s_kt", histName.Data()), 25, 0, 50);
  TH2F *hboost_0j = new TH2F(Form("%s_boost", histName.Data()), Form("%s_boost", histName.Data()), 25, -50, 50, 25, -50, 50);
  
  TH1F *hkx_1j = new TH1F(Form("%s_kx_1jet", histName.Data()), Form("%s_kx_1jet", histName.Data()), 25, -100, 100);
  TH1F *hky_1j = new TH1F(Form("%s_ky_1jet", histName.Data()), Form("%s_ky_1jet", histName.Data()), 25, -100, 100);
  TH1F *hkt_1j = new TH1F(Form("%s_kt_1jet", histName.Data()), Form("%s_kt_1jet", histName.Data()), 25, 0, 100);
  TH2F *hboost_1j = new TH2F(Form("%s_boost_1jet", histName.Data()), Form("%s_boost_1jet", histName.Data()), 25, -50, 50, 25, -50, 50);
  
  TH1F *hkx_2j = new TH1F(Form("%s_kx_2jet", histName.Data()), Form("%s_kx_2jet", histName.Data()), 25, -100, 100);
  TH1F *hky_2j = new TH1F(Form("%s_ky_2jet", histName.Data()), Form("%s_ky_2jet", histName.Data()), 25, -100, 100);
  TH1F *hkt_2j = new TH1F(Form("%s_kt_2jet", histName.Data()), Form("%s_kt_2jet", histName.Data()), 25, 0, 100);
  TH2F *hboost_2j = new TH2F(Form("%s_boost_2jet", histName.Data()), Form("%s_boost_2jet", histName.Data()), 25, -50, 50, 25, -50, 50);

  hkx_0j->Sumw2();
  hky_0j->Sumw2();
  hkt_0j->Sumw2();
  hboost_0j->Sumw2();

  hkx_1j->Sumw2();
  hky_1j->Sumw2();
  hkt_1j->Sumw2();
  hboost_1j->Sumw2();
  
  hkx_2j->Sumw2();
  hky_2j->Sumw2();
  hkt_2j->Sumw2();
  hboost_2j->Sumw2();
    
  // get event based branches..
  unsigned int njets_ = 0;
  unsigned int cuts_ = 0;
  LorentzVector*  dilep_ = 0;
  int type_ = 0; // 0/1/2/3 for mm/me/em/ee
  float pmet_ = 0.0;
  float pTrackMet_ = 0.0;
  LorentzVector*  jet1_ = 0;
  float dPhiDiLepJet1_ = 0.;
  float genmet_ = 0.;
  float genmetPhi_ = 0.;
  int lep1McId_ = 0;
  int lep2McId_ = 0;
  int lep1MotherMcId_ = 0;
  int lep2MotherMcId_ = 0;
  float scale1fb_ = 0.0;  
  
  ch->SetBranchAddress( "njets"         , &njets_         );     
  ch->SetBranchAddress( "cuts"          , &cuts_          );     
  ch->SetBranchAddress( "dilep"         , &dilep_         );   
  ch->SetBranchAddress( "type"          , &type_          );     
  ch->SetBranchAddress( "pmet"          , &pmet_          );     
  ch->SetBranchAddress( "pTrackMet"     , &pTrackMet_     );   
  ch->SetBranchAddress( "jet1"          , &jet1_          );   
  ch->SetBranchAddress( "dPhiDiLepJet1" , &dPhiDiLepJet1_ );     
  ch->SetBranchAddress( "genmet"        , &genmet_        );     
  ch->SetBranchAddress( "genmetPhi"     , &genmetPhi_     );     
  ch->SetBranchAddress( "lep1McId"      , &lep1McId_      );     
  ch->SetBranchAddress( "lep2McId"      , &lep2McId_      );     
  ch->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId_);     
  ch->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId_);     
  ch->SetBranchAddress( "scale1fb"      , &scale1fb_      );   
  
  
  //==========================================
  // Loop All Events
  //==========================================
  std::cout << inputFileName << " has " << ch->GetEntries() << " entries; \n";

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt); 

    // apply ww preselection 
    if ( dilep_->mass() < 12.0) continue;
    if ( dilep_->Pt() < 45) continue;    
    if (min(pmet_,pTrackMet_) < 20.)  continue;
    if ( njets_ > 3 ) continue;
    if ((cuts_ & ww) != ww) continue;
    if ( type_ == 0 || type_ == 3 ) {
      if ( min(pmet_, pTrackMet_) < 45. ) continue;
      if ( jet1_->Pt() > 15. && dPhiDiLepJet1_>(165.*TMath::Pi()/180.) ) continue;  
    } 
    // apply MC truth 
    if ( abs(lep1MotherMcId_) !=  24 || abs(lep2MotherMcId_) != 24 ) continue;
    if ( type_ == 0 ) {
      if ( abs(lep1McId_) != 13) continue;
      if ( abs(lep2McId_) != 13) continue;
    } else if ( type_ == 1 ) {
      if ( abs(lep1McId_) != 13) continue;
      if ( abs(lep1McId_) != 11) continue;
    } else if ( type_ == 2 ) {
      if ( abs(lep1McId_) != 11) continue;
      if ( abs(lep1McId_) != 13) continue;
    } else if ( type_ == 2 ) {
      if ( abs(lep1McId_) != 11) continue;
      if ( abs(lep1McId_) != 11) continue;
    }
    
    
    float kx(-999.), ky(-999.), kt(-999.); 
    kx = dilep_->Px() + genmet_ * cos(genmetPhi_);
    ky = dilep_->Py() + genmet_ * sin(genmetPhi_);
    kt = sqrt(kx*kx+ky*ky);
    
    if ( njets_ == 0 ) {
      hkx_0j->Fill(kx, scale1fb_);
      hky_0j->Fill(ky, scale1fb_);
      hkt_0j->Fill(kt, scale1fb_);
      hboost_0j->Fill(kx, ky, scale1fb_);
    } else if ( njets_ == 1) {
      hkx_1j->Fill(kx, scale1fb_);
      hky_1j->Fill(ky, scale1fb_);
      hkt_1j->Fill(kt, scale1fb_);
      hboost_1j->Fill(kx, ky, scale1fb_);
    } else {
      hkx_2j->Fill(kx, scale1fb_);
      hky_2j->Fill(ky, scale1fb_);
      hkt_2j->Fill(kt, scale1fb_);
      hboost_1j->Fill(kx, ky, scale1fb_);
    }
  }   //nevent
  
  
  outputFile->cd();
  hkx_0j->Write();
  hky_0j->Write();
  hkt_0j->Write();
  hboost_0j->Write();

  hkx_1j->Write();
  hky_1j->Write();
  hkt_1j->Write();
  hboost_1j->Write();

  hkx_2j->Write();
  hky_2j->Write();
  hkt_2j->Write();
  hboost_2j->Write();
}  
