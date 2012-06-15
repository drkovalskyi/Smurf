/*
   This script skims the common smurf ntuples to select 0-jet and ww preselections for the ME calculation  
   smurfproducer(smurfFDir, "ww.root", outputDir,cutstring)
   0 - location of the smurf ntuples
   1 - the sample to run over
   2 - the location of the skimmed smurf files
   3 - the cut string, choose only from those "WW", "ZZ" or "PassFail"
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

// steal from ../../Core/SmurfTree.h
enum Selection {
  BaseLine          = 1UL<<0,  // pt(reco)>20/10, acceptance,!STA muon, mll>12
  ChargeMatch       = 1UL<<1,  // q1*q2<0
  Lep1FullSelection = 1UL<<2,  // full id, isolation, d0, dz etc
  Lep1LooseEleV1    = 1UL<<3,  // electron fakeable object selection is passed V1
  Lep1LooseEleV2    = 1UL<<4,  // electron fakeable object selection is passed V2
  Lep1LooseEleV3    = 1UL<<5,  // electron fakeable object selection is passed V3
  Lep1LooseEleV4    = 1UL<<6,  // electron fakeable object selection is passed V4
  Lep1LooseMuV1     = 1UL<<7,  // muon fakeable object selection (relIso<1.0)
  Lep1LooseMuV2     = 1UL<<8,  // muon fakeable object selection (relIso<0.4)
  Lep2FullSelection = 1UL<<9,  // full id, isolation, d0, dz etc
  Lep2LooseEleV1    = 1UL<<10, // electron fakeable object selection is passed V1
  Lep2LooseEleV2    = 1UL<<11, // electron fakeable object selection is passed V2
  Lep2LooseEleV3    = 1UL<<12, // electron fakeable object selection is passed V3
  Lep2LooseEleV4    = 1UL<<13, // electron fakeable object selection is passed V4
  Lep2LooseMuV1     = 1UL<<14, // muon fakeable object selection (relIso<1.0)
  Lep2LooseMuV2     = 1UL<<15, // muon fakeable object selection (relIso<0.4)
  FullMET           = 1UL<<16, // full met selection
  ZVeto             = 1UL<<17, // event is not in the Z-mass peak for ee/mm final states
  TopTag            = 1UL<<18, // soft muon and b-jet tagging for the whole event regardless of n-jets (non-zero means    tagged)
  TopVeto           = 1UL<<19, // soft muon and b-jet tagging for the whole event regardless of n-jets (zero means tag ged)
  OneBJet           = 1UL<<20, // 1-jet events, where the jet is b-tagged (top control sample with one b-quark missing )
  TopTagNotInJets   = 1UL<<21, // soft muon and b-jet tagging for areas outside primary jets (non-zero means tagged)~
  ExtraLeptonVeto   = 1UL<<22, // extra lepton veto, DR(muon-electron)>=0.3
  Lep3FullSelection = 1UL<<23,  // full id, isolation, d0, dz etc
  Lep3LooseEleV1    = 1UL<<24, // electron fakeable object selection is passed V1
  Lep3LooseEleV2    = 1UL<<25, // electron fakeable object selection is passed V2
  Lep3LooseEleV3    = 1UL<<26, // electron fakeable object selection is passed V3
  Lep3LooseEleV4    = 1UL<<27, // electron fakeable object selection is passed V4
  Lep3LooseMuV1     = 1UL<<28, // muon fakeable object selection (relIso<1.0)
  Lep3LooseMuV2     = 1UL<<29  // muon fakeable object selection (relIso<0.4)
};

// ww-preselection information available through bits (extra cuts are needed)
UInt_t ww = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ZVeto|ExtraLeptonVeto;
// relax met
UInt_t ww_nomet = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ZVeto|ExtraLeptonVeto;
// relax the lepton requirements
UInt_t ww_lepfo = BaseLine|ChargeMatch|ZVeto|ExtraLeptonVeto;
// zz baseline selections using the bits setup in the smurfntuples
// http://www.t2.ucsd.edu/tastwiki/bin/view/Smurf/HZZllvvEventSelections#Reference_selection_V1
UInt_t zz_baseline = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ExtraLeptonVeto;

using namespace std;

//###################
//# main function
//###################
void smurfproducer(TString smurfFDir = "/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/", TString fileName = "zz.root", TString outputDir = "rawsmurfdata/", TString cutstring = "ZZ", int jetbin = 0 ) {

  TFile* fin = new TFile(smurfFDir+fileName);
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  
  TString outputFileName = outputDir + fileName;
  if (cutstring == "PassFail")   outputFileName.ReplaceAll(".root","_PassFail.root");
  if (cutstring == "ZZBTAG")   outputFileName.ReplaceAll(".root","_BTAG.root");
  
  TFile *newfile= new TFile(outputFileName,"recreate");
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // get event based branches..
  unsigned int njets_ = 0;
  unsigned int cuts_ = 0;
  LorentzVector*  lep1_ = 0;
  LorentzVector*  lep2_ = 0;
  LorentzVector*  dilep_ = 0;
  int type_ = 0; // 0/1/2/3 for mm/me/em/ee
  float pmet_ = 0.0;
  float pTrackMet_ = 0.0;
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  float mt_ = 0.;
  float dPhi_ = 0.;
  LorentzVector*  jet1_ = 0;
  float dPhiDiLepJet1_ = 0.;
  float met_ = 0.0;
  float trackMet_ = 0.0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  float jet1Btag_ = 0;
  float jet2Btag_ = 0;
  float jet3Btag_ = 0;
  unsigned int nSoftMuons_ = 0;
  unsigned int nvtx_ = 0;
  
  ch->SetBranchAddress( "njets"         , &njets_     );     
  ch->SetBranchAddress( "cuts"          , &cuts_     );     
  ch->SetBranchAddress( "lep1"          , &lep1_      );   
  ch->SetBranchAddress( "lep2"          , &lep2_      );   
  ch->SetBranchAddress( "dilep"         , &dilep_      );   
  ch->SetBranchAddress( "type"          , &type_     );     
  ch->SetBranchAddress( "pmet"          , &pmet_     );     
  ch->SetBranchAddress( "pTrackMet"     , &pTrackMet_     );   
  ch->SetBranchAddress( "run"           , &run_     );     
  ch->SetBranchAddress( "lumi"          , &lumi_     );     
  ch->SetBranchAddress( "dPhi"          , &dPhi_     );     
  ch->SetBranchAddress( "mt"            , &mt_     );     
  ch->SetBranchAddress( "jet1"          , &jet1_      );   
  ch->SetBranchAddress( "jet2"          , &jet2_      );   
  ch->SetBranchAddress( "jet3"          , &jet3_      );   
  ch->SetBranchAddress( "dPhiDiLepJet1" , &dPhiDiLepJet1_     );     
  ch->SetBranchAddress( "met"           , &met_     );     
  ch->SetBranchAddress( "trackMet"      , &trackMet_     );   
  ch->SetBranchAddress( "jet1Btag"      , &jet1Btag_      );   
  ch->SetBranchAddress( "jet2Btag"      , &jet2Btag_      );   
  ch->SetBranchAddress( "jet3Btag"      , &jet3Btag_      );   
  ch->SetBranchAddress( "nSoftMuons"    , &nSoftMuons_      );   
  ch->SetBranchAddress( "nvtx"          , &nvtx_      );   

  float scale1fb = 0.0;
  ch->SetBranchAddress( "scale1fb"      , &scale1fb     );   
  
  
  //==========================================
  // Loop All Events
  //==========================================
  
  cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; \n";

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt); 

    // select a specific jetbin
    // for 2-jet bin it combines >= 2 jet events

    if ( jetbin < 2 && int(njets_) != jetbin ) continue;
    if ( jetbin == 2 && int(njets_) < 2 ) continue; 
    if ( dilep_->mass() < 12.0) continue;
    if ( dilep_->Pt() < 45) continue;    
    if (TMath::Min(pmet_,pTrackMet_) < 20.)  continue;
    if ( ! (cuts_ & TopVeto) ) continue;
    
    // cuts to select the WW pre-selection
    if (cutstring == "WW") {
      if ((cuts_ & ww) != ww) continue;
    }
    
    // select the PassFail sample for the wjets studies
    else if (cutstring == "PassFail") {  
      if ( (cuts_ & ww_lepfo) != ww_lepfo) continue;
      // skip events with no lepton pass the full selection
      if ( ( (cuts_ & Lep1FullSelection) != Lep1FullSelection)  &&  
	   ( (cuts_ & Lep2FullSelection) != Lep2FullSelection) ) continue;
      
      // skip events with both leptons hat pass the final selection
      if ( ((cuts_ & Lep1FullSelection) == Lep1FullSelection)
	   && ((cuts_ & Lep2FullSelection) == Lep2FullSelection)) continue;
      
      // if lep1 pass full selection, but none of the lep2 pass the FO definition, skip the event
      if ( cuts_ & Lep1FullSelection ) {
	if ( ! ( (cuts_ & Lep2LooseEleV4) || (cuts_ & Lep2LooseMuV2) ) ) continue;
      }
      // if lep2 pass full selection, but none of the lep1 pass the FO definition, skip the event
      if ( cuts_ & Lep2FullSelection ) {
	if ( ! ( (cuts_ & Lep1LooseEleV4) || (cuts_ & Lep1LooseMuV2) ) ) continue;
      }
    } 
    // cuts to select the ZZ pre-selection
    else if (cutstring == "ZZ") {
      if ( TMath::Abs(dilep_->M()-91.1876) > 15.0 ) continue;
      if ( (cuts_ & zz_baseline) != zz_baseline) continue;
      if (lep1_->Pt() < 20.) continue;
      if (lep2_->Pt() < 20.) continue;
      if ( met_ < 50 ) continue;
      if ( dilep_->Pt() < 55.0) continue;
      if (jet1Btag_ > 2.0 && jet1_->Pt() > 30) continue;
      if (jet2Btag_ > 2.0 && jet2_->Pt() > 30) continue;
      if (jet3Btag_ > 2.0 && jet3_->Pt() > 30) continue;
      if( nSoftMuons_ != 0 ) continue;
    }

    else if (cutstring == "ZZBTAG" ) {
      if ( (cuts_ & zz_baseline) != zz_baseline) continue;
      if ( TMath::Abs(dilep_->M()-91.1876) > 15.0 ) continue;
      if (lep1_->Pt() < 20.) continue;
      if (lep2_->Pt() < 20.) continue;
      if ( met_ < 70 ) continue;
      if ( dilep_->Pt() < 55.0) continue;
      bool btag = false;
      if ( jet1Btag_ > 2.0 && jet1_->Pt() > 30) btag = true;
      if ( jet2Btag_ > 2.0 && jet2_->Pt() > 30) btag = true;
      if ( jet3Btag_ > 2.0 && jet3_->Pt() > 30) btag = true;
      if ( nSoftMuons_ != 0 ) btag = true;
      if ( !btag ) continue;
    }
    
    evt_tree->Fill();
  }   //nevent
  
  cout << outputDir + fileName << " has " << evt_tree->GetEntries() << " entries; \n";
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
}  
