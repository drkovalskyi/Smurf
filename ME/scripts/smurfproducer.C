/*
   This script skims the common smurf ntuples to select 0-jet and ww preselections for the ME calculation  
   smurfproducer(smurfFDir, "ww.root", outputDir)
   0 - location of the smurf ntuples
   1 - the sample to run over
   2 - the location of the skimmed smurf files
 */


#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include <utility>

// Steal the code from SmurfTree.h
// 
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
  TopTag            = 1UL<<18, // soft muon and b-jet tagging for the whole event regardless of n-jets (non-zero means tagged)
  TopVeto           = 1UL<<19, // soft muon and b-jet tagging for the whole event regardless of n-jets (zero means tagged)
  OneBJet           = 1UL<<20, // 1-jet events, where the jet is b-tagged (top control sample with one b-quark missing)
  TopTagNotInJets   = 1UL<<21, // soft muon and b-jet tagging for areas outside primary jets (non-zero means tagged)
  ExtraLeptonVeto   = 1UL<<22, // extra lepton veto, DR(muon-electron)>=0.3
  Lep3FullSelection = 1UL<<23,  // full id, isolation, d0, dz etc
  Lep3LooseEleV1    = 1UL<<24, // electron fakeable object selection is passed V1
  Lep3LooseEleV2    = 1UL<<25, // electron fakeable object selection is passed V2
  Lep3LooseEleV3    = 1UL<<26, // electron fakeable object selection is passed V3
  Lep3LooseEleV4    = 1UL<<27, // electron fakeable object selection is passed V4
  Lep3LooseMuV1     = 1UL<<28, // muon fakeable object selection (relIso<1.0)
  Lep3LooseMuV2     = 1UL<<29  // muon fakeable object selection (relIso<0.4)
};


UInt_t ww = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ZVeto|TopVeto|ExtraLeptonVeto;

UInt_t ww_nolepton = BaseLine|ChargeMatch|FullMET|ZVeto|TopVeto|ExtraLeptonVeto;

using namespace std;

//###################
//# main function
//###################
void smurfproducer(TString smurfFDir = "/smurf/data/Run2011_Spring11_SmurfV3/mitf-alljets/", TString fileName = "zz.root", TString outputDir = "/smurf/data/Run2011_Spring11_SmurfV3/mitf-zerojets-wwselection/", TString cutstring = "WW" ) {

  TFile* fin = new TFile(smurfFDir+fileName);
  TFile *newfile = new TFile(outputDir+fileName,"recreate");

  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // Initialize the branches to use to calculate the differential cross-sections
  unsigned int njets_ = 0;
  unsigned int cuts_ = 0;
  int type_ = 0;
  float pmet_ = 0.0;
  float pTrackMet_ = 0.0;


  ch->SetBranchAddress( "njets"      , &njets_     );     
  ch->SetBranchAddress( "cuts"      , &cuts_     );     
  ch->SetBranchAddress( "type"      , &type_     );     
  ch->SetBranchAddress( "pmet"      , &pmet_     );     
  ch->SetBranchAddress( "pTrackMet"      , &pTrackMet_     );     

  //==========================================
  // Loop All Events
  //==========================================
  
  cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; \n";

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt); 
    
    // select only 0 jet final states
    if (njets_ > 0) continue;
    if ( (cuts_ & BaseLine) != BaseLine ) continue;
    if ( (cuts_ & ChargeMatch) != ChargeMatch ) continue;
    if ( (cuts_ & FullMET) != FullMET ) continue;
    if ( (type_ == 0 || type_ == 3) && TMath::Min(pmet_, pTrackMet_) < 35 ) continue;
    if ( (type_ == 1 || type_ == 2) && TMath::Min(pmet_, pTrackMet_) < 20 ) continue;
    if ( (cuts_ & ZVeto) != ZVeto ) continue;
    if ( (cuts_ & TopVeto) != TopVeto ) continue;
    if ( (cuts_ & ExtraLeptonVeto) != ExtraLeptonVeto ) continue;

    /*
    if (cutstring == "WW") {
      if ( (cuts_ & Lep1FullSelection) != Lep1FullSelection) continue;
      if ( (cuts_ & Lep2FullSelection) != Lep2FullSelection) continue;
    }

    if (cutstring == "PassFail") {
      
      // veto the ones that pass the final selection
      if ( ((cuts_ & Lep1FullSelection) == Lep1FullSelection)
	   && ((cuts_ & Lep2FullSelection) == Lep2FullSelection)) 	continue;
      
      // if neither lepton pass the full selection, skip the event
      if ( ( (cuts_ & Lep1FullSelection) != Lep1FullSelection)  &&  
      	   ( (cuts_ & Lep2FullSelection) != Lep2FullSelection) ) continue;
      
   
      // if lep1 pass full selection, but none of the lep2 pass the FO definition, skip the event
      if (     ( (cuts_ & Lep1FullSelection) == Lep1FullSelection) 
	    && ( (cuts_ & Lep2LooseEleV4) != Lep2LooseEleV4)
	       && ( (cuts_ & Lep2LooseMuV1) != Lep2LooseMuV1)) continue;
      
      // if lep1 pass full selection, but none of the lep2 pass the FO definition, skip the event
      if (     ( (cuts_ & Lep2FullSelection) == Lep2FullSelection) 
	       && ( (cuts_ & Lep1LooseEleV4) != Lep1LooseEleV4)
	       && ( (cuts_ & Lep1LooseMuV1) != Lep1LooseMuV1) ) continue;
    }
    */
    
    evt_tree->Fill();
  }   //nevent
  
  cout << outputDir + fileName << " has " << evt_tree->GetEntries() << " entries; \n";
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
}  
