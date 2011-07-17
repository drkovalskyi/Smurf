//root -l Smurf/Analysis/HZZllvv/MakePhotonReweightFactors.C+\(\"/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data_photons.root\",\"/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data-met0-1092ifb-dileptrig.root\"\)

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include "TRandom.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <iomanip>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "TH2D.h"
#include "/home/sixie/CMSSW_analysis/src/Smurf/Core/SmurfTree.h"
#include "/home/sixie/CMSSW_analysis/src/Smurf/HWW/factors.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void MakePhotonReweightFactors
(
  TString PhotonFile    = "/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data_photons.root",
  TString DileptonFile    = "/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data.root"
  )
{

  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("Smurf/Analysis/HZZllvv/goodlumiList.1092ipb.json"); 

  TChain *chPhoton = new TChain("HwwTree0");
  chPhoton->Add(PhotonFile);

  TChain *chDilepton = new TChain("tree");
  chDilepton->Add(DileptonFile);

  TTree *treePhoton  = (TTree*) chPhoton;
  TTree *treeDilepton  = (TTree*) chDilepton;


  Float_t jetbinsarr[5] = {-0.5, 0.5, 1.5, 2.5, 3.5};
  Float_t ptbinsarr[18] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 
                           120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0};
  TH1F *DileptonMass = new TH1F("h1_mass", "h1_mass", 100, 0, 200);
  TH2F *Pt_Photons   = new TH2F("h2_g", "g", 17, ptbinsarr, 4, jetbinsarr);
  TH2F *Pt_Dileptons = new TH2F("h2_z", "z", 17, ptbinsarr, 4, jetbinsarr);
  Pt_Photons->Sumw2();
  Pt_Dileptons->Sumw2();


  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          njets;
  UInt_t          event;
  UInt_t          run;
  UInt_t          lumi;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  jet1  = 0;
  LorentzVector*  jet2  = 0;
  Float_t         dPhi;
  Float_t         dR;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         met;
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
  Float_t         bdt = 0.0;
  Float_t         bdtd = 0.0;
  Float_t         nn = 0.0;
  Float_t         knn = 0.0;
  Float_t         bdtg = 0.0;
  Float_t         higgsPt = -999;



  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------
  //
  //LOOP OVER Photon Data Sample
  //
  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------


  treePhoton->SetBranchAddress( "cuts"         , &cuts         );
  treePhoton->SetBranchAddress( "dstype"       , &dstype       );
  treePhoton->SetBranchAddress( "nvtx"         , &nvtx         );
  treePhoton->SetBranchAddress( "njets"        , &njets        );
  treePhoton->SetBranchAddress( "run"          , &run          );
  treePhoton->SetBranchAddress( "lumi"         , &lumi         );
  treePhoton->SetBranchAddress( "event"        , &event        );
  treePhoton->SetBranchAddress( "scale1fb"     , &scale1fb     );
  treePhoton->SetBranchAddress( "lep1"         , &lep1         );
  treePhoton->SetBranchAddress( "lep2"         , &lep2         );
  treePhoton->SetBranchAddress( "jet1"         , &jet1         );
  treePhoton->SetBranchAddress( "jet2"         , &jet2         );
  treePhoton->SetBranchAddress( "dPhi"         , &dPhi         );
  treePhoton->SetBranchAddress( "dR"           , &dR           );
  treePhoton->SetBranchAddress( "dilep"        , &dilep        );
  treePhoton->SetBranchAddress( "type"         , &type         );
  treePhoton->SetBranchAddress( "pmet"         , &pmet         );
  treePhoton->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  treePhoton->SetBranchAddress( "met"          , &met          );
  treePhoton->SetBranchAddress( "trackMet"     , &trackMet     );
  treePhoton->SetBranchAddress( "mt"           , &mt           );
  treePhoton->SetBranchAddress( "mt1"          , &mt1          );
  treePhoton->SetBranchAddress( "mt2"          , &mt2          );
  treePhoton->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  treePhoton->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  treePhoton->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  treePhoton->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  treePhoton->SetBranchAddress( "lq1"          , &lq1          );
  treePhoton->SetBranchAddress( "lq2"          , &lq2          );
  treePhoton->SetBranchAddress( "lid1"         , &lid1         );
  treePhoton->SetBranchAddress( "lid2"         , &lid2         );
  treePhoton->SetBranchAddress( "lid3"         , &lid3         );
  treePhoton->SetBranchAddress( "processId"    , &processId    );
  treePhoton->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  treePhoton->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  treePhoton->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  treePhoton->SetBranchAddress( "jet2Btag"     , &jet2Btag     );

  Int_t NPhotonEvents = 0;
  for (UInt_t i=0; i<treePhoton->GetEntries(); i++) {
    
    treePhoton->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)treePhoton->GetEntries());

    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

    //****************************************************
    //Photon ID
    //****************************************************
    Bool_t passSpikeKillingEMaxOverE9 = ( (UInt_t(lid1) & 1) == 1);
    Bool_t passSpikeKillingSigmaiEtaiEta = ((UInt_t(lid1) & 2) == 2);
    Bool_t passSpikeKillingSigmaiPhiiPhi = ((UInt_t(lid1) & 4) == 4);
    Bool_t passSpikeKillingSwissCross = ((UInt_t(lid1) & 8) == 8);
    Bool_t passR9Tight = ((UInt_t(lid1) & 16) == 16);
    Bool_t isBarrel = ((UInt_t(lid1) & 32) == 32);
    Bool_t isConverted = ((UInt_t(lid1) & 64) == 64);

    //Apply no cuts for now   
    Bool_t passPhotonCuts = kTRUE;    
    if (!passPhotonCuts) continue;

    //use only barrel photons
    if (fabs(dilep->eta()) > 1.5) continue;
    
    //************************************************************************************************
    // PreSelection + Anti Met Cut
    //************************************************************************************************
    if( njets < 2 && jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( !(TMath::Min(met,trackMet) < 50.0)) continue;

 
    //************************************************************************************************
    // Fill njet, pt
    //************************************************************************************************
    NPhotonEvents++;
    Pt_Photons->Fill(TMath::Min(dilep->pt(), double(299.0)),TMath::Min(Int_t(njets), 3));

  } //End Signal Loop

  
  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------
  //
  //LOOP OVER Dilepton Data Sample
  //
  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------


  treeDilepton->SetBranchAddress( "cuts"         , &cuts         );
  treeDilepton->SetBranchAddress( "dstype"       , &dstype       );
  treeDilepton->SetBranchAddress( "nvtx"         , &nvtx         );
  treeDilepton->SetBranchAddress( "njets"        , &njets        );
  treeDilepton->SetBranchAddress( "run"          , &run          );
  treeDilepton->SetBranchAddress( "lumi"         , &lumi         );
  treeDilepton->SetBranchAddress( "event"        , &event        );
  treeDilepton->SetBranchAddress( "scale1fb"     , &scale1fb     );
  treeDilepton->SetBranchAddress( "lep1"         , &lep1         );
  treeDilepton->SetBranchAddress( "lep2"         , &lep2         );
  treeDilepton->SetBranchAddress( "jet1"         , &jet1         );
  treeDilepton->SetBranchAddress( "jet2"         , &jet2         );
  treeDilepton->SetBranchAddress( "dPhi"         , &dPhi         );
  treeDilepton->SetBranchAddress( "dR"           , &dR           );
  treeDilepton->SetBranchAddress( "dilep"        , &dilep        );
  treeDilepton->SetBranchAddress( "type"         , &type         );
  treeDilepton->SetBranchAddress( "pmet"         , &pmet         );
  treeDilepton->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  treeDilepton->SetBranchAddress( "met"          , &met          );
  treeDilepton->SetBranchAddress( "trackMet"     , &trackMet     );
  treeDilepton->SetBranchAddress( "mt"           , &mt           );
  treeDilepton->SetBranchAddress( "mt1"          , &mt1          );
  treeDilepton->SetBranchAddress( "mt2"          , &mt2          );
  treeDilepton->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  treeDilepton->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  treeDilepton->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  treeDilepton->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  treeDilepton->SetBranchAddress( "lq1"          , &lq1          );
  treeDilepton->SetBranchAddress( "lq2"          , &lq2          );
  treeDilepton->SetBranchAddress( "lid1"         , &lid1         );
  treeDilepton->SetBranchAddress( "lid2"         , &lid2         );
  treeDilepton->SetBranchAddress( "lid3"         , &lid3         );
  treeDilepton->SetBranchAddress( "processId"    , &processId    );
  treeDilepton->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  treeDilepton->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  treeDilepton->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  treeDilepton->SetBranchAddress( "jet2Btag"     , &jet2Btag     );
  treeDilepton->SetBranchAddress( "lep1McId"     , &lep1McId     );
  treeDilepton->SetBranchAddress( "lep2McId"     , &lep2McId     );

  Int_t NZEvents = 0;
  for (UInt_t i=0; i<treeDilepton->GetEntries(); i++) {
    
    treeDilepton->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)treeDilepton->GetEntries());
    
    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

    //************************************************************************************************
    // PreSelection
    //************************************************************************************************
    if( !(lep1->pt() > 20 && lep2->pt() > 20 )) continue;
    if( !((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection )) continue;
    if( !((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection )) continue;
    if( !(type == SmurfTree::mm || type == SmurfTree::ee )) continue;
    if( !(fabs(dilep->mass() - 91.1876) < 15 )) continue;
    if( !(lq1*lq2 < 0)) continue;
    if( !((cuts & SmurfTree::TopVeto) == SmurfTree::TopVeto )) continue;
    if( !((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto )) continue;
    if( !(dilep->pt() > 40)) continue;
    if( njets < 2 && jet1->pt() > 15 && (dPhiDiLepJet1 > 165.0 * TMath::Pi() / 180.0) ) continue; 

    //Fill Mass Distribution before anti Met cut
    DileptonMass->Fill(dilep->mass());

    //Anti Met Cut
    if( !(TMath::Min(met,trackMet) < 50.0)) continue;
    
    //************************************************************************************************
    // Fill njet, pt
    //************************************************************************************************
    NZEvents++;
    Pt_Dileptons->Fill(TMath::Min(dilep->pt(), 299.0),TMath::Min(Int_t(njets), 3));
  
  } //End Signal Loop


  cout << "NPhotonEvents:" <<   NPhotonEvents << endl;
  cout << "NZEvents: " << NZEvents << endl;
  //************************************************************************************************
  // Make Reweight Factor
  //************************************************************************************************

  TH2D *ReweightFactorPtNJet = (TH2D*)Pt_Photons->Clone("h2_weight");
  ReweightFactorPtNJet->SetBinContent(0,0,1.0);
  for(UInt_t a=1; a < ReweightFactorPtNJet->GetXaxis()->GetNbins()+2; ++a) {
    ReweightFactorPtNJet->SetBinContent(a,0,1.0);
    for(UInt_t b=1; b < ReweightFactorPtNJet->GetYaxis()->GetNbins()+2; ++b) {
      if (Pt_Photons->GetBinContent(a,b)>0) {
        Double_t ratio =  Pt_Dileptons->GetBinContent(a,b) / Pt_Photons->GetBinContent(a,b);
        Double_t ratioUncertainty = 0.0;
        if (Pt_Dileptons->GetBinContent(a,b) > 0 && Pt_Photons->GetBinContent(a,b) > 0) {
          ratioUncertainty = ratio * TMath::Sqrt( 1.0 / Pt_Dileptons->GetBinContent(a,b) + 
                                                  1.0 / Pt_Photons->GetBinContent(a,b));
        }
        ReweightFactorPtNJet->SetBinContent(a,b, ratio);
        ReweightFactorPtNJet->SetBinError(a,b, ratioUncertainty);
      } else {
        ReweightFactorPtNJet->SetBinContent(a,b, 0.0);
        ReweightFactorPtNJet->SetBinError(a,b, 0.0);
      }
    }
  }

  //************************************************************************************************
  // Plots
  //************************************************************************************************
  TFile *file = 0;

  file = new TFile("PhotonJetsReweightFactors.root", "UPDATE");
  file->WriteTObject(ReweightFactorPtNJet, ReweightFactorPtNJet->GetName(), "WriteDelete");

  file->WriteTObject(DileptonMass, DileptonMass->GetName(), "WriteDelete");

  file->Close();




  return;
}

