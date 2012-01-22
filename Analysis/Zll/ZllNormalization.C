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
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


int    verboseLevel =   0;
const double sigmaB = 0.35;

//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void ZllNormalization
(
 TString datInputFile    = "/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data.root",
 TString bgdInputFile    = "/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/bkg.root",
 string PUReweightFile   = "/data/smurf/data/Winter11/auxiliar/PileupReweighting.Summer11DYmm_To_Run2011A.root",
 string effFile          = "/data/smurf/data/Winter11/auxiliar/efficiency_results_v7_42x_Run2011A.root",
 string jsonFile         = "/data/smurf/sixie/data/auxiliar/2011a.json",
 Double_t scaleFactorLum = 1.0
 )
{

  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
  chdata->Add(datFile1);

  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;


  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 

  //************************************************************************************************
  // PU Reweighting
  //************************************************************************************************
  TFile *fPUWeightsFile = TFile::Open(PUReweightFile.c_str());
  TH1D *fPUWeights = (TH1D*)(fPUWeightsFile->Get("puWeights"));
  fPUWeights->SetDirectory(0);
  fPUWeightsFile->Close();
  delete fPUWeightsFile;
  TH1D *myPUWeights = fPUWeights;

  //************************************************************************************************
  // Lepton Efficiencies
  //************************************************************************************************
  LeptonScaleLookup leptonEfficiencies(effFile.c_str());


  //----------------------------------------------------------------------------



  double NZmmSelected_Data = 0;
  double NZmmSelected_MC = 0;
  double NZmmSelected_MCUnCorr = 0;
  double NZeeSelected_Data = 0;
  double NZeeSelected_MC = 0;
  double NZeeSelected_MCUnCorr = 0;

  TH1D *ZmmMass_Data = new TH1D("ZmmMass_Data", ";Mass [GeV/c^{2}]; NEvents; ", 200, 0, 200);
  TH1D *ZmmMass_MC = new TH1D("ZmmMass_MC", ";Mass [GeV/c^{2}]; NEvents; ", 200, 0, 200);
  TH1D *ZmmMass_MCUnCorr = new TH1D("ZmmMass_MCUnCorr", ";Mass [GeV/c^{2}]; NEvents; ", 200, 0, 200);
  TH1D *ZeeMass_Data = new TH1D("ZeeMass_Data", ";Mass [GeV/c^{2}]; NEvents; ", 200, 0, 200);
  TH1D *ZeeMass_MC = new TH1D("ZeeMass_MC", ";Mass [GeV/c^{2}]; NEvents; ", 200, 0, 200);
  TH1D *ZeeMass_MCUnCorr = new TH1D("ZeeMass_MCUnCorr", ";Mass [GeV/c^{2}]; NEvents; ", 200, 0, 200);

  double N0Jet_Data = 0;
  double N1Jet_Data = 0;
  double N2Jet_Data = 0;
  double N0Jet_MC = 0;
  double N1Jet_MC = 0;
  double N2Jet_MC = 0;

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
  UInt_t          npu;



  //*****************************************************************************
  // MC Loop
  //*****************************************************************************

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
  background->SetBranchAddress( "npu"           , &npu            );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%10000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    if ( lep1->pt() > 20 && lep2->pt() > 20 && dilep->mass() > 60 && dilep->mass() < 120 
         && ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
         && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
      ) {

      
      //****************************************************************************************
      //Event Weights
      //****************************************************************************************
      Double_t weightFactor = 1.0;

      //PU Reweighting
      double PUWeight = myPUWeights->GetBinContent(myPUWeights->GetXaxis()->FindFixBin(npu));

      //Lepton Efficiency Scale Factor 
      double LeptonEffScaleFactor = leptonEfficiencies.GetExpectedLeptonSF(lep1->eta(),lep1->pt(),  lid1)*leptonEfficiencies.GetExpectedLeptonSF( lep2->eta(),lep2->pt(), lid2);

      //Trigger Efficiency
      double trigEff = leptonEfficiencies.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , fabs(lep2->eta()), lep2->pt(), TMath::Abs(lid1), TMath::Abs(lid2));

      weightFactor = PUWeight * LeptonEffScaleFactor * trigEff;

      Double_t eventWeight = scaleFactorLum*scale1fb*weightFactor;

      if ( type == SmurfTree::mm ) {        
        NZmmSelected_MC += eventWeight;
        NZmmSelected_MCUnCorr += scale1fb*scaleFactorLum;
        ZmmMass_MC->Fill(dilep->mass() , eventWeight);
      }
      
      if ( type == SmurfTree::ee ) {        
        NZeeSelected_MC += eventWeight;
        NZeeSelected_MCUnCorr += scale1fb*scaleFactorLum;
        ZeeMass_MC->Fill(dilep->mass() , eventWeight);
      }

      //NJet Distribution
      if (njets==0) {
        N0Jet_MC += eventWeight;
      } else if (njets ==1) {
        N1Jet_MC += eventWeight;
      } else {
        N2Jet_MC += eventWeight;
      }

    }
  }
  printf("--- Finished Bgdnal loop\n");




  //*****************************************************************************
  // Data Loop
  //*****************************************************************************

  data->SetBranchAddress( "cuts"         , &cuts         );
  data->SetBranchAddress( "nvtx"         , &nvtx         );
  data->SetBranchAddress( "njets"        , &njets        );
  data->SetBranchAddress( "event"        , &event        );
  data->SetBranchAddress( "run"          , &run          );
  data->SetBranchAddress( "lumi"         , &lumi         );
  data->SetBranchAddress( "scale1fb"     , &scale1fb     );
  data->SetBranchAddress( "lep1"         , &lep1         );
  data->SetBranchAddress( "lep2"         , &lep2         );
  data->SetBranchAddress( "jet1"         , &jet1         );
  data->SetBranchAddress( "jet2"         , &jet2         );
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

  NZmmSelected_Data = 0;
  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    //****************************************************************************************
    //Good Lumi Selection
    //****************************************************************************************
    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(!rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
        
    if ( lep1->pt() > 20 && lep2->pt() > 20 && dilep->mass() > 60 && dilep->mass() < 120 
         && ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
         && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
     ) {

      if (type == SmurfTree::mm) {
        NZmmSelected_Data++;
        ZmmMass_Data->Fill(dilep->mass());
      }

      if (type == SmurfTree::ee) {
        NZeeSelected_Data++;
        ZeeMass_Data->Fill(dilep->mass());
      }

      //NJet Distribution
      if (njets==0) {
        N0Jet_Data += 1.0;
      } else if (njets ==1) {
        N1Jet_Data += 1.0;
      } else {
        N2Jet_Data += 1.0;
      }

    }

  }


  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("zyield.root", "UPDATE");

  TCanvas *cv = new TCanvas("cv", "cv", 800, 600);
  TLegend *legend = new TLegend(0.5,0.75,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);



  cout << "Zmm\n";
  cout << "Data: " << NZmmSelected_Data << endl;
  cout << "MC : " << NZmmSelected_MC << endl;
  cout << "MCUnCorr : " << NZmmSelected_MCUnCorr << endl;

  legend->Clear();
  legend->AddEntry(ZmmMass_Data, "Data", "L");
  legend->AddEntry(ZmmMass_MC, "MC @ 187 pb^{-1}", "L");
  ZmmMass_Data->SetLineColor(kBlue);
  ZmmMass_Data->SetLineWidth(2);
  ZmmMass_Data->SetMaximum(ZmmMass_Data->GetMaximum() * 1.2);
  ZmmMass_Data->Draw();
  ZmmMass_MC->SetLineColor(kRed);
  ZmmMass_MC->SetLineWidth(2);
  ZmmMass_MC->Draw("same");
  legend->Draw();
  cv->SaveAs("zmm.gif");
  f->WriteTObject(cv, "zmm", "WriteDelete");
 
  cout << "Zee\n";
  cout << "Data: " << NZeeSelected_Data << endl;
  cout << "MC : " << NZeeSelected_MC << endl;
  cout << "MCUnCorr : " << NZeeSelected_MCUnCorr << endl;

  legend->Clear();
  legend->AddEntry(ZeeMass_Data, "Data", "L");
  legend->AddEntry(ZeeMass_MC, "MC @ 187 pb^{-1}", "L");
  ZeeMass_Data->SetLineColor(kBlue);
  ZeeMass_Data->SetLineWidth(2);
  ZeeMass_Data->SetMaximum(ZeeMass_Data->GetMaximum() * 1.2);
  ZeeMass_Data->Draw();
  ZeeMass_MC->SetLineColor(kRed);
  ZeeMass_MC->SetLineWidth(2);
  ZeeMass_MC->Draw("same");
  legend->Draw(); 
  cv->SaveAs("zee.gif");
  f->WriteTObject(cv, "zee", "WriteDelete");
 
  f->Close();


  cout << "Data: " << N0Jet_Data << " " << N1Jet_Data << " " << N2Jet_Data << endl;
  cout << "MC: " << N0Jet_MC << " " << N1Jet_MC << " " << N2Jet_MC << endl;

}
