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
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"

//root -l -q -b Smurf/Analysis/Zll/ZllGSel.C+'("Ana/nt_scripts/ntuples_zh_53x/lgamma90.root",2012)'
//root -l -q -b Smurf/Analysis/Zll/ZllGSel.C+'("Ana/nt_scripts/ntuples_zh_42x/lgamma90.root",2011)'

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void ZllGSel
(
 TString inputMCFileName  = "inputMC.root",
 int year = 2012
)
{

  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  double lumi = 1.0; double prescale = 1.0;
  if	 (year == 2012){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365; prescale = 0.02;
  }
  else if(year == 2011){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/puWeights_Fall11_42x_True.root";
    lumi     =  4.924; prescale = 0.02;
  }
  else {assert(0);}

  TChain *chinputMC = new TChain("tree");
  chinputMC->Add(inputMCFileName);
  TTree *inputMC = (TTree*) chinputMC;

  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  LeptonScaleLookup trigLookup(effPath.Data());

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  TH1D* hDM3l  = new TH1D("hDM3l", "hDM3l", 120, 0, 120); hDM3l ->Sumw2();
  TH1D* hBM3l  = new TH1D("hBM3l", "hBM3l", 120, 0, 120); hBM3l ->Sumw2();
  TH1D* hDM2l  = new TH1D("hDM2l", "hDM2l", 120, 0, 120); hDM2l ->Sumw2();
  TH1D* hBM2l  = new TH1D("hBM2l", "hBM2l", 120, 0, 120); hBM2l ->Sumw2();
  TH1D* hDPtPh = new TH1D("hDPtPh","hDPtPh",100, 0, 100); hDPtPh->Sumw2();
  TH1D* hBPtPh = new TH1D("hBPtPh","hBPtPh",100, 0, 100); hBPtPh->Sumw2();
  TH1D* hDPt1  = new TH1D("hDPt1" ,"hDPt1", 100, 0, 100); hDPt1 ->Sumw2();
  TH1D* hBPt1  = new TH1D("hBPt1" ,"hBPt1", 100, 0, 100); hBPt1 ->Sumw2();
  TH1D* hDPt2  = new TH1D("hDPt2" ,"hDPt2", 100, 0, 100); hDPt2 ->Sumw2();
  TH1D* hBPt2  = new TH1D("hBPt2" ,"hBPt2", 100, 0, 100); hBPt2 ->Sumw2();

  UInt_t  dstype;
  UInt_t  cuts;
  Int_t   lid1;
  Int_t   lid2;
  Int_t   lid3;
  Int_t   lq1;
  Int_t   lq2;
  Int_t   lq3;
  Float_t npu;
  float scale1fb;
  LorentzVector* dilep = 0;
  LorentzVector* lep1 = 0;
  LorentzVector* lep2 = 0;
  LorentzVector* lep3 = 0;

  inputMC->SetBranchAddress( "dstype"	   , &dstype	  );
  inputMC->SetBranchAddress( "npu"	   , &npu	  );
  inputMC->SetBranchAddress( "lid1"	   , &lid1	  );
  inputMC->SetBranchAddress( "lid2"	   , &lid2	  );
  inputMC->SetBranchAddress( "lid3"	   , &lid3	  );
  inputMC->SetBranchAddress( "scale1fb"	   , &scale1fb	  );
  inputMC->SetBranchAddress( "cuts"        , &cuts        );
  inputMC->SetBranchAddress( "dilep"       , &dilep       );
  inputMC->SetBranchAddress( "lep1"        , &lep1        );
  inputMC->SetBranchAddress( "lep2"        , &lep2        );
  inputMC->SetBranchAddress( "lep3"        , &lep3        );
  inputMC->SetBranchAddress( "lq1"         , &lq1         );
  inputMC->SetBranchAddress( "lq2"         , &lq2         );
  inputMC->SetBranchAddress( "lq3"         , &lq3         );
  for (UInt_t i=0; i<inputMC->GetEntries(); i++) {    
    inputMC->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputMC->GetEntries());

    int FullLid  = (cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection;
	FullLid += (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
	FullLid += (cuts & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection;
    int Photonid  = (cuts & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2;
	Photonid += (cuts & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2;
	Photonid += (cuts & SmurfTree::Lep3LooseEleV2) == SmurfTree::Lep3LooseEleV2;

    if(FullLid != 2 || Photonid != 1) continue;

    if(lep1->pt() <= 10) continue;
    if(lep2->pt() <= 10) continue;
    if(lep3->pt() <= 10) continue;

    bool goodFlavor = false;
    double dilMass  = 0.0; double PtPh = 0.0; double Pt1 = 0.0;  double Pt2 = 0.0; 
    if     ((cuts & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2){
      if(lq2*lq3 < 0 && abs(lid2) == abs(lid3)) {goodFlavor = true; dilMass = (*lep2+*lep3).M(); PtPh = lep1->pt(); Pt1 = lep2->pt(); Pt2 = lep3->pt();}
    }
    else if((cuts & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2){
      if(lq1*lq3 < 0 && abs(lid1) == abs(lid3)) {goodFlavor = true; dilMass = (*lep1+*lep3).M(); PtPh = lep2->pt(); Pt1 = lep1->pt(); Pt2 = lep3->pt();}
    }
    else if((cuts & SmurfTree::Lep3LooseEleV2) == SmurfTree::Lep3LooseEleV2){
      if(lq1*lq2 < 0 && abs(lid1) == abs(lid2)) {goodFlavor = true; dilMass = (*lep1+*lep2).M(); PtPh = lep3->pt(); Pt1 = lep1->pt(); Pt2 = lep2->pt();}
    }
    else {assert(0);}

    if(goodFlavor == false) continue;

    double add = 1.0;
    if(dstype != SmurfTree::data){
      add = add*scale1fb;
      add = add*nPUScaleFactor2012(fhDPU,npu);
      if((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
        add = add*leptonEfficiency(lep1->Pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
      if((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
        add = add*leptonEfficiency(lep2->Pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
      if((cuts & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)
        add = add*leptonEfficiency(lep3->Pt(), lep3->eta(), fhDEffMu, fhDEffEl, lid3);
      add = add*lumi*prescale;
    }

    if(dstype != SmurfTree::data) hBPtPh->Fill(TMath::Min(PtPh,99.999),add); else hDPtPh->Fill(TMath::Min(PtPh,99.999),add);
    if(PtPh <= 90) continue;

    if(dstype != SmurfTree::data) hBM2l->Fill(TMath::Min(dilMass,119.999),add); else hDM2l->Fill(TMath::Min(dilMass,119.999),add);
    if(dilMass <= 12 || dilMass >= 85) continue;

    if(dstype != SmurfTree::data) hBM3l->Fill(TMath::Min((*lep1+*lep2+*lep3).M(),119.999),add); else hDM3l->Fill(TMath::Min((*lep1+*lep2+*lep3).M(),119.999),add);
    if(TMath::Abs((*lep1+*lep2+*lep3).M()-91.1876) > 15.0) continue;
    
    if(dstype != SmurfTree::data) hBPt1->Fill(TMath::Min(Pt1,99.999),add); else hDPt1->Fill(TMath::Min(Pt1,99.999),add);
    if(dstype != SmurfTree::data) hBPt2->Fill(TMath::Min(Pt2,99.999),add); else hDPt2->Fill(TMath::Min(Pt2,99.999),add);
  }

  printf("B/D M2l: %6.2f %3d\n",hBM2l->GetSumOfWeights(),(int)hDM2l->GetSumOfWeights());
  printf("B/D M3l: %6.2f %3d\n",hBM3l->GetSumOfWeights(),(int)hDM3l->GetSumOfWeights());
  printf("B/D Pt1: %6.2f %3d\n",hBPt1->GetSumOfWeights(),(int)hDPt1->GetSumOfWeights());

  TFile *outputFileH0 = new TFile(Form("outputH%d.root",year), "RECREATE");
  outputFileH0->cd();
  hDPtPh->Write();
  hBPtPh->Write();
  hDM3l ->Write();
  hBM3l ->Write();
  hDM2l ->Write();
  hBM2l ->Write();
  hDPt1 ->Write();
  hBPt1 ->Write();
  hDPt2 ->Write();
  hBPt2 ->Write();
  outputFileH0->Close();
}
