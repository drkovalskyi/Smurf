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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void GJScaleFactorB
(
 TString inputZFileName  = "dyll_skim13.root",
 TString inputGFileName  = "gamma.root",
 int year = 2012
)
{

  unsigned int patternTopVeto = SmurfTree::TopVeto;

  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  if	 (year == 2012){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
  }
  else if(year == 2011){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/puWeights_Fall11_42x_True.root";
  }
  else {assert(0);}

  TChain *chinputG = new TChain("tree");
  chinputG->Add(inputGFileName);
  TTree *inputG = (TTree*) chinputG;

  TChain *chinputZ = new TChain("tree");
  chinputZ->Add(inputZFileName);
  TTree *inputZ = (TTree*) chinputZ;

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

  TH1D* hZNvtx   = new TH1D("hZNvtx",  "hZNvtx", 13,  3.5, 29.5); hZNvtx->Sumw2();
  TH1D* hGNvtx   = new TH1D("hGNvtx",  "hGNvtx", 13,  3.5, 29.5); hGNvtx->Sumw2();
  TH1D* hRNvtx   = new TH1D("hRNvtx",  "hRNvtx", 13,  3.5, 29.5); hRNvtx->Sumw2();
  TH1D* hZNjet   = new TH1D("hZNjet",  "hZNjet",  5, -0.5, 4.5);  hZNjet->Sumw2();
  TH1D* hGNjet   = new TH1D("hGNjet",  "hGNjet",  5, -0.5, 4.5);  hGNjet->Sumw2();
  TH1D* hRNjet   = new TH1D("hRNjet",  "hRNjet",  5, -0.5, 4.5);  hRNjet->Sumw2();
  TH1D* hZPtll   = new TH1D("hZPtll",  "hZPtll", 20, 60, 260);    hZPtll->Sumw2();
  TH1D* hGPtll   = new TH1D("hGPtll",  "hGPtll", 20, 60, 260);    hGPtll->Sumw2();
  TH1D* hRPtll   = new TH1D("hRPtll",  "hRPtll", 20, 60, 260);    hRPtll->Sumw2();
  TH1D* hZMet    = new TH1D("hZMet",   "hZMet",  40,  0, 160);    hZMet->Sumw2();
  TH1D* hGMet    = new TH1D("hGMet",   "hGMet",  40,  0, 160);    hGMet->Sumw2();
  TH1D* hRMet    = new TH1D("hRMet",   "hRMet",  40,  0, 160);    hRMet->Sumw2();
  TH1D* hZMT    = new TH1D("hZMT",     "hZMT",   40,  0, 400);    hZMT->Sumw2();
  TH1D* hGMT    = new TH1D("hGMT",     "hGMT",   40,  0, 400);    hGMT->Sumw2();
  TH1D* hRMT    = new TH1D("hRMT",     "hRMT",   40,  0, 400);    hRMT->Sumw2();
  TH1D* hGLid    = new TH1D("hGLid",   "hGLid", 128,-0.5, 127.5); hGLid->Sumw2();
  TH1D* hZEta   = new TH1D("hZEta",  "hZEta",  25, 0, 2.5);  hZEta->Sumw2();
  TH1D* hGEta   = new TH1D("hGEta",  "hGEta",  25, 0, 2.5);  hGEta->Sumw2();
  TH1D* hREta   = new TH1D("hREta",  "hREta",  25, 0, 2.5);  hREta->Sumw2();

  UInt_t  type;
  UInt_t  dstype;
  UInt_t  cuts;
  Int_t   njets;
  Int_t   lid1;
  Int_t   lid2;
  Float_t mt;
  Float_t met;
  Float_t metPhi;
  Float_t npu;
  UInt_t nvtx;
  float scale1fb;
  LorentzVector* dilep = 0;
  LorentzVector* lep1 = 0;
  LorentzVector* lep2 = 0;

  inputZ->SetBranchAddress( "type"	   , &type	  );
  inputZ->SetBranchAddress( "dstype"	   , &dstype	  );
  inputZ->SetBranchAddress( "njets"	   , &njets	  );
  inputZ->SetBranchAddress( "nvtx"	   , &nvtx	  );
  inputZ->SetBranchAddress( "mt"	   , &mt	  );
  inputZ->SetBranchAddress( "met"	   , &met	  );
  inputZ->SetBranchAddress( "metPhi"	   , &metPhi	  );
  inputZ->SetBranchAddress( "npu"	   , &npu	  );
  inputZ->SetBranchAddress( "lid1"	   , &lid1	  );
  inputZ->SetBranchAddress( "lid2"	   , &lid2	  );
  inputZ->SetBranchAddress( "scale1fb"	   , &scale1fb	  );
  inputZ->SetBranchAddress( "cuts"         , &cuts        );
  inputZ->SetBranchAddress( "dilep"        , &dilep       );
  inputZ->SetBranchAddress( "lep1"         , &lep1        );
  inputZ->SetBranchAddress( "lep2"         , &lep2        );
  for (UInt_t i=0; i<inputZ->GetEntries(); i++) {    
    inputZ->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputZ->GetEntries());

    bool lId = (cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
               (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
    if(!lId) continue;
    if(!(lep1->pt() > 20 && lep2->pt() > 20 && fabs(dilep->M()-91.1876) < 15.)) continue;
    if(dilep->pt() <= 60) continue;
    //if(met <= 60) continue;
    //if(TMath::Abs(met-dilep->pt())/dilep->pt() >= 0.25) continue;
    //if(DeltaPhi(dilep->phi(),metPhi)*180.0/TMath::Pi() <= 160) continue;
    if(!((cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto)) continue;
    if(!((cuts & patternTopVeto) == patternTopVeto)) continue;
    if(njets != 0 || TMath::Abs(dilep->Eta()) > 2.3) continue;

    double add = 1.0;
    add = add*scale1fb;
    add = add*nPUScaleFactor2012(fhDPU,npu);

    if(type == SmurfTree::em || type == SmurfTree::me) add = -1.0*add;

    hZNvtx->Fill(TMath::Max(TMath::Min((double)nvtx,29.499),4.0),add);
    hZNjet->Fill(TMath::Min((double)njets,4.499),add);
    hZPtll->Fill(TMath::Min(dilep->pt(),259.999),add);
    hZMet ->Fill(TMath::Min((double)met,159.999),add);
    hZEta ->Fill(TMath::Min(TMath::Abs(dilep->Eta()),2.499),add);
    hZMT  ->Fill(TMath::Min((double)mt,399.999),add);
   }

  inputG->SetBranchAddress( "type"	   , &type	  );
  inputG->SetBranchAddress( "dstype"	   , &dstype	  );
  inputG->SetBranchAddress( "njets"	   , &njets	  );
  inputG->SetBranchAddress( "nvtx"	   , &nvtx	  );
  inputG->SetBranchAddress( "mt"	   , &mt	  );
  inputG->SetBranchAddress( "met"	   , &met	  );
  inputG->SetBranchAddress( "metPhi"	   , &metPhi	  );
  inputG->SetBranchAddress( "npu"	   , &npu	  );
  inputG->SetBranchAddress( "lid1"	   , &lid1	  );
  inputG->SetBranchAddress( "lid2"	   , &lid2	  );
  inputG->SetBranchAddress( "scale1fb"	   , &scale1fb	  );
  inputG->SetBranchAddress( "cuts"         , &cuts        );
  inputG->SetBranchAddress( "dilep"        , &dilep       );
  inputG->SetBranchAddress( "lep1"         , &lep1        );
  inputG->SetBranchAddress( "lep2"         , &lep2        );
  for (UInt_t i=0; i<inputG->GetEntries(); i++) {    
    inputG->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputG->GetEntries());
    if(dilep->pt() <= 60) continue;
    //if(met <= 60) continue;
    //if(TMath::Abs(met-dilep->pt())/dilep->pt() >= 0.25) continue;
    //if(DeltaPhi(dilep->phi(),metPhi)*180.0/TMath::Pi() <= 160) continue;
    if(njets != 0 || (lid1%32) != 31 || TMath::Abs(dilep->Eta()) > 2.3) continue;
    hGNvtx->Fill(TMath::Max(TMath::Min((double)nvtx,29.499),4.0),scale1fb);
  }
  hZNvtx->Scale(1./hZNvtx->GetSumOfWeights());
  hGNvtx->Scale(1./hGNvtx->GetSumOfWeights());
  hRNvtx->Add(hGNvtx);
  hRNvtx->Divide(hZNvtx);

  for (UInt_t i=0; i<inputG->GetEntries(); i++) {    
    inputG->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputG->GetEntries());
    if(dilep->pt() <= 60) continue;
    //if(met <= 60) continue;
    //if(TMath::Abs(met-dilep->pt())/dilep->pt() >= 0.25) continue;
    //if(DeltaPhi(dilep->phi(),metPhi)*180.0/TMath::Pi() <= 160) continue;
    if(njets != 0 || (lid1%32) != 31 || TMath::Abs(dilep->Eta()) > 2.3) continue;

    int bin0 = hRNvtx->FindBin(TMath::Max(TMath::Min((double)nvtx,29.499),4.0));
    double add0 = 1;//hRNvtx->GetBinContent(bin0);
    hGNjet->Fill(TMath::Min((double)njets,4.499),add0*scale1fb);
  }
  hZNjet->Scale(1./hZNjet->GetSumOfWeights());
  hGNjet->Scale(1./hGNjet->GetSumOfWeights());
  hRNjet->Add(hGNjet);
  hRNjet->Divide(hZNjet);

  for (UInt_t i=0; i<inputG->GetEntries(); i++) {    
    inputG->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputG->GetEntries());
    if(dilep->pt() <= 60) continue;
    //if(met <= 60) continue;
    //if(TMath::Abs(met-dilep->pt())/dilep->pt() >= 0.25) continue;
    //if(DeltaPhi(dilep->phi(),metPhi)*180.0/TMath::Pi() <= 160) continue;
    if(njets != 0 || (lid1%32) != 31 || TMath::Abs(dilep->Eta()) > 2.3) continue;

    int bin0 = hRNvtx->FindBin(TMath::Min((double)nvtx,29.499));
    double add0 = 1;//hRNvtx->GetBinContent(bin0);
    int bin1 = hRNjet->FindBin(TMath::Min((double)njets,4.499));
    double add1 = 1;//hRNjet->GetBinContent(bin1);
    hGPtll->Fill(TMath::Min(dilep->pt(),259.999),add0*add1*scale1fb);
  }
  hZPtll->Scale(1./hZPtll->GetSumOfWeights());
  hGPtll->Scale(1./hGPtll->GetSumOfWeights());
  hRPtll->Add(hGPtll);
  hRPtll->Divide(hZPtll);

  for (UInt_t i=0; i<inputG->GetEntries(); i++) {    
    inputG->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputG->GetEntries());
    if(dilep->pt() <= 60) continue;
    //if(met <= 60) continue;
    //if(TMath::Abs(met-dilep->pt())/dilep->pt() >= 0.25) continue;
    //if(DeltaPhi(dilep->phi(),metPhi)*180.0/TMath::Pi() <= 160) continue;
    if(njets != 0 || (lid1%32) != 31 || TMath::Abs(dilep->Eta()) > 2.3) continue;

    int bin0 = hRNvtx->FindBin(TMath::Min((double)nvtx,29.499));
    double add0 = 1;//hRNvtx->GetBinContent(bin0);
    int bin1 = hRNjet->FindBin(TMath::Min((double)njets,4.499));
    double add1 = 1;//hRNjet->GetBinContent(bin1);
    int bin2 = 1;//hRPtll->FindBin(TMath::Min(dilep->pt(),259.999));
    double add2 = hRPtll->GetBinContent(bin2);
    hGMet->Fill(TMath::Min((double)met,159.999),add0*add1*add2*scale1fb);
    hGEta->Fill(TMath::Min(TMath::Abs(dilep->Eta()),2.499),add0*add1*add2*scale1fb);
    hGLid->Fill((double)(lid1%32),add0*add1*add2*scale1fb);
    hGMT->Fill(TMath::Min((double)mt,399.999),add0*add1*add2*scale1fb);
  }
  hZMT->Scale(1./hZMT->GetSumOfWeights());
  hGMT->Scale(1./hGMT->GetSumOfWeights());
  hRMT->Add(hGMT);
  hRMT->Divide(hZMT);
  hZMet->Scale(1./hZMet->GetSumOfWeights());
  hGMet->Scale(1./hGMet->GetSumOfWeights());
  hRMet->Add(hGMet);
  hRMet->Divide(hZMet);
  hZEta->Scale(1./hZEta->GetSumOfWeights());
  hGEta->Scale(1./hGEta->GetSumOfWeights());
  hREta->Add(hGEta);
  hREta->Divide(hZEta);
  hGLid->Scale(1./hGLid->GetSumOfWeights());

/*
  TFile *outputFileG = new TFile(Form("outputG%d.root",year), "RECREATE");
  outputFileG->cd();
  TTree *normalizedTreeG = inputG->CloneTree(0);

  for (UInt_t i=0; i<inputG->GetEntries(); i++) {    
    inputG->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)inputG->GetEntries());
    if(dilep->pt() <= 60) continue;
    if(met <= 60) continue;
    if(njets > 1) continue;
    double add = scale1fb;
      hGMet ->Fill(TMath::Min((double)met,199.999),add);
      scale1fb = scale1fb;
      dstype = dstype;
    normalizedTreeG->Fill();
  }
  normalizedTreeG->Write();
  outputFileG->Close();
*/
  //hZMet ->Scale(1./hZMet ->GetSumOfWeights());
  //hGMet ->Scale(1./hGMet ->GetSumOfWeights());

  TFile *outputFileH0 = new TFile(Form("outputH%d.root",year), "RECREATE");
  outputFileH0->cd();
  hZNvtx  ->Write();
  hGNvtx  ->Write();
  hRNvtx  ->Write();
  hZNjet  ->Write();
  hZNjet  ->Write();
  hRNjet  ->Write();
  hZPtll  ->Write();
  hGPtll  ->Write();
  hRPtll  ->Write();
  hZMT    ->Write();
  hGMT    ->Write();
  hRMT    ->Write();
  hZMet   ->Write();
  hGMet   ->Write();
  hRMet   ->Write();
  hZEta   ->Write();
  hGEta   ->Write();
  hREta   ->Write();
  hGLid   ->Write();
  outputFileH0->Close();
}
