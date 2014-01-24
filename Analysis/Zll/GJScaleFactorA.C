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

// first run Smurf/Analysis/Zll/ZllGSel.C to check prescale
//root -l -q -b Smurf/Analysis/Zll/GJScaleFactorA.C+'("ntuples_zh_53x/mitf-alljets_noweights/gamma50.root","ntuples_zh_53x/mitf-alljets_noweights/gamma75.root","ntuples_zh_53x/mitf-alljets_noweights/gamma90.root",2012)'
//root -l -q -b Smurf/Analysis/Zll/GJScaleFactorA.C+'("ntuples_zh_42x/mitf-alljets_noweights/gamma50.root","ntuples_zh_42x/mitf-alljets_noweights/gamma75.root","ntuples_zh_42x/mitf-alljets_noweights/gamma90.root",2011)'
// hadd gamma.root output?.root
// or
//root -l -q -b Smurf/Analysis/Zll/GJScaleFactorA.C+'("Ana/nt_scripts/ntuples_53x/gamma50.root","Ana/nt_scripts/ntuples_53x/gamma75.root","Ana/nt_scripts/ntuples_53x/gamma90.root",2012)'
//root -l -q -b Smurf/Analysis/Zll/GJScaleFactorA.C+'("Ana/nt_scripts/ntuples_42x_v9/gamma50.root","Ana/nt_scripts/ntuples_42x_v9/gamma75.root","Ana/nt_scripts/ntuples_42x_v9/gamma90.root",2011)'
// hadd gamma.root output?.root

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void GJScaleFactorA
(
 TString input1FileName  = "input1.root",
 TString input2FileName  = "input2.root",
 TString input3FileName  = "input3.root",
 int year = 2012
)
{

  TChain *chinput1 = new TChain("tree");
  chinput1->Add(input1FileName);
  TTree *input1 = (TTree*) chinput1;

  TChain *chinput2 = new TChain("tree");
  chinput2->Add(input2FileName);
  TTree *input2 = (TTree*) chinput2;

  TChain *chinput3 = new TChain("tree");
  chinput3->Add(input3FileName);
  TTree *input3 = (TTree*) chinput3;

  double NOut[3] = {0., 0., 0.};
  float scale1fb;
  LorentzVector* dilep = 0;
  Int_t  dstype = 0;

  input1->SetBranchAddress( "dilep"        , &dilep       );
  input1->SetBranchAddress( "dstype"	   , &dstype	  );
  for (UInt_t i=0; i<input1->GetEntries(); i++) {    
    input1->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)input1->GetEntries());
    if(dilep->pt() > 100 && dilep->pt() < 300 && dstype == SmurfTree::data) NOut[0]++;
  }

  input2->SetBranchAddress( "dilep"        , &dilep       );
  input2->SetBranchAddress( "dstype"	   , &dstype	  );
  for (UInt_t i=0; i<input2->GetEntries(); i++) {    
    input2->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)input2->GetEntries());
    if(dilep->pt() > 100 && dilep->pt() < 300 && dstype == SmurfTree::data) NOut[1]++;
  }

  input3->SetBranchAddress( "dilep"        , &dilep       );
  input3->SetBranchAddress( "dstype"	   , &dstype	  );
  for (UInt_t i=0; i<input3->GetEntries(); i++) {    
    input3->GetEntry(i);
    if (i%1000000 == 0) printf("--- reading event %5d of %5d\n",i,(int)input3->GetEntries());
    if(dilep->pt() > 100 && dilep->pt() < 300 && dstype == SmurfTree::data) NOut[2]++;
  }
  printf("N1: %0.f | N2: %0.f | N3: %0.f\n",NOut[0],NOut[1],NOut[2]);
  printf("N3/N1: %f | N3/N2: %f\n",NOut[2]/NOut[0],NOut[2]/NOut[1]);

  // file 1
  TFile *outputFile1 = new TFile("output1.root", "RECREATE");
  outputFile1->cd();
  TTree *normalizedTree1 = input1->CloneTree(0);

  input1->SetBranchAddress( "scale1fb"	   , &scale1fb	   );
  input1->SetBranchAddress( "dilep"        , &dilep        );
  input1->SetBranchAddress( "dstype"	   , &dstype	   );
  for (UInt_t i=0; i<input1->GetEntries(); i++) {    
    input1->GetEntry(i);
    if(dilep->pt() > 60 && dilep->pt() <= 85 && dstype == SmurfTree::data) {
      scale1fb = NOut[2]/NOut[0];
      normalizedTree1->Fill();
    }
  }

  normalizedTree1->Write();
  outputFile1->Close();

  // file 2
  TFile *outputFile2 = new TFile("output2.root", "RECREATE");
  outputFile2->cd();
  TTree *normalizedTree2 = input2->CloneTree(0);

  input2->SetBranchAddress( "scale1fb"	   , &scale1fb	   );
  input2->SetBranchAddress( "dilep"        , &dilep        );
  input2->SetBranchAddress( "dstype"	   , &dstype	   );
  for (UInt_t i=0; i<input2->GetEntries(); i++) {    
    input2->GetEntry(i);
    if(dilep->pt() > 85 && dilep->pt() <= 100 && dstype == SmurfTree::data) {
      scale1fb = NOut[2]/NOut[1];
      normalizedTree2->Fill();
    }
  }

  normalizedTree2->Write();
  outputFile2->Close();

  // file 3
  TFile *outputFile3 = new TFile("output3.root", "RECREATE");
  outputFile3->cd();
  TTree *normalizedTree3 = input3->CloneTree(0);

  double prescaleMC = 0.02 * 19.365;
  if(year == 2011){
    prescaleMC = 0.02 * 4.924;
  }
  input3->SetBranchAddress( "scale1fb"	   , &scale1fb	   );
  input3->SetBranchAddress( "dilep"        , &dilep        );
  input3->SetBranchAddress( "dstype"	   , &dstype	   );
  for (UInt_t i=0; i<input3->GetEntries(); i++) {    
    input3->GetEntry(i);
    if(dilep->pt() > 100 || (dstype != SmurfTree::data && dilep->pt() > 60)) {
      if(dstype == SmurfTree::data) scale1fb = +1.0;
      else                          scale1fb = -1.0 * scale1fb * prescaleMC;
      dstype = SmurfTree::data;
      normalizedTree3->Fill();
    }
  }

  normalizedTree3->Write();
  outputFile3->Close();
}
