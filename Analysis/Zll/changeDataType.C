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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void changeDataType
(
 TString inputFileName  = "input.root",
 TString outputFileName = "output.root",
 Int_t   newDataType    = SmurfTree::qqbarh
)
{

  TChain *chbackground = new TChain("tree");
  chbackground->Add(inputFileName);
  TTree *background = (TTree*) chbackground;

  TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputFile->cd();

  TTree *normalizedTree = background->CloneTree(0);

  Int_t dstype;

  background->SetBranchAddress( "dstype" , &dstype );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    dstype = newDataType;
    normalizedTree->Fill(); 
  }

  normalizedTree->Write();
  outputFile->Close();

}
