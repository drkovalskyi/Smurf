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

void MakeGoodRunSample
(
 TString inputFileName  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/data.root",
 TString outputFileName = "output.root",
 string jsonFile        = "json/Cert_Current_JSON.txt"
)
{

  TChain *chbackground = new TChain("tree");
  chbackground->Add(inputFileName);
  TTree *background = (TTree*) chbackground;

  TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputFile->cd();

  TTree *normalizedTree = background->CloneTree(0);

  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 

  UInt_t N_all  = background->GetEntries();
  UInt_t N_good = 0;

  //----------------------------------------------------------------------------
  UInt_t          run;
  UInt_t          lumi;

  //*****************************************************************************
  // MC Loop
  //*****************************************************************************

  background->SetBranchAddress( "run"	       , &run	       );
  background->SetBranchAddress( "lumi"         , &lumi         );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    //****************************************************************************************
    //Good Lumi Selection
    //****************************************************************************************
    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(!rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

    N_good++;
    normalizedTree->Fill(); 

  }
  printf("N good/all = %d / %d = %f\n",N_good,N_all,(double)N_good/N_all);

  normalizedTree->Write();
  outputFile->Close();

}
