#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom.h"

using namespace std;

void addBranch(const char* infname, const char* outfile){
  TFile *f = TFile::Open(infname);
  assert(f);
  TTree* t = (TTree*)f->Get("Events");
  assert(t);
  t->SetBranchStatus("*", 1);

  TFile *out = TFile::Open(outfile, "RECREATE");
  TTree *clone;
  clone = t->CloneTree(-1, "fast");
   
  Float_t r(.0);

  TBranch* br = clone->Branch("random", &r, "random/F");
  br->SetTitle("Random number");
  clone->SetAlias("rnd",  "random");

  Int_t nentries = t->GetEntries();
  for(Int_t i = 0; i < nentries; i++) {
    r = gRandom->Rndm();
    br->Fill();
  }

  clone->Write(); 
  out->Close();
  f->Close();
  return;
}
