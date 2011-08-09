//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  3 15:36:55 2011 by ROOT version 5.26/00b
// from TTree h10/h10
// found on file: WWqqbr_tota_mstw8nl_80__80__average_scale.root
//////////////////////////////////////////////////////////

#ifndef mcfm2smurf_h
#define mcfm2smurf_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class SmurfTree;

class mcfm2smurf {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         px3;
   Float_t         py3;
   Float_t         pz3;
   Float_t         E_3;
   Float_t         px4;
   Float_t         py4;
   Float_t         pz4;
   Float_t         E_4;
   Float_t         px5;
   Float_t         py5;
   Float_t         pz5;
   Float_t         E_5;
   Float_t         px6;
   Float_t         py6;
   Float_t         pz6;
   Float_t         E_6;
   Float_t         px7;
   Float_t         py7;
   Float_t         pz7;
   Float_t         E_7;
   Float_t         wt_ALL;
   Float_t         wt_gg;
   Float_t         wt_gq;
   Float_t         wt_qq;
   Float_t         wt_qqb;

   // List of branches
   TBranch        *b_px3;   //!
   TBranch        *b_py3;   //!
   TBranch        *b_pz3;   //!
   TBranch        *b_E_3;   //!
   TBranch        *b_px4;   //!
   TBranch        *b_py4;   //!
   TBranch        *b_pz4;   //!
   TBranch        *b_E_4;   //!
   TBranch        *b_px5;   //!
   TBranch        *b_py5;   //!
   TBranch        *b_pz5;   //!
   TBranch        *b_E_5;   //!
   TBranch        *b_px6;   //!
   TBranch        *b_py6;   //!
   TBranch        *b_pz6;   //!
   TBranch        *b_E_6;   //!
   TBranch        *b_px7;   //!
   TBranch        *b_py7;   //!
   TBranch        *b_pz7;   //!
   TBranch        *b_E_7;   //!
   TBranch        *b_wt_ALL;   //!
   TBranch        *b_wt_gg;   //!
   TBranch        *b_wt_gq;   //!
   TBranch        *b_wt_qq;   //!
   TBranch        *b_wt_qqb;   //!

   mcfm2smurf(TTree *tree=0);
   virtual ~mcfm2smurf();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     MakeSmurfNtuple(const char* filename, int mcfm_process_id);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool     FillProcess(SmurfTree& tree, int mcfm_process_id);
};

#endif

#ifdef mcfm2smurf_cxx
mcfm2smurf::mcfm2smurf(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("WWqqbr_tota_mstw8nl_80__80__average_scale.root");
      if (!f) {
         f = new TFile("WWqqbr_tota_mstw8nl_80__80__average_scale.root");
      }
      tree = (TTree*)gDirectory->Get("h10");

   }
   Init(tree);
}

mcfm2smurf::~mcfm2smurf()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mcfm2smurf::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mcfm2smurf::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mcfm2smurf::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("px3", &px3, &b_px3);
   fChain->SetBranchAddress("py3", &py3, &b_py3);
   fChain->SetBranchAddress("pz3", &pz3, &b_pz3);
   fChain->SetBranchAddress("E_3", &E_3, &b_E_3);
   fChain->SetBranchAddress("px4", &px4, &b_px4);
   fChain->SetBranchAddress("py4", &py4, &b_py4);
   fChain->SetBranchAddress("pz4", &pz4, &b_pz4);
   fChain->SetBranchAddress("E_4", &E_4, &b_E_4);
   fChain->SetBranchAddress("px5", &px5, &b_px5);
   fChain->SetBranchAddress("py5", &py5, &b_py5);
   fChain->SetBranchAddress("pz5", &pz5, &b_pz5);
   fChain->SetBranchAddress("E_5", &E_5, &b_E_5);
   fChain->SetBranchAddress("px6", &px6, &b_px6);
   fChain->SetBranchAddress("py6", &py6, &b_py6);
   fChain->SetBranchAddress("pz6", &pz6, &b_pz6);
   fChain->SetBranchAddress("E_6", &E_6, &b_E_6);
   fChain->SetBranchAddress("px7", &px7, &b_px7);
   fChain->SetBranchAddress("py7", &py7, &b_py7);
   fChain->SetBranchAddress("pz7", &pz7, &b_pz7);
   fChain->SetBranchAddress("E_7", &E_7, &b_E_7);
   fChain->SetBranchAddress("wt_ALL", &wt_ALL, &b_wt_ALL);
   fChain->SetBranchAddress("wt_gg", &wt_gg, &b_wt_gg);
   fChain->SetBranchAddress("wt_gq", &wt_gq, &b_wt_gq);
   fChain->SetBranchAddress("wt_qq", &wt_qq, &b_wt_qq);
   fChain->SetBranchAddress("wt_qqb", &wt_qqb, &b_wt_qqb);
   Notify();
}

Bool_t mcfm2smurf::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mcfm2smurf::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mcfm2smurf::Cut(Long64_t entry)
{
  assert(entry>=0);
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mcfm2smurf_cxx
