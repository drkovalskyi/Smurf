/**********************************************************************************
 * Usage
 *
 * Parameter 1 is the location and leading pattern in the file names. eg. if all files are
 * of the form ntuples_###train_#jets_hww###.root, then the pattern should be "ntuples"
 *
 * Parameter 2 is the name of the output file
 *
 * Example Usage:
 * root -l -b -q MergeMVABranches.C+\(\"/data/smurf/ceballos/test/ntuples\",\"MergeOutput\"\)
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TChainElement.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


//--------------------------------------------------------------------

void MergeMVABranches(
  string inputFilenameHeader = "/data/smurf/ceballos/test/ntuples", 
  string outputFilename = "testMerge.root"
  ) {   

  //*******************************************************************
  //Define MassPoints and Jet Bins
  //*******************************************************************
  vector<string> higgsMassPoints;
  higgsMassPoints.push_back("130");
  vector<string> jetBins;
  jetBins.push_back("0");
  jetBins.push_back("1");

  vector<string> labels;
  vector<string> inputFilenames;
  for(UInt_t i=0; i<higgsMassPoints.size();++i) {
    for(UInt_t j=0; j<jetBins.size();++j) {    
      inputFilenames.push_back(inputFilenameHeader + "_" + higgsMassPoints[i] + "train_" + jetBins[j] + "jets" + "_hww" + higgsMassPoints[i]+".root"); 
      labels.push_back("hww"+higgsMassPoints[i] + "_" + jetBins[j] + "jet_ww");
    }
  }

  //*******************************************************************
  //Define branches to be merged
  //*******************************************************************
  vector<TFile*> files;
  vector<TTree*> trees;
  vector<TBranch*> br_bdt;
  vector<TBranch*> br_bdtd;
  vector<TBranch*> br_nn;
  vector<TBranch*> br_knn;
  vector<TBranch*> br_bdtg;
  vector<TBranch*> br_bdtg_aux0;
  vector<TBranch*> br_bdtg_aux1;
  vector<TBranch*> br_bdtg_aux2;
  vector<Float_t> bdt;
  vector<Float_t> bdtd;
  vector<Float_t> nn;
  vector<Float_t> knn;
  vector<Float_t> bdtg;
  vector<Float_t> bdtg_aux0;
  vector<Float_t> bdtg_aux1;
  vector<Float_t> bdtg_aux2;


  //*******************************************************************
  //Load input trees
  //*******************************************************************
  for(UInt_t i=0; i<labels.size();++i) {
    cout << "file: " << inputFilenames[i] << endl;
    TFile *tmpFile = new TFile( inputFilenames[i].c_str(), "READ");
    assert(tmpFile);
    TTree *tmpTree = (TTree*)tmpFile->Get("tree");
    assert(tmpTree);

    files.push_back(tmpFile);
    trees.push_back(tmpTree);
  }


  //*******************************************************************
  //Create output tree
  //*******************************************************************
  TFile *f = TFile::Open( (inputFilenameHeader + "_" + higgsMassPoints[0] + "train_" + jetBins[0] + "jets" + "_hww" + higgsMassPoints[0]+".root").c_str() , "READ");
  assert(f);
  TTree* t = (TTree*)f->Get("tree");
  assert(t);
  t->SetBranchStatus("*", 1);
  TFile *out = TFile::Open( outputFilename.c_str() ,"RECREATE" );
  TTree *clone;
  clone = t->CloneTree(-1, "fast");


  //*******************************************************************
  //Define branches to be merged
  //*******************************************************************
  vector<TFile*> files;
  vector<TTree*> trees;
  vector<TBranch*> br_bdt;
  vector<TBranch*> br_bdtd;
  vector<TBranch*> br_nn;
  vector<TBranch*> br_knn;
  vector<TBranch*> br_bdtg;
  vector<TBranch*> br_bdtg_aux0;
  vector<TBranch*> br_bdtg_aux1;
  vector<TBranch*> br_bdtg_aux2;
  vector<Float_t> bdt;
  vector<Float_t> bdtd;
  vector<Float_t> nn;
  vector<Float_t> knn;
  vector<Float_t> bdtg;
  vector<Float_t> bdtg_aux0;
  vector<Float_t> bdtg_aux1;
  vector<Float_t> bdtg_aux2;

  for(UInt_t i=0; i<labels.size();++i) {

    bdt.push_back(-1);
    bdtd.push_back(-1);
    nn.push_back(-1);
    knn.push_back(-1);
    bdtg.push_back(-1);
    bdtg_aux0.push_back(-1);
    bdtg_aux1.push_back(-1);
    bdtg_aux2.push_back(-1);
    
    TBranch* tmpbr_bdt        = 0;
    TBranch* tmpbr_bdtd       = 0;
    TBranch* tmpbr_nn         = 0;
    TBranch* tmpbr_knn        = 0;
    TBranch* tmpbr_bdtg       = 0;
    TBranch* tmpbr_bdtg_aux0  = 0;
    TBranch* tmpbr_bdtg_aux1  = 0;
    TBranch* tmpbr_bdtg_aux2  = 0;
    
    tmpbr_bdt        = clone->Branch(("bdt_" +labels[i]).c_str()         , &bdt[i]       , ("bdt_" +labels[i]+"/F").c_str());
    tmpbr_bdtd       = clone->Branch(("bdtd_"+labels[i]).c_str()         , &bdtd[i]      , ("bdtd_"+labels[i]+"/F").c_str()); 
    tmpbr_nn         = clone->Branch(("nn_"  +labels[i]).c_str()         , &nn[i]        , ("nn_"  +labels[i]+"/F").c_str()); 
    tmpbr_knn        = clone->Branch(("knn_" +labels[i]).c_str()         , &knn[i]       , ("knn_" +labels[i]+"/F").c_str()); 
    tmpbr_bdtg       = clone->Branch(("bdtg_"+labels[i]).c_str()         , &bdtg[i]      , ("bdtg_"+labels[i]+"/F").c_str());     
    tmpbr_bdtg_aux0  = clone->Branch(("bdtg_"+labels[i]+"_aux0").c_str() , &bdtg_aux0[i] , ("bdtg_"+labels[i]+"_aux0/F").c_str() ); 
    tmpbr_bdtg_aux1  = clone->Branch(("bdtg_"+labels[i]+"_aux1").c_str() , &bdtg_aux1[i] , ("bdtg_"+labels[i]+"_aux0/F").c_str() ); 
    tmpbr_bdtg_aux2  = clone->Branch(("bdtg_"+labels[i]+"_aux2").c_str() , &bdtg_aux2[i] , ("bdtg_"+labels[i]+"_aux0/F").c_str() ); 
    
    tmpbr_bdt       -> SetTitle(("BDT Output "    + labels[i]).c_str());
    tmpbr_bdt       -> SetTitle(("BDTD Output "   + labels[i]).c_str()); 
    tmpbr_nn        -> SetTitle(("MLPBNN Output " + labels[i]).c_str());
    tmpbr_knn       -> SetTitle(("KNN Output "    + labels[i]).c_str());  
    tmpbr_bdtg      -> SetTitle(("BDTG Output "   + labels[i]).c_str()); 

    br_bdt.push_back(tmpbr_bdt);
    br_bdtd.push_back(tmpbr_bdtd);
    br_nn.push_back(tmpbr_nn);
    br_knn.push_back(tmpbr_knn);
    br_bdtg.push_back(tmpbr_bdtg);
    br_bdtg_aux0.push_back(tmpbr_bdtg_aux0);
    br_bdtg_aux1.push_back(tmpbr_bdtg_aux1);
    br_bdtg_aux2.push_back(tmpbr_bdtg_aux2);  
  }





  //*******************************************************************
  //Event Loop
  //*******************************************************************
  vector<Float_t> var_bdt;
  vector<Float_t> var_bdtd;
  vector<Float_t> var_nn;
  vector<Float_t> var_knn;
  vector<Float_t> var_bdtg;
  vector<Float_t> var_bdtg_aux0;
  vector<Float_t> var_bdtg_aux1;
  vector<Float_t> var_bdtg_aux2;

  for(UInt_t i=0; i<labels.size();++i) {
    var_bdt.push_back(0);
    var_bdtd.push_back(0);
    var_nn.push_back(0);
    var_knn.push_back(0);
    var_bdtg.push_back(0);
    var_bdtg_aux0.push_back(0);
    var_bdtg_aux1.push_back(0);
    var_bdtg_aux2.push_back(0);
    (trees[i])->SetBranchAddress( ("bdt_" +labels[i]).c_str()  , &var_bdt[i]  );
    (trees[i])->SetBranchAddress( ("bdtd_"+labels[i]).c_str()  , &var_bdtd[i] );
    (trees[i])->SetBranchAddress( ("nn_"  +labels[i]).c_str()  , &var_nn[i]   );
    (trees[i])->SetBranchAddress( ("knn_" +labels[i]).c_str()  , &var_knn[i]  );
    (trees[i])->SetBranchAddress( ("bdtg_"+labels[i]).c_str()  , &var_bdtg[i] );
    (trees[i])->SetBranchAddress( ("bdtg_"+labels[i]+"_aux0").c_str() , &var_bdtg_aux0[i] );
    (trees[i])->SetBranchAddress( ("bdtg_"+labels[i]+"_aux1").c_str() , &var_bdtg_aux1[i] );
    (trees[i])->SetBranchAddress( ("bdtg_"+labels[i]+"_aux2").c_str() , &var_bdtg_aux2[i] );
  }

  for (Long64_t ievt=0; ievt<trees[0]->GetEntries();ievt++) {

    if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

    //Get Event Entry for all trees
    for(UInt_t i=0; i<labels.size();++i) {
      (trees[i])->GetEntry(ievt);
    }
      
    //Fill all MVA Branches
    for(UInt_t i=0; i<labels.size();++i) {
      bdt[i] = var_bdt[i];
      bdtd[i] = var_bdtd[i];
      nn[i] = var_nn[i];
      knn[i] = var_knn[i];
      bdtg[i] = var_bdtg[i];
      bdtg_aux0[i] = var_bdtg_aux0[i];
      bdtg_aux1[i] = var_bdtg_aux1[i];
      bdtg_aux2[i] = var_bdtg_aux2[i];        
 
      br_bdt[i]->Fill();
      br_bdtd[i]->Fill();
      br_nn[i]->Fill();
      br_knn[i]->Fill();
      br_bdtg[i]->Fill();
      br_bdtg_aux0[i]->Fill();
      br_bdtg_aux1[i]->Fill();
      br_bdtg_aux2[i]->Fill();

    }
  }

  //*******************************************************************
  // Write output & Close input files
  //*******************************************************************
  clone->Write(); 
  out->Close();
  f->Close();  
  delete out;
  delete f;
  
  for(UInt_t i=0; i<files.size();++i) {
    files[i]->Close();
    delete files[i];
  }
  
}
