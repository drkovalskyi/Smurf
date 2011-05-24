/*
   This script skims the common smurf ntuples to select 0-jet and ww preselections for the ME calculation  
   smurfproducer(smurfFDir, "ww.root", outputDir)
   0 - location of the smurf ntuples
   1 - the sample to run over
   2 - the location of the skimmed smurf files
 */
 
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>

using namespace std;

//###################
//# main function
//###################
void smurfproducer(TString smurfFDir = "/smurf/data/Run2011_Spring11_SmurfV3/mitf-alljets/", TString fileName = "zz.root", TString outputDir = "/smurf/data/Run2011_Spring11_SmurfV3/mitf-zerojets-wwselection/") {

  TFile* fin = new TFile(smurfFDir+fileName);
  TFile *newfile = new TFile(outputDir+fileName,"recreate");

  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // Initialize the branches to use to calculate the differential cross-sections
  unsigned int njets_ = 0;
  unsigned int cuts_ = 0;

  ch->SetBranchAddress( "njets"      , &njets_     );     
  ch->SetBranchAddress( "cuts"      , &cuts_     );     

  //==========================================
  // Loop All Events
  //==========================================
  
  cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; ";

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt);            
    if (njets_ != 0) continue;
    // to select the ww-preselections
    if ( (cuts_ & 4915719) != 4915719 )  continue;
    evt_tree->Fill();
  }   //nevent
  
  cout << outputDir + fileName << " has " << evt_tree->GetEntries() << " entries; ";
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
}  
