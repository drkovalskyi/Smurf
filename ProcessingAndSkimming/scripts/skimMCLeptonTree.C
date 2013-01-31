
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TString.h"
#include "TROOT.h"

// 0 for even
// 1 for odd

void skimMCLeptonTree (unsigned int option, TString inputDir, TString outputFile, Float_t xs, Float_t eventsRead)
{

    TFile* fin = new TFile(inputDir+"/merged.root");
    TTree* ch=(TTree*)fin->Get("leptons"); 
    if (ch==0x0) return; 

    TFile *newfile= new TFile(inputDir+"/"+outputFile,"recreate");
    TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

    unsigned int eventSelection_ = 0;
    ch->SetBranchAddress( "eventSelection"          , &eventSelection_     );
    //unsigned int evt_ = 0;
    //ch->SetBranchAddress( "evt" ,&evt_ );

    // set new scale1fb
    Float_t scale1fb = (xs * 1000.0) / eventsRead;
    std::cout << "scale1fb is " << scale1fb << std::endl;
    evt_tree->Branch("new_scale1fb", &scale1fb, "new_scale1fb/F");
    //evt_tree->SetAlias("scale1fb",  "scale1fb");

    for(int ievt = 0; ievt < ch->GetEntries() ;ievt++) {
        ch->GetEntry(ievt); 
        if ((eventSelection_ & (1<<0)) != (1<<0) && option == 0) continue;
        if ((eventSelection_ & (1<<1)) != (1<<1) && option == 1) continue;
        //if ((odd && (evt_ % 2) == 0) || (!odd && (evt_ %2) != 0)) continue;

        // muon fakes
        if (option == 5) {
            if ((eventSelection_ & (1<<5)) != (1<<5)) continue;
        }

        // electron fakes
        if (option == 4) {
            if ((eventSelection_ & (1<<4)) != (1<<4)) continue;
        }

        evt_tree->Fill();
    }

    newfile->cd(); 
    evt_tree->Write(); 
    newfile->Close();
}  

