
//
// Dave "the one but not the only" Evans 
//

#include "ZFakeLooper.h"

// ROOT includes
#include "TROOT.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2F.h"
#include <cmath>

ZFakeLooper::ZFakeLooper()
{
}

//
// Main function
//

int ZFakeLooper::Loop(bool isData, TChain* chain, const char* name, 
        const std::vector<unsigned int>& ptThresholds)
{

    //
    // setup histograms
    //

    gROOT->cd();
    double binsX[6] = {10, 15, 20, 25, 30, 35};     unsigned int nX = 5;
    double binsY[5] = {0.0, 1.0, 1.479, 2.0, 2.5};  unsigned int nY = 4;
    std::vector<TH2F*> v_denominator;
    std::vector<TH2F*> v_numerator;
    for (unsigned int i = 0; i < ptThresholds.size(); ++i) {
        const char* den_name = Form("den_%s_ptThreshold%u_PtEta", name, ptThresholds[i]);
        const char* num_name = Form("num_%s_ptThreshold%u_PtEta", name, ptThresholds[i]);
        v_denominator.push_back(new TH2F(den_name, Form("%s;p_{T} [GeV]; |#eta|", den_name), nX, binsX, nY, binsY));
        v_numerator.push_back(new TH2F(num_name, Form("%s;p_{T} [GeV]; |#eta|", den_name), nX, binsX, nY, binsY));
        v_denominator[i]->Sumw2();
        v_numerator[i]->Sumw2();
        v_denominator[i]->SetMarkerSize(2);
        v_numerator[i]->SetMarkerSize(2);
    }

    //
    // decide if to loop or not...
    //

    TObjArray *listOfFiles = chain->GetListOfFiles();
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[ZFakeLooper::loop] " << name << " failed to find any files" << std::endl;
        return 1;
    }
    else
        std::cout << "[ZFakeLooper::loop] " << name << std::endl;

    if (!isData) {
        std::cout << "[ZFakeLooper::loop] lumi = " << lumi_ << std::endl;
    }

    unsigned int nEventsChain=0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    TH1F *h1_met = new TH1F(Form("%s_h1_met", name), "met", 20, 0, 100);
    h1_met->Sumw2();

    //
    // loop over content of sample
    //

    TIter fileIter(listOfFiles);
    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

        SmurfTree *smurf = new SmurfTree();
        smurf->LoadTree(currentFile->GetTitle());
        smurf->InitTree();

        // extra variables for the various correction factors
        float sfWeightPU_ = 1.0;
        float sfWeightTrig_ = 1.0;
        float sfWeightEff_ = 1.0;

        if (smurf->tree_->GetBranchStatus("sfWeightPU"))
            smurf->tree_->SetBranchAddress("sfWeightPU", &sfWeightPU_);
        if (smurf->tree_->GetBranchStatus("sfWeightTrig"))
            smurf->tree_->SetBranchAddress("sfWeightTrig", &sfWeightTrig_);
        if (smurf->tree_->GetBranchStatus("sfWeightEff"))
            smurf->tree_->SetBranchAddress("sfWeightEff", &sfWeightEff_);

        //
        // loop over events in file
        //

        ULong64_t nEvents = smurf->tree_->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) 
        {

            smurf->tree_->GetEntry(event);

            ++nEventsTotal;
            int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            } 

            //
            // event selection
            //

            // Baseline
            const unsigned int basic_selection =
                SmurfTree::BaseLine | SmurfTree::ChargeMatch |
                SmurfTree::Lep1FullSelection | SmurfTree::Lep2FullSelection
                    | SmurfTree::TopVeto;
            if ((smurf->cuts_ & basic_selection) != basic_selection) continue;

            // Z event with zero jets
            if (fabs(smurf->dilep_.mass() - 91) > 15.0) continue;
            if (smurf->njets_ > 0) continue;
            if (smurf->lep1_.Pt() < 20. || smurf->lep2_.Pt() < 20.) continue;

            // trigger
            float weight = 1.0;
            if (smurf->dstype_ == SmurfTree::data) {
                if ( (smurf->cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger ) continue;
            }
            else weight = lumi_ * smurf->scale1fb_ * sfWeightTrig_ * sfWeightEff_ * sfWeightPU_ / 1000.0;

            // lep3 is FO
            bool denominatorOnly = false;
            if ((smurf->cuts_ & SmurfTree::Lep3LooseMuV2) == SmurfTree::Lep3LooseMuV2) denominatorOnly = true;            

            // lep3 is full lepton
            bool numerator = false;
            if ((smurf->cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection) numerator = true;

            // reject events that are not denominator
            if (!(denominatorOnly || numerator)) continue;

            // fill met histo
            if (numerator) h1_met->Fill(smurf->met_, weight);

            // reduce diboson backgrounds
            if (smurf->met_ > 20) continue;

            //
            // fill the histograms
            //

            for (unsigned int i = 0; i < ptThresholds.size(); ++i) {
                if (smurf->dilep_.Pt() >= double(ptThresholds[i])) {
                    v_denominator[i]->Fill(smurf->lep3_.Pt(), fabs(smurf->lep3_.Eta()), weight);
                    if (numerator)  v_numerator[i]->Fill(smurf->lep3_.Pt(), fabs(smurf->lep3_.Eta()), weight);
                }
            }

        } // end loop on events


        delete smurf;

    } // end loop on files in chain

    return 0;

}

void ZFakeLooper::setLumi(float lumi)
{
    lumi_   = lumi;
}

