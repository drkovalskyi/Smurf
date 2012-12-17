
//
// Dave "the one but not the only" Evans 
//

#include "FakeLooper.h"

// ROOT includes
#include "TROOT.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2F.h"
#include <cmath>

FakeLooper::FakeLooper(Option option)
{
    option_ = option;
}

//
// Main function
//

int FakeLooper::Loop(bool isData, TChain* chain, const char* name, 
        LeptonTree::EventSelection selection, const std::vector<unsigned int>& ptThresholds)
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

    h1_nvtx_ = new TH1F(Form("%s_nvtx", name),       Form("%s_den_lowPt_mt;N_{Vtx}", name),         50, -0.5, 49.5);

    h1_den_lowPt_mt_           = new TH1F(Form("%s_den_lowPt_mt", name),       Form("%s_den_lowPt_mt;MT [GeV]", name),         20, 0.0, 200.0);
    h1_den_lowPt_mt_lowmet_    = new TH1F(Form("%s_den_lowPt_mt_lowmet", name),Form("%s_den_lowPt_mt_lowmet;MT [GeV]", name),  20, 0.0, 200.0);
    h1_den_lowPt_met_          = new TH1F(Form("%s_den_lowPt_met", name),      Form("%s_den_lowPt_met;MET [GeV]", name),       20, 0.0, 100.0);
    h1_den_lowPt_mt_->Sumw2();
    h1_den_lowPt_mt_lowmet_->Sumw2();
    h1_den_lowPt_met_->Sumw2();

    h1_den_highPt_mt_          = new TH1F(Form("%s_den_highPt_mt", name),       Form("%s_den_highPt_mt;MT [GeV]", name),       20, 0.0, 200.0);
    h1_den_highPt_mt_lowmet_   = new TH1F(Form("%s_den_highPt_mt_lowmet", name),Form("%s_den_highPt_mt_lowmet;MT [GeV]", name),20, 0.0, 200.0);
    h1_den_highPt_met_         = new TH1F(Form("%s_den_highPt_met", name),      Form("%s_den_highPt_met;MET [GeV]", name),     20, 0.0, 100.0);
    h1_den_highPt_mt_->Sumw2();
    h1_den_highPt_mt_lowmet_->Sumw2();
    h1_den_highPt_met_->Sumw2();

    h1_num_lowPt_mt_           = new TH1F(Form("%s_num_lowPt_mt", name),       Form("%s_num_lowPt_mt;MT [GeV]", name),         20, 0.0, 200.0);
    h1_num_lowPt_mt_lowmet_    = new TH1F(Form("%s_num_lowPt_mt_lowmet", name),Form("%s_num_lowPt_mt_lowmet;MT [GeV]", name),  20, 0.0, 200.0);
    h1_num_lowPt_met_          = new TH1F(Form("%s_num_lowPt_met", name),      Form("%s_num_lowPt_met;MET [GeV]", name),       20, 0.0, 100.0);
    h1_num_lowPt_mt_->Sumw2();
    h1_num_lowPt_mt_lowmet_->Sumw2();
    h1_num_lowPt_met_->Sumw2();

    h1_num_highPt_mt_          = new TH1F(Form("%s_num_highPt_mt", name),       Form("%s_num_highPt_mt;MT [GeV]", name),       20, 0.0, 200.0);
    h1_num_highPt_mlj_         = new TH1F(Form("%s_num_highPt_mlj", name),      Form("%s_num_highPt_mlj;M(FO,Jet) [GeV]", name),      32, 40.0, 200.0);
    h1_num_highPt_etaj_        = new TH1F(Form("%s_num_highPt_etaj", name),     Form("%s_num_highPt_etaj;#eta(Jet)", name),     25, -5.0, 5.0);
    h1_num_highPt_mll_  = new TH1F(Form("%s_num_highPt_mll", name),      Form("%s_num_highPt_mll;M(FO,Reco Lepton) [GeV]", name),      25, 40.0, 140.0);
    h1_num_highPt_mllss_  = new TH1F(Form("%s_num_highPt_mllss", name),      Form("%s_num_highPt_mllss;M(FO,Reco Lepton) [GeV]", name),      25, 40.0, 140.0);
    h1_num_highPt_emfcentralj_ = new TH1F(Form("%s_num_highPt_emfcentralj", name),     Form("%s_num_highPt_emfcentralj;f(EM) [GeV]", name), 20, 0.0, 1.0);
    h1_num_highPt_cemfcentralj_ = new TH1F(Form("%s_num_highPt_cemfcentralj", name),     Form("%s_num_highPt_cemfcentralj;f(EM) [GeV]", name), 20, 0.0, 1.0);
    h1_num_highPt_nemfcentralj_ = new TH1F(Form("%s_num_highPt_nemfcentralj", name),     Form("%s_num_highPt_nemfcentralj;f(EM) [GeV]", name), 20, 0.0, 1.0);

    h1_num_highPt_mt_lowmet_   = new TH1F(Form("%s_num_highPt_mt_lowmet", name),Form("%s_num_highPt_mt_lowmet;MT [GeV]", name),20, 0.0, 200.0);
    h1_num_highPt_met_         = new TH1F(Form("%s_num_highPt_met", name),      Form("%s_num_highPt_met;MET [GeV]", name),     20, 0.0, 100.0);
    h1_num_highPt_mt_->Sumw2();
    h1_num_highPt_mlj_->Sumw2();
    h1_num_highPt_etaj_->Sumw2();
    h1_num_highPt_mll_->Sumw2();
    h1_num_highPt_mllss_->Sumw2();
    h1_num_highPt_emfcentralj_->Sumw2();
    h1_num_highPt_cemfcentralj_->Sumw2();
    h1_num_highPt_nemfcentralj_->Sumw2();
    h1_num_highPt_mt_lowmet_->Sumw2();
    h1_num_highPt_met_->Sumw2();

    //
    // decide if to loop or not...
    //

    TObjArray *listOfFiles = chain->GetListOfFiles();
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[FakeLooper::loop] " << name << " failed to find any files" << std::endl;
        return 1;
    }
    else
        std::cout << "[FakeLooper::loop] " << name << std::endl;

    if (!isData) {
        std::cout << "[FakeLooper::loop] HLT_Mu8 = " << lumi_HLT_Mu8_ << std::endl;
        std::cout << "[FakeLooper::loop] HLT_Mu17 = " << lumi_HLT_Mu17_ << std::endl;
        std::cout << "[FakeLooper::loop] HLT_Ele8 = " << lumi_HLT_Ele8_ << std::endl;
        std::cout << "[FakeLooper::loop] HLT_Ele17 = " << lumi_HLT_Ele17_ << std::endl;
    }

    //
    // loop over content of sample
    //

    TIter fileIter(listOfFiles);
    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

        LeptonTree *leptons = new LeptonTree();
        leptons->LoadTree(currentFile->GetTitle());
        leptons->InitTree();

        //
        // get extra variables
        //

        // trigger efficiency assumed
        // to be one in MC
        UInt_t HLT_Mu8_probe = 1;
        UInt_t HLT_Mu17_probe = 1;
        UInt_t HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe = 1;
        UInt_t HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe = 1;
        if (isData) {
            leptons->tree_->SetBranchAddress("HLT_Mu8_probe",    &HLT_Mu8_probe);
            leptons->tree_->SetBranchAddress("HLT_Mu17_probe",   &HLT_Mu17_probe);
            leptons->tree_->SetBranchAddress("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe", 
                    &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe);
            leptons->tree_->SetBranchAddress("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe",                 
                    &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe);
        }

        // scale...
        Float_t weight = 1.0;
        Float_t scale1fb = 1.0;
        if (!isData) {
            leptons->tree_->SetBranchAddress("new_scale1fb", &scale1fb);
        }

        // selections
        const unsigned int muonFO = LeptonTree::PassMuFOICHEP2012;
        const unsigned int muonSelection = (LeptonTree::PassMuIsoICHEP2012 | LeptonTree::PassMuIDICHEP2012);
        const unsigned int eleFO = LeptonTree::PassEleFOICHEP2012;
        const unsigned int eleSelection = (LeptonTree::PassEleIsoICHEP2012 | LeptonTree::PassEleIDICHEP2012);

        //
        // loop over events in file
        //

        ULong64_t nEvents = leptons->tree_->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) 
        {

            leptons->tree_->GetEntry(event);

            // muon event
            // make sure we get muon fake events when running on muon datasets...
            if (selection == LeptonTree::QCDFakeMu 
                    && (leptons->eventSelection_ & selection) == selection) {

                // trigger selection
                if (leptons->probe_.Pt() < 20.0) {
                    if (HLT_Mu8_probe == 0) continue;
                    if (!isData) weight = scale1fb * lumi_HLT_Mu8_ / 1000.0;
                }
                else {
                    if (HLT_Mu17_probe == 0) continue;
                    if (!isData) weight = scale1fb * lumi_HLT_Mu17_ / 1000.0;
                }

                bool numerator = false;
                if ((leptons->leptonSelection_ & muonSelection) == muonSelection) numerator = true;

                // muon selection
                if ((leptons->leptonSelection_ & muonFO) != muonFO)                         continue;
                if (sqrt(pow(leptons->probe_.Eta() - leptons->jet1_.Eta(), 2) 
                            + pow(leptons->probe_.Phi() - leptons->jet1_.Phi(), 2)) < 1.0)  continue;

                // general properties histograms for MC normalisation checks
                // for high pt 
                fillValidationHistograms(leptons, weight, numerator, 15.0, true);

                // reduce ewk contamination
                if (option_ == MET20MT15 && (leptons->mt_ >= 15.0 || leptons->met_ >= 20.0))   continue;
                if (option_ == MET20MT15MLL && (leptons->mt_ >= 15.0 || leptons->met_ >= 20.0
                        || (leptons->qProbe_*leptons->qTag_ < 0 && fabs((leptons->probe_+leptons->tag_).M() - 91.0) < 15.0)))   continue;
                if (option_ == MET20 && leptons->met_ >= 20.0)  continue;

                // fill numerator and denominator according to jet threshold
                for (unsigned int i = 0; i < ptThresholds.size(); ++i) {
                    if (leptons->jet1_.Pt() >= double(ptThresholds[i])) {
                        v_denominator[i]->Fill(leptons->probe_.Pt(), fabs(leptons->probe_.Eta()), weight);
                        if (numerator)  v_numerator[i]->Fill(leptons->probe_.Pt(), fabs(leptons->probe_.Eta()), weight); 
                    }
                }
            }

            // electron event
            // make sure we get electron fake events when running on electron datasets...
            if (selection == LeptonTree::QCDFakeEle 
                    && (leptons->eventSelection_ & selection) == selection) {

                // trigger selection
                if (leptons->probe_.Pt() < 20.0) {
                    if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe == 0) continue;
                    if (!isData) weight = scale1fb * lumi_HLT_Ele8_ / 1000.0;
                }
                else {
                    if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe == 0) continue;
                    if (!isData) weight = scale1fb * lumi_HLT_Ele17_ / 1000.0; 
                }

                bool numerator = false;
                if ((leptons->leptonSelection_ & eleSelection) == eleSelection) numerator = true;

                // muon selection   
                if ((leptons->leptonSelection_ & eleFO) != eleFO)                           continue;
                if (sqrt(pow(leptons->probe_.Eta() - leptons->jet1_.Eta(), 2) 
                            + pow(leptons->probe_.Phi() - leptons->jet1_.Phi(), 2)) < 1.0)  continue;

                // general properties histograms for MC normalisation checks
                // for high pt 
                fillValidationHistograms(leptons, weight, numerator, 35.0, false);

                // reduce ewk contamination
                if (option_ == MET20MT15 && (leptons->mt_ >= 15.0 || leptons->met_ >= 20.0))   continue;
                if (option_ == MET20MT15MLL && (leptons->mt_ >= 15.0 || leptons->met_ >= 20.0 
                        || (leptons->qProbe_*leptons->qTag_ < 0 && fabs((leptons->probe_+leptons->tag_).M() - 91.0) < 15.0)
                        || fabs(leptons->jet1_.Eta()) > 2.50))   continue;
                if (option_ == MET20 && leptons->met_ >= 20.0)  continue;
                // fill numerator and denominator according to jet threshold
                for (unsigned int i = 0; i < ptThresholds.size(); ++i) {
                    if (leptons->jet1_.Pt() >= double(ptThresholds[i])) {
                        v_denominator[i]->Fill(leptons->probe_.Pt(), fabs(leptons->probe_.Eta()), weight);
                        if (numerator)  v_numerator[i]->Fill(leptons->probe_.Pt(), fabs(leptons->probe_.Eta()), weight);
                    }
                }

            }


        } // end loop on events

        delete leptons;

    } // end loop on files in chain

    return 0;

}

void FakeLooper::fillValidationHistograms(const LeptonTree* leptons, const float &weight, bool numerator, float jet, bool muonselection)
{

    h1_nvtx_->Fill(leptons->nvtx_, weight);

    // general properties histograms for MC normalisation checks
    // for low pt 
    if (leptons->jet1_.Pt() >= jet && leptons->probe_.Pt() < 20.0) {
        h1_den_lowPt_met_->Fill(leptons->met_, weight);
        if (numerator) h1_num_lowPt_met_->Fill(leptons->met_, weight);
        if (leptons->met_ >= 30.0) {
            h1_den_lowPt_mt_->Fill(leptons->mt_, weight);
            if (numerator) h1_num_lowPt_mt_->Fill(leptons->mt_, weight);
        }
        if (leptons->met_ < 15.0) {
            h1_den_lowPt_mt_lowmet_->Fill(leptons->mt_, weight);
            if (numerator) h1_num_lowPt_mt_lowmet_->Fill(leptons->mt_, weight);
        }
    }

    // for high pt 
    if (leptons->jet1_.Pt() >= jet && leptons->probe_.Pt() >= 20.0) {
        h1_den_highPt_met_->Fill(leptons->met_, weight);
        if (numerator) h1_num_highPt_met_->Fill(leptons->met_, weight);
        if (leptons->met_ >= 30.0) {
            h1_den_highPt_mt_->Fill(leptons->mt_, weight);
            if (numerator) h1_num_highPt_mt_->Fill(leptons->mt_, weight);
        } 
        if (leptons->met_ < 15.0) {
            h1_den_highPt_mt_lowmet_->Fill(leptons->mt_, weight);
            if (numerator) {
                h1_num_highPt_mt_lowmet_->Fill(leptons->mt_, weight);
                // the FO-leading jet invariant mass at low MT
                if (leptons->mt_ < 15.0) {
                    h1_num_highPt_etaj_->Fill(leptons->jet1_.Eta(), weight);
                    h1_num_highPt_mlj_->Fill((leptons->jet1_ + leptons->probe_).M(), weight);
                    if (fabs(leptons->jet1_.Eta()) < 2.5 || muonselection) {

                        h1_num_highPt_emfcentralj_->Fill(leptons->chargedEmFracJet1_ + leptons->neutralEmFracJet1_, weight);
                        h1_num_highPt_cemfcentralj_->Fill(leptons->chargedEmFracJet1_, weight);
                        h1_num_highPt_nemfcentralj_->Fill(leptons->neutralEmFracJet1_, weight);

                        // FO is "probe" and other reco lepton is "tag"
                        bool passMassVeto = true;
                        float mass = TMath::Max(40.5, TMath::Min((leptons->probe_ + leptons->tag_).M(), 139.5));
                        if (leptons->qTag_ * leptons->qProbe_ < 0) {
                            h1_num_highPt_mll_->Fill(mass, weight);
                        } 
                        else if (leptons->qTag_ * leptons->qProbe_ > 0) {
                            h1_num_highPt_mllss_->Fill(mass, weight);
                        }
                        else {
                            h1_num_highPt_mll_->Fill(40.5, weight);
                        }

                    }
                }
            }
        }
    }


}

void FakeLooper::setMuonLumi(float Mu8, float Mu17)
{
    lumi_HLT_Mu8_   = Mu8;
    lumi_HLT_Mu17_  = Mu17;
}

void FakeLooper::setElectronLumi(float Ele8, float Ele17)
{
    lumi_HLT_Ele8_  = Ele8;
    lumi_HLT_Ele17_ = Ele17; 
}

