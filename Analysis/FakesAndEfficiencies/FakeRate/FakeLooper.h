
#ifndef FAKELOOPER_H
#define FAKELOOPER_H

#include "Enums.h"

// C++ includes
#include <iostream>
#include <vector>

#include "../../../Core/LeptonTree.h"
#include "TChain.h"
#include "TH1F.h"

class FakeLooper {
    public:
        FakeLooper(Option option);
        ~FakeLooper() {};
        int Loop(bool isData, TChain* chain, const char* name, 
            LeptonTree::EventSelection selection, const std::vector<unsigned int>& ptThresholds);

        void setMuonLumi(float Mu8, float Mu17);
        void setElectronLumi(float Ele8, float Ele17);

    private:

void fillValidationHistograms(const LeptonTree* leptons, 
    const float &weight, bool numerator, float jet, bool muonselection);

        Option option_;
        float lumi_HLT_Mu8_;
        float lumi_HLT_Mu17_;
        float lumi_HLT_Ele8_;
        float lumi_HLT_Ele17_;

    // histograms
    TH1F *h1_nvtx_;
    
    TH1F *h1_den_lowPt_mt_;
    TH1F *h1_den_lowPt_mt_lowmet_;
    TH1F *h1_den_lowPt_met_;
    
    TH1F *h1_den_highPt_mt_;
    TH1F *h1_den_highPt_mt_lowmet_;
    TH1F *h1_den_highPt_met_;
        
    TH1F *h1_num_lowPt_mt_;
    TH1F *h1_num_lowPt_mt_lowmet_;
    TH1F *h1_num_lowPt_met_;

    TH1F *h1_num_highPt_mt_;
    TH1F *h1_num_highPt_mlj_;
    TH1F *h1_num_highPt_etaj_;

    TH1F *h1_num_highPt_emfcentralj_;
    TH1F *h1_num_highPt_nemfcentralj_;
    TH1F *h1_num_highPt_cemfcentralj_;

    TH1F *h1_num_highPt_mt_lowmet_;
    TH1F *h1_num_highPt_met_;

    TH1F *h1_num_highPt_mll_;
    TH1F *h1_num_highPt_mllss_;

    TH1F *h1_num_lowPt_mll_;
    TH1F *h1_num_lowPt_mllss_;

};

#endif

