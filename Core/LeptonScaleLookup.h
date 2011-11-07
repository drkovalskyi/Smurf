#ifndef LEPTONSCALELOOKUP_H
#define LEPTONSCALELOOKUP_H

#include "TFile.h"
#include <iostream>
#include <string>
#include "TH2F.h"

class LeptonScaleLookup {

    public:

        LeptonScaleLookup(std::string filename);
        ~LeptonScaleLookup();

        float GetEfficiency(float eta, float pt, TH2F *hist);
        float GetError(float eta, float pt, TH2F *hist);

        float GetExpectedTriggerEfficiency(float eta1, float pt1, float eta2, float pt2, int id1, int id2);
        float GetExpectedLeptonSF(float eta, float pt, int id);
        float GetExpectedLeptonSFErr(float eta, float pt, int id);
        float GetExpectedLeptonEff(float eta, float pt, int id);

        void printTable(std::string name);

    private:
        TFile *file_;
        TH2F *h2_single_e_;
        TH2F *h2_single_m_;
        TH2F *h2_double_e_LeadingLeg_;
        TH2F *h2_double_e_TrailingLeg_;
        TH2F *h2_double_m_LeadingLeg_;
        TH2F *h2_double_m_TrailingLeg_;
        TH2F *h2_cross_e_LeadingLeg_;
        TH2F *h2_cross_e_TrailingLeg_;
        TH2F *h2_cross_m_LeadingLeg_;
        TH2F *h2_cross_m_TrailingLeg_;

        // efficiencies
        TH2F *h2_selection_e_;
        TH2F *h2_selection_m_;
        TH2F *h2_selection_eff_e_;
        TH2F *h2_selection_eff_m_;
};

void LoopupAll(std::string file);

#endif
