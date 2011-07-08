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
        float GetExpectedTriggerEfficiency(float eta1, float pt1, float eta2, float pt2, int id1, int id2);
        float GetExpectedLeptonSF(float eta, float pt, int id);

        void printTable(std::string name);

    private:
        TFile *file_;
        TH2F *h2_double_e_;
        TH2F *h2_single_e_;
        TH2F *h2_double_m_;
        TH2F *h2_single_m_;
        TH2F *h2_cross_m_;
        TH2F *h2_cross_e_;
        TH2F *h2_selection_e_;
        TH2F *h2_selection_m_;
};

void LoopupAll(std::string file);

#endif
