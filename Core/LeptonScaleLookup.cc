
#include "LeptonScaleLookup.h"

LeptonScaleLookup::LeptonScaleLookup(std::string filename)
{
    file_ = new TFile(filename.c_str(), "READ");
    h2_double_e_ = (TH2F*)file_->Get("h2_results_electron_double");
    h2_single_e_ = (TH2F*)file_->Get("h2_results_electron_single");
    h2_double_m_ = (TH2F*)file_->Get("h2_results_muon_double");
    h2_single_m_ = (TH2F*)file_->Get("h2_results_muon_single");
    h2_cross_m_ = (TH2F*)h2_double_m_->Clone("h2_cross_m");
    h2_cross_e_ = (TH2F*)h2_double_e_->Clone("h2_cross_e");

    h2_selection_e_ = (TH2F*)file_->Get("h2_results_electron_selection");
    h2_selection_m_ = (TH2F*)file_->Get("h2_results_muon_selection");

}

LeptonScaleLookup::~LeptonScaleLookup()
{
    delete h2_cross_m_;
    delete h2_cross_e_;
    file_->Close();
    delete file_;
}

float LeptonScaleLookup::GetEfficiency(float eta, float pt, TH2F *hist) 
{
    // make sure pt to look up is in range
    int nbins = hist->GetXaxis()->GetNbins();
    if (pt > (hist->GetXaxis()->GetBinLowEdge(nbins) + hist->GetXaxis()->GetBinWidth(nbins)))
        pt = hist->GetXaxis()->GetBinLowEdge(nbins) + (hist->GetXaxis()->GetBinWidth(nbins)/2.0);

    // make sure eta to look up is in range
    // eta is abs
    eta = fabs(eta);
    nbins = hist->GetYaxis()->GetNbins();
    if (eta > (hist->GetYaxis()->GetBinLowEdge(nbins) + hist->GetYaxis()->GetBinWidth(nbins)))
        eta = hist->GetYaxis()->GetBinLowEdge(nbins) + (hist->GetYaxis()->GetBinWidth(nbins)/2.0);

    // look up the efficiency
    Int_t binX = hist->GetXaxis()->FindBin(pt);
    Int_t binY = hist->GetYaxis()->FindBin(eta);
    return hist->GetBinContent(binX, binY);
}

float LeptonScaleLookup::GetExpectedTriggerEfficiency(float eta1, float pt1, float eta2, float pt2, int id1, int id2)
{

    unsigned int f1 = abs(id1);
    unsigned int f2 = abs(id2);

    float eff_sgl_1, eff_sgl_2;
    float eff_dbl_1, eff_dbl_2;

    // get individual leg efficiencies

    if (f1 == 11 && f2 == 11) {
        eff_sgl_1 = GetEfficiency(eta1, pt1, h2_single_e_);
        eff_sgl_2 = GetEfficiency(eta2, pt2, h2_single_e_);
        eff_dbl_1 = GetEfficiency(eta1, pt1, h2_double_e_);
        eff_dbl_2 = GetEfficiency(eta2, pt2, h2_double_e_);
    }
    else if (f1 == 13 && f2 == 13) {
        eff_sgl_1 = GetEfficiency(eta1, pt1, h2_single_m_);
        eff_sgl_2 = GetEfficiency(eta2, pt2, h2_single_m_);
        eff_dbl_1 = GetEfficiency(eta1, pt1, h2_double_m_);
        eff_dbl_2 = GetEfficiency(eta2, pt2, h2_double_m_);
    }
    else if (f1 == 13 && f2 == 11) {
        eff_sgl_1 = GetEfficiency(eta1, pt1, h2_single_m_);
        eff_sgl_2 = GetEfficiency(eta2, pt2, h2_single_e_);
        eff_dbl_1 = GetEfficiency(eta1, pt1, h2_cross_m_);
        eff_dbl_2 = GetEfficiency(eta2, pt2, h2_cross_e_);
    }
    else if (f1 == 11 && f2 == 13) {
        eff_sgl_1 = GetEfficiency(eta1, pt1, h2_single_e_);
        eff_sgl_2 = GetEfficiency(eta2, pt2, h2_single_m_);
        eff_dbl_1 = GetEfficiency(eta1, pt1, h2_cross_e_);
        eff_dbl_2 = GetEfficiency(eta2, pt2, h2_cross_m_);
    }
    else {
        std::cout << "[LeptonScaleLookup::GetExpectedTriggerEfficiency] ERROR: Invalid flavor combination!" << std::endl;
        return 0.0;
    }

    // calculate event efficiency
    float evt_eff = eff_dbl_1*eff_dbl_2 
                        + eff_sgl_2*(1-eff_dbl_1)
                        + eff_sgl_1*(1-eff_dbl_2);

    // return it
    return evt_eff;        

}

float LeptonScaleLookup::GetExpectedLeptonSF(float eta, float pt, int id)
{

    float sf = 0.0;
    if (abs(id) == 11) {
        sf = GetEfficiency(eta, pt, h2_selection_e_);
    } else if (abs(id) == 13) {
        sf = GetEfficiency(eta, pt, h2_selection_m_);
    } else {
        std::cout << "[LeptonScaleLookup::GetExpectedOfflineSF] ERROR: Invalid flavor!" << std::endl;
    }
    return sf;

}

void LeptonScaleLookup::printTable(std::string name)
{
    TH2F *h2 = (TH2F*)file_->Get(name.c_str());

    printf("\\begin{table}[!ht]\n");
    printf("\\begin{center}\n");

    printf("\\begin{tabular}{c");
    for (int i = 0; i < h2->GetNbinsY(); ++i) printf("|c");
    printf("}\n");

    printf("\\hline\n");

    printf("Measurement ");
    for (Int_t eta = 1; eta < h2->GetYaxis()->GetNbins() + 1; ++eta) {
        Float_t binMin = h2->GetYaxis()->GetBinLowEdge(eta);
        Float_t binMax = binMin + h2->GetYaxis()->GetBinWidth(eta);
        printf(" & $%4.2f<\eta<%4.2f$ ", binMin, binMax);
    }
    printf(" \\\\ \n");

    printf("\\hline\n");
    for (Int_t pt = 1; pt < h2->GetXaxis()->GetNbins() + 1; ++pt) {
        Float_t binMin = h2->GetXaxis()->GetBinLowEdge(pt);
        Float_t binMax = binMin + h2->GetXaxis()->GetBinWidth(pt);
        printf("$%4.0f<p_T<%4.0f$", binMin, binMax);
        for (Int_t eta = 1; eta < h2->GetYaxis()->GetNbins() + 1; ++eta) {
            Float_t sf = h2->GetBinContent(pt, eta);
            Float_t sferr = h2->GetBinError(pt, eta);
            printf(" & %4.2f $\\pm$ %4.2f ", sf, sferr);
        }
        printf(" \\\\ \\hline \n");
    }

    printf("\\end{tabular}\n");
    TString caption = h2->GetTitle();
    caption.ReplaceAll("_", "\\_");
    printf("\\caption{Hello, I'm a table for %s}\n", caption.Data());
    printf("\\label{tab:eff_ele_offline}\n");
    printf("\\end{center}\n");
    printf("\\end{table}\n");

    delete h2;
}

void LoopupAll(std::string file) {

    LeptonScaleLookup lookup(file);

    printf("\\documentclass[article]{revtex4}\n");
    printf("\\usepackage{graphicx}\n");
    printf("\\begin{document}\n");
    printf("\\tableofcontents\n");

    lookup.printTable("h2_results_electron_selection");
    lookup.printTable("h2_results_muon_selection");
    printf("\\clearpage\n");
    lookup.printTable("h2_results_electron_single");
    lookup.printTable("h2_results_electron_double");
    lookup.printTable("h2_results_muon_single");
    lookup.printTable("h2_results_muon_double");

    printf("\\end{document}\n");

}


