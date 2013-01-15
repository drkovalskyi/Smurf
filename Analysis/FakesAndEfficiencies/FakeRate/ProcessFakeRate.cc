
#include "ProcessFakeRate.h"

#include "TCanvas.h"
#include "THStack.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"

#include <algorithm>
#include <iostream>
#include <cmath>

void printline(FILE* fout, TH2F* h2, const char* caption, bool doInt)
{

    fprintf(fout, "\\begin{table}[!ht]\n");
    fprintf(fout, "{\\small\n");
    fprintf(fout, "\\begin{center}\n");
    fprintf(fout, "\\begin{tabular}{c");
    for (int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y) 
        fprintf(fout, "|c");
    fprintf(fout, "}\n");
    fprintf(fout, "\\hline\n");

    fprintf(fout, "bin &");
    for (int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
    {
        Float_t min = h2->GetYaxis()->GetBinLowEdge(y);
        Float_t max = min + h2->GetYaxis()->GetBinWidth(y);
        if (y == h2->GetYaxis()->GetNbins())
            fprintf(fout, "  %4.2f - %4.2f \\\\ ", min, max);
        else
            fprintf(fout, "  %4.2f - %4.2f  & ", min, max);
    }
    fprintf(fout, "\n\\hline \n");

    for (int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {

        Float_t min = h2->GetXaxis()->GetBinLowEdge(x);
        Float_t max = min + h2->GetXaxis()->GetBinWidth(x);
        fprintf(fout, "  %4.1f - %4.1f  & ", min, max);

        for (int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
        {
            Float_t eff = h2->GetBinContent(x, y);
            Float_t err = h2->GetBinError(x, y);
            if (y == h2->GetYaxis()->GetNbins()) {
                if (doInt)  fprintf(fout, "\t%u $\\pm$ %4.2f \\\\", (unsigned int)eff, err);
                else        fprintf(fout, "\t%4.3f $\\pm$ %4.3f \\\\", eff, err);
            }
            else {
                if (doInt)  fprintf(fout, "\t%u $\\pm$ %4.2f & ", (unsigned int)eff, err);
                else        fprintf(fout, "\t%4.3f $\\pm$ %4.3f & ", eff, err);
            }
        }
        fprintf(fout, "\n");
    }
    fprintf(fout, "\\hline \n");
    fprintf(fout, "\\end{tabular}\n");
    fprintf(fout, Form("\\caption{%s}\n", caption));
    fprintf(fout, "\\end{center}}\n");
    fprintf(fout, "\\end{table}\n");

}

void format(TH1F *hist, DataType dataType)
{

    Color_t color = hist->GetLineColor();

    if (dataType == WJETS) {
        color = kGray+1;
    }
    if (dataType == DY) {
        color = kGreen+2;
    }
    if (dataType == DATA) {
        hist->SetMarkerStyle(20);
    }
    if (dataType == WZ) {
        color = kAzure-9;
    }
    if (dataType == ZZ) {
        color = kAzure+2;
    }

    if (dataType != DATA) {
        hist->SetLineColor(color);
        hist->SetFillColor(color);
    }
}

void compareSS(TFile *f, const char* name, const char* fname)
{

    gStyle->SetOptStat(0);

    TH1F* h1_wjets =   (TH1F*)f->Get(Form("WJets_%s", name));
    TH1F* h1_wjets_ss =   (TH1F*)f->Get(Form("WJets_%sss", name));
    TH1F* h1_data_ss =   (TH1F*)f->Get(Form("Data_%sss", name));
    format(h1_wjets, WJETS);
    format(h1_data_ss, DATA);
    h1_wjets_ss->SetMarkerStyle(22);

    TLegend *l1 = new TLegend(0.6, 0.6, 0.85, 0.85);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_data_ss, "SS Data", "lp");
    l1->AddEntry(h1_wjets, "OS W+jets MC", "f");
    l1->AddEntry(h1_wjets_ss, "SS W+jets MC", "lp");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1_data_ss->Draw();
    h1_wjets->Draw("SAME HIST");
    h1_wjets_ss->Draw("SAME E1");
    l1->Draw();
    c1->SaveAs(Form("plots/%s_%s_ssVal.png", name, fname));
    c1->SaveAs(Form("plots/%s_%s_ssVal.pdf", name, fname));
    c1->SetLogy();
    c1->SaveAs(Form("plots/%s_%s_ssVal_log.png", name, fname));
    c1->SaveAs(Form("plots/%s_%s_ssVal_log.pdf", name, fname));

    delete c1;
    delete h1_wjets;
    delete h1_wjets_ss;
    delete h1_data_ss;

}

void printStack(TFile *f, const char* name, const char* fname, Float_t SF_WJets, Float_t Err_WJets, Float_t SF_DY, Float_t Err_DY, bool ss)
{

    TH1F *h1_wjets = 0;
    if (!ss) h1_wjets =   (TH1F*)f->Get(Form("WJets_%s", name));
    else     h1_wjets =   (TH1F*)f->Get(Form("Data_%sss", name));

    TH1F *h1_dy =       (TH1F*)f->Get(Form("DY_%s", name));
    TH1F *h1_data =     (TH1F*)f->Get(Form("Data_%s", name));
    format(h1_wjets, WJETS);
    format(h1_dy, DY);
    format(h1_data, DATA);
    if (!ss) {
        h1_wjets->Scale(SF_WJets);
        setUncertainty(h1_wjets, Err_WJets);
    }
    h1_dy->Scale(SF_DY);
    setUncertainty(h1_dy, Err_DY);

    TH1F *h1_bg = (TH1F*)h1_wjets->Clone("h1_bg");
    h1_bg->Add(h1_dy);
    h1_bg->SetFillColor(kBlack);
    h1_bg->SetFillStyle(3004);
    THStack st;
    st.Add(h1_wjets);
    st.Add(h1_dy);

    TLegend *l1 = new TLegend(0.6, 0.6, 0.85, 0.85);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_data, "Data", "lp");
    if (!ss) l1->AddEntry(h1_wjets, "W+jets", "f");
    else     l1->AddEntry(h1_wjets, "SS", "f");
    l1->AddEntry(h1_dy, "Drell-Yan", "f");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    st.Draw("HIST");
    if (ss) {
        h1_wjets->SetLineColor(kBlack);
        h1_wjets->Draw("SAME E1");
    }
    h1_bg->Draw("E2 SAME");
    h1_data->Draw("SAME");
    st.SetMinimum(0.1);
    st.SetMaximum(std::max(h1_data->GetMaximum(), st.GetMaximum()));
    st.GetXaxis()->SetTitle(h1_dy->GetXaxis()->GetTitle());
    l1->Draw();
    c1->SaveAs(Form("plots/%s_%s.png", name, fname));
    c1->SaveAs(Form("plots/%s_%s.pdf", name, fname));
    c1->SetLogy(1);
    c1->SaveAs(Form("plots/%s_%s_log.png", name, fname));
    c1->SaveAs(Form("plots/%s_%s_log.pdf", name, fname));

    delete h1_wjets;
    delete h1_dy;
    delete h1_data;
    delete h1_bg;
    delete c1;
    delete l1;

}

void printZStack(TFile *f, const char* name, const char* fname, Float_t SF_EWK, Float_t Err_EWK)
{

    TH1F *h1_wz = 0;
    h1_wz =   (TH1F*)f->Get(Form("WZ_%s", name));
    
    TH1F *h1_zz =       (TH1F*)f->Get(Form("ZZ_%s", name));
    TH1F *h1_data =     (TH1F*)f->Get(Form("Data_%s", name));
    format(h1_wz, WZ);
    format(h1_zz, ZZ);
    format(h1_data, DATA);
    h1_wz->Scale(SF_EWK);
    setUncertainty(h1_wz, Err_EWK);
    h1_zz->Scale(SF_EWK);
    setUncertainty(h1_zz, Err_EWK);

    TH1F *h1_bg = (TH1F*)h1_wz->Clone("h1_bg");
    h1_bg->Add(h1_zz);
    h1_bg->SetFillColor(kBlack);
    h1_bg->SetFillStyle(3004);
    THStack st;
    st.Add(h1_zz);
    st.Add(h1_wz);
    
    TLegend *l1 = new TLegend(0.6, 0.6, 0.85, 0.85);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_data, "Data", "lp");
    l1->AddEntry(h1_wz, "WZ", "f");
    l1->AddEntry(h1_zz, "ZZ", "f");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    st.Draw("HIST");
    h1_bg->Draw("E2 SAME");
    h1_data->Draw("SAME");
    st.SetMinimum(0.1);
    st.SetMaximum(std::max(h1_data->GetMaximum() + 2*sqrt(h1_data->GetMaximum()), st.GetMaximum()));
    st.GetXaxis()->SetTitle(h1_zz->GetXaxis()->GetTitle());
    l1->Draw();
    c1->SaveAs(Form("plots/Z%s_%s.png", name, fname));
    c1->SaveAs(Form("plots/Z%s_%s.pdf", name, fname));
    c1->SetLogy(1);
    c1->SaveAs(Form("plots/Z%s_%s_log.png", name, fname));
    c1->SaveAs(Form("plots/Z%s_%s_log.pdf", name, fname));

    delete h1_wz;
    delete h1_zz;
    delete h1_data;
    delete h1_bg;
    delete c1;
    delete l1;

}

void printValidationHistograms(const char *file, const char *name, const char* bin)
{

    TFile *f = new TFile(Form("histos_FakeLooper_%s.root", file), "READ");
    gROOT->cd();

    //TCanvas *c1 = new TCanvas();
    //c1->cd();

    printStack(f, Form("%s_%s_met", name, bin), file, 1.0, 0.0, 1.0, 0.0);
    printStack(f, Form("%s_%s_mt", name, bin), file, 1.0, 0.0, 1.0, 0.0);

    // get W+Jets scale factor
    float SF_WJets = 0.0;
    float Err_WJets = 0.0;
    getWJetsScaleFactor(f, name, SF_WJets, Err_WJets);

    printStack(f, Form("%s_%s_etaj", name, bin), file, SF_WJets, Err_WJets, 1.0, 0.0);
    printStack(f, Form("%s_%s_mll", name, bin), file, SF_WJets, Err_WJets, 1.0, 0.0, true);

    float SF_DY = 0.0;
    float Err_DY = 0.0;
    getDYScaleFactor(f, name, SF_DY, Err_DY);

    // compare the SS, W+jets and data...
    compareSS(f, Form("%s_%s_mll", name, bin), file);

    f->Close();
    delete f;

}

//printFakeRate("met20mt15mll_"+era, "ElectronFakeRate_V4", "Electron FR with MET20, MT15 and Mll Veto", ptThresholds);
//printValidationHistograms("met20mt15mll_"+era, "ElectronFakeRate_V4", "num_highPt");
void printFakeRate(const char* file, const char *name, const char* friendlyName, const std::vector<unsigned int>& ptThresholds)
{

    // root files
    TFile *out = new TFile(Form("FakeRate_Summary_%s_%s.root", name, file), "RECREATE");
    TFile *f = new TFile(Form("histos_FakeLooper_%s.root", file), "READ");
    gROOT->cd();
    gStyle->SetOptStat(0);

    // tex file
    FILE *fout_tex;
    fout_tex = fopen(Form("FakeRate_Summary_%s_%s.tex", name, file), "w");
    fprintf(fout_tex, "\\documentclass{cmspaper}\n");
    fprintf(fout_tex, "\\usepackage{graphicx}\n");
    fprintf(fout_tex, "\\begin{document}\n");
    fprintf(fout_tex, "\\title{%s}\n", friendlyName);
    fprintf(fout_tex, "\\tableofcontents\n");
    fprintf(fout_tex, "\\clearpage\n");

    Int_t palette[7];
    palette[0] = kWhite;
    for (unsigned int i=1;i<7;++i){
        palette[i] = 18-i;
    }
    gStyle->SetPalette(7,palette);
    gStyle->SetPaintTextFormat("4.3f");
    gStyle->SetMarkerSize(2);

    TCanvas *c1 = new TCanvas();
    c1->SetGridx();
    c1->SetGridy();
    c1->cd();

    // get W+Jets scale factor
    float SF_WJets = 0.0;
    float Err_WJets = 0.0;
    getWJetsScaleFactor(f, name, SF_WJets, Err_WJets);
    float SF_DY = 0.0;
    float Err_DY = 0.0;
    getDYScaleFactor(f, name, SF_DY, Err_DY);

    // print the validation histograms
    fprintf(fout_tex, "\\section{Backgrounds}\n");
    fprintf(fout_tex, "\\subsection{W+Jets Background}\n");
    printStack(f, Form("%s_num_highPt_mt", name), file, SF_WJets, Err_WJets, SF_DY, Err_DY);
    fprintf(fout_tex, "\\begin{figure}[!hbtp]\n");
    fprintf(fout_tex, "\\centering\n");
    fprintf(fout_tex, "\\includegraphics[width=.8\\textwidth]{plots/%s_num_highPt_mt_%s.pdf}\n", name, file);
    fprintf(fout_tex, "\\caption{W+Jets control region where FO passes the numerator and MET $>$ 30 GeV. ");
    fprintf(fout_tex, "The default jet pt threshold is applied, and the DY background is subtracted from simulation. ");
    fprintf(fout_tex, "The scale factor is derived in the region $60<\\rm{MT}<90$ GeV to be %4.2f $\\pm$ %4.2f. ", SF_WJets, Err_WJets);
    fprintf(fout_tex, "In the figure, the MC is normalised to the measured scale factors.}\n");
    fprintf(fout_tex, "\\end{figure}\n");

    // print the validation histograms
    fprintf(fout_tex, "\\subsection{DY Background}\n");

    printStack(f, Form("%s_num_highPt_etaj", name), file, SF_WJets, Err_WJets, SF_DY, Err_DY);
    printStack(f, Form("%s_num_highPt_mll", name), file, SF_WJets, Err_WJets, SF_DY, Err_DY, true);
    printStack(f, Form("%s_num_lowPt_mll", name), file, SF_WJets, Err_WJets, SF_DY, Err_DY, true);
    compareSS(f, Form("%s_num_highPt_mll", name), file);

    fprintf(fout_tex, "\\begin{figure}[!hbtp]\n");
    fprintf(fout_tex, "\\centering\n");
    fprintf(fout_tex, "\\includegraphics[width=.8\\textwidth]{plots/%s_num_highPt_etaj_%s.pdf}\n", name, file);
    fprintf(fout_tex, "\\caption{Eta distribution of the leading away jet. ");
    fprintf(fout_tex, "When deriving the electron fake rate, the region where the jet is beyond 2.5 is vetoed}\n");
    fprintf(fout_tex, "\\end{figure}\n");

    fprintf(fout_tex, "\\begin{figure}[!hbtp]\n");
    fprintf(fout_tex, "\\centering\n");
    fprintf(fout_tex, "\\includegraphics[width=.8\\textwidth]{plots/%s_num_highPt_mll_%s_log.pdf}\n", name, file);
    fprintf(fout_tex, "\\caption{Z+Jets control region where FO passes the numerator and MET $<$ 20 GeV and MT $<$ 15 GeV and $76<\\rm{Mll}<106$ GeV. ");
    fprintf(fout_tex, "The default jet pt threshold is applied, and the same sign yield is used to subtract the QCD background. ");
    fprintf(fout_tex, "The scale factor is derived to be %4.2f $\\pm$ %4.2f. ", SF_DY, Err_DY);
    fprintf(fout_tex, "In the figure, the MC is normalised to the measured scale factors.}\n");
    fprintf(fout_tex, "\\end{figure}\n");

    fprintf(fout_tex, "\\begin{figure}[!hbtp]\n");
    fprintf(fout_tex, "\\centering\n");
    fprintf(fout_tex, "\\includegraphics[width=.8\\textwidth]{plots/%s_num_highPt_mll_%s_ssVal.pdf}\n", name, file);
    fprintf(fout_tex, "\\caption{Cross check of the same sign yield in data, and the W+jets same and opposite sign yields in MC.}\n");
    fprintf(fout_tex, "\\end{figure}\n");

    fprintf(fout_tex, "\\section{Fake Rates}\n");

    for (unsigned int i = 0; i < ptThresholds.size(); ++i) 
        //    for (unsigned int i = 0; i < 1; ++i)

    {

        fprintf(fout_tex, "\\clearpage\n");
        fprintf(fout_tex, "\\subsection{Jet pT threshold %u} \n", ptThresholds[i]);

        // get numerator and denominator
        TH2F* h2_den_Data   = (TH2F*)f->Get(Form("den_Data_%s_ptThreshold%u_PtEta", name, ptThresholds[i]));
        TH2F* h2_num_Data   = (TH2F*)f->Get(Form("num_Data_%s_ptThreshold%u_PtEta", name, ptThresholds[i]));
        TH2F* h2_den_WJets  = (TH2F*)f->Get(Form("den_WJets_%s_ptThreshold%u_PtEta", name, ptThresholds[i]));
        TH2F* h2_num_WJets  = (TH2F*)f->Get(Form("num_WJets_%s_ptThreshold%u_PtEta", name, ptThresholds[i]));
        TH2F* h2_den_DY     = (TH2F*)f->Get(Form("den_DY_%s_ptThreshold%u_PtEta", name, ptThresholds[i]));
        TH2F* h2_num_DY     = (TH2F*)f->Get(Form("num_DY_%s_ptThreshold%u_PtEta", name, ptThresholds[i]));

        // scale MC processes
        h2_den_WJets->Scale(SF_WJets);
        h2_num_WJets->Scale(SF_WJets);
        setUncertainty(h2_den_WJets, Err_WJets);
        setUncertainty(h2_num_WJets, Err_WJets);
        h2_den_DY->Scale(SF_DY);
        h2_num_DY->Scale(SF_DY);
        setUncertainty(h2_den_DY, Err_DY);
        setUncertainty(h2_num_DY, Err_DY);

        h2_num_WJets->Draw("TEXT E1");
        c1->SaveAs(Form("%s_num_WJets.png", file));
        h2_num_Data->Draw("TEXT E1");
        c1->SaveAs(Form("%s_num_Data.png", file));
        h2_num_DY->Draw("TEXT E1");
        c1->SaveAs(Form("%s_num_DY.png", file));

        //
        // number of events in N and D in data
        //

        fprintf(fout_tex, "\\subsubsection{Event Yields in data} \n");

        h2_num_Data->SetTitle(Form("%s_ptThreshold%u_PtEta NData (Numer)", name, ptThresholds[i]));
        h2_num_Data->Draw("TEXT");
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_NData_N.png", name, ptThresholds[i], file));
        printline(fout_tex, h2_num_Data, Form("Summary of numerator counts for jet pT $>$ %u", ptThresholds[i]), true);

        h2_den_Data->SetTitle(Form("%s_ptThreshold%u_PtEta NData (Denom)", name, ptThresholds[i]));
        h2_den_Data->Draw("TEXT");
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_NData_D.png", name, ptThresholds[i], file));
        printline(fout_tex, h2_den_Data, Form("Summary of denominator counts for jet pT $>$ %u", ptThresholds[i]), true);

        //
        // compare data to WJets and DY
        //

        fprintf(fout_tex, "\\clearpage\n");
        fprintf(fout_tex, "\\subsubsection{Estimated background fractions} \n");

        TH2F* h2_den_rel_WJets = (TH2F*)h2_den_WJets->Clone("h2_den_rel_WJets");
        h2_den_rel_WJets->SetTitle(Form("%s_ptThreshold%u_PtEta fWJets(Denom)", name, ptThresholds[i]));
        h2_den_rel_WJets->Divide(h2_den_Data);
        h2_den_rel_WJets->Draw("COLZ TEXT E1");
        h2_den_rel_WJets->SetMaximum(1.0);
        h2_den_rel_WJets->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_fWJets_D.png", name, ptThresholds[i], file));

        TH2F* h2_num_rel_WJets = (TH2F*)h2_num_WJets->Clone("h2_num_rel_WJets");
        h2_num_rel_WJets->SetTitle(Form("%s_ptThreshold%u_PtEta fWJets(Numer)", name, ptThresholds[i]));
        h2_num_rel_WJets->Divide(h2_num_Data);
        h2_num_rel_WJets->Draw("COLZ TEXT E1");
        h2_num_rel_WJets->SetMaximum(1.0);
        h2_num_rel_WJets->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_fWJets_N.png", name, ptThresholds[i], file));

        TH2F* h2_den_rel_DY = (TH2F*)h2_den_DY->Clone("h2_den_rel_DY");
        h2_den_rel_DY->SetTitle(Form("%s_ptThreshold%u_PtEta fDY(Denom)", name, ptThresholds[i]));
        h2_den_rel_DY->Divide(h2_den_Data);
        h2_den_rel_DY->Draw("COLZ TEXT E1");
        h2_den_rel_DY->SetMaximum(1.0);
        h2_den_rel_DY->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_fDY_D.png", name, ptThresholds[i], file));

        TH2F* h2_num_rel_DY = (TH2F*)h2_num_DY->Clone("h2_num_rel_DY");
        h2_num_rel_DY->SetTitle(Form("%s_ptThreshold%u_PtEta fDY(Numer)", name, ptThresholds[i]));
        h2_num_rel_DY->Divide(h2_num_Data);
        h2_num_rel_DY->Draw("COLZ TEXT E1");
        h2_num_rel_DY->SetMaximum(1.0);
        h2_num_rel_DY->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_fDY_N.png", name, ptThresholds[i], file));

        printline(fout_tex, h2_den_rel_WJets, Form("n(W+Jets) / n(Data in denominator) for jet pT $>$ %u", ptThresholds[i]), false);
        printline(fout_tex, h2_den_rel_DY, Form("n(DY) / n(Data in denominator) for jet pT $>$ %u", ptThresholds[i]), false);
        printline(fout_tex, h2_num_rel_WJets, Form("n(W+Jets) / n(Data in numerator) for jet pT $>$ %u", ptThresholds[i]), false);
        printline(fout_tex, h2_num_rel_DY, Form("n(DY) / n(Data in numerator) for jet pT $>$ %u", ptThresholds[i]), false);

        // projection of the sum of EWK contamination
        TH2F *h2_num_EWK = (TH2F*)h2_num_WJets->Clone("h2_num_EWK");
        h2_num_EWK->Add(h2_num_DY);
        TH1F *h1_num_EWK_projection = (TH1F*)h2_num_EWK->ProjectionX("h1_num_EWK_projection");
        TH1F *h1_num_Data_projection = (TH1F*)h2_num_Data->ProjectionX("h1_num_Data_projection");
        TH1F *h1_num_fEWK = (TH1F*)h1_num_EWK_projection->Clone("h1_num_fEWK");
        h1_num_fEWK->Divide(h1_num_Data_projection);
        h1_num_fEWK->Draw();
        h1_num_fEWK->GetYaxis()->SetRangeUser(0, 1.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_fEWK_N.png", name, ptThresholds[i], file));        

        TH2F *h2_den_EWK = (TH2F*)h2_den_WJets->Clone("h2_den_EWK");
        h2_den_EWK->Add(h2_den_DY);
        TH1F *h1_den_EWK_projection = (TH1F*)h2_den_EWK->ProjectionX("h1_den_EWK_projection");
        TH1F *h1_den_Data_projection = (TH1F*)h2_den_Data->ProjectionX("h1_den_Data_projection");
        TH1F *h1_den_fEWK = (TH1F*)h1_den_EWK_projection->Clone("h1_den_fEWK");
        h1_den_fEWK->Divide(h1_den_Data_projection);
        h1_den_fEWK->Draw();
        h1_den_fEWK->GetYaxis()->SetRangeUser(0, 1.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_%s_fEWK_D.png", name, ptThresholds[i], file));  
        //
        // make FR
        //

        // projection of corrected and not corrected fake rate
        TH1F *h1_num_Data = (TH1F*)h2_num_Data->ProjectionX("h1_num_Data");
        TH1F *h1_den_Data = (TH1F*)h2_den_Data->ProjectionX("h1_den_Data");
        TH1F *h1_FR = (TH1F*)h1_num_Data->Clone("h1_FR");
        h1_FR->Divide(h1_den_Data);
        h1_FR->SetLineColor(kBlue);

        TH1F *h1_num_Data_bgsub = (TH1F*)h1_num_Data->Clone("h1_num_Data_bgsub");
        h1_num_Data_bgsub->Add(h1_num_EWK_projection, -1.0);
        TH1F *h1_den_Data_bgsub = (TH1F*)h1_den_Data->Clone("h1_den_Data_bgsub");
        h1_den_Data_bgsub->Add(h1_den_EWK_projection, -1.0);
        TH1F *h1_FR_bgsub = (TH1F*)h1_num_Data_bgsub->Clone("h1_FR_bgsub");
        h1_FR_bgsub->Divide(h1_den_Data_bgsub);
        h1_FR_bgsub->SetLineColor(kRed);

        h1_num_Data->Draw();
        h1_num_Data_bgsub->Draw("SAME HIST");
        c1->SaveAs(Form("plots/%s_ptThreshold%u_TEST_%s.png", name, ptThresholds[i], file));

        h1_FR->Draw();
        h1_FR_bgsub->Draw("SAME");
        h1_FR->GetYaxis()->SetRangeUser(0.0, 1.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_FRProjectionX_%s.png", name, ptThresholds[i], file));

        // make FR with no subtraction
        TH2F* h2_FR = (TH2F*)h2_num_Data->Clone(Form("%s_ptThreshold%u_PtEta_raw", name, ptThresholds[i]));
        h2_FR->SetTitle(Form("%s_ptThreshold%u_PtEta RAW", name, ptThresholds[i]));
        h2_FR->Divide(h2_den_Data);
        h2_FR->Draw("COLZ TEXT E1");
        h2_FR->SetMaximum(0.5);
        h2_FR->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_raw_%s.png",name,  ptThresholds[i], file));

        // subtract MC from data numerator and denominator
        h2_den_Data->Add(h2_den_WJets, -1.0);
        h2_num_Data->Add(h2_num_WJets, -1.0);
        h2_den_Data->Add(h2_den_DY, -1.0);
        h2_num_Data->Add(h2_num_DY, -1.0);
        TH2F* h2_FR_bgsub = (TH2F*)h2_num_Data->Clone(Form("%s_ptThreshold%u_PtEta", name, ptThresholds[i]));
        h2_FR_bgsub->SetTitle(Form("%s_ptThreshold%u_PtEta Corrected", name, ptThresholds[i]));
        h2_FR_bgsub->Divide(h2_den_Data);
        h2_FR_bgsub->Draw("COLZ TEXT E1");
        h2_FR_bgsub->SetMaximum(0.5);
        h2_FR_bgsub->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_cor_%s.png", name, ptThresholds[i], file));

        // ratio of fake rate with and without
        // bg subtraction
        TH2F* h2_FR_ratio = (TH2F*)h2_FR_bgsub->Clone(Form("%s_ptThreshold%u_PtEta_ratio", name, ptThresholds[i]));
        h2_FR_ratio->SetTitle(Form("%s_ptThreshold%u_PtEta RAW:Corrected", name, ptThresholds[i]));
        h2_FR_ratio->Divide(h2_FR);
        h2_FR_ratio->Draw("COLZ TEXT E1");
        h2_FR_ratio->SetMaximum(1.0);
        h2_FR_ratio->SetMinimum(0.0);
        c1->SaveAs(Form("plots/%s_ptThreshold%u_PtEta_ratio_%s.png", name, ptThresholds[i], file));

        //fprintf(fout_tex, "\\clearpage\n");
        fprintf(fout_tex, "\\subsubsection{Fake Rate} \n");
        printline(fout_tex, h2_FR, Form("Fake rate before background subtraction for jet pT $>$ %u", ptThresholds[i]), false);
        printline(fout_tex, h2_FR_bgsub, Form("Fake rate after background subtraction for jet pT $>$ %u", ptThresholds[i]), false);


        out->cd();
        h2_FR_bgsub->Write();
        h2_FR->Write();

        delete h2_den_Data;
        delete h2_num_Data;
        delete h2_den_WJets;
        delete h2_num_WJets;
        delete h2_den_DY;
        delete h2_num_DY;
        delete h2_FR;
        delete h2_FR_bgsub;
        delete h2_FR_ratio;

    }

    fprintf(fout_tex, "\\end{document}\n");
    fclose(fout_tex);
    gROOT->ProcessLine(Form(".! pdflatex FakeRate_Summary_%s_%s.tex", name, file));
    gROOT->ProcessLine(Form(".! pdflatex FakeRate_Summary_%s_%s.tex", name, file));

    delete c1;
    out->Close();
    f->Close();
    delete f;
    delete out;

}

void printZFakeRate(const char* file, const char *name, const char* friendlyName)
{

    // root files
    TFile *out = new TFile(Form("ZFakeRate_Summary_%s_%s.root", name, file), "RECREATE");
    TFile *f = new TFile(Form("histos_ZFakeLooper_%s.root", file), "READ");
    gROOT->cd();
    gStyle->SetOptStat(0);

    // tex file
    FILE *fout_tex;
    fout_tex = fopen(Form("ZFakeRate_Summary_%s_%s.tex", name, file), "w");
    fprintf(fout_tex, "\\documentclass{cmspaper}\n");
    fprintf(fout_tex, "\\usepackage{graphicx}\n");
    fprintf(fout_tex, "\\begin{document}\n");
    fprintf(fout_tex, "\\title{%s}\n", friendlyName);
    fprintf(fout_tex, "\\tableofcontents\n");
    fprintf(fout_tex, "\\clearpage\n");

    Int_t palette[7];
    palette[0] = kWhite;
    for (unsigned int i=1;i<7;++i){
        palette[i] = 18-i;
    }
    gStyle->SetPalette(7,palette);
    gStyle->SetPaintTextFormat("4.3f");
    gStyle->SetMarkerSize(2);

    TCanvas *c1 = new TCanvas();
    c1->SetGridx();
    c1->SetGridy();
    c1->cd();

    // set scale factor and error on EWK background
    float SF_EWK = 1.0;
    float Err_EWK = 0.1;

    // print the validation histograms
    fprintf(fout_tex, "\\section{Backgrounds}\n");
    fprintf(fout_tex, "\\subsection{WZ/ZZ Background}\n");
    printZStack(f, Form("%s_h1_met", name), file, SF_EWK, Err_EWK);
    fprintf(fout_tex, "\\begin{figure}[!hbtp]\n");
    fprintf(fout_tex, "\\centering\n");
    fprintf(fout_tex, "\\includegraphics[width=.8\\textwidth]{plots/Z%s_h1_met_%s.pdf}\n", name, file);
    fprintf(fout_tex, "\\caption{The MET distribution where FO passes the numberator.");
    fprintf(fout_tex, "The systematic uncertainty on each of the ZZ and ZZ backgrounds is 10 percent.");
    fprintf(fout_tex, "The default cut to derive the fake rate is MET $<30$~GeV.}\n");
    fprintf(fout_tex, "\\end{figure}\n");

    fprintf(fout_tex, "\\section{Fake Rates}\n");
    // get numerator and denominator
    //unsigned int ptThreshold = 0;
    unsigned int ptThreshold = 0;
    TH2F* h2_den_Data   = (TH2F*)f->Get(Form("den_Data_%s_ptThreshold%u_PtEta", name, ptThreshold));
    TH2F* h2_num_Data   = (TH2F*)f->Get(Form("num_Data_%s_ptThreshold%u_PtEta", name, ptThreshold));
    TH2F* h2_den_ZZ  = (TH2F*)f->Get(Form("den_ZZ_%s_ptThreshold%u_PtEta", name, ptThreshold));
    TH2F* h2_num_ZZ  = (TH2F*)f->Get(Form("num_ZZ_%s_ptThreshold%u_PtEta", name, ptThreshold));
    TH2F* h2_den_WZ     = (TH2F*)f->Get(Form("den_WZ_%s_ptThreshold%u_PtEta", name, ptThreshold));
    TH2F* h2_num_WZ     = (TH2F*)f->Get(Form("num_WZ_%s_ptThreshold%u_PtEta", name, ptThreshold));

    // scale MC processes
    h2_den_ZZ->Scale(SF_EWK);
    h2_num_ZZ->Scale(SF_EWK);
    setUncertainty(h2_den_ZZ, Err_EWK);
    setUncertainty(h2_num_ZZ, Err_EWK);
    h2_den_WZ->Scale(SF_EWK);
    h2_num_WZ->Scale(SF_EWK);
    setUncertainty(h2_den_WZ, Err_EWK);
    setUncertainty(h2_num_WZ, Err_EWK);

    h2_num_ZZ->Draw("TEXT E1");
    c1->SaveAs(Form("Z%s_num_ZZ.png", file));
    h2_num_Data->Draw("TEXT E1");
    c1->SaveAs(Form("Z%s_num_Data.png", file));
    h2_num_WZ->Draw("TEXT E1");
    c1->SaveAs(Form("Z%s_num_WZ.png", file));

        //
        // number of events in N and D in data
        //

        fprintf(fout_tex, "\\subsubsection{Event Yields in data} \n");

        h2_num_Data->SetTitle(Form("Z%s_ptThreshold%u_PtEta NData (Numer)", name, ptThreshold));
        h2_num_Data->Draw("TEXT");
        c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_NData_N.png", name, ptThreshold, file));
        printline(fout_tex, h2_num_Data, Form("Summary of numerator counts for Z pT $>$ %u", ptThreshold), true);

        h2_den_Data->SetTitle(Form("Z%s_ptThreshold%u_PtEta NData (Denom)", name, ptThreshold));
        h2_den_Data->Draw("TEXT");
        c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_NData_D.png", name, ptThreshold, file));
        printline(fout_tex, h2_den_Data, Form("Summary of denominator counts for Z pT $>$ %u", ptThreshold), true);

    //
    // compare data to ZZ and WZ
    //

    fprintf(fout_tex, "\\clearpage\n");
    fprintf(fout_tex, "\\subsubsection{Estimated background fractions} \n");

    TH2F* h2_den_rel_ZZ = (TH2F*)h2_den_ZZ->Clone("h2_den_rel_ZZ");
    h2_den_rel_ZZ->SetTitle(Form("%s_ptThreshold%u_PtEta fZZ(Denom)", name, ptThreshold));
    h2_den_rel_ZZ->Divide(h2_den_Data);
    h2_den_rel_ZZ->Draw("COLZ TEXT E1");
    h2_den_rel_ZZ->SetMaximum(1.0);
    h2_den_rel_ZZ->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_fZZ_D.png", name, ptThreshold, file));

    TH2F* h2_num_rel_ZZ = (TH2F*)h2_num_ZZ->Clone("h2_num_rel_ZZ");
    h2_num_rel_ZZ->SetTitle(Form("%s_ptThreshold%u_PtEta fZZ(Numer)", name, ptThreshold));
    h2_num_rel_ZZ->Divide(h2_num_Data);
    h2_num_rel_ZZ->Draw("COLZ TEXT E1");
    h2_num_rel_ZZ->SetMaximum(1.0);
    h2_num_rel_ZZ->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_fZZ_N.png", name, ptThreshold, file));

    TH2F* h2_den_rel_WZ = (TH2F*)h2_den_WZ->Clone("h2_den_rel_WZ");
    h2_den_rel_WZ->SetTitle(Form("%s_ptThreshold%u_PtEta fWZ(Denom)", name, ptThreshold));
    h2_den_rel_WZ->Divide(h2_den_Data);
    h2_den_rel_WZ->Draw("COLZ TEXT E1");
    h2_den_rel_WZ->SetMaximum(1.0);
    h2_den_rel_WZ->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_fWZ_D.png", name, ptThreshold, file));

    TH2F* h2_num_rel_WZ = (TH2F*)h2_num_WZ->Clone("h2_num_rel_WZ");
    h2_num_rel_WZ->SetTitle(Form("%s_ptThreshold%u_PtEta fWZ(Numer)", name, ptThreshold));
    h2_num_rel_WZ->Divide(h2_num_Data);
    h2_num_rel_WZ->Draw("COLZ TEXT E1");
    h2_num_rel_WZ->SetMaximum(1.0);
    h2_num_rel_WZ->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_fWZ_N.png", name, ptThreshold, file));

    printline(fout_tex, h2_den_rel_ZZ, Form("n(ZZ) / n(Data in denominator) for jet pT $>$ %u", ptThreshold), false);
    printline(fout_tex, h2_den_rel_WZ, Form("n(WZ) / n(Data in denominator) for jet pT $>$ %u", ptThreshold), false);
    printline(fout_tex, h2_num_rel_ZZ, Form("n(ZZ) / n(Data in numerator) for jet pT $>$ %u", ptThreshold), false);
    printline(fout_tex, h2_num_rel_WZ, Form("n(WZ) / n(Data in numerator) for jet pT $>$ %u", ptThreshold), false);

    // projection of the sum of EWK contamination
    TH2F *h2_num_EWK = (TH2F*)h2_num_ZZ->Clone("h2_num_EWK");
    h2_num_EWK->Add(h2_num_WZ);
    TH1F *h1_num_EWK_projection = (TH1F*)h2_num_EWK->ProjectionX("h1_num_EWK_projection");
    TH1F *h1_num_Data_projection = (TH1F*)h2_num_Data->ProjectionX("h1_num_Data_projection");
    TH1F *h1_num_fEWK = (TH1F*)h1_num_EWK_projection->Clone("h1_num_fEWK");
    h1_num_fEWK->Divide(h1_num_Data_projection);
    h1_num_fEWK->Draw();
    h1_num_fEWK->GetYaxis()->SetRangeUser(0, 1.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_fEWK_N.png", name, ptThreshold, file));        

    TH2F *h2_den_EWK = (TH2F*)h2_den_ZZ->Clone("h2_den_EWK");
    h2_den_EWK->Add(h2_den_WZ);
    TH1F *h1_den_EWK_projection = (TH1F*)h2_den_EWK->ProjectionX("h1_den_EWK_projection");
    TH1F *h1_den_Data_projection = (TH1F*)h2_den_Data->ProjectionX("h1_den_Data_projection");
    TH1F *h1_den_fEWK = (TH1F*)h1_den_EWK_projection->Clone("h1_den_fEWK");
    h1_den_fEWK->Divide(h1_den_Data_projection);
    h1_den_fEWK->Draw();
    h1_den_fEWK->GetYaxis()->SetRangeUser(0, 1.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_%s_fEWK_D.png", name, ptThreshold, file));  

    //
    // make FR
    //

    // projection of corrected and not corrected fake rate
    TH1F *h1_num_Data = (TH1F*)h2_num_Data->ProjectionX("h1_num_Data");
    TH1F *h1_den_Data = (TH1F*)h2_den_Data->ProjectionX("h1_den_Data");
    TH1F *h1_FR = (TH1F*)h1_num_Data->Clone("h1_FR");
    h1_FR->Divide(h1_den_Data);
    h1_FR->SetLineColor(kBlue);

    TH1F *h1_num_Data_bgsub = (TH1F*)h1_num_Data->Clone("h1_num_Data_bgsub");
    h1_num_Data_bgsub->Add(h1_num_EWK_projection, -1.0);
    TH1F *h1_den_Data_bgsub = (TH1F*)h1_den_Data->Clone("h1_den_Data_bgsub");
    h1_den_Data_bgsub->Add(h1_den_EWK_projection, -1.0);
    TH1F *h1_FR_bgsub = (TH1F*)h1_num_Data_bgsub->Clone("h1_FR_bgsub");
    h1_FR_bgsub->Divide(h1_den_Data_bgsub);
    h1_FR_bgsub->SetLineColor(kRed);

    h1_num_Data->Draw();
    h1_num_Data_bgsub->Draw("SAME HIST");
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_TEST_%s.png", name, ptThreshold, file));

    h1_FR->Draw();
    h1_FR_bgsub->Draw("SAME");
    h1_FR->GetYaxis()->SetRangeUser(0.0, 1.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_FRProjectionX_%s.png", name, ptThreshold, file));

    // make FR with no subtraction
    TH2F* h2_FR = (TH2F*)h2_num_Data->Clone(Form("%s_ptThreshold%u_PtEta_raw", name, ptThreshold));
    h2_FR->SetTitle(Form("%s_ptThreshold%u_PtEta RAW", name, ptThreshold));
    h2_FR->Divide(h2_den_Data);
    h2_FR->Draw("COLZ TEXT E1");
    h2_FR->SetMaximum(0.5);
    h2_FR->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_raw_%s.png",name,  ptThreshold, file));

    // subtract MC from data numerator and denominator
    h2_den_Data->Add(h2_den_WZ, -1.0);
    h2_num_Data->Add(h2_num_WZ, -1.0);
    h2_den_Data->Add(h2_den_ZZ, -1.0);
    h2_num_Data->Add(h2_num_ZZ, -1.0);
    TH2F* h2_FR_bgsub = (TH2F*)h2_num_Data->Clone(Form("%s_ptThreshold%u_PtEta", name, ptThreshold));
    h2_FR_bgsub->SetTitle(Form("%s_ptThreshold%u_PtEta Corrected", name, ptThreshold));
    h2_FR_bgsub->Divide(h2_den_Data);
    h2_FR_bgsub->Draw("COLZ TEXT E1");
    h2_FR_bgsub->SetMaximum(0.5);
    h2_FR_bgsub->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_cor_%s.png", name, ptThreshold, file));

    // ratio of fake rate with and without
    // bg subtraction
    TH2F* h2_FR_ratio = (TH2F*)h2_FR_bgsub->Clone(Form("%s_ptThreshold%u_PtEta_ratio", name, ptThreshold));
    h2_FR_ratio->SetTitle(Form("%s_ptThreshold%u_PtEta RAW:Corrected", name, ptThreshold));
    h2_FR_ratio->Divide(h2_FR);
    h2_FR_ratio->Draw("COLZ TEXT E1");
    h2_FR_ratio->SetMaximum(1.0);
    h2_FR_ratio->SetMinimum(0.0);
    c1->SaveAs(Form("plots/Z%s_ptThreshold%u_PtEta_ratio_%s.png", name, ptThreshold, file));

    //fprintf(fout_tex, "\\clearpage\n");
    fprintf(fout_tex, "\\subsubsection{Fake Rate} \n");
    printline(fout_tex, h2_FR, Form("Fake rate before background subtraction for Z pT $>$ %u", ptThreshold), false);
    printline(fout_tex, h2_FR_bgsub, Form("Fake rate after background subtraction for Z pT $>$ %u", ptThreshold), false);


    out->cd();
    h2_FR_bgsub->Write();
    h2_FR->Write();

    delete h2_den_Data;
    delete h2_num_Data;
    delete h2_den_WZ;
    delete h2_num_WZ;
    delete h2_den_ZZ;
    delete h2_num_ZZ;
    delete h2_FR;
    delete h2_FR_bgsub;
    delete h2_FR_ratio;

    fprintf(fout_tex, "\\end{document}\n");
    fclose(fout_tex);
    gROOT->ProcessLine(Form(".! pdflatex ZFakeRate_Summary_%s_%s.tex", name, file));
    gROOT->ProcessLine(Form(".! pdflatex ZFakeRate_Summary_%s_%s.tex", name, file));

    delete c1;
    out->Close();
    f->Close();
    delete f;
    delete out;

}

void getWJetsScaleFactor(TFile *f, const char* name, float &sf, float &err)
{

    // get the scale factor for the W+jets
    // and its uncertainty
    TH1F *h1_wjets_ctl_mt =    (TH1F*)f->Get(Form("WJets_%s_num_highPt_mt", name))->Clone("wjets");
    TH1F *h1_dy_ctl_mt =       (TH1F*)f->Get(Form("DY_%s_num_highPt_mt", name))->Clone("dy");
    TH1F *h1_data_ctl_mt =     (TH1F*)f->Get(Form("Data_%s_num_highPt_mt", name))->Clone("data");
    Int_t binLow =  h1_wjets_ctl_mt->GetXaxis()->FindBin(61.0);
    Int_t binHigh  = h1_wjets_ctl_mt->GetXaxis()->FindBin(89.0);
    Float_t nData  = h1_data_ctl_mt->Integral(binLow, binHigh);
    Float_t nDY    = h1_dy_ctl_mt->Integral(binLow, binHigh);
    Float_t nWJets = h1_wjets_ctl_mt->Integral(binLow, binHigh);
    printf("Data  = %u\n", (unsigned int)nData);
    printf("DY    = %4.2f\n", nDY);
    printf("WJets = %4.2f\n", nWJets);
    sf  = (nData - nDY) / nWJets;
    //err = fabs(1-sf) / 2.0;
    err = fabs(1-sf);
    printf("WJets SF = %4.2f +/- %4.2f (syst)\n", sf, err);
    delete h1_wjets_ctl_mt;
    delete h1_dy_ctl_mt;
    delete h1_data_ctl_mt;

}

void getDYScaleFactor(TFile *f, const char* name, float &sf, float &err)
{       

    // get the scale factor for the W+jets
    // and its uncertainty
    TH1F *h1_wjets_ctl_mll =    (TH1F*)f->Get(Form("WJets_%s_num_highPt_mll", name))->Clone("wjets");
    TH1F *h1_dy_ctl_mll =       (TH1F*)f->Get(Form("DY_%s_num_highPt_mll", name))->Clone("dy");
    TH1F *h1_data_ctl_mll =     (TH1F*)f->Get(Form("Data_%s_num_highPt_mll", name))->Clone("data");
    TH1F *h1_data_ctl_mllss =     (TH1F*)f->Get(Form("Data_%s_num_highPt_mllss", name))->Clone("data");

    Int_t binLow =  h1_wjets_ctl_mll->GetXaxis()->FindBin(76.1);
    Int_t binHigh  = h1_wjets_ctl_mll->GetXaxis()->FindBin(105.9);
    Float_t nData  = h1_data_ctl_mll->Integral(binLow, binHigh);
    Float_t nDY    = h1_dy_ctl_mll->Integral(binLow, binHigh);
    Float_t nWJets = h1_wjets_ctl_mll->Integral(binLow, binHigh);
    Float_t nSS    = h1_data_ctl_mllss->Integral(binLow, binHigh);
    printf("Data  = %u\n", (unsigned int)nData);
    printf("DY    = %4.2f\n", nDY);
    printf("WJets = %4.2f\n", nWJets);
    printf("SS = %4.2f\n", nSS);
    sf  = (nData - nSS) / nDY;
    //err = fabs(1-sf) / 2.0;
    err = fabs(1-sf);
    printf("DY SF = %4.2f +/- %4.2f (syst)\n", sf, err);

    delete h1_data_ctl_mllss;
    delete h1_wjets_ctl_mll;
    delete h1_dy_ctl_mll;
    delete h1_data_ctl_mll;

}

void setUncertainty(TH1F *hist, const float &err)
{
    for (Int_t i = 0; i < hist->GetNbinsX()+1; ++i) {
        hist->SetBinError(i, sqrt(pow(hist->GetBinError(i), 2) + pow(hist->GetBinContent(i)*err, 2)) );
    }
}

void setUncertainty(TH2F *hist, const float &err)
{   
    for (Int_t x = 0; x < hist->GetXaxis()->GetNbins()+1; ++x) {
        for (Int_t y = 0; y < hist->GetYaxis()->GetNbins()+1; ++y) {
            hist->SetBinError(x, y, sqrt( pow(hist->GetBinError(x, y), 2) + pow(hist->GetBinContent(x, y)*err, 2)) );
        }
    }
} 

