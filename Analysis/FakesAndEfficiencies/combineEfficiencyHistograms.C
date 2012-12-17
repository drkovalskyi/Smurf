
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"

#include <vector>
#include <iostream>

TH2D *new_histogram(std::vector<TString> fileNames, TString newName, TString hist);

void combineEfficiencyHistograms()
{

    TString effDir = "/smurf/dlevans/Efficiencies/V00-02-09/";
    TString nameEl = "V00-02-09_Moriond_V0";
    TString nameMu = nameEl;
    TFile outfile(effDir + "/summary_Moriond_V0.root", "RECREATE");
    TDirectory *outdir = (TDirectory*)outfile.GetDirectory("");

    gROOT->cd();
    std::vector<TString> fileNames;

    // muon selection
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_ID/extra/smurf.root");
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_Iso/extra/smurf.root");
    TH2D *h2_results_muon_selection = new_histogram(fileNames, "h2_results_muon_selection", "SF");
    h2_results_muon_selection->SetDirectory(outdir);
   
    // electron selection
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_ID/extra/smurf.root");
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_Iso/extra/smurf.root");
    TH2D *h2_results_electron_selection = new_histogram(fileNames, "h2_results_electron_selection", "SF");
    h2_results_electron_selection->SetDirectory(outdir);

    // triggers

    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigSgl/extra/smurf.root");
    TH2D *h2_results_electron_single = new_histogram(fileNames, "h2_results_electron_single", "Data");
    h2_results_electron_single->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigLeadDbl/extra/smurf.root");
    TH2D *h2_results_electron_double_leadingleg = new_histogram(fileNames, "h2_results_electron_double_leadingleg", "Data");
    h2_results_electron_double_leadingleg->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigTrailDbl/extra/smurf.root");
    TH2D *h2_results_electron_double_trailingleg = new_histogram(fileNames, "h2_results_electron_double_trailingleg", "Data");
    h2_results_electron_double_trailingleg->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigDzDbl/extra/smurf.root");
    TH2D *h2_results_electron_double_dz = new_histogram(fileNames, "h2_results_electron_double_dz", "Data");
    h2_results_electron_double_dz->SetDirectory(outdir);

    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigSgl/extra/smurf.root");
    TH2D *h2_results_muon_single = new_histogram(fileNames, "h2_results_muon_single", "Data");
    h2_results_muon_single->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigLeadDbl/extra/smurf.root");
    TH2D *h2_results_muon_double_leadingleg = new_histogram(fileNames, "h2_results_muon_double_leadingleg", "Data");
    h2_results_muon_double_leadingleg->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigTrailDbl/extra/smurf.root");
    TH2D *h2_results_muon_double_trailingleg = new_histogram(fileNames, "h2_results_muon_double_trailingleg", "Data");
    h2_results_muon_double_trailingleg->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigDzDbl/extra/smurf.root");
    TH2D *h2_results_muon_double_dz = new_histogram(fileNames, "h2_results_muon_double_dz", "Data");
    h2_results_muon_double_dz->SetDirectory(outdir);

    // write output
    outfile.Write();
    outfile.Close();

}

TH2D *new_histogram(std::vector<TString> fileNames, TString newName, TString hist)
{

    std::vector<TFile*> files;
    for (unsigned int i = 0; i < fileNames.size(); ++i)
        files.push_back(new TFile(fileNames[i], "READ"));

    gROOT->cd();
    files[0]->ls();

    TH2D *new_hist = 0;
    if (hist == "Data") 
        new_hist = (TH2D*)files[0]->Get("smurf_"+hist+"Eff")->Clone(newName);
    else 
        new_hist = (TH2D*)files[0]->Get("smurf_"+hist)->Clone(newName);
    
    if (new_hist == 0) {
        std::cout << "ERROR" << std::endl;
        return 0;
    }

    for (int ix = 1; ix <= new_hist->GetXaxis()->GetNbins() + 1; ++ix) {
        for (int iy = 1; iy <= new_hist->GetYaxis()->GetNbins() + 1; ++iy) {

            double sf       = 1.0;
            double err2h   = 0.0;
            double err2l    = 0.0;
            for (unsigned int f = 0; f < files.size(); ++f) {
                TH2D *temp_sf   = 0;
                if (hist == "Data") temp_sf = (TH2D*)files[f]->Get("smurf_"+hist+"Eff");
                else temp_sf = (TH2D*)files[f]->Get("smurf_"+hist);
                TH2D *temp_errh = (TH2D*)files[f]->Get("smurf_"+hist+"ErrHigh");
                TH2D *temp_errl = (TH2D*)files[f]->Get("smurf_"+hist+"ErrLow");
                sf      *= temp_sf->GetBinContent(ix, iy);
                err2h   += pow(temp_errh->GetBinContent(ix, iy)/temp_sf->GetBinContent(ix, iy), 2);
                err2l   += pow(temp_errl->GetBinContent(ix, iy)/temp_sf->GetBinContent(ix, iy), 2);
                delete temp_sf;
                delete temp_errl;
                delete temp_errh;
            }

            new_hist->SetBinContent(ix, iy, sf);
            new_hist->SetBinError(ix, iy, sf * (sqrt(err2h) + sqrt(err2l))/2.0);

        }
    }

    for (unsigned int i = 0; i < fileNames.size(); ++i) delete files[i];
    return new_hist;

}

