
{

    TFile f_new("FakeRate_Summary_ElectronFakeRate_V4_met20mt15mll_HCP.root", "READ");
    //TFile f_new("FakeRate_Summary_MuonFakeRate_M2_met20mt15mll_HCP.root", "READ");
    TFile f_old("/smurf/dlevans/FakeRates/V00-02-07_HCP_V0/summary.root", "READ");


    gROOT->cd();

    
    Int_t palette[7];
    palette[0] = kWhite;
    for (unsigned int i=1;i<7;++i){
        palette[i] = 18-i;
    }
    gStyle->SetPalette(7,palette);
    gStyle->SetPaintTextFormat("4.2f");
    gStyle->SetMarkerSize(2);
    gStyle->SetOptStat(0);


    TH2F *fr_new = (TH2F*)f_new.Get("ElectronFakeRate_V4_ptThreshold30_PtEta");
    TH2F *fr_old = (TH2F*)f_old.Get("ElectronFakeRate_V4_ptThreshold30_PtEta");
    //TH2F *fr_new = (TH2F*)f_new.Get("MuonFakeRate_M2_ptThreshold15_PtEta");
    //TH2F *fr_old = (TH2F*)f_old.Get("MuonFakeRate_M2_ptThreshold15_PtEta");

    fr_new->Divide(fr_old);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    fr_new->Draw("COLZ TEXT E1");
    fr_new->SetMinimum(0.3);
    fr_new->SetMaximum(1.4);
    c1->SaveAs("electrons_compareToHCP.png");
    //c1->SaveAs("muons_compareToHCP.png");

}

