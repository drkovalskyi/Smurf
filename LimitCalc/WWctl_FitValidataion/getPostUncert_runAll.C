#include <vector>
#include "TMath.h"

const bool DEBUG = true;

//
// Zoom
//
TH2F* Zoom(TH2F* h2) {

    float mTbins[7]    =   {60,70,80,90,100,110,120};
    float mllbins[4]   =   {12,30,45,60};

    TH2F* h2zoom = new TH2F(Form("%s_zoom", h2->GetName()), Form("%s_zoom", h2->GetTitle()), 6, mTbins, 3, mllbins);

    for(unsigned int ibinX=1; ibinX <= 6; ++ibinX) {
        for(unsigned int ibinY=1; ibinY <= 3; ++ibinY) {
            h2zoom->SetBinContent(    ibinX, ibinY, h2->GetBinContent(ibinX, ibinY) );
            h2zoom->SetBinError(  ibinX, ibinY, h2->GetBinError(ibinX, ibinY) );
        }
    }
    h2zoom->SetStats(0);
    h2zoom->SetXTitle("M_{T} (GeV)");
    h2zoom->SetYTitle("M_{ll} (GeV)");
    h2zoom->SetTitleOffset(1.5,"Y");;
    h2zoom->SetMinimum(0);

    return h2zoom;
}

//
// Rolling back to 2D
//
TH2F* Roll1DTo2D(TH1F *h1, bool doZoom=false) {

    // style
    gStyle->SetPaintTextFormat(".1f");
    gStyle->SetMarkerSize(2.5);

    // palette
    //loadColorPalette();

    unsigned int nbins  = h1->GetXaxis()->GetNbins();

    float mTbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
    float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};
    TH2F* h2 = new TH2F(Form("%s_2D", h1->GetName()), Form("%s_2D",h1->GetTitle()), 14, mTbins, 9, mllbins);
    for(unsigned int ibin=1; ibin <= nbins; ++ibin) {

        int ibinX;
        int ibinY;
        if(ibin<=9*6) {
            ibinX = (int) (ibin-1)/9+1;
            ibinY = (int) (ibin-1)%9+1;
        } else {
            ibinX = (int) (ibin-1)/9+1;
            ibinY = (int) (ibin-1)%9+1;
        }

        h2->SetBinContent(  ibinX, ibinY, h1->GetBinContent(ibin) );
        h2->SetBinError(    ibinX, ibinY, h1->GetBinError(ibin) );
    }
    h2->SetStats(0); 
    h2->SetXTitle("M_{T} (GeV)");
    h2->SetYTitle("M_{ll} (GeV)");
    h2->SetTitleOffset(1.5,"Y");;
    h2->SetMinimum(0);

    // do zoom or not
    TH2F *h2zoom = (TH2F*)Zoom(h2);

    if(doZoom) return h2zoom;
    else return h2;
}

//
// Convert TGraphAsymmErrors to TH1F
//
TH1F* makeHistogram(TGraphAsymmErrors* graph){

    TString histoname = graph->GetName();
    histoname.ReplaceAll("_", "");  
    histoname.ReplaceAll("j0of", "");  
    histoname.ReplaceAll("j1of", "");  
    histoname.ReplaceAll("8tev", "");  
    histoname.ReplaceAll("7tev", "");  

    TH1F* hist = new TH1F(Form("histo_new_%s", histoname.Data()), Form("histo_new_%s",histoname.Data()), 126, -1, 1); // Guillelmo's card 
    //TH1F* hist = new TH1F(Form("histo_new_%s", histoname.Data()), Form("histo_new_%s",histoname.Data()), 126, 0.5, 126.5); // Jae's card
    hist->SetDirectory(0);
    for (int i=0; i<graph->GetN();++i){
        Int_t bin = hist->FindBin(graph->GetX()[i]);
        hist->SetBinContent(bin,graph->GetY()[i]);
    }
    hist->SetStats(0);
    return hist;
}


TH1F* doProjection(TH2F* h2, char* var="mll", char* plotregion="CR1") {
    
    TString histoname = h2->GetName();
    histoname.ReplaceAll("_", "");  
    histoname.ReplaceAll("j0of", "");  
    histoname.ReplaceAll("j1of", "");  
    histoname.ReplaceAll("8tev", "");  
    histoname.ReplaceAll("7tev", "");  
    histoname.ReplaceAll("histo", "");  
    histoname.ReplaceAll("new", "");  
    histoname.ReplaceAll("2D", "");  
    
    // 
    // 
    // 
    if(var=="mll") { // ---------------------- mll ---------------------------------

        if(plotregion=="CR1") { // ----------- CR1 --------------------------------
            float mllbins_CR1[10]  =   {12,30,45,60,75,100,125,150,175,200};
            TH1F* h1_mll_CR1 = new TH1F(Form("h1_mll_CR1_%s", histoname.Data()), Form("h1_mll_CR1_%s", histoname.Data()), 9, mllbins_CR1);
            for(int imll=1; imll<h2->GetYaxis()->GetNbins()+1; imll++){ // 1 - 9
                for(int imt=7; imt<h2->GetXaxis()->GetNbins()+1; imt++){ // 7 - 14
                    h1_mll_CR1->SetBinContent(imll, h2->GetBinContent(imt,imll)+h1_mll_CR1->GetBinContent(imll));
                }
            }
            return h1_mll_CR1;
        } else if (plotregion=="CR2") {  // - CR2 --------------------------------
            float mllbins_CR2[7]   =   {60,75,100,125,150,175,200};
            TH1F* h1_mll_CR2 = new TH1F(Form("h1_mll_CR2_%s", histoname.Data()), Form("h1_mll_CR2_%s", histoname.Data()), 6, mllbins_CR2);
            for(int imll=4; imll<h2->GetYaxis()->GetNbins()+1; imll++){ // 4 - 9
                for(int imt=1; imt<7; imt++){ // 1 - 6
                    h1_mll_CR2->SetBinContent(imll-3, h2->GetBinContent(imt,imll)+h1_mll_CR2->GetBinContent(imll-3));
                }
            }
            return h1_mll_CR2;
        } else {
            cout << " wrong control region name " << endl;
        }

    } else if (var=="mT") { // --------------- mT --------------------------------

        if(plotregion=="CR1") { // ----------- CR1 --------------------------------
            float mTbins_CR1[9]    =   {120,140,160,180,200,220,240,260,280};
            TH1F* h1_mT_CR1  = new TH1F(Form("h1_mT_CR1_%s", histoname.Data()), Form("h1_mT_CR1_%s", histoname.Data()), 8, mTbins_CR1);
            for(int imt=7; imt<h2->GetXaxis()->GetNbins()+1; imt++){ // 7 - 14
                for(int imll=1; imll<h2->GetYaxis()->GetNbins()+1; imll++){ // 1 - 9
                    h1_mT_CR1->SetBinContent(imt-6, h2->GetBinContent(imt,imll)+h1_mT_CR1->GetBinContent(imt-6));
                }
            }
            return h1_mT_CR1;
        } else if (plotregion=="CR2") {  // - CR2 --------------------------------
            float mTbins_CR2[7]    =   {60,70,80,90,100,110,120};
            TH1F* h1_mT_CR2  = new TH1F(Form("h1_mT_CR2_%s", histoname.Data()), Form("h1_mT_CR2_%s", histoname.Data()), 6, mTbins_CR2);
            for(int imt=1; imt<7; imt++){ // 1 - 6
                for(int imll=4; imll<h2->GetYaxis()->GetNbins()+1; imll++){ // 4 - 9
                    h1_mT_CR2->SetBinContent(imt, h2->GetBinContent(imt,imll)+h1_mT_CR2->GetBinContent(imt));
                }
            }
            return h1_mT_CR2;
        } else {
            cout << " wrong control region name " << endl;
        }

    } else {
        cout << " wrong variable name " << endl;
    }



}

//
// Cosmetics for projected histograms
//
void h1cosmetic(TH1F* &h1, char* title="", int linecolor=kRed, int linewidth=1, int fillcolor=0, TString var=""){

    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.);

}
//
// reorganize background processes 
//    
//  sig = ZH, WH, qqH, ggH, 
//  WW = qqWW + ggWW
//  VVWgamma = VV + Wgamma + Wg3l
//  Top = Top
//  Wjets = WjetsE + WjetsM 
//
void  reorganizebkg(  TH1F* &h1_ZH, TH1F* &h1_WH, TH1F* &h1_qqH, TH1F* &h1_ggH, 
                      TH1F* &h1_qqWW, TH1F* &h1_ggWW, TH1F* &h1_VV, TH1F* &h1_Top, 
                      TH1F* &h1_Zjets, TH1F* &h1_WjetsE, TH1F* &h1_WjetsM, 
                      TH1F* &h1_Wgamma, TH1F* &h1_Wg3l, TH1F* &h1_Ztt, 
                      TH1F* &h1_sig, TH1F* &h1_WW, TH1F* &h1_VVWgamma, TH1F* &h1_Wjets,
                      TString var="var") {


    // signal 
    h1_sig = (TH1F*)h1_ZH->Clone("h1_sig");
    h1_sig->Add(h1_WH);
    h1_sig->Add(h1_qqH);
    h1_sig->Add(h1_ggH);

    // WW 
    h1_WW = (TH1F*)h1_qqWW->Clone("h1_WW");
    h1_WW->Add(h1_ggWW);  

    // VVWgamma 
    h1_VVWgamma = (TH1F*)h1_VV->Clone("h1_VVWgamma");
    h1_VVWgamma->Add(h1_Wgamma);  
    h1_VVWgamma->Add(h1_Wg3l); 

    // Top = Top

    // Wjets 
    h1_Wjets = (TH1F*)h1_WjetsE->Clone("h1_Wjets");
    h1_Wjets->Add(h1_WjetsM);  

    // Zjets = Ztt

    // Cosmetic  
    h1cosmetic(h1_sig,          "signal",   kBlack, 1, kRed,        var);
    h1cosmetic(h1_WW,           "WW",       kBlack, 1, kAzure-9,    var);
    h1cosmetic(h1_VVWgamma,     "VV",       kBlack, 1, kAzure-2,    var);
    h1cosmetic(h1_Top,          "Top",      kBlack, 1, kYellow,     var);
    h1cosmetic(h1_Wjets,        "Wjets",    kBlack, 1, kGray+1,     var);
    h1cosmetic(h1_Ztt,          "Zjets",    kBlack, 1, kGreen+2,    var);

}

//
// main  
//
void getPostUncert(bool do_qqWW, TString histofile, int TOTALTOY, char* TESTNAME){ 

    gStyle->SetOptFit(1011); 

    TFile *histo_rootfile; 
    if(do_qqWW) histo_rootfile = new TFile(Form("histoqqww_%s.root",TESTNAME));
    else histo_rootfile = new TFile(Form("histobkg_%s.root",TESTNAME)); 
    // 
    // final histgrams for uncertainty
    // 
    float mllbins_CR1[10]   =   {12,30,45,60,75,100,125,150,175,200};
    float mllbins_CR2[7]    =   {60,75,100,125,150,175,200};
    float mTbins_CR1[9]     =   {120,140,160,180,200,220,240,260,280};
    float mTbins_CR2[7]     =   {60,70,80,90,100,110,120};
    TH1F* h1_mll_CR1_uncert = new TH1F(Form("h1_%s_mll_CR1_uncert",do_qqWW?"qqWW":"bkg"), 
                                       Form("h1_%s_mll_CR1_uncert",do_qqWW?"qqWW":"bkg"), 9, mllbins_CR1);
    TH1F* h1_mll_CR2_uncert = new TH1F(Form("h1_%s_mll_CR2_uncert",do_qqWW?"qqWW":"bkg"), 
                                       Form("h1_%s_mll_CR2_uncert",do_qqWW?"qqWW":"bkg"), 6, mllbins_CR2);
    TH1F* h1_mT_CR1_uncert  = new TH1F(Form("h1_%s_mT_CR1_uncert",do_qqWW?"qqWW":"bkg"), 
                                       Form("h1_%s_mT_CR1_uncert",do_qqWW?"qqWW":"bkg"), 8, mTbins_CR1);
    TH1F* h1_mT_CR2_uncert  = new TH1F(Form("h1_%s_mT_CR2_uncert",do_qqWW?"qqWW":"bkg"), 
                                       Form("h1_%s_mT_CR2_uncert",do_qqWW?"qqWW":"bkg"), 6, mTbins_CR2);
    
    std::vector<float> vec_mll_CR1_uncert_bin[9];
    std::vector<float> vec_mll_CR2_uncert_bin[6];
    std::vector<float> vec_mT_CR1_uncert_bin[8];
    std::vector<float> vec_mT_CR2_uncert_bin[6];
    // loop over toys and do projection 
    for(int itoy=0; itoy<TOTALTOY; itoy++) {
        TH2F* h2_tmp_drawCR1;
        if(do_qqWW) h2_tmp_drawCR1  = (TH2F*)histo_rootfile->Get(Form("h2post_qqWW_CR1_%i",itoy));  
        else h2_tmp_drawCR1         = (TH2F*)histo_rootfile->Get(Form("h2post_bkg_%i",itoy));  
        TH1F* h1_tmp_mll_drawCR1    = doProjection(h2_tmp_drawCR1,         "mll",   "CR1");   
        TH1F* h1_tmp_mT_drawCR1     = doProjection(h2_tmp_drawCR1,         "mT",    "CR1");  
        TH2F* h2_tmp_drawCR2;
        if(do_qqWW) h2_tmp_drawCR2  = (TH2F*)histo_rootfile->Get(Form("h2post_qqWW_CR2_%i",itoy)); 
        else h2_tmp_drawCR2         = (TH2F*)histo_rootfile->Get(Form("h2post_bkg_%i",itoy)); 
        TH1F* h1_tmp_mll_drawCR2    = doProjection(h2_tmp_drawCR2,         "mll",   "CR2");   
        TH1F* h1_tmp_mT_drawCR2     = doProjection(h2_tmp_drawCR2,         "mT",    "CR2");  
        for(int i=0; i<9; i++) vec_mll_CR1_uncert_bin[i].push_back(h1_tmp_mll_drawCR1->GetBinContent(i+1));
        for(int i=0; i<6; i++) vec_mll_CR2_uncert_bin[i].push_back(h1_tmp_mll_drawCR2->GetBinContent(i+1));
        for(int i=0; i<8; i++) vec_mT_CR1_uncert_bin[i].push_back(h1_tmp_mT_drawCR1->GetBinContent(i+1));
        for(int i=0; i<6; i++) vec_mT_CR2_uncert_bin[i].push_back(h1_tmp_mT_drawCR2->GetBinContent(i+1));
        
    }        


    float mean_mll_CR1[9],  uncert_mll_CR1[9];    
    float mean_mll_CR2[6],  uncert_mll_CR2[6];    
    float mean_mT_CR1[8],   uncert_mT_CR1[8];    
    float mean_mT_CR2[6],   uncert_mT_CR2[6];     

    cout << "CR1 mll" << endl;
    for(int i=0; i<9; i++) { 
        mean_mll_CR1[i] = TMath::Mean(TOTALTOY, &vec_mll_CR1_uncert_bin[i][0]);
        uncert_mll_CR1[i] = TMath::RMS(TOTALTOY, &vec_mll_CR1_uncert_bin[i][0]);
        if(DEBUG) cout << i << " ::: " << mean_mll_CR1[i] << " +/- " << uncert_mll_CR1[i] << endl;  
        h1_mll_CR1_uncert->SetBinContent(i+1, mean_mll_CR1[i]);
        h1_mll_CR1_uncert->SetBinError(i+1, uncert_mll_CR1[i]);
    }
    cout << "CR2 mll" << endl;
    for(int i=0; i<6; i++) {
        mean_mll_CR2[i] = TMath::Mean(TOTALTOY, &vec_mll_CR2_uncert_bin[i][0]);
        uncert_mll_CR2[i] = TMath::RMS(TOTALTOY, &vec_mll_CR2_uncert_bin[i][0]);
        if(DEBUG) cout << i << " ::: " << mean_mll_CR2[i] << " +/- " << uncert_mll_CR2[i] << endl;  
        h1_mll_CR2_uncert->SetBinContent(i+1, mean_mll_CR2[i]);
        h1_mll_CR2_uncert->SetBinError(i+1, uncert_mll_CR2[i]);
    }
    cout << "CR1 mT" << endl;
    for(int i=0; i<8; i++) {
        mean_mT_CR1[i] = TMath::Mean(TOTALTOY, &vec_mT_CR1_uncert_bin[i][0]);
        uncert_mT_CR1[i] = TMath::RMS(TOTALTOY, &vec_mT_CR1_uncert_bin[i][0]);
        if(DEBUG) cout << i << " ::: " << mean_mT_CR1[i] << " +/- " << uncert_mT_CR1[i] << endl;  
        h1_mT_CR1_uncert->SetBinContent(i+1, mean_mT_CR1[i]);
        h1_mT_CR1_uncert->SetBinError(i+1, uncert_mT_CR1[i]);
    }
    cout << "CR2 mT" << endl;
    for(int i=0; i<6; i++) {
        mean_mT_CR2[i] = TMath::Mean(TOTALTOY, &vec_mT_CR2_uncert_bin[i][0]);
        uncert_mT_CR2[i] = TMath::RMS(TOTALTOY, &vec_mT_CR2_uncert_bin[i][0]);
        if(DEBUG) cout << i << " ::: " << mean_mT_CR2[i] << " +/- " << uncert_mT_CR2[i] << endl;  
        h1_mT_CR2_uncert->SetBinContent(i+1, mean_mT_CR2[i]);
        h1_mT_CR2_uncert->SetBinError(i+1, uncert_mT_CR2[i]);
    }

    // 
    // histgrams for each bin
    // 
    TH1F* h1_mll_CR1_uncert_bin[9];
    TH1F* h1_mll_CR2_uncert_bin[6];
    TH1F* h1_mT_CR1_uncert_bin[8];
    TH1F* h1_mT_CR2_uncert_bin[6];

    for(int i=0; i<9; i++) h1_mll_CR1_uncert_bin[i] = new TH1F(Form("h1_mll_CR1_uncert_bin%i",i+1), 
                                                               Form("h1_mll_CR1_uncert_bin%i",i+1), 
                                                               100, mean_mll_CR1[i]-5*uncert_mll_CR1[i], mean_mll_CR1[i]+5*uncert_mll_CR1[i]);
    for(int i=0; i<6; i++) h1_mll_CR2_uncert_bin[i] = new TH1F(Form("h1_mll_CR2_uncert_bin%i",i+1), 
                                                               Form("h1_mll_CR2_uncert_bin%i",i+1), 
                                                               100, mean_mll_CR2[i]-5*uncert_mll_CR2[i], mean_mll_CR2[i]+5*uncert_mll_CR2[i]);
    for(int i=0; i<8; i++) h1_mT_CR1_uncert_bin[i]  = new TH1F(Form("h1_mT_CR1_uncert_bin%i",i+1), 
                                                               Form("h1_mT_CR1_uncert_bin%i",i+1), 
                                                               100, mean_mT_CR1[i]-5*uncert_mT_CR1[i], mean_mT_CR1[i]+5*uncert_mT_CR1[i]);
    for(int i=0; i<6; i++) h1_mT_CR2_uncert_bin[i]  = new TH1F(Form("h1_mT_CR2_uncert_bin%i",i+1), 
                                                               Form("h1_mT_CR2_uncert_bin%i",i+1), 
                                                               100, mean_mT_CR2[i]-5*uncert_mT_CR2[i], mean_mT_CR2[i]+5*uncert_mT_CR2[i]);

    // Fill histograms 
    for(int itoy=0; itoy<TOTALTOY; itoy++) {
        for(int i=0; i<9; i++) h1_mll_CR1_uncert_bin[i]->Fill((Double_t)vec_mll_CR1_uncert_bin[i][itoy]);
        for(int i=0; i<6; i++) h1_mll_CR2_uncert_bin[i]->Fill(vec_mll_CR2_uncert_bin[i][itoy]);
        for(int i=0; i<8; i++) h1_mT_CR1_uncert_bin[i]->Fill(vec_mT_CR1_uncert_bin[i][itoy]);
        for(int i=0; i<6; i++) h1_mT_CR2_uncert_bin[i]->Fill(vec_mT_CR2_uncert_bin[i][itoy]);
    } 
    
    if(DEBUG) { // debug 
        // Fit histograms 
        TCanvas *c_fit = new TCanvas("c_fit", "c_fit");
        c_fit->cd(1);
        for(int i=0; i<9; i++) h1_mll_CR1_uncert_bin[i]->Fit("gaus");
        for(int i=0; i<6; i++) h1_mll_CR2_uncert_bin[i]->Fit("gaus");
        for(int i=0; i<8; i++) h1_mT_CR1_uncert_bin[i]->Fit("gaus");
        for(int i=0; i<6; i++) h1_mT_CR2_uncert_bin[i]->Fit("gaus");
        delete c_fit; 

        TCanvas *c_mll_CR2 = new TCanvas("c_mll_CR2", "c_mll_CR2", 900, 900); c_mll_CR2->Divide(3,3); 
        for(int i=0; i<9; i++) { c_mll_CR2->cd(i+1); h1_mll_CR1_uncert_bin[i]->Draw();}
        c_mll_CR2->Print(Form("/home/users/jaehyeok/public_html/misc/c_mll_CR2_%s.pdf",do_qqWW?"qqWW":"bkg"));
    
        TCanvas *c_mll_CR2 = new TCanvas("c_mll_CR2", "c_mll_CR2", 900, 900); c_mll_CR2->Divide(3,3); 
        for(int i=0; i<6; i++) { c_mll_CR2->cd(i+1); h1_mll_CR2_uncert_bin[i]->Draw(); }
        c_mll_CR2->Print(Form("/home/users/jaehyeok/public_html/misc/c_mll_CR2_%s.pdf",do_qqWW?"qqWW":"bkg"));

        TCanvas *c_mT_CR1 = new TCanvas("c_mT_CR1", "c_mT_CR1", 900, 900); c_mT_CR1->Divide(3,3); 
        for(int i=0; i<8; i++) { c_mT_CR1->cd(i+1); h1_mT_CR1_uncert_bin[i]->Draw(); }
        c_mT_CR1->Print(Form("/home/users/jaehyeok/public_html/misc/c_mT_CR1_%s.pdf",do_qqWW?"qqWW":"bkg"));

        TCanvas *c_mT_CR2 = new TCanvas("c_mT_CR2", "c_mT_CR2", 900, 900); c_mT_CR2->Divide(3,3); 
        for(int i=0; i<6; i++) { c_mT_CR2->cd(i+1); h1_mT_CR2_uncert_bin[i]->Draw(); }
        c_mT_CR2->Print(Form("/home/users/jaehyeok/public_html/misc/c_mT_CR2_%s.pdf",do_qqWW?"qqWW":"bkg"));
    } 

    // save uncert histograms in the rootfile
    TFile *outputrootfile = new TFile(histofile, "UPDATE");
    gROOT->cd();
    outputrootfile->cd();
    h1_mll_CR1_uncert->SetDirectory(0);  h1_mll_CR1_uncert->Write();
    h1_mll_CR2_uncert->SetDirectory(0);  h1_mll_CR2_uncert->Write();
    h1_mT_CR1_uncert->SetDirectory(0);  h1_mT_CR1_uncert->Write();
    h1_mT_CR2_uncert->SetDirectory(0);  h1_mT_CR2_uncert->Write(); 
    outputrootfile->Close();
    
    histo_rootfile->Close();
}

void getPostUncert_runAll(int TOTALTOY, char* TESTNAME){
       
    TString uncert_rootfilename = Form("uncert_%s.root",TESTNAME);

    gSystem->Exec(Form("rm %s", uncert_rootfilename.Data()));
    TFile *uncert_rootfile = new TFile(uncert_rootfilename, "CREATE");
    uncert_rootfile->Close();

    getPostUncert(true,     uncert_rootfilename, TOTALTOY, TESTNAME);
    getPostUncert(false,    uncert_rootfilename, TOTALTOY, TESTNAME);

}
