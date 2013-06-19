#include <vector>
#include "TMath.h"

const bool DEBUG = true;

//
// Zoom
//
TH2F* Zoom(TH2F* h2) {

    float mtbins[7]    =   {60,70,80,90,100,110,120};
    float mllbins[4]   =   {12,30,45,60};

    TH2F* h2zoom = new TH2F(Form("%s_zoom", h2->GetName()), Form("%s_zoom", h2->GetTitle()), 6, mtbins, 3, mllbins);

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

    float mtbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
    float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};
    TH2F* h2 = new TH2F(Form("%s_2D", h1->GetName()), Form("%s_2D",h1->GetTitle()), 14, mtbins, 9, mllbins);
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
    //TH2F *h2zoom = (TH2F*)Zoom(h2);

    //if(doZoom) return h2zoom;
    //else return h2;
    return h2;
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

    //float mtbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
    //float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};
  

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
            float mtbins_CR1[9]    =   {120,140,160,180,200,220,240,260,280};
            TH1F* h1_mT_CR1  = new TH1F(Form("h1_mT_CR1_%s", histoname.Data()), Form("h1_mT_CR1_%s", histoname.Data()), 8, mtbins_CR1);
            for(int imt=7; imt<h2->GetXaxis()->GetNbins()+1; imt++){ // 7 - 14
                for(int imll=1; imll<h2->GetYaxis()->GetNbins()+1; imll++){ // 1 - 9
                    h1_mT_CR1->SetBinContent(imt-6, h2->GetBinContent(imt,imll)+h1_mT_CR1->GetBinContent(imt-6));
                }
            }
            return h1_mT_CR1;
        } else if (plotregion=="CR2") {  // - CR2 --------------------------------
            float mtbins_CR2[7]    =   {60,70,80,90,100,110,120};
            TH1F* h1_mT_CR2  = new TH1F(Form("h1_mT_CR2_%s", histoname.Data()), Form("h1_mT_CR2_%s", histoname.Data()), 6, mtbins_CR2);
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

TH2F* CloneTH2F(TH2F* h2, char* title="", int toynumber){
    TH2F* h2_tmp = (TH2F*)h2->Clone(Form("%s_%i", title, toynumber)); 
    h2_tmp->SetTitle(Form("%s_%i", title, toynumber));
    h2_tmp->SetName(Form("%s_%i", title, toynumber));
    return h2_tmp;
}

//
// main function
//
void getpostqqWW(char* histofile, int toynumber=0, char* TESTNAME) {
  
    TFile* postFile_newdefault = TFile::Open(Form("fitoutput_%s/fit_newdefault_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    TFile* postFile_CR1 = TFile::Open(Form("fitoutput_%s/fit_CR1_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    TFile* postFile_CR2 = TFile::Open(Form("fitoutput_%s/fit_CR2_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    
    TH1F *h1post_qqWW_newdefault, *h1post_qqWW_CR1, *h1post_qqWW_CR2;
    TH2F *h2post_qqWW_newdefault, *h2post_qqWW_CR1, *h2post_qqWW_CR2;
    TH2F *h2post_qqWW_newdefault_renormbyCR1, *h2post_qqWW_newdefault_renormbyCR2;

    // get qqWW for nominal and control fit 
    h1post_qqWW_newdefault  = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_qqWW"));   
    h1post_qqWW_CR1         = makeHistogram( (TGraphAsymmErrors*) postFile_CR1->Get("j0of_qqWW"));   
    h1post_qqWW_CR2         = makeHistogram( (TGraphAsymmErrors*) postFile_CR2->Get("j0of_qqWW"));   
    
    // make 2d  
    h2post_qqWW_newdefault  = Roll1DTo2D(h1post_qqWW_newdefault);
    h2post_qqWW_newdefault->SetTitle("h2post_qqWW_newdefault");
    h2post_qqWW_newdefault->SetName("h2post_qqWW_newdefault");
    h2post_qqWW_CR1         = Roll1DTo2D(h1post_qqWW_CR1); 
    h2post_qqWW_CR2         = Roll1DTo2D(h1post_qqWW_CR2); 
    h2post_qqWW_newdefault_renormbyCR1 = CloneTH2F(h2post_qqWW_newdefault,"h2post_qqWW_newdefault_renormbyCR1", toynumber); 
    h2post_qqWW_newdefault_renormbyCR2 = CloneTH2F(h2post_qqWW_newdefault,"h2post_qqWW_newdefault_renormbyCR2", toynumber);

    if(0) { // debug 
        cout << h2post_qqWW_newdefault->Integral(1,6,1,3) << endl;
        cout << h2post_qqWW_CR1->Integral(1,6,1,3) << endl;
        cout << h2post_qqWW_CR2->Integral(1,6,1,3) << endl;

        TCanvas *c_debug = new TCanvas("c_debug", "c_debug", 1200, 400); 
        c_debug->Divide(3,1); 
        c_debug->cd(1);
        h2post_qqWW_newdefault->Draw("colz");
        c_debug->cd(2);
        h2post_qqWW_CR1->Draw("colz");
        c_debug->cd(3);
        h2post_qqWW_CR2->Draw("colz");

        c_debug->Print(Form("/home/users/jaehyeok/public_html/misc/2d_debug_%i.pdf", toynumber));

        delete c_debug; 
    }

    
    //
    // renormalize qqWW using the other control regions  
    //
    float qqww_CR1_normCR1  = h2post_qqWW_CR1->Integral(7,14,1,9); 
    float qqww_newdefault_normCR1  = h2post_qqWW_newdefault_renormbyCR1->Integral(7,14,1,9); 
    h2post_qqWW_newdefault_renormbyCR1->Scale(qqww_CR1_normCR1/qqww_newdefault_normCR1);    

    float qqww_CR2_normCR2  = h2post_qqWW_CR2->Integral(1,6,4,9); 
    float qqww_newdefault_normCR2  = h2post_qqWW_newdefault_renormbyCR2->Integral(1,6,4,9); 
    h2post_qqWW_newdefault_renormbyCR2->Scale(qqww_CR2_normCR2/qqww_newdefault_normCR2);    

    if(0) { // debug 
        cout << "newdefault       : " << h2post_qqWW_newdefault->Integral(1,6,1,3) << endl;
        cout << "CR1              : " << h2post_qqWW_CR1->Integral(1,6,1,3) << endl;
        cout << "CR2              : " << h2post_qqWW_CR2->Integral(1,6,1,3) << endl;
        cout << "renorm using CR1 : " << h2post_qqWW_newdefault_renormbyCR1->Integral(1,6,1,3) << endl;
        cout << "renorm using CR2 : " << h2post_qqWW_newdefault_renormbyCR2->Integral(1,6,1,3) << endl;
    } 
    
    TFile *outputrootfile = new TFile(histofile, "UPDATE");
    gROOT->cd();
    outputrootfile->cd();
    h2post_qqWW_newdefault->SetDirectory(0); h2post_qqWW_newdefault->Write();
    h2post_qqWW_newdefault_renormbyCR1->SetDirectory(0); h2post_qqWW_newdefault_renormbyCR1->Write();
    h2post_qqWW_newdefault_renormbyCR2->SetDirectory(0); h2post_qqWW_newdefault_renormbyCR2->Write();
    outputrootfile->Close();

    postFile_newdefault->Close();
    postFile_CR1->Close(); 
    postFile_CR2->Close();  
    
}

void getpostbkg(char* histofile, int toynumber=0, char* TESTNAME) {

    TFile* postFile_newdefault = TFile::Open(Form("fitoutput_%s/fit_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    
    TH1F *h1post_ZH, *h1post_WH, *h1post_qqH, *h1post_ggH, *h1post_ggWW, *h1post_VV, *h1post_Top,
         *h1post_Zjets, *h1post_WjetsE, *h1post_WjetsM, *h1post_Wgamma, *h1post_Wg3l, *h1post_Ztt;
    TH2F *h2post_ZH, *h2post_WH, *h2post_qqH, *h2post_ggH, *h2post_ggWW, *h2post_VV, *h2post_Top,
         *h2post_Zjets, *h2post_WjetsE, *h2post_WjetsM, *h2post_Wgamma, *h2post_Wg3l, *h2post_Ztt;
    TH2F *h2post_bkg;

    h1post_ZH     = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_ZH"));     
    h1post_WH     = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_WH"));     
    h1post_qqH    = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_qqH"));    
    h1post_ggH    = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_ggH"));    
    h1post_ggWW   = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_ggWW"));   
    h1post_VV     = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_VV"));     
    h1post_Top    = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_Top"));    
    h1post_Zjets  = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_Zjets"));    
    h1post_WjetsE = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_WjetsE")); 
    h1post_WjetsM = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_WjetsM")); 
    h1post_Wgamma = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_Wgamma")); 
    h1post_Wg3l   = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_Wg3l"));   
    h1post_Ztt    = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get("j0of_Ztt"));    

    h2post_ZH       = Roll1DTo2D(h1post_ZH);
    h2post_WH       = Roll1DTo2D(h1post_WH);
    h2post_qqH      = Roll1DTo2D(h1post_qqH);
    h2post_ggH      = Roll1DTo2D(h1post_ggH);
    h2post_ggWW     = Roll1DTo2D(h1post_ggWW);
    h2post_VV       = Roll1DTo2D(h1post_VV);
    h2post_Top      = Roll1DTo2D(h1post_Top);
    h2post_Zjets    = Roll1DTo2D(h1post_Zjets);
    h2post_WjetsE   = Roll1DTo2D(h1post_WjetsE);
    h2post_WjetsM   = Roll1DTo2D(h1post_WjetsM);
    h2post_Wgamma   = Roll1DTo2D(h1post_Wgamma);
    h2post_Wg3l     = Roll1DTo2D(h1post_Wg3l);
    h2post_Ztt      = Roll1DTo2D(h1post_Ztt); 

    h2post_bkg = (TH2F*)h2post_ggWW->Clone(Form("h2post_bkg_%i", toynumber)); 
    h2post_bkg->SetName(Form("h2post_bkg_%i", toynumber)); 
    h2post_bkg->SetTitle("h2post_bkg");
    h2post_bkg->Add(h2post_VV);
    h2post_bkg->Add(h2post_Top);
    h2post_bkg->Add(h2post_Zjets);
    h2post_bkg->Add(h2post_WjetsE);
    h2post_bkg->Add(h2post_WjetsM);
    h2post_bkg->Add(h2post_Wgamma);
    h2post_bkg->Add(h2post_Wg3l);
    h2post_bkg->Add(h2post_Ztt);

    TFile *outputrootfile = new TFile(histofile, "UPDATE");
    gROOT->cd();
    outputrootfile->cd();
    h2post_bkg->SetDirectory(0); h2post_bkg->Write();
    outputrootfile->Close();
}

void getpostone(char* histofile, int toynumber=0, char* TESTNAME, char* process) {

    TFile* postFile_newdefault = TFile::Open(Form("fitoutput_%s/fit_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    
    TH1F *h1post;
    TH2F *h2post; 

    h1post     = makeHistogram( (TGraphAsymmErrors*) postFile_newdefault->Get(Form("j0of_%s",process)));     
    h2post     = Roll1DTo2D(h1post);
    h2post->SetName(Form("h2post_%s_%i",    process, toynumber)); 
    h2post->SetTitle(Form("h2post_%s_%i",   process, toynumber));


    TFile *outputrootfile = new TFile(histofile, "UPDATE");
    gROOT->cd();
    outputrootfile->cd();
    h2post->SetDirectory(0); h2post->Write();
    outputrootfile->Close();

    postFile_newdefault->Close();
}

void getPostUncert(bool do_qqWW, int TOTALTOY, char* TESTNAME, char* process){ 

    gStyle->SetOptFit(1011); 

    TFile *histo_rootfile; 
    histo_rootfile = new TFile(Form("histo_%s.root",process)); 
    // 
    // final histgrams for uncertainty
    // 
    float mllbins_CR1[10]   =   {12,30,45,60,75,100,125,150,175,200};
    float mllbins_CR2[7]    =   {60,75,100,125,150,175,200};
    float mTbins_CR1[9]     =   {120,140,160,180,200,220,240,260,280};
    float mTbins_CR2[7]     =   {60,70,80,90,100,110,120};
    TH1F* h1_mll_CR1_uncert = new TH1F(Form("h1_%s_mll_CR1_uncert",do_qqWW?"qqWW":process), 
                                       Form("h1_%s_mll_CR1_uncert",do_qqWW?"qqWW":process), 9, mllbins_CR1);
    TH1F* h1_mll_CR2_uncert = new TH1F(Form("h1_%s_mll_CR2_uncert",do_qqWW?"qqWW":process), 
                                       Form("h1_%s_mll_CR2_uncert",do_qqWW?"qqWW":process), 6, mllbins_CR2);
    TH1F* h1_mT_CR1_uncert  = new TH1F(Form("h1_%s_mT_CR1_uncert",do_qqWW?"qqWW":process), 
                                       Form("h1_%s_mT_CR1_uncert",do_qqWW?"qqWW":process), 8, mTbins_CR1);
    TH1F* h1_mT_CR2_uncert  = new TH1F(Form("h1_%s_mT_CR2_uncert",do_qqWW?"qqWW":process), 
                                       Form("h1_%s_mT_CR2_uncert",do_qqWW?"qqWW":process), 6, mTbins_CR2);
    
    std::vector<float> vec_mll_CR1_uncert_bin[9];
    std::vector<float> vec_mll_CR2_uncert_bin[6];
    std::vector<float> vec_mT_CR1_uncert_bin[8];
    std::vector<float> vec_mT_CR2_uncert_bin[6];
    // loop over toys and do projection 
    for(int itoy=0; itoy<TOTALTOY; itoy++) {
        TH2F* h2_tmp_drawCR1;
        if(do_qqWW) h2_tmp_drawCR1  = (TH2F*)histo_rootfile->Get(Form("h2post_qqWW_newdefault_renormbyCR2_%i",itoy));  
        else h2_tmp_drawCR1         = (TH2F*)histo_rootfile->Get(Form("h2post_%s_%i", process, itoy));  
        TH1F* h1_tmp_mll_drawCR1    = doProjection(h2_tmp_drawCR1,         "mll",   "CR1");   
        TH1F* h1_tmp_mT_drawCR1     = doProjection(h2_tmp_drawCR1,         "mT",    "CR1");  
        TH2F* h2_tmp_drawCR2;
        if(do_qqWW) h2_tmp_drawCR2  = (TH2F*)histo_rootfile->Get(Form("h2post_qqWW_newdefault_renormbyCR1_%i",itoy)); 
        else h2_tmp_drawCR2         = (TH2F*)histo_rootfile->Get(Form("h2post_%s_%i",process, itoy)); 
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
    
    // Fit histograms 
    TCanvas *c_fit = new TCanvas("c_fit", "c_fit");
    c_fit->cd(1);
    for(int i=0; i<9; i++) h1_mll_CR1_uncert_bin[i]->Fit("gaus");
    for(int i=0; i<6; i++) h1_mll_CR2_uncert_bin[i]->Fit("gaus");
    for(int i=0; i<8; i++) h1_mT_CR1_uncert_bin[i]->Fit("gaus");
    for(int i=0; i<6; i++) h1_mT_CR2_uncert_bin[i]->Fit("gaus");
    delete c_fit; 

    if(DEBUG) { // debug 
        TCanvas *c_debug1 = new TCanvas("c_debug1", "c_debug1", 900, 900); c_debug1->Divide(3,3); 
        for(int i=0; i<9; i++) { c_debug1->cd(i+1); h1_mll_CR1_uncert_bin[i]->Draw();}
        c_debug1->Print(Form("/home/users/jaehyeok/public_html/misc/c_debug1_%s.pdf",do_qqWW?"qqWW":process));
    
        TCanvas *c_debug2 = new TCanvas("c_debug2", "c_debug2", 900, 900); c_debug2->Divide(3,3); 
        for(int i=0; i<6; i++) { c_debug2->cd(i+1); h1_mll_CR2_uncert_bin[i]->Draw(); }
        c_debug2->Print(Form("/home/users/jaehyeok/public_html/misc/c_debug2_%s.pdf",do_qqWW?"qqWW":process));

        TCanvas *c_debug3 = new TCanvas("c_debug3", "c_debug3", 900, 900); c_debug3->Divide(3,3); 
        for(int i=0; i<8; i++) { c_debug3->cd(i+1); h1_mT_CR1_uncert_bin[i]->Draw(); }
        c_debug3->Print(Form("/home/users/jaehyeok/public_html/misc/c_debug3_%s.pdf",do_qqWW?"qqWW":process));

        TCanvas *c_debug4 = new TCanvas("c_debug4", "c_debug4", 900, 900); c_debug4->Divide(3,3); 
        for(int i=0; i<6; i++) { c_debug4->cd(i+1); h1_mT_CR2_uncert_bin[i]->Draw(); }
        c_debug4->Print(Form("/home/users/jaehyeok/public_html/misc/c_debug4_%s.pdf",do_qqWW?"qqWW":process));
    } 

/*
    // save uncert histograms in the rootfile
    TFile *outputrootfile = new TFile(histofile, "UPDATE");
    gROOT->cd();
    outputrootfile->cd();
    h1_mll_CR1_uncert->SetDirectory(0);  h1_mll_CR1_uncert->Write();
    h1_mll_CR2_uncert->SetDirectory(0);  h1_mll_CR2_uncert->Write();
    h1_mT_CR1_uncert->SetDirectory(0);  h1_mT_CR1_uncert->Write();
    h1_mT_CR2_uncert->SetDirectory(0);  h1_mT_CR2_uncert->Write(); 
    outputrootfile->Close();
 */   

//    histo_rootfile->Close();
}

void getPostUncert_One(int TOTALTOY=1000, char* TESTNAME="NOSYST", char* process="WjetsM"){

        gSystem->Exec(Form("rm histo_%s.root",process));

        // create a root file 
        TFile *histo_rootfile = new TFile(Form("histo_%s.root",process), "CREATE");
        histo_rootfile->Close();

        // loop over toys 
        for(int itoy=0; itoy<TOTALTOY; itoy++)  {
            cout << "... reading " << process << " from toy number " << itoy << endl;
            getpostone(Form("histo_%s.root",process),itoy,TESTNAME,process); 
        }
   
        getPostUncert(false, TOTALTOY, TESTNAME, process);
}
