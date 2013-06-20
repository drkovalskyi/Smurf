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
  
    TFile* postFile_CR1 = TFile::Open(Form("fitoutput_%s/fit_newdefault_usingCR1_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    TFile* postFile_CR2 = TFile::Open(Form("fitoutput_%s/fit_newdefault_usingCR2_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    
    TH1F *h1post_qqWW_CR1, *h1post_qqWW_CR2;
    TH2F *h2post_qqWW_CR1, *h2post_qqWW_CR2;

    // get qqWW for nominal and control fit 
    h1post_qqWW_CR1         = makeHistogram( (TGraphAsymmErrors*) postFile_CR1->Get("j0of_qqWW"));   
    h1post_qqWW_CR2         = makeHistogram( (TGraphAsymmErrors*) postFile_CR2->Get("j0of_qqWW"));   
    
    // make 2d  
    h2post_qqWW_CR1         = Roll1DTo2D(h1post_qqWW_CR1); 
    h2post_qqWW_CR1->SetTitle(Form("h2post_qqWW_CR1_%i",toynumber));   h2post_qqWW_CR1->SetName(Form("h2post_qqWW_CR1_%i",toynumber));
    h2post_qqWW_CR2         = Roll1DTo2D(h1post_qqWW_CR2); 
    h2post_qqWW_CR2->SetTitle(Form("h2post_qqWW_CR2_%i",toynumber));   h2post_qqWW_CR2->SetName(Form("h2post_qqWW_CR2_%i",toynumber));

    if(0) { // debug 
        cout << h2post_qqWW_CR1->Integral(1,6,1,3) << endl;
        cout << h2post_qqWW_CR2->Integral(1,6,1,3) << endl;

        TCanvas *c_debug = new TCanvas("c_debug", "c_debug", 800, 400); 
        c_debug->Divide(2,1); 
        c_debug->cd(1);
        h2post_qqWW_CR1->Draw("colz");
        c_debug->cd(2);
        h2post_qqWW_CR2->Draw("colz");

        c_debug->Print(Form("/home/users/jaehyeok/public_html/misc/2d_debug_%i.pdf", toynumber));

        delete c_debug; 
    }

    
    TFile *outputrootfile = new TFile(histofile, "UPDATE");
    gROOT->cd();
    outputrootfile->cd();
    h2post_qqWW_CR1->SetDirectory(0); h2post_qqWW_CR1->Write();
    h2post_qqWW_CR2->SetDirectory(0); h2post_qqWW_CR2->Write();
    outputrootfile->Close();

    postFile_CR1->Close(); 
    postFile_CR2->Close();  
    
}

void getpostbkg(char* histofile, int toynumber=0, char* TESTNAME) {

    TFile* postFile_nominal = TFile::Open(Form("fitoutput_%s/fit_%i_fittedShape_floatMu.root", TESTNAME, toynumber)); 
    
    TH1F *h1post_ZH, *h1post_WH, *h1post_qqH, *h1post_ggH, *h1post_ggWW, *h1post_VV, *h1post_Top,
         *h1post_Zjets, *h1post_WjetsE, *h1post_WjetsM, *h1post_Wgamma, *h1post_Wg3l, *h1post_Ztt;
    TH2F *h2post_ZH, *h2post_WH, *h2post_qqH, *h2post_ggH, *h2post_ggWW, *h2post_VV, *h2post_Top,
         *h2post_Zjets, *h2post_WjetsE, *h2post_WjetsM, *h2post_Wgamma, *h2post_Wg3l, *h2post_Ztt;
    TH2F *h2post_bkg;

    h1post_ZH     = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_ZH"));     
    h1post_WH     = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_WH"));     
    h1post_qqH    = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_qqH"));    
    h1post_ggH    = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_ggH"));    
    h1post_ggWW   = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_ggWW"));   
    h1post_VV     = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_VV"));     
    h1post_Top    = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_Top"));    
    h1post_Zjets  = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_Zjets"));    
    h1post_WjetsE = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_WjetsE")); 
    h1post_WjetsM = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_WjetsM")); 
    h1post_Wgamma = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_Wgamma")); 
    h1post_Wg3l   = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_Wg3l"));   
    h1post_Ztt    = makeHistogram( (TGraphAsymmErrors*) postFile_nominal->Get("j0of_Ztt"));    

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

void getposthisto_runAll(int TOTALTOY, char* TESTNAME){

    bool do_qqWW    = 1;
    bool do_bkg     = 1;

    // 
    // qqWW 
    // 
    if(do_qqWW) {
        // remove previously existing output file 
        gSystem->Exec(Form("rm histoqqww_%s.root",TESTNAME));

        // create a root file 
        TFile *qqwwhisto_rootfile = new TFile(Form("histoqqww_%s.root",TESTNAME), "CREATE");
        qqwwhisto_rootfile->Close();

        // loop over toys 
        for(int itoy=0; itoy<TOTALTOY; itoy++)  {
            cout << "... doing qqWW toy " << itoy << endl;
            getpostqqWW(Form("histoqqww_%s.root",TESTNAME),itoy,TESTNAME); 
        }
    }

    // 
    // Other backgrounds 
    // 
    if(do_bkg) {
        // remove previously existing output file 
        gSystem->Exec(Form("rm histobkg_%s.root",TESTNAME));

        // create a root file 
        TFile *bkghisto_rootfile = new TFile(Form("histobkg_%s.root",TESTNAME), "CREATE");
        bkghisto_rootfile->Close();

        // loop over toys 
        for(int itoy=0; itoy<TOTALTOY; itoy++) {  
            cout << "... doing bkg toy " << itoy << endl;
            getpostbkg(Form("histobkg_%s.root",TESTNAME),itoy,TESTNAME); 
        }
    }

}
