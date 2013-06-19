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

void getUncerBand(TH1F* &h1_uncert, TH1F* &h1_WW, TH1F* &h1_bkg, char* var="mll", char* plotregion="CR1"){

    
    TFile *uncert_rootfile = new TFile("uncert_FullqqWWSyst.root", "READ");

    cout << "... Making Uncertainty band for " << var << " in " << plotregion << endl;
    
    if( var!="mll" && var!="mT" ) { 
        cout << "ERROR : wrong variable name " << endl;
        return; 
    }
    
    if(plotregion!="CR1" && plotregion!="CR2") { 
        cout << "ERROR : wrong control region name " << endl;
        return; 
    }

    for(int i=1; i<h1_WW->GetXaxis()->GetNbins()+1; i++){ 
        
        TH1F *h1_qqWW_tmp   = (TH1F*)uncert_rootfile->Get(Form("h1_qqWW_%s_%s_uncert", var, plotregion));
        TH1F *h1_bkg_tmp    = (TH1F*)uncert_rootfile->Get(Form("h1_bkg_%s_%s_uncert", var, plotregion)); 
        
        if(h1_WW->GetXaxis()->GetNbins() != h1_qqWW_tmp->GetXaxis()->GetNbins()) {
            cout << " ERROR :: something is wrong " << endl;
        }
        // fractional error
        float frac_qqWW_error   = h1_qqWW_tmp->GetBinError(i) / h1_qqWW_tmp->GetBinContent(i); 
        float frac_bkg_error    = h1_bkg_tmp->GetBinError(i) / h1_bkg_tmp->GetBinContent(i);  
        float qqWW_error        = frac_qqWW_error * h1_WW->GetBinContent(i);
        float bkg_error         = frac_bkg_error * h1_bkg->GetBinContent(i);

        h1_uncert->SetBinError(i, TMath::Sqrt(qqWW_error*qqWW_error+bkg_error*bkg_error)); 

        if(1) { // debug
            cout << i << " :: " 
                 << frac_qqWW_error << " " 
                 << frac_bkg_error << " "  
                 << h1_WW->GetBinContent(i) << " " 
                 << h1_bkg->GetBinContent(i) << " " 
                 << qqWW_error << " " 
                 << bkg_error << " " 
                 << endl;
        }
    }

    h1_uncert->SetFillColor(kBlack);
    h1_uncert->SetLineColor(kBlack);
    h1_uncert->SetFillStyle(3005);

    uncert_rootfile->Close();
}

//
//
//
float Chi2(TH1F* h_data, TH1F* h_mc, bool dochi2ndf=true) {
    
    float chi2  = 0.; 
    int NDF     = h_data->GetXaxis()->GetNbins(); 
    for(int i=1; i<h_data->GetXaxis()->GetNbins()+1; i++) {

        float diff2 = (h_data->GetBinContent(i) - h_mc->GetBinContent(i))
                     *(h_data->GetBinContent(i) - h_mc->GetBinContent(i));
        float err2  = h_data->GetBinError(i)*h_data->GetBinError(i)
                     +h_mc->GetBinError(i)*h_mc->GetBinError(i);
        float chi2_tmp = diff2/err2; 

        chi2 = chi2 + chi2_tmp; 
    
        //cout << i << " bin ::: chi2/NDF : " << chi2_tmp << "/" << NDF << endl; 
        
    }

    if(dochi2ndf) return chi2/NDF; 
    else chi2;
}

//
// main function
//
void compareshapes(bool doreweightww=true, char* plotregion="CR1") {
  
    //bool doreweightww = true;
    int njet = 0;
    char* binsize       = "15-25 GeV";
    char* flavor        = "of";
    char* energy        = "8tev"; 
    //char* plotregion    = "CR2";
    char* normregion    = "CR2";
    if(plotregion=="CR2") normregion = "CR1";
   

    cout << "** plot region : " <<  plotregion << endl;
    cout << "** norm region : " <<  normregion << endl;

    //TFile* postFile_default = TFile::Open("output/fit_default_new1-fit-125-all_fittedShape_floatMu.root"); 
    //TFile* postFile_fornorm = TFile::Open(Form("output/fit_%s-fit-125-all_fittedShape_floatMu.root", normregion)); 
    //TFile* rootfile_default = TFile::Open("125/hwwof_0j.input_8TeV.root");
    TFile* postFile_default = TFile::Open("newdefault_fittedShape_floatMu.root"); 
    TFile* postFile_fornorm = TFile::Open(Form("%s_fittedShape_floatMu.root", normregion)); 
    TFile* rootfile_default = TFile::Open("125/hwwof_0j.input_8TeV.root");

    TH1F *h1post_qqWW, *h1post_qqWW_fornorm;
    TH1F *h1_mll_qqWW, *h1_mT_qqWW;
    TH1F *h1post_ZH, *h1post_WH, *h1post_qqH, *h1post_ggH, 
         *h1post_ggWW, *h1post_VV, *h1post_Top, 
         *h1post_Zjets, *h1post_WjetsE, *h1post_WjetsM, 
         *h1post_Wgamma, *h1post_Wg3l, *h1post_Ztt, 
         *h1_data;  
    TH2F *h2post_ZH, *h2post_WH, *h2post_qqH, *h2post_ggH, 
         *h2post_ggWW, *h2post_VV, *h2post_Top, 
         *h2post_Zjets, *h2post_WjetsE, *h2post_WjetsM, 
         *h2post_Wgamma, *h2post_Wg3l, *h2post_Ztt, 
         *h2_data;  
    TH1F *h1_mll_ZH, *h1_mll_WH, *h1_mll_qqH, *h1_mll_ggH, 
         *h1_mll_ggWW, *h1_mll_VV, *h1_mll_Top, 
         *h1_mll_Zjets, *h1_mll_WjetsE, *h1_mll_WjetsM, 
         *h1_mll_Wgamma, *h1_mll_Wg3l, *h1_mll_Ztt, 
         *h1_mll_data;  
    TH1F *h1_mT_ZH, *h1_mT_WH, *h1_mT_qqH, *h1_mT_ggH, 
         *h1_mT_ggWW, *h1_mT_VV, *h1_mT_Top, 
         *h1_mT_Zjets, *h1_mT_WjetsE, *h1_mT_WjetsM, 
         *h1_mT_Wgamma, *h1_mT_Wg3l, *h1_mT_Ztt, 
         *h1_mT_data;  

    // get qqWW for nominal and control fit 
    h1post_qqWW_default     = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_qqWW", njet, flavor, "")));   
    h1post_qqWW_fornorm     = makeHistogram( (TGraphAsymmErrors*) postFile_fornorm->Get(Form("j%i%s%s_qqWW", njet, flavor, "")));   
    
    // get other processes from nominal fit 
    h1post_ZH       = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_ZH", njet, flavor, "")));   
    h1post_WH       = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_WH", njet, flavor, "")));   
    h1post_qqH      = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_qqH", njet, flavor, "")));   
    h1post_ggH      = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_ggH", njet, flavor, "")));   
    h1post_qqWW     = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_qqWW", njet, flavor, "")));   
    h1post_ggWW     = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_ggWW", njet, flavor, "")));   
    h1post_VV       = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_VV", njet, flavor, "")));   
    h1post_Top      = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_Top", njet, flavor, "")));   
    h1post_Zjets    = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_Zjets", njet, flavor, "")));   
    h1post_WjetsE   = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_WjetsE", njet, flavor, "")));   
    h1post_WjetsM   = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_WjetsM", njet, flavor, "")));   
    h1post_Wgamma   = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_Wgamma", njet, flavor, "")));   
    h1post_Wg3l     = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_Wg3l", njet, flavor, "")));   
    h1post_Ztt      = makeHistogram( (TGraphAsymmErrors*) postFile_default->Get(Form("j%i%s%s_Ztt", njet, flavor, "")));   
    
    // get data
    h1_data         = (TH1F*)rootfile_default->Get("histo_Data"); 

    // make 2d  
    h2post_qqWW         = Roll1DTo2D(h1post_qqWW_default); 
    h2post_qqWW_fornorm = Roll1DTo2D(h1post_qqWW_fornorm); 

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
    
    h2_data         = Roll1DTo2D(h1_data); 

    //
    // renormalize qqWW using the other control regions  
    //
    if(plotregion == "CR1") {
        float qqww_normCR2  = h2post_qqWW_fornorm->Integral(1,6,4,9); 
        float qqww          = h2post_qqWW->Integral(1,6,4,9); 
        cout << "scale in CR2 : " << qqww_normCR2 << " / " << qqww << endl; 
        if(doreweightww) h2post_qqWW->Scale(qqww_normCR2/qqww);
    } else { 
        float qqww_normCR1  = h2post_qqWW_fornorm->Integral(7,14,1,9);  
        float qqww          = h2post_qqWW->Integral(7,14,1,9); 
        cout << "scale in CR1 : " << qqww_normCR1 << " / " << qqww << endl; 
        if(doreweightww) h2post_qqWW->Scale(qqww_normCR1/qqww);
    }

    //
    // project all process to mT and mll in the control region
    //
    if(plotregion == "CR1") {
        h1_mll_data     = doProjection(h2_data,         "mll",  "CR1");   
        h1_mT_data      = doProjection(h2_data,         "mT",   "CR1");    
        h1_mll_qqWW     = doProjection(h2post_qqWW,     "mll",  "CR1"); 
        h1_mT_qqWW      = doProjection(h2post_qqWW,     "mT",   "CR1"); 
       
        h1_mll_ZH       = doProjection(h2post_ZH,       "mll",  "CR1"); 
        h1_mT_ZH        = doProjection(h2post_ZH,       "mT",   "CR1"); 
        h1_mll_WH       = doProjection(h2post_WH,       "mll",  "CR1"); 
        h1_mT_WH        = doProjection(h2post_WH,       "mT",   "CR1"); 
        h1_mll_qqH      = doProjection(h2post_qqH,      "mll",  "CR1"); 
        h1_mT_qqH       = doProjection(h2post_qqH,      "mT",   "CR1"); 
        h1_mll_ggH      = doProjection(h2post_ggH,      "mll",  "CR1"); 
        h1_mT_ggH       = doProjection(h2post_ggH,      "mT",   "CR1"); 
        h1_mll_ggWW     = doProjection(h2post_ggWW,     "mll",  "CR1"); 
        h1_mT_ggWW      = doProjection(h2post_ggWW,     "mT",   "CR1"); 
        h1_mll_VV       = doProjection(h2post_VV,       "mll",  "CR1"); 
        h1_mT_VV        = doProjection(h2post_VV,       "mT",   "CR1"); 
        h1_mll_Top      = doProjection(h2post_Top,      "mll",  "CR1"); 
        h1_mT_Top       = doProjection(h2post_Top,      "mT",   "CR1"); 
        h1_mll_Zjets    = doProjection(h2post_Zjets,    "mll",  "CR1"); 
        h1_mT_Zjets     = doProjection(h2post_Zjets,    "mT",   "CR1"); 
        h1_mll_WjetsE   = doProjection(h2post_WjetsE,   "mll",  "CR1"); 
        h1_mT_WjetsE    = doProjection(h2post_WjetsE,   "mT",   "CR1"); 
        h1_mll_WjetsM   = doProjection(h2post_WjetsM,   "mll",  "CR1"); 
        h1_mT_WjetsM    = doProjection(h2post_WjetsM,   "mT",   "CR1"); 
        h1_mll_Wgamma   = doProjection(h2post_Wgamma,   "mll",  "CR1"); 
        h1_mT_Wgamma    = doProjection(h2post_Wgamma,   "mT",   "CR1"); 
        h1_mll_Wg3l     = doProjection(h2post_Wg3l,     "mll",  "CR1"); 
        h1_mT_Wg3l      = doProjection(h2post_Wg3l,     "mT",   "CR1"); 
        h1_mll_Ztt      = doProjection(h2post_Ztt,      "mll",  "CR1"); 
        h1_mT_Ztt       = doProjection(h2post_Ztt,      "mT",   "CR1"); 
    } else { 
        h1_mll_data     = doProjection(h2_data,         "mll",  "CR2");   
        h1_mT_data      = doProjection(h2_data,         "mT",   "CR2");   
        h1_mll_qqWW     = doProjection(h2post_qqWW,     "mll",  "CR2"); 
        h1_mT_qqWW      = doProjection(h2post_qqWW,     "mT",   "CR2"); 
       
        h1_mll_ZH       = doProjection(h2post_ZH,       "mll",  "CR2"); 
        h1_mT_ZH        = doProjection(h2post_ZH,       "mT",   "CR2"); 
        h1_mll_WH       = doProjection(h2post_WH,       "mll",  "CR2"); 
        h1_mT_WH        = doProjection(h2post_WH,       "mT",   "CR2"); 
        h1_mll_qqH      = doProjection(h2post_qqH,      "mll",  "CR2"); 
        h1_mT_qqH       = doProjection(h2post_qqH,      "mT",   "CR2"); 
        h1_mll_ggH      = doProjection(h2post_ggH,      "mll",  "CR2"); 
        h1_mT_ggH       = doProjection(h2post_ggH,      "mT",   "CR2"); 
        h1_mll_ggWW     = doProjection(h2post_ggWW,     "mll",  "CR2"); 
        h1_mT_ggWW      = doProjection(h2post_ggWW,     "mT",   "CR2"); 
        h1_mll_VV       = doProjection(h2post_VV,       "mll",  "CR2"); 
        h1_mT_VV        = doProjection(h2post_VV,       "mT",   "CR2"); 
        h1_mll_Top      = doProjection(h2post_Top,      "mll",  "CR2"); 
        h1_mT_Top       = doProjection(h2post_Top,      "mT",   "CR2"); 
        h1_mll_Zjets    = doProjection(h2post_Zjets,    "mll",  "CR2"); 
        h1_mT_Zjets     = doProjection(h2post_Zjets,    "mT",   "CR2"); 
        h1_mll_WjetsE   = doProjection(h2post_WjetsE,   "mll",  "CR2"); 
        h1_mT_WjetsE    = doProjection(h2post_WjetsE,   "mT",   "CR2"); 
        h1_mll_WjetsM   = doProjection(h2post_WjetsM,   "mll",  "CR2"); 
        h1_mT_WjetsM    = doProjection(h2post_WjetsM,   "mT",   "CR2"); 
        h1_mll_Wgamma   = doProjection(h2post_Wgamma,   "mll",  "CR2"); 
        h1_mT_Wgamma    = doProjection(h2post_Wgamma,   "mT",   "CR2"); 
        h1_mll_Wg3l     = doProjection(h2post_Wg3l,     "mll",  "CR2"); 
        h1_mT_Wg3l      = doProjection(h2post_Wg3l,     "mT",   "CR2"); 
        h1_mll_Ztt      = doProjection(h2post_Ztt,      "mll",  "CR2"); 
        h1_mT_Ztt       = doProjection(h2post_Ztt,      "mT",   "CR2"); 
    } 

    if(0) {
        TCanvas *c_debug_mll = new TCanvas("c_debug_mll", "c_debug_mll", 800, 800);
        c_debug_mll->Divide(4,4);
        c_debug_mll->cd(1); h1_mll_ZH->Draw("histo"); 
        c_debug_mll->cd(2); h1_mll_WH->Draw("histo"); 
        c_debug_mll->cd(3); h1_mll_qqH->Draw("histo"); 
        c_debug_mll->cd(4); h1_mll_ggH->Draw("histo"); 
        c_debug_mll->cd(5); h1_mll_qqWW->Draw("histo"); 
        c_debug_mll->cd(6); h1_mll_ggWW->Draw("histo"); 
        c_debug_mll->cd(7); h1_mll_VV->Draw("histo"); 
        c_debug_mll->cd(8); h1_mll_Top->Draw("histo"); 
        c_debug_mll->cd(9); h1_mll_Zjets->Draw("histo"); 
        c_debug_mll->cd(10); h1_mll_WjetsE->Draw("histo"); 
        c_debug_mll->cd(11); h1_mll_WjetsM->Draw("histo"); 
        c_debug_mll->cd(12); h1_mll_Wgamma->Draw("histo"); 
        c_debug_mll->cd(13); h1_mll_Wg3l->Draw("histo"); 
        c_debug_mll->cd(14); h1_mll_Ztt->Draw("histo"); 
        c_debug_mll->cd(15); h1_mll_data->Draw("e"); 
        TCanvas *c_debug_mT = new TCanvas("c_debug_mT", "c_debug_mT", 800, 800);
        c_debug_mT->Divide(4,4);
        c_debug_mT->cd(1); h1_mT_ZH->Draw("histo"); 
        c_debug_mT->cd(2); h1_mT_WH->Draw("histo"); 
        c_debug_mT->cd(3); h1_mT_qqH->Draw("histo"); 
        c_debug_mT->cd(4); h1_mT_ggH->Draw("histo"); 
        c_debug_mT->cd(5); h1_mT_qqWW->Draw("histo"); 
        c_debug_mT->cd(6); h1_mT_ggWW->Draw("histo"); 
        c_debug_mT->cd(7); h1_mT_VV->Draw("histo"); 
        c_debug_mT->cd(8); h1_mT_Top->Draw("histo"); 
        c_debug_mT->cd(9); h1_mT_Zjets->Draw("histo"); 
        c_debug_mT->cd(10); h1_mT_WjetsE->Draw("histo"); 
        c_debug_mT->cd(11); h1_mT_WjetsM->Draw("histo"); 
        c_debug_mT->cd(12); h1_mT_Wgamma->Draw("histo"); 
        c_debug_mT->cd(13); h1_mT_Wg3l->Draw("histo"); 
        c_debug_mT->cd(14); h1_mT_Ztt->Draw("histo"); 
        c_debug_mT->cd(15); h1_mT_data->Draw("e"); 
    }

    // reorganize background processes 
    // sig = ZH, WH, qqH, ggH, 
    // WW = qqWW + ggWW
    // VVWgamma = VV + Wgamma + Wg3l
    // Top = Top
    // Wjets = WjetsE + WjetsM
    TH1F *h1_mll_sig, *h1_mll_WW, *h1_mll_VVWgamma, *h1_mll_Wjets;
    TH1F *h1_mT_sig, *h1_mT_WW, *h1_mT_VVWgamma, *h1_mT_Wjets;
    reorganizebkg(  h1_mll_ZH, h1_mll_WH, h1_mll_qqH, h1_mll_ggH, 
                    h1_mll_qqWW, h1_mll_ggWW, h1_mll_VV, h1_mll_Top, 
                    h1_mll_Zjets, h1_mll_WjetsE, h1_mll_WjetsM, 
                    h1_mll_Wgamma, h1_mll_Wg3l, h1_mll_Ztt, 
                    h1_mll_sig, h1_mll_WW, h1_mll_VVWgamma, h1_mll_Wjets,
                    "M_{ll} GeV");
    reorganizebkg(  h1_mT_ZH, h1_mT_WH, h1_mT_qqH, h1_mT_ggH, 
                    h1_mT_qqWW, h1_mT_ggWW, h1_mT_VV, h1_mT_Top, 
                    h1_mT_Zjets, h1_mT_WjetsE, h1_mT_WjetsM, 
                    h1_mT_Wgamma, h1_mT_Wg3l, h1_mT_Ztt, 
                    h1_mT_sig, h1_mT_WW, h1_mT_VVWgamma, h1_mT_Wjets,
                    "M_{T} GeV");

    // 
    //  Uncertainty bands
    // 

    TH1F *h1_mll_uncert, *h1_mll_bkg;        // For uncert bands
    TH1F *h1_mT_uncert, *h1_mT_bkg;          // Form uncert bands
    h1_mll_bkg = (TH1F*)h1_mll_Ztt->Clone("h1_mll_bkg");  h1_mll_bkg->SetName("h1_mll_bkg"); h1_mll_bkg->SetTitle("h1_mll_bkg");
    h1_mll_bkg->Add(h1_mll_Top);    h1_mll_bkg->Add(h1_mll_VVWgamma);   h1_mll_bkg->Add(h1_mll_Wjets);
    h1_mT_bkg = (TH1F*)h1_mT_Ztt->Clone("h1_mT_bkg");  h1_mT_bkg->SetName("h1_mT_bkg"); h1_mT_bkg->SetTitle("h1_mT_bkg");
    h1_mT_bkg->Add(h1_mT_Top);  h1_mT_bkg->Add(h1_mT_VVWgamma); h1_mT_bkg->Add(h1_mT_Wjets);
    h1_mll_uncert = (TH1F*)h1_mll_Ztt->Clone("h1_mll_uncert");  h1_mll_uncert->SetName("h1_mll_uncert"); h1_mll_uncert->SetTitle("h1_mll_uncert");
    h1_mll_uncert->Add(h1_mll_Top);  h1_mll_uncert->Add(h1_mll_VVWgamma); h1_mll_uncert->Add(h1_mll_Wjets);  h1_mll_uncert->Add(h1_mll_WW); h1_mll_uncert->Add(h1_mll_sig);
    h1_mT_uncert = (TH1F*)h1_mT_Ztt->Clone("h1_mT_uncert");  h1_mT_uncert->SetName("h1_mT_uncert"); h1_mT_uncert->SetTitle("h1_mT_uncert");
    h1_mT_uncert->Add(h1_mT_Top);  h1_mT_uncert->Add(h1_mT_VVWgamma); h1_mT_uncert->Add(h1_mT_Wjets);  h1_mT_uncert->Add(h1_mT_WW); h1_mT_uncert->Add(h1_mT_sig);

    getUncerBand(h1_mll_uncert, h1_mll_WW, h1_mll_bkg, "mll", plotregion);
    getUncerBand(h1_mT_uncert, h1_mT_WW, h1_mT_bkg, "mT", plotregion);

    //
    // Legend 
    //

    TLegend* l1 = new TLegend (.16,.88-6*.035,.5,.88);
    l1->SetNColumns(2);
    l1->SetBorderSize(0);
    l1->SetFillColor(0);
    l1->SetFillStyle(0);
    l1->SetTextFont(42); 
    l1->SetTextAlign(12);
    l1->SetTextSize(0.03);
    
    l1->AddEntry(h1_mll_data,       "data",         "lp");
    l1->AddEntry(h1_mll_sig,        "M_{H}=125 GeV", "f");
    l1->AddEntry(h1_mll_WW,         "WW",           "f");
    l1->AddEntry(h1_mll_Wjets,      "W+jets",       "f");
    l1->AddEntry(h1_mll_VVWgamma,   "WZ/ZZ",        "f");
    l1->AddEntry(h1_mll_Top,        "Top",          "f");
    l1->AddEntry(h1_mll_Ztt,        "Z+jets",       "f");
    l1->AddEntry(h1_mll_uncert,     "Syst. uncert.","f");
   
    float textSize = 0.03;
    int textfont = 62; // default font is 62
    TLatex *tex_cmsprel = new TLatex(0.5,0.84,"CMS Preliminary");
    tex_cmsprel->SetNDC();
    tex_cmsprel->SetTextSize(textSize+0.02);
    tex_cmsprel->SetLineWidth(2);
    //TLatex *tex_cms = new TLatex(0.5,0.84,"CMS, ");
    //tex_cms->SetNDC();
    //tex_cms->SetTextSize(textSize+0.02);
    //tex_cms->SetLineWidth(2);
    TLatex *tex_8tev = new TLatex(0.6,0.84,Form("#font[%i]{#sqrt{s}=8 TeV, L = %.1f fb^{-1}}", textfont, 19.5));
    tex_8tev->SetNDC();
    tex_8tev->SetTextSize(textSize);
    tex_8tev->SetLineWidth(2);
    TLatex *tex_ch = new TLatex(0.5,0.78, Form("#font[%i]{Category : %s, %s}", textfont, "0jet", "e#mu/#mue"));
    tex_ch->SetNDC();
    tex_ch->SetTextSize(textSize);
    tex_ch->SetLineWidth(2);
    TLatex *tex_plotregion = new TLatex(0.5,0.73, plotregion=="CR1"?
                                                   Form("#font[%i]{120<M_{T}<280 GeV, 12<M_{ll}<200 GeV}",textfont)
                                                  :Form("#font[%i]{60<M_{T}<120 GeV, 60<M_{ll}<200 GeV}",textfont));
    tex_plotregion->SetNDC();
    tex_plotregion->SetTextSize(textSize);
    tex_plotregion->SetLineWidth(2);
    
    TLatex *tex_chi2ndf_mll = new TLatex(0.5,0.68, Form("#chi^{2}/n.d.f. = %.2f/%i = %.2f", 
                                                        Chi2(h1_mll_data, h1_mll_uncert, false), h1_mll_data->GetXaxis()->GetNbins(),
                                                        Chi2(h1_mll_data, h1_mll_uncert, false)/(float)h1_mll_data->GetXaxis()->GetNbins()));
    tex_chi2ndf_mll->SetNDC();
    tex_chi2ndf_mll->SetTextSize(textSize);
    tex_chi2ndf_mll->SetLineWidth(2);
    TLatex *tex_chi2ndf_mT = new TLatex(0.5,0.68, Form("#chi^{2}/n.d.f. = %.2f/%i = %.2f", 
                                                        Chi2(h1_mT_data, h1_mT_uncert, false), h1_mT_data->GetXaxis()->GetNbins(),
                                                        Chi2(h1_mT_data, h1_mT_uncert, false)/(float)h1_mT_data->GetXaxis()->GetNbins()));
    tex_chi2ndf_mT->SetNDC();
    tex_chi2ndf_mT->SetTextSize(textSize);
    tex_chi2ndf_mT->SetLineWidth(2);

    h1_mll_data->SetMarkerStyle(20);   
    h1_mll_data->SetMarkerSize(1.2);   
    h1_mll_data->SetLineWidth(1.2); 
    h1_mll_data->SetLineColor(kBlack);
    h1_mT_data->SetMarkerStyle(20);    
    h1_mT_data->SetMarkerSize(1.2);   
    h1_mT_data->SetLineWidth(1.2);  
    h1_mT_data->SetLineColor(kBlack);

    // make stack of backgrounds   
    TCanvas *c_mll = new TCanvas("mll","mll",600,500);
    c_mll->cd();
    THStack *st_mll = new THStack(Form("st_mll_%s", plotregion), Form(";%s;events / %s", "M_{ll} (GeV)", binsize));
    st_mll->Add(h1_mll_sig);
    st_mll->Add(h1_mll_Ztt);
    st_mll->Add(h1_mll_Top);
    st_mll->Add(h1_mll_VVWgamma);
    st_mll->Add(h1_mll_Wjets);
    st_mll->Add(h1_mll_WW);
    st_mll->Draw("hist");
    h1_mll_data->Draw("e same");  
    h1_mll_uncert->Draw("e2 same");
    st_mll->SetMaximum(h1_mll_data->GetMaximum()*2.0);
    st_mll->GetYaxis()->SetTitleOffset(1.4);
    l1->Draw(); 
    tex_cmsprel->Draw("SAME");
    //tex_cms->Draw("SAME");
    //tex_8tev->Draw("SAME");
    tex_ch->Draw("SAME");
    tex_plotregion->Draw("SAME"); 
    tex_chi2ndf_mll->Draw("SAME");
    c_mll->Print(Form("mll_%s_%s.pdf", plotregion, doreweightww?"weightedqqww":"unweightedqqww")); 

    TCanvas *c_mT = new TCanvas("mT","mT",600,500);
    c_mT->cd();
    if(plotregion == "CR1") binsize = "20 GeV";
    if(plotregion == "CR2") binsize = "10 GeV";
    THStack *st_mT = new THStack(Form("st_mT_%s", plotregion), Form(";%s;events / %s", "M_{T} (GeV)", binsize));
    st_mT->Add(h1_mT_sig);
    st_mT->Add(h1_mT_Ztt);
    st_mT->Add(h1_mT_Top);
    st_mT->Add(h1_mT_VVWgamma);
    st_mT->Add(h1_mT_Wjets);
    st_mT->Add(h1_mT_WW);
    st_mT->Draw("hist");
    h1_mT_data->Draw("e same"); 
    h1_mT_uncert->Draw("e2 same");
    st_mT->SetMaximum(h1_mT_data->GetMaximum()*2.0);
    st_mT->GetYaxis()->SetTitleOffset(1.4);
    l1->Draw(); 
    tex_cmsprel->Draw("SAME");
    //tex_cms->Draw("SAME");
    //tex_8tev->Draw("SAME");
    tex_ch->Draw("SAME");
    tex_plotregion->Draw("SAME");
    tex_chi2ndf_mT->Draw("SAME");

    c_mT->Print(Form("mT_%s_%s.pdf", plotregion, doreweightww?"weightedqqww":"unweightedqqww")); 

//    delete c_mT;
//    delete c_mll; 
}

void compareshapes_runAll() {

    compareshapes(true, "CR1");
    compareshapes(true, "CR2");
//    compareshapes(false, "CR1");
//    compareshapes(false, "CR2");
}
