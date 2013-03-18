//
// print a line 
//
void printtableline(TString process, float Npre, float Npost) {
    
    cout.width(15); cout << process << " & ";   
    cout.width(10); cout << Form("%5.1f", Npre) << " & ";
    cout.width(10); cout << Form("%5.1f", (Npre ? Npost : 0.0)) << " & ";
    cout.width(10); cout << Form("%5.1f", (Npre ? (Npost - Npre) : 0.0)) << " & "; 
    cout.width(10); cout << Form("%5.1f", (Npre ? (Npost - Npre)/Npre*100 : 0.0)); 
    cout.width(10); cout << "\\\\" << endl;
}

//
// Zoom
//
TH2D* Zoom(TH2D* h2) {

    float mtbins[7]    =   {60,70,80,90,100,110,120};
    float mllbins[6]   =   {12,30,45,60,75,100};

    TH2D* h2zoom = new TH2D(Form("%s_zoom", h2->GetName()), Form("%s_zoom", h2->GetTitle()), 6, mtbins, 5, mllbins);

    for(unsigned int ibinX=1; ibinX <= 6; ++ibinX) {
        for(unsigned int ibinY=1; ibinY <= 5; ++ibinY) {
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
// Convert TGraphAsymmErrors to TH1D
// 
TH1D* makeHistogram(TGraphAsymmErrors* graph){
  
    //TH1D* hist = new TH1D(Form("h_%s",graph->GetName()),Form("%s",graph->GetName()),126,-1,1);
    TH1D* hist = new TH1D(Form("h_%s", graph->GetName()), Form("%s",graph->GetName()), 126, 0.5, 126.5);
    hist->SetDirectory(0); 
    for (int i=0; i<graph->GetN();++i){
        Int_t bin = hist->FindBin(graph->GetX()[i]);
        hist->SetBinContent(bin,graph->GetY()[i]);  
    }
    hist->SetStats(0); 
    return hist;
}

//
// B/W color scheme for the color-blind
// from Dave who got from Dima
//
void loadColorPalette() {

    const int nColor = 6;
    Int_t palette[nColor];
    palette[0] = kWhite;
    for (unsigned int i=1;i<nColor;++i){
        palette[i] = 19-i;
    }
    gStyle->SetPalette(nColor,palette);

}

//
// Cosmetics for projected histograms 
//
TH1D * h1cosmetic(TH1D *h1, char* title="", int linecolor=kRed, int linewidth=1, int fillcolor=0){
    
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetStats(0);
    h1->SetMinimum(0.);

    return h1;
}

//
// rolling 1D to 2D 
//
TH2D* Roll1DTo2D(TH1D *h1, bool doZoom=true) {

    // style
    gStyle->SetPaintTextFormat(".1f");
    gStyle->SetMarkerSize(2.5);

    // palette
    loadColorPalette();

    unsigned int nbins  = h1->GetXaxis()->GetNbins();
    
    float mtbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
    float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};
    TH2D* h2 = new TH2D(h1->GetName(), h1->GetTitle(), 14, mtbins, 9, mllbins);
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
    
    return h2;
}

void drawoneprocess(TH2D *h2pre, TH2D *h2post, TString prefix, bool doZoom=false) {
    
    TString title = h2pre->GetTitle();
    title.Remove(0,6);
    
    TH2D *h2pre_draw, *h2post_draw;

    if(doZoom) {
        h2pre_draw  = Zoom(h2pre);
        h2post_draw = Zoom(h2post);
    } else { 
        h2pre_draw  = (TH2D*)h2pre->Clone();
        h2post_draw = (TH2D*)h2post->Clone();
    }

    // redraw missing lines due to white bins    
    TGaxis *axistmp;
    if (!doZoom) axistmp = new TGaxis(60.,200.,280.,200.,0, 2,50510,"+L");
    else if(doZoom) axistmp = new TGaxis(60.,100.,120.,100.,0, 2,50510,"+L");
    axistmp->SetNdivisions(0);
    axistmp->SetLabelFont(0);
   
    // ratio plot
    TH2D *h2ratio = h2post_draw->Clone("h2ratio");
    //h2ratio->Add(h2pre, -1);    
    h2ratio->Divide(h2pre_draw); 
    h2ratio->SetTitle("Ratio"); 

    // projection to mT
    TH1D *h1pre_mT  = h1cosmetic( h2pre_draw->ProjectionX("h1pre_mT",1,1), Form("%s - M_{T} projection", title.Data()), kRed, 2);
    TH1D *h1post_mT = h1cosmetic( h2post_draw->ProjectionX("h1post_mT",1,1), Form("%s - M_{T} projection", title.Data()), kBlue, 2);
    // projection to mll
    TH1D *h1pre_mll  = h1cosmetic( h2pre_draw->ProjectionY("h1pre_mll"), Form("%s - M_{ll} projection", title.Data()), kRed, 2);
    TH1D *h1post_mll = h1cosmetic( h2post_draw->ProjectionY("h1post_mll"), Form("%s - M_{ll} projection", title.Data()), kBlue, 2);

    // Canvas 
    TCanvas *c = new TCanvas("c", "", 1500, 1000);
    c->Divide(3,2);
    
    c->cd(1); 
    h2pre_draw->SetTitle(Form("%s - prefit", title.Data()));
    if(!doZoom) h2pre_draw->Draw("colz");
    else if(doZoom) h2pre_draw->Draw("colz text");
    axistmp->Draw();
    
    c->cd(2); 
    h2post_draw->SetTitle(Form("%s - postfit", title.Data()));
    if(!doZoom) h2post_draw->Draw("colz");
    else if(doZoom) h2post_draw->Draw("colz text");
    axistmp->Draw();
    
    c->cd(3);  // ratio = post / pre
    if(!doZoom) h2ratio->Draw("colz");
    else if(doZoom) h2ratio->Draw("colz text");
    axistmp->Draw();

    c->cd(4); 
    h1pre_mT->SetMaximum( TMath::Max( h1pre_mT->GetBinContent(h1pre_mT->GetMaximumBin()), h1post_mT->GetBinContent(h1post_mT->GetMaximumBin()) )*1.2 );
    h1pre_mT->Draw("histo");
    h1post_mT->Draw("histo same");

    c->cd(5); 
    h1pre_mll->SetMaximum( TMath::Max( h1pre_mll->GetBinContent(h1pre_mll->GetMaximumBin()), h1post_mll->GetBinContent(h1post_mll->GetMaximumBin()) )*1.2 );
    h1pre_mll->Draw("histo");
    h1post_mll->Draw("histo same");

    if(!doZoom) c->Print(Form("prepostplots/%spre-postfit-%s.pdf",      prefix.Data(),  title.Data())); 
    if(doZoom)  c->Print(Form("prepostplots/%spre-postfit-%s_zoom.pdf", prefix.Data(),  title.Data())); 

    delete c;
    delete h2pre_draw; 
    delete h2post_draw;
}

void drawstack (
    TH1D *h1_ZH, TH1D *h1_WH, TH1D *h1_qqH, TH1D *h1_ggH, TH1D *h1_qqWW, TH1D *h1_ggWW, TH1D *h1_VV, TH1D *h1_Top, 
    TH1D *h1_Zjets, TH1D *h1_WjetsE, TH1D *h1_WjetsM, TH1D *h1_Wgamma, TH1D *h1_Wg3l, TH1D *h1_Ztt, TH1D *h1_data, 
    bool drawmll=false
) {

    TH1D *h1tmp_qqWW	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_qqWW)->ProjectionY() : Roll1DTo2D(h1_qqWW)->ProjectionX(),    "qqWW",     1,1,851); 
    TH1D *h1tmp_ggWW	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_ggWW)->ProjectionY() : Roll1DTo2D(h1_ggWW)->ProjectionX(),    "ggWW",     1,1,851); 
    TH1D *h1tmp_VV		=   h1cosmetic( drawmll ? Roll1DTo2D(h1_VV)->ProjectionY() : Roll1DTo2D(h1_VV)->ProjectionX(),    "VV",     1,1,418); //858
    TH1D *h1tmp_Top		=   h1cosmetic( drawmll ? Roll1DTo2D(h1_Top)->ProjectionY() : Roll1DTo2D(h1_Top)->ProjectionX(),    "Top",     1,1,400); 
    TH1D *h1tmp_Wgamma	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_Wgamma)->ProjectionY() : Roll1DTo2D(h1_Wgamma)->ProjectionX(),    "Wgamma",     1,1,858); 
    TH1D *h1tmp_Wg3l	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_Wg3l)->ProjectionY() : Roll1DTo2D(h1_Wg3l)->ProjectionX(),    "Wg3l",     1,1,858); 
    TH1D *h1tmp_WjetsE	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_WjetsE)->ProjectionY() : Roll1DTo2D(h1_WjetsE)->ProjectionX(),    "WjetsE",     1,1,921); 
    TH1D *h1tmp_WjetsM	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_WjetsM)->ProjectionY() : Roll1DTo2D(h1_WjetsM)->ProjectionX(),    "WjetsM",     1,1,921); 
    //TH1D *h1tmp_Zjets	=   h1cosmetic( drawmll ? Roll1DTo2D(h1_Zjets)->ProjectionY() : Roll1DTo2D(h1_Zjets)->ProjectionX(),    "Zjets",     1,1,418); 
    TH1D *h1tmp_Ztt		=   h1cosmetic( drawmll ? Roll1DTo2D(h1_Ztt)->ProjectionY() : Roll1DTo2D(h1_Ztt)->ProjectionX(),    "Ztt",     1,1,418); 
    TH1D *h1tmp_ZH      =   h1cosmetic( drawmll ? Roll1DTo2D(h1_ZH)->ProjectionY() : Roll1DTo2D(h1_ZH)->ProjectionX(),    "ZH",     1,1,2);
    TH1D *h1tmp_WH		=   h1cosmetic( drawmll ? Roll1DTo2D(h1_WH)->ProjectionY() : Roll1DTo2D(h1_WH)->ProjectionX(),    "WH",     1,1,2); 
    TH1D *h1tmp_qqH		=   h1cosmetic( drawmll ? Roll1DTo2D(h1_qqH)->ProjectionY() : Roll1DTo2D(h1_qqH)->ProjectionX(),    "qqH",     1,1,2); 
    TH1D *h1tmp_ggH		=   h1cosmetic( drawmll ? Roll1DTo2D(h1_ggH)->ProjectionY() : Roll1DTo2D(h1_ggH)->ProjectionX(),    "ggH",     1,1,2); 
    TH1D *h1tmp_data    =   h1cosmetic( drawmll ? Roll1DTo2D(h1_data)->ProjectionY() : Roll1DTo2D(h1_data)->ProjectionX(),    "Data",     1,1,2); 


    // Add process 
    TH1D *h1tmp_WW      = (TH1D*)h1tmp_qqWW->Clone("h1tmp_WW");         h1tmp_WW->Add(h1tmp_ggWW);
    TH1D *h1tmp_Wg      = (TH1D*)h1tmp_Wgamma->Clone("h1tmp_Wg");       h1tmp_Wg->Add(h1tmp_Wg3l);
    TH1D *h1tmp_Wjets   = (TH1D*)h1tmp_WjetsE->Clone("h1tmp_Wjets");    h1tmp_Wjets->Add(h1tmp_WjetsM);
    TH1D *h1tmp_sig     = (TH1D*)h1tmp_ggH->Clone("h1tmp_sig");         h1tmp_sig->Add(h1tmp_qqH); h1tmp_sig->Add(h1tmp_WH); h1tmp_sig->Add(h1tmp_ZH);

    THStack* hstack  = new THStack("hstack", "Signal+background fit");
    hstack->Add( h1tmp_WW);
    hstack->Add( h1tmp_VV);
    hstack->Add( h1tmp_Top);
    hstack->Add( h1tmp_Wg);
    hstack->Add( h1tmp_Wjets);
    //hstack->Add( h1tmp_Zjets);
    hstack->Add( h1tmp_Ztt);
    hstack->Add( h1tmp_sig);

   // if(drawmll) hstack->GetXaxis()->SetTitle("M_{ll} (GeV)");
   // if(!drawmll) hstack->GetXaxis()->SetTitle("M_{T} (GeV)");
    
    h1tmp_data->SetMarkerStyle(20); 
    h1tmp_data->SetMarkerSize(1); 

    hstack->SetMaximum(h1tmp_data->GetMaximum()*1.2); 
    hstack->Draw("hist"); 
    h1tmp_data->Draw("PE same"); 

    // Legend 
    //TLegend* legend = new TLegend(x1, y1, x1 + _xoffset, y1 + _yoffset);
    TLegend* legend = new TLegend(0.65, 0.5, 0.88 , 0.85);
    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (    0.03);
    legend->AddEntry(h1tmp_sig,     "MH=125 GeV",   "f");
    legend->AddEntry(h1tmp_WW,      "WW",           "f");
    legend->AddEntry(h1tmp_Top,     "Top",          "f");
    legend->AddEntry(h1tmp_VV,      "VV",           "f");
    legend->AddEntry(h1tmp_Wjets,   "Wjets",        "f");
    legend->AddEntry(h1tmp_Wg,      "W#gamma(*)",   "f");
    legend->Draw(); 

    
}

//
// main function
//
void printNorm_lands(char* carddir="ana_Moriond13_2D_V1", int mh=125, char* flavor="of", int njet=0, char* energy="8TeV", bool doZoom=true) {

    // print signal or not 
    bool printsignal    = false; 
    bool printplots     = true;

    // prefit and postfit yields 
    float pre_ZH      =   -999.;    float post_ZH      =   -999.;
    float pre_WH      =   -999.;    float post_WH      =   -999.;
    float pre_qqH     =   -999.;    float post_qqH     =   -999.;
    float pre_ggH     =   -999.;    float post_ggH     =   -999.;
    float pre_qqWW    =   -999.;    float post_qqWW    =   -999.;
    float pre_ggWW    =   -999.;    float post_ggWW    =   -999.;
    float pre_VV      =   -999.;    float post_VV      =   -999.;
    float pre_Top     =   -999.;    float post_Top     =   -999.;
    float pre_Zjets   =   -999.;    float post_Zjets   =   -999.;
    float pre_WjetsE  =   -999.;    float post_WjetsE  =   -999.;
    float pre_WjetsM  =   -999.;    float post_WjetsM  =   -999.;
    float pre_Wgamma  =   -999.;    float post_Wgamma  =   -999.;
    float pre_Wg3l    =   -999.;    float post_Wg3l    =   -999.;
    float pre_Ztt     =   -999.;    float post_Ztt     =   -999.;
    float pre_Ztt     =   -999.;    float post_Ztt     =   -999.;
    float pre_ZH_SM   =   -999.;    float post_ZH_SM   =   -999.;
    float pre_WH_SM   =   -999.;    float post_WH_SM   =   -999.;
    float pre_qqH_SM  =   -999.;    float post_qqH_SM  =   -999.;
    float pre_ggH_SM  =   -999.;    float post_ggH_SM  =   -999.;

    // -------------------------------------------------------------------------------
    //      PRE-FIT 
    // -------------------------------------------------------------------------------

    // get prefit normalization from a card   
    string line;
    ifstream incard(Form("%s/%i/hww%s_%ij_shape_8TeV.txt", carddir, mh, flavor, njet)); 
    
    char* tmp   =   "";

    if (incard.is_open())
    {
        while ( incard.good() )
        { 
            // get a line from input file
            getline (incard,line);

            // skip the line if it does not contain "rate" 
            if( line.find("rate")==string::npos ) continue;
              
            stringstream stream(line);
            stream >> tmp;
            stream >> pre_ZH;
            stream >> pre_WH;
            stream >> pre_qqH;
            stream >> pre_ggH;
            stream >> pre_qqWW;
            stream >> pre_ggWW;
            stream >> pre_VV;
            stream >> pre_Top;
            stream >> pre_Zjets;
            stream >> pre_WjetsE;
            stream >> pre_Wgamma;
            stream >> pre_Wg3l;
            stream >> pre_Ztt;
            stream >> pre_WjetsM;
            stream >> pre_ZH_SM;
            stream >> pre_WH_SM;
            stream >> pre_qqH_SM;
            stream >> pre_ggH_SM;

            if( !incard.good() ) continue;
        }

    }
    incard.close(); 

    // get pre-fit shapes 
    TFile* f = TFile::Open(Form("%s/%i/hww%s_%ij.input_8TeV.root", carddir, mh, flavor, njet));
    assert(f);
    TH1D *h1pre_ZH, *h1pre_WH, *h1pre_qqH, *h1pre_ggH, *h1pre_qqWW, *h1pre_ggWW, *h1pre_VV, *h1pre_Top, 
         *h1pre_Zjets, *h1pre_WjetsE, *h1pre_WjetsM, *h1pre_Wgamma, *h1pre_Wg3l, *h1pre_Ztt; 
    h1pre_ZH        = (TH1D*)f->Get("histo_ZH");
    h1pre_WH        = (TH1D*)f->Get("histo_WH");
    h1pre_qqH       = (TH1D*)f->Get("histo_qqH");
    h1pre_ggH       = (TH1D*)f->Get("histo_ggH");
    h1pre_qqWW      = (TH1D*)f->Get("histo_qqWW");
    h1pre_ggWW      = (TH1D*)f->Get("histo_ggWW");
    h1pre_VV        = (TH1D*)f->Get("histo_VV");
    h1pre_Top       = (TH1D*)f->Get("histo_Top");
    h1pre_Zjets     = (TH1D*)f->Get("histo_Zjets");
    h1pre_WjetsE    = (TH1D*)f->Get("histo_WjetsE");
    h1pre_WjetsM    = (TH1D*)f->Get("histo_WjetsM"); 
    h1pre_Wgamma    = (TH1D*)f->Get("histo_Wgamma");
    h1pre_Wg3l      = (TH1D*)f->Get("histo_Wg3l");
    h1pre_Ztt       = (TH1D*)f->Get("histo_Ztt");
    h1pre_Data      = (TH1D*)f->Get("histo_Data");

    // -------------------------------------------------------------------------------
    //      POST-FIT 
    // ------------------------------------------------------------------------------- 
    
    // LL_toyGraviton_fitGraviton_seed24280_fittedShape_floatMu.root  
    // LL_toySMHiggs_fitGraviton_seed24280_fittedShape_floatMu.root
    // LL_toyGraviton_fitSMHiggs_seed24280_fittedShape_floatMu.root   
    // LL_toySMHiggs_fitSMHiggs_seed24280_fittedShape_floatMu.root
    
    //TFile* postFile = TFile::Open(Form("%s/output/limits_%ij_shape_%s-fit-%i-all_fittedShape_floatMu.root", carddir, njet, flavor, mh));
    TFile* postFile = TFile::Open("spin2postfit/LL_toySMHiggs_fitSMHiggs_seed24280_fittedShape_floatMu.root");
    TH1D *h1post_ZH, *h1post_WH, *h1post_qqH, *h1post_ggH, *h1post_qqWW, *h1post_ggWW, *h1post_VV, *h1post_Top, 
         *h1post_Zjets, *h1post_WjetsE, *h1post_WjetsM, *h1post_Wgamma, *h1post_Wg3l, *h1post_Ztt; 
    
    RooRealVar *r;
    if(pre_ZH)       { h1post_ZH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_ZH", njet, flavor, energy)));      post_ZH      =   h1post_ZH->Integral();      } 
    if(pre_WH)       { h1post_WH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_WH", njet, flavor, energy)));      post_WH      =   h1post_WH->Integral();      } 
    if(pre_qqH)      { h1post_qqH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_qqH", njet, flavor, energy)));     post_qqH     =   h1post_qqH->Integral();     } 
    if(pre_ggH)      { h1post_ggH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_ggH", njet, flavor, energy)));     post_ggH     =   h1post_ggH->Integral();     } 
    if(pre_qqWW)     { h1post_qqWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_qqWW", njet, flavor, energy)));    post_qqWW    =   h1post_qqWW->Integral();    } 
    if(pre_ggWW)     { h1post_ggWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_ggWW", njet, flavor, energy)));    post_ggWW    =   h1post_ggWW->Integral();    } 
    if(pre_VV)       { h1post_VV     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_VV", njet, flavor, energy)));      post_VV      =   h1post_VV->Integral();      } 
    if(pre_Top)      { h1post_Top    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Top", njet, flavor, energy)));     post_Top     =   h1post_Top->Integral();     } 
    if(pre_Zjets)    { h1post_Zjets  = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Zjets", njet, flavor, energy)));   post_Zjets   =   h1post_Zjets->Integral();   } 
    if(pre_WjetsE)   { h1post_WjetsE = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_WjetsE", njet, flavor, energy)));  post_WjetsE  =   h1post_WjetsE->Integral();  } 
    if(pre_WjetsM)   { h1post_WjetsM = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_WjetsM", njet, flavor, energy)));  post_WjetsM  =   h1post_WjetsM->Integral();  } 
    if(pre_Wgamma)   { h1post_Wgamma = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Wgamma", njet, flavor, energy)));  post_Wgamma  =   h1post_Wgamma->Integral();  } 
    if(pre_Wg3l)     { h1post_Wg3l   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Wg3l", njet, flavor, energy)));    post_Wg3l    =   h1post_Wg3l->Integral();    } 
    if(pre_Ztt)      { h1post_Ztt    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Ztt", njet, flavor, energy)));     post_Ztt     =   h1post_Ztt->Integral();     } 
 
    // -------------------------------------------------------------------------------
    //      print table 
    // -------------------------------------------------------------------------------

    cout << "\\begin{table}[ht!]" << endl;
    cout << "\\begin{center}" << endl;
	cout << "\\begin{tabular}{c|cc|cc}" << endl;
    cout << "\\hline" << endl;   cout << "\\hline" << endl;   
    cout << "        Process &    N(prefit) &   N(postfit) & Difference(raw) &  Difference(\\%)  \\\\  " << endl;   
    cout << "\\hline" << endl;   cout << "\\hline" << endl;   
    if(printsignal) {
        printtableline("ZH",     pre_ZH,     post_ZH);
        printtableline("WH",     pre_WH,     post_WH);
        printtableline("qqH",    pre_qqH,    post_qqH);
        printtableline("ggH",    pre_ggH,    post_ggH);
        cout << "\\hline" << endl;   
    }
    printtableline("qqWW",   pre_qqWW,   post_qqWW);
    printtableline("ggWW",   pre_ggWW,   post_ggWW);
    cout << "\\hline" << endl;   
    printtableline("VV",     pre_VV,     post_VV);
    cout << "\\hline" << endl;   
    printtableline("Top",    pre_Top,    post_Top);
    cout << "\\hline" << endl;   
    printtableline("Zjets",  pre_Zjets,  post_Zjets);
    cout << "\\hline" << endl;   
    printtableline("Wjets($e$)", pre_WjetsE, post_WjetsE);
    printtableline("Wjets($\\mu$)", pre_WjetsM, post_WjetsM);
    cout << "\\hline" << endl;   
    printtableline("W$\\gamma$", pre_Wgamma, post_Wgamma);
    printtableline("W$\\gamma$*",   pre_Wg3l,   post_Wg3l);
    cout << "\\hline" << endl;   
    printtableline("Ztt",    pre_Ztt,    post_Ztt);
    cout << "\\hline" << endl;   cout << "\\hline" << endl;   
	cout << "\\end{tabular}" << endl;
	//cout << "\\caption{ }" << endl;
 	cout << "\\end{center}" << endl;
    cout << "\\end{table}" << endl;

    // -------------------------------------------------------------------------------
    //      draw histograms for each process
    // ------------------------------------------------------------------------------- 
    TString prefix =Form("%s/%i/hww%s_%ij_shape_8TeV.txt", carddir, mh, flavor, njet); 
    prefix.ReplaceAll("/", "_");
    prefix.ReplaceAll(".txt", "_");
    
    if(printplots) {
   /*
        drawoneprocess( Roll1DTo2D(h1pre_ZH),       Roll1DTo2D(h1post_ZH),      prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_WH),       Roll1DTo2D(h1post_WH),      prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_qqH),      Roll1DTo2D(h1post_qqH),     prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_ggH),      Roll1DTo2D(h1post_ggH),     prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_qqWW),     Roll1DTo2D(h1post_qqWW),    prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_ggWW),     Roll1DTo2D(h1post_ggWW),    prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_VV),       Roll1DTo2D(h1post_VV),      prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_Top),      Roll1DTo2D(h1post_Top),     prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_Zjets),    Roll1DTo2D(h1post_Zjets),   prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_WjetsE),   Roll1DTo2D(h1post_WjetsE),  prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_WjetsM),   Roll1DTo2D(h1post_WjetsM),  prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_Wgamma),   Roll1DTo2D(h1post_Wgamma),  prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_Wg3l),     Roll1DTo2D(h1post_Wg3l),    prefix,	doZoom); 
        drawoneprocess( Roll1DTo2D(h1pre_Ztt),      Roll1DTo2D(h1post_Ztt),     prefix,	doZoom); 
  */ 

        drawstack (
            h1post_ZH, h1post_WH, h1post_qqH, h1post_ggH, h1post_qqWW, h1post_ggWW, h1post_VV, h1post_Top, 
            h1post_Zjets, h1post_WjetsE, h1post_WjetsM, h1post_Wgamma, h1post_Wg3l, h1post_Ztt, 
            h1pre_Data, 
            true );

    } 

    if(printsignal) { 
        cout << "Data               : " << h1pre_Data->Integral() << endl;
        cout << "Bkg(prefit)        : " <<  pre_qqWW     + pre_ggWW      + pre_VV        + pre_Top
                                          + pre_Zjets    + pre_WjetsE    + pre_WjetsM 
                                          + pre_Wgamma   + pre_Wg3l      + pre_Ztt 
                                        << endl;
        cout << "Bkg(postfit)       : " <<  post_qqWW    + post_ggWW     + post_VV       + post_Top
                                          + post_Zjets   + post_WjetsE   + post_WjetsM 
                                          + post_Wgamma  + post_Wg3l     + post_Ztt 
                                        << endl;
        cout << "Signal(prefit)     : " << pre_ggH      + pre_qqH       + pre_ZH        + pre_WH    << endl;
        cout << "Signal(postfit)    : " << post_ggH     + post_qqH      + post_ZH       + post_WH   << endl;
    }
}
