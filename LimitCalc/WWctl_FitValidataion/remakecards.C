//
// Zoom
//
TH2F* Zoom(TH2F* h2) {

    float mtbins[7]    =   {60,70,80,90,100,110,120};
    float mllbins[6]   =   {12,30,45,60,75,100};

    TH2F* h2zoom = new TH2F(Form("%s_zoom", h2->GetName()), Form("%s_zoom", h2->GetTitle()), 6, mtbins, 5, mllbins);

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
    TH2F* h2 = new TH2F(h1->GetName(), h1->GetTitle(), 14, mtbins, 9, mllbins);
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
// print a line in a table 
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

void SetBinContentError(TH1F* &h1, int index, float content=0, float error=0){
    h1->SetBinContent(index, content);
    h1->SetBinError(index, error);
}

//
// main function
//
void remakecards(int TEST=1) {
   
    // options
    bool doCR1          = 0;
    bool doCR2          = 0;
    bool doCR1mt200     = 0;
    
    if(TEST==2) doCR1           = 1;
    if(TEST==3) doCR2           = 1; 
  
    if (doCR1+doCR2>1) {
        cout << "!!!! Aborted :: multiple tests selected !!!! "<< endl; 
        return;
    }

    //
    // get post-fit shape from LandS
    //
    int njet = 0;
    char* flavor = "of";
    char* energy = "8tev"; 

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
    float pre_data    =   -999.;    

    TH1F *h1pre_ZH, *h1pre_WH, *h1pre_qqH, *h1pre_ggH,*h1pre_qqWW, *h1pre_ggWW, *h1pre_VV, *h1pre_Top,
         *h1pre_Zjets, *h1pre_WjetsE, *h1pre_WjetsM, *h1pre_Wgamma, *h1pre_Wg3l, *h1pre_Ztt,
         *h1pre_data, *h1pre_qqWWUp, *h1pre_qqWWDown, *h1pre_qqWWNLOUp, *h1pre_qqWWNLODown,
         *h1pre_qqWWJESUp, *h1pre_qqWWJESDown, *h1pre_qqWWLepEffUp, *h1pre_qqWWLepEffDown,
         *h1pre_qqWWLepResUp, *h1pre_qqWWLepResDown, *h1pre_qqWWMETResUp, *h1pre_qqWWMETResDown,
         *h1pre_qqWWStatUp, *h1pre_qqWWStatDown, *h1pre_qqWWPDFUp, *h1pre_qqWWPDFDown;

    TH1F *h1post_ZH, *h1post_WH, *h1post_qqH, *h1post_ggH, *h1post_qqWW, *h1post_ggWW, *h1post_VV, *h1post_Top,
         *h1post_Zjets, *h1post_WjetsE, *h1post_WjetsM, *h1post_Wgamma, *h1post_Wg3l, *h1post_Ztt;

    TFile *preFile = TFile::Open("125/hwwof_0j.input_8TeV.root"); 
    h1pre_data          = (TH1F*)preFile->Get("histo_Data"); 
        h1pre_data->SetTitle("histo_new_Data"); 
        h1pre_data->SetName("histo_new_Data");
    h1pre_ZH          = (TH1F*)preFile->Get("histo_ZH"); 
        h1pre_ZH->SetTitle("histo_new_ZH");
        h1pre_ZH->SetName("histo_new_ZH");
        pre_ZH    =  h1pre_ZH->Integral(); 
    h1pre_WH          = (TH1F*)preFile->Get("histo_WH"); 
        h1pre_WH->SetTitle("histo_new_WH");
        h1pre_WH->SetName("histo_new_WH");
        pre_WH    =  h1pre_WH->Integral(); 
    h1pre_qqH          = (TH1F*)preFile->Get("histo_qqH"); 
        h1pre_qqH->SetTitle("histo_new_qqH");
        h1pre_qqH->SetName("histo_new_qqH");
        pre_qqH    =  h1pre_qqH->Integral(); 
    h1pre_ggH          = (TH1F*)preFile->Get("histo_ggH"); 
        h1pre_ggH->SetTitle("histo_new_ggH");
        h1pre_ggH->SetName("histo_new_ggH");
        pre_ggH    =  h1pre_ggH->Integral(); 
    h1pre_qqWW          = (TH1F*)preFile->Get("histo_qqWW"); 
        h1pre_qqWW->SetTitle("histo_new_qqWW");
        h1pre_qqWW->SetName("histo_new_qqWW");
        pre_qqWW    =  h1pre_qqWW->Integral(); 
    h1pre_ggWW          = (TH1F*)preFile->Get("histo_ggWW"); 
        h1pre_ggWW->SetTitle("histo_new_ggqWW");
        h1pre_ggWW->SetName("histo_new_ggWW");
        pre_ggWW    =  h1pre_ggWW->Integral(); 
    h1pre_VV          = (TH1F*)preFile->Get("histo_VV"); 
        h1pre_VV->SetTitle("histo_new_VV");
        h1pre_VV->SetName("histo_new_VV");
        pre_VV    =  h1pre_VV->Integral(); 
    h1pre_Top          = (TH1F*)preFile->Get("histo_Top"); 
        h1pre_Top->SetTitle("histo_new_Top");
        h1pre_Top->SetName("histo_new_Top");
        pre_Top    =  h1pre_Top->Integral(); 
    h1pre_Zjets          = (TH1F*)preFile->Get("histo_Zjets"); 
        h1pre_Zjets->SetTitle("histo_new_Zjets");
        h1pre_Zjets->SetName("histo_new_Zjets");
        pre_Zjets    =  h1pre_Zjets->Integral(); 
    h1pre_WjetsE          = (TH1F*)preFile->Get("histo_WjetsE"); 
        h1pre_WjetsE->SetTitle("histo_new_WjetsE");
        h1pre_WjetsE->SetName("histo_new_WjetsE");
        pre_WjetsE    =  h1pre_WjetsE->Integral(); 
    h1pre_WjetsM          = (TH1F*)preFile->Get("histo_WjetsM"); 
        h1pre_WjetsM->SetTitle("histo_new_WjetsM");
        h1pre_WjetsM->SetName("histo_new_WjetsM");
        pre_WjetsM    =  h1pre_WjetsM->Integral(); 
    h1pre_Wgamma          = (TH1F*)preFile->Get("histo_Wgamma"); 
        h1pre_Wgamma->SetTitle("histo_new_Wgamma");
        h1pre_Wgamma->SetName("histo_new_Wgamma");
        pre_Wgamma    =  h1pre_Wgamma->Integral(); 
    h1pre_Wg3l          = (TH1F*)preFile->Get("histo_Wg3l"); 
        h1pre_Wg3l->SetTitle("histo_new_Wg3l");
        h1pre_Wg3l->SetName("histo_new_Wg3l");
        pre_Wg3l    =  h1pre_Wg3l->Integral(); 
    h1pre_Ztt          = (TH1F*)preFile->Get("histo_Ztt"); 
        h1pre_Ztt->SetTitle("histo_new_Ztt");
        h1pre_Ztt->SetName("histo_new_Ztt");
        pre_Ztt    =  h1pre_Ztt->Integral(); 

    h1pre_qqWWUp        = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAWWBoundingUp"); 
        h1pre_qqWWUp->SetTitle("histo_new_qqWW_CMS_hww_MVAWWBoundingUp");
        h1pre_qqWWUp->SetName("histo_new_qqWW_CMS_hww_MVAWWBoundingUp");
    h1pre_qqWWDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAWWBoundingDown"); 
        h1pre_qqWWDown->SetTitle("histo_new_qqWW_CMS_hww_MVAWWBoundingDown");
        h1pre_qqWWDown->SetName("histo_new_qqWW_CMS_hww_MVAWWBoundingDown");
    h1pre_qqWWNLOUp     = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAWWNLOBoundingUp"); 
        h1pre_qqWWNLOUp->SetTitle("histo_new_qqWW_CMS_hww_MVAWWNLOBoundingUp");
        h1pre_qqWWNLOUp->SetName("histo_new_qqWW_CMS_hww_MVAWWNLOBoundingUp");
    h1pre_qqWWNLODown   = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAWWNLOBoundingDown"); 
        h1pre_qqWWNLODown->SetTitle("histo_new_qqWW_CMS_hww_MVAWWNLOBoundingDown");
        h1pre_qqWWNLODown->SetName("histo_new_qqWW_CMS_hww_MVAWWNLOBoundingDown");
    h1pre_qqWWJESUp        = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAJESBoundingUp"); 
        h1pre_qqWWJESUp->SetTitle("histo_new_qqWW_CMS_hww_MVAJESBoundingUp");
        h1pre_qqWWJESUp->SetName("histo_new_qqWW_CMS_hww_MVAJESBoundingUp");
    h1pre_qqWWJESDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAJESBoundingDown"); 
        h1pre_qqWWJESDown->SetTitle("histo_new_qqWW_CMS_hww_MVAJESBoundingDown");
        h1pre_qqWWJESDown->SetName("histo_new_qqWW_CMS_hww_MVAJESBoundingDown");
    h1pre_qqWWLepEffUp        = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVALepEffBoundingUp"); 
        h1pre_qqWWLepEffUp->SetTitle("histo_new_qqWW_CMS_hww_MVALepEffBoundingUp");
        h1pre_qqWWLepEffUp->SetName("histo_new_qqWW_CMS_hww_MVALepEffBoundingUp");
    h1pre_qqWWLepEffDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVALepEffBoundingDown"); 
        h1pre_qqWWLepEffDown->SetTitle("histo_new_qqWW_CMS_hww_MVALepEffBoundingDown");
        h1pre_qqWWLepEffDown->SetName("histo_new_qqWW_CMS_hww_MVALepEffBoundingDown");
    h1pre_qqWWLepResUp        = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVALepResBoundingUp"); 
        h1pre_qqWWLepResUp->SetTitle("histo_new_qqWW_CMS_hww_MVALepResBoundingUp");
        h1pre_qqWWLepResUp->SetName("histo_new_qqWW_CMS_hww_MVALepResBoundingUp");
    h1pre_qqWWLepResDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVALepResBoundingDown"); 
        h1pre_qqWWLepResDown->SetTitle("histo_new_qqWW_CMS_hww_MVALepResBoundingDown");
        h1pre_qqWWLepResDown->SetName("histo_new_qqWW_CMS_hww_MVALepResBoundingDown");
    h1pre_qqWWMETResUp        = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAMETResBoundingUp"); 
        h1pre_qqWWMETResUp->SetTitle("histo_new_qqWW_CMS_hww_MVAMETResBoundingUp");
        h1pre_qqWWMETResUp->SetName("histo_new_qqWW_CMS_hww_MVAMETResBoundingUp");
    h1pre_qqWWMETResDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_MVAMETResBoundingDown"); 
        h1pre_qqWWMETResDown->SetTitle("histo_new_qqWW_CMS_hww_MVAMETResBoundingDown");
        h1pre_qqWWMETResDown->SetName("histo_new_qqWW_CMS_hww_MVAMETResBoundingDown");
    h1pre_qqWWStatUp      = (TH1F*)preFile->Get("histo_qqWW_CMS_hwwof_0j_MVAqqWWStatBounding_8TeVUp"); 
        h1pre_qqWWStatUp->SetTitle("histo_new_qqWW_CMS_hwwof_0j_MVAqqWWStatBounding_8TeVUp");
        h1pre_qqWWStatUp->SetName("histo_new_qqWW_CMS_hwwof_0j_MVAqqWWStatBounding_8TeVUp");
    h1pre_qqWWStatDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hwwof_0j_MVAqqWWStatBounding_8TeVDown"); 
        h1pre_qqWWStatDown->SetTitle("histo_new_qqWW_CMS_hwwof_0j_MVAqqWWStatBounding_8TeVDown");
        h1pre_qqWWStatDown->SetName("histo_new_qqWW_CMS_hwwof_0j_MVAqqWWStatBounding_8TeVDown");
    h1pre_qqWWPDFUp      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_PDFqqWWUp"); 
        h1pre_qqWWPDFUp->SetTitle("histo_new_qqWW_CMS_hww_PDFqqWWUp");
        h1pre_qqWWPDFUp->SetName("histo_new_qqWW_CMS_hww_PDFqqWWUp");
    h1pre_qqWWPDFDown      = (TH1F*)preFile->Get("histo_qqWW_CMS_hww_PDFqqWWDown"); 
        h1pre_qqWWPDFDown->SetTitle("histo_new_qqWW_CMS_hww_PDFqqWWDown");
        h1pre_qqWWPDFDown->SetName("histo_new_qqWW_CMS_hww_PDFqqWWDown");

    TFile* postFile = TFile::Open("nominal_fittedShape_floatMu.root");
    h1post_ZH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_ZH", njet, flavor, energy)));     
        post_ZH      =   h1post_ZH->Integral();      
    h1post_WH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_WH", njet, flavor, energy)));     
        post_WH      =   h1post_WH->Integral();      
    h1post_qqH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_qqH", njet, flavor, energy)));    
        post_qqH     =   h1post_qqH->Integral();     
    h1post_ggH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_ggH", njet, flavor, energy)));    
        post_ggH     =   h1post_ggH->Integral();     
    h1post_qqWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_qqWW", njet, flavor, energy)));   
        post_qqWW    =   h1post_qqWW->Integral();    
    h1post_ggWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_ggWW", njet, flavor, energy)));   
        post_ggWW    =   h1post_ggWW->Integral();    
    h1post_VV     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_VV", njet, flavor, energy)));     
        post_VV      =   h1post_VV->Integral();      
    h1post_Top    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Top", njet, flavor, energy)));    
        post_Top     =   h1post_Top->Integral();     
    h1post_Zjets  = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Zjets", njet, flavor, energy)));    
        post_Zjets   =   h1post_Zjets->Integral();   
    h1post_WjetsE = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_WjetsE", njet, flavor, energy))); 
        post_WjetsE  =   h1post_WjetsE->Integral();  
    h1post_WjetsM = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_WjetsM", njet, flavor, energy))); 
        post_WjetsM  =   h1post_WjetsM->Integral();  
    h1post_Wgamma = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Wgamma", njet, flavor, energy))); 
        post_Wgamma  =   h1post_Wgamma->Integral();  
    h1post_Wg3l   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Wg3l", njet, flavor, energy)));   
        post_Wg3l    =   h1post_Wg3l->Integral();    
    h1post_Ztt    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%i%s%s_Ztt", njet, flavor, energy)));    
        post_Ztt     =   h1post_Ztt->Integral();     

        if(0) {  // debug
        //cout << "\\begin{table}[ht!]" << endl;
        //cout << "\\begin{center}" << endl;
        //cout << "\\begin{tabular}{c|cc|cc}" << endl;
        //cout << "\\hline" << endl;   cout << "\\hline" << endl;
        cout << "        Process &    N(prefit) &   N(postfit) & Difference(raw) &  Difference(\\%)  \\\\  " << endl;
        //cout << "\\hline" << endl;   cout << "\\hline" << endl;
        printtableline("ZH",     pre_ZH,     post_ZH);
        printtableline("WH",     pre_WH,     post_WH);
        printtableline("qqH",    pre_qqH,    post_qqH);
        printtableline("ggH",    pre_ggH,    post_ggH);
        //cout << "\\hline" << endl;
        printtableline("qqWW",   pre_qqWW,   post_qqWW);
        printtableline("ggWW",   pre_ggWW,   post_ggWW);
        //cout << "\\hline" << endl;
        printtableline("VV",     pre_VV,     post_VV);
        //cout << "\\hline" << endl;
        printtableline("Top",    pre_Top,    post_Top);
        //cout << "\\hline" << endl;
        printtableline("Zjets",  pre_Zjets,  post_Zjets);
        //cout << "\\hline" << endl;
        printtableline("Wjets($e$)", pre_WjetsE, post_WjetsE);
        printtableline("Wjets($\\mu$)", pre_WjetsM, post_WjetsM);
        //cout << "\\hline" << endl;
        printtableline("W$\\gamma$", pre_Wgamma, post_Wgamma);
        printtableline("W$\\gamma$*",   pre_Wg3l,   post_Wg3l);
        //cout << "\\hline" << endl;
        printtableline("Ztt",    pre_Ztt,    post_Ztt);
        //cout << "\\hline" << endl;   cout << "\\hline" << endl;
        //cout << "\\end{tabular}" << endl;
        //cout << "\\caption{ }" << endl;
        //cout << "\\end{center}" << endl;
        //cout << "\\end{table}" << endl;

        }


    //
    // chop off the unnecessary region
    //
    bool removethisbin = false; 
    for(int i=1; i<127; i++){
        
        // remove CR2 : mll>60 && mT<120
        if ( doCR1 &&  
            ( (i>=4 && i<=9)    || (i>=13 && i<=18) || (i>=22 && i<=27) 
           || (i>=31 && i<=36)  || (i>=40 && i<=45) || (i>=49 && i<=54) )
           )  removethisbin = true; 
        
        // remove CR1 : mT>120
        if ( doCR2 &&  
             i>=55
           )  removethisbin = true; 

        // remove CR1     : mll>60 && mT<120 
        // remove mT>200  
        if ( doCR1mt200 &&  
            ( (i>=4 && i<=9)    || (i>=13 && i<=18) || (i>=22 && i<=27) 
           || (i>=31 && i<=36)  || (i>=40 && i<=45) || (i>=49 && i<=54) 
           || (i>=91) )
           )  removethisbin = true; 

        if(removethisbin) {
            SetBinContentError(h1post_ZH, i);
            SetBinContentError(h1post_WH, i);
            SetBinContentError(h1post_qqH, i);
            SetBinContentError(h1post_ggH, i);
            SetBinContentError(h1post_ggWW, i);
            SetBinContentError(h1post_VV, i);
            SetBinContentError(h1post_Top, i);
            SetBinContentError(h1post_Zjets, i);
            SetBinContentError(h1post_WjetsE, i);
            SetBinContentError(h1post_WjetsM, i);
            SetBinContentError(h1post_Wgamma, i);
            SetBinContentError(h1post_Wg3l, i);
            SetBinContentError(h1post_Ztt, i);

            SetBinContentError(h1pre_data, i);
            SetBinContentError(h1pre_ZH, i);
            SetBinContentError(h1pre_WH, i);
            SetBinContentError(h1pre_qqH, i);
            SetBinContentError(h1pre_ggH, i);
            SetBinContentError(h1pre_qqWW, i);
            SetBinContentError(h1pre_ggWW, i);
            SetBinContentError(h1pre_VV, i);
            SetBinContentError(h1pre_Top, i);
            SetBinContentError(h1pre_Zjets, i);
            SetBinContentError(h1pre_WjetsE, i);
            SetBinContentError(h1pre_WjetsM, i);
            SetBinContentError(h1pre_Wgamma, i);
            SetBinContentError(h1pre_Wg3l, i);
            SetBinContentError(h1pre_Ztt, i);

            SetBinContentError(h1pre_qqWWUp, i);
            SetBinContentError(h1pre_qqWWDown, i);
            SetBinContentError(h1pre_qqWWNLOUp, i);
            SetBinContentError(h1pre_qqWWNLODown, i);
            SetBinContentError(h1pre_qqWWJESUp, i);
            SetBinContentError(h1pre_qqWWJESDown, i);
            SetBinContentError(h1pre_qqWWLepEffUp, i);
            SetBinContentError(h1pre_qqWWLepEffDown, i);
            SetBinContentError(h1pre_qqWWLepResUp, i);
            SetBinContentError(h1pre_qqWWLepResDown, i);
            SetBinContentError(h1pre_qqWWMETResUp, i);
            SetBinContentError(h1pre_qqWWMETResDown, i);
            SetBinContentError(h1pre_qqWWStatUp, i);
            SetBinContentError(h1pre_qqWWStatDown, i);
            SetBinContentError(h1pre_qqWWPDFUp, i);
            SetBinContentError(h1pre_qqWWPDFDown, i);
        }

        removethisbin = false;
    }

    //
    pre_data    =  h1pre_data->Integral(); 
    pre_ZH    =  h1pre_ZH->Integral(); 
    pre_WH    =  h1pre_WH->Integral(); 
    pre_qqH    =  h1pre_qqH->Integral(); 
    pre_ggH    =  h1pre_ggH->Integral(); 
    pre_qqWW    =  h1pre_qqWW->Integral(); 
    pre_ggWW    =  h1pre_ggWW->Integral(); 
    pre_VV    =  h1pre_VV->Integral(); 
    pre_Top    =  h1pre_Top->Integral(); 
    pre_Zjets    =  h1pre_Zjets->Integral(); 
    pre_WjetsE    =  h1pre_WjetsE->Integral(); 
    pre_WjetsM    =  h1pre_WjetsM->Integral(); 
    pre_Wgamma    =  h1pre_Wgamma->Integral(); 
    pre_Wg3l    =  h1pre_Wg3l->Integral(); 
    pre_Ztt    =  h1pre_Ztt->Integral(); 

//    h1pre_qqWW->Draw("histo"); 

    TH2F *h2tmp = Roll1DTo2D(h1pre_data);
    h2tmp->Draw("colz");

    //
    // replace the central shapes
    //
    char* test = "newdefault";
    if(doCR1)       test = "CR1";
    if(doCR2)       test = "CR2";
    if(doCR1mt200)  test = "CR1mt200";
    gSystem->Exec(Form("cp 125/hwwof_0j.input_8TeV.root 125/hwwof_0j.input_8TeV_%s.root", test));
    TFile *inputrootfile = new TFile(Form("125/hwwof_0j.input_8TeV_%s.root", test), "UPDATE");
    gROOT->cd();
    inputrootfile->cd();
  
    h1post_ZH->SetDirectory(0); h1post_ZH->Write();
    h1post_WH->SetDirectory(0); h1post_WH->Write();
    h1post_qqH->SetDirectory(0); h1post_qqH->Write();
    h1post_ggH->SetDirectory(0); h1post_ggH->Write();
    h1post_ggWW->SetDirectory(0); h1post_ggWW->Write();
    h1post_VV->SetDirectory(0); h1post_VV->Write();
    h1post_Top->SetDirectory(0); h1post_Top->Write();
    h1post_Zjets->SetDirectory(0); h1post_Zjets->Write();
    h1post_WjetsE->SetDirectory(0); h1post_WjetsE->Write();
    h1post_WjetsM->SetDirectory(0); h1post_WjetsM->Write();
    h1post_Wgamma->SetDirectory(0); h1post_Wgamma->Write();
    h1post_Wg3l->SetDirectory(0); h1post_Wg3l->Write();
    h1post_Ztt->SetDirectory(0); h1post_Ztt->Write();
    h1pre_data->SetDirectory(0); h1pre_data->Write();
    h1pre_qqWW->SetDirectory(0); h1pre_qqWW->Write();
    //h1post_qqWW->SetDirectory(0); h1post_qqWW->Write(); // save postfit qqWW : fitted with/without the three nuisances
    h1pre_qqWWUp->SetDirectory(0); h1pre_qqWWUp->Write();
    h1pre_qqWWDown->SetDirectory(0); h1pre_qqWWDown->Write();
    h1pre_qqWWNLOUp->SetDirectory(0); h1pre_qqWWNLOUp->Write();
    h1pre_qqWWNLODown->SetDirectory(0); h1pre_qqWWNLODown->Write();
    h1pre_qqWWJESUp->SetDirectory(0); h1pre_qqWWJESUp->Write();
    h1pre_qqWWJESDown->SetDirectory(0); h1pre_qqWWJESDown->Write();
    h1pre_qqWWLepEffUp->SetDirectory(0); h1pre_qqWWLepEffUp->Write();
    h1pre_qqWWLepEffDown->SetDirectory(0); h1pre_qqWWLepEffDown->Write();
    h1pre_qqWWLepResUp->SetDirectory(0); h1pre_qqWWLepResUp->Write();
    h1pre_qqWWLepResDown->SetDirectory(0); h1pre_qqWWLepResDown->Write();
    h1pre_qqWWMETResUp->SetDirectory(0); h1pre_qqWWMETResUp->Write();
    h1pre_qqWWMETResDown->SetDirectory(0); h1pre_qqWWMETResDown->Write();
    h1pre_qqWWStatUp->SetDirectory(0); h1pre_qqWWStatUp->Write();
    h1pre_qqWWStatDown->SetDirectory(0); h1pre_qqWWStatDown->Write();
    h1pre_qqWWPDFUp->SetDirectory(0); h1pre_qqWWPDFUp->Write();
    h1pre_qqWWPDFDown->SetDirectory(0); h1pre_qqWWPDFDown->Write();

    inputrootfile->Close();

    //
    // print out a new card
    // 
    ofstream fout;
    fout.open(Form("125/hwwof_0j_shape_8TeV_%s.txt", test));

    fout << "imax 1 number of channels" << endl; 
    fout << "jmax * number of background" << endl;
    fout << "kmax * number of nuisance parameters" << endl;
    fout << "Observation " << h1pre_data->Integral() << endl; 
    fout << Form("shapes *   *   hwwof_0j.input_8TeV_%s.root  histo_new_$PROCESS histo_new_$PROCESS_$SYSTEMATIC", test) << endl;
    fout << Form("shapes data_obs * hwwof_0j.input_8TeV_%s.root  histo_new_Data", test) << endl;
    fout << "bin                               j0of  j0of  j0of  j0of  j0of  j0of  j0of  j0of  j0of  j0of  j0of  j0of j0of j0of" << endl;
    fout << "process                            ZH    WH   qqH   ggH   qqWW   ggWW  VV    Top  Zjets WjetsE Wgamma Wg3l  Ztt  WjetsM" << endl;
    fout << "process                            -3    -2    -1    0      1     2     3     4     5     6     7    8      9    10" << endl;
    fout << "rate                             " 
         << Form("%.3f ",h1post_ZH->Integral()) 
         << Form("%.3f ",h1post_WH->Integral()) 
         << Form("%.3f ",h1post_qqH->Integral()) 
         << Form("%.3f ",h1post_ggH->Integral()) 
         << Form("%.3f ",h1pre_qqWW->Integral()) 
         << Form("%.3f ",h1post_ggWW->Integral()) 
         << Form("%.3f ",h1post_VV->Integral()) 
         << Form("%.3f ",h1post_Top->Integral()) 
         << Form("%.3f ",h1post_Zjets->Integral()) 
         << Form("%.3f ",h1post_WjetsE->Integral()) 
         << Form("%.3f ",h1post_Wgamma->Integral()) 
         << Form("%.3f ",h1post_Wg3l->Integral()) 
         << Form("%.3f ",h1post_Ztt->Integral()) 
         << Form("%.3f ",h1post_WjetsM->Integral()) 
         << endl;
    fout << "CMS_hww_0j_WW_8TeV_SHAPE               lnU     -     -     -     -   2.000   -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_MVAWWBounding                  shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_MVAWWNLOBounding               shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_MVAJESBounding                 shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_MVALepEffBounding              shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_MVALepResBounding              shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_MVAMETResBounding              shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hww_PDFqqWW                        shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    fout << "CMS_hwwof_0j_MVAqqWWStatBounding_8TeV  shape   -     -     -     -    1.0    -     -     -     -     -     -     -     -    -" << endl;
    
    fout.close();
    
    // ZH     WH      qqH     ggH     qqWW        ggWW    VV      Top     Zjets   WjetsE  Wgamma  Wg3l    Ztt     WjetsM  data  
    // 1.689  5.701   2.930   228.132 3969.634    210.627 132.609 498.701 0.000   282.793 115.578 167.802 45.996  331.797 5729 
   
    if(0) { //debug

        cout << Form("%.3f", 1.689 + 5.701) << "\t" 
            << Form("%.3f", 2.930) << "\t" 
            << Form("%.3f", 228.132) << "\t" 
            << Form("%.3f", 3969.634) << "\t" 
            << Form("%.3f", 210.627) << "\t" 
            << Form("%.3f", 132.609) << "\t" 
            << Form("%.3f", 498.701) << "\t" 
            << Form("%.3f", 282.793) << "\t" 
            << Form("%.3f", 331.797) << "\t" 
            << Form("%.3f", 115.578) << "\t" 
            << Form("%.3f", 167.802) << "\t" 
            << Form("%.3f", 45.996) << "\t"
            << 5729 << endl; 

        cout << Form("%.3f", 1.689 + 5.701 - pre_ZH - pre_WH) << "\t" 
            << Form("%.3f", 2.930 - pre_qqH) << "\t" 
            << Form("%.3f", 228.132 - pre_ggH) << "\t" 
            << Form("%.3f", 3969.634 - pre_qqWW) << "\t" 
            << Form("%.3f", 210.627 - pre_ggWW) << "\t" 
            << Form("%.3f", 132.609 - pre_VV) << "\t" 
            << Form("%.3f", 498.701 - pre_Top) << "\t" 
            << Form("%.3f", 282.793 - pre_WjetsE) << "\t" 
            << Form("%.3f", 331.797 - pre_WjetsM) << "\t" 
            << Form("%.3f", 115.578 - pre_Wgamma) << "\t" 
            << Form("%.3f", 167.802 - pre_Wg3l) << "\t" 
            << Form("%.3f", 45.996 - pre_Ztt) << "\t"
            << 5729 - pre_data << endl; 
    }
}
