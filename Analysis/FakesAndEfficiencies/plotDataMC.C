#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBenchmark.h>             // class to track macro running statistics
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>      // graph class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O

#include "CPlot.hh"	     // helper class for plots
#include "MitStyleRemix.hh"  // style settings for drawing
#endif

void plotDataMC(const TString outdir   = "Data/extra",
        const TString mcfname  = "MC/eff.root",
        const TString datfname = "Data/eff.root"
        ) {
    gBenchmark->Start("plotDataMC");

    //--------------------------------------------------------------------------------------------------------------
    // Settings 
    //============================================================================================================== 

    vector<TString> etalabelv;
    vector<TString> ptlabelv;

    CPlot::sOutDir = outdir;
    TString format = "png";

    // y-axis ranges
    Double_t efflow   = 0.0, effhigh   = 1.20;
    Double_t scalelow = 0.5, scalehigh = 1.5;


    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code 
    //==============================================================================================================   

    TFile *mcfile=0;      
    TFile *datafile=0;
    if (mcfname != "null")     mcfile = new TFile(mcfname);
    if (datfname != "null")    datafile = new TFile(datfname);

    gSystem->mkdir(outdir,true);
    TFile outfile(outdir+"/smurf.root", "RECREATE");

    TH2F *hMCEff=0,   *hMCErrl=0,   *hMCErrh=0;
    TH2F *hDataEff=0, *hDataErrl=0, *hDataErrh=0;

    TGraphAsymmErrors *grMCEffEta=0, *grDataEffEta=0, *grScaleEta=0;
    TGraphAsymmErrors *grMCEffPt=0,  *grDataEffPt=0,  *grScalePt=0;

    vector<TGraphAsymmErrors*> mceff_vs_pt_per_etav;
    vector<TGraphAsymmErrors*> eff_vs_pt_per_etav;
    vector<TGraphAsymmErrors*> mceff_vs_eta_per_ptv;
    vector<TGraphAsymmErrors*> eff_vs_eta_per_ptv;
    vector<TGraphAsymmErrors*> scale_vs_pt_per_etav;
    vector<TGraphAsymmErrors*> scale_vs_eta_per_ptv;

    if (mcfile) {
        grMCEffEta = (TGraphAsymmErrors*)mcfile->Get("grEffEta");
        grMCEffPt  = (TGraphAsymmErrors*)mcfile->Get("grEffPt");
    }

    if(grMCEffPt) {
        if (grMCEffPt->GetX()[grMCEffPt->GetN()-1] > 1000) {
            grMCEffPt->SetPoint(grMCEffPt->GetN()-1, 75, grMCEffPt->GetY()[grMCEffPt->GetN()-1]);
            grMCEffPt->SetPointError(grMCEffPt->GetN()-1,25,25, 
                    grMCEffPt->GetErrorYlow(grMCEffPt->GetN()-1),
                    grMCEffPt->GetErrorYhigh(grMCEffPt->GetN()-1));
        }
    }  			      

    if (mcfile) {         
        hMCEff  = (TH2F*)mcfile->Get("hEffEtaPt");
        hMCErrl = (TH2F*)mcfile->Get("hErrlEtaPt");
        hMCErrh = (TH2F*)mcfile->Get("hErrhEtaPt");  
    }

    if (datafile) {
        grDataEffEta = (TGraphAsymmErrors*)datafile->Get("grEffEta");
        grDataEffPt  = (TGraphAsymmErrors*)datafile->Get("grEffPt");
    }
    if(grDataEffPt) {
        if (grDataEffPt->GetX()[grDataEffPt->GetN()-1] > 1000) {
            grDataEffPt->SetPoint(grDataEffPt->GetN()-1, 75, grDataEffPt->GetY()[grDataEffPt->GetN()-1]);
            grDataEffPt->SetPointError(grDataEffPt->GetN()-1,25,25, 
                grDataEffPt->GetErrorYlow(grDataEffPt->GetN()-1),
                grDataEffPt->GetErrorYhigh(grDataEffPt->GetN()-1));
        }
    }

    if (datafile) {
        hDataEff  = (TH2F*)datafile->Get("hEffEtaPt");
        hDataErrl = (TH2F*)datafile->Get("hErrlEtaPt");
        hDataErrh = (TH2F*)datafile->Get("hErrhEtaPt");
    }

std::cout << "mark 1" << std::endl;

    if(grMCEffEta && grDataEffEta) {
        grScaleEta = new TGraphAsymmErrors(grMCEffEta->GetN());
        for(Int_t i=0; i<grMCEffEta->GetN(); i++) {
            Double_t mcval   = grMCEffEta->GetY()[i];
            Double_t dataval = grDataEffEta->GetY()[i];

            if (mcval > 0) {
                Double_t scale   = dataval/mcval;
                grScaleEta->SetPoint(i,grMCEffEta->GetX()[i],scale);

                Double_t mcerrl   = grMCEffEta->GetErrorYlow(i);
                Double_t mcerrh   = grMCEffEta->GetErrorYhigh(i);
                Double_t dataerrl = grDataEffEta->GetErrorYlow(i);
                Double_t dataerrh = grDataEffEta->GetErrorYhigh(i);
                grScaleEta->SetPointError(i, grMCEffEta->GetErrorXlow(i), grMCEffEta->GetErrorXlow(i),
                        scale*sqrt(mcerrl*mcerrl/mcval/mcval + dataerrl*dataerrl/dataval/dataval),
                        scale*sqrt(mcerrh*mcerrh/mcval/mcval + dataerrh*dataerrh/dataval/dataval));
            } else {
                grScaleEta->SetPoint(i,grMCEffEta->GetX()[i],0);
                grScaleEta->SetPointError(i, 0, 0, 0, 0);
            }
        }
    }

    if(grMCEffPt && grDataEffPt) {
        grScalePt = new TGraphAsymmErrors(grMCEffPt->GetN());
        for(Int_t i=0; i<grMCEffPt->GetN(); i++) {
            Double_t mcval   = grMCEffPt->GetY()[i];
            Double_t dataval = grDataEffPt->GetY()[i];

            if (mcval > 0) {
                Double_t scale   = dataval/mcval;
                grScalePt->SetPoint(i, grMCEffPt->GetX()[i],scale);
                if(i==grMCEffPt->GetN()-1)
                    grScalePt->SetPoint(i,165,scale);

                Double_t mcerrl   = grMCEffPt->GetErrorYlow(i);
                Double_t mcerrh   = grMCEffPt->GetErrorYhigh(i);
                Double_t dataerrl = grDataEffPt->GetErrorYlow(i);
                Double_t dataerrh = grDataEffPt->GetErrorYhigh(i);

                grScalePt->SetPointError(i, grMCEffPt->GetErrorXlow(i), grMCEffPt->GetErrorXlow(i),
                        scale*sqrt(mcerrl*mcerrl/mcval/mcval + dataerrl*dataerrl/dataval/dataval),
                        scale*sqrt(mcerrh*mcerrh/mcval/mcval + dataerrh*dataerrh/dataval/dataval));
            } else {
                grScalePt->SetPoint(i,grMCEffEta->GetX()[i],0);
                grScalePt->SetPointError(i, 0, 0, 0, 0);
            }
        }
    }

std::cout << "mark 2" << std::endl;
std::cout << hDataEff->GetEntries() << std::endl;
    if(!hMCEff && hDataEff && hDataEff->GetEntries()>0) {

        const Int_t nx = hDataEff->GetNbinsX();
        const Int_t ny = hDataEff->GetNbinsY();

        outfile.cd();

        Double_t *arrX = (Double_t*)hDataEff->GetXaxis()->GetXbins()->GetArray();
        Double_t *arrY = (Double_t*)hDataEff->GetYaxis()->GetXbins()->GetArray();
        //arrY[ny] = 100.0;

        // data
        TH2D *h2_smurf_DataEff = new TH2D("smurf_DataEff", "DataEff", ny, arrY, nx, arrX);
        TH2D *h2_smurf_DataErrLow = new TH2D("smurf_DataErrLow", "DataErrLow", ny, arrY, nx, arrX);
        TH2D *h2_smurf_DataErrHigh = new TH2D("smurf_DataErrHigh", "DataErrHigh", ny, arrY, nx, arrX);

        for(Int_t ix=1; ix<=nx; ix++) {

            etalabelv.push_back(Form("%4.3f < |#eta| < %4.3f", hDataEff->GetXaxis()->GetBinLowEdge(ix), hDataEff->GetXaxis()->GetBinLowEdge(ix+1)));

            Double_t xval[ny], xerr[ny];
            Double_t effval[ny],   efferrl[ny],   efferrh[ny];

            for(Int_t iy=1; iy<=ny; iy++) {
                xval[iy-1] = 0.5*(hDataEff->GetYaxis()->GetBinLowEdge(iy+1) + hDataEff->GetYaxis()->GetBinLowEdge(iy));
                xerr[iy-1] = 0.5*(hDataEff->GetYaxis()->GetBinLowEdge(iy+1) - hDataEff->GetYaxis()->GetBinLowEdge(iy));
                if(iy==ny && xval[iy-1] > 1000.0) {
                    xval[iy-1] = 75;
                    xerr[iy-1] = 25;
                }

                Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
                Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
                Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
                effval[iy-1]  = dataeff;
                efferrl[iy-1] = dataerrl;
                efferrh[iy-1] = dataerrh;

                // fill smurf hists
                h2_smurf_DataEff      ->SetBinContent(iy, ix, effval[iy-1]);
                h2_smurf_DataErrLow   ->SetBinContent(iy, ix, efferrl[iy-1]);
                h2_smurf_DataErrHigh  ->SetBinContent(iy, ix, efferrh[iy-1]);

            }

            std::cout << "pushing back" << std::endl;
            eff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,effval,xerr,xerr,efferrl,efferrh));
        }
    }

std::cout << "mark 3" << std::endl;

    if(hMCEff && hDataEff && hMCEff->GetEntries()>0 && hDataEff->GetEntries()>0) {
        const Int_t nx = hMCEff->GetNbinsX();
        const Int_t ny = hMCEff->GetNbinsY(); 

        outfile.cd();

        Double_t *arrX = (Double_t*)hMCEff->GetXaxis()->GetXbins()->GetArray();
        Double_t *arrY = (Double_t*)hMCEff->GetYaxis()->GetXbins()->GetArray();
        arrY[ny] = 100.0;

        // mc
        TH2D *h2_smurf_MCEff = new TH2D("smurf_MCEff", "MCEff", ny, arrY, nx, arrX);
        TH2D *h2_smurf_MCErrLow = new TH2D("smurf_MCErrLow", "MCErrLow", ny, arrY, nx, arrX);
        TH2D *h2_smurf_MCErrHigh = new TH2D("smurf_MCErrHigh", "MCErrHigh", ny, arrY, nx, arrX);

        // data
        TH2D *h2_smurf_DataEff = new TH2D("smurf_DataEff", "DataEff", ny, arrY, nx, arrX);
        TH2D *h2_smurf_DataErrLow = new TH2D("smurf_DataErrLow", "DataErrLow", ny, arrY, nx, arrX);
        TH2D *h2_smurf_DataErrHigh = new TH2D("smurf_DataErrHigh", "DataErrHigh", ny, arrY, nx, arrX);

        // sf
        TH2D *h2_smurf_SF = new TH2D("smurf_SF", "SF", ny, arrY, nx, arrX);
        TH2D *h2_smurf_SFErrLow = new TH2D("smurf_SFErrLow", "SFErrLow", ny, arrY, nx, arrX);
        TH2D *h2_smurf_SFErrHigh = new TH2D("smurf_SFErrHigh", "SFErrHigh", ny, arrY, nx, arrX);

        for(Int_t ix=1; ix<=nx; ix++) {

            etalabelv.push_back(Form("%4.3f < |#eta| < %4.3f", hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1)));

            Double_t xval[ny], xerr[ny];
            Double_t mceffval[ny], mcefferrl[ny], mcefferrh[ny];
            Double_t effval[ny],   efferrl[ny],   efferrh[ny];
            Double_t scaleval[ny], scaleerrl[ny], scaleerrh[ny];

            for(Int_t iy=1; iy<=ny; iy++) {
                xval[iy-1] = 0.5*(hMCEff->GetYaxis()->GetBinLowEdge(iy+1) + hMCEff->GetYaxis()->GetBinLowEdge(iy));
                xerr[iy-1] = 0.5*(hMCEff->GetYaxis()->GetBinLowEdge(iy+1) - hMCEff->GetYaxis()->GetBinLowEdge(iy));
                if(iy==ny && xval[iy] > 1000.0) {
                    xval[iy-1] = 75;
                    xerr[iy-1] = 25;
                }

                Double_t mceff  = hMCEff->GetCellContent(ix,iy);
                Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
                Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
                mceffval[iy-1]  = mceff;
                mcefferrl[iy-1] = mcerrl;
                mcefferrh[iy-1] = mcerrh;

                Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
                Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
                Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
                effval[iy-1]  = dataeff;
                efferrl[iy-1] = dataerrl;
                efferrh[iy-1] = dataerrh;


                Double_t scale = 0.0;
                if (mceff > 0.0) { 
                    scale = dataeff/mceff;
                    scaleval[iy-1]  = scale;
                    scaleerrl[iy-1] = (scale>0) ? scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff) : 0;
                    scaleerrh[iy-1] = (scale>0) ? scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff) : 0;
                } else {
                    scale = 0.0;
                    scaleval[iy-1]  = scale;
                    scaleerrl[iy-1] = 0.0;
                    scaleerrh[iy-1] = 0.0;
                }

                // fill smurf hists
                h2_smurf_MCEff      ->SetBinContent(iy, ix, mceffval[iy-1]);
                h2_smurf_MCErrLow   ->SetBinContent(iy, ix, mcefferrl[iy-1]);
                h2_smurf_MCErrHigh  ->SetBinContent(iy, ix, mcefferrh[iy-1]);

                h2_smurf_DataEff      ->SetBinContent(iy, ix, effval[iy-1]);
                h2_smurf_DataErrLow   ->SetBinContent(iy, ix, efferrl[iy-1]);
                h2_smurf_DataErrHigh  ->SetBinContent(iy, ix, efferrh[iy-1]);

                h2_smurf_SF         ->SetBinContent(iy, ix, scaleval[iy-1]);
                h2_smurf_SFErrLow   ->SetBinContent(iy, ix, scaleerrl[iy-1]);
                h2_smurf_SFErrHigh  ->SetBinContent(iy, ix, scaleerrh[iy-1]);

            }

            mceff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,mceffval,xerr,xerr,mcefferrl,mcefferrh));
            eff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,effval,xerr,xerr,efferrl,efferrh));
            scale_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,scaleval,xerr,xerr,scaleerrl,scaleerrh));
        }  
    }

    //--------------------------------------------------------------------------------------------------------------
    // write smurf file
    //==============================================================================================================

std::cout << "mark 4" << std::endl;

    outfile.Write();
    outfile.Close();

    if (mcfile)     delete mcfile;
    if (datafile)   delete datafile;

std::cout << "mark 5" << std::endl;

    //--------------------------------------------------------------------------------------------------------------
    // Make plots
    //==============================================================================================================
    TCanvas *c = MakeCanvas("c","c",800,600);

    if(grMCEffEta && grDataEffEta) {
        CPlot plotEffEta("effeta","","|#eta|","#varepsilon");
        plotEffEta.AddGraph(grMCEffEta,   "MC","",  kRed, kOpenSquare);
        plotEffEta.AddGraph(grDataEffEta,"data","",kBlue,kFullDotLarge);
        plotEffEta.SetYRange(efflow,effhigh);
        plotEffEta.TransLegend(0,-0.5);
        plotEffEta.Draw(c,kTRUE,format);

        CPlot plotScaleEta("scaleeta","","|#eta|","scale factor");
        plotScaleEta.AddGraph(grScaleEta,"",kBlue,kFullDotLarge);
        plotScaleEta.AddLine(0,1.0,2.7,1.0,kBlack,7);
        plotScaleEta.SetYRange(scalelow,scalehigh);
        plotScaleEta.SetXRange(0,2.7);
        plotScaleEta.Draw(c,kTRUE,format);
    }

    if(grMCEffPt && grDataEffPt) {
        CPlot plotEffPt("effpt","","p_{T} [GeV/c]","#varepsilon");
        plotEffPt.AddGraph(grMCEffPt,   "MC","",  kRed, kOpenSquare);
        plotEffPt.AddGraph(grDataEffPt,"data","",kBlue,kFullDotLarge);
        plotEffPt.SetYRange(efflow,effhigh);
        //plotEffPt.SetXRange(0,100);
        plotEffPt.TransLegend(0,-0.5);
        plotEffPt.Draw(c,kTRUE,format);

        CPlot plotScalePt("scalept","","p_{T} [GeV/c]","scale factor");
        plotScalePt.AddGraph(grScalePt,"",kBlue,kFullDotLarge);
        plotScalePt.AddLine(0,1.0,100,1.0,kBlack,7);
        plotScalePt.SetXRange(0,100);
        plotScalePt.SetYRange(scalelow,scalehigh);
        plotScalePt.Draw(c,kTRUE,format);
    }

    if(mceff_vs_pt_per_etav.size()>0 && eff_vs_pt_per_etav.size()>0) {
        for(UInt_t ig=0; ig<mceff_vs_pt_per_etav.size(); ig++) {
            char pname[100];
            sprintf(pname,"effpt_eta%i",ig);
            CPlot plotEffPt_perEta(pname,"","p_{T} [GeV/c]","#varepsilon");
            plotEffPt_perEta.AddGraph(mceff_vs_pt_per_etav[ig],"MC","",   kRed,  kOpenSquare);
            plotEffPt_perEta.AddGraph(eff_vs_pt_per_etav[ig], "data","", kBlue, kFullDotLarge);
            plotEffPt_perEta.AddTextBox(etalabelv[ig],0.2,0.82,0.4,0.88,0);
            plotEffPt_perEta.TransLegend(0,-0.5); 
            //plotEffPt_perEta.SetXRange(0,100);
            plotEffPt_perEta.SetYRange(efflow,effhigh);
            plotEffPt_perEta.Draw(c,kTRUE,format);
        }    

        for(UInt_t ig=0; ig<scale_vs_pt_per_etav.size(); ig++) {
            char pname[100];
            sprintf(pname,"scalept_eta%i",ig);
            CPlot plotScalePt_perEta(pname,"","p_{T} [GeV/c]","scale factor");
            plotScalePt_perEta.AddGraph(scale_vs_pt_per_etav[ig],"",kBlue,kFullDotLarge);
            plotScalePt_perEta.AddTextBox(etalabelv[ig],0.7,0.82,0.9,0.88,0);
            plotScalePt_perEta.AddLine(0,1,100,1,kBlack,7);
            //plotScalePt_perEta.SetXRange(0,100);
            plotScalePt_perEta.SetYRange(scalelow,scalehigh);
            plotScalePt_perEta.Draw(c,kTRUE,format);
        }
    }

    if(mceff_vs_pt_per_etav.size()==0 && eff_vs_pt_per_etav.size()>0) {
        for(UInt_t ig=0; ig<eff_vs_pt_per_etav.size(); ig++) {
            char pname[100];
            sprintf(pname,"effpt_eta%i",ig);
            CPlot plotEffPt_perEta(pname,"","p_{T} [GeV/c]","#varepsilon");
            plotEffPt_perEta.AddGraph(eff_vs_pt_per_etav[ig], "data","", kBlue, kFullDotLarge);
            plotEffPt_perEta.AddTextBox(etalabelv[ig],0.2,0.82,0.4,0.88,0);
            plotEffPt_perEta.TransLegend(0,-0.5);
            //plotEffPt_perEta.SetXRange(0,100);
            plotEffPt_perEta.SetYRange(efflow,effhigh);
            plotEffPt_perEta.Draw(c,kTRUE,format);
        }
    }

    if(mceff_vs_eta_per_ptv.size()>0 && eff_vs_eta_per_ptv.size()>0) {
        for(UInt_t ig=0; ig<mceff_vs_eta_per_ptv.size(); ig++) {
            char pname[100];
            sprintf(pname,"effeta_pt%i",ig);
            CPlot plotEffEta_perPt(pname,"","|#eta|","#varepsilon");
            plotEffEta_perPt.AddGraph(mceff_vs_eta_per_ptv[ig],"MC","",  kRed,  kOpenSquare);
            plotEffEta_perPt.AddGraph(eff_vs_eta_per_ptv[ig], "data","", kBlue, kFullDotLarge);
            plotEffEta_perPt.AddTextBox(ptlabelv[ig],0.2,0.82,0.4,0.88,0);
            plotEffEta_perPt.SetYRange(efflow,effhigh);
            plotEffEta_perPt.Draw(c,kTRUE,format);
        }

        for(UInt_t ig=0; ig<scale_vs_eta_per_ptv.size(); ig++) {
            char pname[100];
            sprintf(pname,"scaleeta_pt%i",ig);
            CPlot plotScaleEta_perPt(pname,"","|#eta|","scale factor");
            plotScaleEta_perPt.AddGraph(scale_vs_eta_per_ptv[ig],"", kBlue, kFullDotLarge);
            plotScaleEta_perPt.AddLine(-2.7,1.0,2.7,1.0,kBlack,7);
            plotScaleEta_perPt.AddTextBox(ptlabelv[ig],0.7,0.82,0.9,0.88,0);
            plotScaleEta_perPt.SetYRange(scalelow,scalehigh);
            plotScaleEta_perPt.SetXRange(-2.7,2.7);
            plotScaleEta_perPt.Draw(c,kTRUE,format);
        }
    }


    gBenchmark->Show("plotDataMC");
}
