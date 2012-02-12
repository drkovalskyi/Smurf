#include<vector>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THStack.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TArrow.h"
#include <iostream>
#endif

enum samp { iHWW, iWW, iZJets, iTop, iVV, iWJets, iWZ, iZZ, iFakes, nSamples };

float xPos[nSamples+1] = {0.22,0.22,0.22,0.39,0.39,0.39,0.39,0.22,0.39,0.39}; 
float yOff[nSamples+1] = {0,1,2,0,1,2,3,0,1,2};

const Float_t _tsize   = 0.03;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;


//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F* h)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
        TString coordinate,
        TString title)
{
    axis->SetLabelFont  (   42);
    axis->SetLabelOffset(0.015);
    axis->SetLabelSize  (0.050);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (   42);
    axis->SetTitleOffset( 1.15);
    axis->SetTitleSize  (0.050);

    if (coordinate == "y") axis->SetTitleOffset(1.4);
    if (coordinate == "x") axis->SetTitleSize(0.045);
    if (coordinate == "x") axis->SetLabelSize(0.045);

    axis->SetTitle(title);
}


//------------------------------------------------------------------------------
// THStackAxisFonts
//------------------------------------------------------------------------------
void THStackAxisFonts(THStack* h,
        TString  coordinate,
        TString  title)
{
    TAxis* axis = NULL;
    if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
    if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();
    AxisFonts(axis, coordinate, title);
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
        Float_t y1,
        TH1F*   hist,
        TString label,
        TString option)
{
    TLegend* legend = new TLegend(x1,
            y1,
            x1 + _xoffset,
            y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (_tsize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}


class StandardPlot {

    public: 
        StandardPlot() { _hist.resize(nSamples,0); _data = 0; _breakdown = false; _mass = 0; 
          _lowCutValue = -999; _highCutValue = -999; }

        void initHist(TH1F *h) {
          if (h) {
            h->GetXaxis()->SetTitleOffset(1.15);
            h->GetYaxis()->SetTitleOffset(1.4);
          }
        }
        void setMCHist   (const samp &s, TH1F * h)  { initHist(h); _hist[s]       = h;  }
        void setDataHist (TH1F * h)                 { initHist(h); _data          = h;  } 
        void setHWWHist  (TH1F * h)                 { initHist(h); setMCHist(iHWW  ,h); } 
        void setWWHist   (TH1F * h)                 { initHist(h); setMCHist(iWW   ,h); } 
        void setZJetsHist(TH1F * h)                 { initHist(h); setMCHist(iZJets,h); } 
        void setTopHist  (TH1F * h)                 { initHist(h); setMCHist(iTop  ,h); } 
        void setVVHist   (TH1F * h)                 { initHist(h); setMCHist(iVV   ,h); } 
        void setWZHist   (TH1F * h)                 { initHist(h); setMCHist(iWZ   ,h); } 
        void setZZHist   (TH1F * h)                 { initHist(h); setMCHist(iZZ   ,h); } 
        void setFakesHist(TH1F * h)                 { initHist(h); setMCHist(iFakes,h); } 
        void setWJetsHist(TH1F * h)                 { initHist(h); setMCHist(iWJets,h); }
        void setMass(const int &m) {_mass=m;}
        void setLowCutValue (double x) {_lowCutValue=x;}
        void setHighCutValue (double x) {_highCutValue=x;}

        TH1F* getDataHist() { return _data; }
//         TH1* DrawAndRebinTo(const int &rebinTo) {

//             if(rebinTo == 0) return Draw();
//             int rebin = 0, nbins = 0;
//             for (int i=0; i<nSamples; i++) {

//                 // in case the user doesn't set it
//                 if( !_hist[i] ) continue;

//                 nbins = _hist[i]->GetNbinsX();
//             }
//             if (nbins == 0) return Draw();

//             rebin = nbins / rebinTo;
//             while(nbins % rebin != 0) rebin--;
//             return Draw(rebin);

//         }

    void Draw(string outputName, Bool_t isLogY, const int &rebin=1) {
      TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
      if(isLogY == true) c1->SetLogy();


      Color_t _sampleColor[nSamples];
      _sampleColor[iHWW  ] = kRed+1;
      _sampleColor[iWW   ] = kAzure-9;
      _sampleColor[iZJets] = kGreen+2;
      _sampleColor[iTop  ] = kYellow;
      _sampleColor[iVV   ] = kAzure-2;
      _sampleColor[iWJets] = kGray+1;
      _sampleColor[iWZ   ] = kAzure-2;
      _sampleColor[iZZ   ] = kAzure-9;
      _sampleColor[iFakes] = kGray+1;
      //             _sampleColor[iWJets] = kViolet-9;
      //             _sampleColor[iWJets] = kCyan;

      //setUpStyle();
      //if(!gPad) new TCanvas();

      THStack* hstack = new THStack();
      for (int i=0; i<nSamples; i++) {

        // in case the user doesn't set it
        if( !_hist[i] ) continue;

        _hist[i]->Rebin(rebin);
        _hist[i]->SetLineColor(_sampleColor[i]);

        // signal gets overlaid
        if (i == iHWW) continue;

        _hist[i]->SetFillColor(_sampleColor[i]);
        _hist[i]->SetFillStyle(1001);

        initHist(_hist[i]);
        _hist[i]->GetXaxis()->SetTitleOffset(1.15);
        _hist[i]->GetYaxis()->SetTitleOffset(1.4);
        hstack->Add(_hist[i]);
      }

      if(_hist[iHWW]) _hist[iHWW]->SetLineWidth(3);
      if(_data) _data->Rebin(rebin);
      if(_data) _data->SetLineColor  (kBlack);
      if(_data) _data->SetMarkerStyle(kFullCircle);

      hstack->Draw("hist");
//       if(_hist[iHWW]) _hist[iHWW]->Draw("hist");
      hstack->Draw("hist");
      hstack->GetHistogram()->GetXaxis()->SetTitleOffset(1.15);
      hstack->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
      if(_hist[iHWW]) _hist[iHWW]->Draw("hist,same");
      if(_data) _data->Draw("ep,same");

//             hstack->SetTitle("CMS preliminary");
      hstack->SetTitle("");
      hstack->GetXaxis()->SetTitleOffset(1.15);
      hstack->GetYaxis()->SetTitleOffset(1.4);


      //Define Correct YScale
      Float_t theMax = hstack->GetMaximum();
      Float_t theMin = hstack->GetMinimum();
      Float_t maxY = theMax;
      if (_hist[iHWW]) {
        if (_hist[iHWW]->GetMaximum() > theMax) theMax = _hist[iHWW]->GetMaximum();
        if (_hist[iHWW]->GetMinimum() < theMin) theMin = _hist[iHWW]->GetMinimum();
      }            
      if (_data) {              
        Float_t dataMax = GetMaximumIncludingErrors(_data);              
        if (dataMax > theMax) theMax = dataMax;
      }            
      if (gPad->GetLogy()) {
        hstack->SetMaximum(50 * theMax);
        hstack->SetMinimum(0.02);
        maxY = 50 * theMax;
      } else {
        hstack->SetMaximum(1.55 * theMax);
        maxY = 1.55 * theMax;
      }


      //Axis Titles
      if(_breakdown) {
        THStackAxisFonts(hstack, "y", "entries");
        hstack->GetHistogram()->LabelsOption("v");
      } else {
              
        THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
        if(_units.Sizeof() == 1) {
          THStackAxisFonts(hstack, "x", TString::Format("%s",_xLabel.Data()));
          THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f",hstack->GetHistogram()->GetXaxis()->GetBinWidth(1)));
        } else {
          THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
          THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f %s",hstack->GetHistogram()->GetXaxis()->GetBinWidth(1),_units.Data()));
        }
      }
            
      // total mess to get it nice, should be redone
      size_t j=0;
      TString higgsLabel = " HWW";
      if(_mass != 0) higgsLabel.Form(" m_{H}=%d",_mass);

      if(_data        ) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data,         " data",    "lp"); j++; }
      if(_hist[iHWW  ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHWW  ], higgsLabel, "l" ); j++; }
      if(_hist[iWW   ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWW   ], " WW",      "f" ); j++; }
      if(_hist[iZJets]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZJets], " Z+jets",  "f" ); j++; }
      if(_hist[iTop  ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iTop  ], " top",     "f" ); j++; }
      if(_hist[iVV   ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVV   ], " WZ/ZZ",   "f" ); j++; }
      if(_hist[iWJets]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWJets], " W+jets",  "f" ); j++; }
      if(_hist[iWZ   ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWZ   ], " WZ",      "f" ); j++; }
      if(_hist[iZZ   ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZZ   ], " ZZ",      "f" ); j++; }
      if(_hist[iFakes]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iFakes], " fakes",   "f" ); j++; }

      TLatex* luminosity = new TLatex(0.9, 0.815, TString::Format("L = %.2f fb^{-1}",_lumi));
      luminosity->SetNDC();
      luminosity->SetTextAlign(32);
      luminosity->SetTextFont(42);
      luminosity->SetTextSize(_tsize);
      luminosity->Draw("same");
      if(_extraLabel) _extraLabel->Draw("same");

      TLatex* energy = new TLatex(0.89, 0.86, TString::Format("#sqrt{s} = 7 TeV"));
      energy->SetNDC();
      energy->SetTextAlign(32);
      energy->SetTextFont(42);
      energy->SetTextSize(_tsize);
      energy->Draw("same");

      //lines to indicate cuts
      if (_lowCutValue != -999) {
        TLine *lowcutline = new TLine(_lowCutValue,maxY*0.75,_lowCutValue,0.0);
        if (isLogY == true)
          lowcutline = new TLine(_lowCutValue,maxY*0.05,_lowCutValue,0.0);
        lowcutline->SetLineColor(kBlue);
        lowcutline->SetLineStyle(2);
        lowcutline->SetLineWidth(4);
        lowcutline->Draw();

        double xaxisLength = hstack->GetXaxis()->GetXmax()-hstack->GetXaxis()->GetXmin();
        TArrow *lowcutarrow = new TArrow(_lowCutValue,
                                         maxY*0.70,
                                         _lowCutValue + xaxisLength*0.05,
                                         maxY*0.70,0.03,"|>");
        if (isLogY == true)
          lowcutarrow = new TArrow(_lowCutValue,
                                   maxY*0.03,
                                   _lowCutValue + xaxisLength*0.05,
                                   maxY*0.03,0.03,"|>");
        
        lowcutarrow->SetLineColor(kBlue);
        lowcutarrow->SetLineWidth(2);
        lowcutarrow->SetArrowSize(0.01);
        lowcutarrow->SetAngle(30);
        lowcutarrow->Draw();
        

      }

      if (_highCutValue != -999) {
        TLine *highcutline = new TLine(_highCutValue,maxY*0.75,_highCutValue,0.0);
        if (isLogY == true)
          highcutline = new TLine(_highCutValue,maxY*0.05,_highCutValue,0.0);

        highcutline->SetLineColor(kBlue);
        highcutline->SetLineStyle(2);
        highcutline->SetLineWidth(4);
        highcutline->Draw();

        double xaxisLength = hstack->GetXaxis()->GetXmax()-hstack->GetXaxis()->GetXmin();
        TArrow *highcutarrow = new TArrow(_highCutValue,
                                          maxY*0.70,
                                          _highCutValue - xaxisLength*0.05,
                                          maxY*0.70,0.03,"|>");
        if (isLogY == true)
          highcutarrow = new TArrow(_highCutValue,
                                    maxY*0.03,
                                    _highCutValue - xaxisLength*0.05,
                                    maxY*0.03,0.03,"|>");
        
        highcutarrow->SetLineColor(kBlue);
        highcutarrow->SetLineWidth(2);
        highcutarrow->SetArrowSize(0.01);
        highcutarrow->SetAngle(30);
        highcutarrow->Draw();
        

      }


      c1->SaveAs((outputName+".png").c_str());
      c1->SaveAs((outputName+".eps").c_str());
//       c1->SaveAs((outputName+".pdf").c_str());

    }

    void setLumi(const float &l) { _lumi = l; }
    void setLabel(const TString &s) { _xLabel = s; }
    void setUnits(const TString &s) { _units = s; }
    void setBreakdown(const bool &b = true) { _breakdown = b; }
    void addLabel(const std::string &s) {
      _extraLabel = new TLatex(0.9, 0.77, TString(s));
      _extraLabel->SetNDC();
      _extraLabel->SetTextAlign(32);
      _extraLabel->SetTextFont(42);
      _extraLabel->SetTextSize(_tsize);
    }

    private: 
        std::vector<TH1F*> _hist;
        TH1F* _data;

        //MWL
        float    _lumi;
        TString  _xLabel;
        TString  _units;
        TLatex * _extraLabel;
        bool     _breakdown;
        int      _mass;
        double _lowCutValue;
        double _highCutValue;


};


