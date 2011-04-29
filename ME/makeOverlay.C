//
// Plotting script to overlay various Higgs signal with the SM background
// D. Evans and Y. Gao
// ==================
// Instructions
// 1. Create a soft link called data where the full statistics smurfNtuples live
// 2. By default the smurfNtuples containing the ME output is the current directory
// 3. Specify the higgs samples to consider in the main function makeOverlay()
// 4. Do root -l makeOverlay.C+
// 

#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <math.h>
#include <iostream>

// Matrix Element related
#include "TVar.hh"

const unsigned int kNDilep = 4;

class Sample {

    public:

        Sample(TVar::Process process, Color_t color, bool stack) { 

            // initialize data members
            process_ = process;
            stack_ = stack;
            color_ = color;
            file_me_ = TFile::Open(TVar::SmurfProcessName(process)+"_ME.root");
            if (file_me_ == 0x0) {
                std::cout << "Sample (" << TVar::SmurfProcessName(process_) << ") ME file is BAD" << std::endl;
            }
            else {
                tree_me_ = (TTree*)file_me_->Get("tree");
            }

            // get the yields for the normalisations
            InitializeYields();

        }
        ~Sample() {
            file_me_->Close();
        };

        TString GetProcessName() {
            return TVar::SmurfProcessName(process_);
        }

        void SetYield(float yield, unsigned int j) { 
            yield_[j] = yield; 
        }
        const float &GetYield(unsigned int j) { 
            return yield_[j]; 
        }
        float GetTotalYield() {
            float total = 0.0;
            for (unsigned int j = 0; j < kNDilep; ++j) total += yield_[j];
            return total;
        }

        TH1F *GetHistogram(TVar::Process processLR, Int_t nBins, Float_t binMin, Float_t binMax) {
	  TH1F *tmp = new TH1F(Form("h1_%sLR_%s",  TVar::SmurfProcessName(processLR).Data(), GetProcessName().Data()), "LR", nBins, binMin, binMax);
	  tree_me_->Project(Form("h1_%sLR_%s", TVar::SmurfProcessName(processLR).Data(), GetProcessName().Data()), Form("LR[%i]", processLR));
            tmp->Scale(GetTotalYield()/tmp->Integral(0, nBins+1));
            if (stack_) {
                tmp->SetFillColor(color_);
                tmp->SetLineColor(color_);
                // tmp->SetFillStyle(1);
            }
            else {
                tmp->SetLineColor(color_);
                tmp->SetMarkerColor(color_);
                tmp->SetLineWidth(3);
            }
            return tmp;
        }

        bool Stack() { return stack_; }

    private:

        //
        // member functions
        //

        void InitializeYields() 
        {

            TFile *file_norm = TFile::Open("data/"+GetProcessName()+".root");
            TTree *tree_norm = (TTree*)file_norm->Get("tree");

            if (file_norm == 0x0) {
                std::cout << "Sample (" << GetProcessName() << ") NORM file is bad" << std::endl;
                return;
            }

            gROOT->cd();
            for (unsigned int j = 0; j < kNDilep; j++)
            {
                TH1F *tmp = new TH1F("tmp", "tmp", 20, 0, 4);
                if(j==0 || j==2)
                    tree_norm->Project("tmp", "dPhi", Form("scale1fb*(type==%i)", j));
                if(j==1 || j==3)
                    tree_norm->Project("tmp", "dPhi", Form("scale1fb*(type==%i&&lep2.pt()>15)", j));
                if (tmp != 0x0) yield_[j] = tmp->Integral(0, 9999);
                else yield_[j] = 0.0;
                delete tmp;
            }

            file_norm->Close();

        }

        //
        // member data
        //

        TFile *file_me_;
        TTree *tree_me_;
        TVar::Process process_;
        bool stack_;
        Color_t color_;
        float yield_[kNDilep];
};

void drawOverlay(TVar::Process higgsProcess)
{  
  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  std::vector<Sample*> samples;
  samples.push_back(new Sample(higgsProcess, kRed, false));
  samples.push_back(new Sample(TVar::ttbar, kMagenta, true)); 
  samples.push_back(new Sample(TVar::Wp_1jet, kCyan, true));
  samples.push_back(new Sample(TVar::WW, kYellow+2, true)); 
  TH1F *signalHistogram = 0;

  
  THStack *h1_stack = new THStack(); 
  TLegend *lg = new TLegend(0.62, 0.6, 0.88, 0.9);
  for (unsigned int s = 0; s < samples.size(); ++s)
    {
      if (samples[s]->Stack()) {
	TH1F *tmp = samples[s]->GetHistogram(higgsProcess, 20, 0.0, 1.0);
	tmp->SetXTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsProcess).Data()));
	h1_stack->Add(tmp);
	lg->AddEntry(tmp, samples[s]->GetProcessName(), "f");
      } else {
	signalHistogram = samples[s]->GetHistogram(higgsProcess, 20, 0.0, 1.0);
      }
    }
  lg->AddEntry(signalHistogram, TVar::SmurfProcessName(higgsProcess), "lp");
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetShadowColor(0);
  
  signalHistogram->SetXTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsProcess).Data()));
  // now display them
  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1_stack->Draw("HIST");
  h1_stack->GetXaxis()->SetTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsProcess).Data()));
  signalHistogram->Draw("SAMEHIST");
  lg->Draw();
  c1->Print( TVar::SmurfProcessName(higgsProcess)+"_LR.eps");
  c1->Print( TVar::SmurfProcessName(higgsProcess)+"_LR.png");
  
  // tidy up
  for (unsigned int s = 0; s < samples.size(); s++) delete samples[s];
  delete signalHistogram;
  delete h1_stack;
  delete c1;
  delete lg;
}

void makeOverlay()
{
  std::vector<TVar::Process> HiggsSamples;
  HiggsSamples.push_back(TVar::HWW160); 
  HiggsSamples.push_back(TVar::HWW120);
  HiggsSamples.push_back(TVar::HWW200); 
  
  for(unsigned int s = 0; s < HiggsSamples.size(); ++s)
    drawOverlay(HiggsSamples[s]);
      
}

