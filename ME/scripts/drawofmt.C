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
#include "TArrow.h"
#include "TStyle.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TF1.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void drawsingle(TChain *& ch, int mH, TString type, float metcut, float min, float max);
void drawdatamc(TString dirName, int mH, float metcut, float min, float max);

void drawofmt()
{
  TString dirName = "data/";
  //TString dirName = "/smurf/ksung/ntuples/";

  int mH[6] = {250, 300, 350, 400, 500, 600};

  // draw comparison of ee/mm vs em

   // draw comparison of ee/mm vs em
  TChain *ch_top = new TChain("tree");

  ch_top->Add(dirName + "/ttbar2l-powheg.root");
  ch_top->Add(dirName + "/wtop-powheg.root");
  ch_top->Add(dirName + "/wtopb-powheg.root");
  ch_top->Add(dirName + "/ww-madgraph.root");


  if ( ch_top == 0x0 ) {
    std::cout << " ch_top is not filled exiting\n";
    return;
  }
  
  drawsingle(ch_top, 250, "OF", 70, 230, 300);
  drawdatamc(dirName, 250, 70, 230, 300);  
  /*
  drawsingle(ch_top, 300, "OF", 80, 250, 350);
  drawsingle(ch_top, 350, "OF", 80, 250, 400);
  drawsingle(ch_top, 400, "OF", 80, 250, 450);
  drawsingle(ch_top, 500, "OF", 80, 250, 600);
  drawsingle(ch_top, 600, "OF", 80, 300, 750);


  drawdatamc(dirName, 300, 80, 250, 350);
  drawdatamc(dirName, 350, 80, 250, 400);
  drawdatamc(dirName, 400, 80, 250, 450);
  drawdatamc(dirName, 500, 80, 250, 600);
  drawdatamc(dirName, 600, 80, 250, 700);
  */
  delete ch_top;
  
}
    
void drawsingle(TChain *& ch, int mH, TString type, float metcut, float min, float max)
{

  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();

  TString plotDir = "tempplots/ofmtvar/nonbtag/";
  TString fileName;
  TString modeName;

  float met_ = 0.;
  float scale1fb_ = 0.;
  LorentzVector*  dilep_ = 0;
  int type_ = 0; // 0/1/2/3 for mm/me/em/ee
  float metPhi_ = 0.;
  LorentzVector*  jet1_ = 0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  unsigned int nSoftMuons_ = 0;
  float jet1Btag_ = 0;
  float jet2Btag_ = 0;
  float jet3Btag_ = 0;


  ch->SetBranchAddress( "met",  &met_);
  ch->SetBranchAddress( "metPhi",  &metPhi_);
  ch->SetBranchAddress( "scale1fb",  &scale1fb_);
  ch->SetBranchAddress( "dilep", &dilep_); 
  ch->SetBranchAddress( "type"      , &type_     );     
  ch->SetBranchAddress( "jet1", &jet1_);
  ch->SetBranchAddress( "jet2", &jet2_);
  ch->SetBranchAddress( "jet3", &jet3_);
  ch->SetBranchAddress( "nSoftMuons"      , &nSoftMuons_     );     
  ch->SetBranchAddress( "jet1Btag", &jet1Btag_);
  ch->SetBranchAddress( "jet2Btag", &jet2Btag_);
  ch->SetBranchAddress( "jet3Btag", &jet3Btag_);
  
  TH1F *hist_ee = new TH1F("hist_ee", "hist_ee", 10, min, max);
  TH1F *hist_mm = new TH1F("hist_mm", "hist_mm", 10, min, max);
  TH1F *hist_em = new TH1F("hist_em", "hist_em", 10, min, max);
  
  hist_ee->Sumw2();
  hist_mm->Sumw2();
  hist_em->Sumw2();
  
   // loop over the events to fill histograms
   for(int ievt = 0; ievt< ch->GetEntries(); ievt++){  
     ch->GetEntry(ievt);  
     if ( met_ < metcut) continue;
     if ( TMath::Abs(dilep_->M() - 91.1876) > 15) continue;
     if ( dilep_->Pt() < 55.0) continue;
     
     if (jet1Btag_ > 2.0 && jet1_->Pt() > 30) continue;
     if (jet2Btag_ > 2.0 && jet2_->Pt() > 30) continue;
     if (jet3Btag_ > 2.0 && jet3_->Pt() > 30) continue;
     if ( nSoftMuons_ != 0) continue;     
     
     
     // dphi cut
     
     float dphiCut = 0.5;
     
     float dPhi1, dPhi2, dPhi3;
     jet1_->Pt() > 30 ? dPhi1 = acos(cos(metPhi_ - jet1_->Phi())) : dPhi1 = 999.9;
     jet2_->Pt() > 30 ? dPhi2 = acos(cos(metPhi_ - jet2_->Phi())) : dPhi2 = 999.9;
     jet3_->Pt() > 30 ? dPhi3 = acos(cos(metPhi_ - jet3_->Phi())) : dPhi3 = 999.9;
     
     // apply the dphi Cut
     if ( TMath::Min(dPhi1, TMath::Min(dPhi2, dPhi3)) < dphiCut ) continue;

     // calculate the mT
     float termA = sqrt( dilep_->Pt()*dilep_->Pt() + dilep_->M() * dilep_->M() );
     float termB = sqrt( met_*met_ + dilep_->M() * dilep_->M() );
     float newX = dilep_->Px() + met_ * cos(metPhi_);
     float newY = dilep_->Py() + met_ * sin(metPhi_);
     float termC = newX*newX + newY*newY;
     float mt = sqrt( pow(termA + termB, 2) - termC );
     if ( mt > max || mt < min) continue;
     
     if ( type_ == 0 )   hist_mm->Fill(mt, scale1fb_);
     if ( type_ == 3 )    hist_ee->Fill(mt, scale1fb_);
     if ( type_ == 1 || type_ == 2) hist_em->Fill(mt, scale1fb_);
   }
   
   // normalise all plots to 1.
   
   hist_ee->Scale(1./hist_ee->Integral(0, 9999));
   hist_mm->Scale(1./hist_mm->Integral(0, 9999));
   hist_em->Scale(1./hist_em->Integral(0, 9999));

   hist_ee->SetLineColor(kBlue);
   hist_ee->SetMarkerColor(kBlue);

   hist_mm->SetLineColor(kRed);
   hist_mm->SetMarkerColor(kRed);
   
   
   float yMax = hist_ee->GetMaximum();
   yMax = hist_mm->GetMaximum() > yMax ? hist_mm->GetMaximum() : yMax;
   yMax = hist_em->GetMaximum() > yMax ? hist_em->GetMaximum() : yMax;
   
   hist_ee->SetMaximum( yMax * 1.2);
   hist_mm->SetMaximum( yMax * 1.2);
   hist_em->SetMaximum( yMax * 1.2);

   hist_ee->SetMinimum( 0.001);
   hist_mm->SetMinimum( 0.001);
   hist_em->SetMinimum( 0.001);
   
   hist_ee->SetXTitle("mT (GeV)");
   hist_mm->SetXTitle("mT (GeV)");
   hist_em->SetXTitle("mT (GeV)");
   
   // get legend
   TLegend *leg = new TLegend(0.55, 0.7, 0.85, 0.9);
   leg->SetFillColor(0);
   leg->AddEntry(hist_ee, "ee");
   leg->AddEntry(hist_mm, "#mu#mu");
   leg->AddEntry(hist_em, "e#mu");


   // draw plots..
   TCanvas *c1 = new TCanvas();

   hist_mm->Draw("HISTE0");
   hist_ee->Draw("samehiste0");
   hist_em->Draw("samehiste0");
   leg->Draw("same");
   
   c1->SaveAs(Form("%s/%s_mT_mH%i_lin.png", plotDir.Data(), type.Data(), mH));
   c1->SaveAs(Form("%s/%s_mT_mH%i_lin.eps", plotDir.Data(), type.Data(), mH));
   
   c1->SetLogy();
   
   c1->SaveAs(Form("%s/%s_mT_mH%i_log.png", plotDir.Data(), type.Data(), mH));
   c1->SaveAs(Form("%s/%s_mT_mH%i_log.eps", plotDir.Data(), type.Data(), mH));

   // tidy up
   delete c1;
   delete hist_ee;
   delete hist_mm;
   delete hist_em;
   
}

void drawdatamc(TString dirName, int mH, float metcut, float min, float max)
{  
  
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();

  TString plotDir = "tempplots/ofmtvar/nonbtag/";

  TChain *ch = new TChain("tree");

  ch->Add(dirName + "data_2l.goodlumiRun2011AB.root");
  ch->Add(dirName + "/ttbar2l-powheg.root");
  ch->Add(dirName + "/wtop-powheg.root");
  ch->Add(dirName + "/wtopb-powheg.root");
  ch->Add(dirName + "/ww-madgraph.root");

  if ( ch == 0x0 ) {
    std::cout << " ch is not filled exiting\n";
    return;
  }
  
  TH1F *hist_data = new TH1F("hist_data", "hist_data", 10, min, max);
  hist_data->Sumw2();
  TH1F *hist_mc = new TH1F("hist_mc", "hist_mc", 10, min, max);
  hist_mc->Sumw2();
  
  float met_ = 0.;
  float scale1fb_ = 0.;
  LorentzVector*  dilep_ = 0;
  int type_ = 0; // 0/1/2/3 for mm/me/em/ee
  float metPhi_ = 0.;
  LorentzVector*  jet1_ = 0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  unsigned int nSoftMuons_ = 0;
  float jet1Btag_ = 0;
  float jet2Btag_ = 0;
  float jet3Btag_ = 0;
  int dstype_ = 0; // 0 for data

  ch->SetBranchAddress( "met",  &met_);
  ch->SetBranchAddress( "metPhi",  &metPhi_);
  ch->SetBranchAddress( "scale1fb",  &scale1fb_);
  ch->SetBranchAddress( "dilep", &dilep_); 
  ch->SetBranchAddress( "type"      , &type_     );     
  ch->SetBranchAddress( "jet1", &jet1_);
  ch->SetBranchAddress( "jet2", &jet2_);
  ch->SetBranchAddress( "jet3", &jet3_);
  ch->SetBranchAddress( "nSoftMuons"      , &nSoftMuons_     );     
  ch->SetBranchAddress( "jet1Btag", &jet1Btag_);
  ch->SetBranchAddress( "jet2Btag", &jet2Btag_);
  ch->SetBranchAddress( "jet3Btag", &jet3Btag_);
  ch->SetBranchAddress( "dstype"      , &dstype_     );     
  
  // loop over the events to fill histograms
  for(int ievt = 0; ievt< ch->GetEntries(); ievt++){  
     ch->GetEntry(ievt);  
     if ( met_ < metcut) continue;
     if ( TMath::Abs(dilep_->M() - 91.1876) > 15) continue;
     if ( TMath::Abs(dilep_->Pt() < 55)) continue;
     if ( type_ == 0 ) continue;
     if ( type_ == 3) continue;
 
     // dphi cut
     float dphiCut = 0.5;
     
     float dPhi1, dPhi2, dPhi3;
     jet1_->Pt() > 30 ? dPhi1 = acos(cos(metPhi_ - jet1_->Phi())) : dPhi1 = 999.9;
     jet2_->Pt() > 30 ? dPhi2 = acos(cos(metPhi_ - jet2_->Phi())) : dPhi2 = 999.9;
     jet3_->Pt() > 30 ? dPhi3 = acos(cos(metPhi_ - jet3_->Phi())) : dPhi3 = 999.9;
     
     // apply the dphi Cut
     if ( TMath::Min(dPhi1, TMath::Min(dPhi2, dPhi3)) < dphiCut ) continue;

     // calculate the mT
     float termA = sqrt( dilep_->Pt()*dilep_->Pt() + dilep_->M() * dilep_->M() );
     float termB = sqrt( met_*met_ + dilep_->M() * dilep_->M() );
     float newX = dilep_->Px() + met_ * cos(metPhi_);
     float newY = dilep_->Py() + met_ * sin(metPhi_);
     float termC = newX*newX + newY*newY;
     float mt = sqrt( pow(termA + termB, 2) - termC );
     if ( mt > max || mt < min) continue;
   
     if ( dstype_ == 0)  {
       hist_data->Fill(mt);
     }
     else {
       hist_mc->Fill(mt, scale1fb_);
     }
  }

  
  hist_mc->Scale(hist_data->Integral(0, 1000) / hist_mc->Integral(0, 1000) );
  

  // draw overlay 
  float yMax = hist_mc->GetMaximum() > hist_data->GetMaximum() ? hist_mc->GetMaximum() : hist_data->GetMaximum();
  hist_mc->SetMaximum(yMax * 1.1);
  hist_data->SetMaximum(yMax * 1.1);
  
  
  // draw plots..
  hist_data->SetLineColor(kBlue);
  hist_data->SetMarkerColor(kBlue);
  hist_mc->SetLineColor(kRed);
  hist_mc->SetMarkerColor(kRed);
  
  hist_data->SetXTitle("mT (GeV)");
  hist_mc->SetXTitle("mT (GeV)");

  
  // get legend
  TLegend *leg = new TLegend(0.55, 0.7, 0.85, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(hist_data, "Data (e#mu)");
  leg->AddEntry(hist_mc, "MC (e#mu)");

  TH1F* hist_ratio = (TH1F*) hist_data->Clone("hist_ratio");
  hist_ratio->Divide(hist_ratio, hist_mc, 1, 1);
  hist_ratio->SetYTitle("Data/MC");
  hist_ratio->GetXaxis()->SetTitleSize(0.15);
  hist_ratio->GetYaxis()->SetTitleSize(0.12);
  hist_ratio->GetYaxis()->SetTitleOffset(0.5);
  hist_ratio->GetXaxis()->SetLabelSize(0.1);
  hist_ratio->GetYaxis()->SetLabelSize(0.1);
  hist_ratio->GetYaxis()->SetRangeUser(0.0, 3.0);
  hist_ratio->GetYaxis()->SetNdivisions(504);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 800);
  c1->cd();
  
  TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.05);
  pad1->SetRightMargin(0.07);
  pad1->SetLeftMargin(0.18);
  pad1->Draw();
  
  c1->cd();
  TPad *pad2 = new TPad("p_leg", "p_leg", 0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.18);
  pad2->SetRightMargin(0.07);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();

  pad1->cd();
  hist_data->Draw("HISTE0");
  hist_mc->Draw("SAMEE0");
  leg->Draw("same");

  pad2->cd();  
  TLine *line = new TLine(min, 1, max, 1);
  line->SetLineWidth(2);
  TLine *line_up = new TLine(min, 1.25, max, 1.25);
  line_up->SetLineWidth(2);
  line_up->SetLineStyle(2);
  TLine *line_down = new TLine(min, 0.75, max, 0.75);
  line_down->SetLineWidth(2);
  line_down->SetLineStyle(2);
  
  hist_ratio->Draw();
  line->Draw("same");
  line_up->Draw("same");
  line_down->Draw("same"); 

  c1->SaveAs(Form("%s/OF_mT_mH%i_datamc_lin.png", plotDir.Data(), mH));
  c1->SaveAs(Form("%s/OF_mT_mH%i_datamc_lin.eps", plotDir.Data(), mH));
  
  pad1->SetLogy();
  c1->SaveAs(Form("%s/OF_mT_mH%i_datamc_log.png", plotDir.Data(), mH));
  c1->SaveAs(Form("%s/OF_mT_mH%i_datamc_log.eps", plotDir.Data(), mH));

  

  // tidy up
  delete ch;
  delete c1;
  delete hist_data;
  delete hist_mc;
  delete hist_ratio;
}


