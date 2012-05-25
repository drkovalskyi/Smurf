#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->CompileMacro("tauify.cc","k");
  //  gSystem->CompileMacro("../../Core/SmurfTree.h","k");
  gSystem->CompileMacro("dytt.C","k");
  // dytt();
}
#endif 

/*
void dytt(){
  TauData* taudata = new TauData();
  double pt = 10;
  TH1D* h = new TH1D("h","h",100,0,pt);
  LorentzVector v(pt,0,0,sqrt(pt*pt+1.777*1.777));
  for ( unsigned int i=0; i<1e5; ++i ){
    taudata->Tauify(v,0);
    h->Fill(taudata->GetLeptonFromTau().pt());
    // cout << "Lepton E: " << taudata->GetLeptonFromTau().E() <<endl;
    // cout << "\tLepton Pt: " << taudata->GetLeptonFromTau().pt() <<endl;
    // cout << "Missing E: " << taudata->GetMetFromTau().E() <<endl;
    // cout << "Missing Et: " << taudata->GetMetFromTau().Et() <<endl;
    // cout << "\tMissing Pt: " << taudata->GetMetFromTau().pt() <<endl;
    h->Draw();
  }
}
*/
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "../../Core/SmurfTree.h"
// #include "smurfAnalysis.h"
#include <fstream>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"

#ifndef __CINT__
#include <algorithm>
#include "tauify.h"
double projectedMet(const LorentzVector& met, 
		    const LorentzVector& lep1, 
		    const LorentzVector& lep2)
{
  double phi = met.phi();
  double dPhi1 = fabs(lep1.phi() - phi);
  dPhi1 = std::min(2*M_PI - dPhi1, dPhi1);
  double dPhi2 = fabs(lep2.phi() - phi);
  dPhi2 = std::min(2*M_PI - dPhi2, dPhi2);
  double dPhi = std::min(dPhi1,dPhi2);
  if (dPhi < M_PI/2) return met.pt()*sin(dPhi);
  return met.pt();
}
double pow2(const double& a){
  return a*a;
}
double mt(LorentzVector dilep, LorentzVector met){
  return sqrt(pow2(dilep.pt()+met.pt())-pow2(dilep.px()+met.px())-pow2(dilep.py()+met.py()));
}
double mt(double pt, double met, double dphi)
{
  return sqrt(2*pt*met*(1-cos(dphi)));
}
double integralError(const TH1* h, int firstBin, int lastBin){
  double sum2(0);
  for (int i=firstBin; i<=lastBin; ++i)
    sum2 += h->GetBinError(i)*h->GetBinError(i);
  return sqrt(sum2);
}


#endif


void dytt(const unsigned int nToys = 100)
{
  std::cout << "nToys: " << nToys << std::endl;
  const unsigned int njets = 0;
  const unsigned int prescale = 1;
  // double bfscale = 0.0620/2*prescale;
  double bfscale = pow(0.1739+0.1782,2)/4*prescale; //Zmm to Ztt(em - one channel)
  const double lumi = 1.55;
  const int spin1 = 0;
  const int spin2 = 0;
  const bool scaleMCByArea = false;
  const double dymm_cor = 1.30;
  const double dyee_cor = 0.86;
  const double minMT = -1;
  const bool useMCWeights = true;
  const bool useMassCut = true;
  // const double dymm_cor = 1.35;
  // const double dyee_cor = 0.90;

  TH1D* hztt_mass_mc = new TH1D("hztt_mass_mc","Drell-Yan(#tau#tau) Mass",200,0,200);
  hztt_mass_mc->SetDirectory(0);
  hztt_mass_mc->Sumw2();
  hztt_mass_mc->SetMarkerColor(kRed);
  hztt_mass_mc->SetMarkerStyle(20);
  TH1D* hztt_met0_mc = new TH1D("hztt_met0_mc","Drell-Yan(#tau#tau) projected(pfMET)",50,0,100);
  hztt_met0_mc->SetDirectory(0);
  hztt_met0_mc->Sumw2();
  hztt_met0_mc->SetMarkerColor(kRed);
  hztt_met0_mc->SetMarkerStyle(20);
  TH1D* hztt_met1_mc = new TH1D("hztt_met1_mc","Drell-Yan(#tau#tau) projected(trkMET)",50,0,100);
  hztt_met1_mc->SetDirectory(0);
  hztt_met1_mc->Sumw2();
  hztt_met1_mc->SetMarkerColor(kRed);
  hztt_met1_mc->SetMarkerStyle(20);
  TH1D* hztt_met2_mc = new TH1D("hztt_met2_mc","Drell-Yan(#tau#tau) minProjectedMET",50,0,100);
  hztt_met2_mc->SetDirectory(0);
  hztt_met2_mc->Sumw2();
  hztt_met2_mc->SetMarkerColor(kRed);
  hztt_met2_mc->SetMarkerStyle(20);
  TH1D* hztt_mt = new TH1D("hztt_mt","Drell-Yan(#tau#tau) MT",20,0,200);
  hztt_mt->SetDirectory(0);
  hztt_mt->Sumw2();
  TH1D* hztt_mt_of = new TH1D("hztt_mt_of","Drell-Yan(#tau#tau) MT",20,0,200);
  hztt_mt_of->SetDirectory(0);
  hztt_mt_of->Sumw2();
  hztt_mt_of->SetFillColor(kMagenta);
  TH1D* hztt_mt_mc = new TH1D("hztt_mt_mc","Drell-Yan(#tau#tau) MT",20,0,200);
  hztt_mt_mc->SetDirectory(0);
  hztt_mt_mc->Sumw2();
  hztt_mt_mc->SetMarkerColor(kRed);
  hztt_mt_mc->SetMarkerStyle(20);
  {
    SmurfTree tree;
    tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dytt.root");
    tree.InitTree();
    Long64_t nDataEntries = tree.tree_->GetEntries();
    for (Long64_t i = 0; i < nDataEntries; i++){
      tree.tree_->GetEntry(i);
      if ( (tree.cuts_&518)!=518 ) continue;
      // if (tree.type_==0 || tree.type_==3) continue;
      if ( tree.njets_!=njets ) continue;
      if ( tree.mt_<minMT ) continue;
      hztt_mass_mc->Fill(tree.dilep_.mass());
      hztt_met0_mc->Fill(tree.pmet_,tree.scale1fb_*lumi/2);
      hztt_met1_mc->Fill(tree.pTrackMet_,tree.scale1fb_*lumi/2);
      double minmet = std::min(tree.pmet_,tree.pTrackMet_);
      hztt_met2_mc->Fill(minmet,tree.scale1fb_*lumi/2);
      LorentzVector met(tree.met_*cos(tree.metPhi_),tree.met_*sin(tree.metPhi_),0,tree.met_);
      LorentzVector dilep(tree.dilep_);
      // if (minmet>20) hztt_mt_mc->Fill(mt(tree.dilep_.pt(),tree.met_,tree.dPhiDiLepMET_),tree.scale1fb_*lumi/2);
      if (minmet>20) hztt_mt_mc->Fill(mt(dilep,met),tree.scale1fb_*lumi/2);
    }
  }

  TauData taudata;
  SmurfTree tree;
  // tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/data_2l.root");
  // tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dymm.root");
  tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dyee.root");
  tree.InitTree();
  Long64_t nDataEntries = tree.tree_->GetEntries();
  TH1D* hdy = new TH1D("hdy","Drell-Yan(ee/mm) Mass",200,0,200);
  hdy->SetDirectory(0);
  hdy->Sumw2();
  TH1D* hzll_met0 = new TH1D("hzll_met0","Drell-Yan(ee/mm) projected(pfMET)",50,0,100);
  TH1D* hzll_met1 = new TH1D("hzll_met1","Drell-Yan(ee/mm) projected(trkMET)",50,0,100);
  TH1D* hzll_met2 = new TH1D("hzll_met2","Drell-Yan(ee/mm) minProjectedMET",50,0,100);
  TH1D* hzll_met0_of = new TH1D("hzll_met0_of","hzll_met0_of",50,0,100);
  TH1D* hzll_met1_of = new TH1D("hzll_met1_of","hzll_met1_of",50,0,100);
  TH1D* hzll_met2_of = new TH1D("hzll_met2_of","hzll_met2_of",50,0,100);
  hzll_met0->SetDirectory(0);
  hzll_met1->SetDirectory(0);
  hzll_met2->SetDirectory(0);
  hzll_met0_of->SetDirectory(0);
  hzll_met1_of->SetDirectory(0);
  hzll_met2_of->SetDirectory(0);
  hzll_met0->Sumw2();
  hzll_met1->Sumw2();
  hzll_met2->Sumw2();
  hzll_met0_of->Sumw2();
  hzll_met1_of->Sumw2();
  hzll_met2_of->Sumw2();

  TH1D* hztt_mass = new TH1D("hztt_mass","Drell-Yan(#tau#tau) Mass",200,0,200);
  hztt_mass->SetDirectory(0);
  hztt_mass->Sumw2();
  TH1D* hztt_met0 = new TH1D("hztt_met0","Drell-Yan(#tau#tau) projected(pfMET)",50,0,100);
  TH1D* hztt_met1 = new TH1D("hztt_met1","Drell-Yan(#tau#tau) projected(trkMET)",50,0,100);
  TH1D* hztt_met2 = new TH1D("hztt_met2","Drell-Yan(#tau#tau) minProjectedMET",50,0,100);
  TH1D* hztt_met0_of = new TH1D("hztt_met0_of","hztt_met0_of",50,0,100);
  TH1D* hztt_met1_of = new TH1D("hztt_met1_of","hztt_met1_of",50,0,100);
  TH1D* hztt_met2_of = new TH1D("hztt_met2_of","hztt_met2_of",50,0,100);
  hztt_met0->Sumw2();
  hztt_met0->SetDirectory(0);
  hztt_met1->Sumw2();
  hztt_met1->SetDirectory(0);
  hztt_met2->Sumw2();
  hztt_met2->SetDirectory(0);
  hztt_met0_of->Sumw2();
  hztt_met0_of->SetDirectory(0);
  hztt_met1_of->Sumw2();
  hztt_met1_of->SetDirectory(0);
  hztt_met2_of->Sumw2();
  hztt_met2_of->SetDirectory(0);

  for (Long64_t i = 0; i < nDataEntries; i++){
    if (i%prescale != 0 ) continue;
    tree.tree_->GetEntry(i);
    double weight = 1;
    if (useMCWeights) weight = 2*tree.scale1fb_*lumi;
    if ( (tree.cuts_&518)!=518 ) continue;
    if ( tree.njets_!=njets ) continue;
    if ( tree.type_==0 || tree.type_==3 ) hdy->Fill(tree.dilep_.mass(),weight);
    if ( tree.dilep_.mass()>40 &&
	 (!useMassCut || fabs(tree.dilep_.mass()-91)<10) ){
      if ( tree.type_==0 ){ // mm
	hzll_met0->Fill(tree.pmet_,weight/dymm_cor);
	hzll_met1->Fill(tree.pTrackMet_,weight/dymm_cor);
	hzll_met2->Fill(std::min(tree.pmet_,tree.pTrackMet_),weight/dymm_cor);
      }
      if ( tree.type_==3 ){ // ee
	hzll_met0->Fill(tree.pmet_,weight/dyee_cor);
	hzll_met1->Fill(tree.pTrackMet_,weight/dyee_cor);
	hzll_met2->Fill(std::min(tree.pmet_,tree.pTrackMet_),weight/dyee_cor);
      }
      if ( tree.type_==1 || tree.type_==2 ){ // em/me
	hzll_met0_of->Fill(tree.pmet_,weight);
	hzll_met1_of->Fill(tree.pTrackMet_,weight);
	hzll_met2_of->Fill(std::min(tree.pmet_,tree.pTrackMet_),weight);
      }
      LorentzVector t1(tree.lep1_.px(), tree.lep1_.py(), tree.lep1_.pz(), sqrt(tree.lep1_.P()*tree.lep1_.P()+1.777*1.777));
      LorentzVector t2(tree.lep2_.px(), tree.lep2_.py(), tree.lep2_.pz(), sqrt(tree.lep2_.P()*tree.lep2_.P()+1.777*1.777));
      for (unsigned int j=0; j<nToys; ++j){
	LorentzVector met(tree.met_*cos(tree.metPhi_),tree.met_*sin(tree.metPhi_),0,tree.met_);
	LorentzVector trkmet(tree.trackMet_*cos(tree.trackMetPhi_),
			     tree.trackMet_*sin(tree.trackMetPhi_),0,tree.trackMet_);
	taudata.Tauify(t1,spin1);
	LorentzVector lep1 = taudata.GetLeptonFromTau();
	met += taudata.GetMetFromTau(); 
	taudata.Tauify(t2,spin2);
	LorentzVector lep2 = taudata.GetLeptonFromTau();
	met += taudata.GetMetFromTau(); 
	LorentzVector dilep = lep1 + lep2;
	trkmet += tree.lep1_+tree.lep2_;
	trkmet -= lep1+lep2;
	
	if ( mt(dilep,met)<minMT ) continue;

	if ( dilep.mass()>12 && ( (lep1.pt()>10 && lep2.pt()>20) ||
				  (lep1.pt()>20 && lep2.pt()>10) ) &&
	     fabs(lep1.eta())<2.5 && fabs(lep2.eta())<2.5 ){
	  double minmet = std::min(projectedMet(met,lep1,lep2),
				   projectedMet(trkmet,lep1,lep2));
	  if ( tree.type_==0 ){
	    hztt_met0->Fill(projectedMet(met,lep1,lep2),1.0/nToys*bfscale*weight/dymm_cor);
	    hztt_met1->Fill(projectedMet(trkmet,lep1,lep2),1.0/nToys*bfscale*weight/dymm_cor);
	    hztt_met2->Fill(minmet,1.0/nToys*bfscale*weight/dymm_cor);
	    hztt_mass->Fill(dilep.mass(),1.0/nToys*bfscale*weight/dymm_cor);
	    if (minmet>20) hztt_mt->Fill(mt(dilep,met),1.0/nToys*bfscale*weight/dymm_cor);
	  }
	  if ( tree.type_==3 ){
	    hztt_met0->Fill(projectedMet(met,lep1,lep2),1.0/nToys*bfscale*weight/dyee_cor);
	    hztt_met1->Fill(projectedMet(trkmet,lep1,lep2),1.0/nToys*bfscale*weight/dyee_cor);
	    hztt_met2->Fill(minmet,1.0/nToys*bfscale*weight/dyee_cor);
	    hztt_mass->Fill(dilep.mass(),1.0/nToys*bfscale*weight/dyee_cor);
	    if (minmet>20) hztt_mt->Fill(mt(dilep,met),1.0/nToys*bfscale*weight/dyee_cor);
	  }
	  if (tree.type_==1 || tree.type_==2){
	    hztt_met0_of->Fill(projectedMet(met,lep1,lep2),1.0/nToys*bfscale*weight);
	    hztt_met1_of->Fill(projectedMet(trkmet,lep1,lep2),1.0/nToys*bfscale*weight);
	    hztt_met2_of->Fill(minmet,1.0/nToys*bfscale*weight);
	    if (minmet>20) hztt_mt_of->Fill(mt(dilep,met),1.0/nToys*bfscale*weight);
	  }
	}
      }
    }
  }

  std::cout << "data-driven event yield: " << hztt_met1->Integral(0,51) << std::endl;
  std::cout << "MC event yield: "   << hztt_met1_mc->Integral(0,51) << std::endl;
  double ratio = hztt_met1->Integral(0,51)/hztt_met1_mc->Integral(0,51);
  std::cout << "data/mc ratio: "   << ratio << std::endl;

  TCanvas* c0 = new TCanvas("c0","c0",800,400);
  c0->Divide(2,1);
  c0->cd(1);
  gPad->SetLogy(1);
  hdy->Draw();
  c0->cd(2);
  gPad->SetLogy(0);
  hztt_mass->Draw("hist");
  // hztt_mass_mc->Draw("hist e");

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(3,2);
  c1->cd(1);
  gPad->SetLogy(1);
  hzll_met0->Draw();
  hzll_met0_of->SetFillColor(kMagenta);
  hzll_met0_of->Draw("same hist");
  c1->cd(2);
  gPad->SetLogy(1);
  hzll_met1->Draw();
  hzll_met1_of->SetFillColor(kMagenta);
  hzll_met1_of->Draw("same hist");
  c1->cd(3);
  gPad->SetLogy(1);
  hzll_met2->Draw();
  hzll_met2_of->SetFillColor(kMagenta);
  hzll_met2_of->Draw("same hist");
  c1->cd(4);
  gPad->SetLogy(1);
  hztt_met0->Draw("hist");
  hztt_met0_of->SetFillColor(kMagenta);
  hztt_met0_of->Draw("same hist");
  if (scaleMCByArea) hztt_met0_mc->Scale(ratio);
  hztt_met0_mc->Draw("same e");
  std::cout << "pmet>20 (data): " << hztt_met0->Integral(11,51) << "+/-" << integralError(hztt_met0,11,51) << std::endl;
  std::cout << "pmet>20 (MC): "   << hztt_met0_mc->Integral(11,51) << "+/-" << integralError(hztt_met0_mc,11,51) << std::endl;
  c1->cd(5);
  gPad->SetLogy(1);
  hztt_met1->Draw("hist");
  hztt_met1_of->SetFillColor(kMagenta);
  hztt_met1_of->Draw("same hist");
  if (scaleMCByArea) hztt_met1_mc->Scale(ratio);
  hztt_met1_mc->Draw("same e");
  std::cout << "ptrkmet>20 (data): " << hztt_met1->Integral(11,51) << "+/-" << integralError(hztt_met1,11,51) << std::endl;
  std::cout << "ptrkmet>20 (MC): "   << hztt_met1_mc->Integral(11,51) << "+/-" << integralError(hztt_met0_mc,11,51) << std::endl;
  c1->cd(6);
  gPad->SetLogy(1);
  hztt_met2->Draw("hist");
  hztt_met2_of->SetFillColor(kMagenta);
  hztt_met2_of->Draw("same hist");
  if (scaleMCByArea) hztt_met2_mc->Scale(ratio);
  hztt_met2_mc->Draw("same e");
  std::cout << "minPMet>20 (data): " << hztt_met2->Integral(11,51) << "+/-" << integralError(hztt_met2,11,51) << std::endl;
  std::cout << "minPMet>20 (MC): " << hztt_met2_mc->Integral(11,51) << "+/-" << integralError(hztt_met2_mc,11,51) << std::endl;

  new TCanvas("c2","c2",500,500);
  gPad->SetLogy(1);
  hztt_mt->SetLineWidth(2);
  hztt_mt->Draw("hist");
  hztt_mt_of->Draw("hist same");
  // hztt_mt_mc->SetLineColor(kRed);
  if (scaleMCByArea) hztt_mt_mc->Scale(ratio);
  hztt_mt_mc->Draw("hist same e");
  std::cout << "minPMet>20 && MT>80 (data): " << hztt_mt->Integral(9,21) << "+/-" << integralError(hztt_mt,9,21) << std::endl;
  std::cout << "minPMet>20 && MT>80 (MC): " << hztt_mt_mc->Integral(9,21) << "+/-" << integralError(hztt_mt_mc,9,21) << std::endl;
  
}
