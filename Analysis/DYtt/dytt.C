#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->CompileMacro("tauify.cc","k");
  //  gSystem->CompileMacro("../../Core/SmurfTree.h","k");
  gSystem->CompileMacro("dytt.C","k");
  // dytt();
}
#endif 

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
double projectedMet(const SmurfTree::LorentzVector& met, 
		    const SmurfTree::LorentzVector& lep1, 
		    const SmurfTree::LorentzVector& lep2)
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
double mt(SmurfTree::LorentzVector dilep, SmurfTree::LorentzVector met){
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

TH1D* hztt_mass_mc[3];
TH1D* hztt_met0_mc[3];
TH1D* hztt_met1_mc[3];
TH1D* hztt_met2_mc[3];
TH1D* hztt_mt[3];
TH1D* hztt_mt_of[3];
TH1D* hztt_mt_mc[3];
TH1D* hdy[3];
TH1D* hzll_met0[3];
TH1D* hzll_met1[3];
TH1D* hzll_met2[3];
TH1D* hzll_met0_of[3];
TH1D* hzll_met1_of[3];
TH1D* hzll_met2_of[3];
TH1D* hztt_mass[3];
TH1D* hztt_met0[3];
TH1D* hztt_met1[3]; 
TH1D* hztt_met2[3]; 
TH1D* hztt_met0_of[3];
TH1D* hztt_met1_of[3];
TH1D* hztt_met2_of[3];

TH1D* makeTH1D(const char* iname, const char* ititle, Int_t nbinsx, Double_t xlow, Double_t xup, unsigned int i){
  assert(i<3);
  const char* names[3] = {"0-jets","1-jet", "2 and more jets"};
  std::string name(Form("%s_%u",iname,i));
  std::string title(Form("%s (%s)",ititle,names[i]));
  return new TH1D(name.c_str(),title.c_str(),nbinsx,xlow,xup);
}

void init(){
  for(unsigned int i=0; i<3; ++i){
    hztt_mass_mc[i] = makeTH1D("hztt_mass_mc_%u","Drell-Yan(#tau#tau) Mass",200,0,200,i);
    hztt_mass_mc[i]->SetDirectory(0);
    hztt_mass_mc[i]->Sumw2();
    hztt_mass_mc[i]->SetMarkerColor(kRed);
    hztt_mass_mc[i]->SetMarkerStyle(20);
    
    hztt_met0_mc[i] = makeTH1D("hztt_met0_mc","Drell-Yan(#tau#tau) projected(pfMET)",50,0,100,i);
    hztt_met0_mc[i]->SetDirectory(0);
    hztt_met0_mc[i]->Sumw2();
    hztt_met0_mc[i]->SetMarkerColor(kRed);
    hztt_met0_mc[i]->SetMarkerStyle(20);
    hztt_met1_mc[i] = makeTH1D("hztt_met1_mc","Drell-Yan(#tau#tau) projected(trkMET)",50,0,100,i);
    hztt_met1_mc[i]->SetDirectory(0);
    hztt_met1_mc[i]->Sumw2();
    hztt_met1_mc[i]->SetMarkerColor(kRed);
    hztt_met1_mc[i]->SetMarkerStyle(20);
    hztt_met2_mc[i] = makeTH1D("hztt_met2_mc","Drell-Yan(#tau#tau) minProjectedMET",50,0,100,i);
    hztt_met2_mc[i]->SetDirectory(0);
    hztt_met2_mc[i]->Sumw2();
    hztt_met2_mc[i]->SetMarkerColor(kRed);
    hztt_met2_mc[i]->SetMarkerStyle(20);
    hztt_mt[i] = makeTH1D("hztt_mt","Drell-Yan(#tau#tau) MT",20,0,200,i);
    hztt_mt[i]->SetDirectory(0);
    hztt_mt[i]->Sumw2();
    hztt_mt_of[i] = makeTH1D("hztt_mt_of","Drell-Yan(#tau#tau) MT",20,0,200,i);
    hztt_mt_of[i]->SetDirectory(0);
    hztt_mt_of[i]->Sumw2();
    hztt_mt_of[i]->SetFillColor(kMagenta);
    hztt_mt_mc[i] = makeTH1D("hztt_mt_mc","Drell-Yan(#tau#tau) MT",20,0,200,i);
    hztt_mt_mc[i]->SetDirectory(0);
    hztt_mt_mc[i]->Sumw2();
    hztt_mt_mc[i]->SetMarkerColor(kRed);
    hztt_mt_mc[i]->SetMarkerStyle(20);
    
    hdy[i] = makeTH1D("hdy","Drell-Yan(ee/mm) Mass",200,0,200,i);
    hdy[i]->SetDirectory(0);
    hdy[i]->Sumw2();
    hzll_met0[i] = makeTH1D("hzll_met0","Drell-Yan(ee/mm) projected(pfMET)",50,0,100,i);
    hzll_met1[i] = makeTH1D("hzll_met1","Drell-Yan(ee/mm) projected(trkMET)",50,0,100,i);
    hzll_met2[i] = makeTH1D("hzll_met2","Drell-Yan(ee/mm) minProjectedMET",50,0,100,i);
    hzll_met0_of[i] = makeTH1D("hzll_met0_of","hzll_met0_of",50,0,100,i);
    hzll_met1_of[i] = makeTH1D("hzll_met1_of","hzll_met1_of",50,0,100,i);
    hzll_met2_of[i] = makeTH1D("hzll_met2_of","hzll_met2_of",50,0,100,i);
    hzll_met0[i]->SetDirectory(0);
    hzll_met1[i]->SetDirectory(0);
    hzll_met2[i]->SetDirectory(0);
    hzll_met0_of[i]->SetDirectory(0);
    hzll_met1_of[i]->SetDirectory(0);
    hzll_met2_of[i]->SetDirectory(0);
    hzll_met0[i]->Sumw2();
    hzll_met1[i]->Sumw2();
    hzll_met2[i]->Sumw2();
    hzll_met0_of[i]->Sumw2();
    hzll_met1_of[i]->Sumw2();
    hzll_met2_of[i]->Sumw2();

    hztt_mass[i] = makeTH1D("hztt_mass","Drell-Yan(#tau#tau) Mass",200,0,200,i);
    hztt_mass[i]->SetDirectory(0);
    hztt_mass[i]->Sumw2();
    hztt_met0[i] = makeTH1D("hztt_met0","Drell-Yan(#tau#tau) projected(pfMET)",50,0,100,i);
    hztt_met1[i] = makeTH1D("hztt_met1","Drell-Yan(#tau#tau) projected(trkMET)",50,0,100,i);
    // cout << "hztt_met1[i]: " << hztt_met1[i] << endl;
    hztt_met2[i] = makeTH1D("hztt_met2","Drell-Yan(#tau#tau) minProjectedMET",50,0,100,i);
    hztt_met0_of[i] = makeTH1D("hztt_met0_of","hztt_met0_of",50,0,100,i);
    hztt_met1_of[i] = makeTH1D("hztt_met1_of","hztt_met1_of",50,0,100,i);
    hztt_met2_of[i] = makeTH1D("hztt_met2_of","hztt_met2_of",50,0,100,i);
    hztt_met0[i]->Sumw2();
    hztt_met0[i]->SetDirectory(0);
    hztt_met1[i]->Sumw2();
    hztt_met1[i]->SetDirectory(0);
    hztt_met2[i]->Sumw2();
    hztt_met2[i]->SetDirectory(0);
    hztt_met0_of[i]->Sumw2();
    hztt_met0_of[i]->SetDirectory(0);
    hztt_met1_of[i]->Sumw2();
    hztt_met1_of[i]->SetDirectory(0);
    hztt_met2_of[i]->Sumw2();
    hztt_met2_of[i]->SetDirectory(0);
  }
}


#endif

/////////////////////////////////////////////////////////////////////////////////////////////////

namespace dy{ enum EventType {MCee,MCmm,MCall,Data};}

void dytt(const unsigned int nToys, const double lumi,
	  const char* zll, dy::EventType zll_type,
	  const char* mc_dytt)
{
  const unsigned int prescale = 10;
  // Assume that we deal with data, where we have both Zee and Zmm events
  // In this case we get additional branching fraction from tau to e and 
  // tau to mu decays. Out of all 4 final states, we need only 2: em and me.
  // We also should take into account that we have 2 times more data due to 
  // presence of Zee and Zmm
  double bfscale = pow(0.1739+0.1782,2)/2/2*prescale;
  const int spin1 = 0;
  const int spin2 = 0;
  const bool scaleMCByArea = false;
  const double dymm_cor = 1.30; 
  const double dyee_cor = 0.86;
  const double minMT = -1;
  const bool useMassCut = true;
  // const double dymm_cor = 1.35;
  // const double dyee_cor = 0.90;
  init();
  if (mc_dytt){
    printf("Processing DYtt Monte Carlo\n\tPlots are scaled to %0.1f/fb\n",lumi);
    SmurfTree tree;
    tree.LoadTree(mc_dytt);
    tree.InitTree();
    Long64_t nDataEntries = tree.tree_->GetEntries();
    int i_permille_old = 0;
    for (Long64_t i = 0; i < nDataEntries; i++){
      int i_permille = (int)floor(1000 * i / double(nDataEntries));
      if (i_permille != i_permille_old) {
	// xterm magic from L. Vacavant and A. Cerri
	printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
	       "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	fflush(stdout);
	i_permille_old = i_permille;
      }
      tree.tree_->GetEntry(i);
      if ( (tree.cuts_&518)!=518 ) continue;
      // if (tree.type_==0 || tree.type_==3) continue;
      unsigned int njets = tree.njets_;
      if (njets>2) njets = 2;
      if ( tree.mt_<minMT ) continue;
      hztt_mass_mc[njets]->Fill(tree.dilep_.mass());
      hztt_met0_mc[njets]->Fill(tree.pmet_,tree.scale1fb_*lumi/2);
      hztt_met1_mc[njets]->Fill(tree.pTrackMet_,tree.scale1fb_*lumi/2);
      double minmet = std::min(tree.pmet_,tree.pTrackMet_);
      hztt_met2_mc[njets]->Fill(minmet,tree.scale1fb_*lumi/2);
      SmurfTree::LorentzVector met(tree.met_*cos(tree.metPhi_),tree.met_*sin(tree.metPhi_),0,tree.met_);
      SmurfTree::LorentzVector dilep(tree.dilep_);
      // if (minmet>20) hztt_mt_mc[njets]->Fill(mt(tree.dilep_.pt(),tree.met_,tree.dPhiDiLepMET_),tree.scale1fb_*lumi/2);
      if (minmet>20) hztt_mt_mc[njets]->Fill(mt(dilep,met),tree.scale1fb_*lumi/2);
    }
    printf("Done processing DYtt Monte Carlo\n");
  }

  TFile* fOut = TFile::Open("dytt_from_z.root","RECREATE");
  assert(fOut);
  TauData taudata;
  SmurfTree tree;
  tree.LoadTree(zll);
  tree.InitTree();
  TTree* outTree = tree.tree_->CloneTree(0);
  // Variables that need to be updated:
  //   + lep1
  //   + lep2
  //   + dilep
  //   + met
  //   + metPhi
  //   + pmet
  //   + trackMet
  //   + trackMetPhi
  //   + pTrackMet
  //   + mt
  //   + scale1fb
  //   - cuts
  //   + dPhi
  //   + dR
  //   - dPhiLep1Jet1
  //   - dRLep1Jet1
  //   - dPhiLep2Jet1
  //   - dRLep2Jet1
  //   - dPhiLep1MET
  //   - dPhiLep2MET
  //   + dPhiDiLepMET
  //   + dPhiDiLepJet1
  
  outTree->SetDirectory(fOut);
  SmurfTree::LorentzVector lep1;
  SmurfTree::LorentzVector* lepPtr1(&lep1);
  outTree->SetBranchAddress("lep1", &lepPtr1);
  SmurfTree::LorentzVector lep2;
  SmurfTree::LorentzVector* lepPtr2(&lep2);
  outTree->SetBranchAddress("lep2", &lepPtr2);
  SmurfTree::LorentzVector dilep;
  SmurfTree::LorentzVector* dilepPtr(&dilep);
  outTree->SetBranchAddress("dilep", &dilepPtr);
  float met;           outTree->SetBranchAddress("met", &met);
  float pmet;          outTree->SetBranchAddress("pmet", &pmet);
  float metPhi;        outTree->SetBranchAddress("metPhi", &metPhi);
  float trackMet;      outTree->SetBranchAddress("trackMet", &trackMet);
  float trackMetPhi;   outTree->SetBranchAddress("trackMetPhi", &trackMetPhi);
  float pTrackMet;     outTree->SetBranchAddress("pTrackMet", &pTrackMet);
  float mtVar;         outTree->SetBranchAddress("mt", &mtVar);
  Int_t type;          outTree->SetBranchAddress("type", &type);
  float scale1fb;      outTree->SetBranchAddress("scale1fb", &scale1fb);
  float dR;            outTree->SetBranchAddress("dR", &dR);
  float dPhi;          outTree->SetBranchAddress("dPhi", &dPhi);
  float dPhiDiLepMET;  outTree->SetBranchAddress("dPhiDiLepMET", &dPhiDiLepMET);
  float dPhiDiLepJet1; outTree->SetBranchAddress("dPhiDiLepJet1", &dPhiDiLepJet1);

  
  Long64_t nDataEntries = tree.tree_->GetEntries();
  printf("Processing Zll data\n\tNumber of pseudo tau decays per single Zll event: %u\n\tPrescale: %u\n",nToys,prescale);
  printf("\tDYmm correction to get DYtt(emu) rate: %0.2f\n",dymm_cor);
  printf("\tDYee correction to get DYtt(emu) rate: %0.2f\n",dyee_cor);
  bool useMCWeights = false;
  if (zll_type == dy::MCee ){
    printf("\tInput is Zee Monte Carlo.\n");
    bfscale *= 2;
    useMCWeights = true;
  }
  if (zll_type == dy::MCmm ){
    printf("\tInput is Zmm Monte Carlo.\n");
    bfscale *= 2;
    useMCWeights = true;
  }
  if (zll_type == dy::MCall ){
    printf("\tInput is Zll Monte Carlo.\n");
    useMCWeights = true;
  }

  int i_permille_old = 0;
  for (Long64_t i = 0; i < nDataEntries; i++){
    if (i%prescale != 0 ) continue;
      int i_permille = (int)floor(1000 * i / double(nDataEntries));
      if (i_permille != i_permille_old) {
	// xterm magic from L. Vacavant and A. Cerri
	printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
	       "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	fflush(stdout);
	i_permille_old = i_permille;
      }
    tree.tree_->GetEntry(i);
    double weight = 1;
    if (useMCWeights) weight = tree.scale1fb_*lumi;
    if ( (tree.cuts_&518)!=518 ) continue;
    unsigned int njets = tree.njets_;
    if (njets>2) njets = 2;
    if ( tree.type_==0 || tree.type_==3 ) hdy[njets]->Fill(tree.dilep_.mass(),weight);
    if ( tree.dilep_.mass()>40 &&
	 (!useMassCut || fabs(tree.dilep_.mass()-91)<10) ){
      // copy event
      if ( tree.type_==0 ){ // mm
	hzll_met0[njets]->Fill(tree.pmet_,weight/dymm_cor);
	hzll_met1[njets]->Fill(tree.pTrackMet_,weight/dymm_cor);
	hzll_met2[njets]->Fill(std::min(tree.pmet_,tree.pTrackMet_),weight/dymm_cor);
      }
      if ( tree.type_==3 ){ // ee
	hzll_met0[njets]->Fill(tree.pmet_,weight/dyee_cor);
	hzll_met1[njets]->Fill(tree.pTrackMet_,weight/dyee_cor);
	hzll_met2[njets]->Fill(std::min(tree.pmet_,tree.pTrackMet_),weight/dyee_cor);
      }
      if ( tree.type_==1 || tree.type_==2 ){ // em/me
	hzll_met0_of[njets]->Fill(tree.pmet_,weight);
	hzll_met1_of[njets]->Fill(tree.pTrackMet_,weight);
	hzll_met2_of[njets]->Fill(std::min(tree.pmet_,tree.pTrackMet_),weight);
      }
      SmurfTree::LorentzVector t1(tree.lep1_.px(), tree.lep1_.py(), tree.lep1_.pz(), sqrt(tree.lep1_.P()*tree.lep1_.P()+1.777*1.777));
      SmurfTree::LorentzVector t2(tree.lep2_.px(), tree.lep2_.py(), tree.lep2_.pz(), sqrt(tree.lep2_.P()*tree.lep2_.P()+1.777*1.777));
      for (unsigned int j=0; j<nToys; ++j){
	SmurfTree::LorentzVector met_p4(tree.met_*cos(tree.metPhi_),tree.met_*sin(tree.metPhi_),0,tree.met_);
	SmurfTree::LorentzVector trkmet_p4(tree.trackMet_*cos(tree.trackMetPhi_),
					   tree.trackMet_*sin(tree.trackMetPhi_),0,tree.trackMet_);
	taudata.Tauify(t1,spin1);
	lep1 = taudata.GetLeptonFromTau();
	met_p4 += taudata.GetMetFromTau(); 
	taudata.Tauify(t2,spin2);
	lep2 = taudata.GetLeptonFromTau();
	met_p4 += taudata.GetMetFromTau(); 
	if (lep1.pt()<lep2.pt()) std::swap(lep1,lep2);
	dilep = lep1 + lep2;
	trkmet_p4 += tree.lep1_+tree.lep2_;
	trkmet_p4 -= lep1+lep2;
	
	if ( mt(dilep,met_p4)<minMT ) continue;

	if ( dilep.mass()>12 && ( (lep1.pt()>10 && lep2.pt()>20) ||
				  (lep1.pt()>20 && lep2.pt()>10) ) &&
	     fabs(lep1.eta())<2.5 && fabs(lep2.eta())<2.5 ){
	  double minmet = std::min(projectedMet(met_p4,lep1,lep2),
				   projectedMet(trkmet_p4,lep1,lep2));
	  met = met_p4.pt();
	  pmet = projectedMet(met_p4,lep1,lep2);
	  metPhi = met_p4.phi();
	  trackMet = trkmet_p4.pt();
	  pTrackMet = projectedMet(trkmet_p4,lep1,lep2);
	  trackMetPhi = trkmet_p4.phi();
	  mtVar = mt(dilep,met_p4);
	  type = SmurfTree::me;
	  dPhi = acos(cos(lep1.phi()-lep2.phi()));
	  dR = sqrt(dPhi*dPhi+pow(lep1.eta()-lep2.eta(),2));
	  dPhiDiLepJet1 = acos(cos(dilep.phi()-tree.jet1_.phi()));
	  dPhiDiLepMET = acos(cos(dilep.phi()-metPhi));

	  if ( tree.type_==0 ){
	    scale1fb = 1.0/nToys*bfscale*weight/dymm_cor/lumi;
	    hztt_met0[njets]->Fill(projectedMet(met_p4,lep1,lep2),1.0/nToys*bfscale*weight/dymm_cor);
	    hztt_met1[njets]->Fill(projectedMet(trkmet_p4,lep1,lep2),1.0/nToys*bfscale*weight/dymm_cor);
	    hztt_met2[njets]->Fill(minmet,1.0/nToys*bfscale*weight/dymm_cor);
	    hztt_mass[njets]->Fill(dilep.mass(),1.0/nToys*bfscale*weight/dymm_cor);
	    if (minmet>20) hztt_mt[njets]->Fill(mt(dilep,met_p4),1.0/nToys*bfscale*weight/dymm_cor);
	    outTree->Fill();
	  }
	  if ( tree.type_==3 ){
	    scale1fb = 1.0/nToys*bfscale*weight/dyee_cor/lumi;
	    hztt_met0[njets]->Fill(projectedMet(met_p4,lep1,lep2),1.0/nToys*bfscale*weight/dyee_cor);
	    hztt_met1[njets]->Fill(projectedMet(trkmet_p4,lep1,lep2),1.0/nToys*bfscale*weight/dyee_cor);
	    hztt_met2[njets]->Fill(minmet,1.0/nToys*bfscale*weight/dyee_cor);
	    hztt_mass[njets]->Fill(dilep.mass(),1.0/nToys*bfscale*weight/dyee_cor);
	    if (minmet>20) hztt_mt[njets]->Fill(mt(dilep,met_p4),1.0/nToys*bfscale*weight/dyee_cor);
	    outTree->Fill();
	  }
	  if (tree.type_==1 || tree.type_==2){
	    hztt_met0_of[njets]->Fill(projectedMet(met_p4,lep1,lep2),1.0/nToys*bfscale*weight);
	    hztt_met1_of[njets]->Fill(projectedMet(trkmet_p4,lep1,lep2),1.0/nToys*bfscale*weight);
	    hztt_met2_of[njets]->Fill(minmet,1.0/nToys*bfscale*weight);
	    if (minmet>20) hztt_mt_of[njets]->Fill(mt(dilep,met_p4),1.0/nToys*bfscale*weight);
	  }
	}
      }
    }
  }
  fOut->cd();
  outTree->Write();
  fOut->Close();

  printf("Done processing Zll data\n");

  const char* name[3] = {"0-jets","1-jet", "2 and more jets"};
  for (unsigned int i=0; i<3; ++i){
    std::cout << name[i] << std::endl;
    std::cout << "\tdata-driven event yield: " << hztt_met1[i]->Integral(0,51) << std::endl;
    std::cout << "\tMC event yield: "   << hztt_met1_mc[i]->Integral(0,51) << std::endl;
    double ratio = hztt_met1[i]->Integral(0,51)/hztt_met1_mc[i]->Integral(0,51);
    std::cout << "\tdata/mc ratio: "   << ratio << std::endl;

    TCanvas* c0 = new TCanvas(Form("c0_%u",i),"Canvas",800,400);
    c0->Divide(2,1);
    c0->cd(1);
    gPad->SetLogy(1);
    hdy[i]->Draw();
    c0->cd(2);
    gPad->SetLogy(0);
    hztt_mass[i]->Draw("hist");
    // hztt_mass_mc[i]->Draw("hist e");

    TCanvas* c1 = new TCanvas(Form("c1_%u",i),"Canvas",1200,800);
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogy(1);
    hzll_met0[i]->Draw();
    hzll_met0_of[i]->SetFillColor(kMagenta);
    hzll_met0_of[i]->Draw("same hist");
    c1->cd(2);
    gPad->SetLogy(1);
    hzll_met1[i]->Draw();
    hzll_met1_of[i]->SetFillColor(kMagenta);
    hzll_met1_of[i]->Draw("same hist");
    c1->cd(3);
    gPad->SetLogy(1);
    hzll_met2[i]->Draw();
    hzll_met2_of[i]->SetFillColor(kMagenta);
    hzll_met2_of[i]->Draw("same hist");
    c1->cd(4);
    gPad->SetLogy(1);
    hztt_met0[i]->Draw("hist");
    hztt_met0_of[i]->SetFillColor(kMagenta);
    hztt_met0_of[i]->Draw("same hist");
    if (scaleMCByArea) hztt_met0_mc[i]->Scale(ratio);
    hztt_met0_mc[i]->Draw("same e");
    std::cout << "\tpmet>20 (data): " << hztt_met0[i]->Integral(11,51) << "+/-" << integralError(hztt_met0[i],11,51) << std::endl;
    std::cout << "\tpmet>20 (MC): "   << hztt_met0_mc[i]->Integral(11,51) << "+/-" << integralError(hztt_met0_mc[i],11,51) << std::endl;
    c1->cd(5);
    gPad->SetLogy(1);
    hztt_met1[i]->Draw("hist");
    hztt_met1_of[i]->SetFillColor(kMagenta);
    hztt_met1_of[i]->Draw("same hist");
    if (scaleMCByArea) hztt_met1_mc[i]->Scale(ratio);
    hztt_met1_mc[i]->Draw("same e");
    std::cout << "\tptrkmet>20 (data): " << hztt_met1[i]->Integral(11,51) << "+/-" << integralError(hztt_met1[i],11,51) << std::endl;
    std::cout << "\tptrkmet>20 (MC): "   << hztt_met1_mc[i]->Integral(11,51) << "+/-" << integralError(hztt_met0_mc[i],11,51) << std::endl;
    c1->cd(6);
    gPad->SetLogy(1);
    hztt_met2[i]->Draw("hist");
    hztt_met2_of[i]->SetFillColor(kMagenta);
    hztt_met2_of[i]->Draw("same hist");
    if (scaleMCByArea) hztt_met2_mc[i]->Scale(ratio);
    hztt_met2_mc[i]->Draw("same e");
    std::cout << "\tminPMet>20 (data): " << hztt_met2[i]->Integral(11,51) << "+/-" << integralError(hztt_met2[i],11,51) << std::endl;
    std::cout << "\tminPMet>20 (MC): " << hztt_met2_mc[i]->Integral(11,51) << "+/-" << integralError(hztt_met2_mc[i],11,51) << std::endl;

    new TCanvas(Form("c2_%u",i),"c2",500,500);
    gPad->SetLogy(1);
    hztt_mt[i]->SetLineWidth(2);
    hztt_mt[i]->Draw("hist");
    hztt_mt_of[i]->Draw("hist same");
    // hztt_mt_mc[i]->SetLineColor(kRed);
    if (scaleMCByArea) hztt_mt_mc[i]->Scale(ratio);
    hztt_mt_mc[i]->Draw("hist same e");
    std::cout << "\tminPMet>20 && MT>80 (data): " << hztt_mt[i]->Integral(9,21) << "+/-" << integralError(hztt_mt[i],9,21) << std::endl;
    std::cout << "\tminPMet>20 && MT>80 (MC): " << hztt_mt_mc[i]->Integral(9,21) << "+/-" << integralError(hztt_mt_mc[i],9,21) << std::endl;
  }
}
