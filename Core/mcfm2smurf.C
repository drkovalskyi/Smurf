#define mcfm2smurf_cxx
#include "mcfm2smurf.h"
#include "SmurfTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include "TRandom3.h"
#include "TMath.h"

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

bool passedWWSelection(const SmurfTree& tree){
  if (tree.lep1_.pt()<20||abs(tree.lep1_.eta())>2.4) return false;
  if (tree.lep2_.pt()<10||abs(tree.lep2_.eta())>2.4) return false;
  if (tree.dilep_.mass()<12) return false;
  if ( tree.type_==SmurfTree::em || tree.type_==SmurfTree::me ){
    if (tree.pmet_<20) return false;
  } else {
    if (tree.pmet_<40) return false;
    if (tree.jet1_.pt()>15 && tree.dPhiDiLepJet1_>M_PI/180*165) return false;
  }
  return true;
}

double nearestDeltaPhi(const SmurfTree& tree, double phi)
{
  double dPhi1 = fabs(tree.lep1_.phi() - phi);
  dPhi1 = std::min(2*TMath::Pi() - dPhi1, dPhi1);
  double dPhi2 = fabs(tree.lep2_.phi() - phi);
  dPhi2 = std::min(2*TMath::Pi() - dPhi2, dPhi2);
  return std::min(dPhi1, dPhi2);
}

double projectedMet(const SmurfTree& tree, double met, double phi)
{
  double dPhi = nearestDeltaPhi(tree,phi);
  if (dPhi < TMath::Pi()/2) return met*sin(dPhi);
  return met;
}

void mcfm2smurf::MakeSmurfNtuple(const char* filename, int mcfm_process_id)
{
  SmurfTree tree;
  tree.CreateTree();
  tree.tree_->SetDirectory(0);
  TRandom3 gen;
  
  if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);
      tree.InitVariables();

      switch (mcfm_process_id){
      case 61:
	{
	  tree.lep1_ = SmurfTree::LorentzVector(px4,py4,pz4,E_4);
	  tree.lep2_ = SmurfTree::LorentzVector(px5,py5,pz5,E_5);
	  SmurfTree::LorentzVector nu1(px3,py3,pz3,E_3);
	  SmurfTree::LorentzVector nu2(px6,py6,pz6,E_6);
	  tree.jet1_ = -(tree.lep1_+tree.lep2_+nu1+nu2);
	  if (tree.jet1_.pt()>30) 
	    tree.njets_ = 1;
	  else
	    tree.njets_ = 0;
	  if (tree.lep1_.pt() < tree.lep2_.pt()) std::swap(tree.lep1_,tree.lep2_);
	  tree.met_ = (nu1+nu2).pt();
	  tree.metPhi_ = (nu1+nu2).phi();
	}
	break;
      default:
	std::cout << "Unsupported process id: " << mcfm_process_id << std::endl;
	return;
      }
	
      tree.dilep_ = tree.lep1_ + tree.lep2_;
      tree.dPhi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(tree.lep1_,tree.lep2_));
      tree.dR_ = ROOT::Math::VectorUtil::DeltaR(tree.lep1_,tree.lep2_);
      tree.weight_ = wt_ALL;
      // tree.type_ = SmurfTree::Type(gen.Uniform(4));
      tree.type_ = SmurfTree::Type(0);
      
      tree.trackMet_ = tree.met_;
      tree.trackMetPhi_ = tree.metPhi_;
      tree.pmet_ = projectedMet(tree,tree.met_,tree.metPhi_);
      tree.pTrackMet_ = tree.pmet_;
      
      tree.mt_ = mt(tree.dilep_.pt(),tree.met_,tree.dPhiDiLepMET_);
      tree.mt1_ = mt(tree.lep1_.pt(),tree.met_,tree.dPhiLep1MET_);
      tree.mt2_ = mt(tree.lep2_.pt(),tree.met_,tree.dPhiLep2MET_);

      tree.dPhiDiLepMET_  = acos(cos(tree.dilep_.phi()-tree.metPhi_));
      tree.dPhiLep1MET_   = acos(cos(tree.lep1_.phi()-tree.metPhi_));
      tree.dPhiLep2MET_   = acos(cos(tree.lep2_.phi()-tree.metPhi_));
      tree.dPhiDiLepJet1_ = acos(cos(tree.dilep_.phi()-tree.jet1_.phi()));
      
      if (passedWWSelection(tree)) tree.tree_->Fill();
   }
   TFile* fSmurf = TFile::Open(filename,"RECREATE");
   assert(fSmurf);
   tree.tree_->Write();
   fSmurf->Close();
   
}
