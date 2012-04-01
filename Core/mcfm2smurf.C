#define mcfm2smurf_cxx
#include "mcfm2smurf.h"
#include "SmurfTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include "TRandom3.h"
#include "TMath.h"

typedef std::pair<SmurfTree::LorentzVector,int> mpair;

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

bool accept(const SmurfTree::LorentzVector& p4){
  return fabs(p4.eta())<2.4;
}

bool passedWWSelection(const SmurfTree& tree){
  if (tree.lep1_.pt()<20 || !accept(tree.lep1_)) return false;
  if (tree.lep2_.pt()<10 || !accept(tree.lep2_)) return false;
  if (tree.dilep_.mass()<12) return false;
  if (tree.pmet_<20) return false;
//   if ( tree.type_==SmurfTree::ee || tree.type_==SmurfTree::mm ){
//     if (tree.pmet_<40) return false;
//     if (tree.jet1_.pt()>15 && tree.dPhiDiLepJet1_>M_PI/180*165) return false;
//   }
  // isolation
  if (ROOT::Math::VectorUtil::DeltaR(tree.lep1_,tree.lep2_)<0.3) return false;
  if (tree.lep3_.pt()>0 && ROOT::Math::VectorUtil::DeltaR(tree.lep1_,tree.lep3_)<0.3) return false;
  if (tree.lep3_.pt()>0 && ROOT::Math::VectorUtil::DeltaR(tree.lep2_,tree.lep3_)<0.3) return false;
  
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

bool sortVectors(const mpair& v1, const mpair& v2)
{
  return v1.first.pt() > v2.first.pt();
}

bool mcfm2smurf::FillProcess(SmurfTree& tree, int mcfm_process_id){
  switch (mcfm_process_id){
  case 61:
    {
      tree.lep1_ = SmurfTree::LorentzVector(px4,py4,pz4,E_4);
      tree.lep2_ = SmurfTree::LorentzVector(px5,py5,pz5,E_5);
      SmurfTree::LorentzVector nu1(px3,py3,pz3,E_3);
      SmurfTree::LorentzVector nu2(px6,py6,pz6,E_6);
      tree.quadlep_ = tree.lep1_+tree.lep2_+nu1+nu2; 
      tree.jet1_ = -tree.quadlep_;
      if (tree.jet1_.pt()>30) 
	tree.njets_ = 1;
      else
	tree.njets_ = 0;
      if (tree.lep1_.pt() < tree.lep2_.pt()) std::swap(tree.lep1_,tree.lep2_);
      tree.met_ = (nu1+nu2).pt();
      tree.metPhi_ = (nu1+nu2).phi();
      tree.dstype_ = SmurfTree::qqww;
    }
    break;
  case 71:
    {
      SmurfTree::LorentzVector l1(px4,py4,pz4,E_4);
      SmurfTree::LorentzVector l2(px5,py5,pz5,E_5);
      SmurfTree::LorentzVector l3(px6,py6,pz6,E_6);
      SmurfTree::LorentzVector nu1(px3,py3,pz3,E_3);
      SmurfTree::LorentzVector met = nu1;
      std::vector<mpair> leps;
      if ( accept(l1) )
	leps.push_back(mpair(l1,24));
      else
	met += l1;
      if ( accept(l2) )
	leps.push_back(mpair(l2,23));
      else
	met += l2;
      if ( accept(l3) )
	leps.push_back(mpair(l3,23));
      else
	met += l3;
      if ( leps.size()<2 ) return false;
      std::sort(leps.begin(),leps.end(),sortVectors);
      tree.lep1_ = leps.at(0).first;
      tree.lep1MotherMcId_ = leps.at(0).second;
      tree.lep2_ = leps.at(1).first;
      tree.lep2MotherMcId_ = leps.at(1).second;
      if ( leps.size()==3 ){
	tree.lep3_ = leps.at(2).first;
	tree.lep3MotherMcId_ = leps.at(2).second;
      }

      tree.quadlep_ = l1+l2+l3+nu1;
      tree.jet1_ = -tree.quadlep_;
      if (tree.jet1_.pt()>30) 
	tree.njets_ = 1;
      else
	tree.njets_ = 0;
      tree.met_ = met.pt();
      tree.metPhi_ = met.phi();
      tree.dstype_ = SmurfTree::wz;
    }
    break;
  default:
    std::cout << "Unsupported process id: " << mcfm_process_id << std::endl;
    exit(1);
  }
  return true;
}

void mcfm2smurf::MakeSmurfNtuple(const char* filename, int mcfm_process_id)
{
  TFile* fSmurf = TFile::Open(filename,"RECREATE");
  assert(fSmurf);
  SmurfTree tree;
  tree.CreateTree();
  TRandom3 gen;
  
  if (fChain == 0) {
    std::cout << "Tree is invalid" << std::endl;
    return;
  }
  Long64_t nentries = fChain->GetEntriesFast();
  int i_permille_old = 0;
  std::cout << "Number of entries to process: " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries; ++jentry) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);
      
      int i_permille = (int)floor(1000 * jentry / float(nentries));
      if (i_permille != i_permille_old) {
	printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
	       "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	fflush(stdout);
	i_permille_old = i_permille;
      }

      tree.InitVariables();
      
      if (!FillProcess(tree,mcfm_process_id)) continue;
      tree.dilep_ = tree.lep1_ + tree.lep2_;
      tree.dPhi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(tree.lep1_,tree.lep2_));
      tree.dR_ = ROOT::Math::VectorUtil::DeltaR(tree.lep1_,tree.lep2_);
      tree.scale1fb_ = wt_ALL;
      // tree.type_ = SmurfTree::Type(gen.Uniform(4));
      tree.type_ = SmurfTree::Type(0);
      tree.cuts_ = SmurfTree::FullSelection;
      
      tree.trackMet_ = tree.met_;
      tree.trackMetPhi_ = tree.metPhi_;
      tree.pmet_ = projectedMet(tree,tree.met_,tree.metPhi_);
      tree.pTrackMet_ = tree.pmet_;
      
      tree.dPhiDiLepMET_  = acos(cos(tree.dilep_.phi()-tree.metPhi_));
      tree.dPhiLep1MET_   = acos(cos(tree.lep1_.phi()-tree.metPhi_));
      tree.dPhiLep2MET_   = acos(cos(tree.lep2_.phi()-tree.metPhi_));
      tree.dPhiDiLepJet1_ = acos(cos(tree.dilep_.phi()-tree.jet1_.phi()));
      
      tree.mt_ = mt(tree.dilep_.pt(),tree.met_,tree.dPhiDiLepMET_);
      tree.mt1_ = mt(tree.lep1_.pt(),tree.met_,tree.dPhiLep1MET_);
      tree.mt2_ = mt(tree.lep2_.pt(),tree.met_,tree.dPhiLep2MET_);

      if (passedWWSelection(tree)) tree.tree_->Fill();
   }
   tree.tree_->Write();
   fSmurf->Close();
   
}
