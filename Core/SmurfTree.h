#ifndef SmurfTree_H
#define SmurfTree_H
#include <vector>
#include <set>
#include "TFile.h"
#include "TTree.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

//
// Ntiple structure:
//  * plain ntuples without vectors
//  * hypothesis based (not event based). Upto to a user to decide if
//    he allows multiple hypothesis per event or not.
//
// Ntuple content:
//  * event info: run, lumi, event, nVertex, weight (for MC), event type (data, ttbar, wjets etc)
//  * hypothesis type (leading, trailing): ee, em, me, mm
//  * lepton1 (leading)
//  * lepton2 (sub-leading)
//  * jet1 (leading)
//  * jet2 (sub-leading)
//  * jetEventType (0,1,VBF)
//  * nJets
//  * dilepton p4
//  * met
//  * metPhi
//  * sumEt
//  * projectedMet
//  * MT - transverse higgs mass
//  * MT1 - transverse mass of lepton1 and met
//  * MT2 - transverse mass of lepton2 and met
//  * dPhi - delta phi between leptons
//  * dR -   delta R   between leptons   
//  * dPhiLep1Jet1 - delta phi between leading jet and leading lepton
//  * dRLep1Jet1 -   delta R   between leading jet and leading lepton
//  * dPhiLep2Jet1 - delta phi between leading jet and sub-leading lepton
//  * dRLep2Jet1 -   delta R   between leading jet and sub-leading lepton
//  * dPhiLep1Met - delta phi between leading lepton and MET
//  * dPhiLep2Met - delta phi between sub-leading lepton and MET
//  * dPhiDiLepMet - delta phi between di-lepton and MET
//  * dPhiDiLepJet1 - delta phi between di-lepton and leading jet
//  * mc info
//    * lep1 mc_id
//    * lep2 mc_id
//    * jet1 hard-scatter match (flavor)
//    * jet2 hard-scatter match (flavor)
//    * gen met
//
// Lepton variables:
//  * p4
//  * charge
//  * selection type
//  
// Jet variables:
//  * p4
//  * b-tag
// 

struct SmurfTree {
  // float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 
  // bit map
  enum Selection {
    BaseLine    = 0x01,
    Loose       = 0x02,
    Alternative = 0x04
  };
  // first is leading lepton
  enum Type {
    mm, 
    me, 
    em, 
    ee
  };
  enum EventType {
    ZeroJet, 
    OneJet, 
    VBF, 
    AllJets
  };
  enum DataType {
    data,
    qqww,
    ggww,
    hww120,
    hww130,
    hww140,
    hww150,
    hww160,
    hww170,
    hww180,
    hww190,
    hww200,
    hww210,
    hww220,
    hww230,
    hww250,
    hww300,
    hww350,
    hww400,
    hww450,
    hww500,
    hww550,
    hww600,
    vbfhww120,
    vbfhww130,
    vbfhww140,
    vbfhww150,
    vbfhww160,
    vbfhww170,
    vbfhww180,
    vbfhww190,
    vbfhww200,
    vbfhww210,
    vbfhww220,
    vbfhww230,
    vbfhww250,
    vbfhww300,
    vbfhww350,
    vbfhww400,
    vbfhww450,
    vbfhww500,
    vbfhww550,
    vbfhww600,
    ttbar,
    tw,
    dyee,
    dymm,
    dytt,
    wjets,
    wz,
    zz,
    wgamma,
    qcd,
    other
};

  // variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   nvtx_;
  float          scale1fb_;
  float          met_;
  float          metPhi_;
  float          sumet_;
  float          genmet_;
  float          genmetPhi_;
  Type           type_;
  DataType       dstype_;
  LorentzVector  lep1_;
  int            lq1_;
  unsigned int   lid1_;
  LorentzVector  lep2_;
  int            lq2_;
  unsigned int   lid2_;
  LorentzVector  jet1_;
  float          jet1Btag_;
  LorentzVector  jet2_;
  float          jet2Btag_;
  unsigned int   njets_;
  EventType      evtype_;
  LorentzVector  dilep_;
  float          pmet_;
  float          mt_;
  float          mt1_;
  float          mt2_;
  float          dPhi_;
  float          dR_;
  float          dPhiLep1Jet1_;
  float          dRLep1Jet1_;
  float          dPhiLep2Jet1_;
  float          dRLep2Jet1_;
  float          dPhiLep1MET_;
  float          dPhiLep2MET_;
  float          dPhiDiLepMET_;
  float          dPhiDiLepJet1_;
  int            lep1McId_;
  int            lep2McId_;
  int            jet1McId_;
  int            jet2McId_;

  TTree *tree_;
  SmurfTree():
    lepPtr1_(&lep1_),
    lepPtr2_(&lep2_),
    jetPtr1_(&jet1_),
    jetPtr2_(&jet2_),
    dilepPtr_(&dilep_)
  {}

  void InitVariables(){
    event_        = 0;
    run_          = 0;
    lumi_         = 0;
    nvtx_         = 0;
    scale1fb_     = 0;
    met_          = -999.;
    metPhi_       = -999.;
    sumet_        = -999.;
    genmet_       = 0;
    genmetPhi_    = 0;
    type_         = mm;
    dstype_       = data;
    lq1_          = 0;
    lid1_         = 0;
    lq2_          = 0;
    lid2_         = 0;
    jet1Btag_     = -999.;
    jet2Btag_     = -999.;
    njets_        = 0;
    evtype_       = ZeroJet;
    pmet_         = -999.;
    mt_           = -999.;
    mt1_          = -999.;
    mt2_          = -999.;
    dPhi_         = -999.;
    dR_           = -999.;
    dPhiLep1Jet1_ = -999.;
    dRLep1Jet1_   = -999.;
    dPhiLep2Jet1_ = -999.;
    dRLep2Jet1_   = -999.;
    dPhiLep2MET_  = -999.;
    dPhiLep2MET_  = -999.;
    dPhiDiLepMET_ = -999.;
    dPhiDiLepJet1_= -999.;
    lep1McId_   = 0;
    lep2McId_   = 0;
    jet1McId_   = 0;
    jet1McId_   = 0;
  }
  void LoadTree(const char* file){
    TFile* f = TFile::Open(file);
    assert(f);
    tree_ = dynamic_cast<TTree*>(f->Get("tree"));
    assert(tree_);
  }
  void CreateTree(){
    tree_ = new TTree("tree","Smurf ntuples");
    InitVariables();

    //book the branches
    tree_->Branch("event",         &event_,      "event/i");
    tree_->Branch("run",           &run_,        "run/i");
    tree_->Branch("lumi",          &lumi_,       "lumi/i");
    tree_->Branch("nvtx",          &nvtx_,       "nvtx/i");
    tree_->Branch("scale1fb",      &scale1fb_,   "scale1fb/F");
    tree_->Branch("met",           &met_,        "met/F");
    tree_->Branch("metPhi",        &metPhi_,     "metPhi/F");
    tree_->Branch("sumet",         &sumet_,      "sumet/F");
    tree_->Branch("genmet",        &genmet_,     "genmet/F");
    tree_->Branch("genmetPhi",     &genmetPhi_,  "genmetPhi/F");
    tree_->Branch("type",          &type_,       "type/I");
    tree_->Branch("dstype",        &dstype_,     "dstype/I");
    tree_->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lq1",           &lq1_,        "lq1/I");
    tree_->Branch("lid1",          &lid1_,       "lid1/I");
    tree_->Branch("lep2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lq2",           &lq2_,        "lq2/I");
    tree_->Branch("lid2",          &lid2_,       "lid2/I");
    tree_->Branch("jet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1Btag",      &jet1Btag_,   "jet1Btag/F");
    tree_->Branch("jet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet2Btag",      &jet2Btag_,   "jet2Btag/F");
    tree_->Branch("njets",         &njets_,      "njets/I");
    tree_->Branch("evtype",        &evtype_,     "evtype/I");
    tree_->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &dilepPtr_);
    tree_->Branch("pmet",          &pmet_,       "pmet/F");
    tree_->Branch("mt",            &mt_,         "mt/F");
    tree_->Branch("mt1",           &mt1_,        "mt1/F");
    tree_->Branch("mt2",           &mt2_,        "mt2/F");
    tree_->Branch("dPhi",          &dPhi_,       "dPhi/F");
    tree_->Branch("dR",            &dR_,         "dR/F");
    tree_->Branch("dPhiLep1Jet1",  &dPhiLep1Jet1_,   "dPhiLep1Jet1/F");
    tree_->Branch("dRLep1Jet1",    &dRLep1Jet1_,     "dRLep1Jet1/F");
    tree_->Branch("dPhiLep2Jet1",  &dPhiLep2Jet1_,   "dPhiLep2Jet1/F");
    tree_->Branch("dRLep2Jet1",    &dRLep2Jet1_,     "dRLep2Jet1/F");
    tree_->Branch("dPhiDiLepMET",  &dPhiDiLepMET_,   "dPhiDiLepMET/F");
    tree_->Branch("dPhiLep1MET",   &dPhiLep1MET_,    "dPhiLep1MET/F");
    tree_->Branch("dPhiLep1MET",   &dPhiLep2MET_,    "dPhiLep2MET/F");
    tree_->Branch("dPhiDiLepJet1", &dPhiDiLepJet1_,  "dPhiDiLepJet1/F");
    tree_->Branch("lep1McId",      &lep1McId_,   "lep1McId/I");
    tree_->Branch("lep2McId",      &lep2McId_,   "lep2McId/I");
    tree_->Branch("jet1McId",      &jet1McId_,   "jet1McId/I");
    tree_->Branch("jet2McId",      &jet2McId_,   "jet2McId/I");
  }
  private:
  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* dilepPtr_;

}; 
#endif
