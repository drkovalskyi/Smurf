#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TGraphErrors.h>           // Graph class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O

#include "Smurf/Core/RecoilCorrector.hh"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h"
#include "Smurf/Analysis/HWWlvlv/HWWCuts.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#endif

static const unsigned int basic_selection = SmurfTree::BaseLine | 
                                            SmurfTree::ChargeMatch | 
                                            SmurfTree::Lep1FullSelection |
                                            SmurfTree::Lep2FullSelection |
                                            SmurfTree::ExtraLeptonVeto |
                                            SmurfTree::TopVeto;

static const Bool_t verbose = kFALSE;

//=== FUNCTION DECLARATIONS ======================================================================================

// compute projected MET
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi);

// compute systematic uncertainty
Double_t computeSyst(const TH1D *hout, const TH1D *hin, Int_t binUsed, Double_t rErrorMax, Bool_t isDebug = false);

//=== MAIN MACRO =================================================================================================
// Options:
//
// --- ZWindowSubtractionMethod == 0 : Opposite Flavor Background Subtraction
// --- ZWindowSubtractionMethod == 1 : Data Corrected MC Background Subtraction
// 
// --- MassCutLow  : lower  mass cut with respect to z pole mass
// --- MassCutHigh : higher mass cut with respect to z pole mass
//

void ComputeDYBkgScaleFactor(Int_t period = 0, Bool_t useRecoilModel = kFALSE, Int_t ZWindowSubtractionMethod = 0, 
                             Double_t MassCutLow = 7.5, Double_t MassCutHigh = 7.5)
{

  //*******************************************************
  // Settings 
  //*******************************************************
  bool WWXSSel = false;
  double ptLepMin = 10.0;
  if(WWXSSel == true) ptLepMin = 20.;

  bool forBDTAna = false;
  bool useRFromData[3] = {0, 0, 1};

  Double_t lumi = 1;
  TString filesPath   = "dummy";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  enum { kMuMu, kEleEle, kEleMu, kMuEle };
  
  Int_t nmet = 1;
  if(useRecoilModel) nmet = 100;

  if     (period == 0){ // Full2011-Fall11-V9-3500ipb
    lumi = 3.553;minRun =      0;maxRun = 999999;
    filesPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_3500ipb";
  }
  else if(period == -1){ // Full2011-Fall11-V9 NoJetId
    lumi = 3.553;minRun =      0;maxRun = 999999;
    filesPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_NoJetId";
  }
  else if(period == 1){ // Full2011-Fall11-V9-5000ipb-2nd
    lumi = 5.098;minRun =      0;maxRun = 999999;
    //filesPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_5000ipb";
    filesPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets";
  }
  else if(period == 2){ // Full2011-Fall11-V9-5000ipb-newId
    lumi = 5.098;minRun =      0;maxRun = 999999;
    filesPath  = "/data/smurf/data/Run2012_Summer12_SmurfV10_52X/mitf-alljets";
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  

  infilenamev.push_back(Form("%s/data.root",filesPath.Data()));

  if (ZWindowSubtractionMethod == 0) {
   infilenamev.push_back(Form("%s/dyll.root",filesPath.Data()));
   infilenamev.push_back(Form("%s/wz.root",filesPath.Data()));
   infilenamev.push_back(Form("%s/zz.root",filesPath.Data()));
  } else if (ZWindowSubtractionMethod == 1) {
    infilenamev.push_back(Form("%s/backgroundA.root",filesPath.Data()));
  }

  //*******************************************************
  //Met Corrections
  //*******************************************************
  vector <RecoilCorrector*> metCorrections[2];
  metCorrections[0].push_back(new RecoilCorrector("/data/smurf/sixie/data/auxiliar/recoilFits/recoilfit_datamm_0jet.root"));
  metCorrections[0].push_back(new RecoilCorrector("/data/smurf/sixie/data/auxiliar/recoilFits/recoilfit_datamm_1jet.root"));
  metCorrections[0].push_back(new RecoilCorrector("/data/smurf/sixie/data/auxiliar/recoilFits/recoilfit_datamm_2jet.root"));
  metCorrections[1].push_back(new RecoilCorrector("/data/smurf/sixie/data/auxiliar/recoilFits/recoilfit_dataee_0jet.root"));
  metCorrections[1].push_back(new RecoilCorrector("/data/smurf/sixie/data/auxiliar/recoilFits/recoilfit_dataee_1jet.root"));
  metCorrections[1].push_back(new RecoilCorrector("/data/smurf/sixie/data/auxiliar/recoilFits/recoilfit_dataee_2jet.root"));
  
  //*******************************************************
  //mH dependent Cuts
  //*******************************************************
  const Double_t mZ = 91.1876;
  
  const Int_t nbins = 4;  
  Float_t bins[nbins+1];

  const Int_t nmass = 22;
  const Double_t mH[nmass]     = {0,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};  
        Bool_t useDYMVA[nmass] = {0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
  
  //*******************************************************
  //Yields and  histograms
  //*******************************************************
  vector<Double_t> nin_kee_data, nin_kmm_data;
  
  vector<vector<TH1D*> > hNin_ree_mc,   hNout_ree_mc,   hNin_rmm_mc,   hNout_rmm_mc;
  vector<vector<TH1D*> > hNin_ree_da,   hNout_ree_da,   hNin_rmm_da,   hNout_rmm_da;

  vector<vector<Double_t> > nin_ee_dy, nout_ee_dy, nin_ee_vz, nout_ee_vz, nin_ee_data;  
  vector<vector<Double_t> > nin_mm_dy, nout_mm_dy, nin_mm_vz, nout_mm_vz, nin_mm_data;
  vector<vector<Double_t> > varin_ee_dy, varout_ee_dy, varin_ee_vz, varout_ee_vz;
  vector<vector<Double_t> > varin_mm_dy, varout_mm_dy, varin_mm_vz, varout_mm_vz;
  vector<vector<Double_t> > nin_em_data, nin_me_data;
  vector<vector<Double_t> > nin_ee_bkg, nin_mm_bkg;
  vector<vector<Double_t> > varin_ee_bkg, varin_mm_bkg;
  vector<vector<Double_t> > nin_ee_ww, nin_mm_ww;
  vector<vector<Double_t> > varin_ee_ww, varin_mm_ww;
  vector<vector<Double_t> > nin_ee_wjets, nin_mm_wjets;
  vector<vector<Double_t> > varin_ee_wjets, varin_mm_wjets;
  vector<vector<Double_t> > nin_ee_top, nin_mm_top;
  vector<vector<Double_t> > varin_ee_top, varin_mm_top;
  vector<vector<Double_t> > nin_ee_dytt, nin_mm_dytt;
  vector<vector<Double_t> > varin_ee_dytt, varin_mm_dytt;

  
  for(UInt_t jetIndex = 0; jetIndex < 3; ++jetIndex) {
    Double_t tmp_nin_kee_data=0, tmp_nin_kmm_data=0;
    vector<TH1D*> tmp_hNin_ree_mc,   tmp_hNout_ree_mc,   tmp_hNin_rmm_mc,   tmp_hNout_rmm_mc;
    vector<TH1D*> tmp_hNin_ree_da,   tmp_hNout_ree_da,   tmp_hNin_rmm_da,   tmp_hNout_rmm_da;
    
    vector<Double_t> tmp_nin_ee_dy, tmp_nout_ee_dy, tmp_nin_ee_vz, tmp_nout_ee_vz, tmp_nin_ee_data;  
    vector<Double_t> tmp_nin_mm_dy, tmp_nout_mm_dy, tmp_nin_mm_vz, tmp_nout_mm_vz, tmp_nin_mm_data;
    vector<Double_t> tmp_varin_ee_dy, tmp_varout_ee_dy, tmp_varin_ee_vz, tmp_varout_ee_vz;
    vector<Double_t> tmp_varin_mm_dy, tmp_varout_mm_dy, tmp_varin_mm_vz, tmp_varout_mm_vz;
    vector<Double_t> tmp_nin_em_data, tmp_nin_me_data;
    vector<Double_t> tmp_nin_ee_bkg, tmp_nin_mm_bkg;
    vector<Double_t> tmp_varin_ee_bkg, tmp_varin_mm_bkg;
    vector<Double_t> tmp_nin_ee_ww, tmp_nin_mm_ww;
    vector<Double_t> tmp_varin_ee_ww, tmp_varin_mm_ww;
    vector<Double_t> tmp_nin_ee_wjets, tmp_nin_mm_wjets;
    vector<Double_t> tmp_varin_ee_wjets, tmp_varin_mm_wjets;
    vector<Double_t> tmp_nin_ee_top, tmp_nin_mm_top;
    vector<Double_t> tmp_varin_ee_top, tmp_varin_mm_top;
    vector<Double_t> tmp_nin_ee_dytt, tmp_nin_mm_dytt;
    vector<Double_t> tmp_varin_ee_dytt, tmp_varin_mm_dytt;

    char hname[50];
    for(Int_t imass=0; imass<nmass; imass++) {
      if     (useDYMVA[imass] == kTRUE){
        if     (jetIndex == 0){
	  bins[0] = -0.90; bins[1] = -0.85; bins[2] = -0.60; bins[3] =  0.60; bins[4] = 1.0; 
	}
	else if(jetIndex == 1){
	  bins[0] = -0.90; bins[1] = -0.85; bins[2] = -0.60; bins[3] =  0.30; bins[4] = 1.0; 
	}
	else if(jetIndex == 2){
	  bins[0] = 20; bins[1] = 25; bins[2] = 30; bins[3] =  45; bins[4] = 50; 
	}
      } else {
	bins[0] = 20; bins[1] = 25; bins[2] = 30; bins[3] =  45; bins[4] = 50; 
      }

      sprintf(hname,"hNin_%iJet_ree_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_ree_mc.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_ree_mc[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_ree_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_ree_mc.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_ree_mc[imass]->Sumw2();
      sprintf(hname,"hNin_%iJet_rmm_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_rmm_mc.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_rmm_mc[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_rmm_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_rmm_mc.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_rmm_mc[imass]->Sumw2();
      
      sprintf(hname,"hNin_%iJet_ree_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_ree_da.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_ree_da[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_ree_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_ree_da.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_ree_da[imass]->Sumw2();
      sprintf(hname,"hNin_%iJet_rmm_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_rmm_da.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_rmm_da[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_rmm_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_rmm_da.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_rmm_da[imass]->Sumw2();
      
      tmp_nin_ee_dy.push_back(0), tmp_nout_ee_dy.push_back(0), tmp_nin_ee_vz.push_back(0), tmp_nout_ee_vz.push_back(0), tmp_nin_ee_data.push_back(0);
      tmp_nin_mm_dy.push_back(0), tmp_nout_mm_dy.push_back(0), tmp_nin_mm_vz.push_back(0), tmp_nout_mm_vz.push_back(0), tmp_nin_mm_data.push_back(0);
      tmp_varin_ee_dy.push_back(0), tmp_varout_ee_dy.push_back(0), tmp_varin_ee_vz.push_back(0), tmp_varout_ee_vz.push_back(0);    
      tmp_varin_mm_dy.push_back(0), tmp_varout_mm_dy.push_back(0), tmp_varin_mm_vz.push_back(0), tmp_varout_mm_vz.push_back(0);
      tmp_nin_em_data.push_back(0), tmp_nin_me_data.push_back(0);
      tmp_nin_ee_bkg.push_back(0), tmp_nin_mm_bkg.push_back(0);
      tmp_varin_ee_bkg.push_back(0), tmp_varin_mm_bkg.push_back(0);
      tmp_nin_ee_ww.push_back(0), tmp_nin_mm_ww.push_back(0);
      tmp_varin_ee_ww.push_back(0), tmp_varin_mm_ww.push_back(0);
      tmp_nin_ee_wjets.push_back(0), tmp_nin_mm_wjets.push_back(0);
      tmp_varin_ee_wjets.push_back(0), tmp_varin_mm_wjets.push_back(0);
      tmp_nin_ee_top.push_back(0), tmp_nin_mm_top.push_back(0);
      tmp_varin_ee_top.push_back(0), tmp_varin_mm_top.push_back(0);
      tmp_nin_ee_dytt.push_back(0), tmp_nin_mm_dytt.push_back(0);
      tmp_varin_ee_dytt.push_back(0), tmp_varin_mm_dytt.push_back(0);
     }

    nin_kee_data.push_back(tmp_nin_kee_data);
    nin_kmm_data.push_back(tmp_nin_kmm_data);

    hNin_ree_mc.push_back(tmp_hNin_ree_mc); 
    hNout_ree_mc.push_back(tmp_hNout_ree_mc); 
    hNin_rmm_mc.push_back(tmp_hNin_rmm_mc); 
    hNout_rmm_mc.push_back(tmp_hNout_rmm_mc);

    hNin_ree_da.push_back(tmp_hNin_ree_da); 
    hNout_ree_da.push_back(tmp_hNout_ree_da); 
    hNin_rmm_da.push_back(tmp_hNin_rmm_da); 
    hNout_rmm_da.push_back(tmp_hNout_rmm_da);

    nin_ee_dy.push_back(tmp_nin_ee_dy);
    nout_ee_dy.push_back(tmp_nout_ee_dy);
    nin_ee_vz.push_back(tmp_nin_ee_vz);
    nout_ee_vz.push_back(tmp_nout_ee_vz);
    nin_ee_data.push_back(tmp_nin_ee_data);
    nin_mm_dy.push_back(tmp_nin_mm_dy);
    nout_mm_dy.push_back(tmp_nout_mm_dy);
    nin_mm_vz.push_back(tmp_nin_mm_vz);
    nout_mm_vz.push_back(tmp_nout_mm_vz);
    nin_mm_data.push_back(tmp_nin_mm_data);
    varin_ee_dy.push_back(tmp_varin_ee_dy);
    varout_ee_dy.push_back(tmp_varout_ee_dy);
    varin_ee_vz.push_back(tmp_varin_ee_vz);
    varout_ee_vz.push_back(tmp_varout_ee_vz);
    varin_mm_dy.push_back(tmp_varin_mm_dy);
    varout_mm_dy.push_back(tmp_varout_mm_dy);
    varin_mm_vz.push_back(tmp_varin_mm_vz);
    varout_mm_vz.push_back(tmp_varout_mm_vz);
    nin_em_data.push_back(tmp_nin_em_data);
    nin_me_data.push_back(tmp_nin_me_data);

    nin_ee_bkg.push_back(tmp_nin_ee_bkg);
    nin_mm_bkg.push_back(tmp_nin_mm_bkg);
    varin_ee_bkg.push_back(tmp_varin_ee_bkg);
    varin_mm_bkg.push_back(tmp_varin_mm_bkg);
    nin_ee_ww.push_back(tmp_nin_ee_ww);
    nin_mm_ww.push_back(tmp_nin_mm_ww);
    varin_ee_ww.push_back(tmp_varin_ee_ww);
    varin_mm_ww.push_back(tmp_varin_mm_ww);
    nin_ee_wjets.push_back(tmp_nin_ee_wjets);
    nin_mm_wjets.push_back(tmp_nin_mm_wjets);
    varin_ee_wjets.push_back(tmp_varin_ee_wjets);
    varin_mm_wjets.push_back(tmp_varin_mm_wjets);
    nin_ee_top.push_back(tmp_nin_ee_top);
    nin_mm_top.push_back(tmp_nin_mm_top);
    varin_ee_top.push_back(tmp_varin_ee_top);
    varin_mm_top.push_back(tmp_varin_mm_top);
    nin_ee_dytt.push_back(tmp_nin_ee_dytt);
    nin_mm_dytt.push_back(tmp_nin_mm_dytt);
    varin_ee_dytt.push_back(tmp_varin_ee_dytt);
    varin_mm_dytt.push_back(tmp_varin_mm_dytt);

  }

  //*******************************************************
  //Systematic Error on VZ Normalization
  //*******************************************************
  const Double_t vzNormSystematic = 0.10; //10% systematic on the VZ cross section inside Z window

  //*******************************************************
  //Event Loop
  //*******************************************************
  SmurfTree tree;
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    cout << "Processing " << infilenamev[ifile] << "..." << endl;
    tree.LoadTree(infilenamev[ifile]);
    tree.InitTree(0);
    Float_t dymva = -100.0;
    tree.tree_->SetBranchAddress(Form("dymva"), &dymva);
    
    for(UInt_t ientry = 0; ientry <tree.tree_->GetEntries(); ientry++){
      tree.tree_->GetEntry(ientry);
        
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      if(tree.dstype_ == SmurfTree::data && tree.run_ <  minRun) continue;
      if(tree.dstype_ == SmurfTree::data && tree.run_ >  maxRun) continue;

      //apply trigger requirement
      if(tree.dstype_ == SmurfTree::data)
        if ( (tree.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger ) continue; 
        
      if(!((tree.cuts_ & SmurfTree::BaseLine) == SmurfTree::BaseLine)            ) continue;
      if(!((tree.cuts_ & SmurfTree::ChargeMatch) == SmurfTree::ChargeMatch)      ) continue;
      if( (tree.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
      if( !((tree.cuts_ & SmurfTree::TopVeto) == SmurfTree::TopVeto)             ) continue; // cut on btagging
         
      Int_t ijet = tree.njets_;
      if(ijet >= 2){
        if(tree.jet3_.Pt() <= 30)					                       ijet = 2;
        else if(tree.jet3_.Pt() > 30 && (
    	  (tree.jet1_.Eta()-tree.jet3_.Eta() > 0 && tree.jet2_.Eta()-tree.jet3_.Eta() < 0) ||
    	  (tree.jet2_.Eta()-tree.jet3_.Eta() > 0 && tree.jet1_.Eta()-tree.jet3_.Eta() < 0)))   ijet = 3;
        else							                               ijet = 2;
        if(tree.njets_ < 2 || tree.njets_ > 3)                                                 ijet = 3;
      }
      if(ijet>2) continue;

      if(tree.lep1_.Pt() < 20) continue;
      if(tree.lep2_.Pt() < ptLepMin) continue;
 
      //*******************************************************
      //Classify background type
      //*******************************************************          
      Int_t ProcessType = -99;
      if (tree.dstype_==SmurfTree::data ) {
        if ( ifile == 0)                               ProcessType = -1; //2l data
      }
      if ((tree.dstype_==SmurfTree::wz || 
           tree.dstype_==SmurfTree::zz)
          && (tree.lep1MotherMcId_==23 && 
              tree.lep2MotherMcId_==23))               ProcessType = 1; //WZ,ZZ with both lepton from same Z
      if(tree.dstype_ == SmurfTree::tw              )  ProcessType = 2; //top
      if(tree.dstype_ == SmurfTree::ttbar           )  ProcessType = 2; //top
      if(tree.dstype_ == SmurfTree::qqww            )  ProcessType = 0; //WW
      if(tree.dstype_ == SmurfTree::ggww            )  ProcessType = 0; //WW
      if(tree.dstype_ == SmurfTree::dytt            )  ProcessType = 5; //DYtt
      if(tree.dstype_ == SmurfTree::dyttDataDriven  )  ProcessType = 5; //DYtt
      if(tree.dstype_ == SmurfTree::qcd             )  ProcessType = 5; //DYtt
      if(tree.dstype_ == SmurfTree::wgamma          )  ProcessType = 4; //wgamma
            
      Double_t weight = 1;
      //*******************************************************
      //Handle Fake Bkg
      //*******************************************************          
      int nFake = 0;
      if(((tree.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (tree.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((tree.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (tree.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((tree.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (tree.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((tree.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (tree.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      
      if(nFake > 1){
        weight = 0.0;
      }
      else if(nFake == 1){
        if(tree.dstype_ == SmurfTree::data){
          weight = tree.sfWeightFR_;
          ProcessType = 4;
        }
        // MC real lepton fake contamination
        else if(TMath::Abs(tree.lep1McId_)*TMath::Abs(tree.lep2McId_) > 0 
                || tree.dstype_ == SmurfTree::wgamma
          ) {
	  // negative weight is already included
          weight = lumi*tree.scale1fb_*tree.sfWeightPU_*tree.sfWeightEff_*tree.sfWeightTrig_*tree.sfWeightFR_;
          ProcessType = 4;
	  if(weight > 0) printf("warning, weight > 0! (%f)\n",weight);
        }
        else {
          weight = 0.0;
        }
      }
      else if(tree.dstype_ != SmurfTree::data) {
        //normal MC
        //require both leptons pass full selection
        if (!(((tree.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
           && ((tree.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;

        weight *= lumi*tree.scale1fb_*tree.sfWeightPU_*tree.sfWeightEff_*tree.sfWeightTrig_;
      } else {
        //tight+tight data
        if (!(((tree.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
           && ((tree.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
      }
      //*******************************************************
      //Data Driven Normalization Scale Factors
      //*******************************************************          
      //Top Scale Factors
      if(ProcessType == 2) {
    	if(tree.njets_ == 0) { weight=weight*TopBkgScaleFactor(0); }
    	if(tree.njets_ == 1) { weight=weight*TopBkgScaleFactor(1); }
    	if(tree.njets_ >= 2) { weight=weight*TopBkgScaleFactor(2); }
      }
      //WW Scale Factors
      if(ProcessType == 0) {
    	if(tree.njets_ == 0) { weight=weight*WWBkgScaleFactorCutBased(0,0); }
        if(tree.njets_ >= 1) { weight=weight*WWBkgScaleFactorCutBased(0,1); }
      }

      Double_t pfmet     = tree.met_;
      Double_t pfmetphi  = tree.metPhi_;
      Double_t trkmet    = tree.trackMet_;
      Double_t trkmetphi = tree.trackMetPhi_;
    
      Int_t finalState=-1;
      if(tree.type_==SmurfTree::mm) finalState=kMuMu;
      if(tree.type_==SmurfTree::ee) finalState=kEleEle;
      if(tree.type_==SmurfTree::me) finalState=kMuEle;
      if(tree.type_==SmurfTree::em) finalState=kEleMu;
      assert(finalState > -1);
    
      if(tree.dstype_==SmurfTree::data && mZ - tree.dilep_.M() < MassCutLow && tree.dilep_.M() - mZ < MassCutHigh) {
        if(finalState==kEleEle) nin_kee_data[ijet]++;
        if(finalState==kMuMu)	nin_kmm_data[ijet]++;
      }
      
      //For Z->mm/ee MC
      if((tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm) && weight >0) {

        //sample Met model many times (nmet times)
        for(Int_t imet=0; imet<nmet; imet++) {
          
          //Apply Met Model Prediction
          if (useRecoilModel) {
            if(finalState==kMuMu || finalState==kEleEle) {
              metCorrections[finalState][ijet] ->Correct(pfmet,pfmetphi,trkmet,trkmetphi,tree.dilep_.Pt(),tree.dilep_.Phi(),tree.dilep_.Pt(),tree.dilep_.Phi());
            }
          }

	  Double_t minpfmet  = TMath::Min(projectedMET(pfmet,pfmetphi,tree.lep1_.Phi()),  projectedMET(pfmet,pfmetphi,tree.lep2_.Phi()));
          Double_t mintrkmet = TMath::Min(projectedMET(trkmet,trkmetphi,tree.lep1_.Phi()),projectedMET(trkmet,trkmetphi,tree.lep2_.Phi()));
          Double_t minmet    = TMath::Min(minpfmet,mintrkmet);
          if(minmet>=49) minmet=49;
	  Double_t mt = sqrt( 2.0 * (tree.dilep_.Pt()) * pfmet * (1.0-cos( acos(cos(tree.dilep_.Phi() - pfmetphi)) ) ));

          //loop over Higgs masses
	  for(Int_t imass=0; imass<nmass; imass++) {
	    int whichIMass = imass;
	    double theCutMassHigh = cutMassHigh(mH[imass]);
            bool passMtCut = mt > cutMTLow(mH[imass]) && mt < cutMTHigh(mH[imass]);

            if(forBDTAna == true && imass != 0){
	      whichIMass = 0;
	      theCutMassHigh = DileptonMassPreselectionCut(mH[imass]);
	      passMtCut = mt > 80.0 && mt < mH[imass];
	    }

	    bool passKinCuts = tree.lep1_.Pt() > cutPtMaxLow(mH[whichIMass])   &&
	                       tree.lep2_.Pt() > cutPtMinLow(mH[whichIMass],0) &&
			       tree.dPhi_ < cutDeltaphiHigh(mH[whichIMass])*TMath::Pi()/180.;
	    bool passWBFSel = true;
            if(tree.njets_ >= 2 && imass != 0){
	      theCutMassHigh = DileptonMassPreselectionCut(mH[imass]);
	      passMtCut = mt > 30.0 && mt < mH[imass];  
              int centrality = 0;
	      if(((tree.jet1_.Eta()-tree.lep1_.Eta() > 0 && tree.jet2_.Eta()-tree.lep1_.Eta() < 0) ||
        	  (tree.jet2_.Eta()-tree.lep1_.Eta() > 0 && tree.jet1_.Eta()-tree.lep1_.Eta() < 0)) &&
        	 ((tree.jet1_.Eta()-tree.lep2_.Eta() > 0 && tree.jet2_.Eta()-tree.lep2_.Eta() < 0) ||
        	  (tree.jet2_.Eta()-tree.lep2_.Eta() > 0 && tree.jet1_.Eta()-tree.lep2_.Eta() < 0))) centrality = 1; 
	      passWBFSel = (tree.jet1_+tree.jet2_).M() > 450. && TMath::Abs(tree.jet1_.Eta()-tree.jet2_.Eta()) > 3.5 && 
                            centrality == 1;
	      passKinCuts = true;
	    }
            if(!passWBFSel) continue;

	    if(!passKinCuts) continue;
	    if(tree.dilep_.Pt() <= 45.0) continue;
	    if(minmet <= 20.0) continue;

            bool dPhiDiLepJetCut = kTRUE;
            if(useDYMVA[imass] == kFALSE){
              if(tree.njets_ <= 1) dPhiDiLepJetCut = tree.jet1_.Pt() <= 15. || tree.dPhiDiLepJet1_*180.0/TMath::Pi() < 165.;
              else                 dPhiDiLepJetCut = DeltaPhi((tree.jet1_+tree.jet2_).Phi(),tree.dilep_.Phi())*180.0/TMath::Pi() < 165.;
            }
            if(tree.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((tree.jet1_+tree.jet2_).Phi(),tree.dilep_.Phi())*180.0/TMath::Pi() < 165.;
	    if(dPhiDiLepJetCut == kFALSE) continue;

	    Double_t varMet = minmet;
	    if     (ijet == 2)                varMet = tree.met_;
	    else if(useDYMVA[imass] == kTRUE) varMet = dymva;

	    if(mZ - tree.dilep_.M() < MassCutLow && tree.dilep_.M() - mZ < MassCutHigh) {
	      if(finalState==kEleEle && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) {
                hNin_ree_mc[ijet][imass]->Fill(varMet,weight/(Double_t)nmet);                
              }
	      if(finalState==kMuMu   && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) {
                hNin_rmm_mc[ijet][imass]->Fill(varMet,weight/(Double_t)nmet);              	 
              }
            } else if(fabs(tree.dilep_.M() - mZ) >= 15 && tree.dilep_.M() < theCutMassHigh) {
	      if(finalState==kEleEle && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) {
                hNout_ree_mc[ijet][imass]->Fill(varMet,weight/(Double_t)nmet);	
              }
	      if(finalState==kMuMu   && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) {
                hNout_rmm_mc[ijet][imass]->Fill(varMet,weight/(Double_t)nmet);
              }
            }
	
 	    if(!passMtCut) continue;

            bool passMET = minmet > 20.;
            if(useDYMVA[imass] == kTRUE){
              if     (tree.njets_ == 0) passMET = passMET && dymva >  0.60;
              else if(tree.njets_ == 1) passMET = passMET && dymva >  0.30;
              else                      passMET = passMET && tree.met_ > 45.0;
            } else {
              if     (tree.njets_ == 0) passMET = passMET && minmet > 45.;
              else if(tree.njets_ == 1) passMET = passMET && minmet > 45.;
              else                      passMET = passMET && tree.met_ > 45.0;
            }

            if(passMET == kFALSE) continue;

	    if(mZ - tree.dilep_.M() < MassCutLow && tree.dilep_.M() - mZ < MassCutHigh) {
	      if(finalState==kEleEle && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) { 
	        nin_ee_dy[ijet][imass]+=weight/(Double_t)nmet; 
	        varin_ee_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet; 
	      }
	
	      if(finalState==kMuMu   && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) { 
	        nin_mm_dy[ijet][imass]+=weight/(Double_t)nmet; 
	        varin_mm_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
	      }

             } else if(fabs(tree.dilep_.M() - mZ) >= 15 && tree.dilep_.M() < theCutMassHigh) {

	      if(finalState==kEleEle && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) { 
	        nout_ee_dy[ijet][imass]+=weight/(Double_t)nmet; 
	        varout_ee_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
	      }
	  
	      if(finalState==kMuMu   && (tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm)) { 
	        nout_mm_dy[ijet][imass]+=weight/(Double_t)nmet;
	        varout_mm_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
	      }
            }
          }                
	}
      } else {
      //For Other Bkg MC

        Double_t minpfmet  = TMath::Min(projectedMET(pfmet,pfmetphi,tree.lep1_.Phi()),  projectedMET(pfmet,pfmetphi,tree.lep2_.Phi()));
        Double_t mintrkmet = TMath::Min(projectedMET(trkmet,trkmetphi,tree.lep1_.Phi()),projectedMET(trkmet,trkmetphi,tree.lep2_.Phi()));
        Double_t minmet    = TMath::Min(minpfmet,mintrkmet);
        if(minmet>50) minmet=49;
	
	Double_t mt = sqrt( 2.0 * (tree.dilep_.Pt()) * pfmet * (1.0-cos(acos(cos(tree.dilep_.Phi()-pfmetphi)))) );
	
	for(Int_t imass=0; imass<nmass; imass++) {
	  int whichIMass = imass;
	  double theCutMassHigh = cutMassHigh(mH[imass]);
          bool passMtCut = mt > cutMTLow(mH[imass]) && mt < cutMTHigh(mH[imass]);

          if(forBDTAna == true && imass != 0){
	    whichIMass = 0;
	    theCutMassHigh = DileptonMassPreselectionCut(mH[imass]);
	    passMtCut = mt > 80.0 && mt < mH[imass];
	  }

	  bool passKinCuts = tree.lep1_.Pt() > cutPtMaxLow(mH[whichIMass])   &&
	  		     tree.lep2_.Pt() > cutPtMinLow(mH[whichIMass],0) &&
	        	     tree.dPhi_ < cutDeltaphiHigh(mH[whichIMass])*TMath::Pi()/180.;
	  bool passWBFSel = true;
          if(tree.njets_ >= 2 && imass != 0){
	    theCutMassHigh = DileptonMassPreselectionCut(mH[imass]);
	    passMtCut = mt > 30.0 && mt < mH[imass];  
            int centrality = 0;
	    if(((tree.jet1_.Eta()-tree.lep1_.Eta() > 0 && tree.jet2_.Eta()-tree.lep1_.Eta() < 0) ||
        	(tree.jet2_.Eta()-tree.lep1_.Eta() > 0 && tree.jet1_.Eta()-tree.lep1_.Eta() < 0)) &&
               ((tree.jet1_.Eta()-tree.lep2_.Eta() > 0 && tree.jet2_.Eta()-tree.lep2_.Eta() < 0) ||
        	(tree.jet2_.Eta()-tree.lep2_.Eta() > 0 && tree.jet1_.Eta()-tree.lep2_.Eta() < 0))) centrality = 1; 
	    passWBFSel = (tree.jet1_+tree.jet2_).M() > 450. && TMath::Abs(tree.jet1_.Eta()-tree.jet2_.Eta()) > 3.5 && 
                          centrality == 1;
	    passKinCuts = true;
	  }
          if(!passWBFSel) continue;

	  if(!passKinCuts) continue;
          if(tree.dilep_.Pt() <= 45.0) continue;
	  if(minmet <= 20.0) continue;

          bool dPhiDiLepJetCut = kTRUE;
          if(useDYMVA[imass] == kFALSE){
            if(tree.njets_ <= 1) dPhiDiLepJetCut = tree.jet1_.Pt() <= 15. || tree.dPhiDiLepJet1_*180.0/TMath::Pi() < 165.;
            else		 dPhiDiLepJetCut = DeltaPhi((tree.jet1_+tree.jet2_).Phi(),tree.dilep_.Phi())*180.0/TMath::Pi() < 165.;
          }
          if(tree.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((tree.jet1_+tree.jet2_).Phi(),tree.dilep_.Phi())*180.0/TMath::Pi() < 165.;
	  if(dPhiDiLepJetCut == kFALSE) continue;

	  Double_t varMet = minmet;
	  if	 (ijet == 2)		    varMet = tree.met_;
	  else if(useDYMVA[imass] == kTRUE) varMet = dymva;
          double theDataWeight = 0.0;
          if     (ProcessType == -1 && (finalState==kEleEle||finalState==kMuMu)) theDataWeight =  1.0;
	  else if(ProcessType == -1                                            ) theDataWeight = -1.0;
	  else if(ProcessType ==  1 && weight > 0                              ) theDataWeight = -1.0 * weight;
	  // In reality ee == ee+mm for data
	  if(theDataWeight != 0.0){
	    if(mZ - tree.dilep_.M() < MassCutLow && tree.dilep_.M() - mZ < MassCutHigh) {
              hNin_ree_da[ijet][imass]->Fill(varMet,theDataWeight/(Double_t)nmet);                
            } else if(fabs(tree.dilep_.M() - mZ) >= 15 && tree.dilep_.M() < theCutMassHigh) {
              hNout_ree_da[ijet][imass]->Fill(varMet,theDataWeight/(Double_t)nmet);	
            }
          }
	  
          if(!passMtCut) continue;

          bool passMET = minmet > 20.;
          if(useDYMVA[imass] == kTRUE){
            if     (tree.njets_ == 0) passMET = passMET && dymva >  0.60;
            else if(tree.njets_ == 1) passMET = passMET && dymva >  0.30;
            else		      passMET = passMET && tree.met_ > 45.0;
          } else {
            if     (tree.njets_ == 0) passMET = passMET && minmet > 45.;
            else if(tree.njets_ == 1) passMET = passMET && minmet > 45.;
            else		      passMET = passMET && tree.met_ > 45.0;
          }
          if(passMET == kFALSE) continue;

          //*********************************************************************
          //Count Yields
          //*********************************************************************
          //In Z peak region
	  if(mZ - tree.dilep_.M() < MassCutLow && tree.dilep_.M() - mZ < MassCutHigh) {
            //if(tree.dstype_==SmurfTree::data && ProcessType == -1 && imass==0 && ijet==2) printf("SSS %d %d\n",tree.event_,tree.type_);
	  
            //first one 
            if(finalState==kEleEle) {
              if(ProcessType == -1) { 
                nin_ee_data[ijet][imass]++; 
              }
              if(ProcessType == 1)  {
                nin_ee_vz[ijet][imass]+=weight;
                varin_ee_vz[ijet][imass]+=weight*weight;
              }
              if(ProcessType == 0) {
                nin_ee_bkg[ijet][imass]+=weight;
                varin_ee_bkg[ijet][imass]+=weight*weight;
                nin_ee_ww[ijet][imass]+=weight;
                varin_ee_ww[ijet][imass]+=weight*weight;
              }
              if(ProcessType == 4) {
                nin_ee_bkg[ijet][imass]+=weight;
                varin_ee_bkg[ijet][imass]+=weight*weight;
                nin_ee_wjets[ijet][imass]+=weight;
                varin_ee_wjets[ijet][imass]+=weight*weight;
              }
              if(ProcessType == 2) {
                nin_ee_bkg[ijet][imass]+=weight;
                varin_ee_bkg[ijet][imass]+=weight*weight;
                nin_ee_top[ijet][imass]+=weight;
                varin_ee_top[ijet][imass]+=weight*weight;
              }
              if(ProcessType == 5) {
                nin_ee_bkg[ijet][imass]+=weight;
                varin_ee_bkg[ijet][imass]+=weight*weight;
                nin_ee_dytt[ijet][imass]+=weight;
                varin_ee_dytt[ijet][imass]+=weight*weight;
              }
            }
	
            if(finalState==kMuMu) {
              if(ProcessType == -1 ) { 
                nin_mm_data[ijet][imass]++;
              }	  
              if(ProcessType == 1 ) {
                nin_mm_vz[ijet][imass]+=weight;
                varin_mm_vz[ijet][imass]+=weight*weight;
              }              
              if(ProcessType == 0 ) {
                nin_mm_bkg[ijet][imass]+=weight;
                varin_mm_bkg[ijet][imass]+=weight*weight;
                nin_mm_ww[ijet][imass]+=weight;
                varin_mm_ww[ijet][imass]+=weight*weight;
              }              
              if(ProcessType == 4 ) {
                nin_mm_bkg[ijet][imass]+=weight;
                varin_mm_bkg[ijet][imass]+=weight*weight;
                nin_mm_wjets[ijet][imass]+=weight;
                varin_mm_wjets[ijet][imass]+=weight*weight;
              }              
              if(ProcessType == 2 ) {
                nin_mm_bkg[ijet][imass]+=weight;
                varin_mm_bkg[ijet][imass]+=weight*weight;
                nin_mm_top[ijet][imass]+=weight;
                varin_mm_top[ijet][imass]+=weight*weight;
              }              
              if(ProcessType == 5 ) {
                nin_mm_bkg[ijet][imass]+=weight;
                varin_mm_bkg[ijet][imass]+=weight*weight;
                nin_mm_dytt[ijet][imass]+=weight;
                varin_mm_dytt[ijet][imass]+=weight*weight;
              }  

            }

	    if(finalState==kMuEle && tree.dstype_==SmurfTree::data && ProcessType == -1) {
              nin_me_data[ijet][imass]++;           
            }
	    if(finalState==kEleMu && tree.dstype_==SmurfTree::data && ProcessType == -1) {
              nin_em_data[ijet][imass]++; 
            }

          } 
          // Out of Z peak region
          else if(fabs(tree.dilep_.M() - mZ) >= 15 && tree.dilep_.M() < theCutMassHigh) {

	    if(finalState==kEleEle) {
              if (ProcessType == 1) {
                nout_ee_vz[ijet][imass]+=weight;
                varout_ee_vz[ijet][imass]+=weight*weight;
              }
	    }
	  
	    if(finalState==kMuMu) {
              if ( ProcessType == 1) {
                nout_mm_vz[ijet][imass]+=weight;
                varout_mm_vz[ijet][imass]+=weight*weight;
              }
	    }
          }

        }
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  char DYEstimateTableName[200];
  sprintf(DYEstimateTableName,"DYEstimateTable.txt");
  if(forBDTAna == true) sprintf(DYEstimateTableName,"DYEstimateTableBDT.txt");
  ofstream fout(DYEstimateTableName);

  vector<vector<Double_t> > DYBkgScaleFactorHiggsSelection;
  vector<vector<Double_t> > DYBkgScaleFactorHiggsSelectionErr;
  vector<Double_t> DYBkgScaleFactorWWPreselection;
  vector<Double_t> DYBkgScaleFactorWWPreselectionErr;

  for(UInt_t jetIndex = 0; jetIndex < 3; ++jetIndex) {
    
    vector<Double_t> tmpDYBkgScaleFactorHiggsSelection;
    vector<Double_t> tmpDYBkgScaleFactorHiggsSelectionErr;
    Double_t tmpDYBkgScaleFactorWWPreselection;
    Double_t tmpDYBkgScaleFactorWWPreselectionErr;

    fout << "************************************************************************\n";
    fout << jetIndex << "-Jet Bin DY Bkg Scale Factor Computation\n";
    fout << "************************************************************************\n";

    Double_t k    = sqrt(nin_kee_data[jetIndex]/nin_kmm_data[jetIndex]);
    Double_t kerr = 0.5*k*sqrt(1.0/nin_kee_data[jetIndex] + 1.0/nin_kmm_data[jetIndex]);
    char buffer[200];
  
    fout << "Electron to muon efficiency ratio is " << k << " +/- " << kerr << endl;

    fout << endl;
    fout << jetIndex << "-jet bin summary:" << endl;
    fout << setw(4) << "sel" << "   ";
    fout << setw(25) << "R_out/in       ";
    fout << setw(15) << "mm/me/em/ee";
    fout << setw(20) << "N_in (OF,VZ sub)";
    fout << setw(20) << "N_in (data/MC)";
    fout << setw(20) << "   N_out data   ";
    fout << setw(20) << "N_out MC  ";
    fout << setw(20) << "N_out (data/MC)" << endl;


    for(Int_t imass=0; imass<nmass; imass++) {
      fout << setw(4) << mH[imass] << "   ";
    
      //
      // compute Routin from MC
      //

      Int_t MetBinToComputeRoutin = nbins-1;

      Double_t nout_ee   = 0;
      Double_t errout_ee = 0;
      Double_t nout_mm   = 0;
      Double_t errout_mm = 0;
      Double_t nout_ll   = 0;
      Double_t errout_ll = 0;      
      Double_t nin_ee   = 0;
      Double_t errin_ee = 0;
      Double_t nin_mm   = 0;
      Double_t errin_mm = 0;
      Double_t nin_ll   = 0;
      Double_t errin_ll = 0;      
      Double_t Ree        = 0;
      Double_t ReeErrStat = 0;
      Double_t ReeErrSyst = 0;      
      Double_t Rmm        = 0;
      Double_t RmmErrStat = 0;
      Double_t RmmErrSyst = 0;
      Double_t Rll = 0;
      Double_t RllErrStat = 0;
      Double_t RllErrSyst = 0;

      if(useRFromData[jetIndex] == false){
	nout_ee   = hNout_ree_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errout_ee = hNout_ree_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nout_mm   = hNout_rmm_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errout_mm = hNout_rmm_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nout_ll   = nout_ee+nout_mm;
	errout_ll = sqrt(errout_ee*errout_ee + errout_mm*errout_mm);

	nin_ee   = hNin_ree_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errin_ee = hNin_ree_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nin_mm   = hNin_rmm_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errin_mm = hNin_rmm_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nin_ll   = nin_ee+nin_mm;
	errin_ll = sqrt(errin_ee*errin_ee + errin_mm*errin_mm);

	Ree        = nout_ee/nin_ee;
	ReeErrStat = Ree*sqrt(errin_ee*errin_ee/nin_ee/nin_ee + errout_ee*errout_ee/nout_ee/nout_ee);
	ReeErrSyst = computeSyst(hNout_ree_mc[jetIndex][imass],hNin_ree_mc[jetIndex][imass], MetBinToComputeRoutin, 0.4);    

	Rmm        = nout_mm/nin_mm;
	RmmErrStat = Rmm*sqrt(errin_mm*errin_mm/nin_mm/nin_mm + errout_mm*errout_mm/nout_mm/nout_mm);
	RmmErrSyst = computeSyst(hNout_rmm_mc[jetIndex][imass],hNin_rmm_mc[jetIndex][imass], MetBinToComputeRoutin, 0.4);

	TH1D *hout = (TH1D*)hNout_ree_mc[jetIndex][imass]->Clone("hout");
	hout->Add(hNout_rmm_mc[jetIndex][imass]);
	TH1D *hin = (TH1D*)hNin_ree_mc[jetIndex][imass]->Clone("hin");
	hin->Add(hNin_rmm_mc[jetIndex][imass]);
	Rll        = nout_ll/nin_ll;
	RllErrStat = Rll*sqrt(errin_ll*errin_ll/nin_ll/nin_ll + errout_ll*errout_ll/nout_ll/nout_ll);
	RllErrSyst = computeSyst(hout,hin, MetBinToComputeRoutin, 0.4);

	delete hout;
	delete hin;
      } else {
	nout_ll   = hNout_ree_da[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errout_ll = hNout_ree_da[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);

	nin_ll   = hNin_ree_da[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errin_ll = hNin_ree_da[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);

	TH1D *hout = (TH1D*)hNout_ree_da[jetIndex][imass]->Clone("hout");
	TH1D *hin = (TH1D*)hNin_ree_da[jetIndex][imass]->Clone("hin");

	Rll        = nout_ll/nin_ll;
	RllErrStat = Rll*sqrt(errin_ll*errin_ll/nin_ll/nin_ll + errout_ll*errout_ll/nout_ll/nout_ll);
	Double_t rErrorMax = 100000.4;
	cout << "M: " << mH[imass] << endl;
	RllErrSyst = computeSyst(hout,hin,MetBinToComputeRoutin,rErrorMax,true);
        //while(RllErrSyst <= 0.0 && rErrorMax < 1.0){
	//  rErrorMax += 0.1;
	//  RllErrSyst = computeSyst(hout,hin,MetBinToComputeRoutin,rErrorMax);
	//}  
	delete hout;
	delete hin;
      }

      sprintf(buffer,"%.2f +/- %.2f +/- %.2f",Rll,RllErrStat,RllErrSyst);
//     sprintf(buffer,"%.3f +/- %.3f +/- %.3f",Ree,ReeErrStat,ReeErrSyst);
//     sprintf(buffer,"%.3f +/- %.3f +/- %.3f",Rmm,RmmErrStat,RmmErrSyst);
      fout << setw(25) << buffer; 

      //
      // raw in-yields in data
      //
      sprintf(buffer,"%i/%i/%i/%i",Int_t(nin_mm_data[jetIndex][imass]), Int_t(nin_me_data[jetIndex][imass]), Int_t(nin_em_data[jetIndex][imass]), Int_t(nin_ee_data[jetIndex][imass]));
      fout << setw(15) << buffer;
    
      //
      // in-yields in data after OF and VZ subtraction
      //
      Double_t nof = nin_em_data[jetIndex][imass] + nin_me_data[jetIndex][imass];
      Double_t nin_ee_sub, err_ee_sub,nin_mm_sub,err_mm_sub,nin_ll_sub,err_ll_sub;

      if (ZWindowSubtractionMethod == 0) {
        //use Opposite Flavor data for bkg subtraction
        nin_ee_sub = nin_ee_data[jetIndex][imass] - 0.5*k*nof - nin_ee_vz[jetIndex][imass];
        err_ee_sub = sqrt(nin_ee_data[jetIndex][imass] + 0.5*0.5*k*k*nof*nof*(kerr*kerr/k/k + 1.0/nof) + varin_ee_vz[jetIndex][imass] + pow(nin_ee_vz[jetIndex][imass]*vzNormSystematic,2) );
        nin_mm_sub = nin_mm_data[jetIndex][imass] - 0.5/k*nof - nin_mm_vz[jetIndex][imass];
        err_mm_sub = sqrt(nin_mm_data[jetIndex][imass] + 0.5*0.5/k/k*nof*nof*(kerr*kerr/k/k + 1.0/nof) + varin_mm_vz[jetIndex][imass] + pow(nin_mm_vz[jetIndex][imass]*vzNormSystematic,2));
        nin_ll_sub = nin_ee_sub + nin_mm_sub;
        err_ll_sub = sqrt(nin_mm_data[jetIndex][imass] + nin_ee_data[jetIndex][imass] 
                                   + 0.5*0.5*( (k+1.0/k)*(k+1.0/k)*nof*nof*kerr*kerr + (k+1.0/k)*(k+1.0/k)*nof ) 
                                   + varin_mm_vz[jetIndex][imass] + varin_ee_vz[jetIndex][imass] + pow((nin_ee_vz[jetIndex][imass]+nin_mm_vz[jetIndex][imass])*vzNormSystematic,2));
        if(nin_ee_sub <= 0) nin_ee_sub = 1;
        if(nin_mm_sub <= 0) nin_mm_sub = 1;
        if(nin_ll_sub <= 0) nin_ll_sub = 1;

        if (verbose) {
          sprintf(buffer,"%.2f +/- %.2f : %.2f : %.2f + %.2f  ",nin_ll_sub,err_ll_sub, 0.5*k*nof + 0.5/k*nof, nin_ee_vz[jetIndex][imass], nin_mm_vz[jetIndex][imass]);
        }

      } else if (ZWindowSubtractionMethod == 1) {
        //Use data corrected Monte Carlo for bkg subtraction
        nin_ee_sub = nin_ee_data[jetIndex][imass] - nin_ee_bkg[jetIndex][imass] - nin_ee_vz[jetIndex][imass];
        err_ee_sub = sqrt(nin_ee_data[jetIndex][imass] +  
                          varin_ee_bkg[jetIndex][imass] + varin_ee_vz[jetIndex][imass] +
                          pow(nin_ee_ww[jetIndex][imass]*(1.0 - WWBkgScaleFactorKappaCutBased(115,min(int(jetIndex),1))),2) + 
                          pow(nin_ee_wjets[jetIndex][imass]*0.36,2) + 
                          pow(nin_ee_top[jetIndex][imass]*(1.0 - TopBkgScaleFactorKappa(jetIndex)),2) + 
                          pow(nin_ee_dytt[jetIndex][imass]*0.20,2) + 
                          pow(nin_ee_vz[jetIndex][imass]*vzNormSystematic,2)
          );
        nin_mm_sub = nin_mm_data[jetIndex][imass] - nin_mm_bkg[jetIndex][imass] - nin_mm_vz[jetIndex][imass];
        err_mm_sub = sqrt(nin_mm_data[jetIndex][imass] +  
                          varin_mm_bkg[jetIndex][imass] + varin_mm_vz[jetIndex][imass] +
                          pow(nin_mm_ww[jetIndex][imass]*(1.0 - WWBkgScaleFactorKappaCutBased(115,min(int(jetIndex),1))),2) + 
                          pow(nin_mm_wjets[jetIndex][imass]*0.36,2) + 
                          pow(nin_mm_top[jetIndex][imass]*(1.0 - TopBkgScaleFactorKappa(jetIndex)),2) + 
                          pow(nin_mm_dytt[jetIndex][imass]*0.20,2) + 
                          pow(nin_mm_vz[jetIndex][imass]*vzNormSystematic,2)
          );
        nin_ll_sub = nin_ee_sub + nin_mm_sub;
        err_ll_sub = sqrt(nin_ee_data[jetIndex][imass] + nin_mm_data[jetIndex][imass] + 
                          varin_ee_bkg[jetIndex][imass] +  varin_mm_bkg[jetIndex][imass] + 
                          varin_ee_vz[jetIndex][imass] + varin_mm_vz[jetIndex][imass] + 
                          pow((nin_ee_ww[jetIndex][imass]+nin_mm_ww[jetIndex][imass])*(1.0 - WWBkgScaleFactorKappaCutBased(115,min(int(jetIndex),1))),2) + 
                          pow((nin_ee_wjets[jetIndex][imass]+nin_mm_wjets[jetIndex][imass])*0.36,2) + 
                          pow((nin_ee_top[jetIndex][imass]+nin_mm_top[jetIndex][imass])*(1.0 - TopBkgScaleFactorKappa(jetIndex)),2) + 
                          pow((nin_ee_dytt[jetIndex][imass]+nin_mm_dytt[jetIndex][imass])*0.20,2) + 
                          pow((nin_ee_vz[jetIndex][imass] + nin_mm_vz[jetIndex][imass])*vzNormSystematic , 2)
          );

        if (verbose) {
          sprintf(buffer," %.2f +/- %.2f (%.2f +/- %.2f [bkg] : %.2f +/- %.2f [ww]   %.2f +/- %.2f [wjets]  %.2f +/- %.2f [top] | %.2f +/- %.2f [vz] )",nin_ll_sub,err_ll_sub, 
                  nin_ee_bkg[jetIndex][imass]+nin_mm_bkg[jetIndex][imass] , 
                  sqrt(varin_ee_bkg[jetIndex][imass] +  varin_mm_bkg[jetIndex][imass] + 
                       pow((nin_ee_ww[jetIndex][imass]+nin_mm_ww[jetIndex][imass])*(1.0 - WWBkgScaleFactorKappaCutBased(115,min(int(jetIndex),1))),2) + 
                       pow((nin_ee_wjets[jetIndex][imass]+nin_mm_wjets[jetIndex][imass])*0.36,2) + 
                       pow((nin_ee_top[jetIndex][imass]+nin_mm_top[jetIndex][imass])*(1.0 - TopBkgScaleFactorKappa(jetIndex)),2) + 
                       pow((nin_ee_dytt[jetIndex][imass]+nin_mm_dytt[jetIndex][imass])*0.20,2)
                    ),
                  nin_ee_ww[jetIndex][imass]+nin_mm_ww[jetIndex][imass] , 
                  sqrt(varin_ee_ww[jetIndex][imass] + varin_mm_ww[jetIndex][imass] + 
                       pow((nin_ee_ww[jetIndex][imass]+nin_mm_ww[jetIndex][imass])*(1.0 - WWBkgScaleFactorKappaCutBased(115,min(int(jetIndex),1))),2) 
                    ),
                  nin_ee_wjets[jetIndex][imass]+nin_mm_wjets[jetIndex][imass] , 
                  sqrt(varin_ee_wjets[jetIndex][imass] + varin_mm_wjets[jetIndex][imass] +
                       pow((nin_ee_wjets[jetIndex][imass]+nin_mm_wjets[jetIndex][imass])*0.36,2)
                    ),
                  nin_ee_top[jetIndex][imass]+nin_mm_top[jetIndex][imass] , 
                  sqrt(varin_ee_top[jetIndex][imass] + varin_mm_top[jetIndex][imass] +
                       pow((nin_ee_top[jetIndex][imass]+nin_mm_top[jetIndex][imass])*(1.0 - TopBkgScaleFactorKappa(jetIndex)),2)
                    ),
                  nin_ee_vz[jetIndex][imass]+nin_mm_vz[jetIndex][imass] ,
                  sqrt(varin_ee_vz[jetIndex][imass] + varin_mm_vz[jetIndex][imass] +                           
                       pow((nin_ee_vz[jetIndex][imass] + nin_mm_vz[jetIndex][imass])*vzNormSystematic , 2)
                       
                    )
            );
        }

      }
    
      if (!verbose) { sprintf(buffer,"%.2f +/- %.2f ",nin_ll_sub,err_ll_sub);}
//    sprintf(buffer,"%.2f +/- %.2f",nin_ee_sub,err_ee_sub);
//    sprintf(buffer,"%.2f +/- %.2f",nin_mm_sub,err_mm_sub);
      fout << setw(20) << buffer;
    
      //
      // in-yield data/MC scale factor
      //

      Double_t sfin_ee     = nin_ee_sub/nin_ee_dy[jetIndex][imass];
      Double_t sfin_ee_err = sfin_ee*sqrt(err_ee_sub*err_ee_sub/nin_ee_sub/nin_ee_sub + (varin_ee_dy[jetIndex][imass])*(varin_ee_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]));
      Double_t sfin_mm     = nin_mm_sub/nin_mm_dy[jetIndex][imass];
      Double_t sfin_mm_err = sfin_mm*sqrt(err_mm_sub*err_mm_sub/nin_mm_sub/nin_mm_sub + (varin_mm_dy[jetIndex][imass])*(varin_mm_dy[jetIndex][imass])/(nin_mm_dy[jetIndex][imass])/(nin_mm_dy[jetIndex][imass]));
      Double_t sfin_ll     = nin_ll_sub/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass]);
//       Double_t sfin_ll_err = sfin_ll*sqrt(err_ll_sub*err_ll_sub/nin_ll_sub/nin_ll_sub 
//                                           + (varin_ee_dy[jetIndex][imass]+varin_mm_dy[jetIndex][imass])*(varin_ee_dy[jetIndex][imass]+varin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass]));
      Double_t sfin_ll_err = sfin_ll*sqrt(err_ll_sub*err_ll_sub/nin_ll_sub/nin_ll_sub 
                                          + (varin_ee_dy[jetIndex][imass]+varin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass]));
   
      sprintf(buffer,"%.2f +/- %.2f",sfin_ll,sfin_ll_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfin_ee,sfin_ee_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfin_mm,sfin_mm_err);
      fout << setw(20) << buffer;
    
      //
      // out-yield prediction
      //
      Double_t nout_ee_pre = Ree*nin_ee_sub;
      Double_t nout_ee_sta = nout_ee_pre*sqrt(ReeErrStat*ReeErrStat/Ree/Ree + err_ee_sub*err_ee_sub/nin_ee_sub/nin_ee_sub);
      Double_t nout_ee_sys = nout_ee_pre*ReeErrSyst/Ree;
      Double_t nout_ee_err = sqrt(nout_ee_sys*nout_ee_sys + nout_ee_sta*nout_ee_sta);
      Double_t nout_mm_pre = Rmm*nin_mm_sub;
      Double_t nout_mm_sta = nout_mm_pre*sqrt(RmmErrStat*RmmErrStat/Rmm/Rmm + err_mm_sub*err_mm_sub/nin_mm_sub/nin_mm_sub);
      Double_t nout_mm_sys = nout_mm_pre*RmmErrSyst/Rmm;
      Double_t nout_mm_err = sqrt(nout_mm_sys*nout_mm_sys + nout_mm_sta*nout_mm_sta);
      Double_t nout_ll_pre = Rll*nin_ll_sub;
      Double_t nout_ll_sta = nout_ll_pre*sqrt((RllErrStat*RllErrStat)/Rll/Rll + err_ll_sub*err_ll_sub/nin_ll_sub/nin_ll_sub);
      Double_t nout_ll_sys = nout_ll_pre*RllErrSyst/Rll;
      Double_t nout_ll_err = sqrt(nout_ll_sys*nout_ll_sys + nout_ll_sta*nout_ll_sta);

      sprintf(buffer,"%.2f +/- %.2f +/- %.2f",nout_ll_pre,nout_ll_sta,nout_ll_sys);
//    sprintf(buffer,"%.2f +/- %.2f +/- %.2f",nout_ee_pre,nout_ee_sta,nout_ee_sys);
//    sprintf(buffer,"%.2f +/- %.2f +/- %.2f",nout_mm_pre,nout_mm_sta,nout_mm_sys);
      fout << "   " << setw(20) << buffer;
    
      //
      // out-yield in MC
      //
      sprintf(buffer,"%.2f +/- %.2f",nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass],sqrt(varout_ee_dy[jetIndex][imass] + varout_mm_dy[jetIndex][imass]));
//    sprintf(buffer,"%.2f +/- %.2f",nout_ee_dy[jetIndex][imass],sqrt(varout_ee_dy[jetIndex][imass]));
//    sprintf(buffer,"%.2f +/- %.2f",nout_mm_dy[jetIndex][imass],sqrt(varout_mm_dy[jetIndex][imass]));
      fout << setw(20) << buffer;
    
      //
      // out-yield data/MC scale factor
      //
      Double_t sfout_ee     = nout_ee_pre/nout_ee_dy[jetIndex][imass];
      Double_t sfout_ee_err = sfout_ee*sqrt(nout_ee_err*nout_ee_err/nout_ee_pre/nout_ee_pre );
      Double_t sfout_mm     = nout_mm_pre/nout_mm_dy[jetIndex][imass];
      Double_t sfout_mm_err = sfout_mm*sqrt(nout_mm_err*nout_mm_err/nout_mm_pre/nout_mm_pre );
      Double_t sfout_ll     = nout_ll_pre/(nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass]);
      Double_t sfout_ll_err = sfout_ll*sqrt(nout_ll_err*nout_ll_err/nout_ll_pre/nout_ll_pre);
    
      sprintf(buffer,"%.2f +/- %.2f",sfout_ll,sfout_ll_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfout_ee,sfout_ee_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfout_mm,sfout_mm_err);
      fout << setw(20) << buffer;        
      fout << endl;
      if (imass == 0) {
        tmpDYBkgScaleFactorWWPreselection    = sfout_ll;
        tmpDYBkgScaleFactorWWPreselectionErr = sfout_ll_err;
      } else {
        // BIG CHANGE, we quote
        //tmpDYBkgScaleFactorHiggsSelection.push_back(sfout_ll);
        //tmpDYBkgScaleFactorHiggsSelectionErr.push_back(sfout_ll_err);
        tmpDYBkgScaleFactorHiggsSelection.push_back(nout_ll_pre);
        tmpDYBkgScaleFactorHiggsSelectionErr.push_back(sqrt(nout_ll_sta*nout_ll_sta+nout_ll_sys*nout_ll_sys));
      }
    }

    DYBkgScaleFactorHiggsSelection.push_back(tmpDYBkgScaleFactorHiggsSelection);
    DYBkgScaleFactorHiggsSelectionErr.push_back(tmpDYBkgScaleFactorHiggsSelectionErr);
    DYBkgScaleFactorWWPreselection.push_back(tmpDYBkgScaleFactorWWPreselection);
    DYBkgScaleFactorWWPreselectionErr.push_back(tmpDYBkgScaleFactorWWPreselectionErr);
  }

  fout.close();


  //***************************************************************************
  // Generate DY Scale Factor and Systematics Code for card creation
  //***************************************************************************
  char DYBkgScaleFactorsName[200];
  sprintf(DYBkgScaleFactorsName,"DYBkgScaleFactors.h");
  if(forBDTAna == true) sprintf(DYBkgScaleFactorsName,"DYBkgScaleFactorsBDT.h");
  ofstream outf(DYBkgScaleFactorsName);

  if(forBDTAna == true) outf << "static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {" << endl;
  else                  outf << "static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {" << endl;
  
  outf << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << mH[i+1];
    if (i < nmass-1-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t DYBkgScaleFactorWWPreselection[3] = { " 
       << DYBkgScaleFactorWWPreselection[0] << ", "
       << DYBkgScaleFactorWWPreselection[1] << ", "
       << DYBkgScaleFactorWWPreselection[2] << " "
       << " };" << endl;

  outf << "  Double_t DYBkgScaleFactorHiggsSelection[3][" << nmass-1 << "] = { " << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << DYBkgScaleFactorHiggsSelection[0][i];    
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << DYBkgScaleFactorHiggsSelection[1][i];    
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << DYBkgScaleFactorHiggsSelection[2][i];    
    if (i < nmass-1-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];" << endl;
  
  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return DYBkgScaleFactorWWPreselection[jetBin];" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;


  if(forBDTAna == true) outf << "static Double_t DYBkgScaleFactorBDTKappa(Int_t mH, Int_t jetBin) {" << endl;
  else                  outf << "static Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {" << endl;
  
  outf << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << mH[i+1];
    if (i < nmass-1-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { " 
       << (1.0 + DYBkgScaleFactorWWPreselectionErr[0]/DYBkgScaleFactorWWPreselection[0]) << ", "
       << (1.0 + DYBkgScaleFactorWWPreselectionErr[1]/DYBkgScaleFactorWWPreselection[1]) << ", "
       << (1.0 + DYBkgScaleFactorWWPreselectionErr[2]/DYBkgScaleFactorWWPreselection[2]) << " "
       << " };" << endl;

  outf << "  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][" << nmass-1 << "] = { " << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << (1.0 + DYBkgScaleFactorHiggsSelectionErr[0][i] / DYBkgScaleFactorHiggsSelection[0][i]);
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << (1.0 + DYBkgScaleFactorHiggsSelectionErr[1][i] / DYBkgScaleFactorHiggsSelection[1][i]);
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << (1.0 + DYBkgScaleFactorHiggsSelectionErr[2][i] / DYBkgScaleFactorHiggsSelection[2][i]);
    if (i < nmass-1-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];" << endl;
  
  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return DYBkgScaleFactorWWPreselectionKappa[jetBin];" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;

  outf.close();



}



//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
    Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi) 
{
  Double_t dphi = acos(cos(lepPhi-metPhi));
  if(dphi > 0.5*TMath::Pi())
    return met;
    
  return met*sin(dphi);
}

//--------------------------------------------------------------------------------------------------
Double_t computeSyst(const TH1D *hout, const TH1D *hin, Int_t binUsed, Double_t rErrorMax, Bool_t isDebug)
{
  const Int_t nbins = hout->GetXaxis()->GetNbins();
  Double_t r0=(hout->GetBinContent(binUsed))/(hin->GetBinContent(binUsed));
  Double_t dr=0;
  for(Int_t ibin=1; ibin<=nbins; ibin++) {
    Double_t r = (hout->GetBinContent(ibin))/(hin->GetBinContent(ibin));

    //only consider last bin for systematics if the uncertainty is less than rErrorMax to protect against statistical fluctuations
    Double_t rError = sqrt( pow(hout->GetBinError(ibin)/hout->GetBinContent(ibin),2) + pow(hin->GetBinError(ibin)/hin->GetBinContent(ibin),2));

    if(r < 0.0) r = 0.0;
    if(isDebug == true) printf("bin(%1d) Nin = %8.3f Nout = %8.3f R0 = %5.3f R = %5.3f dr = %5.3f drMax = %5.3f rError = %5.3f\n",ibin,hin->GetBinContent(ibin),hout->GetBinContent(ibin),r0,r,fabs(r-r0),dr,rError);

    if(rError > rErrorMax) continue;

    if(fabs(r-r0) > dr) dr = fabs(r-r0);
  }
  return dr;
}

TGraphErrors* MakeRoutinGraph(const TH1D *hout, const TH1D *hin, string graphname) {

  Int_t nbins = hout->GetNbinsX();
  assert(nbins <= 200);

  Double_t x[200];
  Double_t y[200];
  Double_t xErr[200];
  Double_t yErr[200];

  for (int i=0;i < 200; i++) {
    x[i] = 0;
    y[i] = 0;
    xErr[i] = 0;
    yErr[i] = 0;
  }

  //don't take the overflow bins
  for (int b=0; b<nbins ; ++b) {

    x[b] = hout->GetXaxis()->GetBinCenter(b+1);    
    xErr[b] = 0.0;

    Double_t ratio;
    Double_t err;

    Double_t n1 = hout->GetBinContent(b+1);
    Double_t n2 = hin->GetBinContent(b+1);
    if (n2 > 0) ratio = n1 / n2;
    else ratio = 0;
    if (n1 > 0 && n2 > 0) err = ratio * sqrt( pow(hout->GetBinError(b+1) / n1 , 2) + pow( hin->GetBinError(b+1) / n2,2) );    
    else err = 0;    

    //cerr << " done bin " << b << " " << x[b] << " : " << n1 << " +/-" << hout->GetBinError(b+1) << " / " << n2 << " +/- " << hin->GetBinError(b+1) << " = " << ratio << " " << err << endl;
    y[b] = ratio;
    yErr[b] = err;
   }

  TGraphErrors *graphRatio = new TGraphErrors(nbins, x, y, xErr, yErr );
  graphRatio->SetName(graphname.c_str());
  graphRatio->SetTitle(graphname.c_str());
  graphRatio->GetXaxis()->SetTitle(hout->GetXaxis()->GetTitle());
  graphRatio->GetYaxis()->SetTitle("R_{out/in}");

  graphRatio->SetMarkerSize(1);
  graphRatio->SetLineWidth(2);

  return graphRatio;
}



TGraphErrors* MakeRoutinGraphDataDrivenOFSubtraction(const TH1D *hout_ee_data, const TH1D *hin_ee_data, 
                                                     const TH1D *hout_mm_data, const TH1D *hin_mm_data, 
                                                     const TH1D *hout_ee_vz, const TH1D *hin_ee_vz, 
                                                     const TH1D *hout_mm_vz, const TH1D *hin_mm_vz, 
                                                     const TH1D *hout_OF_data, const TH1D *hin_OF_data, 
                                                     Double_t EleToMuEffRatio, Double_t EleToMuEffRatioErr,
                                                     Int_t FinalState, string graphname) {
  Double_t vzNormSystematic = 0.10;

  Int_t nbins = hout_ee_data->GetNbinsX();
  assert(nbins <= 200);

  Double_t x[200];
  Double_t y[200];
  Double_t xErr[200];
  Double_t yErr[200];

  for (int i=0;i < 200; i++) {
    x[i] = 0;
    y[i] = 0;
    xErr[i] = 0;
    yErr[i] = 0;
  }

  //don't take the overflow bins
  for (int b=0; b<nbins ; ++b) {

    x[b] = hout_ee_data->GetXaxis()->GetBinCenter(b+1);    
    xErr[b] = 0.0;

    Double_t Nout_data = 0;
    Double_t Nin_data = 0;
    Double_t Nout_OFBkg = 0;
    Double_t Nin_OFBkg = 0;
    Double_t Nout_VZ = 0;
    Double_t Nin_VZ = 0;
    Double_t Nout_BkgSubtracted = 0;
    Double_t Nin_BkgSubtracted = 0;
    Double_t NErrSqr_out_BkgSubtracted = 0;
    Double_t NErrSqr_in_BkgSubtracted = 0;
    
    if (FinalState == 0) {
//       cout << "Graph : " << graphname << endl;
      Nout_data = hout_ee_data->GetBinContent(b+1) + hout_mm_data->GetBinContent(b+1);
      Nin_data = hin_ee_data->GetBinContent(b+1) + hin_mm_data->GetBinContent(b+1);
      Nout_OFBkg = 0.5 * EleToMuEffRatio * hout_OF_data->GetBinContent(b+1) + 0.5 / EleToMuEffRatio * hout_OF_data->GetBinContent(b+1) ;
      Nin_OFBkg = 0.5 * EleToMuEffRatio * hin_OF_data->GetBinContent(b+1) + 0.5 / EleToMuEffRatio * hin_OF_data->GetBinContent(b+1);
      Nout_VZ = hout_ee_vz->GetBinContent(b+1) + hout_mm_vz->GetBinContent(b+1);
      Nin_VZ = hin_ee_vz->GetBinContent(b+1) + hout_mm_vz->GetBinContent(b+1);

      Nout_BkgSubtracted = Nout_data - Nout_OFBkg - Nout_VZ;
      Nin_BkgSubtracted  = Nin_data  - Nin_OFBkg  - Nin_VZ;

      NErrSqr_out_BkgSubtracted = sqrt(hout_ee_data->GetBinContent(b+1) + hout_mm_data->GetBinContent(b+1) 
                                       + 0.5*0.5*( (EleToMuEffRatio+1.0/EleToMuEffRatio)*(EleToMuEffRatio+1.0/EleToMuEffRatio)*
                                                   hout_OF_data->GetBinContent(b+1)*hout_OF_data->GetBinContent(b+1)*
                                                   EleToMuEffRatioErr*EleToMuEffRatioErr + 
                                                   (EleToMuEffRatio+1.0/EleToMuEffRatio)*(EleToMuEffRatio+1.0/EleToMuEffRatio)*hout_OF_data->GetBinContent(b+1) ) + 
                                       pow(hout_mm_vz->GetBinError(b+1),2) + pow(hout_ee_vz->GetBinError(b+1),2) + 
                                       pow((Nout_VZ)*vzNormSystematic,2));
      NErrSqr_in_BkgSubtracted = sqrt(hin_ee_data->GetBinContent(b+1) + hin_mm_data->GetBinContent(b+1) 
                                       + 0.5*0.5*( (EleToMuEffRatio+1.0/EleToMuEffRatio)*(EleToMuEffRatio+1.0/EleToMuEffRatio)*
                                                   hin_OF_data->GetBinContent(b+1)*hin_OF_data->GetBinContent(b+1)*
                                                   EleToMuEffRatioErr*EleToMuEffRatioErr + 
                                                   (EleToMuEffRatio+1.0/EleToMuEffRatio)*(EleToMuEffRatio+1.0/EleToMuEffRatio)*hin_OF_data->GetBinContent(b+1) ) + 
                                       pow(hin_mm_vz->GetBinError(b+1),2) + pow(hin_ee_vz->GetBinError(b+1),2) + 
                                       pow((Nin_VZ)*vzNormSystematic,2));
//       cout << EleToMuEffRatio << " : " << hout_OF_data->GetBinContent(b+1) << endl;
//       cout << EleToMuEffRatio << " : " << hin_OF_data->GetBinContent(b+1) << endl;
//       cout << Nout_data << " - " << Nout_OFBkg << " - " << Nout_VZ << " = " << Nout_BkgSubtracted << " +/- " << NErrSqr_out_BkgSubtracted << endl;
//       cout << Nin_data << " - " << Nin_OFBkg << " - " << Nin_VZ << " = " << Nin_BkgSubtracted << " +/- " << NErrSqr_in_BkgSubtracted << endl;
      
    }
    if (FinalState == 1) {
      Nout_data = hout_ee_data->GetBinContent(b+1);
      Nin_data = hin_ee_data->GetBinContent(b+1);

      Nout_OFBkg = 0.5 * EleToMuEffRatio * hout_OF_data->GetBinContent(b+1);
      Nin_OFBkg = 0.5 * EleToMuEffRatio * hin_OF_data->GetBinContent(b+1);
      Nout_VZ = hout_ee_vz->GetBinContent(b+1);
      Nin_VZ = hin_ee_vz->GetBinContent(b+1);

      Nout_BkgSubtracted = Nout_data - Nout_OFBkg - Nout_VZ;
      Nin_BkgSubtracted  = Nin_data  - Nin_OFBkg  - Nin_VZ;

      NErrSqr_out_BkgSubtracted = sqrt(hout_ee_data->GetBinContent(b+1) 
                                       + 0.5*0.5*EleToMuEffRatio*EleToMuEffRatio*
                                       hout_OF_data->GetBinContent(b+1)*hout_OF_data->GetBinContent(b+1)*
                                       ( EleToMuEffRatioErr*EleToMuEffRatioErr/EleToMuEffRatio/EleToMuEffRatio + 
                                         1.0 / hout_OF_data->GetBinContent(b+1) )  + 
                                       pow(hout_ee_vz->GetBinError(b+1),2) + pow((Nout_VZ)*vzNormSystematic,2));
      NErrSqr_in_BkgSubtracted = sqrt(hin_ee_data->GetBinContent(b+1) 
                                      + 0.5*0.5*EleToMuEffRatio*EleToMuEffRatio*
                                      hin_OF_data->GetBinContent(b+1)*hin_OF_data->GetBinContent(b+1)*
                                      ( EleToMuEffRatioErr*EleToMuEffRatioErr/EleToMuEffRatio/EleToMuEffRatio + 
                                        1.0 / hin_OF_data->GetBinContent(b+1) )  + 
                                      pow(hin_ee_vz->GetBinError(b+1),2) + pow((Nin_VZ)*vzNormSystematic,2));


    }
    if (FinalState == 2) {
      Nout_data = hout_mm_data->GetBinContent(b+1);
      Nin_data = hin_mm_data->GetBinContent(b+1);

      Nout_OFBkg = 0.5 / EleToMuEffRatio * hout_OF_data->GetBinContent(b+1) ;
      Nin_OFBkg = 0.5 / EleToMuEffRatio * hin_OF_data->GetBinContent(b+1);
      Nout_VZ = hout_mm_vz->GetBinContent(b+1);
      Nin_VZ = hout_mm_vz->GetBinContent(b+1);

      Nout_BkgSubtracted = Nout_data - Nout_OFBkg - Nout_VZ;
      Nin_BkgSubtracted  = Nin_data  - Nin_OFBkg  - Nin_VZ;

      NErrSqr_out_BkgSubtracted = sqrt(hout_mm_data->GetBinContent(b+1) 
                                       + 0.5*0.5/EleToMuEffRatio/EleToMuEffRatio*
                                       hout_OF_data->GetBinContent(b+1)*hout_OF_data->GetBinContent(b+1)*
                                       (EleToMuEffRatioErr*EleToMuEffRatioErr/EleToMuEffRatio/EleToMuEffRatio + 
                                        1.0/hout_OF_data->GetBinContent(b+1) ) + 
                                       pow(hout_mm_vz->GetBinError(b+1),2) + pow((Nout_VZ)*vzNormSystematic,2));
      NErrSqr_in_BkgSubtracted = sqrt(hin_mm_data->GetBinContent(b+1) 
                                       + 0.5*0.5/EleToMuEffRatio/EleToMuEffRatio*
                                       hin_OF_data->GetBinContent(b+1)*hin_OF_data->GetBinContent(b+1)*
                                       (EleToMuEffRatioErr*EleToMuEffRatioErr/EleToMuEffRatio/EleToMuEffRatio + 
                                        1.0/hin_OF_data->GetBinContent(b+1) ) + 
                                       pow(hin_mm_vz->GetBinError(b+1),2) + pow((Nin_VZ)*vzNormSystematic,2));
    }


    Double_t ratio;
    Double_t err;

    Double_t n1 = Nout_BkgSubtracted;
    Double_t n2 = Nin_BkgSubtracted;
    if (n2 > 0) ratio = n1 / n2;
    else ratio = 0;
    if (n1 > 0 && n2 > 0) err = ratio * sqrt( NErrSqr_out_BkgSubtracted / pow(n1,2) + NErrSqr_in_BkgSubtracted / pow(n2,2) );    
    else err = 0;    

    if (FinalState == 0) {
//       cout << " done bin " << b << " " << x[b] << " : " << n1 << " +/-" << sqrt(NErrSqr_out_BkgSubtracted) 
//            << " / " << n2 << " +/- " << sqrt(NErrSqr_in_BkgSubtracted) << " = " 
//            << ratio << " " << err << endl;
    }
    y[b] = ratio;
    yErr[b] = err;
   }

  TGraphErrors *graphRatio = new TGraphErrors(nbins, x, y, xErr, yErr );
  graphRatio->SetName(graphname.c_str());
  graphRatio->SetTitle(graphname.c_str());
  graphRatio->GetXaxis()->SetTitle(hout_ee_data->GetXaxis()->GetTitle());
  graphRatio->GetYaxis()->SetTitle("R_{out/in}");

  graphRatio->SetMarkerSize(1);
  graphRatio->SetLineWidth(2);

  return graphRatio;
}


TGraphErrors* MakeRoutinGraphDataDrivenMCSubtraction(const TH1D *hout_ee_data, const TH1D *hin_ee_data, 
                                                     const TH1D *hout_mm_data, const TH1D *hin_mm_data, 
                                                     const TH1D *hout_ee_vz, const TH1D *hin_ee_vz, 
                                                     const TH1D *hout_mm_vz, const TH1D *hin_mm_vz, 
                                                     const TH1D *hout_ee_ww, const TH1D *hin_ee_ww, 
                                                     const TH1D *hout_mm_ww, const TH1D *hin_mm_ww, 
                                                     const TH1D *hout_ee_wjets, const TH1D *hin_ee_wjets, 
                                                     const TH1D *hout_mm_wjets, const TH1D *hin_mm_wjets, 
                                                     const TH1D *hout_ee_top, const TH1D *hin_ee_top, 
                                                     const TH1D *hout_mm_top, const TH1D *hin_mm_top, 
                                                     const TH1D *hout_ee_dytt, const TH1D *hin_ee_dytt, 
                                                     const TH1D *hout_mm_dytt, const TH1D *hin_mm_dytt, 
                                                     Double_t WWSystematicError, Double_t TopSystematicError,
                                                     Int_t FinalState, string graphname) {

  Double_t vzNormSystematic = 0.10;

  Int_t nbins = hout_ee_data->GetNbinsX();
  assert(nbins <= 200);

  Double_t x[200];
  Double_t y[200];
  Double_t xErr[200];
  Double_t yErr[200];

  for (int i=0;i < 200; i++) {
    x[i] = 0;
    y[i] = 0;
    xErr[i] = 0;
    yErr[i] = 0;
  }

  //don't take the overflow bins
  for (int b=0; b<nbins ; ++b) {

    x[b] = hout_ee_data->GetXaxis()->GetBinCenter(b+1);    
    xErr[b] = 0.0;

    Double_t Nout_data = 0;
    Double_t Nin_data = 0;
    Double_t Nout_BkgSubtracted = 0;
    Double_t Nin_BkgSubtracted = 0;
    Double_t NErrSqr_out_BkgSubtracted = 0;
    Double_t NErrSqr_in_BkgSubtracted = 0;
    
    if (FinalState == 0) {
      Nout_data = hout_ee_data->GetBinContent(b+1) + hout_mm_data->GetBinContent(b+1);
      Nin_data = hin_ee_data->GetBinContent(b+1) + hin_mm_data->GetBinContent(b+1);
 
      Nout_BkgSubtracted = Nout_data - ( hout_ee_vz->GetBinContent(b+1) + hout_mm_vz->GetBinContent(b+1) +
                                         hout_ee_ww->GetBinContent(b+1) + hout_mm_ww->GetBinContent(b+1) +
                                         hout_ee_wjets->GetBinContent(b+1) + hout_mm_wjets->GetBinContent(b+1) +
                                         hout_ee_top->GetBinContent(b+1) + hout_mm_top->GetBinContent(b+1) +
                                         hout_ee_dytt->GetBinContent(b+1) + hout_mm_dytt->GetBinContent(b+1)
        );
      Nin_BkgSubtracted = Nin_data - ( hin_ee_vz->GetBinContent(b+1) + hin_mm_vz->GetBinContent(b+1) +
                                       hin_ee_ww->GetBinContent(b+1) + hin_mm_ww->GetBinContent(b+1) +
                                       hin_ee_wjets->GetBinContent(b+1) + hin_mm_wjets->GetBinContent(b+1) +
                                       hin_ee_top->GetBinContent(b+1) + hin_mm_top->GetBinContent(b+1) +
                                       hin_ee_dytt->GetBinContent(b+1) + hin_mm_dytt->GetBinContent(b+1)
        );


      NErrSqr_out_BkgSubtracted = sqrt(Nout_data +  
                                       hout_ee_vz->GetBinError(b+1) +  hout_mm_vz->GetBinError(b+1) + 
                                       hout_ee_ww->GetBinError(b+1) +  hout_mm_ww->GetBinError(b+1) + 
                                       hout_ee_wjets->GetBinError(b+1) +  hout_mm_wjets->GetBinError(b+1) + 
                                       hout_ee_top->GetBinError(b+1) +  hout_mm_top->GetBinError(b+1) + 
                                       hout_ee_dytt->GetBinError(b+1) +  hout_mm_dytt->GetBinError(b+1) + 
                                       pow((hout_ee_ww->GetBinContent(b+1)+hout_mm_ww->GetBinContent(b+1))*WWSystematicError,2) + 
                                       pow((hout_ee_wjets->GetBinContent(b+1)+hout_mm_wjets->GetBinContent(b+1))*0.36,2) + 
                                       pow((hout_ee_top->GetBinContent(b+1)+hout_mm_top->GetBinContent(b+1))*TopSystematicError,2) + 
                                       pow((hout_ee_dytt->GetBinContent(b+1)+hout_mm_dytt->GetBinContent(b+1))*0.20,2) + 
                                       pow((hout_ee_vz->GetBinContent(b+1)+hout_mm_vz->GetBinContent(b+1))*vzNormSystematic , 2)
        );
      NErrSqr_in_BkgSubtracted = sqrt(Nin_data +  
                                       hin_ee_vz->GetBinError(b+1) +  hin_mm_vz->GetBinError(b+1) + 
                                       hin_ee_ww->GetBinError(b+1) +  hin_mm_ww->GetBinError(b+1) + 
                                       hin_ee_wjets->GetBinError(b+1) +  hin_mm_wjets->GetBinError(b+1) + 
                                       hin_ee_top->GetBinError(b+1) +  hin_mm_top->GetBinError(b+1) + 
                                       hin_ee_dytt->GetBinError(b+1) +  hin_mm_dytt->GetBinError(b+1) + 
                                       pow((hin_ee_ww->GetBinContent(b+1)+hin_mm_ww->GetBinContent(b+1))*WWSystematicError,2) + 
                                       pow((hin_ee_wjets->GetBinContent(b+1)+hin_mm_wjets->GetBinContent(b+1))*0.36,2) + 
                                       pow((hin_ee_top->GetBinContent(b+1)+hin_mm_top->GetBinContent(b+1))*TopSystematicError,2) + 
                                       pow((hin_ee_dytt->GetBinContent(b+1)+hin_mm_dytt->GetBinContent(b+1))*0.20,2) + 
                                       pow((hin_ee_vz->GetBinContent(b+1)+hin_mm_vz->GetBinContent(b+1))*vzNormSystematic , 2)
        );
      
      
    }
    if (FinalState == 1) {
      Nout_data = hout_ee_data->GetBinContent(b+1);
      Nin_data = hin_ee_data->GetBinContent(b+1);
 
      Nout_BkgSubtracted = Nout_data - ( hout_ee_vz->GetBinContent(b+1) +
                                         hout_ee_ww->GetBinContent(b+1) + 
                                         hout_ee_wjets->GetBinContent(b+1) + 
                                         hout_ee_top->GetBinContent(b+1) + 
                                         hout_ee_dytt->GetBinContent(b+1)  
        );
      Nin_BkgSubtracted = Nin_data - ( hin_ee_vz->GetBinContent(b+1) + 
                                       hin_ee_ww->GetBinContent(b+1) + 
                                       hin_ee_wjets->GetBinContent(b+1) +
                                       hin_ee_top->GetBinContent(b+1) + 
                                       hin_ee_dytt->GetBinContent(b+1) 
        );


      NErrSqr_out_BkgSubtracted = sqrt(Nout_data +  
                                       hout_ee_vz->GetBinError(b+1) + 
                                       hout_ee_ww->GetBinError(b+1) + 
                                       hout_ee_wjets->GetBinError(b+1) + 
                                       hout_ee_top->GetBinError(b+1) + 
                                       hout_ee_dytt->GetBinError(b+1) + 
                                       pow((hout_ee_ww->GetBinContent(b+1))*WWSystematicError,2) + 
                                       pow((hout_ee_wjets->GetBinContent(b+1))*0.36,2) + 
                                       pow((hout_ee_top->GetBinContent(b+1))*TopSystematicError,2) + 
                                       pow((hout_ee_dytt->GetBinContent(b+1))*0.20,2) + 
                                       pow((hout_ee_vz->GetBinContent(b+1))*vzNormSystematic,2)
        );
      NErrSqr_in_BkgSubtracted = sqrt(Nin_data +  
                                       hin_ee_vz->GetBinError(b+1) +
                                       hin_ee_ww->GetBinError(b+1) +
                                       hin_ee_wjets->GetBinError(b+1) +
                                       hin_ee_top->GetBinError(b+1) +
                                       hin_ee_dytt->GetBinError(b+1) +
                                       pow((hin_ee_ww->GetBinContent(b+1))*WWSystematicError,2) + 
                                       pow((hin_ee_wjets->GetBinContent(b+1))*0.36,2) + 
                                       pow((hin_ee_top->GetBinContent(b+1))*TopSystematicError,2) + 
                                       pow((hin_ee_dytt->GetBinContent(b+1))*0.20,2) + 
                                       pow((hin_ee_vz->GetBinContent(b+1))*vzNormSystematic , 2)
        );
    }
    if (FinalState == 2) {
      Nout_data = hout_mm_data->GetBinContent(b+1);
      Nin_data = hin_mm_data->GetBinContent(b+1);
 
      Nout_BkgSubtracted = Nout_data - ( hout_mm_vz->GetBinContent(b+1) +
                                         hout_mm_ww->GetBinContent(b+1) + 
                                         hout_mm_wjets->GetBinContent(b+1) + 
                                         hout_mm_top->GetBinContent(b+1) + 
                                         hout_mm_dytt->GetBinContent(b+1)  
        );
      Nin_BkgSubtracted = Nin_data - ( hin_mm_vz->GetBinContent(b+1) + 
                                       hin_mm_ww->GetBinContent(b+1) + 
                                       hin_mm_wjets->GetBinContent(b+1) +
                                       hin_mm_top->GetBinContent(b+1) + 
                                       hin_mm_dytt->GetBinContent(b+1) 
        );


      NErrSqr_out_BkgSubtracted = sqrt(Nout_data +  
                                       hout_mm_vz->GetBinError(b+1) + 
                                       hout_mm_ww->GetBinError(b+1) + 
                                       hout_mm_wjets->GetBinError(b+1) + 
                                       hout_mm_top->GetBinError(b+1) + 
                                       hout_mm_dytt->GetBinError(b+1) + 
                                       pow((hout_mm_ww->GetBinContent(b+1))*WWSystematicError,2) + 
                                       pow((hout_mm_wjets->GetBinContent(b+1))*0.36,2) + 
                                       pow((hout_mm_top->GetBinContent(b+1))*TopSystematicError,2) + 
                                       pow((hout_mm_dytt->GetBinContent(b+1))*0.20,2) + 
                                       pow((hout_mm_vz->GetBinContent(b+1))*vzNormSystematic,2)
        );
      NErrSqr_in_BkgSubtracted = sqrt(Nin_data +  
                                       hin_mm_vz->GetBinError(b+1) +
                                       hin_mm_ww->GetBinError(b+1) +
                                       hin_mm_wjets->GetBinError(b+1) +
                                       hin_mm_top->GetBinError(b+1) +
                                       hin_mm_dytt->GetBinError(b+1) +
                                       pow((hin_mm_ww->GetBinContent(b+1))*WWSystematicError,2) + 
                                       pow((hin_mm_wjets->GetBinContent(b+1))*0.36,2) + 
                                       pow((hin_mm_top->GetBinContent(b+1))*TopSystematicError,2) + 
                                       pow((hin_mm_dytt->GetBinContent(b+1))*0.20,2) + 
                                       pow((hin_mm_vz->GetBinContent(b+1))*vzNormSystematic , 2)
        );
    }


    Double_t ratio;
    Double_t err;

    Double_t n1 = Nout_BkgSubtracted;
    Double_t n2 = Nin_BkgSubtracted;
    if (n2 > 0) ratio = n1 / n2;
    else ratio = 0;
    if (n1 > 0 && n2 > 0) err = ratio * sqrt( NErrSqr_out_BkgSubtracted / pow(n1,2) + NErrSqr_in_BkgSubtracted / pow(n2,2) );    
    else err = 0;    

//     cerr << " done bin " << b << " " << x[b] << " : " << n1 << " +/-" << hout->GetBinError(b+1) << " / " << n2 << " +/- " << hin->GetBinError(b+1) << " = " << ratio << " " << err << endl;
    y[b] = ratio;
    yErr[b] = err;
   }

  TGraphErrors *graphRatio = new TGraphErrors(nbins, x, y, xErr, yErr );
  graphRatio->SetName(graphname.c_str());
  graphRatio->SetTitle(graphname.c_str());
  graphRatio->GetXaxis()->SetTitle(hout_ee_data->GetXaxis()->GetTitle());
  graphRatio->GetYaxis()->SetTitle("R_{out/in}");

  graphRatio->SetMarkerSize(1);
  graphRatio->SetLineWidth(2);

  return graphRatio;
}
