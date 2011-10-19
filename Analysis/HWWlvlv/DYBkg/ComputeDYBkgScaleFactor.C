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
#include "Smurf/Analysis/HZZllvv/PileupReweighting.h"
#include "Smurf/Core/SmurfTree.h"
#endif

static const unsigned int basic_selection = SmurfTree::BaseLine | 
                                            SmurfTree::ChargeMatch | 
                                            SmurfTree::Lep1FullSelection |
                                            SmurfTree::Lep2FullSelection |
                                            SmurfTree::ExtraLeptonVeto |
                                            SmurfTree::TopVeto;


//=== FUNCTION DECLARATIONS ======================================================================================

// compute projected MET
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi);

// compute systematic uncertainty
Double_t computeSyst(const TH1F *hout, const TH1F *hin, Int_t binUsed);
TGraphErrors* MakeRoutinGraph(const TH1F *hout, const TH1F *hin, string graphname, vector<Double_t> bins);

//=== MAIN MACRO =================================================================================================

void ComputeDYBkgScaleFactor(Bool_t useRecoilModel = kFALSE)
{

  //*******************************************************
  // Settings 
  //*******************************************************
  const TString format("png");
//  const Double_t lumi = 1.545;
  const Double_t lumi = 2.121;
  enum { kMuMu, kEleEle, kEleMu, kMuEle };
  
  Int_t nmet = 1;
  if(useRecoilModel) nmet = 100;

  bool isOldAna = true;
  TString filesPath  = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets";
  if(isOldAna == false){
    filesPath  = "/data/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets";
  }

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  

  infilenamev.push_back(Form("%s/data_2l.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/dyee.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/dymm.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/wz.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/zz.root",filesPath.Data()));

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
  
  const Int_t nbins = 5;  
  const Float_t bins[nbins+1] = {20, 22, 26, 31, 40, 50};
  vector<Double_t> binEdges;
  for (UInt_t k=0; k<nbins+1; ++k) binEdges.push_back(bins[k]);
  
  const Int_t nmass = 13;
  const Double_t mH[nmass] = {0,115,120,130,140,150,160,170,180,190,200,250,300};  

  Double_t _Pt1[nmass]    = {	20,  20,  20,  25,  25,  27,  30,  34,  36,  38,  40,  55,  70 };
  Double_t _Pt2[nmass]    = {	10,  10,  10,  10,  15,  25,  25,  25,  25,  25,  25,  25,  25 };
  Double_t _Mll[nmass]    = { 7000,  40,  40,  45,  45,  50,  50,  50,  60,  80,  90, 150, 200 };
  Double_t _DPhi[nmass]   = {  180, 115, 115,  90,  90,  90,  60,  60,  70,  90, 100, 140, 175 };
  Double_t _MtLow[nmass]  = {	 0,  70,  70,  75,  80,  80,  90, 110, 120, 120, 120, 120, 120 };
  Double_t _MtHigh[nmass] = { 7000, 110, 120, 125, 130, 150, 160, 170, 180, 190, 200, 250, 300 };
    
  //*******************************************************
  //Yields and  histograms
  //*******************************************************
  vector<Double_t> nin_kee_data, nin_kmm_data;
  
  vector<vector<TH1F*> > hNin_ree_mc,   hNout_ree_mc,   hNin_rmm_mc,   hNout_rmm_mc;
  vector<vector<TH1F*> > hNin_ree_data, hNout_ree_data, hNin_rmm_data, hNout_rmm_data;

  vector<vector<Double_t> > nin_ee_dy, nout_ee_dy, nin_ee_vz, nout_ee_vz, nin_ee_data;  
  vector<vector<Double_t> > nin_mm_dy, nout_mm_dy, nin_mm_vz, nout_mm_vz, nin_mm_data;
  vector<vector<Double_t> > varin_ee_dy, varout_ee_dy, varin_ee_vz, varout_ee_vz;
  vector<vector<Double_t> > varin_mm_dy, varout_mm_dy, varin_mm_vz, varout_mm_vz;
  vector<vector<Double_t> > nin_em_data, nin_me_data;

  
  for(UInt_t jetIndex = 0; jetIndex < 3; ++jetIndex) {
    Double_t tmp_nin_kee_data=0, tmp_nin_kmm_data=0;
    vector<TH1F*> tmp_hNin_ree_mc,   tmp_hNout_ree_mc,   tmp_hNin_rmm_mc,   tmp_hNout_rmm_mc;
    vector<TH1F*> tmp_hNin_ree_data, tmp_hNout_ree_data, tmp_hNin_rmm_data, tmp_hNout_rmm_data;
    
    vector<Double_t> tmp_nin_ee_dy, tmp_nout_ee_dy, tmp_nin_ee_vz, tmp_nout_ee_vz, tmp_nin_ee_data;  
    vector<Double_t> tmp_nin_mm_dy, tmp_nout_mm_dy, tmp_nin_mm_vz, tmp_nout_mm_vz, tmp_nin_mm_data;
    vector<Double_t> tmp_varin_ee_dy, tmp_varout_ee_dy, tmp_varin_ee_vz, tmp_varout_ee_vz;
    vector<Double_t> tmp_varin_mm_dy, tmp_varout_mm_dy, tmp_varin_mm_vz, tmp_varout_mm_vz;
    vector<Double_t> tmp_nin_em_data, tmp_nin_me_data;

    char hname[50];
    for(Int_t imass=0; imass<nmass; imass++) {
      sprintf(hname,"hNin_%iJet_ree_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_ree_mc.push_back(new TH1F(hname,"",nbins,bins));  tmp_hNin_ree_mc[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_ree_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_ree_mc.push_back(new TH1F(hname,"",nbins,bins)); tmp_hNout_ree_mc[imass]->Sumw2();
      sprintf(hname,"hNin_%iJet_rmm_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_rmm_mc.push_back(new TH1F(hname,"",nbins,bins));  tmp_hNin_rmm_mc[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_rmm_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_rmm_mc.push_back(new TH1F(hname,"",nbins,bins)); tmp_hNout_rmm_mc[imass]->Sumw2();
      
      sprintf(hname,"hNin_%iJet_ree_data_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_ree_data.push_back(new TH1F(hname,"",nbins-1,bins));  tmp_hNin_ree_data[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_ree_data_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_ree_data.push_back(new TH1F(hname,"",nbins-1,bins)); tmp_hNout_ree_data[imass]->Sumw2();
      sprintf(hname,"hNin_%iJet_rmm_data_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_rmm_data.push_back(new TH1F(hname,"",nbins-1,bins));  tmp_hNin_rmm_data[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_rmm_data_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_rmm_data.push_back(new TH1F(hname,"",nbins-1,bins)); tmp_hNout_rmm_data[imass]->Sumw2();
      
      tmp_nin_ee_dy.push_back(0), tmp_nout_ee_dy.push_back(0), tmp_nin_ee_vz.push_back(0), tmp_nout_ee_vz.push_back(0), tmp_nin_ee_data.push_back(0);
      tmp_nin_mm_dy.push_back(0), tmp_nout_mm_dy.push_back(0), tmp_nin_mm_vz.push_back(0), tmp_nout_mm_vz.push_back(0), tmp_nin_mm_data.push_back(0);
      tmp_varin_ee_dy.push_back(0), tmp_varout_ee_dy.push_back(0), tmp_varin_ee_vz.push_back(0), tmp_varout_ee_vz.push_back(0);    
      tmp_varin_mm_dy.push_back(0), tmp_varout_mm_dy.push_back(0), tmp_varin_mm_vz.push_back(0), tmp_varout_mm_vz.push_back(0);
      tmp_nin_em_data.push_back(0), tmp_nin_me_data.push_back(0);
    }

    nin_kee_data.push_back(tmp_nin_kee_data);
    nin_kmm_data.push_back(tmp_nin_kmm_data);

    hNin_ree_mc.push_back(tmp_hNin_ree_mc); 
    hNout_ree_mc.push_back(tmp_hNout_ree_mc); 
    hNin_rmm_mc.push_back(tmp_hNin_rmm_mc); 
    hNout_rmm_mc.push_back(tmp_hNout_rmm_mc);
    hNin_ree_data.push_back(tmp_hNin_ree_data); 
    hNout_ree_data.push_back(tmp_hNout_ree_data); 
    hNin_rmm_data.push_back(tmp_hNin_rmm_data); 
    hNout_rmm_data.push_back(tmp_hNout_rmm_data);
    
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
    
    for(UInt_t ientry = 0; ientry <tree.tree_->GetEntries(); ientry++){
      tree.tree_->GetEntry(ientry);
        
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      if(!((tree.cuts_ & basic_selection) == basic_selection)) continue;
      if(tree.dilep_.M() < 12) continue;
        
      Int_t ijet = tree.njets_;
      if(ijet>2) ijet=2;
     
      if(tree.lep1_.Pt() < 20) continue;
      if(tree.lep2_.Pt() < 10) continue;

      Double_t weight = 1;
      if(tree.dstype_!=SmurfTree::data) {
        weight *= lumi*tree.scale1fb_*tree.sfWeightPU_*tree.sfWeightEff_*tree.sfWeightTrig_;
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
    
      if(tree.dstype_==SmurfTree::data && fabs(tree.dilep_.M() - mZ)<15) {
        if(finalState==kEleEle) nin_kee_data[ijet]++;
        if(finalState==kMuMu)	nin_kmm_data[ijet]++;
      }
      
      if(tree.jet1_.Pt()>15 && tree.dPhiDiLepJet1_>165.*TMath::Pi()/180.) continue;

      //For Z->mm/ee MC
      if(tree.dstype_==SmurfTree::dyee || tree.dstype_==SmurfTree::dymm) {
        

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
    
	    if(tree.lep1_.Pt() < _Pt1[imass]) continue;
	    if(tree.lep2_.Pt() < _Pt2[imass]) continue;
	
	    if(tree.dPhi_ > _DPhi[imass]*TMath::Pi()/180.) continue;

	    if(fabs(tree.dilep_.M() - mZ) < 15) {
	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dyee) {
                hNin_ree_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);                
               }
	      if(finalState==kMuMu   && tree.dstype_==SmurfTree::dymm) {
                hNin_rmm_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);              	 
              }
            } else if(tree.dilep_.M() < _Mll[imass]) {
	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dyee) {
                hNout_ree_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);	
              }
	      if(finalState==kMuMu   && tree.dstype_==SmurfTree::dymm) {
                hNout_rmm_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);
              }
            }
	
	    if(mt < _MtLow[imass] || mt > _MtHigh[imass]) continue;
	    if(minmet < 40) continue;
	
	    if(fabs(tree.dilep_.M() - mZ) < 15) {
	  
	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dyee) { 
	        nin_ee_dy[ijet][imass]+=weight/(Double_t)nmet; 
	        varin_ee_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet; 
	      }	  
	
	      if(finalState==kMuMu && tree.dstype_==SmurfTree::dymm) { 
	        nin_mm_dy[ijet][imass]+=weight/(Double_t)nmet; 
	        varin_mm_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
	      }	  
	 
            } else if(tree.dilep_.M() < _Mll[imass]) {
	  
	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dyee) { 
	        nout_ee_dy[ijet][imass]+=weight/(Double_t)nmet; 
	        varout_ee_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
	      }
	  
	      if(finalState==kMuMu && tree.dstype_==SmurfTree::dymm) { 
	        nout_mm_dy[ijet][imass]+=weight/(Double_t)nmet;
	        varout_mm_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
	      }
            }

//*************************************************************************************************
//*************************************************************************************************
//          For OLD LP2011 NTUPLES
// 	    if(fabs(tree.dilep_.M() - mZ) < 15) {
// 	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dymm) {
//                 hNin_ree_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);                
//                }
// 	      if(finalState==kMuMu   && tree.dstype_==SmurfTree::dyee) {
//                 hNin_rmm_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);              	 
//               }
//             } else if(tree.dilep_.M() < _Mll[imass]) {
// 	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dymm) {
//                 hNout_ree_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);	
//               }
// 	      if(finalState==kMuMu   && tree.dstype_==SmurfTree::dyee) {
//                 hNout_rmm_mc[ijet][imass]->Fill(minmet,weight/(Double_t)nmet);
//               }
//             }
	
// 	    if(mt < _MtLow[imass] || mt > _MtHigh[imass]) continue;
// 	    if(minmet < 40) continue;
	
// 	    if(fabs(tree.dilep_.M() - mZ) < 15) {
	  
// 	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dymm) { 
// 	        nin_ee_dy[ijet][imass]+=weight/(Double_t)nmet; 
// 	        varin_ee_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet; 
// 	      }	  
	
// 	      if(finalState==kMuMu && tree.dstype_==SmurfTree::dyee) { 
// 	        nin_mm_dy[ijet][imass]+=weight/(Double_t)nmet; 
// 	        varin_mm_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
// 	      }	  
	 
//             } else if(tree.dilep_.M() < _Mll[imass]) {
	  
// 	      if(finalState==kEleEle && tree.dstype_==SmurfTree::dymm) { 
// 	        nout_ee_dy[ijet][imass]+=weight/(Double_t)nmet; 
// 	        varout_ee_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
// 	      }
	  
// 	      if(finalState==kMuMu && tree.dstype_==SmurfTree::dyee) { 
// 	        nout_mm_dy[ijet][imass]+=weight/(Double_t)nmet;
// 	        varout_mm_dy[ijet][imass]+=weight*weight/(Double_t)nmet/(Double_t)nmet;
// 	      }
//             }
//*************************************************************************************************
//*************************************************************************************************



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
    
	  if(tree.lep1_.Pt() < _Pt1[imass]) continue;
	  if(tree.lep2_.Pt() < _Pt2[imass]) continue;
	
	  if(tree.dPhi_ > _DPhi[imass]*TMath::Pi()/180.) continue;

	  if(mt < _MtLow[imass] || mt > _MtHigh[imass]) continue;
	  if(minmet < 40) continue;
	
	  if(fabs(tree.dilep_.M() - mZ) < 15) {
	  
	    if(finalState==kEleEle && tree.dstype_==SmurfTree::data) { nin_ee_data[ijet][imass]++; }	  
	    if(finalState==kEleEle 
	       && (tree.dstype_==SmurfTree::wz || tree.dstype_==SmurfTree::zz)
	       && (tree.lep1MotherMcId_==23 && tree.lep2MotherMcId_==23)) {
	      nin_ee_vz[ijet][imass]+=weight;
	      varin_ee_vz[ijet][imass]+=weight*weight;
	    }
	
	    if(finalState==kMuMu && tree.dstype_==SmurfTree::data) { nin_mm_data[ijet][imass]++; }	  
	    if(finalState==kMuMu 
	       && (tree.dstype_==SmurfTree::wz || tree.dstype_==SmurfTree::zz)
	       && (tree.lep1MotherMcId_==23 && tree.lep2MotherMcId_==23)) {
	      nin_mm_vz[ijet][imass]+=weight;
	      varin_mm_vz[ijet][imass]+=weight*weight;
	    }	  
	    if(finalState==kMuEle && tree.dstype_==SmurfTree::data) nin_me_data[ijet][imass]++;           
	    if(finalState==kEleMu && tree.dstype_==SmurfTree::data) nin_em_data[ijet][imass]++; 

          } else if(tree.dilep_.M() < _Mll[imass]) {

	    if(finalState==kEleEle 
	       && (tree.dstype_==SmurfTree::wz || tree.dstype_==SmurfTree::zz)
	       && (tree.lep1MotherMcId_==23 && tree.lep2MotherMcId_==23)) {
	      nout_ee_vz[ijet][imass]+=weight;
	      varout_ee_vz[ijet][imass]+=weight*weight;
	    }
	  
	    if(finalState==kMuMu 
	       && (tree.dstype_==SmurfTree::wz || tree.dstype_==SmurfTree::zz)
	       && (tree.lep1MotherMcId_==23 && tree.lep2MotherMcId_==23)) {
	      nout_mm_vz[ijet][imass]+=weight;
	      varout_mm_vz[ijet][imass]+=weight*weight;
	    }
          }
        }    
      }
    }
  }

  
  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  vector<vector<Double_t> > DYBkgScaleFactorHiggsSelection;
  vector<vector<Double_t> > DYBkgScaleFactorHiggsSelectionErr;
  vector<Double_t> DYBkgScaleFactorWWPreselection;
  vector<Double_t> DYBkgScaleFactorWWPreselectionErr;



  for(UInt_t jetIndex = 0; jetIndex < 3; ++jetIndex) {
    vector<Double_t> tmpDYBkgScaleFactorHiggsSelection;
    vector<Double_t> tmpDYBkgScaleFactorHiggsSelectionErr;
    Double_t tmpDYBkgScaleFactorWWPreselection;
    Double_t tmpDYBkgScaleFactorWWPreselectionErr;

    cout << "************************************************************************\n";
    cout << jetIndex << "-Jet Bin DY Bkg Scale Factor Computation\n";
    cout << "************************************************************************\n";

    Double_t k    = sqrt(nin_kee_data[jetIndex]/nin_kmm_data[jetIndex]);
    Double_t kerr = 0.5*k*sqrt(1.0/nin_kee_data[jetIndex] + 1.0/nin_kmm_data[jetIndex]);
    char buffer[200];
  
    cout << "Electron to muon efficiency ratio is " << k << " +/- " << kerr << endl;

    cout << endl;
    cout << jetIndex << "-jet bin summary:" << endl;
    cout << setw(4) << "sel" << "   ";
    cout << setw(25) << "R_out/in       ";
    cout << setw(15) << "mm/me/em/ee";
    cout << setw(20) << "N_in (OF,VZ sub)";
    cout << setw(20) << "N_in (data/MC)";
    cout << setw(20) << "   N_out data   ";
    cout << setw(20) << "N_out MC  ";
    cout << setw(20) << "N_out (data/MC)" << endl;


    for(Int_t imass=0; imass<nmass; imass++) {
      cout << setw(4) << mH[imass] << "   ";
    
      //
      // compute Routin from MC
      //

      Bool_t passProtectionAgainstLargeStatisticalFluctuations = kFALSE;
      Int_t MetBinToComputeRoutin = nbins;

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

      //Protect against huge statistical uncertainty for Rout/in
      //If the stat uncertainty is > 65% , then take R from previous MET bin
      while (!passProtectionAgainstLargeStatisticalFluctuations) {
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
        ReeErrSyst = computeSyst(hNout_ree_mc[jetIndex][imass],hNin_ree_mc[jetIndex][imass], MetBinToComputeRoutin);    
    
        Rmm        = nout_mm/nin_mm;
        RmmErrStat = Rmm*sqrt(errin_mm*errin_mm/nin_mm/nin_mm + errout_mm*errout_mm/nout_mm/nout_mm);
        RmmErrSyst = computeSyst(hNout_rmm_mc[jetIndex][imass],hNin_rmm_mc[jetIndex][imass], MetBinToComputeRoutin);
    

        TH1F *hout = (TH1F*)hNout_ree_mc[jetIndex][imass]->Clone("hout");
        hout->Add(hNout_rmm_mc[jetIndex][imass]);
        TH1F *hin = (TH1F*)hNin_ree_mc[jetIndex][imass]->Clone("hin");
        hin->Add(hNin_rmm_mc[jetIndex][imass]);
        Rll        = nout_ll/nin_ll;
        RllErrStat = Rll*sqrt(errin_ll*errin_ll/nin_ll/nin_ll + errout_ll*errout_ll/nout_ll/nout_ll);
        RllErrSyst = computeSyst(hout,hin, MetBinToComputeRoutin);
        delete hout;
        delete hin;

        if (RllErrStat / Rll <= 0.65 || MetBinToComputeRoutin == nbins - 1) {
          passProtectionAgainstLargeStatisticalFluctuations = kTRUE;
        } else {
          MetBinToComputeRoutin -= 1;
          
        }
      }
      
//       cout << "\n" << nin_ee << " " << nout_ee << " : " << nin_mm << " " << nout_mm << " "
//            << endl;

      sprintf(buffer,"%.2f +/- %.2f +/- %.2f",Rll,RllErrStat,RllErrSyst);
//     sprintf(buffer,"%.3f +/- %.3f +/- %.3f",Ree,ReeErrStat,ReeErrSyst);
//     sprintf(buffer,"%.3f +/- %.3f +/- %.3f",Rmm,RmmErrStat,RmmErrSyst);
      cout << setw(25) << buffer; 

      //
      // raw in-yields in data
      //
      sprintf(buffer,"%i/%i/%i/%i",Int_t(nin_mm_data[jetIndex][imass]), Int_t(nin_me_data[jetIndex][imass]), Int_t(nin_em_data[jetIndex][imass]), Int_t(nin_ee_data[jetIndex][imass]));
      cout << setw(15) << buffer;
    
      //
      // in-yields in data after OF and VZ subtraction
      //
      Double_t nof = nin_em_data[jetIndex][imass] + nin_me_data[jetIndex][imass];
      Double_t nin_ee_sub = nin_ee_data[jetIndex][imass] - 0.5*k*nof - nin_ee_vz[jetIndex][imass];
      Double_t err_ee_sub = sqrt(nin_ee_data[jetIndex][imass] + 0.5*0.5*k*k*nof*nof*(kerr*kerr/k/k + 1.0/nof) + varin_ee_vz[jetIndex][imass] + pow(nin_ee_vz[jetIndex][imass]*vzNormSystematic,2) );
      Double_t nin_mm_sub = nin_mm_data[jetIndex][imass] - 0.5/k*nof - nin_mm_vz[jetIndex][imass];
      Double_t err_mm_sub = sqrt(nin_mm_data[jetIndex][imass] + 0.5*0.5/k/k*nof*nof*(kerr*kerr/k/k + 1.0/nof) + varin_mm_vz[jetIndex][imass] + pow(nin_mm_vz[jetIndex][imass]*vzNormSystematic,2));
      Double_t nin_ll_sub = nin_ee_sub + nin_mm_sub;
      Double_t err_ll_sub = sqrt(nin_mm_data[jetIndex][imass] + nin_ee_data[jetIndex][imass] 
                                 + 0.5*0.5*( (k+1.0/k)*(k+1.0/k)*nof*nof*kerr*kerr + (k+1.0/k)*(k+1.0/k)*nof ) 
                                 + varin_mm_vz[jetIndex][imass] + varin_ee_vz[jetIndex][imass] + pow((nin_ee_vz[jetIndex][imass]+nin_mm_vz[jetIndex][imass])*vzNormSystematic,2));
      if(nin_ee_sub <= 0) nin_ee_sub = 1;
      if(nin_mm_sub <= 0) nin_mm_sub = 1;
      if(nin_ll_sub <= 0) nin_ll_sub = 1;
    
      sprintf(buffer,"%.2f +/- %.2f",nin_ll_sub,err_ll_sub);
//    sprintf(buffer,"%.2f +/- %.2f",nin_ee_sub,err_ee_sub);
//    sprintf(buffer,"%.2f +/- %.2f",nin_mm_sub,err_mm_sub);
      cout << setw(20) << buffer;
    
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
      cout << setw(20) << buffer;
    
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
      cout << "   " << setw(20) << buffer;
    
      //
      // out-yield in MC
      //
      sprintf(buffer,"%.2f +/- %.2f",nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass],sqrt(varout_ee_dy[jetIndex][imass] + varout_mm_dy[jetIndex][imass]));
//    sprintf(buffer,"%.2f +/- %.2f",nout_ee_dy[jetIndex][imass],sqrt(varout_ee_dy[jetIndex][imass]));
//    sprintf(buffer,"%.2f +/- %.2f",nout_mm_dy[jetIndex][imass],sqrt(varout_mm_dy[jetIndex][imass]));
      cout << setw(20) << buffer;
    
      //
      // out-yield data/MC scale factor
      //
      Double_t sfout_ee     = nout_ee_pre/nout_ee_dy[jetIndex][imass];
      Double_t sfout_ee_err = sfout_ee*sqrt(nout_ee_err*nout_ee_err/nout_ee_pre/nout_ee_pre + varout_ee_dy[jetIndex][imass]/nout_ee_dy[jetIndex][imass]/nout_ee_dy[jetIndex][imass]);
      Double_t sfout_mm     = nout_mm_pre/nout_mm_dy[jetIndex][imass];
      Double_t sfout_mm_err = sfout_mm*sqrt(nout_mm_err*nout_mm_err/nout_mm_pre/nout_mm_pre + varout_mm_dy[jetIndex][imass]/nout_mm_dy[jetIndex][imass]/nout_mm_dy[jetIndex][imass]);
      Double_t sfout_ll     = nout_ll_pre/(nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass]);
      Double_t sfout_ll_err = sfout_ll*sqrt(nout_ll_err*nout_ll_err/nout_ll_pre/nout_ll_pre 
                                            + (varout_ee_dy[jetIndex][imass]+varout_mm_dy[jetIndex][imass])/(nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass])/(nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass]));
    
      sprintf(buffer,"%.2f +/- %.2f",sfout_ll,sfout_ll_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfout_ee,sfout_ee_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfout_mm,sfout_mm_err);
      cout << setw(20) << buffer;        
      cout << endl;
      if (imass == 0) {
        tmpDYBkgScaleFactorWWPreselection = sfout_ll;
        tmpDYBkgScaleFactorWWPreselectionErr = sfout_ll_err;
      } else {
        tmpDYBkgScaleFactorHiggsSelection.push_back(sfout_ll);
        tmpDYBkgScaleFactorHiggsSelectionErr.push_back(sfout_ll_err);
      }

    }

    DYBkgScaleFactorHiggsSelection.push_back(tmpDYBkgScaleFactorHiggsSelection);
    DYBkgScaleFactorHiggsSelectionErr.push_back(tmpDYBkgScaleFactorHiggsSelectionErr);
    DYBkgScaleFactorWWPreselection.push_back(tmpDYBkgScaleFactorWWPreselection);
    DYBkgScaleFactorWWPreselectionErr.push_back(tmpDYBkgScaleFactorWWPreselectionErr);

  }

  //***************************************************************************
  // Generate DY Scale Factor and Systematics Code for card creation
  //***************************************************************************
  ofstream outf("DYBkgScaleFactors.h");

  outf << "Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {" << endl;
  
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


  outf << "Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {" << endl;
  
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
  const Double_t pi = 3.14159265358979;
  Double_t dphi = acos(cos(lepPhi-metPhi));
  if(dphi > 0.5*pi)
    return met;
    
  return met*sin(dphi);
}

//--------------------------------------------------------------------------------------------------
Double_t computeSyst(const TH1F *hout, const TH1F *hin, Int_t binUsed)
{
  const Int_t nbins = binUsed;
  Double_t r0=(hout->GetBinContent(nbins))/(hin->GetBinContent(nbins));
  Double_t dr=0;
//cout << "R0 = " << r0 << endl;
  for(Int_t ibin=1; ibin<=nbins; ibin++) {
    Double_t r = (hout->GetBinContent(ibin))/(hin->GetBinContent(ibin));
    if(fabs(r-r0) > dr) dr = fabs(r-r0);
//cout << "Nin = " << hin->GetBinContent(ibin) << ", Nout = " << hout->GetBinContent(ibin) << ", R = " << r << ", dr = " << dr << endl;
  }  
  return dr;
}

TGraphErrors* MakeRoutinGraph(const TH1F *hout, const TH1F *hin, string graphname, vector<Double_t> bins) {

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

//     cerr << " done bin " << b << " " << x[b] << " : " << n1 << " +/-" << hout->GetBinError(b+1) << " / " << n2 << " +/- " << hin->GetBinError(b+1) << " = " << ratio << " " << err << endl;
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
