#ifndef _TEVTPROB_HH_
#define _TEVTPROB_HH_
//-----------------------------------------------------------------------------
// Description: Class TEvtProb: EvtProb base class
// ------------
//
//      Event Probability Density Calculation
//
// Feb 21 2011
// Sergo Jindariani
// Yanyan Gao
// Kevin Burkett
//-----------------------------------------------------------------------------
#include <sstream>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"

#include "TMCFM.hh"
#include "TVar.hh"
#include "TUtil.hh"
#include "PhaseSpace.hh"
#include "foam/TFoamIntegrand.hh"

#include "TH2F.h"
#include "TH1F.h"
#include "assert.h"
#include "TROOT.h"



//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject {

public:
  //--------------------
  // Variables
  //--------------------
  int _seed;
  int _ncalls;
  int _smearLevel;
  int _isApplyFake;

  TVar::Process _process;
  TVar::MatrixElement _matrixElement;
  TVar::HWWPhaseSpace _hwwPhaseSpace;
  EffHist _effhist;
  FRHist _FRhist;
  BoostHist _boosthist;

  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb();
 ~TEvtProb();

  //----------------------
  // Function
  //----------------------
  void SetSeed(int tmp) { _seed = tmp; }
  void SetNcalls(int tmp) { _ncalls= tmp; }
  void SetSmearLevel(int tmp) { _smearLevel=tmp; }
  void SetApplyFake (int tmp) { _isApplyFake=tmp; }
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void SetHWWPhaseSpace(TVar::HWWPhaseSpace tmp){ _hwwPhaseSpace = tmp; }

  void SetMCHist(int proc, TVar::VerbosityLevel verbosity) {
    if (verbosity >= TVar::DEBUG) std::cout << "TEvtProb::SetMCHist for process " << TVar::ProcessName(proc) << std::endl;
    TFile *fUtil = TFile::Open("Util.root", "READ");
    assert(fUtil);
    gROOT->cd();
    // lepton efficiencies are taken from the WW MC
    _effhist.els_eff_mc = (TH2F*) fUtil->Get("ww_heleEff")->Clone();
    _effhist.mus_eff_mc = (TH2F*) fUtil->Get("ww_hmuEff")->Clone();

    if (TVar::SmurfProcessName(proc)=="wjets") {
      _boosthist.kx = (TH1F*) fUtil->Get(TString("ww_kx"))->Clone();
      _boosthist.ky = (TH1F*) fUtil->Get(TString("ww_ky"))->Clone();
    } 
   else
     {
       _boosthist.kx = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_kx"))->Clone();
       _boosthist.ky = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_ky"))->Clone();
     }
    
    // fake-rate histograms
    _FRhist.els_fr = (TH2F*) fUtil->Get("wjets_heleFR")->Clone();
    _FRhist.mus_fr = (TH2F*) fUtil->Get("wjets_hmuFR")->Clone();
    _FRhist.els_part_fo = (TH2F*) fUtil->Get("wjets_heleGenFR")->Clone();
    _FRhist.mus_part_fo = (TH2F*) fUtil->Get("wjets_hmuGenFR")->Clone();
    _FRhist.els_ptres = (TH2F*) fUtil->Get("wjets_heleFOResponse")->Clone();
    _FRhist.mus_ptres = (TH2F*) fUtil->Get("wjets_hmuFOResponse")->Clone();
    
    fUtil->Close();
  }

 void SetFRHist() {
    TFile *elFRfile = TFile::Open("ww_el_fr_EG.root", "READ");
    assert(elFRfile);
    gROOT->cd();
    _FRhist.els_fr = (TH2F*) elFRfile->Get("el_fr_v4_wwV1")->Clone();
    elFRfile->Close();

    TFile *muFRfile = TFile::Open("ww_mu_fr_Mu.root", "READ");
    assert(muFRfile);
    gROOT->cd();
    _FRhist.mus_fr = (TH2F*) muFRfile->Get("mu_fr_fo_wwV1_10")->Clone();
    muFRfile->Close();
 }
  
  void LOXsec(double* Xsec,double* Err);
  void LOXsec_foam(double* Xsec,double* Err); 
 
  void NeutrinoIntegrate(TVar::Process proc,
                         cdf_event_type cdf_event,
                         double *Xsec,
                         double *XsecErr, TVar::VerbosityLevel verbosity);

  void SetHiggsMass(double mass){
    masses_mcfm_.hmass=mass;
    masses_mcfm_.hwidth=HiggsWidth(mass);
}

 
  
  ClassDef(TEvtProb,0);
};


//----------------------------------------
// Global Function
//----------------------------------------
double Integrand_LOXsec(double * arg, unsigned int dim, void * param);

class  Integrand_LOXsec_foam : public TFoamIntegrand{
 public:
   Double_t Density(Int_t ndim, Double_t* );
};


double Integrand_NeutrinoIntegration(double * r, unsigned int NDim, void * param, BoostHist boosthist);
  
#endif
