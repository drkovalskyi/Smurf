#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include "TRandom.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <iomanip>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include "Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


int    verboseLevel =   0;
const double sigmaB = 0.35;

//------------------------------------------------------------------------------
// PlotHiggsRes
//------------------------------------------------------------------------------
void ZllNormalization
(
 TString datInputFile    = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/data.root",
 TString bgdInputFile    = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/backgroundA.root",
 Double_t scaleFactorLum = 12.1,
 Int_t period = 2
 )
{

  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
  chdata->Add(datFile1);

  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;

  TString effPath      = "";
  TString fakePath     = "";
  TString puPath       = "";
  if	 (period == 0){ // Full2012-Summer12-V9-3500ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
  }
  else if(period == -1){ // Full2012-Summer12-V7
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-02_V3/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-02_V3/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X//auxiliar/puWeights_Summer12.root";
  }
  else if(period == 1){ // Full2012-Summer12-V9-5000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-06_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_5000ipb_71mb.root";
  }
  else if(period == 2){ // Full2012-Summer12-V9-12000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-07_HCP_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True.root";
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  //----------------------------------------------------------------------------
  // Lepton Efficiency Scale Factors, Fake Rates, 
  // Trigger Efficiencies
  // Pileup Weights
  //----------------------------------------------------------------------------
  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open(fakePath.Data());
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open(fakePath.Data());
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  LeptonScaleLookup trigLookup(effPath.Data());

  TFile *fPU2012File = TFile::Open(puPath.Data());
  TH1D *fhDPU2012 = (TH1D*)(fPU2012File->Get("puWeights"));
  assert(fhDPU2012);
  fhDPU2012->SetDirectory(0);
  delete fPU2012File;

  double massCut[2] = {82,100};
  double metCut[2] = {-30,30};
  //----------------------------------------------------------------------------

  double NZmmSelected_MC = 0;
  double NZmmSelected_DA = 0;
  double NZeeSelected_MC = 0;
  double NZeeSelected_DA = 0;

  TH1D *ZmmMass_MC    	 = new TH1D("ZmmMass_MC"      , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_MC_BB 	 = new TH1D("ZmmMass_MC_BB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_MC_EE 	 = new TH1D("ZmmMass_MC_EE"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_MC_EB 	 = new TH1D("ZmmMass_MC_EB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_MC    	 = new TH1D("ZeeMass_MC"      , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_MC_BB 	 = new TH1D("ZeeMass_MC_BB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_MC_EE 	 = new TH1D("ZeeMass_MC_EE"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_MC_EB 	 = new TH1D("ZeeMass_MC_EB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_DA    	 = new TH1D("ZmmMass_DA"      , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_DA_BB 	 = new TH1D("ZmmMass_DA_BB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_DA_EE 	 = new TH1D("ZmmMass_DA_EE"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_DA_EB 	 = new TH1D("ZmmMass_DA_EB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_DA    	 = new TH1D("ZeeMass_DA"      , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_DA_BB 	 = new TH1D("ZeeMass_DA_BB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_DA_EE 	 = new TH1D("ZeeMass_DA_EE"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_DA_EB 	 = new TH1D("ZeeMass_DA_EB"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_CorMC    = new TH1D("ZmmMass_CorMC"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_CorMC_BB = new TH1D("ZmmMass_CorMC_BB", ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_CorMC_EE = new TH1D("ZmmMass_CorMC_EE", ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZmmMass_CorMC_EB = new TH1D("ZmmMass_CorMC_EB", ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_CorMC    = new TH1D("ZeeMass_CorMC"   , ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_CorMC_BB = new TH1D("ZeeMass_CorMC_BB", ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_CorMC_EE = new TH1D("ZeeMass_CorMC_EE", ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);
  TH1D *ZeeMass_CorMC_EB = new TH1D("ZeeMass_CorMC_EB", ";Mass [GeV/c^{2}]; NEvents; ", 200, massCut[0],massCut[1]);

  TH1D *Zll0Met_MC_X         = new TH1D("Zll0Met_MC_X"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0Met_MC_Y         = new TH1D("Zll0Met_MC_Y"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0TrackMet_MC_X    = new TH1D("Zll0TrackMet_MC_X"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0TrackMet_MC_Y    = new TH1D("Zll0TrackMet_MC_Y"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0Met_DA_X         = new TH1D("Zll0Met_DA_X"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0Met_DA_Y         = new TH1D("Zll0Met_DA_Y"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0TrackMet_DA_X    = new TH1D("Zll0TrackMet_DA_X"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0TrackMet_DA_Y    = new TH1D("Zll0TrackMet_DA_Y"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0Met_CorMC_X      = new TH1D("Zll0Met_CorMC_X"      , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0Met_CorMC_Y      = new TH1D("Zll0Met_CorMC_Y"      , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0TrackMet_CorMC_X = new TH1D("Zll0TrackMet_CorMC_X" , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll0TrackMet_CorMC_Y = new TH1D("Zll0TrackMet_CorMC_Y" , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1Met_MC_X         = new TH1D("Zll1Met_MC_X"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1Met_MC_Y         = new TH1D("Zll1Met_MC_Y"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1TrackMet_MC_X    = new TH1D("Zll1TrackMet_MC_X"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1TrackMet_MC_Y    = new TH1D("Zll1TrackMet_MC_Y"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1Met_DA_X         = new TH1D("Zll1Met_DA_X"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1Met_DA_Y         = new TH1D("Zll1Met_DA_Y"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1TrackMet_DA_X    = new TH1D("Zll1TrackMet_DA_X"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1TrackMet_DA_Y    = new TH1D("Zll1TrackMet_DA_Y"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1Met_CorMC_X      = new TH1D("Zll1Met_CorMC_X"      , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1Met_CorMC_Y      = new TH1D("Zll1Met_CorMC_Y"      , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1TrackMet_CorMC_X = new TH1D("Zll1TrackMet_CorMC_X" , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll1TrackMet_CorMC_Y = new TH1D("Zll1TrackMet_CorMC_Y" , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2Met_MC_X         = new TH1D("Zll2Met_MC_X"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2Met_MC_Y         = new TH1D("Zll2Met_MC_Y"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2TrackMet_MC_X    = new TH1D("Zll2TrackMet_MC_X"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2TrackMet_MC_Y    = new TH1D("Zll2TrackMet_MC_Y"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2Met_DA_X         = new TH1D("Zll2Met_DA_X"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2Met_DA_Y         = new TH1D("Zll2Met_DA_Y"         , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2TrackMet_DA_X    = new TH1D("Zll2TrackMet_DA_X"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2TrackMet_DA_Y    = new TH1D("Zll2TrackMet_DA_Y"    , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2Met_CorMC_X      = new TH1D("Zll2Met_CorMC_X"      , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2Met_CorMC_Y      = new TH1D("Zll2Met_CorMC_Y"      , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2TrackMet_CorMC_X = new TH1D("Zll2TrackMet_CorMC_X" , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);
  TH1D *Zll2TrackMet_CorMC_Y = new TH1D("Zll2TrackMet_CorMC_Y" , ";Missing Energy [GeV]; NEvents; ", 200, metCut[0],metCut[1]);

  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          njets;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         met;
  Float_t         trackMet;
  Float_t         metPhi;
  Float_t         trackMetPhi;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;
  Int_t           lep1McId;
  Int_t           lep2McId;
  Int_t           processId;
  Float_t         sfWeightPU;
  Float_t         sfWeightEff; 
  Float_t         sfWeightTrig;
  Float_t         sfWeightFR;
  Float_t         npu;

  //*****************************************************************************
  // MC Loop
  //*****************************************************************************

  background->SetBranchAddress( "cuts"          , &cuts 	  );
  background->SetBranchAddress( "dstype"        , &dstype	  );
  background->SetBranchAddress( "nvtx"          , &nvtx 	  );
  background->SetBranchAddress( "njets"         , &njets	  );
  background->SetBranchAddress( "scale1fb"      , &scale1fb	  );
  background->SetBranchAddress( "lep1"          , &lep1 	  );
  background->SetBranchAddress( "lep2"          , &lep2 	  );
  background->SetBranchAddress( "dilep"         , &dilep	  );
  background->SetBranchAddress( "type"          , &type 	  );
  background->SetBranchAddress( "met"	     	, &met            );
  background->SetBranchAddress( "trackMet"   	, &trackMet       );
  background->SetBranchAddress( "metPhi"	, &metPhi         );
  background->SetBranchAddress( "trackMetPhi"	, &trackMetPhi    );
  background->SetBranchAddress( "lq1"           , &lq1  	  );
  background->SetBranchAddress( "lq2"           , &lq2  	  );
  background->SetBranchAddress( "lid1"          , &lid1 	  );
  background->SetBranchAddress( "lid2"          , &lid2 	  );
  background->SetBranchAddress( "lid3"          , &lid3 	  );
  background->SetBranchAddress( "lep1McId"	, &lep1McId	  );
  background->SetBranchAddress( "lep2McId"	, &lep2McId	  );
  background->SetBranchAddress( "processId"     , &processId	  );
  background->SetBranchAddress( "sfWeightPU"	, &sfWeightPU	  );
  background->SetBranchAddress( "sfWeightEff"	, &sfWeightEff    );
  background->SetBranchAddress( "sfWeightTrig"  , &sfWeightTrig   );
  background->SetBranchAddress( "sfWeightFR"	, &sfWeightFR	  );
  background->SetBranchAddress( "npu"           , &npu            );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    if ( lep1->pt() > 20 && lep2->pt() > 10 && dilep->mass() > massCut[0] && dilep->mass() < massCut[1]
         && (cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && lq1*lq2 < 0
	 && (cuts & SmurfTree::TopTag) != SmurfTree::TopTag
      ) {

      //----------------------------------------------------------------------------
      // Define event weights    
      //----------------------------------------------------------------------------
      double myWeight = 1.0;
      double add      = 1.0;

      //----------------------------------------------------------------------------
      // First classify event into tight+tight, tight+fail, fail+fail 
      //----------------------------------------------------------------------------
      int nFake = 0;
      if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(nFake < 0) assert(0);

      bool isRealLepton = false;
      if((TMath::Abs(lep1McId) == 11 || TMath::Abs(lep1McId) == 13) &&
	 (TMath::Abs(lep2McId) == 11 || TMath::Abs(lep2McId) == 13)) isRealLepton = true;
      double addLepEff = 1.0;
      double addFR     = 1.0;

      //----------------------------------------------------------------------------
      // Explicitly neglect fail+fail events
      //----------------------------------------------------------------------------
      if(nFake > 1){
	myWeight = 0.0;
      }
      //----------------------------------------------------------------------------
      // Tight+Fail events
      //----------------------------------------------------------------------------
      else if(nFake == 1){
	//----------------------------------------------------------------------------
	// For data, apply fake rate to generate bkg prediction
	//----------------------------------------------------------------------------
	if(dstype == SmurfTree::data){
          addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        							            (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	  addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        							            (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = addFR;
          myWeight	      = add;
	}
	//----------------------------------------------------------------------------
	// For real lepton or w+gamma:
	// apply fake rates, and give negative weight to subtract non jet-fake
	// contamination in the tight+fail sample
	//----------------------------------------------------------------------------
	else if(isRealLepton == true || dstype == SmurfTree::wgamma){
    	  addFR =       fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    		        						    (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	  addFR = addFR*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    									    (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

          add = addFR;
	  add = add*nPUScaleFactor2012(fhDPU2012,npu);

          addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
                      leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
          add = add*addLepEff;

    	  add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							    fabs(lep2->eta()), lep2->pt(), 
    							    TMath::Abs( lid1), TMath::Abs(lid2));
          myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
	}
	//----------------------------------------------------------------------------
	// Neglect fake lepton events in MC
	//----------------------------------------------------------------------------
	else {
          myWeight = 0.0;
	}
      }

      //----------------------------------------------------------------------------
      // Neglect any tight+tight events from data
      //----------------------------------------------------------------------------
      else if(dstype == SmurfTree::data) myWeight = 0.0;

      //----------------------------------------------------------------------------
      // DY->TauTau Embedded sample events: apply normalization weight.
      //----------------------------------------------------------------------------
      else if(dstype== SmurfTree::dyttDataDriven || dstype == SmurfTree::qcd) {
	myWeight = ZttScaleFactor(nvtx,period,scale1fb)*scaleFactorLum;
      }

      //----------------------------------------------------------------------------
      // Regular Tight+Tight events from Monte Carlo
      //----------------------------------------------------------------------------
      else if(dstype != SmurfTree::data){
	//----------------------------------------------------------------------------      
	// Apply lepton efficiency scale factors, trigger efficiencies,
	// Pileup weights
	//----------------------------------------------------------------------------
	add = 1.0;
	add = add*nPUScaleFactor2012(fhDPU2012,npu);

	addLepEff = leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1)*
        	    leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
	add = add*addLepEff;
	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
	myWeight = scale1fb*scaleFactorLum*add;
      }
      
      Double_t eventWeight = scaleFactorLum*scale1fb*sfWeightPU*sfWeightEff*sfWeightTrig*sfWeightFR;
      if(dstype == SmurfTree::data) eventWeight = sfWeightFR;

      if(fabs(eventWeight-myWeight)/eventWeight > 0.00001) printf("%f %f %d -> %f %f %f %f %f\n",eventWeight,myWeight,dstype,
                                                                  scale1fb,sfWeightPU,sfWeightEff,sfWeightTrig,sfWeightFR);

      double corr[2] = {1.0, 1.0};
      if     (TMath::Abs(lid1) == 13 && TMath::Abs(lep1->eta()) <  1.479){
        corr[0] = 0.99920 + gRandom->Gaus(0.00,0.010);
      }
      else if(TMath::Abs(lid1) == 13 && TMath::Abs(lep1->eta()) >= 1.479){
        corr[0] = 0.99934 + gRandom->Gaus(0.00,0.017);
      }
      else if(TMath::Abs(lid1) == 11 && TMath::Abs(lep1->eta()) <  1.479){
        corr[0] = 0.99807 + gRandom->Gaus(0.00,0.015);
      }
      else if(TMath::Abs(lid1) == 11 && TMath::Abs(lep1->eta()) >= 1.479){
        corr[0] = 0.99952 + gRandom->Gaus(0.00,0.030);
      }
      if     (TMath::Abs(lid2) == 13 && TMath::Abs(lep2->eta()) <  1.479){
        corr[1] = 0.99920 + gRandom->Gaus(0.00,0.010);
      }
      else if(TMath::Abs(lid2) == 13 && TMath::Abs(lep2->eta()) >= 1.479){
        corr[1] = 0.99934 + gRandom->Gaus(0.00,0.017);
      }
      else if(TMath::Abs(lid2) == 11 && TMath::Abs(lep2->eta()) <  1.479){
        corr[1] = 0.99807 + gRandom->Gaus(0.00,0.015);
      }
      else if(TMath::Abs(lid2) == 11 && TMath::Abs(lep2->eta()) >= 1.479){
        corr[1] = 0.99952 + gRandom->Gaus(0.00,0.030);
      }
      double pllx  = lep1->px()*corr[0]+lep2->px()*corr[1];
      double plly  = lep1->py()*corr[0]+lep2->py()*corr[1];
      double pllz  = lep1->pz()*corr[0]+lep2->pz()*corr[1];
      double ell   = lep1->E()*corr[0] +lep2->E() *corr[1];
      double dilmass = ell*ell -pllx*pllx -plly*plly -pllz*pllz;
      if(dilmass >=0) dilmass = sqrt(dilmass); else dilmass = 0.0;

      if ( type == SmurfTree::mm ) {        
        ZmmMass_MC->Fill(dilep->mass() , eventWeight);
	if     (TMath::Abs(lep1->eta()) < 1.479 && TMath::Abs(lep2->eta()) < 1.479) ZmmMass_MC_BB->Fill(dilep->mass() , eventWeight);
	else if(TMath::Abs(lep1->eta()) > 1.479 && TMath::Abs(lep2->eta()) > 1.479) ZmmMass_MC_EE->Fill(dilep->mass() , eventWeight);
	else                                                                        ZmmMass_MC_EB->Fill(dilep->mass() , eventWeight);
        ZmmMass_CorMC->Fill(dilmass , eventWeight);
	if     (TMath::Abs(lep1->eta()) < 1.479 && TMath::Abs(lep2->eta()) < 1.479) ZmmMass_CorMC_BB->Fill(dilmass , eventWeight);
	else if(TMath::Abs(lep1->eta()) > 1.479 && TMath::Abs(lep2->eta()) > 1.479) ZmmMass_CorMC_EE->Fill(dilmass , eventWeight);
	else                                                                        ZmmMass_CorMC_EB->Fill(dilmass , eventWeight);
        NZmmSelected_MC += eventWeight;
      }
      
      if ( type == SmurfTree::ee ) {        
        ZeeMass_MC->Fill(dilep->mass() , eventWeight);
	if     (TMath::Abs(lep1->eta()) < 1.479 && TMath::Abs(lep2->eta()) < 1.479) ZeeMass_MC_BB->Fill(dilep->mass() , eventWeight);
	else if(TMath::Abs(lep1->eta()) > 1.479 && TMath::Abs(lep2->eta()) > 1.479) ZeeMass_MC_EE->Fill(dilep->mass() , eventWeight);
	else                                                                        ZeeMass_MC_EB->Fill(dilep->mass() , eventWeight);
        ZeeMass_CorMC->Fill(dilmass , eventWeight);
	if     (TMath::Abs(lep1->eta()) < 1.479 && TMath::Abs(lep2->eta()) < 1.479) ZeeMass_CorMC_BB->Fill(dilmass , eventWeight);
	else if(TMath::Abs(lep1->eta()) > 1.479 && TMath::Abs(lep2->eta()) > 1.479) ZeeMass_CorMC_EE->Fill(dilmass , eventWeight);
	else                                                                        ZeeMass_CorMC_EB->Fill(dilmass , eventWeight);
        NZeeSelected_MC += eventWeight;
      }

      double metx=0.0;double mety=0.0;double trkmetx=0.0;double trkmety=0.0;
      if    (njets == 0){
        metx	= met*cos(metPhi)	   -0.45189+gRandom->Gaus(0.0,3.2);
        mety	= met*sin(metPhi)	   -0.20148+gRandom->Gaus(0.0,3.2);
        trkmetx = trackMet*cos(trackMetPhi)+0.12580+gRandom->Gaus(0.0,1.0);
        trkmety = trackMet*sin(trackMetPhi)+0.02615+gRandom->Gaus(0.0,1.0);
      }
      else if(njets == 1){
        metx	= met*cos(metPhi)	   -0.39040+gRandom->Gaus(0.0,3.6);
        mety	= met*sin(metPhi)	   -0.20427+gRandom->Gaus(0.0,3.6);
        trkmetx = trackMet*cos(trackMetPhi)+0.07639+gRandom->Gaus(0.0,4.5);
        trkmety = trackMet*sin(trackMetPhi)+0.01167+gRandom->Gaus(0.0,4.5);
      }
      else if(njets >= 2){
        metx	= met*cos(metPhi)	   -0.27127+gRandom->Gaus(0.0,4.3);
        mety	= met*sin(metPhi)	   -0.18935+gRandom->Gaus(0.0,4.3);
        trkmetx = trackMet*cos(trackMetPhi)+0.13328+gRandom->Gaus(0.0,6.0);
        trkmety = trackMet*sin(trackMetPhi)-0.01351+gRandom->Gaus(0.0,6.0);
      }
      if     (njets == 0){
        Zll0Met_MC_X     ->Fill(met*cos(metPhi)  	  , eventWeight);
	Zll0Met_MC_Y     ->Fill(met*sin(metPhi) 	  , eventWeight);
	Zll0TrackMet_MC_X->Fill(trackMet*cos(trackMetPhi) , eventWeight);
	Zll0TrackMet_MC_Y->Fill(trackMet*sin(trackMetPhi) , eventWeight);
        Zll0Met_CorMC_X     ->Fill(metx  	     	  , eventWeight);
	Zll0Met_CorMC_Y     ->Fill(mety 	     	  , eventWeight);
	Zll0TrackMet_CorMC_X->Fill(trkmetx 		  , eventWeight);
	Zll0TrackMet_CorMC_Y->Fill(trkmety 		  , eventWeight);
      } 
      else if(njets == 1){
        Zll1Met_MC_X     ->Fill(met*cos(metPhi)  	  , eventWeight);
	Zll1Met_MC_Y     ->Fill(met*sin(metPhi) 	  , eventWeight);
	Zll1TrackMet_MC_X->Fill(trackMet*cos(trackMetPhi) , eventWeight);
	Zll1TrackMet_MC_Y->Fill(trackMet*sin(trackMetPhi) , eventWeight);
        Zll1Met_CorMC_X     ->Fill(metx  	     	  , eventWeight);
	Zll1Met_CorMC_Y     ->Fill(mety 	     	  , eventWeight);
	Zll1TrackMet_CorMC_X->Fill(trkmetx 		  , eventWeight);
	Zll1TrackMet_CorMC_Y->Fill(trkmety 		  , eventWeight);
      } 
      else if(njets >= 2){
        Zll2Met_MC_X     ->Fill(met*cos(metPhi)  	  , eventWeight);
	Zll2Met_MC_Y     ->Fill(met*sin(metPhi) 	  , eventWeight);
	Zll2TrackMet_MC_X->Fill(trackMet*cos(trackMetPhi) , eventWeight);
	Zll2TrackMet_MC_Y->Fill(trackMet*sin(trackMetPhi) , eventWeight);
        Zll2Met_CorMC_X     ->Fill(metx  	     	  , eventWeight);
	Zll2Met_CorMC_Y     ->Fill(mety 	     	  , eventWeight);
	Zll2TrackMet_CorMC_X->Fill(trkmetx 		  , eventWeight);
	Zll2TrackMet_CorMC_Y->Fill(trkmety 		  , eventWeight);
      } 

    }
  }
  printf("--- Finished Bgdnal loop\n");

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************

  data->SetBranchAddress( "cuts"         , &cuts         );
  data->SetBranchAddress( "nvtx"         , &nvtx         );
  data->SetBranchAddress( "njets"        , &njets        );
  data->SetBranchAddress( "scale1fb"     , &scale1fb     );
  data->SetBranchAddress( "lep1"         , &lep1         );
  data->SetBranchAddress( "lep2"         , &lep2         );
  data->SetBranchAddress( "dilep"        , &dilep        );
  data->SetBranchAddress( "type"         , &type         );
  data->SetBranchAddress( "met" 	 , &met 	 );
  data->SetBranchAddress( "trackMet"	 , &trackMet	 );
  data->SetBranchAddress( "metPhi"	 , &metPhi	 );
  data->SetBranchAddress( "trackMetPhi"  , &trackMetPhi  );
  data->SetBranchAddress( "lq1"          , &lq1          );
  data->SetBranchAddress( "lq2"          , &lq2          );
  data->SetBranchAddress( "lid1"         , &lid1         );
  data->SetBranchAddress( "lid2"         , &lid2         );
  data->SetBranchAddress( "lid3"         , &lid3         );
  data->SetBranchAddress( "processId"    , &processId    );

  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    if ( lep1->pt() > 20 && lep2->pt() > 10 && dilep->mass() > massCut[0] && dilep->mass() < massCut[1]
         && (cuts & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && lq1*lq2 < 0
	 && (cuts & SmurfTree::TopTag) != SmurfTree::TopTag
         && ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
         && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
     ) {

      Double_t eventWeight = 1.0;

      if ( type == SmurfTree::mm ) {        
        ZmmMass_DA->Fill(dilep->mass() , eventWeight);
	if     (TMath::Abs(lep1->eta()) < 1.479 && TMath::Abs(lep2->eta()) < 1.479) ZmmMass_DA_BB->Fill(dilep->mass() , eventWeight);
	else if(TMath::Abs(lep1->eta()) > 1.479 && TMath::Abs(lep2->eta()) > 1.479) ZmmMass_DA_EE->Fill(dilep->mass() , eventWeight);
	else                                                                        ZmmMass_DA_EB->Fill(dilep->mass() , eventWeight);
        NZmmSelected_DA += eventWeight;
      }
      
      if ( type == SmurfTree::ee ) {        
        ZeeMass_DA->Fill(dilep->mass() , eventWeight);
	if     (TMath::Abs(lep1->eta()) < 1.479 && TMath::Abs(lep2->eta()) < 1.479) ZeeMass_DA_BB->Fill(dilep->mass() , eventWeight);
	else if(TMath::Abs(lep1->eta()) > 1.479 && TMath::Abs(lep2->eta()) > 1.479) ZeeMass_DA_EE->Fill(dilep->mass() , eventWeight);
	else                                                                        ZeeMass_DA_EB->Fill(dilep->mass() , eventWeight);
        NZeeSelected_DA += eventWeight;
      }

      if     (njets == 0){
        Zll0Met_DA_X     ->Fill(met*cos(metPhi) , eventWeight); 
	Zll0Met_DA_Y     ->Fill(met*sin(metPhi) , eventWeight);
	Zll0TrackMet_DA_X->Fill(trackMet*cos(trackMetPhi) , eventWeight);
	Zll0TrackMet_DA_Y->Fill(trackMet*sin(trackMetPhi) , eventWeight);
      } 
      else if(njets == 1){
        Zll1Met_DA_X     ->Fill(met*cos(metPhi) , eventWeight); 
	Zll1Met_DA_Y     ->Fill(met*sin(metPhi) , eventWeight);
	Zll1TrackMet_DA_X->Fill(trackMet*cos(trackMetPhi) , eventWeight);
	Zll1TrackMet_DA_Y->Fill(trackMet*sin(trackMetPhi) , eventWeight);
      } 
      else if(njets >= 2){
        Zll2Met_DA_X     ->Fill(met*cos(metPhi) , eventWeight); 
	Zll2Met_DA_Y     ->Fill(met*sin(metPhi) , eventWeight);
	Zll2TrackMet_DA_X->Fill(trackMet*cos(trackMetPhi) , eventWeight);
	Zll2TrackMet_DA_Y->Fill(trackMet*sin(trackMetPhi) , eventWeight);
      } 

    }
  }

  ZmmMass_MC      ->Scale(ZmmMass_DA   ->GetSumOfWeights()/ZmmMass_MC	->GetSumOfWeights());
  ZmmMass_MC_BB   ->Scale(ZmmMass_DA_BB->GetSumOfWeights()/ZmmMass_MC_BB->GetSumOfWeights());
  ZmmMass_MC_EE   ->Scale(ZmmMass_DA_EE->GetSumOfWeights()/ZmmMass_MC_EE->GetSumOfWeights());
  ZmmMass_MC_EB   ->Scale(ZmmMass_DA_EB->GetSumOfWeights()/ZmmMass_MC_EB->GetSumOfWeights());
  ZeeMass_MC      ->Scale(ZeeMass_DA   ->GetSumOfWeights()/ZeeMass_MC	->GetSumOfWeights());
  ZeeMass_MC_BB   ->Scale(ZeeMass_DA_BB->GetSumOfWeights()/ZeeMass_MC_BB->GetSumOfWeights());
  ZeeMass_MC_EE   ->Scale(ZeeMass_DA_EE->GetSumOfWeights()/ZeeMass_MC_EE->GetSumOfWeights());
  ZeeMass_MC_EB   ->Scale(ZeeMass_DA_EB->GetSumOfWeights()/ZeeMass_MC_EB->GetSumOfWeights());

  ZmmMass_CorMC   ->Scale(ZmmMass_DA   ->GetSumOfWeights()/ZmmMass_MC	->GetSumOfWeights());
  ZmmMass_CorMC_BB->Scale(ZmmMass_DA_BB->GetSumOfWeights()/ZmmMass_MC_BB->GetSumOfWeights());
  ZmmMass_CorMC_EE->Scale(ZmmMass_DA_EE->GetSumOfWeights()/ZmmMass_MC_EE->GetSumOfWeights());
  ZmmMass_CorMC_EB->Scale(ZmmMass_DA_EB->GetSumOfWeights()/ZmmMass_MC_EB->GetSumOfWeights());
  ZeeMass_CorMC   ->Scale(ZeeMass_DA   ->GetSumOfWeights()/ZeeMass_MC	->GetSumOfWeights());
  ZeeMass_CorMC_BB->Scale(ZeeMass_DA_BB->GetSumOfWeights()/ZeeMass_MC_BB->GetSumOfWeights());
  ZeeMass_CorMC_EE->Scale(ZeeMass_DA_EE->GetSumOfWeights()/ZeeMass_MC_EE->GetSumOfWeights());
  ZeeMass_CorMC_EB->Scale(ZeeMass_DA_EB->GetSumOfWeights()/ZeeMass_MC_EB->GetSumOfWeights());

  cout << "Zmm -->  Data/MC: " << NZmmSelected_DA << " - " << NZmmSelected_MC << endl;
  printf("MeanMC: %8.5f  %8.5f  %8.5f  %8.5f MeanDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZmmMass_MC->GetMean(),ZmmMass_MC_BB->GetMean(),ZmmMass_MC_EE->GetMean(),ZmmMass_MC_EB->GetMean(),
	 ZmmMass_DA->GetMean(),ZmmMass_DA_BB->GetMean(),ZmmMass_DA_EE->GetMean(),ZmmMass_DA_EB->GetMean());
  printf("RMSMC: %8.5f  %8.5f  %8.5f  %8.5f RMSDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZmmMass_MC->GetRMS(),ZmmMass_MC_BB->GetRMS(),ZmmMass_MC_EE->GetRMS(),ZmmMass_MC_EB->GetRMS(),
	 ZmmMass_DA->GetRMS(),ZmmMass_DA_BB->GetRMS(),ZmmMass_DA_EE->GetRMS(),ZmmMass_DA_EB->GetRMS());
  cout << "Zee -->  Data/MC: " << NZeeSelected_DA << " - " << NZeeSelected_MC << endl;
  printf("MeanMC: %8.5f  %8.5f  %8.5f  %8.5f MeanDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZeeMass_MC->GetMean(),ZeeMass_MC_BB->GetMean(),ZeeMass_MC_EE->GetMean(),ZeeMass_MC_EB->GetMean(),
	 ZeeMass_DA->GetMean(),ZeeMass_DA_BB->GetMean(),ZeeMass_DA_EE->GetMean(),ZeeMass_DA_EB->GetMean());
  printf("RMSMC: %8.5f  %8.5f  %8.5f  %8.5f RMSDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZeeMass_MC->GetRMS(),ZeeMass_MC_BB->GetRMS(),ZeeMass_MC_EE->GetRMS(),ZeeMass_MC_EB->GetRMS(),
	 ZeeMass_DA->GetRMS(),ZeeMass_DA_BB->GetRMS(),ZeeMass_DA_EE->GetRMS(),ZeeMass_DA_EB->GetRMS());

  double bias[4],smear[4];
  bias[0]  = ZmmMass_DA_BB->GetMean()/ZmmMass_MC_BB->GetMean();
  bias[1]  = ZmmMass_DA_EE->GetMean()/ZmmMass_MC_EE->GetMean();
  bias[2]  = ZeeMass_DA_BB->GetMean()/ZeeMass_MC_BB->GetMean();
  bias[3]  = ZeeMass_DA_EE->GetMean()/ZeeMass_MC_EE->GetMean();
  smear[0] = sqrt(TMath::Max(ZmmMass_DA_BB->GetRMS()*ZmmMass_DA_BB->GetRMS()-ZmmMass_MC_BB->GetRMS()*ZmmMass_MC_BB->GetRMS(),0.))/ZmmMass_DA_BB->GetMean();
  smear[1] = sqrt(TMath::Max(ZmmMass_DA_EE->GetRMS()*ZmmMass_DA_EE->GetRMS()-ZmmMass_MC_EE->GetRMS()*ZmmMass_MC_EE->GetRMS(),0.))/ZmmMass_DA_EE->GetMean();
  smear[2] = sqrt(TMath::Max(ZeeMass_DA_BB->GetRMS()*ZeeMass_DA_BB->GetRMS()-ZeeMass_MC_BB->GetRMS()*ZeeMass_MC_BB->GetRMS(),0.))/ZeeMass_DA_BB->GetMean();
  smear[3] = sqrt(TMath::Max(ZeeMass_DA_EE->GetRMS()*ZeeMass_DA_EE->GetRMS()-ZeeMass_MC_EE->GetRMS()*ZeeMass_MC_EE->GetRMS(),0.))/ZeeMass_DA_EE->GetMean();
  printf("bias:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smear: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);

  printf("******************************\n");
  printf("MeanMC: %8.5f  %8.5f  %8.5f  %8.5f MeanDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZmmMass_CorMC->GetMean(),ZmmMass_CorMC_BB->GetMean(),ZmmMass_CorMC_EE->GetMean(),ZmmMass_CorMC_EB->GetMean(),
	 ZmmMass_DA->GetMean(),ZmmMass_DA_BB->GetMean(),ZmmMass_DA_EE->GetMean(),ZmmMass_DA_EB->GetMean());
  printf("RMSMC: %8.5f  %8.5f  %8.5f  %8.5f RMSDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZmmMass_CorMC->GetRMS(),ZmmMass_CorMC_BB->GetRMS(),ZmmMass_CorMC_EE->GetRMS(),ZmmMass_CorMC_EB->GetRMS(),
	 ZmmMass_DA->GetRMS(),ZmmMass_DA_BB->GetRMS(),ZmmMass_DA_EE->GetRMS(),ZmmMass_DA_EB->GetRMS());
  printf("MeanMC: %8.5f  %8.5f  %8.5f  %8.5f MeanDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZeeMass_CorMC->GetMean(),ZeeMass_CorMC_BB->GetMean(),ZeeMass_CorMC_EE->GetMean(),ZeeMass_CorMC_EB->GetMean(),
	 ZeeMass_DA->GetMean(),ZeeMass_DA_BB->GetMean(),ZeeMass_DA_EE->GetMean(),ZeeMass_DA_EB->GetMean());
  printf("RMSMC: %8.5f  %8.5f  %8.5f  %8.5f RMSDA: %8.5f  %8.5f  %8.5f  %8.5f\n",
         ZeeMass_CorMC->GetRMS(),ZeeMass_CorMC_BB->GetRMS(),ZeeMass_CorMC_EE->GetRMS(),ZeeMass_CorMC_EB->GetRMS(),
	 ZeeMass_DA->GetRMS(),ZeeMass_DA_BB->GetRMS(),ZeeMass_DA_EE->GetRMS(),ZeeMass_DA_EB->GetRMS());
  bias[0]  = ZmmMass_DA_BB->GetMean()/ZmmMass_CorMC_BB->GetMean();
  bias[1]  = ZmmMass_DA_EE->GetMean()/ZmmMass_CorMC_EE->GetMean();
  bias[2]  = ZeeMass_DA_BB->GetMean()/ZeeMass_CorMC_BB->GetMean();
  bias[3]  = ZeeMass_DA_EE->GetMean()/ZeeMass_CorMC_EE->GetMean();
  smear[0] = sqrt(TMath::Max(ZmmMass_DA_BB->GetRMS()*ZmmMass_DA_BB->GetRMS()-ZmmMass_CorMC_BB->GetRMS()*ZmmMass_CorMC_BB->GetRMS(),0.))/ZmmMass_DA_BB->GetMean();
  smear[1] = sqrt(TMath::Max(ZmmMass_DA_EE->GetRMS()*ZmmMass_DA_EE->GetRMS()-ZmmMass_CorMC_EE->GetRMS()*ZmmMass_CorMC_EE->GetRMS(),0.))/ZmmMass_DA_EE->GetMean();
  smear[2] = sqrt(TMath::Max(ZeeMass_DA_BB->GetRMS()*ZeeMass_DA_BB->GetRMS()-ZeeMass_CorMC_BB->GetRMS()*ZeeMass_CorMC_BB->GetRMS(),0.))/ZeeMass_DA_BB->GetMean();
  smear[3] = sqrt(TMath::Max(ZeeMass_DA_EE->GetRMS()*ZeeMass_DA_EE->GetRMS()-ZeeMass_CorMC_EE->GetRMS()*ZeeMass_CorMC_EE->GetRMS(),0.))/ZeeMass_DA_EE->GetMean();
  printf("bias:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smear: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);

  Zll0Met_MC_X     ->Scale(Zll0Met_DA_X     ->GetSumOfWeights()/Zll0Met_MC_X	 ->GetSumOfWeights());
  Zll0Met_MC_Y     ->Scale(Zll0Met_DA_Y     ->GetSumOfWeights()/Zll0Met_MC_Y	 ->GetSumOfWeights());
  Zll0TrackMet_MC_X->Scale(Zll0TrackMet_DA_X->GetSumOfWeights()/Zll0TrackMet_MC_X->GetSumOfWeights());
  Zll0TrackMet_MC_Y->Scale(Zll0TrackMet_DA_Y->GetSumOfWeights()/Zll0TrackMet_MC_Y->GetSumOfWeights());
  Zll1Met_MC_X     ->Scale(Zll1Met_DA_X     ->GetSumOfWeights()/Zll1Met_MC_X	 ->GetSumOfWeights());
  Zll1Met_MC_Y     ->Scale(Zll1Met_DA_Y     ->GetSumOfWeights()/Zll1Met_MC_Y	 ->GetSumOfWeights());
  Zll1TrackMet_MC_X->Scale(Zll1TrackMet_DA_X->GetSumOfWeights()/Zll1TrackMet_MC_X->GetSumOfWeights());
  Zll1TrackMet_MC_Y->Scale(Zll1TrackMet_DA_Y->GetSumOfWeights()/Zll1TrackMet_MC_Y->GetSumOfWeights());
  Zll2Met_MC_X     ->Scale(Zll2Met_DA_X     ->GetSumOfWeights()/Zll2Met_MC_X	 ->GetSumOfWeights());
  Zll2Met_MC_Y     ->Scale(Zll2Met_DA_Y     ->GetSumOfWeights()/Zll2Met_MC_Y	 ->GetSumOfWeights());
  Zll2TrackMet_MC_X->Scale(Zll2TrackMet_DA_X->GetSumOfWeights()/Zll2TrackMet_MC_X->GetSumOfWeights());
  Zll2TrackMet_MC_Y->Scale(Zll2TrackMet_DA_Y->GetSumOfWeights()/Zll2TrackMet_MC_Y->GetSumOfWeights());
  Zll0Met_CorMC_X     ->Scale(Zll0Met_DA_X     ->GetSumOfWeights()/Zll0Met_MC_X	    ->GetSumOfWeights());
  Zll0Met_CorMC_Y     ->Scale(Zll0Met_DA_Y     ->GetSumOfWeights()/Zll0Met_MC_Y	    ->GetSumOfWeights());
  Zll0TrackMet_CorMC_X->Scale(Zll0TrackMet_DA_X->GetSumOfWeights()/Zll0TrackMet_MC_X->GetSumOfWeights());
  Zll0TrackMet_CorMC_Y->Scale(Zll0TrackMet_DA_Y->GetSumOfWeights()/Zll0TrackMet_MC_Y->GetSumOfWeights());
  Zll1Met_CorMC_X     ->Scale(Zll1Met_DA_X     ->GetSumOfWeights()/Zll1Met_MC_X	    ->GetSumOfWeights());
  Zll1Met_CorMC_Y     ->Scale(Zll1Met_DA_Y     ->GetSumOfWeights()/Zll1Met_MC_Y	    ->GetSumOfWeights());
  Zll1TrackMet_CorMC_X->Scale(Zll1TrackMet_DA_X->GetSumOfWeights()/Zll1TrackMet_MC_X->GetSumOfWeights());
  Zll1TrackMet_CorMC_Y->Scale(Zll1TrackMet_DA_Y->GetSumOfWeights()/Zll1TrackMet_MC_Y->GetSumOfWeights());
  Zll2Met_CorMC_X     ->Scale(Zll2Met_DA_X     ->GetSumOfWeights()/Zll2Met_MC_X	    ->GetSumOfWeights());
  Zll2Met_CorMC_Y     ->Scale(Zll2Met_DA_Y     ->GetSumOfWeights()/Zll2Met_MC_Y	    ->GetSumOfWeights());
  Zll2TrackMet_CorMC_X->Scale(Zll2TrackMet_DA_X->GetSumOfWeights()/Zll2TrackMet_MC_X->GetSumOfWeights());
  Zll2TrackMet_CorMC_Y->Scale(Zll2TrackMet_DA_Y->GetSumOfWeights()/Zll2TrackMet_MC_Y->GetSumOfWeights());

  printf("********MET******************\n");
  printf("   MeanMC: %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f\n",
         Zll0Met_MC_X->GetMean(),Zll0Met_MC_Y->GetMean(),Zll0TrackMet_MC_X->GetMean(),Zll0TrackMet_MC_Y->GetMean(),
         Zll1Met_MC_X->GetMean(),Zll1Met_MC_Y->GetMean(),Zll1TrackMet_MC_X->GetMean(),Zll1TrackMet_MC_Y->GetMean(),
         Zll2Met_MC_X->GetMean(),Zll2Met_MC_Y->GetMean(),Zll2TrackMet_MC_X->GetMean(),Zll2TrackMet_MC_Y->GetMean());
  printf("   MeanDA: %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f\n",
         Zll0Met_DA_X->GetMean(),Zll0Met_DA_Y->GetMean(),Zll0TrackMet_DA_X->GetMean(),Zll0TrackMet_DA_Y->GetMean(),
         Zll1Met_DA_X->GetMean(),Zll1Met_DA_Y->GetMean(),Zll1TrackMet_DA_X->GetMean(),Zll1TrackMet_DA_Y->GetMean(),
         Zll2Met_DA_X->GetMean(),Zll2Met_DA_Y->GetMean(),Zll2TrackMet_DA_X->GetMean(),Zll2TrackMet_DA_Y->GetMean());
  printf("MeanCorMC: %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f\n",
         Zll0Met_CorMC_X->GetMean(),Zll0Met_CorMC_Y->GetMean(),Zll0TrackMet_CorMC_X->GetMean(),Zll0TrackMet_CorMC_Y->GetMean(),
         Zll1Met_CorMC_X->GetMean(),Zll1Met_CorMC_Y->GetMean(),Zll1TrackMet_CorMC_X->GetMean(),Zll1TrackMet_CorMC_Y->GetMean(),
         Zll2Met_CorMC_X->GetMean(),Zll2Met_CorMC_Y->GetMean(),Zll2TrackMet_CorMC_X->GetMean(),Zll2TrackMet_CorMC_Y->GetMean());
  printf("   RMSMC: %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f\n",
         Zll0Met_MC_X->GetRMS(),Zll0Met_MC_Y->GetRMS(),Zll0TrackMet_MC_X->GetRMS(),Zll0TrackMet_MC_Y->GetRMS(),
         Zll1Met_MC_X->GetRMS(),Zll1Met_MC_Y->GetRMS(),Zll1TrackMet_MC_X->GetRMS(),Zll1TrackMet_MC_Y->GetRMS(),
         Zll2Met_MC_X->GetRMS(),Zll2Met_MC_Y->GetRMS(),Zll2TrackMet_MC_X->GetRMS(),Zll2TrackMet_MC_Y->GetRMS());
  printf("   RMSDA: %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f\n",
         Zll0Met_DA_X->GetRMS(),Zll0Met_DA_Y->GetRMS(),Zll0TrackMet_DA_X->GetRMS(),Zll0TrackMet_DA_Y->GetRMS(),
         Zll1Met_DA_X->GetRMS(),Zll1Met_DA_Y->GetRMS(),Zll1TrackMet_DA_X->GetRMS(),Zll1TrackMet_DA_Y->GetRMS(),
         Zll2Met_DA_X->GetRMS(),Zll2Met_DA_Y->GetRMS(),Zll2TrackMet_DA_X->GetRMS(),Zll2TrackMet_DA_Y->GetRMS());
  printf("RMSCorMC: %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f | %8.5f  %8.5f  %8.5f  %8.5f\n",
         Zll0Met_CorMC_X->GetRMS(),Zll0Met_CorMC_Y->GetRMS(),Zll0TrackMet_CorMC_X->GetRMS(),Zll0TrackMet_CorMC_Y->GetRMS(),
         Zll1Met_CorMC_X->GetRMS(),Zll1Met_CorMC_Y->GetRMS(),Zll1TrackMet_CorMC_X->GetRMS(),Zll1TrackMet_CorMC_Y->GetRMS(),
         Zll2Met_CorMC_X->GetRMS(),Zll2Met_CorMC_Y->GetRMS(),Zll2TrackMet_CorMC_X->GetRMS(),Zll2TrackMet_CorMC_Y->GetRMS());
  bias[0]  = Zll0Met_DA_X     ->GetMean()-Zll0Met_MC_X     ->GetMean();
  bias[1]  = Zll0Met_DA_Y     ->GetMean()-Zll0Met_MC_Y     ->GetMean();
  bias[2]  = Zll0TrackMet_DA_X->GetMean()-Zll0TrackMet_MC_X->GetMean();
  bias[3]  = Zll0TrackMet_DA_Y->GetMean()-Zll0TrackMet_MC_Y->GetMean();
  smear[0] = sqrt(TMath::Max(Zll0Met_DA_X     ->GetRMS()*Zll0Met_DA_X	  ->GetRMS()-Zll0Met_MC_X     ->GetRMS()*Zll0Met_MC_X	  ->GetRMS(),0.));
  smear[1] = sqrt(TMath::Max(Zll0Met_DA_Y     ->GetRMS()*Zll0Met_DA_Y	  ->GetRMS()-Zll0Met_MC_Y     ->GetRMS()*Zll0Met_MC_Y	  ->GetRMS(),0.));
  smear[2] = sqrt(TMath::Max(Zll0TrackMet_DA_X->GetRMS()*Zll0TrackMet_DA_X->GetRMS()-Zll0TrackMet_MC_X->GetRMS()*Zll0TrackMet_MC_X->GetRMS(),0.));
  smear[3] = sqrt(TMath::Max(Zll0TrackMet_DA_Y->GetRMS()*Zll0TrackMet_DA_Y->GetRMS()-Zll0TrackMet_MC_Y->GetRMS()*Zll0TrackMet_MC_Y->GetRMS(),0.));
  printf(" bias0:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smear0: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);
  bias[0]  = Zll0Met_DA_X     ->GetMean()-Zll0Met_CorMC_X     ->GetMean();
  bias[1]  = Zll0Met_DA_Y     ->GetMean()-Zll0Met_CorMC_Y     ->GetMean();
  bias[2]  = Zll0TrackMet_DA_X->GetMean()-Zll0TrackMet_CorMC_X->GetMean();
  bias[3]  = Zll0TrackMet_DA_Y->GetMean()-Zll0TrackMet_CorMC_Y->GetMean();
  smear[0] = sqrt(TMath::Max(Zll0Met_CorMC_X     ->GetRMS()*Zll0Met_CorMC_X     ->GetRMS()-Zll0Met_DA_X     ->GetRMS()*Zll0Met_DA_X     ->GetRMS(),0.));
  smear[1] = sqrt(TMath::Max(Zll0Met_CorMC_Y     ->GetRMS()*Zll0Met_CorMC_Y     ->GetRMS()-Zll0Met_DA_Y     ->GetRMS()*Zll0Met_DA_Y     ->GetRMS(),0.));
  smear[2] = sqrt(TMath::Max(Zll0TrackMet_CorMC_X->GetRMS()*Zll0TrackMet_CorMC_X->GetRMS()-Zll0TrackMet_DA_X->GetRMS()*Zll0TrackMet_DA_X->GetRMS(),0.));
  smear[3] = sqrt(TMath::Max(Zll0TrackMet_CorMC_Y->GetRMS()*Zll0TrackMet_CorMC_Y->GetRMS()-Zll0TrackMet_DA_Y->GetRMS()*Zll0TrackMet_DA_Y->GetRMS(),0.));
  printf(" biasCor0:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smearCor0: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);
  bias[0]  = Zll1Met_DA_X     ->GetMean()-Zll1Met_MC_X     ->GetMean();
  bias[1]  = Zll1Met_DA_Y     ->GetMean()-Zll1Met_MC_Y     ->GetMean();
  bias[2]  = Zll1TrackMet_DA_X->GetMean()-Zll1TrackMet_MC_X->GetMean();
  bias[3]  = Zll1TrackMet_DA_Y->GetMean()-Zll1TrackMet_MC_Y->GetMean();
  smear[0] = sqrt(TMath::Max(Zll1Met_DA_X     ->GetRMS()*Zll1Met_DA_X	  ->GetRMS()-Zll1Met_MC_X     ->GetRMS()*Zll1Met_MC_X	  ->GetRMS(),0.));
  smear[1] = sqrt(TMath::Max(Zll1Met_DA_Y     ->GetRMS()*Zll1Met_DA_Y	  ->GetRMS()-Zll1Met_MC_Y     ->GetRMS()*Zll1Met_MC_Y	  ->GetRMS(),0.));
  smear[2] = sqrt(TMath::Max(Zll1TrackMet_DA_X->GetRMS()*Zll1TrackMet_DA_X->GetRMS()-Zll1TrackMet_MC_X->GetRMS()*Zll1TrackMet_MC_X->GetRMS(),0.));
  smear[3] = sqrt(TMath::Max(Zll1TrackMet_DA_Y->GetRMS()*Zll1TrackMet_DA_Y->GetRMS()-Zll1TrackMet_MC_Y->GetRMS()*Zll1TrackMet_MC_Y->GetRMS(),0.));
  printf(" bias1:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smear1: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);
  bias[0]  = Zll1Met_DA_X     ->GetMean()-Zll1Met_CorMC_X     ->GetMean();
  bias[1]  = Zll1Met_DA_Y     ->GetMean()-Zll1Met_CorMC_Y     ->GetMean();
  bias[2]  = Zll1TrackMet_DA_X->GetMean()-Zll1TrackMet_CorMC_X->GetMean();
  bias[3]  = Zll1TrackMet_DA_Y->GetMean()-Zll1TrackMet_CorMC_Y->GetMean();
  smear[0] = sqrt(TMath::Max(Zll1Met_CorMC_X     ->GetRMS()*Zll1Met_CorMC_X     ->GetRMS()-Zll1Met_DA_X     ->GetRMS()*Zll1Met_DA_X     ->GetRMS(),0.));
  smear[1] = sqrt(TMath::Max(Zll1Met_CorMC_Y     ->GetRMS()*Zll1Met_CorMC_Y     ->GetRMS()-Zll1Met_DA_Y     ->GetRMS()*Zll1Met_DA_Y     ->GetRMS(),0.));
  smear[2] = sqrt(TMath::Max(Zll1TrackMet_CorMC_X->GetRMS()*Zll1TrackMet_CorMC_X->GetRMS()-Zll1TrackMet_DA_X->GetRMS()*Zll1TrackMet_DA_X->GetRMS(),0.));
  smear[3] = sqrt(TMath::Max(Zll1TrackMet_CorMC_Y->GetRMS()*Zll1TrackMet_CorMC_Y->GetRMS()-Zll1TrackMet_DA_Y->GetRMS()*Zll1TrackMet_DA_Y->GetRMS(),0.));
  printf(" biasCor1:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smearCor1: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);
  bias[0]  = Zll2Met_DA_X     ->GetMean()-Zll2Met_MC_X     ->GetMean();
  bias[1]  = Zll2Met_DA_Y     ->GetMean()-Zll2Met_MC_Y     ->GetMean();
  bias[2]  = Zll2TrackMet_DA_X->GetMean()-Zll2TrackMet_MC_X->GetMean();
  bias[3]  = Zll2TrackMet_DA_Y->GetMean()-Zll2TrackMet_MC_Y->GetMean();
  smear[0] = sqrt(TMath::Max(Zll2Met_DA_X     ->GetRMS()*Zll2Met_DA_X	  ->GetRMS()-Zll2Met_MC_X     ->GetRMS()*Zll2Met_MC_X	  ->GetRMS(),0.));
  smear[1] = sqrt(TMath::Max(Zll2Met_DA_Y     ->GetRMS()*Zll2Met_DA_Y	  ->GetRMS()-Zll2Met_MC_Y     ->GetRMS()*Zll2Met_MC_Y	  ->GetRMS(),0.));
  smear[2] = sqrt(TMath::Max(Zll2TrackMet_DA_X->GetRMS()*Zll2TrackMet_DA_X->GetRMS()-Zll2TrackMet_MC_X->GetRMS()*Zll2TrackMet_MC_X->GetRMS(),0.));
  smear[3] = sqrt(TMath::Max(Zll2TrackMet_DA_Y->GetRMS()*Zll2TrackMet_DA_Y->GetRMS()-Zll2TrackMet_MC_Y->GetRMS()*Zll2TrackMet_MC_Y->GetRMS(),0.));
  printf(" bias2:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smear2: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);
  bias[0]  = Zll2Met_DA_X     ->GetMean()-Zll2Met_CorMC_X     ->GetMean();
  bias[1]  = Zll2Met_DA_Y     ->GetMean()-Zll2Met_CorMC_Y     ->GetMean();
  bias[2]  = Zll2TrackMet_DA_X->GetMean()-Zll2TrackMet_CorMC_X->GetMean();
  bias[3]  = Zll2TrackMet_DA_Y->GetMean()-Zll2TrackMet_CorMC_Y->GetMean();
  smear[0] = sqrt(TMath::Max(Zll2Met_CorMC_X     ->GetRMS()*Zll2Met_CorMC_X     ->GetRMS()-Zll2Met_DA_X     ->GetRMS()*Zll2Met_DA_X     ->GetRMS(),0.));
  smear[1] = sqrt(TMath::Max(Zll2Met_CorMC_Y     ->GetRMS()*Zll2Met_CorMC_Y     ->GetRMS()-Zll2Met_DA_Y     ->GetRMS()*Zll2Met_DA_Y     ->GetRMS(),0.));
  smear[2] = sqrt(TMath::Max(Zll2TrackMet_CorMC_X->GetRMS()*Zll2TrackMet_CorMC_X->GetRMS()-Zll2TrackMet_DA_X->GetRMS()*Zll2TrackMet_DA_X->GetRMS(),0.));
  smear[3] = sqrt(TMath::Max(Zll2TrackMet_CorMC_Y->GetRMS()*Zll2TrackMet_CorMC_Y->GetRMS()-Zll2TrackMet_DA_Y->GetRMS()*Zll2TrackMet_DA_Y->GetRMS(),0.));
  printf(" biasCor2:  %8.5f  %8.5f  %8.5f  %8.5f\n",bias[0],bias[1],bias[2],bias[3]);
  printf("smearCor2: %8.5f  %8.5f  %8.5f  %8.5f\n",smear[0],smear[1],smear[2],smear[3]);

  TFile *f = new TFile("zyield.root", "UPDATE");
  f->cd();
  ZmmMass_MC   ->Write();
  ZmmMass_MC_BB->Write();
  ZmmMass_MC_EE->Write();
  ZmmMass_MC_EB->Write();
  ZeeMass_MC   ->Write();
  ZeeMass_MC_BB->Write();
  ZeeMass_MC_EE->Write();
  ZeeMass_MC_EB->Write();
  ZmmMass_DA   ->Write();
  ZmmMass_DA_BB->Write();
  ZmmMass_DA_EE->Write();
  ZmmMass_DA_EB->Write();
  ZeeMass_DA   ->Write();
  ZeeMass_DA_BB->Write();
  ZeeMass_DA_EE->Write();
  ZeeMass_DA_EB->Write();
  ZmmMass_CorMC   ->Write();
  ZmmMass_CorMC_BB->Write();
  ZmmMass_CorMC_EE->Write();
  ZmmMass_CorMC_EB->Write();
  ZeeMass_CorMC   ->Write();
  ZeeMass_CorMC_BB->Write();
  ZeeMass_CorMC_EE->Write();
  ZeeMass_CorMC_EB->Write();
  Zll0Met_MC_X         ->Write();
  Zll0Met_MC_Y         ->Write();
  Zll0TrackMet_MC_X    ->Write();
  Zll0TrackMet_MC_Y    ->Write();
  Zll0Met_DA_X         ->Write();
  Zll0Met_DA_Y         ->Write();
  Zll0TrackMet_DA_X    ->Write();
  Zll0TrackMet_DA_Y    ->Write();
  Zll0Met_CorMC_X      ->Write();
  Zll0Met_CorMC_Y      ->Write();
  Zll0TrackMet_CorMC_X ->Write();
  Zll0TrackMet_CorMC_Y ->Write();
  Zll1Met_MC_X         ->Write();
  Zll1Met_MC_Y         ->Write();
  Zll1TrackMet_MC_X    ->Write();
  Zll1TrackMet_MC_Y    ->Write();
  Zll1Met_DA_X         ->Write();
  Zll1Met_DA_Y         ->Write();
  Zll1TrackMet_DA_X    ->Write();
  Zll1TrackMet_DA_Y    ->Write();
  Zll1Met_CorMC_X      ->Write();
  Zll1Met_CorMC_Y      ->Write();
  Zll1TrackMet_CorMC_X ->Write();
  Zll1TrackMet_CorMC_Y ->Write();
  Zll2Met_MC_X         ->Write();
  Zll2Met_MC_Y         ->Write();
  Zll2TrackMet_MC_X    ->Write();
  Zll2TrackMet_MC_Y    ->Write();
  Zll2Met_DA_X         ->Write();
  Zll2Met_DA_Y         ->Write();
  Zll2TrackMet_DA_X    ->Write();
  Zll2TrackMet_DA_Y    ->Write();
  Zll2Met_CorMC_X      ->Write();
  Zll2Met_CorMC_Y      ->Write();
  Zll2TrackMet_CorMC_X ->Write();
  Zll2TrackMet_CorMC_Y ->Write();
  f->Close();
}
