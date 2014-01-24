// @(#)root/tmva $Id: trainMVA_smurf_MultiClass.C,v 1.5 2012/11/26 09:48:34 ceballos Exp $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include "TRandom.h"

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChainElement.h"

#include "TMVAGui.C"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "../Core/SmurfTree.h"
#include "../Analysis/HWWlvlv/factors.h"
#include "../Analysis/HWWlvlv/HWWCuts.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

/// first is leading lepton
enum Type {
  mm, 
  me, 
  em, 
  ee
};

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void trainMVA_smurf_MultiClass(
 UInt_t  nJetsType       = 0,			   
 TString sigInputFile    = "data/histo_hww130_std_pu11_randomized.training.root",
 TString bgdInputFile    = "data/histo_w10-ggww-z2-v8-pu11_all_noskim_normalized.root",
 TString outTag          = "default",
 TString myMethodList    = "Fisher,BDT",
 UInt_t  mH              = 130
 )
{
  // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
  // if you use your private .rootrc, or run from a different directory, please copy the
  // corresponding lines from .rootrc

  // methods to be processed can be given as an argument; use format:
  //
  // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
  //
  // if you like to use a method via the plugin mechanism, we recommend using
  //
  // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
  // (an example is given for using the BDT as plugin (see below),
  // but of course the real application is when you write your own
  // method based)
 
  //-----------------------------------------------------
  // choose higgs mass
  //-----------------------------------------------------

  cout << "Training MVA Higgs mass " << mH << endl;
   
  //-------------------------------------------------------------------------------------
  // define event selection
  //
  // current selection is: WW selection, except anti b-tagging and soft muon veto
  //
  // The event selection is applied below, see: 
  // SIGNAL EVENT SELECTION and BACKGROUND EVENT SELECTION
  // if you want to apply any additional selection this needs to be implemented below
  //-------------------------------------------------------------------------------------

  //-----------------------------------------------------
  // choose which variables to include in MVA training
  //-----------------------------------------------------
  std::map<std::string,int> mvaVar;
  if(nJetsType == 0 || nJetsType == 1){
    mvaVar[ "lep1pt" ]            = 1;  //pt of leading lepton
    mvaVar[ "lep2pt" ]            = 1;  //pt of sub-leading lepton
    mvaVar[ "dPhi" ]              = 1;  //delta phi btw leptons
    mvaVar[ "dR" ]                = 1;  //delta R btw leptons
    mvaVar[ "dilmass" ]           = 1;  //dilepton mass
    mvaVar[ "pmet" ]              = 0;  //projected met
    mvaVar[ "met" ]               = 0;  //met
    mvaVar[ "mt" ]                = 1;  //transverse higgs mass
    mvaVar[ "dPhiLep1MET" ]       = 0;  //delta phi btw leading lepton and met
    mvaVar[ "dPhiLep2MET" ]       = 0;  //delta phi btw leading sub-lepton and met
    mvaVar[ "razor" ]             = 0;  //razor
    mvaVar[ "mt1" ]		  = 0;  //transverse mass of leading lepton and met
    mvaVar[ "mt2" ]		  = 0;  //transverse mass of sub-leading lepton and met
    mvaVar[ "dilpt" ]		  = 1;  //dilepton pt
    mvaVar[ "type" ]		  = 0;  //dilepton flavor type
  }
  if(nJetsType == 1){
    mvaVar[ "dPhiDiLepMET" ]   = 0; //delta phi btw dilepton and met
    mvaVar[ "dPhiDiLepJet1" ]  = 0; //delta phi btw dilepton and jet1
  }
  if (nJetsType == 2) {
    mvaVar[ "mjj" ]             = 1;  //invariant mass of the dijet
    mvaVar[ "detajj" ]          = 1;  //deta of the dijet
    mvaVar[ "dphijj" ]          = 1;  //dphi of the dijet
    mvaVar[ "ptjj" ]            = 1;  //pt of the dijet
    mvaVar[ "higgspt" ]         = 1;  //higgs pt
    mvaVar[ "dphihjj" ]         = 1;  //dphi between the Higgs and the jj system
  }
  TCut sel = "";
  
  float dilmass_cut = DileptonMassPreselectionCut(mH);
  if     (mH <= 160) dilmass_cut = 100;
  else if(mH <= 250) dilmass_cut = 150;
  else               dilmass_cut = 600;
  float mt_cut = mH;
  if(mH <= 250) mt_cut = 280;
  else          mt_cut = 600;

  cout << "Using dilepton mass < " << dilmass_cut << endl;
  cout << "Using mt mass       < " << mt_cut      << endl;

  double rndLim[3] = {1.0,1.0,1.0};
  if     (mH==125 && nJetsType == 0) {rndLim[0] = 0.35;rndLim[1] = 1.00;rndLim[2] = 0.07;}
  else if(mH==125 && nJetsType == 1) {rndLim[0] = 0.35;rndLim[1] = 1.00;rndLim[2] = 0.10;}
  else if(mH==200 && nJetsType == 0) {rndLim[0] = 0.25;rndLim[1] = 1.00;rndLim[2] = 0.11;}
  else if(mH==200 && nJetsType == 1) {rndLim[0] = 0.25;rndLim[1] = 1.00;rndLim[2] = 0.17;}
  else if(mH==500 && nJetsType == 0) {rndLim[0] = 0.50;rndLim[1] = 1.00;rndLim[2] = 0.15;}
  else if(mH==500 && nJetsType == 1) {rndLim[0] = 0.50;rndLim[1] = 1.00;rndLim[2] = 0.20;}

  //---------------------------------
  //choose bkg samples to include
  //---------------------------------
  
  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdInputFile);

  //---------------------------------
  //choose signal sample to include
  //---------------------------------

  TChain *chsignal = new TChain("tree");
  chsignal->Add(sigInputFile);

  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // --- Cut optimisation
  Use["Cuts"]            = 1;
  Use["CutsD"]           = 1;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 1;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 1;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 1;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 1; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 1; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 1;
  // 
  // --- Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 1;
  //
  // --- multi-output MVA's
  Use["multi_BDTG"]      = 1;
  Use["multi_MLP"]       = 1;
  Use["multi_FDA_GA"]    = 0;
  //
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
        std::cout << std::endl;
        return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  // --- Here the preparation phase begins

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName = outTag + ".root";
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is 
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TMVA::Factory *factory = new TMVA::Factory( outTag.Data(), outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  //factory->AddVariable( "myvar1 := var1+var2", 'F' );
  //factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
  //factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
  //factory->AddVariable( "var4",                "Variable 4", "units", 'F' );

  if (mvaVar["lep1pt"])       factory->AddVariable( "lep1pt",        "1st lepton pt",       "GeV", 'F' );
  if (mvaVar["lep2pt"])       factory->AddVariable( "lep2pt",        "2nd lepton pt",       "GeV", 'F' );
  if (mvaVar["dPhi"])         factory->AddVariable( "dPhi",          "dphi(ll)",            "",    'F' );
  if (mvaVar["dR"])           factory->AddVariable( "dR",            "dr(ll)",              "",    'F' );
  if (mvaVar["dilmass"])      factory->AddVariable( "dilmass",       "M(ll)",               "GeV", 'F' );
  if (mvaVar["type"])         factory->AddVariable( "type",          "Dil Flavor Type",     "",    'F' );
  if (mvaVar["pmet"])         factory->AddVariable( "pmet",          "Proj. MET",           "GeV", 'F' );
  if (mvaVar["met"])          factory->AddVariable( "met",           "MET",                 "GeV", 'F' );
  if (mvaVar["mt"])           factory->AddVariable( "mt",            "MT(Higgs)",           "GeV", 'F' );
  if (mvaVar["mt1"])          factory->AddVariable( "mt1",           "MT(lep1,MET)",        "GeV", 'F' );
  if (mvaVar["mt2"])          factory->AddVariable( "mt2",           "MT(lep2,MET)",        "GeV", 'F' );
  if (mvaVar["dPhiLep1MET"])  factory->AddVariable( "dPhiLep1MET",   "dphi(lep1,MET)",      "GeV", 'F' );
  if (mvaVar["dPhiLep2MET"])  factory->AddVariable( "dPhiLep2MET",   "dphi(lep2,MET)",      "GeV", 'F' );
  if (mvaVar["dilpt"])        factory->AddVariable( "dilpt",         "dilepton pt",         "GeV", 'F' );
  if (mvaVar["razor"])        factory->AddVariable( "razor",         "razor",               "",    'F' );
  if (mvaVar["dPhiDiLepMET"]) factory->AddVariable( "dPhiDiLepMET",  "dphi(dilep,MET)",     "",    'F' );
  if (mvaVar["dPhiDiLepJet1"])factory->AddVariable( "dPhiDiLepJet1", "dphi(dilep,jet)",     "",    'F' );
  if (mvaVar["mjj"])	      factory->AddVariable( "mjj",	     "mjj",		    "",	   'F' );
  if (mvaVar["detajj"])       factory->AddVariable( "detajj",	     "detajj",		    "",	   'F' );
  if (mvaVar["dphijj"])       factory->AddVariable( "dphijj",	     "dphijj",		    "",	   'F' );
  if (mvaVar["ptjj"])	      factory->AddVariable( "ptjj",	     "ptjj",		    "",	   'F' );
  if (mvaVar["higgspt"])      factory->AddVariable( "higgspt",	     "higgspt",		    "",	   'F' );
  if (mvaVar["dphihjj"])      factory->AddVariable( "dphihjj",	     "dphihjj",		    "",	   'F' );

  int nVariablesTemp = 0;

  if (mvaVar["lep1pt"])       { cout << "Adding variable to MVA training: lep1.pt()"      << endl; nVariablesTemp++; }
  if (mvaVar["lep2pt"])       { cout << "Adding variable to MVA training: lep2.pt()"      << endl; nVariablesTemp++; }
  if (mvaVar["dPhi"])         { cout << "Adding variable to MVA training: dPhi"           << endl; nVariablesTemp++; }
  if (mvaVar["dR"])           { cout << "Adding variable to MVA training: dR"             << endl; nVariablesTemp++; }
  if (mvaVar["dilmass"])      { cout << "Adding variable to MVA training: dilep.mass()"   << endl; nVariablesTemp++; }
  if (mvaVar["type"])         { cout << "Adding variable to MVA training: type"           << endl; nVariablesTemp++; }
  if (mvaVar["pmet"])         { cout << "Adding variable to MVA training: pmet"           << endl; nVariablesTemp++; }
  if (mvaVar["met"])          { cout << "Adding variable to MVA training: met"            << endl; nVariablesTemp++; }
  if (mvaVar["mt"])           { cout << "Adding variable to MVA training: mt"             << endl; nVariablesTemp++; }
  if (mvaVar["mt1"])          { cout << "Adding variable to MVA training: mt1"            << endl; nVariablesTemp++; }
  if (mvaVar["mt2"])          { cout << "Adding variable to MVA training: mt2"            << endl; nVariablesTemp++; }
  if (mvaVar["dPhiLep1MET"])  { cout << "Adding variable to MVA training: dPhiLep1MET"    << endl; nVariablesTemp++; }
  if (mvaVar["dPhiLep2MET"])  { cout << "Adding variable to MVA training: dPhiLep2MET"    << endl; nVariablesTemp++; }
  if (mvaVar["dilpt"])        { cout << "Adding variable to MVA training: dil.pt()"       << endl; nVariablesTemp++; }
  if (mvaVar["razor"])        { cout << "Adding variable to MVA training: razor"          << endl; nVariablesTemp++; }
  if (mvaVar["dPhiDiLepMET"]) { cout << "Adding variable to MVA training: dPhiDiLepMET"   << endl; nVariablesTemp++; }
  if (mvaVar["dPhiDiLepJet1"]){ cout << "Adding variable to MVA training: dPhiDiLepJet1"  << endl; nVariablesTemp++; }
  if (mvaVar["mjj"])	      { cout << "Adding variable to MVA training: mjj"  	  << endl; nVariablesTemp++; }
  if (mvaVar["detajj"])       { cout << "Adding variable to MVA training: detajj"	  << endl; nVariablesTemp++; }
  if (mvaVar["dphijj"])       { cout << "Adding variable to MVA training: dphijj"	  << endl; nVariablesTemp++; }
  if (mvaVar["ptjj"])	      { cout << "Adding variable to MVA training: ptjj"  	  << endl; nVariablesTemp++; }
  if (mvaVar["higgspt"])      { cout << "Adding variable to MVA training: higgspt"	  << endl; nVariablesTemp++; }
  if (mvaVar["dphihjj"])      { cout << "Adding variable to MVA training: dphihjj"	  << endl; nVariablesTemp++; }

  const unsigned int nVariables = nVariablesTemp;
  cout << "Using " << nVariables << " variables for MVA training" << endl;

  // You can add so-called "Spectator variables", which are not used in the MVA training,
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
  // input variables, the response values of all trained MVAs, and the spectator variables
  //factory->AddSpectator( "njets",  "Event Type", "units", 'F' );
  //factory->AddSpectator( "mydilmass := ",  "Spectator 2", "units", 'F' );

  TTree *signal     = (TTree*) chsignal;
  TTree *background = (TTree*) chbackground;
   
  std::cout << "--- TMVAClassification       : Using bkg input files: -------------------" <<  std::endl;

  TObjArray *listOfBkgFiles = chbackground->GetListOfFiles();
  TIter bkgFileIter(listOfBkgFiles);
  TChainElement* currentBkgFile = 0;

  while((currentBkgFile = (TChainElement*)bkgFileIter.Next())) {
    std::cout << currentBkgFile->GetTitle() << std::endl;
  }

  std::cout << "--- TMVAClassification       : Using sig input files: -------------------" <<  std::endl;
   
  TObjArray *listOfSigFiles = chsignal->GetListOfFiles();
  TIter sigFileIter(listOfSigFiles);
  TChainElement* currentSigFile = 0;

  while((currentSigFile = (TChainElement*)sigFileIter.Next())) {
    std::cout << currentSigFile->GetTitle() << std::endl;
  }

  // global event weights per tree (see below for setting event-wise weights)
  //Double_t signalWeight     = 1.0;
  //Double_t backgroundWeight = 1.0;
   
  // You can add an arbitrary number of signal or background trees
  //factory->AddSignalTree    ( signal,     signalWeight     );
  //factory->AddBackgroundTree( background, backgroundWeight );
      
  // To give different trees for training and testing, do as follows:
  //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
  //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

  // Use the following code instead of the above two or four lines to add signal and background
  // training and test events "by hand"
  // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
  //      variable definition, but simply compute the expression before adding the event
  //
  
  std::vector<Double_t> vars( nVariables );

  UInt_t          processId;
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          njets;
  UInt_t          event;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  jet1  = 0;
  LorentzVector*  jet2  = 0;
  LorentzVector*  jet3  = 0;
  Float_t         dPhi;
  Float_t         dR;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         met;
  Float_t         metPhi;
  Float_t         mt;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         dPhiLep1MET;
  Float_t         dPhiLep2MET;
  Float_t         dPhiDiLepMET;
  Float_t         dPhiDiLepJet1;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid3;
  Float_t         jetLowBtag;
  UInt_t          nSoftMuons;
  Float_t         jet1Btag;

  unsigned int patternTopTag = SmurfTree::TopTag;

  signal->SetBranchAddress( "processId"    , &processId    );
  signal->SetBranchAddress( "cuts"         , &cuts         );
  signal->SetBranchAddress( "dstype"       , &dstype       );
  signal->SetBranchAddress( "njets"        , &njets        );
  signal->SetBranchAddress( "event"        , &event        );
  signal->SetBranchAddress( "scale1fb"     , &scale1fb     );
  signal->SetBranchAddress( "lep1"         , &lep1         );
  signal->SetBranchAddress( "lep2"         , &lep2         );
  signal->SetBranchAddress( "jet1"         , &jet1         );
  signal->SetBranchAddress( "jet2"         , &jet2         );
  signal->SetBranchAddress( "jet3"         , &jet3         );
  signal->SetBranchAddress( "dPhi"         , &dPhi         );
  signal->SetBranchAddress( "dR"           , &dR           );
  signal->SetBranchAddress( "dilep"        , &dilep        );
  signal->SetBranchAddress( "type"         , &type         );
  signal->SetBranchAddress( "pmet"         , &pmet         );
  signal->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  signal->SetBranchAddress( "met"          , &met          );
  signal->SetBranchAddress( "metPhi"       , &metPhi       );
  signal->SetBranchAddress( "mt"           , &mt           );
  signal->SetBranchAddress( "mt1"          , &mt1          );
  signal->SetBranchAddress( "mt2"          , &mt2          );
  signal->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  signal->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  signal->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  signal->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  signal->SetBranchAddress( "lq1"          , &lq1          );
  signal->SetBranchAddress( "lq2"          , &lq2          );
  signal->SetBranchAddress( "lid3"         , &lid3         );
  signal->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  signal->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  signal->SetBranchAddress( "jet1Btag"     , &jet1Btag     );

  int nsigtrain = 0;
  int nsigtest  = 0;
  int nbkgtrain[3] = {0,0,0};
  int nbkgtest[3]  = {0,0,0};

  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);

    //--------------------------------------------------
    // SIGNAL EVENT SELECTION
    //--------------------------------------------------

    int classType = 0;
    unsigned int Njet3 = njets;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(jet3->pt() <= 30)					           Njet3 = 2;
      else if(jet3->pt() > 30 && (
    	(jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    	(jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 3;
      else							           Njet3 = 2;
      if(njets < 2 || njets > 3)                                           Njet3 = 3;
      if(processId != 10001) classType = 1;
    }

    if(!((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
         (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)) continue; // two good leptons

    if( Njet3 != nJetsType                      ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut             ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                             ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0                   ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	                ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	                ) continue; // cut on trailing lepton pt
    if( TMath::Min(pmet,pTrackMet) <= 20        ) continue; // cut on pmet for all lepton-pair flavors
    if( lid3 != 0	                        ) continue; // cut on dileptons
    if( dilep->pt() <= 45                       ) continue; // cut on dilepton pt
    if( (cuts & patternTopTag) == patternTopTag ) continue; // cut on btagging
    if( njets >= 1 && dPhiDiLepJet1 > 100       ) continue; // cut on impossible dPhiDiLepJet1
    if( mt <= 80                                ) continue;
    if( mt >= mt_cut                            ) continue;

    // VBF Cuts
    if ( nJetsType == 2 ) {
      // centrality cut
      float largestEta = jet1->Eta();
      float smallestEta = jet2->Eta();
      if (jet2->Eta() > jet1->Eta()) {
	largestEta = jet2->Eta();
	smallestEta = jet1->Eta();
      }
      if( !(lep1->Eta() > smallestEta && lep1->Eta() < largestEta)) continue;
      if(! (lep2->Eta() > smallestEta && lep2->Eta() < largestEta)) continue;
      if( TMath::Abs(jet1->Eta() - jet2->Eta()) <= 2.0            ) continue;
      if( (*jet1+*jet2).M() <= 200                                ) continue;
    }

    int varCounter = 0;
    
    double higgsV[2] = {lep1->px()+lep2->px()+met*cos(metPhi),
                        lep1->py()+lep2->py()+met*sin(metPhi)};
    if (mvaVar["lep1pt"])        vars[varCounter++] = lep1->pt();
    if (mvaVar["lep2pt"])        vars[varCounter++] = lep2->pt();    
    if (mvaVar["dPhi"])          vars[varCounter++] = dPhi;
    if (mvaVar["dR"])            vars[varCounter++] = dR;
    if (mvaVar["dilmass"])       vars[varCounter++] = dilep->mass();
    if (mvaVar["type"])          vars[varCounter++] = type;
    if (mvaVar["pmet"])          vars[varCounter++] = pmet;
    if (mvaVar["met"])           vars[varCounter++] = met;
    if (mvaVar["mt"])            vars[varCounter++] = mt;
    if (mvaVar["mt1"])           vars[varCounter++] = mt1;
    if (mvaVar["mt2"])           vars[varCounter++] = mt2;
    if (mvaVar["dPhiLep1MET"])   vars[varCounter++] = dPhiLep1MET;
    if (mvaVar["dPhiLep2MET"])   vars[varCounter++] = dPhiLep2MET;
    if (mvaVar["dilpt"])         vars[varCounter++] = dilep->pt();
    if (mvaVar["razor"])         vars[varCounter++] = CalcGammaMRstar(*lep1,*lep2);
    if (mvaVar["dPhiDiLepMET"])  vars[varCounter++] = dPhiDiLepMET;
    if (mvaVar["dPhiDiLepJet1"]) vars[varCounter++] = dPhiDiLepJet1;
    if (mvaVar["mjj"])    	 vars[varCounter++] = (*jet1+*jet2).M();
    if (mvaVar["detajj"]) 	 vars[varCounter++] = TMath::Abs(jet1->Eta() - jet2->Eta());
    if (mvaVar["dphijj"])  	 vars[varCounter++] = DeltaPhi(jet1->Phi(),jet2->Phi());
    if (mvaVar["ptjj"])   	 vars[varCounter++] = (*jet1+*jet2).Pt();
    if (mvaVar["higgspt"])	 vars[varCounter++] = sqrt(higgsV[0]*higgsV[0]+higgsV[1]*higgsV[1]);
    if (mvaVar["dphihjj"]) 	 vars[varCounter++] = DeltaPhi((*jet1+*jet2).Phi(),TMath::ATan2(higgsV[1],higgsV[0]));

    if ( gRandom->Uniform(0,1) < 0.5 ){
      if     (classType == 0 && gRandom->Uniform(0,1) < rndLim[0]) {factory->AddTestEvent( "sig" , vars, scale1fb ); nsigtrain++;}
      else if(classType == 1 && gRandom->Uniform(0,1) < rndLim[1]) {factory->AddTestEvent( "bkg0", vars, scale1fb ); assert(0);}
      else if(classType == 2 && gRandom->Uniform(0,1) < rndLim[2]) {factory->AddTestEvent( "bkg1", vars, scale1fb ); assert(0);}
    }
    else{
      if     (classType == 0 && gRandom->Uniform(0,1) < rndLim[0]) {factory->AddTrainingEvent( "sig" , vars, scale1fb ); nsigtest++;}
      else if(classType == 1 && gRandom->Uniform(0,1) < rndLim[1]) {factory->AddTrainingEvent( "bkg0", vars, scale1fb ); assert(0);}
      else if(classType == 2 && gRandom->Uniform(0,1) < rndLim[2]) {factory->AddTrainingEvent( "bkg1", vars, scale1fb ); assert(0);}
    }
  }

  background->SetBranchAddress( "cuts"         , &cuts         );
  background->SetBranchAddress( "dstype"       , &dstype       );
  background->SetBranchAddress( "njets"        , &njets        );  
  background->SetBranchAddress( "event"        , &event        );
  background->SetBranchAddress( "scale1fb"     , &scale1fb     );
  background->SetBranchAddress( "lep1"         , &lep1         );
  background->SetBranchAddress( "lep2"         , &lep2         );
  background->SetBranchAddress( "jet1"         , &jet1         );
  background->SetBranchAddress( "jet2"         , &jet2         );
  background->SetBranchAddress( "jet3"         , &jet3         );
  background->SetBranchAddress( "dPhi"         , &dPhi         );
  background->SetBranchAddress( "dR"           , &dR           );
  background->SetBranchAddress( "dilep"        , &dilep        );
  background->SetBranchAddress( "type"         , &type         );
  background->SetBranchAddress( "pmet"         , &pmet         );
  background->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  background->SetBranchAddress( "met"          , &met          );
  background->SetBranchAddress( "metPhi"       , &metPhi       );
  background->SetBranchAddress( "mt"           , &mt           );
  background->SetBranchAddress( "mt1"          , &mt1          );
  background->SetBranchAddress( "mt2"          , &mt2          );
  background->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  background->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  background->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  background->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  background->SetBranchAddress( "lq1"	       , &lq1	       );
  background->SetBranchAddress( "lq2"	       , &lq2	       );
  background->SetBranchAddress( "lid3"         , &lid3         );
  background->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  background->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  background->SetBranchAddress( "jet1Btag"     , &jet1Btag     );

  cout << "Add signal events" << endl;
  cout << "Added " << nsigtrain << " training events" << endl;
  cout << "Added " << nsigtest  << " test events" << endl;
  
  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);

    //--------------------------------------------------
    // BACKGROUND EVENT SELECTION
    //--------------------------------------------------

    unsigned int Njet3 = njets;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(jet3->pt() <= 30)					           Njet3 = 2;
      else if(jet3->pt() > 30 && (
    	(jet1->eta()-jet3->eta() > 0 && jet2->eta()-jet3->eta() < 0) ||
    	(jet2->eta()-jet3->eta() > 0 && jet1->eta()-jet3->eta() < 0)))     Njet3 = 3;
      else							           Njet3 = 2;
      if(njets < 2 || njets > 3)                                           Njet3 = 3;
    }

    bool lId = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
                (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
		dstype == SmurfTree::wgamma;
    if(!lId) continue; // two good leptons || wgamma

    if( Njet3 != nJetsType                      ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut             ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                             ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0                   ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	                ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	                ) continue; // cut on trailing lepton pt
    if( TMath::Min(pmet,pTrackMet) <= 20        ) continue; // cut on pmet for all lepton-pair flavors
    if( lid3 != 0	                        ) continue; // cut on dileptons
    if( dilep->pt() <= 45                       ) continue; // cut on dilepton pt
    if( (cuts & patternTopTag) == patternTopTag ) continue; // cut on btagging
    if( njets >= 1 && dPhiDiLepJet1 > 100       ) continue; // cut on impossible dPhiDiLepJet1
    if( mt <= 80                                ) continue;
    if( mt >= mt_cut                            ) continue;

    // VBF Cuts
    int classType = 1;
    if ( nJetsType == 2 ) {
      // centrality cut
      float largestEta = jet1->Eta();
      float smallestEta = jet2->Eta();
      if (jet2->Eta() > jet1->Eta()) {
	largestEta = jet2->Eta();
	smallestEta = jet1->Eta();
      }
      if( !(lep1->Eta() > smallestEta && lep1->Eta() < largestEta)) continue;
      if(! (lep2->Eta() > smallestEta && lep2->Eta() < largestEta)) continue;
      if( TMath::Abs(jet1->Eta() - jet2->Eta()) <= 2.0            ) continue;
      if( (*jet1+*jet2).M() <= 200                                ) continue;
      classType = 2;
    }
    
    if     (dstype == SmurfTree::wgamma) classType = 1;
    else if(dstype == SmurfTree::ttbar)  classType = 2;
    else if(dstype == SmurfTree::tw)	 classType = 2;
    else if(dstype == SmurfTree::qqww)   classType = 2;
    else if(dstype == SmurfTree::ggww)   classType = 2;
    else				 classType = 3;
    
    if(classType == 3) continue;

    int varCounter = 0;
    
    double higgsV[2] = {lep1->px()+lep2->px()+met*cos(metPhi),
                        lep1->py()+lep2->py()+met*sin(metPhi)};
    if (mvaVar["lep1pt"])        vars[varCounter++] = lep1->pt();
    if (mvaVar["lep2pt"])        vars[varCounter++] = lep2->pt();    
    if (mvaVar["dPhi"])          vars[varCounter++] = dPhi;
    if (mvaVar["dR"])            vars[varCounter++] = dR;
    if (mvaVar["dilmass"])       vars[varCounter++] = dilep->mass();
    if (mvaVar["type"])          vars[varCounter++] = type;
    if (mvaVar["pmet"])          vars[varCounter++] = pmet;
    if (mvaVar["met"])           vars[varCounter++] = met;
    if (mvaVar["mt"])            vars[varCounter++] = mt;
    if (mvaVar["mt1"])           vars[varCounter++] = mt1;
    if (mvaVar["mt2"])           vars[varCounter++] = mt2;
    if (mvaVar["dPhiLep1MET"])   vars[varCounter++] = dPhiLep1MET;
    if (mvaVar["dPhiLep2MET"])   vars[varCounter++] = dPhiLep2MET;
    if (mvaVar["dilpt"])         vars[varCounter++] = dilep->pt();
    if (mvaVar["razor"])         vars[varCounter++] = CalcGammaMRstar(*lep1,*lep2);
    if (mvaVar["dPhiDiLepMET"])  vars[varCounter++] = dPhiDiLepMET;
    if (mvaVar["dPhiDiLepJet1"]) vars[varCounter++] = dPhiDiLepJet1;
    if (mvaVar["mjj"])    	 vars[varCounter++] = (*jet1+*jet2).M();
    if (mvaVar["detajj"]) 	 vars[varCounter++] = TMath::Abs(jet1->Eta() - jet2->Eta());
    if (mvaVar["dphijj"])  	 vars[varCounter++] = DeltaPhi(jet1->Phi(),jet2->Phi());
    if (mvaVar["ptjj"])   	 vars[varCounter++] = (*jet1+*jet2).Pt();
    if (mvaVar["higgspt"])	 vars[varCounter++] = sqrt(higgsV[0]*higgsV[0]+higgsV[1]*higgsV[1]);
    if (mvaVar["dphihjj"]) 	 vars[varCounter++] = DeltaPhi((*jet1+*jet2).Phi(),TMath::ATan2(higgsV[1],higgsV[0]));

    if ( gRandom->Uniform(0,1) < 0.5 ){
      if     (classType == 0 && gRandom->Uniform(0,1) < rndLim[0]) {factory->AddTestEvent( "sig" , vars, scale1fb );nbkgtrain[0]++;}
      else if(classType == 1 && gRandom->Uniform(0,1) < rndLim[1]) {factory->AddTestEvent( "bkg0", vars, scale1fb );nbkgtrain[1]++;}
      else if(classType == 2 && gRandom->Uniform(0,1) < rndLim[2]) {factory->AddTestEvent( "bkg1", vars, scale1fb );nbkgtrain[2]++;}
      
    }
    else{
      if     (classType == 0 && gRandom->Uniform(0,1) < rndLim[0]) {factory->AddTrainingEvent( "sig" , vars, scale1fb );nbkgtest[0]++;}
      else if(classType == 1 && gRandom->Uniform(0,1) < rndLim[1]) {factory->AddTrainingEvent( "bkg0", vars, scale1fb );nbkgtest[1]++;}
      else if(classType == 2 && gRandom->Uniform(0,1) < rndLim[2]) {factory->AddTrainingEvent( "bkg1", vars, scale1fb );nbkgtest[2]++;}
    }
  }
  
  cout << "Done adding background" << endl;
  cout << "Added0 " << nbkgtrain[0] << " training events" << endl;
  cout << "Added1 " << nbkgtrain[1] << " training events" << endl;
  cout << "Added2 " << nbkgtrain[2] << " training events" << endl;
  cout << "Added0 " << nbkgtest[0]  << " test events" << endl;
  cout << "Added1 " << nbkgtest[1]  << " test events" << endl;
  cout << "Added2 " << nbkgtest[2]  << " test events" << endl;
  

  // --- end ------------------------------------------------------------
  //
  // --- end of tree registration 
   
  // Set individual event weights (the variables must exist in the original TTree)
  //factory->SetSignalWeightExpression    ("scale1fb");
  //factory->SetBackgroundWeightExpression("scale1fb");
  cout << "Done setting weights" << endl;

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = sel; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = sel; // for example: TCut mycutb = "abs(var1)<0.5";

  // Tell the factory how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used 
  // for training, and the other half for testing:
  //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  // To also specify the number of testing events, use:
  //    factory->PrepareTrainingAndTestTree( mycut,
  //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   
  //Use random splitting
  factory->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );

  //Use alternate splitting 
  //(this is preferable since its easier to track which events were used for training, but the job crashes! need to fix this...)
  //factory->PrepareTrainingAndTestTree( mycuts, mycutb,
  //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Alternate:NormMode=NumEvents:!V" );

  // ---- Book MVA methods
  //
  // Please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

  if (Use["CutsD"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

  if (Use["CutsPCA"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

  if (Use["CutsGA"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                         "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

  if (Use["CutsSA"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                         "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
    factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                         "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=30:NSmoothBkg[0]=30:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

  // Decorrelated likelihood
  if (Use["LikelihoodD"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                         "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=40:VarTransform=Decorrelate" );

  // PCA-transformed likelihood
  if (Use["LikelihoodPCA"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                         "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

  // Use a kernel density estimator to approximate the PDFs
  if (Use["LikelihoodKDE"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                         "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

  // Use a variable-dependent mix of splines and kernel density estimator
  if (Use["LikelihoodMIX"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                         "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

  // Test the multi-dimensional probability density estimator
  // here are the options strings for the MinMax and RMS methods, respectively:
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
  if (Use["PDERS"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                         "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

  if (Use["PDERSD"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                         "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

  if (Use["PDERSPCA"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                         "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

  // Multi-dimensional likelihood estimator using self-adapting phase-space binning
  if (Use["PDEFoam"])
    factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                         "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

  if (Use["PDEFoamBoost"])
    factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                         "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"])
    factory->BookMethod( TMVA::Types::kKNN, "KNN",
                         "H:nkNN=31:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

  // H-Matrix (chi2-squared) method
  if (Use["HMatrix"])
    factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

  // Linear discriminant (same as Fisher discriminant)
  if (Use["LD"])
    factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=40:NsmoothMVAPdf=10" );

  // Fisher with Gauss-transformed input variables
  if (Use["FisherG"])
    factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

  // Composite classifier: ensemble (tree) of boosted Fisher classifiers
  if (Use["BoostedFisher"])
    factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                         "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" );

  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

  if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

  if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  if (Use["FDA_MT"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

  if (Use["FDA_GAMT"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

  if (Use["FDA_MCMT"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  if (Use["MLPBFGS"])
    factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

  if (Use["MLPBNN"])
    factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+3:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

  // CF(Clermont-Ferrand)ANN
  if (Use["CFMlpANN"])
    factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

  // Tmlp(Root)ANN
  if (Use["TMlpANN"])
    factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

  // Support Vector Machine
  if (Use["SVM"])
    factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                         "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=1000:NNodesMax=5" );

  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod(TMVA::Types::kBDT,"BDT",
			 "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=2000:NNodesMax=5:VarTransform=Decorrelate");

  if (Use["BDTB"]) // Bagging
    factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                         "!H:!V:NTrees=500:BoostType=Bagging:SeparationType=GiniIndex:nCuts=2000:PruneMethod=NoPruning" );

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod(TMVA::Types::kBDT,"BDTD",
                    "!H:!V:NTrees=500:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=2000:PruneMethod=CostComplexity:PruneStrength=25.0:VarTransform=Decorrelate");

  // RuleFit -- TMVA implementation of Friedman's method
  if (Use["RuleFit"])
    factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                         "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   
  // For an example of the category classifier usage, see: TMVAClassificationCategory

  // --------------------------------------------------------------------------------------------------

  // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

  // factory->OptimizeAllMethods("SigEffAt001","Scan");
  // factory->OptimizeAllMethods("ROCIntegral","GA");

  // --------------------------------------------------------------------------------------------------

  // ---- Now you can tell the factory to train, test, and evaluate the MVAs
  
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  
  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
