// @(#)root/tmva $Id: trainMVA_smurf_wh3l.C,v 1.1 2012/03/26 09:02:40 ceballos Exp $
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
#include "../../Ana/nt_scripts/trilepton.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void trainMVA_smurf_wh3l(
 UInt_t  njetsType       = 999,			   
 TString sigInputFile    = "data/histo_hww130_std_pu11_randomized.training.root",
 TString bgdInputFile    = "data/histo_w10-ggww-z2-v8-pu11_all_noskim_normalized.root",
 TString outTag          = "default",
 TString myMethodList    = "Fisher,BDT",
 UInt_t  mH              = 130
 )
{
  if(njetsType != 999) return;
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
  mvaVar[ "lep1pt" ]   = 1;
  mvaVar[ "lep2pt" ]   = 1;
  mvaVar[ "lep3pt" ]   = 1;
  mvaVar[ "massZMin" ] = 1;
  mvaVar[ "massMin" ]  = 1;
  mvaVar[ "dRMin" ]    = 1;
  mvaVar[ "met" ]      = 1;
  mvaVar[ "mt1" ]      = 1;
  mvaVar[ "mt2" ]      = 1;
  mvaVar[ "mt3" ]      = 1;
  mvaVar[ "mass3l" ]   = 0;

  TCut sel = "";

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
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

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
  if (mvaVar["lep3pt"])       factory->AddVariable( "lep3pt",        "3rd lepton pt",       "GeV", 'F' );
  if (mvaVar["massZMin"])     factory->AddVariable( "massZMin",      "|M(ll)-mZ|",          "",    'F' );
  if (mvaVar["massMin"])      factory->AddVariable( "massMin",       "M(ll) Min",           "GeV", 'F' );
  if (mvaVar["dRMin"])        factory->AddVariable( "dRMin",         "DeltaR Min",          "",    'F' );
  if (mvaVar["met"])          factory->AddVariable( "met",           "min-MET",             "GeV", 'F' );
  if (mvaVar["mt1"])          factory->AddVariable( "mt1",           "MT(lep1,MET)",        "GeV", 'F' );
  if (mvaVar["mt2"])          factory->AddVariable( "mt2",           "MT(lep2,MET)",        "GeV", 'F' );
  if (mvaVar["mt3"])          factory->AddVariable( "mt3",           "MT(lep3,MET)",        "GeV", 'F' );
  if (mvaVar["mass3l"])       factory->AddVariable( "mass3l",        "MT(lll)",             "GeV", 'F' );

  int nVariablesTemp = 0;

  if (mvaVar["lep1pt"])    { cout << "Adding variable to MVA training: lep1.pt()"      << endl; nVariablesTemp++; }
  if (mvaVar["lep2pt"])    { cout << "Adding variable to MVA training: lep2.pt()"      << endl; nVariablesTemp++; }
  if (mvaVar["lep3pt"])    { cout << "Adding variable to MVA training: lep3.pt()"      << endl; nVariablesTemp++; }
  if (mvaVar["massZMin"])  { cout << "Adding variable to MVA training: mass Z min"     << endl; nVariablesTemp++; }
  if (mvaVar["massMin"])   { cout << "Adding variable to MVA training: dil mass min"   << endl; nVariablesTemp++; }
  if (mvaVar["dRMin"])     { cout << "Adding variable to MVA training: deltaR min"     << endl; nVariablesTemp++; }
  if (mvaVar["met"])	   { cout << "Adding variable to MVA training: min-met"	       << endl; nVariablesTemp++; }
  if (mvaVar["mt1"])	   { cout << "Adding variable to MVA training: mt1"	       << endl; nVariablesTemp++; }
  if (mvaVar["mt2"])	   { cout << "Adding variable to MVA training: mt2"	       << endl; nVariablesTemp++; }
  if (mvaVar["mt3"])	   { cout << "Adding variable to MVA training: mt3"	       << endl; nVariablesTemp++; }
  if (mvaVar["mass3l"])    { cout << "Adding variable to MVA training: trilep mass"    << endl; nVariablesTemp++; }

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

  unsigned int patternTopTag = SmurfTree::TopTag;

  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          event;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  lep3  = 0;
  LorentzVector*  jet1  = 0;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         mt3;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lq3;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;

  signal->SetBranchAddress( "cuts"         , &cuts         );
  signal->SetBranchAddress( "dstype"       , &dstype       );
  signal->SetBranchAddress( "event"        , &event        );
  signal->SetBranchAddress( "scale1fb"     , &scale1fb     );
  signal->SetBranchAddress( "lep1"         , &lep1         );
  signal->SetBranchAddress( "lep2"         , &lep2         );
  signal->SetBranchAddress( "lep3"         , &lep3         );
  signal->SetBranchAddress( "jet1"         , &jet1         );
  signal->SetBranchAddress( "pmet"         , &pmet         );
  signal->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  signal->SetBranchAddress( "mt1"          , &mt1          );
  signal->SetBranchAddress( "mt2"          , &mt2          );
  signal->SetBranchAddress( "mt3"          , &mt3          );
  signal->SetBranchAddress( "lq1"          , &lq1          );
  signal->SetBranchAddress( "lq2"          , &lq2          );
  signal->SetBranchAddress( "lq3"          , &lq3          );
  signal->SetBranchAddress( "lid1"         , &lid1         );
  signal->SetBranchAddress( "lid2"         , &lid2         );
  signal->SetBranchAddress( "lid3"         , &lid3         );

  int nsigtrain = 0;
  int nsigtest  = 0;
  int nbkgtrain = 0;
  int nbkgtest  = 0;

  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);

    //--------------------------------------------------
    // SIGNAL EVENT SELECTION
    //--------------------------------------------------

    if(!((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
         (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && 
         (cuts & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)) continue; // three good leptons

    if( lid3 == 0	                        ) continue;
    if( TMath::Abs(lq1+lq2+lq3) != 1            ) continue;
    if( lep1->pt() <= 20	                ) continue;
    if( lep2->pt() <= 10	                ) continue;
    if( lep3->pt() <= 10	                ) continue;
    if( jet1->pt() >  40	                ) continue;
    if( TMath::Min(pmet,pTrackMet) <= 30        ) continue;
    if( (cuts & patternTopTag) == patternTopTag ) continue;

    double  massZMin = trilepton_info(0,*lep1,*lep2,*lep3,
                                  lq1 ,lq2 ,lq3,
		                  lid1,lid2,lid3,
				  mt1 ,mt2 ,mt3);
    double  massMin  = trilepton_info(1,*lep1,*lep2,*lep3,
                                  lq1 ,lq2 ,lq3,
		                  lid1,lid2,lid3,
				  mt1 ,mt2 ,mt3);
    double  dRMin = trilepton_info(2,*lep1,*lep2,*lep3,
                                  lq1 ,lq2 ,lq3,
		                  lid1,lid2,lid3,
				  mt1 ,mt2 ,mt3);
    if( massMin <= 12	                        ) continue;
    if( massMin >= 100	                        ) continue;
    if( massZMin <= 15	                        ) continue;

    int varCounter = 0;
    
    if (mvaVar["lep1pt"])        vars[varCounter++] = lep1->pt();
    if (mvaVar["lep2pt"])        vars[varCounter++] = lep2->pt();    
    if (mvaVar["lep3pt"])        vars[varCounter++] = lep3->pt();    
    if (mvaVar["massZMin"])      vars[varCounter++] = TMath::Min(massZMin,99.999);
    if (mvaVar["massMin"])       vars[varCounter++] = massMin;
    if (mvaVar["dRMin"])         vars[varCounter++] = dRMin;
    if (mvaVar["met"])	         vars[varCounter++] = TMath::Min(pmet,pTrackMet);
    if (mvaVar["mt1"])	         vars[varCounter++] = mt1;
    if (mvaVar["mt2"])	         vars[varCounter++] = mt2;
    if (mvaVar["mt3"])	         vars[varCounter++] = mt3;
    if (mvaVar["mass3l"])        vars[varCounter++] = (*lep1+*lep2+*lep3).mass();
 
    if ( event%2 != 0 ){
      factory->AddSignalTrainingEvent( vars, scale1fb );
      nsigtrain++;
    }
    else{
      factory->AddSignalTestEvent    ( vars, scale1fb );
      nsigtest++;
    }
  }

  background->SetBranchAddress( "cuts"         , &cuts         );
  background->SetBranchAddress( "dstype"       , &dstype       );
  background->SetBranchAddress( "event"        , &event        );  
  background->SetBranchAddress( "scale1fb"     , &scale1fb     );
  background->SetBranchAddress( "lep1"         , &lep1         );
  background->SetBranchAddress( "lep2"         , &lep2         );
  background->SetBranchAddress( "lep3"         , &lep3         );
  background->SetBranchAddress( "jet1"         , &jet1         );
  background->SetBranchAddress( "pmet"         , &pmet         );
  background->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  background->SetBranchAddress( "mt1"	       , &mt1	       );
  background->SetBranchAddress( "mt2"	       , &mt2	       );
  background->SetBranchAddress( "mt3"	       , &mt3	       );
  background->SetBranchAddress( "lq1"	       , &lq1	       );
  background->SetBranchAddress( "lq2"	       , &lq2	       );
  background->SetBranchAddress( "lq3"	       , &lq3	       );
  background->SetBranchAddress( "lid1"         , &lid1         );
  background->SetBranchAddress( "lid2"         , &lid2         );
  background->SetBranchAddress( "lid3"         , &lid3         );

  cout << "Add background events" << endl;
  cout << "Added " << nsigtrain << " training events" << endl;
  cout << "Added " << nsigtest  << " test events" << endl;
  
  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);

    //--------------------------------------------------
    // BACKGROUND EVENT SELECTION
    //--------------------------------------------------

    if(dstype == SmurfTree::data)  continue; // cut on dstype

    if(!((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
         (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && 
         (cuts & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)) continue; // three good leptons

    if( lid3 == 0	                        ) continue;
    if( TMath::Abs(lq1+lq2+lq3) != 1            ) continue;
    if( lep1->pt() <= 20	                ) continue;
    if( lep2->pt() <= 10	                ) continue;
    if( lep3->pt() <= 10	                ) continue;
    if( jet1->pt() >  40	                ) continue;
    if( TMath::Min(pmet,pTrackMet) <= 30        ) continue;
    if( (cuts & patternTopTag) == patternTopTag ) continue;

    double  massZMin = trilepton_info(0,*lep1,*lep2,*lep3,
                                  lq1 ,lq2 ,lq3,
		                  lid1,lid2,lid3,
				  mt1 ,mt2 ,mt3);
    double  massMin  = trilepton_info(1,*lep1,*lep2,*lep3,
                                  lq1 ,lq2 ,lq3,
		                  lid1,lid2,lid3,
				  mt1 ,mt2 ,mt3);
    double  dRMin = trilepton_info(2,*lep1,*lep2,*lep3,
                                  lq1 ,lq2 ,lq3,
		                  lid1,lid2,lid3,
				  mt1 ,mt2 ,mt3);
    if( massMin <= 12	                        ) continue;
    if( massMin >= 100	                        ) continue;
    if( massZMin <= 15	                        ) continue;

    int varCounter = 0;
    
    if (mvaVar["lep1pt"])        vars[varCounter++] = lep1->pt();
    if (mvaVar["lep2pt"])        vars[varCounter++] = lep2->pt();    
    if (mvaVar["lep3pt"])        vars[varCounter++] = lep3->pt();    
    if (mvaVar["massZMin"])      vars[varCounter++] = TMath::Min(massZMin,99.999);
    if (mvaVar["massMin"])       vars[varCounter++] = massMin;
    if (mvaVar["dRMin"])         vars[varCounter++] = dRMin;
    if (mvaVar["met"])	         vars[varCounter++] = TMath::Min(pmet,pTrackMet);
    if (mvaVar["mt1"])	         vars[varCounter++] = mt1;
    if (mvaVar["mt2"])	         vars[varCounter++] = mt2;
    if (mvaVar["mt3"])	         vars[varCounter++] = mt3;
    if (mvaVar["mass3l"])        vars[varCounter++] = (*lep1+*lep2+*lep3).mass();

    if ( event%2 != 0 ){
      factory->AddBackgroundTrainingEvent( vars, scale1fb );
      nbkgtrain++;
    }
    else{
      factory->AddBackgroundTestEvent    ( vars, scale1fb );
      nbkgtest++;
    }
  }
  
  cout << "Done adding background" << endl;
  cout << "Added " << nbkgtrain << " training events" << endl;
  cout << "Added " << nbkgtest  << " test events" << endl;
  

  // --- end ------------------------------------------------------------
  //
  // --- end of tree registration 
   
  // Set individual event weights (the variables must exist in the original TTree)
  factory->SetSignalWeightExpression    ("scale1fb");
  factory->SetBackgroundWeightExpression("scale1fb");
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
  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                       "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

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
    factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=1000:HiddenLayers=N+1:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

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
			 "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=1000:NNodesMax=5:VarTransform=Decorrelate");

  if (Use["BDTB"]) // Bagging
    factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                         "!H:!V:NTrees=500:BoostType=Bagging:SeparationType=GiniIndex:nCuts=1000:PruneMethod=NoPruning" );

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod(TMVA::Types::kBDT,"BDTD",
                    "!H:!V:NTrees=500:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=1000:PruneMethod=CostComplexity:PruneStrength=50:VarTransform=Decorrelate");

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
