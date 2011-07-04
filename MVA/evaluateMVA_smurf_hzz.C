/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TMath.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void evaluateMVA_smurf_hzz(
char* inputFile      = "data/histo_hww130_std_pu11_randomized.training.root", 
int mH               = 130,
TString myMethodList = "Fisher,BDT",
TString outTag       = "default",
TString path         = "") 
{   
#ifdef __CINT__
  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

  //--------------------------------------------------------------------
  // path to weights dir (this is where MVA training info is stored)
  // output root file will be stored at [path]/output
  //--------------------------------------------------------------------

  cout << "Looking for weights dir at " << path << endl;
  
  vector<char*> samples;
  samples.push_back(inputFile);

  //--------------------------------------------------------------------------------
  // IMPORTANT: set the following variables to the same set used for MVA training!!!
  //--------------------------------------------------------------------------------

  std::map<std::string,int> mvaVar;
  mvaVar[ "lep1pt" ]            = 1;  //pt of leading lepton
  mvaVar[ "lep2pt" ]            = 1;  //pt of sub-leading lepton
  mvaVar[ "dPhi" ]              = 1;  //delta phi btw leptons
  mvaVar[ "dR" ]                = 0;  //delta R btw leptons
  mvaVar[ "dilmass" ]           = 0;  //dilepton mass
  mvaVar[ "type" ]              = 0;  //dilepton flavor type
  mvaVar[ "pmet" ]              = 0;  //projected met
  mvaVar[ "met" ]               = 1;  //met
  mvaVar[ "mt" ]                = 0;  //transverse higgs mass
  mvaVar[ "mt1" ]               = 1;  //transverse mass of leading lepton and met
  mvaVar[ "mt2" ]               = 1;  //transverse mass of sub-leading lepton and met
  mvaVar[ "dPhiLep1MET" ]       = 1;  //delta phi btw leading lepton and met
  mvaVar[ "dPhiLep2MET" ]       = 1;  //delta phi btw leading sub-lepton and met
  mvaVar[ "dPhiDiLepMET" ]	= 0;  //delta phi btw dilepton and met
  mvaVar[ "dPhiDiLepJet1" ]	= 0;  //delta phi btw dilepton and jet1 (only for njet>0)
  mvaVar[ "mtHZZ" ]             = 1;  //tranvserse higgs mass for HZZ

  //---------------------------------------------------------------
  // specifies the selection applied to events in the training
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
  // ---------------------------------------------------------------
  Use["Plugin"]          = 0;
  Use["Category"]        = 0;
  Use["SVM_Gauss"]       = 0;
  Use["SVM_Poly"]        = 0;
  Use["SVM_Lin"]         = 0;

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassificationApplication" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod 
                  << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
          std::cout << it->first << " ";
        }
        std::cout << std::endl;
        return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  const unsigned int nsamples = samples.size();
  
  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    // --- Create the Reader object

    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    //    Float_t var1, var2;
    //    Float_t var3, var4;
    //    reader->AddVariable( "myvar1 := var1+var2", &var1 );
    //    reader->AddVariable( "myvar2 := var1-var2", &var2 );
    //    reader->AddVariable( "var3",                &var3 );
    //    reader->AddVariable( "var4",                &var4 );

    Float_t lep1pt;
    Float_t lep2pt;
    Float_t dPhi;
    Float_t dR;
    Float_t dilmass;
    Float_t type;
    Float_t pmet;
    Float_t met;
    Float_t metPhi;
    Float_t mt;
    Float_t mt1;
    Float_t mt2;
    Float_t dPhiLep1MET;
    Float_t dPhiLep2MET;
    Float_t dPhiDiLepMET;
    Float_t dPhiDiLepJet1;
    Float_t mtHZZ;

    if (mvaVar["lep1pt"])        reader->AddVariable( "lep1pt",        &lep1pt       );
    if (mvaVar["lep2pt"])        reader->AddVariable( "lep2pt",        &lep2pt       );
    if (mvaVar["dPhi"])          reader->AddVariable( "dPhi",          &dPhi         );
    if (mvaVar["dR"])            reader->AddVariable( "dR",            &dR           );
    if (mvaVar["dilmass"])       reader->AddVariable( "dilmass",       &dilmass      );
    if (mvaVar["type"])          reader->AddVariable( "type",          &type         );
    if (mvaVar["pmet"])          reader->AddVariable( "pmet",          &pmet         );
    if (mvaVar["met"])           reader->AddVariable( "met",           &met          );
    if (mvaVar["metPhi"])        reader->AddVariable( "metPhi",        &metPhi       );
    if (mvaVar["mt"])            reader->AddVariable( "mt",            &mt           );
    if (mvaVar["mt1"])           reader->AddVariable( "mt1",           &mt1          );
    if (mvaVar["mt2"])           reader->AddVariable( "mt2",           &mt2          );
    if (mvaVar["dPhiLep1MET"])   reader->AddVariable( "dPhiLep1MET",   &dPhiLep1MET  );
    if (mvaVar["dPhiLep2MET"])   reader->AddVariable( "dPhiLep2MET",   &dPhiLep2MET  );
    if (mvaVar["dPhiDiLepMET"])  reader->AddVariable( "dPhiDiLepMET",  &dPhiDiLepMET );
    if (mvaVar["dPhiDiLepJet1"]) reader->AddVariable( "dPhiDiLepJet1", &dPhiDiLepJet1);
    if (mvaVar["mtHZZ"])         reader->AddVariable( "mtHZZ",         &mtHZZ        );
 
    // Spectator variables declared in the training have to be added to the reader, too
    //    Float_t spec1,spec2;
    //    reader->AddSpectator( "spec1 := var1*2",   &spec1 );
    //    reader->AddSpectator( "spec2 := var1*3",   &spec2 );

    //Float_t Category_cat1, Category_cat2, Category_cat3;
    if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      //       reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      //       reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      //       reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
    }

    // --- Book the MVA methods

    //--------------------------------------------------------------------------------------
    // tell evaluateMVA_smurf_hzz where to find the weights dir, which contains the trained MVA's. 
    // In this example, the weights dir is located at [path]/[dir]
    // and the output root file is written to [path]/[output]
    //--------------------------------------------------------------------------------------

    TString dir    = path + "weights/";
    TString outdir = path + "output/";
    TString prefix = outTag;

    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
        TString methodName = TString(it->first) + TString(" method");
        TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
        reader->BookMVA( methodName, weightfile ); 
      }
    }
   
    // Book output histograms
    UInt_t nbin = 1000;
    TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
    TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
    TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
    TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
    TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

    if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );               
    if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
    if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
    if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
    if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
    if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
    if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
    if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
    if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
    if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
    if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
    if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
    if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
    if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
    if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
    if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
    if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
    if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
    if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
    if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -1. , 1. );
    if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
    if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
    if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
    if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
    if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
    if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
    if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
    if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
    if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
    if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

    if (Use["Likelihood"])    histLk      ->Sumw2();
    if (Use["LikelihoodD"])   histLkD     ->Sumw2();
    if (Use["LikelihoodPCA"]) histLkPCA   ->Sumw2();
    if (Use["LikelihoodKDE"]) histLkKDE   ->Sumw2();
    if (Use["LikelihoodMIX"]) histLkMIX   ->Sumw2();
    if (Use["PDERS"])         histPD      ->Sumw2();
    if (Use["PDERSD"])        histPDD     ->Sumw2();
    if (Use["PDERSPCA"])      histPDPCA   ->Sumw2();
    if (Use["KNN"])           histKNN     ->Sumw2();
    if (Use["HMatrix"])       histHm      ->Sumw2();
    if (Use["Fisher"])        histFi      ->Sumw2();
    if (Use["FisherG"])       histFiG     ->Sumw2();
    if (Use["BoostedFisher"]) histFiB     ->Sumw2();
    if (Use["LD"])            histLD      ->Sumw2();
    if (Use["MLP"])           histNn      ->Sumw2();
    if (Use["MLPBFGS"])       histNnbfgs  ->Sumw2();
    if (Use["MLPBNN"])        histNnbnn   ->Sumw2();
    if (Use["CFMlpANN"])      histNnC     ->Sumw2();
    if (Use["TMlpANN"])       histNnT     ->Sumw2();
    if (Use["BDT"])           histBdt     ->Sumw2();
    if (Use["BDTD"])          histBdtD    ->Sumw2();
    if (Use["BDTG"])          histBdtG    ->Sumw2();
    if (Use["RuleFit"])       histRf      ->Sumw2();
    if (Use["SVM_Gauss"])     histSVMG    ->Sumw2();
    if (Use["SVM_Poly"])      histSVMP    ->Sumw2();
    if (Use["SVM_Lin"])       histSVML    ->Sumw2();
    if (Use["FDA_MT"])        histFDAMT   ->Sumw2();
    if (Use["FDA_GA"])        histFDAGA   ->Sumw2();
    if (Use["Category"])      histCat     ->Sumw2();
    if (Use["Plugin"])        histPBdt    ->Sumw2();

    // PDEFoam also returns per-event error, fill in histogram, and also fill significance
    if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
    }

    // Book example histogram for probability (the other methods are done similarly)
    TH1F *probHistFi(0), *rarityHistFi(0);
    if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
    }

    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //   

    TChain *ch = new TChain("tree");
    ch -> Add( Form("%s" , samples.at(i)) );

    TFile *f = TFile::Open( Form("%s" , samples.at(i)) );
    assert(f);
    TTree* t = (TTree*)f->Get("tree");
    assert(t);
    t->SetBranchStatus("*", 1);

    TString ofn(samples.at(i));
    ofn.ReplaceAll("data/","");
    const char* mydir = outdir;
    TFile *out = TFile::Open( Form("%s/%s_%s" , mydir, outTag.Data(), ofn.Data() ) ,"RECREATE" );
    TTree *clone;
    clone = t->CloneTree(-1, "fast");
   
    Float_t bdt;
    Float_t bdtd;
    Float_t nn;
    Float_t knn;
    Float_t bdtg;
    Int_t   test;

    TBranch* br_bdt        = 0;
    TBranch* br_bdtd       = 0;
    TBranch* br_nn         = 0;
    TBranch* br_knn        = 0;
    TBranch* br_bdtg       = 0;
    TBranch* br_test = clone->Branch(Form("test_hww%i_ww",mH) , &test , Form("test_hww%i_ww/I" , mH) );

    if(Use["BDT"])         br_bdt        = clone->Branch(Form("bdt_hww%i_ww"        ,mH) , &bdt        , Form("bdt_hww%i_ww/F"        ,mH) );
    if(Use["BDTD"])        br_bdtd       = clone->Branch(Form("bdtd_hww%i_ww"       ,mH) , &bdtd       , Form("bdtd_hww%i_ww/F"       ,mH) );
    if(Use["MLPBNN"])      br_nn         = clone->Branch(Form("nn_hww%i_ww"         ,mH) , &nn         , Form("nn_hww%i_ww/F"         ,mH) );
    if(Use["KNN"])         br_knn        = clone->Branch(Form("knn_hww%i_ww"        ,mH) , &knn        , Form("knn_hww%i_ww/F"        ,mH) );
    if(Use["BDTG"])        br_bdtg       = clone->Branch(Form("bdtg_hww%i_ww"       ,mH) , &bdtg       , Form("bdtg_hww%i_ww/F"       ,mH) );

    if(Use["BDT"])         br_bdt       -> SetTitle(Form("BDT Output H%i"    , mH));
    if(Use["BDTD"])        br_bdt       -> SetTitle(Form("BDTD Output H%i"   , mH));
    if(Use["MLPBNN"])      br_nn        -> SetTitle(Form("MLPBNN Output H%i" , mH));
    if(Use["KNN"])         br_knn       -> SetTitle(Form("KNN Output H%i"    , mH));
    if(Use["BDTG"])        br_bdtg      -> SetTitle(Form("BDTG Output H%i"   , mH));

    // --- Event loop

    // Prepare the event tree
    // - here the variable names have to corresponds to your tree
    // - you can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
  
    TTree *theTree     = (TTree*) ch;

    std::cout << "--- Using input files: -------------------" <<  std::endl;

    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TChainElement* currentFile = 0;
    
    while((currentFile = (TChainElement*)fileIter.Next())) {
      std::cout << currentFile->GetTitle() << std::endl;
    }

    UInt_t          event_;
    Float_t         scale1fb_;
    LorentzVector*  lep1_   = 0;
    LorentzVector*  lep2_   = 0;
    Float_t         dPhi_;
    Float_t         dR_;
    LorentzVector*  dilep_  = 0;
    UInt_t          type_;
    Float_t         pmet_;
    Float_t         met_;
    Float_t         metPhi_;
    Float_t         mt_;
    Float_t         mt1_;
    Float_t         mt2_;
    Float_t         dPhiLep1MET_;
    Float_t         dPhiLep2MET_;
    Float_t         dPhiDiLepMET_;
    Float_t         dPhiDiLepJet1_;

    theTree->SetBranchAddress( "event"        , &event_        );
    theTree->SetBranchAddress( "scale1fb"     , &scale1fb_     );
    theTree->SetBranchAddress( "lep1"         , &lep1_         );
    theTree->SetBranchAddress( "lep2"         , &lep2_         );
    theTree->SetBranchAddress( "dPhi"         , &dPhi_         );
    theTree->SetBranchAddress( "dR"           , &dR_           );
    theTree->SetBranchAddress( "dilep"        , &dilep_        );
    theTree->SetBranchAddress( "type"         , &type_         );
    theTree->SetBranchAddress( "pmet"         , &pmet_         );
    theTree->SetBranchAddress( "met"          , &met_          );
    theTree->SetBranchAddress( "metPhi"       , &metPhi_       );
    theTree->SetBranchAddress( "mt"           , &mt_           );
    theTree->SetBranchAddress( "mt1"          , &mt1_          );
    theTree->SetBranchAddress( "mt2"          , &mt2_          );
    theTree->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET_  );
    theTree->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET_  );
    theTree->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET_ );
    theTree->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1_);

    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;

    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    int npass   = 0;
    float yield = 0.;

    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      //--------------------------------------------------------
      // important: here we associate branches to MVA variables
      //--------------------------------------------------------

      lep1pt         = lep1_->pt();
      lep2pt         = lep2_->pt();
      dPhi           = dPhi_;
      dR             = dR_;
      dilmass        = dilep_->mass();
      type           = type_;
      pmet           = pmet_;
      met            = met_;
      mt             = mt_;
      mt1            = mt1_;
      mt2            = mt2_;
      dPhiLep1MET    = dPhiLep1MET_;
      dPhiLep2MET    = dPhiLep2MET_;
      dPhiDiLepMET   = dPhiDiLepMET_;
      dPhiDiLepJet1  = dPhiDiLepJet1_;

      double pxzll = dilep_->px() + met*cos(metPhi);
      double pyzll = dilep_->py() + met*sin(metPhi);
      mtHZZ = TMath::Power(sqrt(dilep_->pt()*dilep_->pt()+91.1876*91.1876)+
                           sqrt(met        *met          +91.1876*91.1876),2)
		    -pxzll*pxzll-pyzll*pyzll;
      if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;

      npass++;
      yield+=scale1fb_;

      // --- Return the MVA outputs and fill into histograms

      if (Use["BDT"]){
        bdt = reader->EvaluateMVA( "BDT method" );
        br_bdt->Fill();
      }
      if (Use["BDTD"]){
        bdtd = reader->EvaluateMVA( "BDTD method" );
        br_bdtd->Fill();
      }
      if (Use["MLPBNN"]){
        nn  = reader->EvaluateMVA( "MLPBNN method" );
        br_nn->Fill();
      }
      if (Use["KNN"]){
        knn  = reader->EvaluateMVA( "KNN method" );
        br_knn->Fill();
      }
      if (Use["BDTG"]){
        bdtg  = reader->EvaluateMVA( "BDTG method" );
        br_bdtg->Fill();
      }

      test = 0;
      if(event_ %2 == 1 ) test = 1; 
      br_test->Fill();


      if (Use["CutsGA"]) {
        // Cuts is a special case: give the desired signal efficienciy
        Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
        if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) , scale1fb_);
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) , scale1fb_);
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) , scale1fb_);
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) , scale1fb_);
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) , scale1fb_);
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) , scale1fb_);
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) , scale1fb_);
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) , scale1fb_);
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) , scale1fb_);
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) , scale1fb_);
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) , scale1fb_);
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) , scale1fb_);
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) , scale1fb_);
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) , scale1fb_);
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) , scale1fb_);
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) , scale1fb_);
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) , scale1fb_);
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) , scale1fb_);
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) , scale1fb_);
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) , scale1fb_);
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) , scale1fb_);
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) , scale1fb_);
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) , scale1fb_);
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) , scale1fb_);
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) , scale1fb_);
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) , scale1fb_);
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) , scale1fb_);
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) , scale1fb_);
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) , scale1fb_);
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) , scale1fb_);

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
        Double_t val = reader->EvaluateMVA( "PDEFoam method" );
        Double_t err = reader->GetMVAError();
        histPDEFoam   ->Fill( val );
        histPDEFoamErr->Fill( err );         
        if (err>1.e-50) histPDEFoamSig->Fill( val/err , scale1fb_);
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
        probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) , scale1fb_ );
        rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) , scale1fb_ );
      }
    }

    std::cout << npass << " events passing selection, yield " << yield << std::endl;
 
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    // Get efficiency for cuts classifier
    if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                 << " (for a required signal efficiency of " << effS << ")" << std::endl;

    if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
        std::vector<Double_t> cutsMin;
        std::vector<Double_t> cutsMax;
        mcuts->GetCuts( 0.7, cutsMin, cutsMax );
        std::cout << "--- -------------------------------------------------------------" << std::endl;
        std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
        for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
          std::cout << "... Cut: " 
                    << cutsMin[ivar] 
                    << " < \"" 
                    << mcuts->GetInputVar(ivar)
                    << "\" <= " 
                    << cutsMax[ivar] << std::endl;
        }
        std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
    }

    // --- write output baby

    clone->Write(); 
    out->Close();
    f->Close();

    // --- Write histograms
    cout << "dir " << dir << endl;

    TFile *target  =   new TFile( Form("%s/%s_histos_%s",mydir, outTag.Data(), ofn.Data() ) ,"RECREATE" );
    cout << "Writing to file " << Form("%s/%s_histos_%s",mydir, outTag.Data(), ofn.Data() ) << endl;

    if (Use["Likelihood"   ])   histLk     ->Write();
    if (Use["LikelihoodD"  ])   histLkD    ->Write();
    if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
    if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
    if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
    if (Use["PDERS"        ])   histPD     ->Write();
    if (Use["PDERSD"       ])   histPDD    ->Write();
    if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
    if (Use["KNN"          ])   histKNN    ->Write();
    if (Use["HMatrix"      ])   histHm     ->Write();
    if (Use["Fisher"       ])   histFi     ->Write();
    if (Use["FisherG"      ])   histFiG    ->Write();
    if (Use["BoostedFisher"])   histFiB    ->Write();
    if (Use["LD"           ])   histLD     ->Write();
    if (Use["MLP"          ])   histNn     ->Write();
    if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
    if (Use["MLPBNN"       ])   histNnbnn  ->Write();
    if (Use["CFMlpANN"     ])   histNnC    ->Write();
    if (Use["TMlpANN"      ])   histNnT    ->Write();
    if (Use["BDT"          ])   histBdt    ->Write();
    if (Use["BDTD"         ])   histBdtD   ->Write();
    if (Use["BDTG"         ])   histBdtG   ->Write(); 
    if (Use["RuleFit"      ])   histRf     ->Write();
    if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
    if (Use["SVM_Poly"     ])   histSVMP   ->Write();
    if (Use["SVM_Lin"      ])   histSVML   ->Write();
    if (Use["FDA_MT"       ])   histFDAMT  ->Write();
    if (Use["FDA_GA"       ])   histFDAGA  ->Write();
    if (Use["Category"     ])   histCat    ->Write();
    if (Use["Plugin"       ])   histPBdt   ->Write();

    // Write also error and significance histos
    if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

    // Write also probability hists
    if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
    target->Close();

    delete reader;
    
    std::cout << "==> TMVAClassificationApplication is done with sample " << samples.at(i) << endl << std::endl;
  } 
}

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom.h"

using namespace std;


