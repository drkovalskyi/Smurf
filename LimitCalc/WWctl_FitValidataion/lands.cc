#include <TString.h>
#include <sstream>
#include <map>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <TSystem.h>
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "CRandom.h"
#include "UtilsROOT.h"
#include "PlotUtilities.h"
#include "CLsLimit.h"
#include "SignificanceBands.h"
#include "MLLxsBands.h"
#include "RooRandom.h"
#include "TStopwatch.h"
#include "RooDataSet.h"
#include "TFile.h"
#include "RooMsgService.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "RooEllipse.h"
#include "TMatrixD.h"
#include "TH2D.h"
#include "RooMultiVarGaussian.h"

#ifdef LANDS_PROFILE
#include <google/profiler.h>
#endif

using namespace std;
using namespace lands;
using namespace RooFit;

typedef std::map<TString, vector<TString> > TStrMap;
typedef std::pair<TString, vector<TString> > TStrPair;
TStrMap options;
vector<double> GetListToEval(TString parname, vector<TString> tmpv);
vector<double> GetListToEvalForNuis(TString& parname, vector<TString> tmpv);
void processParameters(int argc, const char* argv[]);
void PrintHelpMessage();
double Err_Bisection(CountingModel*cms, TString spar, double * bestFitPars, double fmin, double bestfitPOI, int hilo);

// parameters :
vector<TString> datacards; // enable multiple cards combination, and allow "* ?" in the input
TString method; // ProfiledLikelihood, Bayesian, Hybrid, FeldmanCousins 
TString jobname; // default constructed as "fist_card_name"+"method"
TString dataset; // default: data_obs;  asimov_b, asimov_sb
int debug; // default 0
bool doExpectation; // default false; 
int toys; // number of toys to do expectation, e.g. bands and mean/median limits.   default 1000,  you must also turn on doExpectation
int toysHybrid; // number of toys used to build CLsb, CLb, CLs, as well as number of toys per point in FeldmanCousins construction
int toysBayesian; // number of toys used to average out nuisance paramereters in Bayesian method 
int toysPreBayesian; // number of toys used to average out nuisance paramereters in Bayesian method, for pre estimation 
int seed;  // default 1234;   seed of random number generator engine
int systematics; // default 1, will determin if use systematics according to data card; if 0, then will not use systematics 
float CL;  // default 0.95;  Confidence Level...
lands::PRIOR prior; //default flat;   prior on signal.  other options:  corr,  1/sqrt(r)
int calcsignificance; // default 0;  if 1, then calc significance instead of limit
bool lowerLimit; // in FeldmanCousins,  0 by default for upper limit,  if 1 then calc lower limit  
int toysFactor; // in FeldmanCousins, Increase the toys per point by this factor w.r.t. the minimum from adaptive sampling,   default =1
bool adaptiveSampling; // currently only implemented in FeldmanCousins,  turn on (=1) by default.  =0 off
int testStat; //default LEP.  other options: TEV, Atlas, AtlasAllowMuHatNeg ...
int rule; //default CLs.  other options: CLsb
float clsAcc; // 0.001; Absolute accuracy on CLs to reach to terminate the scan
double rAbsAcc;// 0.01; Absolute accuracy on r to reach to terminate the scan
double rRelAcc;// 0.01; Relative accuracy on r to reach to terminate the scan
bool singlePoint;
double testR;
int scanRs;
int scanMs;
int nSteps;
int nMsteps;
int bPlots = 0;
int tossToyConvention;
int tossPseudoDataConvention;
int UseBestEstimateToCalcQ;

// for FeldmanCousins
bool bQuickEstimateInitialLimit = true; 
double initialRmin = 1., initialRmax = 21;// only when bQuickEstimateInitialLimit==false

int oneside = 2; //for PLR limit     for default
TString PLalgorithm; // algorithm for PL approximation method (    Minos,    Migrad )

bool bReadPars = false;
bool bWritePars = false;
TString fileFittedPars="";
bool bNotCalcQdata = false;
int nToysForCLb = -1;
int nToysForCLsb = -1;
bool bNotCalcCLssbb = false;
bool bSaveM2lnQ = false;
bool bSkipM2lnQ = false; 
bool bOnlyEvaluateQdata= false; 

bool bM2lnQGridPreComputed = false;
TString sFileM2lnQGrid = "";


bool bCalcObservedLimit = true;

TString sFileLimitsDistribution = "";

vector<TString> librariesToBeLoaded;

bool bOnlyEvalCL_forVR = false;
vector<double> vR_toEval;
vector<double> vM_toEval;

int makeCLsBands = 0;

TString sFileBonlyM2lnQ = "";

double HiggsMass = -1;

int nTries = 1;

int maximumFunctionCallsInAFit = 5000000;
int minuitSTRATEGY = 0;
int minuitPrintLevel= -1;
double minuitTolerance = 0.001;
double nuisancesRange = 5.;

double flatPriorCropNsigma = 3;

TString sPhysicsModel = "StandardModelHiggs";

bool bRunOnlyWithBestFittedNuisances_bayesian = false;
double inputMu_bayesian = 0;

bool bWriteToys = false;

TString ExpectationHints = "";
bool runAsymptoticLimits();
bool runAsymptoticCLs();
double runProfileLikelihoodApproximation(double neg2_llr=-9999);
bool runBayesian();
bool saveExpectation();
bool runMaxLikelihoodFit(bool printInfo, bool makePlots, TString smakePlots="");
bool runPValue(double &pvalue, double& significance, double& muhat, bool printInfo=false);

double preHints_obs=0, preHints_m2s=0, preHints_m1s=0, preHints_p1s=0, preHints_p2s=0, preHints_median=0;
// common results
double rmean;
vector<double> difflimits;


bool bConstrainMuPositive = false;  // in the max_ll fit of ScanningMuFit
bool bAlsoExtract2SigmErrorBar = false;

// for LHC-CLs method
bool bFixNuisancsAtNominal = false;

// custom, constrain mu to be [ rMin, rMax ]  during fit
double customRMin=0, customRMax=0;

bool bForceSymmetryError = false;

bool bMultiSigProcShareSamePDF = false;

TString scanRAroundBand = "all";


bool bConstructAsimovbFromNominal = true; 
bool bConstructAsimovsbFromFit = false;

bool bRedefineObservableRange = false;
double ObservableRangeMin = 0;
double ObservableRangeMax = 0;
int ObservableBins = 100;

// for parametric shapes in RooFit
vector<TString> rebin_observables ;
vector<int>     rebin_observables_nbins;
vector<double>  rebin_observables_xmin;
vector<double>  rebin_observables_xmax;

bool bDumpFitResults = false;
TString FileNuisanceFeeding = ""; // FIXME jaehyeok
bool    ManualNuisanceFeeding = false; // FIXME jaehyeok
int maxsets_forcaching = 0;
int PrintParameterFrom = -1;
int PrintParameterTo = -1;
bool bRunMinuitContour = false;

vector<structPOI> addtionalPOIs;

vector<TString> vCouplingsDef;

TString sSMBrFile;
vector<TString> CVsetting, CFsetting;
bool scanCvCf=false;
vector<double> vCv_toEval;
vector<double> vCf_toEval;
TString sbinningCV, sbinningCF;
TString sbinningMH, sbinningMU;

double InjectingSignalStrength = 1;
int ToysAtDifferentSignalStrength = 0; // for MLLxsBands toys // b-only,    1 to use injectingSignalStrength,  2 to use best fitted signal strength	


vector<TString> CvvSetting, CggSetting, CttSetting, CbbSetting, CglglSetting;
bool scanCX=false;
vector<double> vCvv_toEval, vCgg_toEval, vCtt_toEval, vCbb_toEval, vCglgl_toEval;

TString ErrEstAlgo = "Minos";
double PrintPdfEvlCycle = -1;

bool doMemoryCheck = false;

bool DoNotRunGlobalFit = false;

int NumDegreeOfFreedom = 1; // for 68% CL extraction in  n dimensions  

bool GenerateToysAtBestFitSB = false;

bool NoErrorEstimate = false;

//more generic option
TString mainPOI="";
bool scanPOI=false;
vector<double> scanPOIrange;
TString sbinningPOI;

bool bFixNuisancesAtBestFit = false; // when doing scanning mass without systematics --> if remove all sys from cards --> bias,  --> so need to fit to data and then fix to that 

int intRemoveBins = -1; // only remove// 0: remove only bins with total bkg<=0,  1: remove bins with total sig<=0,  2: both;   <0 don't remove anything 

bool useHist = false;

bool newExpectedAsymBand = false; // for new expected asymptotic bands, see https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/297.html

TString sToysFromFile = "";

TString sMinimizingApproach = "Normal";


// 
double fitBestMu, fitBestMuErrHi, fitBestMuErrLo, fitBestMH, fitBestMHErrLo, fitBestMHErrHi;
vector<int> viPseudoPOIs; 
vector<TString> vsPseudoPOIs; 
vector<double> vdPseudoPOIs; 

vector<TString> vsParToBeFixedWhenGenToys; 
vector<double>  vdParToBeFixedWhenGenToys; 

TString scanNuisance = "";
vector<double> vNuis_toEval;

int main(int argc, const char*argv[]){
    processParameters(argc, argv);


    if(debug<2)RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
    if(debug<10)RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
    if(debug<10)RooMsgService::instance().getStream(1).removeTopic(Caching) ;
    if(debug<2)RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    TString slib = gSystem->GetFromPipe("ls $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so");
    if(!slib.BeginsWith("ls")) gSystem->Load(slib);
    for(int i=0; i<librariesToBeLoaded.size(); i++){
        gSystem->Load(librariesToBeLoaded[i]);
    }

    CountingModel *cms; // this instance will be the one combining all sub cards
    cms = new CountingModel();
    cms->SetDebug(debug);
    cms->ForceSymmetryError(bForceSymmetryError);
    if(sPhysicsModel=="ChargedHiggs") cms->SetPhysicsModel(typeChargedHiggs);
    if(bRedefineObservableRange) cms->SetObservableRange(ObservableBins, ObservableRangeMin, ObservableRangeMax);

    /*
     * combining at most 100 datacards
     */
    CountingModel *tmp1[100];
    CountingModel *tmp[100];// = new CountingModel(); 
    if(datacards.size()>100) {cout<<"too many datacards "<<datacards.size()<<endl; exit(0);}
    for(int i=0; i<datacards.size(); i++){
        tmp[i] = new CountingModel();
        tmp[i] -> SetDebug(debug);
        tmp[i] -> ForceSymmetryError(bForceSymmetryError);
        if(bRedefineObservableRange) tmp[i]->SetObservableRange(ObservableBins, ObservableRangeMin, ObservableRangeMax);
        if(sPhysicsModel=="ChargedHiggs") tmp[i]->SetPhysicsModel(typeChargedHiggs);
        ConfigureModel(tmp[i], HiggsMass, datacards[i].Data(), useHist, debug);
        tmp[i]->SetUseSystematicErrors(true);
    }
    if(debug)cout<<"totally "<<datacards.size()<<" data cards processed"<<endl;
    if(datacards.size()==1) cms=tmp[0];
    else if(datacards.size()>=2){
        tmp1[1] = CombineModels(tmp[0], tmp[1]);
        tmp1[1]->SetUseSystematicErrors(true);
        if(debug)cout<<"2 data cards have been combined"<<endl;
        for(int i=2; i<datacards.size(); i++){
            if(debug)cout<<"Adding "<<i+1<<" data cards have been combined"<<endl;
            tmp1[i] = CombineModels(tmp1[i-1], tmp[i]);
            tmp1[i]->SetDebug(debug);
            if(sPhysicsModel=="ChargedHiggs") tmp1[i]->SetPhysicsModel(typeChargedHiggs);
            tmp1[i]->SetUseSystematicErrors(true);
            if(debug)cout<<i+1<<" data cards have been combined"<<endl;
        }	
        cms = tmp1[datacards.size()-1];
    }else{
        if(sFileLimitsDistribution=="" && makeCLsBands<2 && sFileBonlyM2lnQ=="") exit(0);
    }
    cout<<"totally "<<datacards.size()<<" data cards combined"<<endl;

    for(int i=0; i<addtionalPOIs.size(); i++){
        cms->AddFlatParam(addtionalPOIs[i].name.Data(), addtionalPOIs[i].value, addtionalPOIs[i].minV, addtionalPOIs[i].maxV );
    }

    for(int i=0; i<vCouplingsDef.size(); i++){
        cms->AddCouplingParameter(vCouplingsDef[i]);
    }


    RooArgSet ras = cms->GetWorkSpaceVaried()->allVars();
    std::auto_ptr<TIterator> itertmp(ras.createIterator());
    for(RooRealVar * obs= (RooRealVar*)itertmp->Next(); obs!=0; obs=(RooRealVar*)itertmp->Next()){
        if(obs->isConstant()) continue;
        else { 
            TString sobs = obs->GetName();

            bool alreayIn = false;
            for(int i=0; i<cms->Get_v_uncname().size(); i++) {
                TString sunc = cms->Get_v_uncname()[i];				
                if(sobs==sunc || sobs==(sunc+"_x") || sobs==(sunc+"_mean"))
                { alreayIn=true; break;}
            }
            if(!alreayIn){
                double errlo = obs->getErrorLo();
                double errhi = obs->getErrorHi();
                double xmin=0, xmax=0;
                errlo=fabs(errlo);
                errhi=fabs(errhi);

                if( ( obs->getVal() - errlo*7 > obs->getMin()) and errlo>0 ) xmin = obs->getVal()-errlo*7;
                else xmin = obs->getMin();
                if( ( obs->getVal() + errhi*7 < obs->getMax()) and errhi>0 ) xmax = obs->getVal()+errhi*7;
                else xmax = obs->getMax();
                //cms->AddFlatParam(obs->GetName(), obs->getVal(), xmin, xmax);
                //cout<<"Adding flatParam ******** Non-constant variable: "<<obs->GetName()<<" "<<obs->getVal()<<" ["<<obs->getMin()<<","<<obs->getMax()<<"] --> ["<<xmin<<","<<xmax<<"]"<<endl;
                cout<<" *Non-constant variable: "<<obs->GetName()<<"  flatParam             "<<endl;//<<obs->getVal()<<" ["<<obs->getMin()<<","<<obs->getMax()<<"] --> ["<<xmin<<","<<xmax<<"]"<<endl;
            }
        }
    }


    vector<structPOI> vpoi_cx;
    if(sPhysicsModel=="ChargedHiggs") cms->SetPhysicsModel(typeChargedHiggs);
    if(sPhysicsModel=="C5Higgs") {cms->SetPhysicsModel(typeC5Higgs); vpoi_cx=cms->AddCX(CvvSetting, CggSetting, CttSetting, CbbSetting, CglglSetting); vCouplingsDef.push_back("Cvv"); vCouplingsDef.push_back("Cgg"); vCouplingsDef.push_back("Ctt"); vCouplingsDef.push_back("Cbb"); vCouplingsDef.push_back("Cglgl");}
    if(sPhysicsModel=="CvCfHiggs") {cms->SetPhysicsModel(typeCvCfHiggs); vpoi_cx= cms->AddCvCf( CVsetting, CFsetting); vCouplingsDef.push_back("CV"); vCouplingsDef.push_back("CF");}
    if(sPhysicsModel=="CvCfHiggs" or sPhysicsModel=="C5Higgs") { cms->GetSMHiggsBuilder()->readSMBr(sSMBrFile);}	
    cms->SetUseSystematicErrors(systematics);

    cms->Set_printPdfEvlCycle(PrintPdfEvlCycle);

    if(vCouplingsDef.size()>0){
        cms->CheckCouplingSet();
    }
    // done combination

    cms->MultiSigProcShareSamePDF(bMultiSigProcShareSamePDF);

    CRandom *rdm = new CRandom(seed);  //initilize a random generator
    if(seed==0) cout<<" TrueRandom Seed = "<<rdm->GetSeed()<<endl;
    RooRandom::randomGenerator()->SetSeed(rdm->GetSeed());

    // common operations
    if(debug)cms->Print();
    cms->SetRdm(rdm);
    //cms->RemoveChannelsWithExpectedSignal0orBkg0();
    int nch_removed  = cms->RemoveChannelsWithExpectedSignal0orBkg0(intRemoveBins); // 0: remove only bins with total bkg<=0,  1: remove bins with total sig<=0,  2: both
    // FIXME need add an option for speficying if to remove some channels 
    if(debug and nch_removed )cms->Print();
    if(nch_removed)cms->SetUseSystematicErrors(systematics);//need reconfig,  otherwise crash
    //cms->SetAllowNegativeSignalStrength(false);

    cms->SetMoveUpShapeUncertainties(1);


    cms->SetMass(HiggsMass);

    cms->SetMinimzingApproach(sMinimizingApproach);


    if(debug) cms->GetWorkSpace()->Print();
    //**********************  setValueDirty() takes a lot of timing ************************
    //RooAbsArg::setDirtyInhibit(0);

    cms->PrintParametricChannelDataEntries();

    TStopwatch watch;  
    watch.Start();

#ifdef LANDS_PROFILE
    if (gSystem->Getenv("CPUPROFILE"))
    {
        printf("Starting profiling into '%s'\n", gSystem->Getenv("CPUPROFILE"));
        ProfilerStart(gSystem->Getenv("CPUPROFILE"));
    }
#endif


    cms_global= cms;
    cms_global->SetDebug(debug);
    cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);
    cms_global->Set_minuitPrintLevel(minuitPrintLevel);
    cms_global->Set_minuitTolerance(minuitTolerance);
    cms_global->Set_nuisancesRange(nuisancesRange);
    cms_global->Set_maximumFunctionCallsInAFit(maximumFunctionCallsInAFit);

    cms_global->Set_maxSetsCaching(maxsets_forcaching);
    cms_global->Set_PrintParameter(PrintParameterFrom, PrintParameterTo);

    if(rebin_observables.size()>0){
        for(int i=0; i<rebin_observables.size(); i++){
            if(cms_global->GetWorkSpace()->var(rebin_observables[i])) {
                RooRealVar *x = cms_global->GetWorkSpace()->var(rebin_observables[i]);
                x->setBins(rebin_observables_nbins[i]);
                if(rebin_observables_xmin[i]<rebin_observables_xmax[i]){
                    x->setMin(rebin_observables_xmin[i]);
                    x->setMax(rebin_observables_xmax[i]);
                }

                x = cms_global->GetWorkSpaceVaried()->var(rebin_observables[i]);
                x->setBins(rebin_observables_nbins[i]);
                if(rebin_observables_xmin[i]<rebin_observables_xmax[i]){
                    x->setMin(rebin_observables_xmin[i]);
                    x->setMax(rebin_observables_xmax[i]);
                }


                for(int j=0; j<cms->Get_v_pdfs_roodataset_real().size(); j++){
                    RooArgSet * observable = (RooArgSet*)cms->Get_v_pdfs_roodataset_real()[j]->get();
                    std::auto_ptr<TIterator> iter(observable->createIterator());
                    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
                        if(TString(obs->GetName()) == rebin_observables[i]){
                            obs->setBins(rebin_observables_nbins[i]);
                            if(rebin_observables_xmin[i]<rebin_observables_xmax[i]){
                                obs->setMin(rebin_observables_xmin[i]);
                                obs->setMax(rebin_observables_xmax[i]);
                            }
                        }
                    }
                }
            }else{
                cout<< " WARNING potential error ***** : --RebinObservables has a observable \""<<rebin_observables[i]<<"\" which is not in workspace "<<endl;
            }
        }
    }


    structPOI poiMU("signal_strength", 0, 0, 0, 0, 0);
    structPOI poiM("MH", 0, 0, 0, 0, 0);
    cms_global->addPOI(poiMU);
    if(cms_global->Get_MH_i()>0)cms_global->addPOI(poiM);
    for(int i=0; i<addtionalPOIs.size(); i++) cms_global->addPOI(addtionalPOIs[i]);
    if(sPhysicsModel=="CvCfHiggs" or sPhysicsModel=="C5Higgs"){
        for(int i=0; i<vpoi_cx.size(); i++) {
            cms_global->addPOI(vpoi_cx[i]);
        }
    } 

    TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
    if(customRMax!=customRMin) {_customRMin=customRMin; _customRMax = customRMax;}

    _startNuisances = cms->Get_norminalPars();
    _inputNuisances = cms->Get_norminalPars();
    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
    //bConstructAsimovbFromNominal = true;
    //if(method == "Asymptotic") bConstructAsimovbFromNominal = false;
    if(dataset=="asimov_b" or method=="Asymptotic" or ExpectationHints=="Asymptotic")cms->ConstructAsimovData(0, bConstructAsimovbFromNominal); // b-only asimov 
    if(dataset=="asimov_sb")cms->ConstructAsimovData(1, !bConstructAsimovsbFromFit, InjectingSignalStrength); // sb asimov
    if(dataset == "asimov_b")cms->UseAsimovData(0);
    else if(dataset == "asimov_sb")cms->UseAsimovData(1);

    cms->Set_minuitPrintLevel(minuitPrintLevel);

    cms->SetErrEstAlgo(ErrEstAlgo);

    if(doMemoryCheck) {delete cms; return 0;}

    if(cms->GetPhysicsModel()==typeCvCfHiggs and debug) {
        for(int i=0; i< cms->Get_vv_productionMode().size(); i++){
            for(int j=0; j< cms->Get_vv_productionMode()[i].size(); j++){
                cout<<" channel "<<cms->GetChannelName(i)<<" proc "<<cms->GetProcessNames(i)[j]<<" pm = "<<cms->Get_vv_productionMode()[i][j]<<" dm= "<<cms->Get_vv_channelDecayMode()[i][j]<<endl;
            }
        }
        for(int i=0; i< cms->Get_vv_pdfs_productionMode().size(); i++){
            for(int j=0; j< cms->Get_vv_pdfs_productionMode()[i].size(); j++){
                cout<<" channel "<<cms->Get_v_pdfs_channelname()[i]<<" proc "<<cms->Get_vv_pdfs_procname()[i][j]<<" pm = "<<cms->Get_vv_pdfs_productionMode()[i][j]<<" dm= "<<cms->Get_vv_pdfs_channelDecayMode()[i][j]<<endl;
            }
        }
    }

    if(calcsignificance==0){// calc limits
        if(method == "Bayesian"){
            runBayesian();
        }else if(method == "Hybrid"){
            if(ExpectationHints == "Asymptotic"){
                TString tmp_dataset = dataset;
                dataset= "data_obs";
                vdata_global=cms_global->Get_v_data();
                vdata_global_th=cms_global->Get_v_dataTH();
                cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
                runAsymptoticLimits();
                if(scanRAroundBand=="all"){
                    dataset= "asimov_b";
                    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
                    vdata_global = cms->Get_AsimovData(0);
                    runAsymptoticLimits();
                }
                dataset = tmp_dataset;
                SaveResults(jobname+"Asymptotic"+"_limitbands", HiggsMass, preHints_obs, 0, 0, 0, preHints_m2s, preHints_m1s, preHints_median, 0, preHints_p1s, preHints_p2s);
            }else if(ExpectationHints == "Bayesian"){
                if(!doExpectation) {
                    if(scanRAroundBand=="all"){
                        doExpectation = true;
                        toysBayesian = 100;
                    }
                }
                runBayesian();
                cms->SetData(cms->Get_v_data_real());
                cms->SetDataForUnbinned(cms->Get_v_pdfs_roodataset_real());
            }else if(ExpectationHints!=""){
                TFile *fhints = new TFile(ExpectationHints);
                if(fhints->IsZombie()) {cout<<" File supposing to contain bands hints: "<<ExpectationHints<<" is not exist or it's bad, exit"<<endl; exit(1);};

                TTree *fChain = (TTree*) fhints->Get("T");
                if(fChain == NULL) {cout<<" TTree \"T\" is not exist in the file "<<ExpectationHints<<", exit"<<endl; exit(1);};
                // Declaration of leaf types
                Double_t        mH;
                Double_t        limit;
                Double_t        rm2s;
                Double_t        rm1s;
                Double_t        rmedian;
                Double_t        rp1s;
                Double_t        rp2s;

                // List of branches
                TBranch        *b_mH;   //!
                TBranch        *b_limit;   //!
                TBranch        *b_rm2s;   //!
                TBranch        *b_rm1s;   //!
                TBranch        *b_rmedian;   //!
                TBranch        *b_rp1s;   //!
                TBranch        *b_rp2s;   //!

                if(fChain->GetBranch("mH")){
                    fChain->SetBranchAddress("mH", &mH, &b_mH);
                }
                if(fChain->GetBranch("mh")){
                    fChain->SetBranchAddress("mh", &mH, &b_mH);
                }

                fChain->SetBranchAddress("limit", &limit, &b_limit);
                fChain->SetBranchAddress("rm2s", &rm2s, &b_rm2s);
                fChain->SetBranchAddress("rm1s", &rm1s, &b_rm1s);
                fChain->SetBranchAddress("rmedian", &rmedian, &b_rmedian);
                fChain->SetBranchAddress("rp1s", &rp1s, &b_rp1s);
                fChain->SetBranchAddress("rp2s", &rp2s, &b_rp2s);

                Long64_t nentries = fChain->GetEntries();
                bool only_one_entry = true;
                if(nentries != 1) only_one_entry = false;

                //Long64_t nbytes = 0, nb = 0;
                for (Long64_t jentry=0; jentry<nentries;jentry++) {
                    Long64_t ientry = fChain->GetEntry(jentry);
                    if (ientry < 0) break;
                    if((mH==HiggsMass && only_one_entry==false) or only_one_entry) {
                        preHints_median = rmedian;	
                        preHints_p2s = rp2s;
                        preHints_p1s = rp1s;
                        preHints_m1s = rm1s;
                        preHints_m2s = rm2s;
                        preHints_obs = limit;
                    }
                }
            }
            if(ExpectationHints!=""){
                cout<<" hints:  "<<endl;
                cout<<" observed limit = "<<preHints_obs<<endl;
                cout<<" bands :  "<<preHints_m2s<<" "<<preHints_m1s<<" "<<preHints_median<<" "<<preHints_p1s<<" "<<preHints_p2s<<endl;
                if(scanRAroundBand=="all"){
                    double tmpstep = (1.1*preHints_p2s - preHints_m2s*0.8)/20.;
                    if(tmpstep > 0.001) for(double tmpr = preHints_m2s*0.8; tmpr<=preHints_p2s*1.1; tmpr+=tmpstep){
                        vR_toEval.push_back(tmpr);
                    }

                    tmpstep = (1.2*preHints_m2s - preHints_m2s*0.8)/5.;
                    if(tmpstep > 0.001)for(double tmpr = preHints_m2s*0.75; tmpr<=preHints_m2s*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);

                    if(preHints_obs < preHints_m2s) {
                        tmpstep = (preHints_obs*1.2 - preHints_obs*0.8)/5.;
                        if(tmpstep>0.01)for(double tmpr = preHints_obs*0.8; tmpr<=preHints_obs*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                    }
                    if(preHints_obs > preHints_p2s) {
                        tmpstep = (preHints_obs*1.2 - preHints_obs*0.8)/5.;
                        if(tmpstep>0.01)for(double tmpr = preHints_obs*0.8; tmpr<=preHints_obs*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                    }
                }else if(scanRAroundBand=="obs"){
                    double tmpstep = (preHints_obs*1.2 - preHints_obs*0.8)/5.;
                    if(tmpstep>0.001)for(double tmpr = preHints_obs*0.8; tmpr<=preHints_obs*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                }else if(scanRAroundBand=="-2"){
                    double tmpstep = (preHints_m2s*1.2 - preHints_m2s*0.8)/5.;
                    if(tmpstep>0.001)for(double tmpr = preHints_m2s*0.8; tmpr<=preHints_m2s*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                }else if(scanRAroundBand=="-1"){
                    double tmpstep = (preHints_m1s*1.2 - preHints_m1s*0.8)/5.;
                    if(tmpstep>0.001)for(double tmpr = preHints_m1s*0.8; tmpr<=preHints_m1s*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                }else if(scanRAroundBand=="0"){
                    double tmpstep = (preHints_median*1.2 - preHints_median*0.8)/5.;
                    if(tmpstep>0.001)for(double tmpr = preHints_median*0.8; tmpr<=preHints_median*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                }else if(scanRAroundBand=="1"){
                    double tmpstep = (preHints_p1s*1.2 - preHints_p1s*0.8)/5.;
                    if(tmpstep>0.001)for(double tmpr = preHints_p1s*0.8; tmpr<=preHints_p1s*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                }else if(scanRAroundBand=="2"){
                    double tmpstep = (preHints_p2s*1.2 - preHints_p2s*0.8)/5.;
                    if(tmpstep>0.001)for(double tmpr = preHints_p2s*0.8; tmpr<=preHints_p2s*1.2; tmpr+=tmpstep) vR_toEval.push_back(tmpr);
                }
                std::sort(vR_toEval.begin(), vR_toEval.end());
                vector<double>::iterator it;
                it = std::unique(vR_toEval.begin(), vR_toEval.end());
                vR_toEval.resize(it - vR_toEval.begin());
            }


            cms->SetTossToyConvention(tossToyConvention);
            cms->SetUseBestEstimateToCalcQ(UseBestEstimateToCalcQ);

            // initialize the calculator
            CLsBase frequentist;
            frequentist.SetDebug(debug);
            frequentist.SetRdm(rdm);
            frequentist.SetTestStatistics(testStat);
            cms_global= cms;
            vdata_global=cms->Get_v_data();
            vdata_global_th=cms_global->Get_v_dataTH();

            CLsLimit clsr95;
            clsr95.SetDebug(debug);
            clsr95.SetRule(rule);
            double rtmp;
            clsr95.SetAlpha(1-CL);

            if(singlePoint){ cms->SetSignalScaleFactor(testR); }
            else { cms->SetSignalScaleFactor(1.); testR=1.;}

            frequentist.SetModel(cms);

            if(bOnlyEvaluateQdata){
                if(testStat==1)frequentist.prepareLogNoverB();
                frequentist.BuildM2lnQ_data();
                watch.Print();
                delete cms;
                return 0;
            }

            if(GenerateToysAtBestFitSB){
                //DoAfit(vR_toEval[i], cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
                double tmpr = 1, tmperr=1, ErrorDef=0;
                int success[1]={0};
                cms->SetNoErrorEstimation(1);
                double *pars = new double[cms->Get_max_uncorrelation()+1];
                MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //202, 201, 2:  allow mu to be negative
                cms->SetSignalScaleFactor(pars[0]);
                cms->Set_fittedParsInData_sb(pars);
                cout<<" ***  gonna generateToys with s+b best fit :  mu="<<pars[0]<<endl;
                //frequentist.BuildM2lnQ_sb(toysHybrid); // read from "-tH"
                frequentist.BuildM2lnQ_sb(toysHybrid, false, bWriteToys, jobname);
                watch.Print();
                delete cms;
                return 0;
            }

            //frequentist.checkFittedParsInData(true, false, "fittedPars.root");
            if(tossToyConvention==1 && !bFixNuisancsAtNominal )frequentist.checkFittedParsInData(bReadPars, bWritePars, fileFittedPars);
            if(tossToyConvention==1 && bFixNuisancsAtNominal){
                cms->Set_fittedParsInData_sb(cms->Get_norminalPars());
                cms->Set_fittedParsInData_b(cms->Get_norminalPars());
            }

            if(debug>=1000){
                cout<<" *** DELETEME checking 1"<<endl;
                int npars = cms_global->Get_max_uncorrelation();
                for(int i=0; i<=npars; i++) 
                    printf("DELETEME  par %30s      %.6f \n", i>0?cms->Get_v_uncname()[i-1].c_str():"signal_strength", cms->Get_fittedParsInData_b()[i]);
            }

            vector<double> vb, vsb;
            if(scanRs && vR_toEval.size()>0){
                for(int i=0; i<vR_toEval.size(); i++){
                    cout<<" running with signal stength = "<<vR_toEval[i]<<endl;
                    cms->SetSignalScaleFactor(vR_toEval[i]);
                    if(tossToyConvention==1 && bFixNuisancsAtNominal){
                        cms->Set_fittedParsInData_sb(cms->Get_norminalPars());
                        cms->Set_fittedParsInData_b(cms->Get_norminalPars());
                    }
                    if(tossToyConvention==1 && !bFixNuisancsAtNominal){
                        DoAfit(vR_toEval[i], cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
                        cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
                    }

                    if(testStat==1)frequentist.prepareLogNoverB();
                    if(nToysForCLsb<=0) nToysForCLsb=toysHybrid;
                    if(nToysForCLb<=0) nToysForCLb=toysHybrid;

                    int nTimes = 1; 
                    if(vR_toEval[i]<=1.2*preHints_m2s) nTimes = 2;
                    frequentist.BuildM2lnQ_sb(nToysForCLsb*nTimes);
                    vsb = frequentist.Get_m2logQ_sb();
                    frequentist.BuildM2lnQ_b(nToysForCLb*nTimes);
                    vb = frequentist.Get_m2logQ_b();

                    frequentist.BuildM2lnQ_data();
                    double qdata = frequentist.Get_m2lnQ_data();
                    cout<<" Q[mu="<<vR_toEval[i]<<"] ="<<qdata <<endl;
                    TString s_r; s_r.Form("TESTED_R%.5f", cms->GetSignalScaleFactor());
                    TString s_qdata; s_qdata.Form("DATA_R%.5f_Q%.5f", cms->GetSignalScaleFactor(), qdata);
                    TString s_sb = "SAMPLING_SB_"; s_sb+=s_r;
                    TString s_b = "SAMPLING_B_"; s_b+=s_r;
                    TString fileM2lnQ = jobname; fileM2lnQ+="_m2lnQ_deprecated.root";
                    TString option = (i==0?"RECREATE":"UPDATE");
                    //FillTree(fileM2lnQ, cms->GetSignalScaleFactor(), qdata, vsb, vb, s_r, s_qdata, s_sb, s_b, option);
                    fileM2lnQ = jobname; fileM2lnQ+="_m2lnQ.root";
                    FillTree2(fileM2lnQ, cms->GetSignalScaleFactor(), qdata, vsb, vb, s_r, s_sb, s_b, option);
                }
                watch.Print();
                delete cms;
                return 0;
            }

            //frequentist.BuildM2lnQ(toysHybrid);
            if(!bSkipM2lnQ){
                if(testStat==1)frequentist.prepareLogNoverB();
                if(!bNotCalcQdata){
                    frequentist.BuildM2lnQ_data();
                    cout<<" Q_data = " << frequentist.Get_m2lnQ_data()<<endl;
                }
                if(nToysForCLsb<=0) nToysForCLsb=toysHybrid;
                if(nToysForCLb<=0) nToysForCLb=toysHybrid;
                frequentist.BuildM2lnQ_sb(nToysForCLsb, false, bWriteToys, jobname);
                vsb = frequentist.Get_m2logQ_sb();

                if(debug>=1000){
                    cout<<" *** DELETEME checking 2"<<endl;
                    int npars = cms_global->Get_max_uncorrelation();
                    for(int ii=0; ii<=npars; ii++) 
                        printf("DELETEME  par %30s      %.6f \n", ii>0?cms->Get_v_uncname()[ii-1].c_str():"signal_strength", cms->Get_fittedParsInData_b()[ii]);
                }
                frequentist.BuildM2lnQ_b(nToysForCLb, false, bWriteToys, jobname);

                vb = frequentist.Get_m2logQ_b();
                if(bWriteToys){
                    watch.Print();
                    delete cms;
                    return 0;
                }
            }

            if(makeCLsBands>=2){
                if(debug) cout<<"MakeCLsValues from precomputed file = "<<sFileM2lnQGrid<<endl;
                std::map<double, TTree*> gridCLsb; //r, <clsb, clserr>
                std::map<double, TTree*> gridCLb; //r, <clsb, clserr>
                std::map<double, double> gridQdata; //r, q_data
                ReadM2lnQGridFromFile(sFileM2lnQGrid, gridCLsb, gridCLb, debug);
                ReadM2lnQGridFromFile(sFileM2lnQGrid, gridQdata, debug);

                TTree *tsb = gridCLsb[testR];
                TTree *tb = gridCLb[testR];
                double qdata = gridQdata[testR];

                if(tsb==0 || tb==0) {
                    cout<<"ERROR: grid of R="<<testR<<" not in the file"<<endl; 
                    cout<<"   all available R's : "<<endl;
                    for (std::map<double, TTree *>::iterator itg = gridCLsb.begin(), edg = gridCLsb.end(); itg != edg; ++itg) {
                        if(itg->first == testR) continue;
                        cout<<"  "<<itg->first<<"  ";
                    }
                    cout<<endl;
                    exit(1);
                }

                GetM2lnQ(tsb, tb, vsb, vb, debug);
                clsr95.DoingStatisticalBandsForCLs(vsb, vb);	
                double cls, errs;
                GetCLs(qdata, tsb, tb, cls, errs);
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;

                double * clsbands = clsr95.Get_CLsProjected();
                SaveResults(jobname+"_clsbands", HiggsMass, cls, errs, 0, 0, clsbands[0], clsbands[1], clsbands[2], clsbands[5], clsbands[3], clsbands[4]);
                delete cms;
                return 0;
            }

            if(!bSkipM2lnQ && makeCLsBands==1){
                clsr95.DoingStatisticalBandsForCLs(vsb, vb);	
                double errs;
                double cls = frequentist.CLs(errs);
                double * clsbands = clsr95.Get_CLsProjected();
                SaveResults(jobname+"_clsbands", HiggsMass, cls, errs, 0, 0, clsbands[0], clsbands[1], clsbands[2], clsbands[5], clsbands[3], clsbands[4]);
            }

            if(!bNotCalcCLssbb && !bSkipM2lnQ){
                double errs, errb, errsb;
                double cls = frequentist.CLs(errs);
                double clsb = frequentist.CLsb(errsb);
                double clb = frequentist.CLb(errb);
                cout<<"---------------testing at signal strength r = "<<cms->GetSignalScaleFactor()<<"-------------------------"<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<endl;
                cout<<"------------------------------------------------------------"<<endl;
            }

            if(bSaveM2lnQ && !bSkipM2lnQ){
                double qdata = bNotCalcQdata?0:frequentist.Get_m2lnQ_data();
                TString s_r; s_r.Form("TESTED_R%.5f", cms->GetSignalScaleFactor());
                TString s_qdata; s_qdata.Form("DATA_R%.5f_Q%.10f", cms->GetSignalScaleFactor(), qdata);
                TString s_sb = "SAMPLING_SB_"; s_sb+=s_r;
                TString s_b = "SAMPLING_B_"; s_b+=s_r;
                TString fileM2lnQ = jobname; fileM2lnQ+="_m2lnQ.root";
                FillTree(fileM2lnQ, cms->GetSignalScaleFactor(), qdata, vsb, vb, s_r, s_qdata, s_sb, s_b);
                fileM2lnQ = jobname; fileM2lnQ+="_m2lnQ2.root";
                FillTree2(fileM2lnQ, cms->GetSignalScaleFactor(), qdata, vsb, vb, s_r, s_sb, s_b);
            }

            if(0){// throw pseudo data
                double *pars = NULL;
                if(tossPseudoDataConvention==1){
                    pars = cms->Get_fittedParsInData_b();
                    if(pars==NULL) DoAfit(0, cms->Get_v_data(), cms->Get_v_pdfs_roodataset(), pars);
                }
                VIChannel vpseudodata; 
                vpseudodata=cms->GetToyData_H0(pars);
                vector<RooDataSet*> vrds; vrds.clear();
                for(int c=0; c<cms->Get_vv_pdfs().size(); c++){
                    RooDataSet *rds = new RooDataSet(*((RooDataSet*)cms->Get_v_pdfs_roodataset_toy()[c]));
                    vrds.push_back(rds);
                }
                TFile *f = new TFile("data.root", "RECREATE");
                TH1D * hobs = new TH1D("hobs","hobs", vpseudodata.size(), 0, vpseudodata.size());
                for(int i=0; i<vpseudodata.size(); i++) hobs->SetBinContent(i+1, vpseudodata[i]);
                f->WriteTObject(hobs);
                for(int c=0; c<cms->Get_vv_pdfs().size(); c++){
                    vrds[c]->SetName(cms->Get_v_pdfs_channelname()[c].c_str());
                    f->WriteTObject(vrds[c]);
                }
                f->Close();

            }

            //delete cms;
            for(int hello=0; hello<-1; hello++){
                RooRandom::randomGenerator()->SetSeed(seed);
                CRandom *rdm1 = new CRandom(seed);  //initilize a random generator
                CountingModel *cmsClone; // this instance will be the one combining all sub cards
                cmsClone = new CountingModel();
                CountingModel *tmpn1[100];
                CountingModel *tmpn[100];// = new CountingModel(); 
                if(datacards.size()>100) {cout<<"too many datacards "<<datacards.size()<<endl; exit(0);}
                for(int i=0; i<datacards.size(); i++){
                    tmpn[i] = new CountingModel();
                    ConfigureModel(tmpn[i], HiggsMass, datacards[i].Data());
                    tmpn[i]->SetUseSystematicErrors(true);
                }
                if(datacards.size()==1) cmsClone=tmpn[0];
                else if(datacards.size()>=2){
                    tmpn1[1] = CombineModels(tmpn[0], tmpn[1]);
                    tmpn1[1]->SetUseSystematicErrors(true);
                    for(int i=2; i<datacards.size(); i++){
                        tmpn1[i] = CombineModels(tmpn1[i-1], tmpn[i]);
                        tmpn1[i]->SetUseSystematicErrors(true);
                    }	
                    cmsClone = tmpn1[datacards.size()-1];
                }else{exit(0);}
                cmsClone->SetUseSystematicErrors(systematics);
                // done combination

                cmsClone->SetRdm(rdm1);
                cmsClone->SetUseSystematicErrors(systematics);
                int nch_removed = cmsClone->RemoveChannelsWithExpectedSignal0orBkg0(0); // 0: remove only bins with total bkg<=0,  1: remove bins with total sig<=0,  2: both
                cmsClone->SetMoveUpShapeUncertainties(1);
                cmsClone->SetTossToyConvention(tossToyConvention);

                cmsClone->SetDebug(debug);
                cmsClone->SetUseBestEstimateToCalcQ(UseBestEstimateToCalcQ);
                frequentist.SetModel(cmsClone);

                frequentist.checkFittedParsInData(true, false, "fittedPars.root");

                if(singlePoint){ cmsClone->SetSignalScaleFactor(testR); }
                else { cmsClone->SetSignalScaleFactor(1.); }
                frequentist.BuildM2lnQ(toysHybrid);
                double errs, errb, errsb;
                double cls = frequentist.CLs(errs);
                double clsb = frequentist.CLsb(errsb);
                double clb = frequentist.CLb(errb);
                cout<<"---------------testing at signal strength r = "<<cmsClone->GetSignalScaleFactor()<<"-------------------------"<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<endl;
                cout<<"------------------------------------------------------------"<<endl;


                cmsClone->~CountingModel();
                for(int i=0; i<datacards.size()-1; i++) {
                    if(tmpn[i]) tmpn[i]->~CountingModel();
                    if(tmpn1[i]) tmpn1[i]->~CountingModel();
                }
            }

            if(bPlots && !bSkipM2lnQ){
                DrawPdfM2logQ pdfM2logQ(frequentist.Get_m2logQ_sb(),frequentist.Get_m2logQ_b(), frequentist.Get_m2lnQ_data(), 
                        "-2lnQ on data", "; -2lnQ; entries", "lnq", pt);
                pdfM2logQ.draw();
                double m2lnqdata = frequentist.Get_m2lnQ_data();
                //cout<<m2lnqdata<<endl;
                vector<double> vsb = frequentist.Get_m2logQ_sb();
                for(int i=0; i<vsb.size(); i++){
                    //	if(vsb[i]>=m2lnqdata)cout<<i<<" "<<vsb[i]<<endl;
                }
                FillTree("m2lnQ_b.root", frequentist.Get_m2logQ_b());
                FillTree("m2lnQ_sb.root", frequentist.Get_m2logQ_sb());
            }

            if(debug>=10 && !bSkipM2lnQ) {
                cout<<"-2lnQ on data = "<<frequentist.Get_m2lnQ_data()<<endl;
                vector<double> qsb = frequentist.Get_m2logQ_sb();
                vector<double> qb = frequentist.Get_m2logQ_b();
                int n10 = int(qsb.size()/10); if(n10<=0) n10=1;
                if(debug>=100) n10=1;
                cout<<"-2lnQ for SB"<<endl;
                for(int i=0; i<qsb.size(); i+=n10) cout<<qsb[i]<<endl;
                n10 = int(qb.size()/10); if(n10<=0) n10=1;
                if(debug>=100) n10=1;
                cout<<"-2lnQ for B"<<endl;
                for(int i=0; i<qb.size(); i+=n10) cout<<qb[i]<<endl;
            }

            if(singlePoint){ 
                watch.Print();
                delete cms;
                return 1;
            }

            if(scanRs){
                vector<double> vr, vc; vr.clear(); vc.clear();
                if(nSteps<=0) { cout<<"ERROR: nSteps must be > 0"<<endl; return 0; }
                for(int i=0; i<nSteps+1; i++){
                    double testr = initialRmin + i*(initialRmax - initialRmin)/(float)nSteps;
                    cms->SetSignalScaleFactor(testr);
                    frequentist.BuildM2lnQ(cms, toysHybrid);
                    double errs;
                    double cls = frequentist.CLs(errs);
                    vr.push_back(testr);
                    vc.push_back(cls);
                }
                printf("\n results of scanned r vs. CLs: \n");
                for(int i=0; i<vr.size(); i++){
                    printf("   r=%10.3f  CLs=%7.5f\n", vr[i], vc[i]);
                }
                if(bPlots){
                    DrawEvolution2D d2d(vr, vc, "; r ; CLs", (jobname+"_scanned_r_vs_cl").Data(), pt);
                    d2d.draw();
                    d2d.save();
                }

                watch.Print();
                delete cms;
                return 0;
            }


            if(bCalcObservedLimit){
                if(bQuickEstimateInitialLimit) rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist, toysHybrid);
                else rtmp = clsr95.LimitOnSignalScaleFactor(cms, initialRmin, initialRmax, &frequentist, toysHybrid, 3);

                if(bPlots){
                    DrawEvolution2D d2d(clsr95.GetvTestedScaleFactors(), clsr95.GetvTestedCLs(), "; r ; CLs", (jobname+"_r_vs_cl").Data(), pt);
                    d2d.draw();
                    d2d.save();

                }
                cout<<"------------------------------------------------------------"<<endl;
                if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
                cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<endl;
                cout<<"------------------------------------------------------------"<<endl;

                SaveResults(jobname+"_freqObsLimit", HiggsMass, rtmp, clsr95.LimitErr(), 0, 0, 0, 0, 0, 0, 0, 0);
            }


            if(doExpectation && sFileLimitsDistribution==""){
                cms->SetSignalScaleFactor(1.);
                LimitBands lb(&clsr95, &frequentist, cms);	
                lb.SetQuickEstimateInitialLimit(bQuickEstimateInitialLimit);
                lb.SetInitialR(initialRmin, initialRmax);
                lb.SetTossPseudoDataConvention(tossPseudoDataConvention);
                lb.SetDebug(debug);
                lb.SetPlotLevel(bPlots);
                lb.IsM2lnQGridPreComputed(bM2lnQGridPreComputed, sFileM2lnQGrid);
                int noutcomes = toys;
                if(bOnlyEvalCL_forVR){
                    lb.Set_bOnlyEvalCL_forVR(true);	
                    lb.Set_vR_toEval(vR_toEval);
                    lb.CLsLimitBands(1-CL, noutcomes, toysHybrid);
                    vector< vector<double> > vvCL_forVR = lb.Get_vvCL_forVR();
                    cout<<endl<<"VR: ";
                    for(int jj=0; jj<vR_toEval.size(); jj++){
                        printf(" %8.5f ",vR_toEval[jj]);
                    }
                    cout<<endl;
                    for(int ii=0; ii<vvCL_forVR.size(); ii++){
                        cout<<"CL: ";
                        for(int jj=0; jj<vvCL_forVR[ii].size(); jj++){
                            printf(" %8.5f ",vvCL_forVR[ii][jj]);
                        }
                        cout<<endl;
                    }
                    watch.Print();
                    delete cms;
                    return 0;
                }
                lb.CLsLimitBands(1-CL, noutcomes, toysHybrid);
                rmean=lb.GetCLsLimitMean();
                difflimits = lb.GetDifferentialLimitsCLs();
                saveExpectation();
            }

        }else if(method == "FeldmanCousins"){

            CLsBase frequentist;
            frequentist.SetDebug(debug);
            frequentist.SetRdm(rdm);
            cms_global= cms;
            //vdata_global=cms->Get_v_data();

            double tmp;
            vdata_global=cms->Get_v_data();

            // FeldmanCousins must use specified test statisics:  profile likelihood ratio, which allow mu hat > probed mu
            frequentist.SetTestStatistics(31);

            cms->SetAllowNegativeSignalStrength(false);


            CLsLimit clsr95;
            clsr95.SetDebug(debug);
            double rtmp;
            clsr95.SetAlpha(1-CL);


            clsr95.SetAdaptiveSampling(adaptiveSampling);
            bool lowerLimit_ = false;

            double r95_fc;

            double fcMid = 0, fcErr = 0; 
            double rmin = initialRmin,  rmax = initialRmax;

            if(bQuickEstimateInitialLimit){
                BayesianBase bys(cms, 0.05, 1.e-2);
                bys.SetNumToys(1000);
                rtmp= bys.Limit();
                rmin = rtmp/10.,  rmax = rtmp*2;
            }

            int nsteps = 10;
            vector< vector<double> > vv_all;
            vector< vector<double> > vv;
            // the following iterating algothm comes from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HiggsAnalysis/CombinedLimit/src/FeldmanCousins.cc
            do { 
                if (debug) std::cout << "scan in range [" << rmin << ", " << rmax << "]" << std::endl;
                if(fcMid!=0){
                    double bin = (rmax-rmin)/10.;
                    rmin = rmin+0.5*bin;
                    rmax = rmax-0.5*bin;
                }
                r95_fc = clsr95.FeldmanCousins(cms, rmin, rmax, &frequentist, toysHybrid, nsteps);
                vv = clsr95.GetFCconstruction();
                for(int i=0; i<vv.size(); i++){
                    vv_all.push_back(vv[i]);
                }
                int found = -1; 
                if (lowerLimit) {
                    for(int i=0; i<nsteps; ++i){
                        int n = vv[i].size();
                        if(vv[i][n-1]>vv[i][n-2]) found=i;
                        else break;
                    }
                } else {
                    for(int i=0; i<nsteps; ++i){
                        int n = vv[i].size();
                        bool inside = (vv[i][n-1]<=vv[i][n-2]) ; // important ...
                        if (inside) found = i;
                        else if (found != -1) break;
                    }
                    if (found == -1) {
                        std::cout << "Points are either all inside or all outside the bound." << std::endl;
                        delete cms;
                        return false;
                    }
                }
                double fcBefore = (found > -1 ? (rmin+(rmax-rmin)/double(nsteps)*found) : rmin);
                double fcAfter  = (found < nsteps-1 ? (rmin+(rmax-rmin)/double(nsteps)*(found+1)) : rmax);
                fcMid = 0.5*(fcAfter+fcBefore);
                fcErr = 0.5*(fcAfter-fcBefore);
                if (debug) std::cout << "  would be r < " << fcMid << " +/- "<<fcErr << std::endl;
                rmin=(std::max(rmin, fcMid-3*fcErr)); 
                rmax=(std::min(rmax, fcMid+3*fcErr));
                if (fcErr < 4*std::max(rAbsAcc, rRelAcc* fcMid)) { // make last scan more precise
                    clsr95.SetAdditionalNToysFactor(4*toysFactor);
                }
            } while (fcErr > std::max(rAbsAcc, rRelAcc* fcMid));

            r95_fc = fcMid;
            if (0) {
                std::cout << "\n -- FeldmanCousins++ -- \n";
                std::cout << "Limit: r " << (lowerLimit_ ? "> " : "< ") << fcMid << " +/- "<<fcErr << "\n";
            }


            if(debug) cout <<"95\% CL upper limit by FC: "<<r95_fc<<",   use bys: "<<rtmp<<endl;

            cout<<"------------------------------------------------------------"<<endl;
            if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
            cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<r95_fc<<endl;
            cout<<"------------------------------------------------------------"<<endl;

            SaveResults(jobname+"_fcObsLimit", HiggsMass, r95_fc, 0, 0, 0, 0, 0, 0, 0, 0, 0);

            if(bPlots){
                // for plots...
                TCanvas *can = new TCanvas("c","c");	
                TString ssave = "fc";
                TString stmp = "hframe"; stmp += ssave;
                double xmax = 0;
                for(int i=0; i<vv_all.size(); i++){
                    int n = vv_all[i].size();
                    if(vv_all[i][n-2]>xmax) xmax= vv_all[i][n-2];
                }
                TH1F *hframe= new TH1F(stmp, "; q_{#mu}; #mu", 1000, 0, 2*xmax);
                hframe->SetMinimum(rtmp/10.);
                hframe->SetMaximum(rtmp*2.);
                hframe->SetStats(0);
                hframe->SetFillStyle(1);
                hframe->Draw(" ");

                float *r_tested = new float[vv_all.size()];
                float *q_data = new float[vv_all.size()];

                float q95=0;
                for(int i=0; i<vv_all.size(); i++){
                    int n = vv_all[i].size();
                    double x[2]={0, vv_all[i][n-2]};	// q_up	
                    double y[2]={vv_all[i][n-3], vv_all[i][n-3]}; // r
                    TGraph *gr = new TGraph(2, x, y);
                    gr->Draw("l same");

                    r_tested[i]=vv_all[i][n-3];
                    q_data[i]=vv_all[i][n-1];

                    //if(r95_fc==vv_all[i][n-3]) q95 = q_data[i];
                    if(fabs(r95_fc-r_tested[i])/r95_fc<0.00001) q95 = q_data[i];
                }

                TGraph gr(vv_all.size(), q_data, r_tested);
                gr.Sort(&TGraph::CompareY);
                gr.Draw("l same");

                TArrow *arrow95 = new TArrow(0, r95_fc, q95, r95_fc, 0.03, "|>");
                arrow95->SetLineWidth(3);
                arrow95->SetLineColor(kRed);
                //arrow95->SetFillColor(kRed);
                arrow95->SetFillStyle(0);
                arrow95->Draw();
                Save(can, "fc");
                delete r_tested;
                delete q_data;
            }

            watch.Print();
            delete cms;
            return 1;
        }else if(method=="ProfiledLikelihood" or method=="ProfileLikelihood"){
            vdata_global=cms->Get_v_data();
            vdata_global_th=cms_global->Get_v_dataTH();
            runProfileLikelihoodApproximation();
            watch.Print();
            delete cms;
            return 0;

        }else if(method=="AsymptoticCLs"){
            runAsymptoticCLs();
            watch.Print();
            delete cms;
            return 0;
        }else if(method=="Asymptotic"){
            runAsymptoticLimits();
            watch.Print();
            delete cms;
            return 0;
        }else if(method == "MultiDimFit"){

            cms_global= cms;
            int idMH=-1;
            vector<string> v_uncname = cms_global->Get_v_uncname();
            for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]=="MH") idMH=i;
            for(int k=0; k<cms_global->POIs().size(); k++){
                for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]== cms_global->POIs()[k].name ) 
                { cout<<" ***** psPOI:  "<<v_uncname[i-1]<<endl; viPseudoPOIs.push_back(i); vsPseudoPOIs.push_back(v_uncname[i-1]); }
            }

            // data fit
            vdata_global=cms->Get_v_data();
            vdata_global_th=cms->Get_v_dataTH();
            if(cms->hasParametricShape()){
                cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
            }
            runMaxLikelihoodFit(true, true);


            // toy fit  s+b 
            if(doExpectation){
                //FIXME SetSignalStrength(ToysAtDifferentSignalStrength, InjectingSignalStrength);
                vector<double> vmuhat, vmuerrhi, vmuerrlo, vmhhat, vmherrhi, vmherrlo;
                vector< vector<double> > vvdPseudoPOIs; vvdPseudoPOIs.resize(viPseudoPOIs.size()); 
                TFile * fin;
                if(sToysFromFile!=""){
                    toys = cms->NumberOfToysFromFile(sToysFromFile);
                    fin = new TFile(sToysFromFile);
                }
                for(int t=0; t<toys; t++)
                {
                    if(debug>=0)cout<<" Toy "<<t<<endl;
                    if(sToysFromFile!=""){
                        vdata_global = (VDChannel)cms->GetToyDataFromFile(fin, t);
                        if(cms->hasParametricShape()){
                            cms->SetTmpDataForUnbinned(cms->GetToyUnbinnedDataFromFile(fin, t));
                        }
                    }
                    else{
                        for(int p=0; p<vsParToBeFixedWhenGenToys.size(); p++){
                            int poi_i =-1;
                            for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]==vsParToBeFixedWhenGenToys[p].Data()) poi_i=i;
                            if( poi_i >=0 ){
                                if(cms->Get_v_pdftype()[poi_i] == typeFlat){
                                    vector< vector<double> > v_Pars = cms->Get_v_Pars(); 
                                    double tmp = ( vdParToBeFixedWhenGenToys[p] - v_Pars[poi_i][1] )/ ( v_Pars[poi_i][2] - v_Pars[poi_i][1]);
                                    cms->SetParForGenToy(vsParToBeFixedWhenGenToys[p], tmp); // set to init value from command line
                                }else	cms->SetParForGenToy(vsParToBeFixedWhenGenToys[p], vdParToBeFixedWhenGenToys[p]); // set to init value from command line
                            }
                            else {cout<<"ERROR:  fixing a parameter which is not exist in the list : "<<vsParToBeFixedWhenGenToys[p]<<endl; exit(1);}
                        }
                        vdata_global = (VDChannel)cms->GetToyData_H1();//_cms->Get_fittedParsInData_sb());
                        if(cms->hasParametricShape()){
                            cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_toy());
                        }
                        cms->ClearParForGenToy();

                    }
                    runMaxLikelihoodFit(debug>0?true:false,debug>0?true:false, TString(t));
                    vmuhat.push_back(fitBestMu);
                    vmuerrhi.push_back(fitBestMuErrHi);
                    vmuerrlo.push_back(fitBestMuErrLo);
                    vmhhat.push_back(fitBestMH);
                    vmherrhi.push_back(fitBestMHErrHi);
                    vmherrlo.push_back(fitBestMHErrLo);

                    for(int p=0; p<viPseudoPOIs.size(); p++){
                        if(cms->Get_v_pdftype()[viPseudoPOIs[p]] == typeFlat){
                            int poi_i= viPseudoPOIs[p];
                            vector< vector<double> > v_Pars = cms->Get_v_Pars(); 
                            double tmp = vdPseudoPOIs[p] *(v_Pars[poi_i][2]-v_Pars[poi_i][1])+v_Pars[poi_i][1];
                            vvdPseudoPOIs[p].push_back( tmp );
                        }else vvdPseudoPOIs[p].push_back(vdPseudoPOIs[p]);
                    }
                }

                vector< vector<double> > vvd; vvd.push_back(vmuhat); vvd.push_back(vmuerrlo); vvd.push_back(vmuerrhi); 
                //vvd.push_back(vmhhat); 
                vvd.push_back(vmherrlo); vvd.push_back(vmherrhi); 
                vector< TString > vs; vs.push_back("fitBestMu"); vs.push_back("fitBestMuErrLo"); vs.push_back("fitBestMuErrHi"); 
                //vs.push_back("fitBestMH"); 
                vs.push_back("fitBestMHErrLo"); vs.push_back("fitBestMHErrHi"); 
                for(int p=0; p<viPseudoPOIs.size(); p++){
                    vvd.push_back(vvdPseudoPOIs[p]);
                    vs.push_back(vsPseudoPOIs[p]);
                }
                TString ts=jobname; ts+="_maxlls";
                FillTree(ts, vvd, vs); 
                //CopyTrees(jobname+"_maxllfit.root", "T", ts+"_tree.root", "T_data");
            }

            watch.Print();

        }else if(method == "ScanningMuFit" or method=="MaxLikelihoodFit"){
            if(debug) cout<<" YouAreHere MaxLikelihoodFit bFixNuisancesAtBestFit"<<bFixNuisancesAtBestFit<<endl;
            int idMH=-1;
            vector<string> v_uncname = cms_global->Get_v_uncname();
            for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]=="MH") idMH=i;
            vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();

            cms_global= cms;
            vdata_global=cms->Get_v_data();
            vdata_global_th=cms_global->Get_v_dataTH();
            cms_global->SetDebug(debug);
            _inputNuisances = cms->Get_norminalPars();
            _startNuisances = cms->Get_norminalPars();

            if(bDumpFitResults)cms->DumpFitResults(_inputNuisances, jobname+"_nominalShape");

            double *pars = new double[cms->Get_max_uncorrelation()+1]; // nsys + r
            for(int i=0; i<=cms->Get_max_uncorrelation(); i++) pars[i]=0;

            //double bestFitPars[cms->Get_maximumFunctionCallsInAFit()+1];  //  not static ... not allowed 
            double *bestFitPars= new double[cms->Get_max_uncorrelation()+1];
            for(int i=0; i<=cms->Get_max_uncorrelation(); i++) bestFitPars[i]=0;

            double tmp=0, tmpe=0, tmpr = 0, tmperr=0;
            double ErrorDef = TMath::ChisquareQuantile(0.68 , NumDegreeOfFreedom);// (confidenceLevel, ndf)
            _ErrorDef = ErrorDef;
            int success[1]={0};
            _countPdfEvaluation = 0;
            if(bDumpFitResults){
                lands::_bDumpFinalFitResults = true;
            }
            double y0_2=0;
            // for CvCf   ==>   y0_2,  CV, CF and their errors 

            if(mainPOI != "") {
                cms->keepOnlyPOI(mainPOI);	
            }
            if(DoNotRunGlobalFit==false){
                if(NoErrorEstimate) cms->SetNoErrorEstimation(1);
                if(vCouplingsDef.size()==0) y0_2 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //202, 201, 2:  allow mu to be negative
                else { y0_2 =  MinuitFit(3, tmp, tmperr, 1.0, pars, false, debug, 0, 0); }
                if(debug) cout<<" _countPdfEvaluation = "<<_countPdfEvaluation<<endl;
                cout<<y0_2<<" fitter u="<<pars[0]<<", from minos fit asymmetry 68% CL:  ["<<tmperr<<","<<tmpr<<"]"<<endl; // upperL, lowerL
                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                    myMinuit->GetParameter(i, tmp, tmpe);
                    bestFitPars[i]=tmp;
                }
                cout<<" fAmin = "<<myMinuit->fAmin<<endl;
            }
            //double y0_2 =  MinuitFit(102, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //102, 101, 21:  don't allow mu to be negative
            double mu_hat = pars[0];
            double mu_hat_up = tmpr;
            double mu_hat_low = tmperr;
            watch.Print();


            if(cms->GetPhysicsModel()==typeCvCfHiggs and !DoNotRunGlobalFit) cms->ShowCvCfHiggsScales(bestFitPars);
            if(cms->GetPhysicsModel()==typeCvCfHiggs and 0) {
                // check 10000 calls of GammaTot  timing 
                TStopwatch watch2;
                watch2.Start();
                for(int i=0; i<500000; i++){
                    cms->SetFlatPars(bestFitPars);
                    cms->CalcGammaTot();
                }
                watch2.Stop();
                cout<<"500000 GammaTot calls take: "<<endl;
                watch2.Print();
            }
            if(cms->GetPhysicsModel()==typeC5Higgs ){
                if(ErrEstAlgo =="Minos"){
                    cout<<" 	68% CL : "<<endl;
                    for(int k=1; k<cms_global->POIs().size(); k++)cout<<" POI: "<<cms_global->POIs()[k].name<<"= "<<cms_global->POIs()[k].value<<" + "<<cms_global->POIs()[k].errUp<<" - "<<fabs(cms_global->POIs()[k].errDown)<<endl;
                }
                if(ErrEstAlgo=="Bisect"){
                    double glbminCgg = cms->Get_v_Pars()[cms->Get_par_i("Cgg")][0];
                    double glbminCbb = cms->Get_v_Pars()[cms->Get_par_i("Cbb")][0];
                    double glbminCtt = cms->Get_v_Pars()[cms->Get_par_i("Ctt")][0];
                    double glbminCvv = cms->Get_v_Pars()[cms->Get_par_i("Cvv")][0];
                    double glbminCglgl = cms->Get_v_Pars()[cms->Get_par_i("Cglgl")][0];
                    double errhiCgg = Err_Bisection(cms, "Cgg", bestFitPars, y0_2, glbminCgg, 1);
                    double errloCgg = Err_Bisection(cms, "Cgg", bestFitPars, y0_2, glbminCgg, 0);
                    double errhiCbb = Err_Bisection(cms, "Cbb", bestFitPars, y0_2, glbminCbb, 1);
                    double errloCbb = Err_Bisection(cms, "Cbb", bestFitPars, y0_2, glbminCbb, 0);
                    double errhiCtt = Err_Bisection(cms, "Ctt", bestFitPars, y0_2, glbminCtt, 1);
                    double errloCtt = Err_Bisection(cms, "Ctt", bestFitPars, y0_2, glbminCtt, 0);
                    double errhiCglgl = Err_Bisection(cms, "Cglgl", bestFitPars, y0_2, glbminCglgl, 1);
                    double errloCglgl = Err_Bisection(cms, "Cglgl", bestFitPars, y0_2, glbminCglgl, 0);
                    double errhiCvv = Err_Bisection(cms, "Cvv", bestFitPars, y0_2, glbminCvv, 1);
                    double errloCvv = Err_Bisection(cms, "Cvv", bestFitPars, y0_2, glbminCvv, 0);

                    cout<<" POI Cgg: "<<glbminCgg<<" +"<<errhiCgg<<" -"<<errloCgg<<endl;
                    cout<<" POI Cvv: "<<glbminCvv<<" +"<<errhiCvv<<" -"<<errloCvv<<endl;
                    cout<<" POI Cbb: "<<glbminCbb<<" +"<<errhiCbb<<" -"<<errloCbb<<endl;
                    cout<<" POI Ctt: "<<glbminCtt<<" +"<<errhiCtt<<" -"<<errloCtt<<endl;
                    cout<<" POI Cglgl: "<<glbminCglgl<<" +"<<errhiCglgl<<" -"<<errloCglgl<<endl;
                }
                watch.Print();
                delete pars;
                delete cms;
                return 0;

            }

            if(bRunMinuitContour && idMH>0 ){

                TFile f(jobname+"_contourRvsM.root","RECREATE");

                // using mncont to extract points around contour = 1 or 2 sigma   
                myMinuit->SetErrorDef(TMath::ChisquareQuantile(0.68, NumDegreeOfFreedom));
                TGraph * g1sigma = (TGraph*)myMinuit->Contour(100, 0, idMH);
                if(g1sigma){
                    g1sigma->SetName("ContourRvsM_1");
                    for(int i=0; i<g1sigma->GetN(); i++){
                        double r, m; g1sigma->GetPoint(i, r,m);
                        g1sigma-> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
                    }
                    f.WriteTObject(g1sigma); 
                } else {cout<<" WARNING: minuit doesn't give 1 sigma contour (return 0 points) "<<endl;}


                myMinuit->SetErrorDef(TMath::ChisquareQuantile(0.95, NumDegreeOfFreedom));
                TGraph * g2sigma = (TGraph*)myMinuit->Contour(100, 0, idMH);
                if(g2sigma){
                    g2sigma->SetName("ContourRvsM_2");
                    for(int i=0; i<g2sigma->GetN(); i++){
                        double r, m; g2sigma->GetPoint(i, r,m);
                        g2sigma-> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
                    }
                    f.WriteTObject(g2sigma);
                } else {cout<<" WARNING: minuit doesn't give 2 sigma contour (return 0 points) "<<endl;}
                f.Close();
            }


            {  // check covariance matrix and plot the ellipse if MH variable present in the floating parameter list 
                int npar = cms_global->Get_max_uncorrelation()+1;
                TMatrixDSym mat(npar); 
                if(0) {
                    myMinuit->mnemat( mat.GetMatrixArray(), npar);
                    mat.Print();
                }
                if(0){

                    TMatrixDSym *cormat = (TMatrixDSym*)mat.Clone();
                    for (Int_t i=0 ; i<cormat->GetNrows() ; i++) {
                        for (Int_t j=0 ; j<cormat->GetNcols() ; j++) {
                            if (i!=j) {
                                (*cormat)(i,j) = (*cormat)(i,j) / sqrt((*cormat)(i,i)*(*cormat)(j,j)) ;
                            }
                        }
                    }
                    for (Int_t i=0 ; i<cormat->GetNrows() ; i++) {
                        (*cormat)(i,i) = 1.0 ;
                    }

                    if(debug)cormat->Print();

                    if(vCouplingsDef.size() > 0){
                        //  the signal_strength is fixed at 1.
                        // it's put in the last par in the matrix
                        vector<int> map1, map2 ;

                        for(int i=0; i<cms->Get_max_uncorrelation();  i++){
                            map1.push_back(i);
                        } 
                        map2.push_back(cms->Get_max_uncorrelation()); // the last par is swtiched to be signal_strength 

                        TMatrixDSym S11, S22 ;
                        TMatrixD S12, S21 ;
                        RooMultiVarGaussian::blockDecompose(mat,map1,map2,S11,S12,S21,S22) ;
                        if(debug>=10)S11.Print();

                        map1.clear(); map2.clear();
                        vector<TString> scpPars; 
                        for(int i=0; i<cms->Get_max_uncorrelation();  i++){
                            if(TString(cms->Get_v_uncname()[i]).BeginsWith("Coupling_")){
                                map1.push_back(i);
                                scpPars.push_back(TString(cms->Get_v_uncname()[i]).ReplaceAll("Coupling_",""));
                            }
                            else map2.push_back(i);
                        } 
                        TMatrixDSym Scp;
                        RooMultiVarGaussian::blockDecompose(S11,map1,map2,Scp,S12,S21,S22) ;
                        Scp.Print();

                        TMatrixDSym *cormat_cp = (TMatrixDSym*)Scp.Clone();
                        for (Int_t i=0 ; i<cormat_cp->GetNrows() ; i++) {
                            for (Int_t j=0 ; j<cormat_cp->GetNcols() ; j++) {
                                if (i!=j) {
                                    (*cormat_cp)(i,j) = (*cormat_cp)(i,j) / sqrt((*cormat_cp)(i,i)*(*cormat_cp)(j,j)) ;
                                }
                            }
                        }
                        for (Int_t i=0 ; i<cormat_cp->GetNrows() ; i++) {
                            (*cormat_cp)(i,i) = 1.0 ;
                        }

                        cormat_cp->Print();

                        // Construct TH2D of correlation matrix 
                        Int_t n = cormat_cp->GetNcols() ;

                        TH2D* hh = new TH2D("CouplingCorrelationMatrix","Coupling Correlation Matrix",n,0,n,n,0,n) ;

                        for (Int_t i = 0 ; i<n ; i++) {
                            for (Int_t j = 0 ; j<n; j++) {
                                hh->Fill(i+0.5,n-j-0.5,(*cormat_cp)(i,j)) ;
                            }
                            hh->GetXaxis()->SetBinLabel(i+1,scpPars[i]);
                            hh->GetYaxis()->SetBinLabel(n-i,scpPars[i]);
                        }
                        hh->SetMinimum(-1) ;
                        hh->SetMaximum(+1) ;
                        TFile f(jobname+"_coupling.root","RECREATE");
                        f.WriteTObject(hh); f.Close();

                    }

                    /*




                    //_____________________________________________________________________________
                    TMatrixDSym RooFitResult::reducedCovarianceMatrix(const RooArgList& params) const 
                    {
                    // Return a reduced covariance matrix (Note that Vred _is_ a simple sub-matrix of V,
                    // row/columns are ordered to matched the convention given in input argument 'params'

                    const TMatrixDSym& V = covarianceMatrix() ;

                    // Handle case where V==Vred here
                    if (V.GetNcols()==params.getSize()) {
                    return V ;
                    }


                    // Make sure that all given params were floating parameters in the represented fit
                    RooArgList params2 ;
                    TIterator* iter = params.createIterator() ;
                    RooAbsArg* arg ;
                    while((arg=(RooAbsArg*)iter->Next())) {
                    if (_finalPars->find(arg->GetName())) {
                    params2.add(*arg) ;
                    } else {
                    coutW(InputArguments) << "RooFitResult::reducedCovarianceMatrix(" << GetName() << ") WARNING input variable " 
                    << arg->GetName() << " was not a floating parameters in fit result and is ignored" << endl ;
                    }
                    }
                    delete iter ;

                    // Need to order params in vector in same order as in covariance matrix
                    RooArgList params3 ;
                    iter = _finalPars->createIterator() ;
                    while((arg=(RooAbsArg*)iter->Next())) {
                    if (params2.find(arg->GetName())) {
                    params3.add(*arg) ;
                    }
                    }
                    delete iter ;

                    // Find (subset) of parameters that are stored in the covariance matrix
                    vector<int> map1, map2 ;
                    for (int i=0 ; i<_finalPars->getSize() ; i++) {
                    if (params3.find(_finalPars->at(i)->GetName())) {
                    map1.push_back(i) ;
                    } else {
                    map2.push_back(i) ;
                    }
                    }

                    TMatrixDSym S11, S22 ;
                    TMatrixD S12, S21 ;
                    RooMultiVarGaussian::blockDecompose(V,map1,map2,S11,S12,S21,S22) ;

                    return S11 ;
                    }
                     */

                }	
                if(idMH>0 && !DoNotRunGlobalFit){
                    vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();
                    double corr_ij = mat(0, idMH)/sqrt(mat(0,0)*mat(idMH,idMH));
                    cout<<corr_ij<<endl;
                    double x1, s1, x2, s2; 
                    myMinuit->GetParameter(0, x1, s1);
                    myMinuit->GetParameter(idMH, x2, s2);
                    RooEllipse *contour= new RooEllipse("contour",x2,x1,s2,s1,corr_ij);
                    contour->SetLineWidth(2) ;
                    for(int i=0; i<contour->GetN(); i++){
                        double r, m; contour->GetPoint(i, m,r);
                        contour -> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
                    }
                    TCanvas c("ellipse","ellipse", 600, 600);
                    contour->SetTitle("1 #sigma contour ;mass [GeV]; #mu");
                    contour->Draw("AC");
                    c.Print(jobname+"_ellipse_RM1sigma.pdf");
                    c.Print(jobname+"_ellipse_RM1sigma.root");

                    TFile f(jobname+"_ellipse_RM1sigma.root", "UPDATE");
                    double m[1] = { x2 * v_paramsUnc[idMH][1]+ v_paramsUnc[idMH][3] };
                    double merr[1] = { s2 * v_paramsUnc[idMH][1] };
                    double mu[1] = { x1 };
                    double muerr[1] = { s1 };
                    TGraphErrors tge(1, m, mu, merr, muerr);
                    tge.SetName("pointWithError");
                    f.WriteTObject(&tge);

                    double merrup[1] = { cms_global->POIs()[1].errUp};
                    double merrdown[1] = { fabs(cms_global->POIs()[1].errDown)};
                    double muerrup[1] = { cms_global->POIs()[0].errUp};
                    double muerrdown[1] = { fabs(cms_global->POIs()[0].errDown)};
                    TGraphAsymmErrors tgae(1, m, mu, merrdown, merrup, muerrdown, muerrup);
                    tgae.SetName("pointWithAsymError");
                    f.WriteTObject(&tgae);

                    TString snumbers5  = "corr_"; snumbers5+=corr_ij;
                    TObjString tsnumbers5(snumbers5);
                    tsnumbers5.Write("corr");
                    //f.WriteTObject(&tsnumbers5);
                    TString snumbers6  = ""; snumbers6+=y0_2;
                    TObjString tsnumbers6(snumbers6);
                    tsnumbers6.Write("fAmin");
                    f.Close();

                }

                cout<<" 	68% CL : "<<endl;
                cout<<" POI: mu = "<<cms_global->POIs()[0].value<<" + "<<cms_global->POIs()[0].errUp<<" - "<<fabs(cms_global->POIs()[0].errDown)<<endl;
                //if(idMH>0)cout<<" POI: MH = "<<cms_global->POIs()[1].value<<" + "<<cms_global->POIs()[1].errUp<<" - "<<fabs(cms_global->POIs()[1].errDown)<<endl;
                for(int k=1; k<cms_global->POIs().size(); k++)cout<<" POI: "<<cms_global->POIs()[k].name<<"= "<<cms_global->POIs()[k].value<<" + "<<cms_global->POIs()[k].errUp<<" - "<<fabs(cms_global->POIs()[k].errDown)<<endl;
                watch.Print();
            }

            if(debug) cout<<" YouAreHere debug 1"<<endl;

            if(cms->GetPhysicsModel()==typeCvCfHiggs)SaveResults(jobname+"_maxllfit", HiggsMass, 0, 0, 0, 0, 0, mu_hat_low, 0, mu_hat, mu_hat_up, 0);
            else if(cms->GetPhysicsModel()==typeC5Higgs)SaveResults(jobname+"_maxllfit", HiggsMass, 0, 0, 0, 0, 0, mu_hat_low, 0, mu_hat, mu_hat_up, 0);
            else SaveResults(jobname+"_maxllfit", HiggsMass, y0_2, 0, 0, 0, 0, mu_hat_low, 0, mu_hat, mu_hat_up, 0);

            // FIXME : jaehyeok     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            if(!ManualNuisanceFeeding) { 
                if(bDumpFitResults)cms->DumpFitResults(pars, jobname+"_fittedShape_floatMu"); 
                double *pars_maxll = new double[cms->Get_max_uncorrelation()+1]; 
                for(int i=0; i<=cms->Get_max_uncorrelation(); i++) pars_maxll[i]=pars[i];
                cms->Set_fittedParsInData_sb(pars_maxll);
            } else {  
                
                double *pars_new = new double[cms->Get_max_uncorrelation()+1];
                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++) pars_new[i] = pars[i]; 
                
                // get the file of nuisances   
                TString tmp[11];
                double central;
                bool SBfit = true;
                string line;
                cout << "FileNuisanceFeeding : " << FileNuisanceFeeding << endl; 
                ifstream inFile(FileNuisanceFeeding);
                int i_nuisance=1;
                if (inFile.is_open())
                {
                    while ( inFile.good() )
                    {
                        // get a line from the input file
                        getline (inFile,line);
                        if( line.find("Real")!=string::npos ) SBfit=false;
                        if(!SBfit) continue;

                        // skip line if it does not contain "rate"
                        if( line.find("par")==string::npos ) continue;
                        if( line.find("name")!=string::npos ) continue;
                        if( line.find("signal_strength")!=string::npos ) continue;
                        if( line.find(".....")!=string::npos ) continue;

                        cout << "************************* " << line << endl;

                        stringstream stream(line);
                        for(int i=0; i<11; i++) { 
                            if (i==2) stream >> central;
                            else stream >> tmp[i];
                        }

                        pars_new[i_nuisance] = central;
                        i_nuisance++;

                        if( !inFile.good() ) continue;
                    }
                }
                inFile.close();

                if(bDumpFitResults)cms->DumpFitResults(pars_new, jobname+"_fittedShape_floatMu");

                double *pars_maxll = new double[cms->Get_max_uncorrelation()+1]; 
                for(int i=0; i<=cms->Get_max_uncorrelation(); i++) pars_maxll[i]=pars_new[i];
                cms->Set_fittedParsInData_sb(pars_maxll);
            }
            // FIXME jaehyeok ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


            vector< vector<double> > v_Pars = cms_global->Get_v_Pars();	
            if(scanCvCf && cms_global->GetPhysicsModel()==typeCvCfHiggs){
                cms_global->SetNoErrorEstimation(1);
                TStopwatch watch2;
                watch2.Start();
                vector< std::pair<double, double> > vrm; vrm.clear(); 
                vector<double> vc; vc.clear();


                int _Cv_i = cms_global->Get_Cv_i();
                int _Cf_i = cms_global->Get_Cf_i();
                if(_Cv_i<0 or _Cf_i<0) {cout<<" no CV or CF parameter, exit"<<endl;  exit(1);}



                if(DoNotRunGlobalFit == false) { vrm.push_back(make_pair(pars[_Cv_i]*(v_Pars[_Cv_i][2]-v_Pars[_Cv_i][1])+v_Pars[_Cv_i][1],  pars[_Cf_i]*(v_Pars[_Cf_i][2]-v_Pars[_Cf_i][1])+v_Pars[_Cf_i][1]) );vc.push_back(y0_2-y0_2); }

                bool best = false;
                for(int i=0; i<vCv_toEval.size(); i++){
                    double testr = vCv_toEval[i];
                    for(int j=0; j<vCf_toEval.size(); j++){
                        if(debug)watch.Start();
                        double testm = vCf_toEval[j];
                        if(testm==0 and testr==0) continue;
                        _countPdfEvaluation = 0;
                        cms_global->SetPOItoBeFixed("CV",testr);
                        cms_global->SetPOItoBeFixed("CF",testm);


                        double y0_1 =0; 
                        if(bFixNuisancesAtBestFit){
                            bestFitPars[_Cv_i] =  (testr - v_Pars[_Cv_i][1])/(v_Pars[_Cv_i][2]-v_Pars[_Cv_i][1]);
                            bestFitPars[_Cf_i] =  (testr - v_Pars[_Cf_i][1])/(v_Pars[_Cf_i][2]-v_Pars[_Cf_i][1]);
                            y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);

                        }else y0_1 =  MinuitFit(3, tmp, tmperr, 1.0, pars, best?true:false, debug, 0, best?bestFitPars:0);

                        if(!bFixNuisancesAtBestFit){
                            for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                                myMinuit->GetParameter(i, tmp, tmpe);
                                bestFitPars[i]=tmp;
                            }
                            best = true;
                        }


                        if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                        if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                        vrm.push_back(make_pair(testr, testm));
                        vc.push_back(y0_1-y0_2);
                        TString sj = jobname; sj+="_fittedShape_cv"; sj+=testr; sj+="_cf"; sj+=testm;
                        if(bDumpFitResults)cms->DumpFitResults(pars, sj);
                        if(debug){
                            watch.Stop();
                            watch.Print();
                        }
                    }
                }
                if(debug){
                    printf("\n results of scanned Cv,Cf vs. q: \n");
                    for(int i=0; i<vrm.size(); i++){
                        printf("  [ Cv =%10.3f, Cf=%10.3f],  delta_q=%7.5f\n", vrm[i].first, vrm[i].second, vc[i]);
                    }
                }
                if(bPlots){
                    DrawTH2D d2h(sbinningCV, sbinningCF, vrm, vc, "-2lnQ; CV; CF", (jobname+"_scanned_cvcf_vs_q_TH2").Data(), pt);
                    d2h.draw();
                    d2h.save();
                    DrawTGraph2D d2d(vrm, vc, "-2lnQ; CV; CF", (jobname+"_scanned_cvcf_vs_q").Data(), pt);
                    d2d.draw();
                    d2d.save();
                    //TFile f(jobname+"_scanned_cvcf_vs_q.root", "UPDATE");
                    //d2d.getGraph()->SetName("Graph2D_cvcf_vs_q");					
                    //f.WriteTObject(d2d.getGraph());
                    //f.Close();
                }

                watch2.Stop();
                cout<<"Scanning total takes: "<<endl;
                watch2.Print();
                delete pars;
                delete cms;
                return 0;
            }
            if(debug) cout<<" YouAreHere debug 2  scanPOI="<<scanPOI<<endl;

            if(scanPOI && mainPOI!=""){
                if(debug) cout<<" YouAreHere scanPOI bFixNuisancesAtBestFit"<<bFixNuisancesAtBestFit<<endl;
                cms_global->SetNoErrorEstimation(1);
                TStopwatch watch2;
                watch2.Start();
                vector<double> vc; vc.clear();
                vector<double> vr; vr.clear();

                int poi_i = -1;
                for(int i=0; i<cms_global->Get_v_uncname().size(); i++) {
                    if(mainPOI==(TString)cms_global->Get_v_uncname()[i]) poi_i=i+1;
                }


                if(poi_i<0 ) {cout<<" no "<<mainPOI<<" parameter, exit"<<endl;  exit(1);}


                vector< vector<double> > v_Pars = cms_global->Get_v_Pars();	

                if(DoNotRunGlobalFit == false) { 
                    vr.push_back(pars[poi_i]*(v_Pars[poi_i][2]-v_Pars[poi_i][1])+v_Pars[poi_i][1]);
                    vc.push_back(y0_2-y0_2); 
                }

                bool best = false;
                for(int j=0; j<scanPOIrange.size(); j++){
                    if(debug)watch.Start();
                    double testm = scanPOIrange[j];
                    _countPdfEvaluation = 0;
                    cms_global->SetPOItoBeFixed(mainPOI,testm);


                    double y0_1 =0; 
                    if(bFixNuisancesAtBestFit){
                        if(debug) cout<<" MakeSureYouAreHere "<<endl;
                        bestFitPars[poi_i] =  (testm - v_Pars[poi_i][1])/(v_Pars[poi_i][2]-v_Pars[poi_i][1]);
                        y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);

                    }else y0_1 =  MinuitFit(3, tmp, tmperr, 1.0, pars, best?true:false, debug, 0, best?bestFitPars:0);

                    if(!bFixNuisancesAtBestFit){
                        for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                            myMinuit->GetParameter(i, tmp, tmpe);
                            bestFitPars[i]=tmp;
                        }
                        //best = true;
                    }

                    if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                    if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                    vr.push_back(testm);
                    vc.push_back(y0_1-y0_2);
                    if(debug){
                        watch.Stop();
                        watch.Print();
                    }
                }
                if(debug){
                    printf("\n results of scanned POI vs. q: \n");
                    for(int i=0; i<vr.size(); i++){
                        printf("  poi=%10.3f],  delta_q=%7.5f\n", vr[i], vc[i]);
                    }
                }
                if(bPlots){
                    TString st1 = "-2lnQ;"; st1+=mainPOI; st1+=";q"; 
                    TString st2 = jobname; st2+="_scanned_"; st2+=mainPOI; st2+="_vs_q_TH1";
                    DrawTH1D d1h(sbinningPOI, vr, vc, st1.Data(), st2.Data(), pt);
                    d1h.draw();
                    d1h.save();
                    st2 = jobname; st2+="_scanned_"; st2+=mainPOI; st2+="_vs_q";
                    DrawEvolution2D d2d(vr, vc, st1.Data(), st2.Data(), pt);
                    d2d.draw();
                    d2d.save();
                }

                watch2.Stop();
                cout<<"Scanning total takes: "<<endl;
                watch2.Print();
                delete pars;
                delete cms;
                return 0;
            }


            if(vCouplingsDef.size() and !scanMs) {delete pars;
                delete cms;
                return 0;}

                double y0_3, mu_hat2, mu_hat_up2, mu_hat_low2;

                if(bAlsoExtract2SigmErrorBar){
                    ErrorDef = TMath::ChisquareQuantile(0.95 , NumDegreeOfFreedom);// (confidenceLevel, ndf)
                    _ErrorDef = ErrorDef;
                    y0_3 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //102, 101, 21:  don't allow mu to be negative
                    //double y0_3 =  MinuitFit(102, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //102, 101, 21:  don't allow mu to be negative
                    cout<<y0_3<<" fitter u="<<pars[0]<<", from minos fit asymmetry 95% CL:  ["<<tmperr<<","<<tmpr<<"]"<<endl; // upperL, lowerL
                    mu_hat2 = pars[0];
                    mu_hat_up2 = tmpr;
                    mu_hat_low2 = tmperr;

                    {  // check covariance matrix and plot the ellipse if MH variable present in the floating parameter list 
                        int npar = cms_global->Get_max_uncorrelation()+1;
                        TMatrixDSym mat(npar); 
                        myMinuit->mnemat( mat.GetMatrixArray(), npar);
                        if(debug)mat.Print();

                        if(idMH>0){
                            vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();
                            double corr_ij = mat(0, idMH)/sqrt(mat(0,0)*mat(idMH,idMH));
                            cout<<corr_ij<<endl;
                            double x1, s1, x2, s2; 
                            myMinuit->GetParameter(0, x1, s1);
                            myMinuit->GetParameter(idMH, x2, s2);
                            RooEllipse *contour= new RooEllipse("contour",x2,x1,s2,s1,corr_ij);
                            contour->SetLineWidth(2) ;
                            for(int i=0; i<contour->GetN(); i++){
                                double r, m; contour->GetPoint(i, m, r);
                                contour -> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
                            }
                            TCanvas c("ellipse","ellipse", 600, 600);
                            contour->SetTitle("1 #sigma contour ;mass [GeV]; #mu");
                            contour->Draw("AC");
                            c.Print(jobname+"_ellipse_RM2sigma.pdf");
                            c.Print(jobname+"_ellipse_RM2sigma.root");

                            TFile f(jobname+"_ellipse_RM2sigma.root", "UPDATE");
                            double m[1] = { x2 * v_paramsUnc[idMH][1]+ v_paramsUnc[idMH][3] };
                            double merr[1] = { s2 * v_paramsUnc[idMH][1] };
                            double mu[1] = { x1 };
                            double muerr[1] = { s1 };
                            TGraphErrors tge(1, m, mu, merr, muerr);
                            tge.SetName("pointWithError");
                            f.WriteTObject(&tge);

                            double merrup[1] = { cms_global->POIs()[1].errUp};
                            double merrdown[1] = { fabs(cms_global->POIs()[1].errDown)};
                            double muerrup[1] = { cms_global->POIs()[0].errUp};
                            double muerrdown[1] = { fabs(cms_global->POIs()[0].errDown)};
                            TGraphAsymmErrors tgae(1, m, mu, merrdown, merrup, muerrdown, muerrup);
                            tgae.SetName("pointWithAsymError");
                            f.WriteTObject(&tgae);

                            TString snumbers5  = "corr_"; snumbers5+=corr_ij;
                            TObjString tsnumbers5(snumbers5);
                            tsnumbers5.Write("corr");
                            f.Close();
                        }
                        cout<<" 	95% CL : "<<endl;
                        cout<<" POI: mu = "<<cms_global->POIs()[0].value<<" + "<<cms_global->POIs()[0].errUp<<" - "<<fabs(cms_global->POIs()[0].errDown)<<endl;
                        if(idMH>0)cout<<" POI: MH = "<<cms_global->POIs()[1].value<<" + "<<cms_global->POIs()[1].errUp<<" - "<<fabs(cms_global->POIs()[1].errDown)<<endl;
                        watch.Print();
                    }

                }

                if(debug) cout<<" YouAreHere debug 3"<<endl;

                ErrorDef = TMath::ChisquareQuantile(0.68 , NumDegreeOfFreedom);// (confidenceLevel, ndf)
                _ErrorDef = ErrorDef;
                bool b_global_fit = true;
                // FIXME also need to check the 2 sigma fit is fine 
                // if not, need to run brute-force approach
                if( (success[0]!=0 /*or tmpr==tmperr*/) and !DoNotRunGlobalFit){
                    b_global_fit = false;
                    double tmpr = 0;
                    _initialRforFit = 0;
                    y0_2 =  MinuitFit(21, tmpr, tmperr, 0, pars, false, debug) ;
                    _initialRforFit = 1;
                    if(debug)	cout<<y0_2<<" fitter u="<<tmpr<<" +/- "<<tmperr<<endl;
                    mu_hat = pars[0];
                }


                if(debug) cout<<" YouAreHere debug 4"<<endl;
                if(b_global_fit==false){
                    double x1 =0, x2 =1;
                    double y0 = y0_2;  
                    x1 = mu_hat; x2=2*mu_hat;   // has to be positive, otherwise mess up


                    if(debug) cout<<" y0="<<y0<<" x1="<<x1<<" x2="<<x2<<endl;

                    double CI = ErrorDef/2.;
                    double precision = 0.001;
                    int nsearched = 3;
                    double y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                    if(fabs((y-y0)/2. - CI) > precision ){
                        //first, got a number > 1.921,  otherwise increase it by a factor of 10 ...
                        while( (y-y0)/2. < CI ){
                            x1 =  x2;
                            x2 *=10.;
                            y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                            if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                            nsearched++;
                        }
                        y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                        if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                        if(debug) cout<< " r="<<(x1+x2)/2.<<" NLL = "<<y<<endl;
                        while( fabs((y-y0)/2. - CI)/CI > precision ){
                            double  tmpx = (x1+x2)/2.;
                            if( (y-y0)/2. < CI ){
                                x1 = tmpx; 
                            }else{
                                x2 = tmpx;
                            }	
                            y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                            if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                            if(debug) printf(" r=%.6f,  NLL=%.10f\n", (x1+x2)/2., y);
                            nsearched++;
                        }
                        mu_hat_up= (x2+x1)/2.;
                    }else{
                        mu_hat_up = x2;
                    }

                }
                // look for asymmetry errors on muhat
                if(b_global_fit==false){
                    double x1 =0, x2 =1;
                    double y0 = y0_2;  
                    x1 = mu_hat; x2=0.5*mu_hat;  // has to be positive, otherwise mess up
                    double precision = 0.001;

                    double CI = ErrorDef/2.;
                    double y =  MinuitFit(3, tmp, tmp, 0, pars, false, debug );
                    if (fabs(y-y0_2)/2. <= CI+precision ){
                        mu_hat_low = 0.;
                    }else{
                        if(debug) cout<<" y0="<<y0<<" x1="<<x1<<" x2="<<x2<<endl;
                        int nsearched = 3;
                        y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                        if(fabs((y-y0)/2. - CI) > precision ){
                            //first, got a number > 1.921,  otherwise increase it by a factor of 10 ...
                            while( (y-y0)/2. < CI ){
                                x1 =  x2;
                                x2 /=10.;
                                y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                                if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                                nsearched++;
                            }
                            y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                            if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                            if(debug) cout<< " r="<<(x1+x2)/2.<<" NLL = "<<y<<endl;
                            while( fabs((y-y0)/2. - CI)/CI > precision ){
                                double  tmpx = (x1+x2)/2.;
                                if( (y-y0)/2. < CI ){
                                    x1 = tmpx; 
                                }else{
                                    x2 = tmpx;
                                }	
                                y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                                if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                                if(debug) printf(" r=%.6f,  NLL=%.10f\n", (x1+x2)/2., y);
                                nsearched++;
                            }
                            mu_hat_low= (x2+x1)/2.;
                        }else{
                            mu_hat_low = x2;
                        }
                    }

                }

                cout<<" muhat = "<<mu_hat<<" and 68% CL: ["<<mu_hat_low<<","<<mu_hat_up<<"]"<<endl;
                watch.Print();
                if(bPlots){
                    TGraphAsymmErrors *ge = new TGraphAsymmErrors(1);
                    ge->SetName("muhat");
                    ge->SetPoint(0, mu_hat, y0_2);
                    ge->SetPointError(0, mu_hat_low, mu_hat_up, 0, 0);
                    TFile f((jobname+"_muhat.root").Data(),"RECREATE");
                    f.WriteTObject(ge);

                    if(bAlsoExtract2SigmErrorBar){
                        TGraphAsymmErrors *ge2 = new TGraphAsymmErrors(1);
                        ge2->SetName("muhat_2sigma");
                        ge2->SetPoint(0, mu_hat2, y0_3);
                        ge2->SetPointError(0, mu_hat_low2, mu_hat_up2, 0, 0);
                        f.WriteTObject(ge2);
                    }
                }
                if(debug) cout<<" YouAreHere debug 5"<<endl;
                if(scanRs and scanMs){
                    if(debug) cout<<" YouAreHere debug 6"<<endl;
                    cms_global->SetNoErrorEstimation(1);
                    TStopwatch watch2;
                    watch2.Start();
                    vector< std::pair<double, double> > vrm; vrm.clear(); 
                    vector<double> vc; vc.clear();


                    if(idMH<0) {cout<<" no MH parameter, exit"<<endl;  exit(1);}


                    if(!DoNotRunGlobalFit){vrm.push_back(make_pair(pars[idMH]*v_paramsUnc[idMH][1]+v_paramsUnc[idMH][3], pars[0]) );vc.push_back(y0_2-y0_2);}

                    bool best = false;

                    for(int i=0; i<vR_toEval.size(); i++){
                        double testr = vR_toEval[i];
                        for(int j=0; j<vM_toEval.size(); j++){
                            if(debug)watch.Start();
                            double testm = vM_toEval[j];
                            _countPdfEvaluation = 0;
                            cms_global->SetPOItoBeFixed("MH",testm);
                            double y0_1 = 0;
                            if(bFixNuisancesAtBestFit){
                                bestFitPars[0] =  testr; 
                                bestFitPars[idMH] =  (testr - v_paramsUnc[idMH][3])/v_paramsUnc[idMH][1];
                                y0_1 =  MinuitFit(10, tmp, tmperr, testr, bestFitPars);

                            }else 
                                y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, best?true:false, debug, 0, best?bestFitPars:0);

                            if(!bFixNuisancesAtBestFit){
                                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                                    myMinuit->GetParameter(i, tmp, tmpe);
                                    bestFitPars[i]=tmp;
                                }
                                best = true;
                            }


                            if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                            if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                            vrm.push_back(make_pair(testm, testr));
                            vc.push_back(y0_1-y0_2);
                            TString sj = jobname; sj+="_fittedShape_r"; sj+=testr; sj+="_m"; sj+=testm;
                            if(bDumpFitResults)cms->DumpFitResults(pars, sj);
                            if(debug){
                                watch.Stop();
                                watch.Print();
                            }
                        }
                    }
                    if(debug){
                        printf("\n results of scanned r,m vs. q: \n");
                        for(int i=0; i<vrm.size(); i++){
                            printf("  [ m=%10.3f, r=%10.3f],  delta_q=%7.5f\n", vrm[i].first, vrm[i].second, vc[i]);
                        }
                    }
                    if(bPlots){
                        DrawTH2D d2h(sbinningMH, sbinningMU, vrm, vc, "#delta(-2lnQ); mass [GeV]; #mu", (jobname+"_scanned_rm_vs_q_TH2").Data(), pt);
                        d2h.draw();
                        d2h.save();

                        DrawTGraph2D d2d(vrm, vc, "#delta(2lnQ); mass [GeV]; #mu", (jobname+"_scanned_rm_vs_q").Data(), pt);
                        d2d.draw();
                        d2d.save();
                        TFile f(jobname+"_scanned_rm_vs_q.root", "UPDATE");
                        d2d.getGraph()->SetName("Graph2D_rm_vs_q");					
                        f.WriteTObject(d2d.getGraph());
                        f.Close();
                    }

                    watch2.Stop();
                    cout<<"Scanning total takes: "<<endl;
                    watch2.Print();
                    delete pars;
                    delete cms;
                    return 0;
                }
                else {
                    if(debug) cout<<" YouAreHere debug 7"<<endl;
                    cms_global->SetNoErrorEstimation(1);
                    watch.Start();
                    if(scanRs){
                        if(debug) cout<<" YouAreHere debug 8"<<endl;
                        vector<double> vr, vc; vr.clear(); vc.clear();
                        if(!DoNotRunGlobalFit){
                            vr.push_back(pars[0]);vc.push_back(y0_2-y0_2);
                        }
                        bool best = false;
                        for(int i=0; i<vR_toEval.size(); i++){
                            double testr = vR_toEval[i];
                            _countPdfEvaluation = 0;
                            double y0_1 = 0;
                            if(bFixNuisancesAtBestFit){
                                bestFitPars[0] =  testr; 
                                y0_1 =  MinuitFit(10, tmp, tmperr, testr, bestFitPars);

                            }else 
                                y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, best?true:false, debug, 0, best?bestFitPars:0);
                            if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                            if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                            vr.push_back(testr);
                            vc.push_back(y0_1-y0_2);
                            TString sj = jobname; sj+="_fittedShape_mu"; sj+=testr;
                            if(bDumpFitResults)cms->DumpFitResults(pars, sj);

                            if(!bFixNuisancesAtBestFit){
                                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                                    myMinuit->GetParameter(i, tmp, tmpe);
                                    bestFitPars[i]=tmp;
                                }
                                best = true;
                            }

                        }
                        if( scanRs >= 2) {
                            if(mu_hat<vR_toEval[0]){
                                double dr = (vR_toEval[0] - mu_hat)/10.;
                                if(dr<0.2) dr=0.2;
                                for(double testr = mu_hat+0.2; testr<vR_toEval[0]; testr+=dr){
                                    double y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, false, debug);
                                    if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                                    vr.push_back(testr);
                                    vc.push_back(y0_1-y0_2);
                                }
                            }
                            if(mu_hat_up>vR_toEval.back()){
                                double dr = (mu_hat_up - vR_toEval.back())/10.;
                                if(dr<0.2) dr=0.2;
                                for(double testr = vR_toEval.back(); testr<=mu_hat_up+0.2; testr+=dr){
                                    double y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, false, debug);
                                    if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                                    vr.push_back(testr);
                                    vc.push_back(y0_1-y0_2);
                                }
                            }
                        }

                        printf("\n results of scanned r vs. q: \n");
                        for(int i=0; i<vr.size(); i++){
                            printf("   r=%10.3f  delta_q=%7.5f\n", vr[i], vc[i]);
                        }
                        if(bPlots){
                            DrawEvolution2D d2d(vr, vc, "; r ; q", (jobname+"_scanned_r_vs_q").Data(), pt);
                            d2d.draw();
                            d2d.save();
                        }
                        watch.Stop();
                        watch.Print();
                        delete pars;
                        delete cms;
                        return 0;
                    }
                    if(scanMs){
                        if(debug) cout<<" MakeSureYouAreHere scanMs"<<endl;
                        bool best = false;
                        vector<double> vr, vc; vr.clear(); vc.clear();
                        if(!DoNotRunGlobalFit){
                            vr.push_back(pars[0]);vc.push_back(y0_2-y0_2);
                        }
                        if(1)for(int i=0; i<vM_toEval.size(); i++){
                            double testr = vM_toEval[i];
                            _countPdfEvaluation = 0;
                            cms_global->SetPOItoBeFixed("MH",testr);

                            //double y0_1 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, 1./*common mu*/, pars, best?true:false, debug, success, best?bestFitPars:0) ;  //202, 201, 2:  allow mu to be negative
                            double y0_1 =0;
                            if(bFixNuisancesAtBestFit){
                                //bestFitPars[idMH] =  (testr - v_paramsUnc[idMH][3])/v_paramsUnc[idMH][1];
                                //y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);
                                best = true;
                                cms->SetFastScan(true);
                                y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
                            }else
                                y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);

                            if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                            vr.push_back(testr);
                            vc.push_back(y0_1-y0_2);
                            TString sj = jobname; sj+="_fittedShape_m"; sj+=testr;
                            if(bDumpFitResults)cms->DumpFitResults(pars, sj);

                            //if(!bFixNuisancesAtBestFit)
                            {
                                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                                    myMinuit->GetParameter(i, tmp, tmpe);
                                    bestFitPars[i]=tmp;
                                }
                                best = true;
                            }

                        }

                        if(0){
                            double *bestFitParsBak= new double[cms->Get_max_uncorrelation()+1];
                            for(int i=0; i<=cms->Get_max_uncorrelation(); i++) bestFitParsBak[i]=bestFitPars[i];
                            for(int i=0; i<100; i++){
                                double testr = 126.3+0.01*i;
                                cout<<" testm = "<<testr<<endl;
                                _countPdfEvaluation = 0;
                                cms_global->SetPOItoBeFixed("MH",testr);

                                //double y0_1 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, 1./*common mu*/, pars, best?true:false, debug, success, best?bestFitPars:0) ;  //202, 201, 2:  allow mu to be negative
                                double y0_1 =0; 
                                if(bFixNuisancesAtBestFit){
                                    //bestFitPars[idMH] =  (testr - v_paramsUnc[idMH][3])/v_paramsUnc[idMH][1];
                                    //y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);
                                    best = true;
                                    cms->SetFastScan(true);
                                    y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
                                }else 
                                    y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);

                                if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                                vr.push_back(testr);
                                vc.push_back(y0_1-y0_2);
                                TString sj = jobname; sj+="_fittedShape_m"; sj+=testr;
                                if(bDumpFitResults)cms->DumpFitResults(pars, sj);

                                //if(!bFixNuisancesAtBestFit)
                                {
                                    for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                                        myMinuit->GetParameter(i, tmp, tmpe);
                                        bestFitPars[i]=tmp;
                                    }
                                    best = true;
                                }

                            }
                            for(int i=0; i<=cms->Get_max_uncorrelation(); i++) bestFitPars[i]=bestFitParsBak[i];
                            for(int i=0; i<100; i++){
                                double testr = 126.3-0.01*i;
                                _countPdfEvaluation = 0;                                    
                                cms_global->SetPOItoBeFixed("MH",testr);

                                //double y0_1 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, 1./*common mu*/, pars, best?true:false, debug, success, best?bestFitPars:0) ;  //202, 201, 2:  allow mu to be negative                                             
                                double y0_1 =0;                 
                                if(bFixNuisancesAtBestFit){     
                                    //bestFitPars[idMH] =  (testr - v_paramsUnc[idMH][3])/v_paramsUnc[idMH][1];
                                    //y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);
                                    best = true;            
                                    cms->SetFastScan(true); 
                                    y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
                                }else                           
                                    y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);

                                if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                                vr.push_back(testr);            
                                vc.push_back(y0_1-y0_2);        
                                TString sj = jobname; sj+="_fittedShape_m"; sj+=testr;
                                if(bDumpFitResults)cms->DumpFitResults(pars, sj);

                                //if(!bFixNuisancesAtBestFit)   
                                {                               
                                    for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                                        myMinuit->GetParameter(i, tmp, tmpe);
                                        bestFitPars[i]=tmp;
                                    }
                                    best = true;
                                }

                            }
                        }

                        printf("\n results of scanned m vs. q: \n");
                        for(int i=0; i<vr.size(); i++){
                            printf("   r=%10.3f  delta_q=%7.5f\n", vr[i], vc[i]);
                        }
                        if(bPlots){
                            TString st1 = "-2lnQ;"; st1+="Higgs boson mass (GeV)"; st1+=";q"; 
                            TString st2 = jobname; st2+="_scanned_m_vs_q_TH1";
                            DrawTH1D d1h(sbinningMH, vr, vc, st1.Data(), st2.Data(), pt);
                            d1h.draw();
                            d1h.save();
                            DrawEvolution2D d2d(vr, vc, "; m ; q", (jobname+"_scanned_m_vs_q").Data(), pt);
                            d2d.draw();
                            d2d.save();
                        }

                        watch.Stop();
                        watch.Print();
                        delete pars;
                        delete cms;
                        return 0;
                    }
                }
                if(doExpectation){
                    MLLxsBands lb(cms);	
                    lb.SetDebug(debug);
                    lb.SetSignalStrength(ToysAtDifferentSignalStrength, InjectingSignalStrength);
                    int noutcomes = toys;
                    if(sToysFromFile!="")lb.ToysFromFile(sToysFromFile);
                    lb.Bands(noutcomes);
                    double rmean, rm1s, rm2s, rp1s, rp2s;	
                    rmean=lb.GetMLLxsMean();
                    rm1s=lb.GetMLLxs(-1);rm2s=lb.GetMLLxs(-2);rp1s=lb.GetMLLxs(1);rp2s=lb.GetMLLxs(2);
                    difflimits=lb.GetDifferentialMLLxs();

                    TString ts=jobname; ts+="_maxlls";
                    FillTree(ts, lb.Get_mu_unsorted(), "muhat", lb.Get_lh_unsorted(), "likelihood");
                    CopyTrees(jobname+"_maxllfit.root", "T", ts+"_tree.root", "T_data");

                    TCanvas *c=new TCanvas("csig","cSig");
                    c->SetLogy(1);
                    TH1F *h=new TH1F("h","; #hat{#mu}; entries", 200, 0, 8);	
                    for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
                    h->Draw();
                    Save(c, (jobname+"differential_muhat").Data());

                    vector<double> all_calculated_R95s;
                    vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
                    SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

                    pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
                    pt->AddText("#hat{#mu} statistical bands");
                    string ssave=(jobname+"_cummulativeMuhat").Data();
                    string stitle="; #hat{#mu}; cumulative probability;";
                    PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
                            rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
                    plotRvsP.draw();

                    cout<<" Expected MLLxs:   median       = "<<rmean<<endl;
                    cout<<" Expected MLLxs:   68\% coverage = ["<<rm1s<< "," <<rp1s<<"]"<<endl;
                    cout<<" Expected MLLxs:   95\% coverage = ["<<rm2s<< "," <<rp2s<<"]"<<endl;
                }


        }else{
            if(doExpectation && sFileLimitsDistribution!="") saveExpectation();
        }



    }else { // calc significances 

        // initialize the calculator
        CLsBase frequentist;
        frequentist.SetDebug(debug);
        frequentist.SetRdm(rdm);
        frequentist.SetTestStatistics(testStat);

        double tmp, tmperr;
        double rmean, rm1s, rm2s, rp1s, rp2s;	
        vector<double> difflimits; 
        vector<double> muhats; 
        if(method == "ProfiledLikelihood" or method=="ProfileLikelihood"){
            cms_global= cms;
            vdata_global=cms->Get_v_data();
            vdata_global_th=cms_global->Get_v_dataTH();
            if(cms->hasParametricShape()){
                cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
            }

            double pvalue, significance, muhat;
            runPValue(pvalue, significance, muhat, true);
            if(doExpectation){
                TFile *fin;
                if(sToysFromFile!=""){
                    toys = cms->NumberOfToysFromFile(sToysFromFile);
                    fin = new TFile(sToysFromFile);
                }
                for(int t=0; t<toys; t++)
                {
                    if(debug>=0)cout<<" Toy "<<t<<endl;
                    if(sToysFromFile!=""){
                        vdata_global = (VDChannel)cms->GetToyDataFromFile(fin, t);
                        if(cms->hasParametricShape()){
                            cms->SetTmpDataForUnbinned(cms->GetToyUnbinnedDataFromFile(fin, t));
                        }
                    }
                    else{
                        vdata_global = (VDChannel)cms->GetToyData_H1();//_cms->Get_fittedParsInData_sb());
                        if(cms->hasParametricShape()){
                            cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_toy());
                        }
                    }
                    runPValue(pvalue, significance, muhat, false);
                    difflimits.push_back(significance);
                    muhats.push_back(muhat);
                }

                /*
                   SignificanceBands lb(&frequentist, cms);	
                   lb.SetDebug(debug);
                   int noutcomes = toys;
                   lb.Bands(noutcomes);
                   rmean=lb.GetSignificanceMean();
                   rm1s=lb.GetSignificance(-1);rm2s=lb.GetSignificance(-2);rp1s=lb.GetSignificance(1);rp2s=lb.GetSignificance(2);
                   difflimits=lb.GetDifferentialSignificances();
                   muhats=lb.GetDifferentialMuhats();
                 */
            }
        }else if(method == "Hybrid"){
            cout<<" calc p-value ... "<<endl;

            cms->SetTossToyConvention(tossToyConvention);
            cms->SetSignalScaleFactor(1.);
            cms->SetUseBestEstimateToCalcQ(UseBestEstimateToCalcQ);
            frequentist.SetTestStatistics(testStat);
            frequentist.SetModel(cms);
            vector<double> vclb;

            double  signi, pvalue;
            if(sFileBonlyM2lnQ==""){
                int ntoysToDoSignificance = toysHybrid; //10000000;
                if(tossToyConvention==1)frequentist.checkFittedParsInData(bReadPars, bWritePars, fileFittedPars);
                cout<<"\n Running "<<ntoysToDoSignificance<<" toys to evaluate significance for data "<<endl;
                signi = frequentist.SignificanceForData(ntoysToDoSignificance);

                //cout<<"Q_b_data = "<<frequentist.Get_m2lnQ_data()<<endl;	
                printf("Q_b_data = %.10f\n",frequentist.Get_m2lnQ_data());	
                vclb = frequentist.Get_m2logQ_b();
                TString  s = jobname; 
                s+="_hybridSig_ts"; s+=testStat;
                s+="_seed"; s+=seed;
                TString s_qdata = TString::Format("%.10f",frequentist.Get_m2lnQ_data());
                FillTree(s, vclb, s_qdata);
                //signi = frequentist.SignificanceForData(-s_qdata.Atof(), vclb);
                //cout<<"  "<<s_qdata.Atof()<<"  "<<frequentist.Get_m2lnQ_data()<<endl;

                //for(int i=0; i<vclb.size(); i++){
                //printf("Q_b_toy = %.10f\n", vclb[i]);
                //}
            }else{
                cout<<"\n reading toys to evaluate significance for data "<<endl;
                TString s_qdata = "";
                TTree * tb = (TTree*)LoadTreeBonly(sFileBonlyM2lnQ, s_qdata);				
                vclb = GetVectorFrom(tb, "brT");
                signi = frequentist.SignificanceForData(-s_qdata.Atof(), vclb);
            }

            cout<<"------------------------------------------------------------"<<endl;
            if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
            cout<<" Observed Significance for the data = "<<signi<<endl;
            cout<<"------------------------------------------------------------"<<endl;

            pvalue = ROOT::Math::normal_cdf_c(signi);
            SaveResults(jobname+"_frqObsPvalue", HiggsMass, 0, 0, signi, pvalue, 0, 0, 0, 0, 0, 0);
            if(doExpectation){
                // FIXME
            }
        }
        if(doExpectation){

            TString ts=jobname; ts+="_significances";
            FillTree(ts, difflimits, "significance", muhats, "muhat");

            TCanvas *c=new TCanvas("csig","cSig");
            c->SetLogy(1);
            TH1F *h=new TH1F("h",";ProfiledLikelihood significance; entries", 200, 0, 8);	
            for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
            h->Draw();
            Save(c, (jobname+"differential_sig").Data());

            vector<double> all_calculated_R95s;
            vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
            SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

            pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
            pt->AddText("PL significance statistical bands");
            string ssave=(jobname+"_cummulativeS").Data();
            string stitle="; ProfiledLikelihood significance; cumulative probability;";
            PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
                    rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
            plotRvsP.draw();

            cout<<" Expected Significance:   median       = "<<rmean<<endl;
            cout<<" Expected Significance:   68\% coverage = ["<<rm1s<< "," <<rp1s<<"]"<<endl;
            cout<<" Expected Significance:   95\% coverage = ["<<rm2s<< "," <<rp2s<<"]"<<endl;
        }
    }

#ifdef LANDS_PROFILE
    if (gSystem->Getenv("CPUPROFILE"))
    {
        ProfilerStop();
        printf("Finished profiling into '%s'\n", gSystem->Getenv("CPUPROFILE"));
    }
#endif

    watch.Print();
    delete cms;
    return 1;
}
vector<double> GetListToEvalForNuis(TString& parname, vector<TString> tmpv){
    //for(int i=0; i<tmpv.size(); i++) cout<<tmpv[i]<<endl;
    parname = tmpv[0]; 
    tmpv.erase(tmpv.begin());
    return GetListToEval(parname, tmpv);
}
vector<double> GetListToEval(TString parname, vector<TString> tmpv){
    // tmpv[0] for binning,  tmpv[1] and so on for actual scan if tmpv.size>1
    // cout<<"tmpv.size="<<tmpv.size()<<endl;
    for(int i=0; i<tmpv.size(); i++) cout<<tmpv[i]<<endl;
    vector<double> vtoEval;vtoEval.clear();
    for(int i= (tmpv.size()>1?1:0); i<tmpv.size(); i++){
        if(tmpv[i].IsFloat()) vtoEval.push_back(tmpv[i].Atof());
        else if(tmpv[i].BeginsWith("[") and tmpv[i].EndsWith("]")){
            tmpv[i].ReplaceAll("[","");tmpv[i].ReplaceAll("]","");	
            vector<string> vstr;
            StringSplit(vstr, tmpv[i].Data(), ",");
            if(vstr.size()==3 and (TString(vstr[2]).IsFloat() or TString(vstr[2]).BeginsWith("x")) ){
                double r0 = TString(vstr[0]).Atof(), r1=TString(vstr[1]).Atof();	
                //if(r0>r1 or r0<0) continue;
                if(TString(vstr[2]).BeginsWith("x")) {
                    //if(r0==0) continue;
                    double step = TString(vstr[2]).ReplaceAll("x", "").Atof();
                    if(step<=1) continue;
                    for(double r=r0; r<=r1; r*=step) vtoEval.push_back(r);
                }else{
                    double step = TString(vstr[2]).Atof();
                    if(step<=0) continue;
                    for(double r=r0; r<=r1; r+=step) vtoEval.push_back(r);
                };
            }else{
                cout<<"ERROR: wrong format of "<<parname<<", should be sth like [1.2,2.0,x1.05] or [1.2,2.0,0.05]"<<endl; exit(1);
            };
        }else {cout<<"ERROR: wrong format of args of "<<parname<<":  \""<<tmpv[i]<<"\""<<endl; exit(1);}
    }
    std::sort(vtoEval.begin(), vtoEval.end());
    vector<double>::iterator it;
    it = std::unique(vtoEval.begin(), vtoEval.end());
    vtoEval.resize(it - vtoEval.begin());
    cout<<parname<<"_toEval"<<vtoEval.size()<<endl;
    return vtoEval;
}
double Err_Bisection(CountingModel*cms, TString spar, double * bestFitPars, double y0_2, double glbmin, int hilo){
    double *pars = new double[cms->Get_max_uncorrelation()+1]; // nsys + r
    double tmp, tmpe;
    vector< vector<double> > v_Pars = cms_global->Get_v_Pars();	
    int par_i= cms->Get_par_i(spar);
    int iteration=1, maxIteration=50;	
    double precision = 1e-3;
    double deltaq = 9999.;
    bool best = true;
    if(hilo==1){
        double err=v_Pars[par_i][2]-glbmin;
        if(err<=1e-6) { cout<<spar<<" up error at limit !"<<endl; delete pars; return err;}
        cms_global->SetPOItoBeFixed(spar, glbmin+err);
        double y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
        deltaq = y0_1 - (y0_2+1);
        double xmin = 0, xmax=err;
        if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
        if(deltaq>precision){
            while(deltaq > precision and iteration<=maxIteration) {
                xmax/=2.;
                err = xmax;
                cms_global->SetPOItoBeFixed(spar, glbmin+err);
                y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
                deltaq = y0_1 - (y0_2+1);
                iteration++;	
                if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
            }
            xmin = err; xmax*=2;
        }
        else if(deltaq<-precision){
            cout<<spar<<" up error at limit !"<<endl; delete pars; return err;
            /*while(deltaq<-precision and iteration<=maxIteration) {
              xmin = xmax;xmax = xmax*2;
              err = xmax;
              cms_global->SetPOItoBeFixed(spar, glbmin+err);
              y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
              deltaq = y0_1 - (y0_2+1);
              iteration++;	
              if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
              }
              xmin=xmax/2.;xmax = err;*/ 
        }
        if(debug)cout<<spar<<" err xmin="<<xmin<<" "<<xmax<<endl;

        err = xmin + (xmax-xmin)/2.;
        cms_global->SetPOItoBeFixed(spar, glbmin+err);
        y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
        deltaq = y0_1 - (y0_2+1);
        if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
        while(fabs(deltaq)>precision and iteration<=maxIteration) {
            if(deltaq>0) { xmax=err; err=(xmax-xmin)/2.+xmin; }
            else {xmin=err; err= (xmax-xmin)/2.+xmin; }
            cms_global->SetPOItoBeFixed(spar, glbmin+err);
            y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
            deltaq = y0_1 - (y0_2+1);
            iteration++;	
            if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
        }
        if(debug)cout<<spar<<" err final = "<<err<<endl;
        delete pars;
        return err;
    }
    else {
        double err=glbmin-v_Pars[par_i][1];
        if(err<=1e-6) { cout<<spar<<" low error at limit !"<<endl;delete pars; return err;}
        cms_global->SetPOItoBeFixed(spar, glbmin-err);
        double y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
        deltaq = y0_1 - (y0_2+1);
        double xmin = 0, xmax=err;
        if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
        if(deltaq>precision){
            while(deltaq > precision and iteration<=maxIteration) {
                xmax/=2.;
                err = xmax;
                cms_global->SetPOItoBeFixed(spar, glbmin-err);
                y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
                deltaq = y0_1 - (y0_2+1);
                iteration++;	
                if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
            }
            xmin = err; xmax*=2;
        }
        else if(deltaq<-precision){
            cout<<spar<<" low error at limit !"<<endl;delete pars; return err;
        }
        if(debug)cout<<spar<<" err xmin="<<xmin<<" "<<xmax<<endl;

        err = xmin + (xmax-xmin)/2.;
        cms_global->SetPOItoBeFixed(spar, glbmin-err);
        y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
        deltaq = y0_1 - (y0_2+1);
        if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
        while(fabs(deltaq)>precision and iteration<=maxIteration) {
            if(deltaq>0) { xmax=err; err=(xmax-xmin)/2.+xmin; }
            else {xmin=err; err= (xmax-xmin)/2.+xmin; }
            cms_global->SetPOItoBeFixed(spar, glbmin-err);
            y0_1 =  MinuitFit(3, tmp, tmpe, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);
            deltaq = y0_1 - (y0_2+1);
            iteration++;	
            if(debug)cout<<spar<<" err = "<<err<<" deltaQ = "<<deltaq<<endl;
        }
        if(debug)cout<<spar<<" err final = "<<err<<endl;
        delete pars;
        return err;
    }
}

bool runPValue(double &pvalue, double& significance, double& muhat, bool printInfo){
    CountingModel*cms= cms_global;
    _inputNuisances = cms->Get_norminalPars();
    _startNuisances = cms->Get_norminalPars();

    double tmp, tmperr;
    double rmean, rm1s, rm2s, rp1s, rp2s;	
    double sig_data;
    double m2lnQ;
    double x2;
    /*
       x2 =  MinuitFit(2, tmp, tmperr); // allow mu<0
       cout<<"fitted r = "<<tmp<<endl;
       m2lnQ = MinuitFit(3,tmp, tmp) - x2; // mu=0
       double sig_data = sqrt(fabs(m2lnQ));
       cout<<"Observed significance using PLR method = "<<sig_data<<endl;
       cout<<"Observed p-value = "<<1-TMath::Erf(sig_data)<<endl;
     */
    //x2 =  MinuitFit(21, tmp, tmperr); // mu>0
    //x2 =  MinuitFit(2, tmp, tmperr); // allow mu<0
    double * fittedPars = new double[cms->Get_max_uncorrelation()+1];
    int success[1];success[0]=0;
    _countPdfEvaluation = 0;
    x2 = MinuitFit(2, tmp, tmperr, 1, fittedPars, false, debug, success);  // MinuitFit(mode, r, err_r)
    if(debug) cout<<" _countPdfEvaluation ="<<_countPdfEvaluation<<endl;
    /*
       int ntmp = 0;
       while(success[0]!=0 && ntmp < 10){
       cms->FluctuatedNumbers(); // toss nuisance around fitted b_hat_0 in data  
       _startNuisances = cms->Get_randomizedPars();	
       x2 = MinuitFit(2, tmp, tmperr, 0, fittedPars, false, 0, success);  // MinuitFit(mode, r, err_r)
       ntmp++;
       }
       if(ntmp>0)cout<<" after "<<ntmp+1<<" tries: fit "<<(success[0]?"fails":"successful")<<endl;
     */
    _initialRforFit = 0;
    if(success[0])x2 =  MinuitFit(21, tmp, tmperr, 0, fittedPars, false, debug, success); // mu>0
    _initialRforFit = 1;
    if(printInfo)cout<<"fitted r = "<<tmp<<endl;
    double mu_hat = tmp;

    double tmp_muhat, tmp_muhaterr ;
    myMinuit->GetParameter(0, tmp_muhat, tmp_muhaterr);

    if(debug)cout<<" mu_hat = "<<mu_hat<<" tmp_muhat="<<tmp_muhat<<" tmp_muhaterr="<<tmp_muhaterr<<endl;

    if(mu_hat<=0) {
        if(printInfo)cout<<" fitted r <=0 at maximum likelihood ratio ,  signicance = 0"<<endl;
        m2lnQ = 0;
    }else{
        _countPdfEvaluation = 0;
        double * fittedPars2 = new double[cms->Get_max_uncorrelation()+1];
        m2lnQ = MinuitFit(3,tmp, tmp, 0, fittedPars2, true, debug, success, fittedPars) - x2; // mu=0
        if(debug)cout<<"   x2 = "<<x2<<",  m2lnQ = "<<m2lnQ<<endl;
        if(debug) cout<<" _countPdfEvaluation ="<<_countPdfEvaluation<<endl;
        delete fittedPars2;
    }
    sig_data = sqrt(fabs(m2lnQ));
    pvalue = ROOT::Math::normal_cdf_c(sig_data);
    significance = sig_data;
    muhat = mu_hat;
    if(printInfo){ 
        if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
        cout<<"Observed significance using PLR method = "<<sig_data<<endl;
        if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
        cout<<"Observed p-value = "<<pvalue<<endl;
        if(sig_data>=4 ) cout<<"WARNING: please contact mschen@cern.ch for a potential issue on the nuisances' range. Needs change from [-5,5] to larger one"<<endl;
        if(printInfo)SaveResults(jobname+"_plrObsPvalue", HiggsMass, 0, 0, sig_data, pvalue, 0, tmp_muhat-tmp_muhaterr, tmp_muhat, mu_hat, tmp_muhat+tmp_muhaterr, 0);
    }
    if(debug>=10) { // show a plot for   -log(Lambda(mu)) vs. mu ...
        for(double r=0; r<2; r+=0.1){
            m2lnQ = MinuitFit(3,tmp, tmp, r) -x2; 
        }
    }
}
bool runMaxLikelihoodFit(bool printInfo, bool makePlots, TString smakePlots){

    TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
    TStopwatch watch;  
    watch.Start();
    CountingModel* cms = cms_global;
    int idMH=-1;
    vector<string> v_uncname = cms_global->Get_v_uncname();
    for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]=="MH") idMH=i;
    vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();

    if(scanNuisance!=""){
        int idNuis=-1;
        vector<string> v_uncname = cms_global->Get_v_uncname();
        for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]==scanNuisance.Data()) idNuis=i;
        if(idNuis<0) {cout<<"ERROR:  intend to scan "<<scanNuisance<<", but it's not in the list, exit"<<endl; return 1;}
    }
    //vdata_global=cms->Get_v_data();
    //vdata_global_th=cms_global->Get_v_dataTH();

    cms_global->SetDebug(debug);
    _inputNuisances = cms->Get_norminalPars();
    _startNuisances = cms->Get_norminalPars();

    if(bDumpFitResults)cms->DumpFitResults(_inputNuisances, jobname+"_nominalShape");

    double *pars = new double[cms->Get_max_uncorrelation()+1]; // nsys + r
    for(int i=0; i<=cms->Get_max_uncorrelation(); i++) pars[i]=0;

    //double bestFitPars[cms->Get_maximumFunctionCallsInAFit()+1];  //  not static ... not allowed 
    double *bestFitPars= new double[cms->Get_max_uncorrelation()+1];
    for(int i=0; i<=cms->Get_max_uncorrelation(); i++) bestFitPars[i]=0;

    double tmp=0, tmpe=0, tmpr = 0, tmperr=0;
    double ErrorDef = TMath::ChisquareQuantile(0.68 , NumDegreeOfFreedom);// (confidenceLevel, ndf)
    _ErrorDef = ErrorDef;
    int success[1]={0};
    _countPdfEvaluation = 0;
    if(bDumpFitResults){
        lands::_bDumpFinalFitResults = true;
    }
    double y0_2=0;
    // for CvCf   ==>   y0_2,  CV, CF and their errors 


    if(mainPOI != "") {
        cms->keepOnlyPOI(mainPOI);	
    }
    if(DoNotRunGlobalFit==false){
        if(NoErrorEstimate or ErrEstAlgo=="Bisect") cms->SetNoErrorEstimation(1);
        if(vCouplingsDef.size()==0) y0_2 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //202, 201, 2:  allow mu to be negative
        else { y0_2 =  MinuitFit(3, tmp, tmperr, 1.0, pars, false, debug, 0, 0); }
        if(debug) cout<<" _countPdfEvaluation = "<<_countPdfEvaluation<<endl;
        if(printInfo) cout<<y0_2<<" fitter u="<<pars[0]<<", from minos fit asymmetry 68% CL:  ["<<tmperr<<","<<tmpr<<"]"<<endl; // upperL, lowerL
        fitBestMu = pars[0]; fitBestMuErrLo=pars[0]-tmperr; fitBestMuErrHi = tmpr-pars[0];
        for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
            myMinuit->GetParameter(i, tmp, tmpe);
            bestFitPars[i]=tmp;
        }
        if(printInfo)cout<<" fAmin = "<<myMinuit->fAmin<<endl;
        vdPseudoPOIs.clear();
        for(int p=0; p<viPseudoPOIs.size(); p++){
            for(int i=0; i<cms->Get_max_uncorrelation()+1; i++) {
                if(i==viPseudoPOIs[p]) vdPseudoPOIs.push_back(pars[i]); 
            }
        }
    }
    //double y0_2 =  MinuitFit(102, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //102, 101, 21:  don't allow mu to be negative
    double mu_hat = pars[0];
    double mu_hat_up = tmpr;
    double mu_hat_low = tmperr;
    if(printInfo)watch.Print();

    if(scanNuisance!=""){
        int idNuis=-1;
        vector<string> v_uncname = cms_global->Get_v_uncname();
        for(int i=1; i<=v_uncname.size(); i++) if(v_uncname[i-1]==scanNuisance.Data()) idNuis=i;
        if(idNuis<0) return 1;

        bool best = false;
        vector<double> vr, vc; vr.clear(); vc.clear();
        if(!DoNotRunGlobalFit){
            vr.push_back(pars[idNuis]);vc.push_back(y0_2-y0_2);
        }
        for(int i=0; i<vNuis_toEval.size(); i++){
            double testr = vNuis_toEval[i];
            _countPdfEvaluation = 0;
            double y0_1 =0; 

            best=true;
            bestFitPars[idNuis] =  testr ;
            y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);

            if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
            vr.push_back(testr);
            vc.push_back(y0_1-y0_2);
            TString sj = jobname; sj+="_fittedShape_m"; sj+=testr;
            if(bDumpFitResults)cms->DumpFitResults(pars, sj);

            if(!bFixNuisancesAtBestFit){
                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                    myMinuit->GetParameter(i, tmp, tmpe);
                    bestFitPars[i]=tmp;
                }
                best = true;
            }

        }
        if(printInfo){
            printf("\n results of scanned %s vs. q: \n", scanNuisance.Data());
            for(int i=0; i<vr.size(); i++){
                printf("   r=%10.3f  delta_q=%7.5f\n", vr[i], vc[i]);
            }
        }
        if(bPlots and makePlots){
            TString st1 = "-2lnQ;"; st1+=scanNuisance; st1+=";q"; 
            TString st2 = jobname; st2+="_scanned_"; st2+=scanNuisance; st2+="_vs_q_TH1"; st2+=smakePlots;
            //DrawTH1D d1h(sbinningMH, vr, vc, st1.Data(), st2.Data(), pt);
            //d1h.draw();
            //d1h.save();
            DrawEvolution2D d2d(vr, vc, TString("; "+scanNuisance+" ; q").Data(), (jobname+"_scanned_"+scanNuisance+"_vs_q"+smakePlots).Data(), pt);
            d2d.draw();
            d2d.save();
        }

        watch.Stop();
        if(printInfo)watch.Print();
        delete pars;
        return 0;
    }


    // mass error
    if(idMH>0){
        if(ErrEstAlgo =="Minos" and printInfo and !NoErrorEstimate){
            cout<<" 	68% CL : "<<endl;
            for(int k=1; k<cms_global->POIs().size(); k++){
                cout<<" POI: "<<cms_global->POIs()[k].name<<"= "<<cms_global->POIs()[k].value<<" + "<<cms_global->POIs()[k].errUp<<" - "<<fabs(cms_global->POIs()[k].errDown)<<endl;
                if(cms_global->POIs()[k].name== "MH"){ fitBestMH =cms_global->POIs()[k].value ; fitBestMHErrLo=fabs(cms_global->POIs()[k].errDown); fitBestMHErrHi =cms_global->POIs()[k].errUp ;}
            }
        }
        if(ErrEstAlgo=="Bisect"){
            double glbminMH = cms->Get_v_Pars()[cms->Get_par_i("MH")][0];
            double errhiMH= Err_Bisection(cms, "MH", bestFitPars, y0_2, glbminMH, 1);
            double errloMH= Err_Bisection(cms, "MH", bestFitPars, y0_2, glbminMH, 0);
            if(printInfo){
                cout<<" 	68% CL : "<<endl;
                cout<<" POI MH: "<<glbminMH<<" +"<<errhiMH<<" -"<<errloMH<<endl;
            }
            fitBestMH = glbminMH; fitBestMHErrLo=errloMH; fitBestMHErrHi = errhiMH;
        }
    }

    if(cms->GetPhysicsModel()==typeCvCfHiggs and !DoNotRunGlobalFit) cms->ShowCvCfHiggsScales(bestFitPars);
    if(cms->GetPhysicsModel()==typeC5Higgs ){
        if(ErrEstAlgo =="Minos" and printInfo){
            cout<<" 	68% CL : "<<endl;
            for(int k=1; k<cms_global->POIs().size(); k++)cout<<" POI: "<<cms_global->POIs()[k].name<<"= "<<cms_global->POIs()[k].value<<" + "<<cms_global->POIs()[k].errUp<<" - "<<fabs(cms_global->POIs()[k].errDown)<<endl;
        }
        if(ErrEstAlgo=="Bisect"){
            double glbminCgg = cms->Get_v_Pars()[cms->Get_par_i("Cgg")][0];
            double glbminCbb = cms->Get_v_Pars()[cms->Get_par_i("Cbb")][0];
            double glbminCtt = cms->Get_v_Pars()[cms->Get_par_i("Ctt")][0];
            double glbminCvv = cms->Get_v_Pars()[cms->Get_par_i("Cvv")][0];
            double glbminCglgl = cms->Get_v_Pars()[cms->Get_par_i("Cglgl")][0];
            double errhiCgg = Err_Bisection(cms, "Cgg", bestFitPars, y0_2, glbminCgg, 1);
            double errloCgg = Err_Bisection(cms, "Cgg", bestFitPars, y0_2, glbminCgg, 0);
            double errhiCbb = Err_Bisection(cms, "Cbb", bestFitPars, y0_2, glbminCbb, 1);
            double errloCbb = Err_Bisection(cms, "Cbb", bestFitPars, y0_2, glbminCbb, 0);
            double errhiCtt = Err_Bisection(cms, "Ctt", bestFitPars, y0_2, glbminCtt, 1);
            double errloCtt = Err_Bisection(cms, "Ctt", bestFitPars, y0_2, glbminCtt, 0);
            double errhiCglgl = Err_Bisection(cms, "Cglgl", bestFitPars, y0_2, glbminCglgl, 1);
            double errloCglgl = Err_Bisection(cms, "Cglgl", bestFitPars, y0_2, glbminCglgl, 0);
            double errhiCvv = Err_Bisection(cms, "Cvv", bestFitPars, y0_2, glbminCvv, 1);
            double errloCvv = Err_Bisection(cms, "Cvv", bestFitPars, y0_2, glbminCvv, 0);

            if(printInfo){
                cout<<" POI Cgg: "<<glbminCgg<<" +"<<errhiCgg<<" -"<<errloCgg<<endl;
                cout<<" POI Cvv: "<<glbminCvv<<" +"<<errhiCvv<<" -"<<errloCvv<<endl;
                cout<<" POI Cbb: "<<glbminCbb<<" +"<<errhiCbb<<" -"<<errloCbb<<endl;
                cout<<" POI Ctt: "<<glbminCtt<<" +"<<errhiCtt<<" -"<<errloCtt<<endl;
                cout<<" POI Cglgl: "<<glbminCglgl<<" +"<<errhiCglgl<<" -"<<errloCglgl<<endl;
            }
        }
        if(printInfo)watch.Print();
        delete pars;
        return 0;

    }

    if(bRunMinuitContour && idMH>0 ){

        TFile f(jobname+"_contourRvsM.root","RECREATE");

        // using mncont to extract points around contour = 1 or 2 sigma   
        myMinuit->SetErrorDef(TMath::ChisquareQuantile(0.68, NumDegreeOfFreedom));
        TGraph * g1sigma = (TGraph*)myMinuit->Contour(100, 0, idMH);
        if(g1sigma){
            g1sigma->SetName("ContourRvsM_1");
            for(int i=0; i<g1sigma->GetN(); i++){
                double r, m; g1sigma->GetPoint(i, r,m);
                g1sigma-> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
            }
            f.WriteTObject(g1sigma); 
        } else {cout<<" WARNING: minuit doesn't give 1 sigma contour (return 0 points) "<<endl;}


        myMinuit->SetErrorDef(TMath::ChisquareQuantile(0.95, NumDegreeOfFreedom));
        TGraph * g2sigma = (TGraph*)myMinuit->Contour(100, 0, idMH);
        if(g2sigma){
            g2sigma->SetName("ContourRvsM_2");
            for(int i=0; i<g2sigma->GetN(); i++){
                double r, m; g2sigma->GetPoint(i, r,m);
                g2sigma-> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
            }
            f.WriteTObject(g2sigma);
        } else {cout<<" WARNING: minuit doesn't give 2 sigma contour (return 0 points) "<<endl;}
        f.Close();
    }


    if(ErrEstAlgo == "Minos"){  // check covariance matrix and plot the ellipse if MH variable present in the floating parameter list 
        int npar = cms_global->Get_max_uncorrelation()+1;
        TMatrixDSym mat(npar); 
        if(idMH>0 && !DoNotRunGlobalFit){
            vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();
            double corr_ij = mat(0, idMH)/sqrt(mat(0,0)*mat(idMH,idMH));
            cout<<corr_ij<<endl;
            double x1, s1, x2, s2; 
            myMinuit->GetParameter(0, x1, s1);
            myMinuit->GetParameter(idMH, x2, s2);
            RooEllipse *contour= new RooEllipse("contour",x2,x1,s2,s1,corr_ij);
            contour->SetLineWidth(2) ;
            for(int i=0; i<contour->GetN(); i++){
                double r, m; contour->GetPoint(i, m,r);
                contour -> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
            }
            TCanvas c("ellipse","ellipse", 600, 600);
            contour->SetTitle("1 #sigma contour ;mass [GeV]; #mu");
            contour->Draw("AC");
            c.Print(jobname+"_ellipse_RM1sigma.pdf");
            c.Print(jobname+"_ellipse_RM1sigma.root");

            TFile f(jobname+"_ellipse_RM1sigma.root", "UPDATE");
            double m[1] = { x2 * v_paramsUnc[idMH][1]+ v_paramsUnc[idMH][3] };
            double merr[1] = { s2 * v_paramsUnc[idMH][1] };
            double mu[1] = { x1 };
            double muerr[1] = { s1 };
            TGraphErrors tge(1, m, mu, merr, muerr);
            tge.SetName("pointWithError");
            f.WriteTObject(&tge);

            double merrup[1] = { cms_global->POIs()[1].errUp};
            double merrdown[1] = { fabs(cms_global->POIs()[1].errDown)};
            double muerrup[1] = { cms_global->POIs()[0].errUp};
            double muerrdown[1] = { fabs(cms_global->POIs()[0].errDown)};
            TGraphAsymmErrors tgae(1, m, mu, merrdown, merrup, muerrdown, muerrup);
            tgae.SetName("pointWithAsymError");
            f.WriteTObject(&tgae);

            TString snumbers5  = "corr_"; snumbers5+=corr_ij;
            TObjString tsnumbers5(snumbers5);
            tsnumbers5.Write("corr");
            //f.WriteTObject(&tsnumbers5);
            TString snumbers6  = ""; snumbers6+=y0_2;
            TObjString tsnumbers6(snumbers6);
            tsnumbers6.Write("fAmin");
            f.Close();

        }

        if(printInfo){
            cout<<" 	68% CL : "<<endl;
            cout<<" POI: mu = "<<cms_global->POIs()[0].value<<" + "<<cms_global->POIs()[0].errUp<<" - "<<fabs(cms_global->POIs()[0].errDown)<<endl;
            for(int k=1; k<cms_global->POIs().size(); k++)cout<<" POI: "<<cms_global->POIs()[k].name<<"= "<<cms_global->POIs()[k].value<<" + "<<cms_global->POIs()[k].errUp<<" - "<<fabs(cms_global->POIs()[k].errDown)<<endl;
            watch.Print();
        }
    }


    if(makePlots){
        if(cms->GetPhysicsModel()==typeCvCfHiggs)SaveResults(jobname+"_maxllfit"+smakePlots, HiggsMass, 0, 0, 0, 0, 0, mu_hat_low, 0, mu_hat, mu_hat_up, 0);
        else if(cms->GetPhysicsModel()==typeC5Higgs)SaveResults(jobname+"_maxllfit"+smakePlots, HiggsMass, 0, 0, 0, 0, 0, mu_hat_low, 0, mu_hat, mu_hat_up, 0);
        else SaveResults(jobname+"_maxllfit"+smakePlots, HiggsMass, y0_2, 0, 0, 0, 0, mu_hat_low, 0, mu_hat, mu_hat_up, 0);
    }

    if(bDumpFitResults)cms->DumpFitResults(pars, jobname+"_fittedShape_floatMu");

    double *pars_maxll = new double[cms->Get_max_uncorrelation()+1]; 
    for(int i=0; i<=cms->Get_max_uncorrelation(); i++) pars_maxll[i]=pars[i];
    cms->Set_fittedParsInData_sb(pars_maxll);


    vector< vector<double> > v_Pars = cms_global->Get_v_Pars();	
    if(scanCvCf && cms_global->GetPhysicsModel()==typeCvCfHiggs){
        cms_global->SetNoErrorEstimation(1);
        TStopwatch watch2;
        watch2.Start();
        vector< std::pair<double, double> > vrm; vrm.clear(); 
        vector<double> vc; vc.clear();


        int _Cv_i = cms_global->Get_Cv_i();
        int _Cf_i = cms_global->Get_Cf_i();
        if(_Cv_i<0 or _Cf_i<0) {cout<<" no CV or CF parameter, exit"<<endl;  exit(1);}



        if(DoNotRunGlobalFit == false) { vrm.push_back(make_pair(pars[_Cv_i]*(v_Pars[_Cv_i][2]-v_Pars[_Cv_i][1])+v_Pars[_Cv_i][1],  pars[_Cf_i]*(v_Pars[_Cf_i][2]-v_Pars[_Cf_i][1])+v_Pars[_Cf_i][1]) );vc.push_back(y0_2-y0_2); }

        bool best = false;
        for(int i=0; i<vCv_toEval.size(); i++){
            double testr = vCv_toEval[i];
            for(int j=0; j<vCf_toEval.size(); j++){
                if(debug)watch.Start();
                double testm = vCf_toEval[j];
                if(testm==0 and testr==0) continue;
                _countPdfEvaluation = 0;
                cms_global->SetPOItoBeFixed("CV",testr);
                cms_global->SetPOItoBeFixed("CF",testm);


                double y0_1 =0; 
                if(bFixNuisancesAtBestFit){
                    bestFitPars[_Cv_i] =  (testr - v_Pars[_Cv_i][1])/(v_Pars[_Cv_i][2]-v_Pars[_Cv_i][1]);
                    bestFitPars[_Cf_i] =  (testr - v_Pars[_Cf_i][1])/(v_Pars[_Cf_i][2]-v_Pars[_Cf_i][1]);
                    y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);

                }else y0_1 =  MinuitFit(3, tmp, tmperr, 1.0, pars, best?true:false, debug, 0, best?bestFitPars:0);

                if(!bFixNuisancesAtBestFit){
                    for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                        myMinuit->GetParameter(i, tmp, tmpe);
                        bestFitPars[i]=tmp;
                    }
                    best = true;
                }


                if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                vrm.push_back(make_pair(testr, testm));
                vc.push_back(y0_1-y0_2);
                TString sj = jobname; sj+="_fittedShape_cv"; sj+=testr; sj+="_cf"; sj+=testm;
                if(bDumpFitResults)cms->DumpFitResults(pars, sj);
                if(debug){
                    watch.Stop();
                    watch.Print();
                }
            }
        }
        if(debug){
            printf("\n results of scanned Cv,Cf vs. q: \n");
            for(int i=0; i<vrm.size(); i++){
                printf("  [ Cv =%10.3f, Cf=%10.3f],  delta_q=%7.5f\n", vrm[i].first, vrm[i].second, vc[i]);
            }
        }
        if(bPlots and makePlots){
            DrawTH2D d2h(sbinningCV, sbinningCF, vrm, vc, "-2lnQ; CV; CF", (jobname+"_scanned_cvcf_vs_q_TH2"+smakePlots).Data(), pt);
            d2h.draw();
            d2h.save();
            DrawTGraph2D d2d(vrm, vc, "-2lnQ; CV; CF", (jobname+"_scanned_cvcf_vs_q"+smakePlots).Data(), pt);
            d2d.draw();
            d2d.save();
            //TFile f(jobname+"_scanned_cvcf_vs_q.root", "UPDATE");
            //d2d.getGraph()->SetName("Graph2D_cvcf_vs_q");					
            //f.WriteTObject(d2d.getGraph());
            //f.Close();
        }

        watch2.Stop();
        if(printInfo){
            cout<<"Scanning total takes: "<<endl;
            watch2.Print();
        }
        delete pars;
        return 0;
    }

    if(scanPOI && mainPOI!=""){
        cms_global->SetNoErrorEstimation(1);
        TStopwatch watch2;
        watch2.Start();
        vector<double> vc; vc.clear();
        vector<double> vr; vr.clear();

        int poi_i = -1;
        for(int i=0; i<cms_global->Get_v_uncname().size(); i++) {
            if(mainPOI==(TString)cms_global->Get_v_uncname()[i]) poi_i=i+1;
        }


        if(poi_i<0 ) {cout<<" no "<<mainPOI<<" parameter, exit"<<endl;  exit(1);}


        vector< vector<double> > v_Pars = cms_global->Get_v_Pars();	

        if(DoNotRunGlobalFit == false) { 
            vr.push_back(pars[poi_i]*(v_Pars[poi_i][2]-v_Pars[poi_i][1])+v_Pars[poi_i][1]);
            vc.push_back(y0_2-y0_2); 
        }

        bool best = false;
        for(int j=0; j<scanPOIrange.size(); j++){
            if(debug)watch.Start();
            double testm = scanPOIrange[j];
            _countPdfEvaluation = 0;
            cms_global->SetPOItoBeFixed(mainPOI,testm);


            double y0_1 =0; 
            if(bFixNuisancesAtBestFit){
                bestFitPars[poi_i] =  (testm - v_Pars[poi_i][1])/(v_Pars[poi_i][2]-v_Pars[poi_i][1]);
                y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);

            }else y0_1 =  MinuitFit(3, tmp, tmperr, 1.0, pars, best?true:false, debug, 0, best?bestFitPars:0);

            if(!bFixNuisancesAtBestFit){
                for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                    myMinuit->GetParameter(i, tmp, tmpe);
                    bestFitPars[i]=tmp;
                }
                best = true;
            }

            if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
            if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
            vr.push_back(testm);
            vc.push_back(y0_1-y0_2);
            if(debug){
                watch.Stop();
                watch.Print();
            }
        }
        if(debug){
            printf("\n results of scanned POI vs. q: \n");
            for(int i=0; i<vr.size(); i++){
                printf("  poi=%10.3f],  delta_q=%7.5f\n", vr[i], vc[i]);
            }
        }
        if(bPlots and makePlots){
            TString st1 = "-2lnQ;"; st1+=mainPOI; st1+=";q"; 
            TString st2 = jobname; st2+="_scanned_"; st2+=mainPOI; st2+="_vs_q_TH1"; st2+=smakePlots;
            DrawTH1D d1h(sbinningPOI, vr, vc, st1.Data(), st2.Data(), pt);
            d1h.draw();
            d1h.save();
            st2 = jobname; st2+="_scanned_"; st2+=mainPOI; st2+="_vs_q"; st2+=smakePlots;
            DrawEvolution2D d2d(vr, vc, st1.Data(), st2.Data(), pt);
            d2d.draw();
            d2d.save();
        }

        watch2.Stop();
        if(printInfo){
            cout<<"Scanning total takes: "<<endl;
            watch2.Print();
        }
        delete pars;
        return 0;
    }


    if(vCouplingsDef.size() and !scanMs) {delete pars;
        return 0;}

        double y0_3, mu_hat2, mu_hat_up2, mu_hat_low2;

        if(bAlsoExtract2SigmErrorBar){
            ErrorDef = TMath::ChisquareQuantile(0.95 , NumDegreeOfFreedom);// (confidenceLevel, ndf)
            _ErrorDef = ErrorDef;
            y0_3 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //102, 101, 21:  don't allow mu to be negative
            //double y0_3 =  MinuitFit(102, tmpr, tmperr, ErrorDef, pars, false, debug, success) ;  //102, 101, 21:  don't allow mu to be negative
            if(printInfo)cout<<y0_3<<" fitter u="<<pars[0]<<", from minos fit asymmetry 95% CL:  ["<<tmperr<<","<<tmpr<<"]"<<endl; // upperL, lowerL
            mu_hat2 = pars[0];
            mu_hat_up2 = tmpr;
            mu_hat_low2 = tmperr;

            {  // check covariance matrix and plot the ellipse if MH variable present in the floating parameter list 
                int npar = cms_global->Get_max_uncorrelation()+1;
                TMatrixDSym mat(npar); 
                myMinuit->mnemat( mat.GetMatrixArray(), npar);
                if(debug)mat.Print();

                if(idMH>0){
                    vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();
                    double corr_ij = mat(0, idMH)/sqrt(mat(0,0)*mat(idMH,idMH));
                    cout<<corr_ij<<endl;
                    double x1, s1, x2, s2; 
                    myMinuit->GetParameter(0, x1, s1);
                    myMinuit->GetParameter(idMH, x2, s2);
                    RooEllipse *contour= new RooEllipse("contour",x2,x1,s2,s1,corr_ij);
                    contour->SetLineWidth(2) ;
                    for(int i=0; i<contour->GetN(); i++){
                        double r, m; contour->GetPoint(i, m, r);
                        contour -> SetPoint(i, m*v_paramsUnc[idMH][1] + v_paramsUnc[idMH][3], r);
                    }
                    TCanvas c("ellipse","ellipse", 600, 600);
                    contour->SetTitle("1 #sigma contour ;mass [GeV]; #mu");
                    contour->Draw("AC");
                    c.Print(jobname+"_ellipse_RM2sigma.pdf");
                    c.Print(jobname+"_ellipse_RM2sigma.root");

                    TFile f(jobname+"_ellipse_RM2sigma.root", "UPDATE");
                    double m[1] = { x2 * v_paramsUnc[idMH][1]+ v_paramsUnc[idMH][3] };
                    double merr[1] = { s2 * v_paramsUnc[idMH][1] };
                    double mu[1] = { x1 };
                    double muerr[1] = { s1 };
                    TGraphErrors tge(1, m, mu, merr, muerr);
                    tge.SetName("pointWithError");
                    f.WriteTObject(&tge);

                    double merrup[1] = { cms_global->POIs()[1].errUp};
                    double merrdown[1] = { fabs(cms_global->POIs()[1].errDown)};
                    double muerrup[1] = { cms_global->POIs()[0].errUp};
                    double muerrdown[1] = { fabs(cms_global->POIs()[0].errDown)};
                    TGraphAsymmErrors tgae(1, m, mu, merrdown, merrup, muerrdown, muerrup);
                    tgae.SetName("pointWithAsymError");
                    f.WriteTObject(&tgae);

                    TString snumbers5  = "corr_"; snumbers5+=corr_ij;
                    TObjString tsnumbers5(snumbers5);
                    tsnumbers5.Write("corr");
                    f.Close();
                }
                if(printInfo){
                    cout<<" 	95% CL : "<<endl;
                    cout<<" POI: mu = "<<cms_global->POIs()[0].value<<" + "<<cms_global->POIs()[0].errUp<<" - "<<fabs(cms_global->POIs()[0].errDown)<<endl;
                    if(idMH>0)cout<<" POI: MH = "<<cms_global->POIs()[1].value<<" + "<<cms_global->POIs()[1].errUp<<" - "<<fabs(cms_global->POIs()[1].errDown)<<endl;
                    watch.Print();
                }
            }

        }


        ErrorDef = TMath::ChisquareQuantile(0.68 , NumDegreeOfFreedom);// (confidenceLevel, ndf)
        _ErrorDef = ErrorDef;
        bool b_global_fit = true;
        // FIXME also need to check the 2 sigma fit is fine 
        // if not, need to run brute-force approach
        if( (success[0]!=0 or tmpr==tmperr) and !DoNotRunGlobalFit){
            b_global_fit = false;
            double tmpr = 0;
            _initialRforFit = 0;
            y0_2 =  MinuitFit(21, tmpr, tmperr, 0, pars, false, debug) ;
            _initialRforFit = 1;
            if(debug)	cout<<y0_2<<" fitter u="<<tmpr<<" +/- "<<tmperr<<endl;
            mu_hat = pars[0];
        }


        if(b_global_fit==false){
            double x1 =0, x2 =1;
            double y0 = y0_2;  
            x1 = mu_hat; x2=2*mu_hat;   // has to be positive, otherwise mess up


            if(debug) cout<<" y0="<<y0<<" x1="<<x1<<" x2="<<x2<<endl;

            double CI = ErrorDef/2.;
            double precision = 0.001;
            int nsearched = 3;
            double y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
            if(fabs((y-y0)/2. - CI) > precision ){
                //first, got a number > 1.921,  otherwise increase it by a factor of 10 ...
                while( (y-y0)/2. < CI ){
                    x1 =  x2;
                    x2 *=10.;
                    y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                    if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                    nsearched++;
                }
                y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                if(debug) cout<< " r="<<(x1+x2)/2.<<" NLL = "<<y<<endl;
                while( fabs((y-y0)/2. - CI)/CI > precision ){
                    double  tmpx = (x1+x2)/2.;
                    if( (y-y0)/2. < CI ){
                        x1 = tmpx; 
                    }else{
                        x2 = tmpx;
                    }	
                    y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                    if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                    if(debug) printf(" r=%.6f,  NLL=%.10f\n", (x1+x2)/2., y);
                    nsearched++;
                }
                mu_hat_up= (x2+x1)/2.;
            }else{
                mu_hat_up = x2;
            }

        }
        // look for asymmetry errors on muhat
        if(b_global_fit==false){
            double x1 =0, x2 =1;
            double y0 = y0_2;  
            x1 = mu_hat; x2=0.5*mu_hat;  // has to be positive, otherwise mess up
            double precision = 0.001;

            double CI = ErrorDef/2.;
            double y =  MinuitFit(3, tmp, tmp, 0, pars, false, debug );
            if (fabs(y-y0_2)/2. <= CI+precision ){
                mu_hat_low = 0.;
            }else{
                if(debug) cout<<" y0="<<y0<<" x1="<<x1<<" x2="<<x2<<endl;
                int nsearched = 3;
                y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                if(fabs((y-y0)/2. - CI) > precision ){
                    //first, got a number > 1.921,  otherwise increase it by a factor of 10 ...
                    while( (y-y0)/2. < CI ){
                        x1 =  x2;
                        x2 /=10.;
                        y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                        if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                        nsearched++;
                    }
                    y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                    if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                    if(debug) cout<< " r="<<(x1+x2)/2.<<" NLL = "<<y<<endl;
                    while( fabs((y-y0)/2. - CI)/CI > precision ){
                        double  tmpx = (x1+x2)/2.;
                        if( (y-y0)/2. < CI ){
                            x1 = tmpx; 
                        }else{
                            x2 = tmpx;
                        }	
                        y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                        if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                        if(debug) printf(" r=%.6f,  NLL=%.10f\n", (x1+x2)/2., y);
                        nsearched++;
                    }
                    mu_hat_low= (x2+x1)/2.;
                }else{
                    mu_hat_low = x2;
                }
            }

        }

        if(printInfo){
            cout<<" muhat = "<<mu_hat<<" and 68% CL: ["<<mu_hat_low<<","<<mu_hat_up<<"]"<<endl;
            watch.Print();
        }
        if(bPlots and makePlots){
            TGraphAsymmErrors *ge = new TGraphAsymmErrors(1);
            ge->SetName("muhat");
            ge->SetPoint(0, mu_hat, y0_2);
            ge->SetPointError(0, mu_hat_low, mu_hat_up, 0, 0);
            TFile f((jobname+"_muhat.root"+smakePlots).Data(),"RECREATE");
            f.WriteTObject(ge);

            if(bAlsoExtract2SigmErrorBar){
                TGraphAsymmErrors *ge2 = new TGraphAsymmErrors(1);
                ge2->SetName("muhat_2sigma");
                ge2->SetPoint(0, mu_hat2, y0_3);
                ge2->SetPointError(0, mu_hat_low2, mu_hat_up2, 0, 0);
                f.WriteTObject(ge2);
            }
        }
        if(scanRs and scanMs){
            cms_global->SetNoErrorEstimation(1);
            TStopwatch watch2;
            watch2.Start();
            vector< std::pair<double, double> > vrm; vrm.clear(); 
            vector<double> vc; vc.clear();


            if(idMH<0) {cout<<" no MH parameter, exit"<<endl;  exit(1);}


            if(!DoNotRunGlobalFit){vrm.push_back(make_pair(pars[idMH]*v_paramsUnc[idMH][1]+v_paramsUnc[idMH][3], pars[0]) );vc.push_back(y0_2-y0_2);}

            bool best = false;

            for(int i=0; i<vR_toEval.size(); i++){
                double testr = vR_toEval[i];
                for(int j=0; j<vM_toEval.size(); j++){
                    if(debug)watch.Start();
                    double testm = vM_toEval[j];
                    _countPdfEvaluation = 0;
                    cms_global->SetPOItoBeFixed("MH",testm);
                    double y0_1 = 0;
                    if(bFixNuisancesAtBestFit){
                        bestFitPars[0] =  testr; 
                        bestFitPars[idMH] =  (testr - v_paramsUnc[idMH][3])/v_paramsUnc[idMH][1];
                        y0_1 =  MinuitFit(10, tmp, tmperr, testr, bestFitPars);

                    }else 
                        y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, best?true:false, debug, 0, best?bestFitPars:0);

                    if(!bFixNuisancesAtBestFit){
                        for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                            myMinuit->GetParameter(i, tmp, tmpe);
                            bestFitPars[i]=tmp;
                        }
                        best = true;
                    }


                    if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                    if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                    vrm.push_back(make_pair(testm, testr));
                    vc.push_back(y0_1-y0_2);
                    TString sj = jobname; sj+="_fittedShape_r"; sj+=testr; sj+="_m"; sj+=testm;
                    if(bDumpFitResults)cms->DumpFitResults(pars, sj);
                    if(debug){
                        watch.Stop();
                        watch.Print();
                    }
                }
            }
            if(debug){
                printf("\n results of scanned r,m vs. q: \n");
                for(int i=0; i<vrm.size(); i++){
                    printf("  [ m=%10.3f, r=%10.3f],  delta_q=%7.5f\n", vrm[i].first, vrm[i].second, vc[i]);
                }
            }
            if(bPlots and makePlots){
                DrawTH2D d2h(sbinningMH, sbinningMU, vrm, vc, "#delta(-2lnQ); mass [GeV]; #mu", (jobname+"_scanned_rm_vs_q_TH2"+smakePlots).Data(), pt);
                d2h.draw();
                d2h.save();

                DrawTGraph2D d2d(vrm, vc, "#delta(2lnQ); mass [GeV]; #mu", (jobname+"_scanned_rm_vs_q"+smakePlots).Data(), pt);
                d2d.draw();
                d2d.save();
                TFile f(jobname+"_scanned_rm_vs_q"+smakePlots+".root", "UPDATE");
                d2d.getGraph()->SetName("Graph2D_rm_vs_q");					
                f.WriteTObject(d2d.getGraph());
                f.Close();
            }

            watch2.Stop();
            if(printInfo){
                cout<<"Scanning total takes: "<<endl;
                watch2.Print();
            }
            delete pars;
            return 0;
        }
        else {
            cms_global->SetNoErrorEstimation(1);
            watch.Start();
            if(scanRs){
                vector<double> vr, vc; vr.clear(); vc.clear();
                if(!DoNotRunGlobalFit){
                    vr.push_back(pars[0]);vc.push_back(y0_2-y0_2);
                }
                bool best = false;
                for(int i=0; i<vR_toEval.size(); i++){
                    double testr = vR_toEval[i];
                    _countPdfEvaluation = 0;
                    double y0_1 = 0;
                    if(bFixNuisancesAtBestFit){
                        bestFitPars[0] =  testr; 
                        y0_1 =  MinuitFit(10, tmp, tmperr, testr, bestFitPars);

                    }else 
                        y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, best?true:false, debug, 0, best?bestFitPars:0);
                    if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                    if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                    vr.push_back(testr);
                    vc.push_back(y0_1-y0_2);
                    TString sj = jobname; sj+="_fittedShape_mu"; sj+=testr;
                    if(bDumpFitResults)cms->DumpFitResults(pars, sj);

                    if(!bFixNuisancesAtBestFit){
                        for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                            myMinuit->GetParameter(i, tmp, tmpe);
                            bestFitPars[i]=tmp;
                        }
                        best = true;
                    }

                }
                if( scanRs >= 2) {
                    if(mu_hat<vR_toEval[0]){
                        double dr = (vR_toEval[0] - mu_hat)/10.;
                        if(dr<0.2) dr=0.2;
                        for(double testr = mu_hat+0.2; testr<vR_toEval[0]; testr+=dr){
                            double y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, false, debug);
                            if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                            vr.push_back(testr);
                            vc.push_back(y0_1-y0_2);
                        }
                    }
                    if(mu_hat_up>vR_toEval.back()){
                        double dr = (mu_hat_up - vR_toEval.back())/10.;
                        if(dr<0.2) dr=0.2;
                        for(double testr = vR_toEval.back(); testr<=mu_hat_up+0.2; testr+=dr){
                            double y0_1 =  MinuitFit(3, tmp, tmperr, testr, pars, false, debug);
                            if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
                            vr.push_back(testr);
                            vc.push_back(y0_1-y0_2);
                        }
                    }
                }

                if(printInfo){
                    printf("\n results of scanned r vs. q: \n");
                    for(int i=0; i<vr.size(); i++){
                        printf("   r=%10.3f  delta_q=%7.5f\n", vr[i], vc[i]);
                    }
                }
                if(bPlots and makePlots){
                    DrawEvolution2D d2d(vr, vc, "; r ; q", (jobname+"_scanned_r_vs_q"+smakePlots).Data(), pt);
                    d2d.draw();
                    d2d.save();
                }
                watch.Stop();
                if(printInfo)watch.Print();
                delete pars;
                return 0;
            }
            if(scanMs){
                bool best = false;
                vector<double> vr, vc; vr.clear(); vc.clear();
                if(!DoNotRunGlobalFit){
                    vr.push_back(pars[0]);vc.push_back(y0_2-y0_2);
                }
                for(int i=0; i<vM_toEval.size(); i++){
                    double testr = vM_toEval[i];
                    _countPdfEvaluation = 0;
                    cms_global->SetPOItoBeFixed("MH",testr);

                    //double y0_1 =  MinuitFit(bConstrainMuPositive?102:202, tmpr, tmperr, 1./*common mu*/, pars, best?true:false, debug, success, best?bestFitPars:0) ;  //202, 201, 2:  allow mu to be negative
                    double y0_1 =0; 
                    if(bFixNuisancesAtBestFit){
                        bestFitPars[idMH] =  (testr - v_paramsUnc[idMH][3])/v_paramsUnc[idMH][1];
                        y0_1 =  MinuitFit(10, tmp, tmperr, 1., bestFitPars);

                    }else y0_1 =  MinuitFit(3, tmp, tmperr, 1., pars, best?true:false, debug, 0, best?bestFitPars:0);

                    if(debug) cout<<"_countPdfEvaluation="<<_countPdfEvaluation<<endl;
                    vr.push_back(testr);
                    vc.push_back(y0_1-y0_2);
                    TString sj = jobname; sj+="_fittedShape_m"; sj+=testr;
                    if(bDumpFitResults)cms->DumpFitResults(pars, sj);

                    if(!bFixNuisancesAtBestFit){
                        for(int i=0; i<cms->Get_max_uncorrelation()+1; i++){
                            myMinuit->GetParameter(i, tmp, tmpe);
                            bestFitPars[i]=tmp;
                        }
                        best = true;
                    }

                }
                if(printInfo){
                    printf("\n results of scanned m vs. q: \n");
                    for(int i=0; i<vr.size(); i++){
                        printf("   r=%10.3f  delta_q=%7.5f\n", vr[i], vc[i]);
                    }
                }
                if(bPlots and makePlots){
                    TString st1 = "-2lnQ;"; st1+="Higgs boson mass (GeV)"; st1+=";q"; 
                    TString st2 = jobname; st2+="_scanned_m_vs_q_TH1"; st2+=smakePlots;
                    DrawTH1D d1h(sbinningMH, vr, vc, st1.Data(), st2.Data(), pt);
                    d1h.draw();
                    d1h.save();
                    DrawEvolution2D d2d(vr, vc, "; m ; q", (jobname+"_scanned_m_vs_q"+smakePlots).Data(), pt);
                    d2d.draw();
                    d2d.save();
                }

                watch.Stop();
                if(printInfo)watch.Print();
                delete pars;
                return 0;
            }
        }

        return 0;

}


double runProfileLikelihoodApproximation(double neg2_llr){
    // change to use parabora approximation ....

    CountingModel * cms = cms_global;
    cms_global->SetDebug(debug);

    double r95;
    double tmp;
    double tmperr;
    double *pars = new double[cms->Get_max_uncorrelation()+1]; // nsys + r

    if(_inputNuisances == NULL) _inputNuisances = cms->Get_norminalPars();
    if(_startNuisances == NULL) _startNuisances = cms->Get_norminalPars();


    double ErrorDef = TMath::ChisquareQuantile(CL , 1);// (confidenceLevel, ndf)
    if(neg2_llr>=0) ErrorDef = neg2_llr;
    if(PLalgorithm == "Minos"){
        double upperL=1, lowerL=0; 
        //if(dataset=="asimov_b") upperL = 0.000001; // the starting mu value in the global fit
        // set it to be close to 0, because of it's asimov_b dataset, which may speed up fits during b-only toys
        // but it cause some problems ,  instabilities  FIXME
        double y0_2 =  MinuitFit(102, upperL, lowerL, ErrorDef, pars, false, debug) ;

        cout<<" best fitted r, upper bound = "<<upperL<<", nllmin="<<y0_2<<endl;
        if(upperL==lowerL){
            cout<<"WARNING: First Attempt Fails, try one more time with different set of starting values"<<endl;
            y0_2 =  MinuitFit(102, upperL, lowerL, ErrorDef, pars, true, debug, 0, pars) ;
            if(upperL==lowerL){
                cout<<"ERROR: need to be investigate --> two attempts fails "<<endl;
                cout<<" -----> trying Migrad"<<endl;
                PLalgorithm = "Migrad";
            }
        }
        r95 = upperL;
    }
    if(PLalgorithm == "Migrad"){
        double tmpr = 0;
        double y0_2 =  MinuitFit(2, tmpr, tmperr, 0, pars, false, debug) ;
        if(debug)	cout<<y0_2<<" fitter u="<<tmpr<<" +/- "<<tmperr<<endl;

        double y0_1 =  MinuitFit(3, tmp, tmperr, 0, pars, false, debug);
        if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;

        double x1 =0, x2 =1;
        double y0 = y0_1;  
        if ( (y0_1 > y0_2) && tmpr>0) {
            y0 = y0_2;
            x1 = tmpr; x2=2*tmpr;
        }


        //y0 = y0_2;  x1 = tmpr; x2=2*tmpr;
        //if(x2<1) x2 = 1.;
        if(debug) cout<<" y0="<<y0<<" x1="<<x1<<" x2="<<x2<<endl;


        //------------
        //If we have a background as a Gaussian distribution G(x|b) with mean b,
        //sigma sqrt(b), and observe x = b, then the upper limit on signal is
        //mu = 1.64*sqrt(b), which gives a 5% chance for P(obs < b).
        //
        //Therefore:
        //
        //1) muhat = 0, given the observation x=b
        //
        //2) -2 * ln (lambda(mu)) = -2 * ln( G(x|b+mu) / G(x|b+muhat) ) = 1.64^2
        //
        //If we use the 1.96-rule, it may *artificially* improve, but not cure,
        //the coverage for low statistics case, but would now give wrong coverage
        //in asymptotic.

        //http://en.wikipedia.org/wiki/Chi-square_distribution
        //double CI = 1.921;  // = 1.96**2/2 ,  probably for two sided 
        //double CI = 1.64*1.64/2.;
        double CI = ErrorDef/2.;
        /*
           if (oneside==2) CI= 1.921;  // = 1.96**2/2 , two sided 95% CL --->  one sided 97.5%
           else CI = 1.64*1.64/2.; // two sided 90% CL
         */
        double precision = 0.001;
        int nsearched = 3;
        //1.925 for 95% CL,   ....   
        //              assume the profile likelihood function to be an increasing function
        //              bisection search ...
        //              y1 always > y0

        double y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
        if(fabs((y-y0)/2. - CI) > precision ){
            //first, got a number > 1.921,  otherwise increase it by a factor of 10 ...
            while( (y-y0)/2. < CI ){
                x1 =  x2;
                x2 *=10.;
                y =  MinuitFit(3, tmp, tmp, x2, pars, false, debug );
                if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                nsearched++;
            }
            y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
            if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
            if(debug) cout<< " r="<<(x1+x2)/2.<<" NLL = "<<y<<endl;
            while( fabs((y-y0)/2. - CI)/CI > precision ){
                double  tmpx = (x1+x2)/2.;
                if( (y-y0)/2. < CI ){
                    x1 = tmpx; 
                }else{
                    x2 = tmpx;
                }	
                y = MinuitFit(3, tmp, tmp, (x1+x2)/2., pars, false, debug );
                if(debug) cout<<" NLL = "<<y<<" r="<<pars[0]<<endl;
                if(debug) printf(" r=%.6f,  NLL=%.10f\n", (x1+x2)/2., y);
                nsearched++;
            }
            r95 = (x2+x1)/2.;
        }else{
            r95 = x2;
        }


        if(debug)cout<<"r95 = "<<r95<<",  "<<nsearched<<" steps"<<endl;

    }

    cout<<"------------------------------------------------------------"<<endl;
    if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
    cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<r95<<endl;
    cout<<"------------------------------------------------------------"<<endl;

    SaveResults(jobname+"_plrObsLimit", HiggsMass, r95, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    if(bPlots){ // show a plot for  -log(Lambda(mu)) vs. mu 
        double x0 =  MinuitFit(21, tmp, tmperr, 0, pars, true, debug, 0, pars);
        //double x0 = y0;
        cout<<" x0 = "<<x0<<" optimized r="<<tmp<<endl;
        vector<double> vrxsec, vmlnq; 
        vrxsec.clear(); vmlnq.clear();
        for(double r=0; r<=2*r95; r+=r95/10.){
            vrxsec.push_back(r);
            double x3 = MinuitFit(3,tmp, tmp, r, pars, true, debug, 0, pars) ;
            cout<<"r= "<<r<<", x3 = "<<x3<<endl;
            double m2lnQ = x3 - x0;
            cout<<"Profiled:  r ="<<r<<"	m2lnQ="<<m2lnQ<<endl;
            vmlnq.push_back(m2lnQ/2.0);
            if(m2lnQ/2.>3) break;
        }
        TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
        DrawEvolution2D profiledL(vrxsec, vmlnq, ";r = #sigma/#sigma_{SM}; -log#lambda(r)", (jobname+"_profiledL").Data(), pt);
        profiledL.setLogY(0);
        profiledL.draw();
    }
    if(pars) delete pars;
    return r95;
}
bool runAsymptoticCLs(){
    // change to use parabora approximation ....
    CountingModel * cms = cms_global;
    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
    vdata_global = cms->Get_AsimovData(0);
    double r95;
    double tmp;
    double tmpr ;
    double tmperr;
    double *pars = new double[cms->Get_max_uncorrelation()+1]; // nsys + r

    _inputNuisances = cms->Get_norminalPars();
    _startNuisances = cms->Get_norminalPars();

    double starting_mu_for_globalFit = 1; // for data

    // for expected median limit,  asimov data set
    //double ErrorDef = TMath::NormQuantile(0.975)*TMath::NormQuantile(0.975); 
    double ErrorDef = TMath::ChisquareQuantile(CL , 1);// (confidenceLevel, ndf)

    int nsteps_in_asymptotic = 0;
    // running asymptotic limit from the hint 
    int success[1] = {0};
    int minuitSTRATEGY_old=minuitSTRATEGY;
    minuitSTRATEGY = minuitSTRATEGY_old;
    cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);

    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
    vdata_global = cms->Get_AsimovData(0);
    if(debug) cout<<" fitting asimovb with MinuitFit(2)"<<endl;
    double L_asimovb_glbmin =  MinuitFit(2, tmpr, tmperr, 0, pars, false, debug, success) ;
    minuitSTRATEGY_old=minuitSTRATEGY;
    while(success[0] and minuitSTRATEGY<2) {
        cout<<" failed global fit of asimovb with minuitSTRATEGY = "<<minuitSTRATEGY<<endl;
        //minuitSTRATEGY +=2 ; 
        cout<<" try minuitSTRATEGY = "<<minuitSTRATEGY<<endl;
        cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);
        //cms_global->FluctuatedNumbers();
        //L_asimovb_glbmin =  MinuitFit(201, tmpr, tmperr, 0.1, cms_global->Get_randomizedPars(), true, debug, success) ;

        cms_global->SetNoErrorEstimation(1);
        L_asimovb_glbmin =  MinuitFit(201, tmpr, tmperr, 0.1, pars, false, debug, success) ;
        cms_global->SetNoErrorEstimation(0);
    }
    minuitSTRATEGY = minuitSTRATEGY_old;
    cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);

    double mu =  1.;

    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
    vdata_global = cms->Get_AsimovData(0);
    if(debug) cout<<" fitting asimovb with MinuitFit(3)"<<endl;
    double L_asimovb_condmin =  MinuitFit(3, tmp, tmperr, mu, pars, false, debug);

    nsteps_in_asymptotic +=4;
    if(L_asimovb_condmin - L_asimovb_glbmin <0 ) {
        cout<<"2 L_asimovb_condmin = "<<L_asimovb_condmin<<"  < "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
        exit(1);
    }
    double L_data_condmin = L_asimovb_condmin; 
    double L_data_glbmin = L_asimovb_glbmin;
    double tmpcls = (  1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) ) / ( ROOT::Math::normal_cdf(
                sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin)
                ) ); 


    if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
    cout<<"Observed Asymptotic CLs  "<<tmpcls<<endl;
    cout<<"------------------------------------------------------------"<<endl;

    SaveResults(jobname+"_asymptoticObsCLs", HiggsMass, tmpcls, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    if(pars) delete pars;
    return true;
}
bool runAsymptoticLimits(){
    // change to use parabora approximation ....
    CountingModel * cms = cms_global;
    vdata_global=cms->Get_v_data();
    cms_global->SetDebug(debug);

    if(dataset=="asimov_b"){
        cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
        vdata_global = cms->Get_AsimovData(0);
    }
    double r95;
    double tmp;
    double tmpr ;
    double tmperr;
    double *pars = new double[cms->Get_max_uncorrelation()+1]; // nsys + r

    _inputNuisances = cms->Get_norminalPars();
    _startNuisances = cms->Get_norminalPars();

    double starting_mu_for_globalFit = 1; // for data   1. 

    // for expected median limit,  asimov data set
    //double ErrorDef = TMath::NormQuantile(0.975)*TMath::NormQuantile(0.975); 
    double ErrorDef = TMath::ChisquareQuantile(CL , 1);// (confidenceLevel, ndf)
    if(1){ // first get hint of the limit from  PLR method, wrong asymptotic formula
        //cout<<"  DELETEME :   _inputNuisances[2] = "<<_inputNuisances[2]<<endl;
        if(dataset == "asimov_b" and bConstructAsimovbFromNominal==false){
            _inputNuisances = cms->Get_fittedParsInData_b();
            _startNuisances = cms->Get_fittedParsInData_b();
        }
        //cout<<"  DELETEME :   _inputNuisances[2] = "<<_inputNuisances[2]<<endl;
        r95 = runProfileLikelihoodApproximation(ErrorDef);
        if(dataset=="asimov_b"){
            double N = 0;
            double m2s, m1s, p1s, p2s;
            if(scanRAroundBand=="all" || scanRAroundBand=="0"){
                cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
                vdata_global = cms->Get_AsimovData(0);
                preHints_median = r95;
            }
            if(newExpectedAsymBand){
                if(scanRAroundBand=="all" || scanRAroundBand=="-1"){
                    N = -1;
                    ErrorDef = TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ; 
                    ErrorDef *= ErrorDef;	
                    if(debug)cout<<" ErrorDef = "<<ErrorDef<<endl;
                    m1s = runProfileLikelihoodApproximation(ErrorDef);
                    preHints_m1s=m1s;
                }
                if(scanRAroundBand=="all" || scanRAroundBand=="1"){
                    N = 1;
                    ErrorDef = TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ; 
                    ErrorDef *= ErrorDef;	
                    if(debug)cout<<" ErrorDef = "<<ErrorDef<<endl;
                    p1s = runProfileLikelihoodApproximation(ErrorDef);
                    preHints_p1s=p1s;
                }
                if(scanRAroundBand=="all" || scanRAroundBand=="-2"){
                    N = -2;
                    ErrorDef = TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ; 
                    ErrorDef *= ErrorDef;	
                    if(debug)cout<<" ErrorDef = "<<ErrorDef<<endl;
                    m2s = runProfileLikelihoodApproximation(ErrorDef);
                    preHints_m2s=m2s;
                }
                if(scanRAroundBand=="all" || scanRAroundBand=="1"){
                    N = 2;
                    ErrorDef = TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ; 
                    ErrorDef *= ErrorDef;	
                    if(debug)cout<<" ErrorDef = "<<ErrorDef<<endl;
                    p2s = runProfileLikelihoodApproximation(ErrorDef);
                    preHints_p2s=p2s;
                }
            }else{


                // old Asymptoic ******* 
                // http://cms.cern.ch/iCMS/jsp/analysis/admin/analysismanagement.jsp?ancode=HIG-11-011
                // AN 11-298
                // equation (38) : mu_median = sigma * normal_quantile(1-0.5*alpha) 
                double sigma = preHints_median / TMath::NormQuantile(1-0.5*(1-CL));

                cout<<" median = "<<r95<<", sigma = "<<sigma<<endl;

                // equation(37)  : mu_N = sigma*(normal_quantile(1 - quantile*alpha, 1.0) + normal_quantile(quantile))
                // equation(37)  : mu_N = sigma*(normal_quantile(1 - normal_cdf(N)*alpha) + N ) 
                N = -2;
                m2s = sigma *  ( TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ) ; 
                N = -1;
                m1s = sigma *  ( TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ) ; 
                N = 1 ;
                p1s = sigma *  ( TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ) ; 
                N = 2;
                p2s = sigma *  ( TMath::NormQuantile(1-(1-CL)*ROOT::Math::normal_cdf(N))+N ) ; 

                preHints_p2s = p2s; preHints_p1s=p1s; preHints_m2s=m2s; preHints_m1s=m1s;
            }
            cout<<"BandsFromAsymptotic R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
            if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
            printf("BANDS %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", m2s, m1s, 0.00, p1s, p2s, r95);
            cout<<"------------------------------------------------------------"<<endl;
            SaveResults(jobname+"Asymptotic"+"_limitbands", HiggsMass, preHints_median, 0, 0, 0, preHints_m2s, preHints_m1s, preHints_median, 0, preHints_p1s, preHints_p2s);
        }
    }
    int nsteps_in_asymptotic = 0;
    // running asymptotic limit from the hint 
    if(dataset=="data_obs"){
        vdata_global = cms->Get_v_data();
        cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
        if(debug) cout<<" fitting data with MinuitFit(2)"<<endl;
        int success[1] = {0};

        _inputNuisances = cms->Get_norminalPars();
        double L_data_glbmin =  MinuitFit(2, tmpr, tmperr, starting_mu_for_globalFit, pars, false, debug, success) ;
        int minuitSTRATEGY_old=minuitSTRATEGY;
        while(success[0] and minuitSTRATEGY<2) {
            cout<<" failed global fit of data with minuitSTRATEGY = "<<minuitSTRATEGY<<endl;
            minuitSTRATEGY +=2 ; 
            cout<<" try minuitSTRATEGY = "<<minuitSTRATEGY<<endl;
            cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);
            L_data_glbmin =  MinuitFit(2, tmpr, tmperr, starting_mu_for_globalFit+0.1, pars, false, debug, success) ;
        }
        if(tmpr<0){
            L_data_glbmin =  MinuitFit(3, tmp, tmperr, 0, pars, false, debug);
            tmpr=0;
        }
        double rfitted_d = tmpr;
        minuitSTRATEGY = minuitSTRATEGY_old;
        cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);

        cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
        vdata_global = cms->Get_AsimovData(0);
        if(!bConstructAsimovbFromNominal)_inputNuisances = cms->Get_fittedParsInData_b();
        if(debug) cout<<" fitting asimovb with MinuitFit(2)"<<endl;
        double L_asimovb_glbmin =  MinuitFit(2, tmpr, tmperr, 0, pars, false, debug, success) ;
        cout<<" fitted r = "<<tmpr<<" nllmin = "<<L_asimovb_glbmin<<endl;
        minuitSTRATEGY_old=minuitSTRATEGY;

        while(success[0] and minuitSTRATEGY<2) {
            cout<<" failed global fit of asimovb with minuitSTRATEGY = "<<minuitSTRATEGY<<endl;
            //minuitSTRATEGY +=2 ; 
            cout<<" try minuitSTRATEGY = "<<minuitSTRATEGY<<endl;
            cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);
            //cms_global->FluctuatedNumbers();
            //L_asimovb_glbmin =  MinuitFit(201, tmpr, tmperr, 0.1, cms_global->Get_randomizedPars(), true, debug, success) ;
            cms_global->SetNoErrorEstimation(1);
            L_asimovb_glbmin =  MinuitFit(201, tmpr, tmperr, 0.1, pars, false, debug, success) ;
            cms_global->SetNoErrorEstimation(0);
            cout<<" fitted r = "<<tmpr<<" nllmin = "<<L_asimovb_glbmin<<endl;
        }
        if(tmpr<0) {
            L_asimovb_glbmin = MinuitFit(3, tmp, tmperr, 0,  pars, false, debug);
            tmpr = 0 ;
        }
        double rfitted_a = tmpr;
        if(1)cout<<" data_global min ="<<L_data_glbmin<<", best fitted r = "<<rfitted_d<<endl;
        if(1)cout<<" asimov_global min ="<<L_asimovb_glbmin<<", best fitted r = "<<rfitted_a<<endl;

        minuitSTRATEGY = minuitSTRATEGY_old;
        cms_global->Set_minuitSTRATEGY(minuitSTRATEGY);

        double mu =  r95;

        vdata_global = cms->Get_v_data();
        cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
        _inputNuisances = cms->Get_norminalPars();
        if(debug) cout<<" fitting data with MinuitFit(3)"<<endl;
        double L_data_condmin =  MinuitFit(3, tmp, tmperr, mu, pars, false, debug);

        cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
        vdata_global = cms->Get_AsimovData(0);
        if(!bConstructAsimovbFromNominal)_inputNuisances = cms->Get_fittedParsInData_b();
        if(debug) cout<<" fitting asimovb with MinuitFit(3)"<<endl;
        double L_asimovb_condmin =  MinuitFit(3, tmp, tmperr, mu, pars, false, debug);

        nsteps_in_asymptotic +=4;
        if(L_data_condmin - L_data_glbmin <0 ) {
            cout<<"1 L_data_condmin = "<<L_data_condmin<<"  < "<<"L_data_glbmin="<<L_data_glbmin<<", exit"<<endl;
            exit(1);
        }
        if(L_asimovb_condmin - L_asimovb_glbmin <0 ) {
            cout<<"2 L_asimovb_condmin = "<<L_asimovb_condmin<<"  < "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
            exit(1);
        }
        double tmpcls = (  1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) ) / ( ROOT::Math::normal_cdf(
                    sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin)
                    ) ); 

        double precision = 0.001;
        double mu_up = mu, mu_down=mu;

        if(debug or tmpcls==0)cout<<"tmpcls = "<<tmpcls<<endl;
        if(debug or tmpcls==0)cout<<"mu = "<<mu<<" mu_up="<<mu_up<<" mu_down="<<mu_down<<endl;
        if(tmpcls==0) {
            cout<<" L_data_condmin = "<<L_data_condmin<<", L_data_glbmin="<<L_data_glbmin<<"; L_asimovb_condmin="<<L_asimovb_condmin<<", L_asimovb_glbmin="<<L_asimovb_glbmin<<endl;
            cout<<"ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) = "<<ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) <<endl;
        }
        if(singlePoint){
            vdata_global = cms->Get_v_data();
            cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
            _inputNuisances = cms->Get_norminalPars();
            L_data_condmin =  MinuitFit(3, tmp, tmperr, testR, pars, false, debug);
            vdata_global = cms->Get_AsimovData(0);
            cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
            if(!bConstructAsimovbFromNominal)_inputNuisances = cms->Get_fittedParsInData_b();
            L_asimovb_condmin =  MinuitFit(3, tmp, tmperr, testR, pars, false, debug);

            if(L_data_condmin - L_data_glbmin <0 ) {
                cout<<"3 L_data_condmin = "<<L_data_condmin<<"  < "<<"L_data_glbmin="<<L_data_glbmin<<", exit"<<endl;
                exit(1);
            }
            if(L_asimovb_condmin - L_asimovb_glbmin <0 ) {
                cout<<"4 L_asimovb_condmin = "<<L_asimovb_condmin<<"  < "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
                exit(1);
            }
            tmpcls = (  1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) ) / ( ROOT::Math::normal_cdf(
                        sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin)
                        ) ); 

            cout<<" L_asimovb_condmin = "<<L_asimovb_condmin<<" , "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
            cout<<" qmu = "<<L_data_condmin - L_data_glbmin<<" qasimov = "<<L_asimovb_condmin - L_asimovb_glbmin<<endl;
            cout<<" RESULTS:  mu  = "<<testR;
            cout<<" clsb = "<<1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) );
            cout<<" clb = "<<ROOT::Math::normal_cdf(sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin));
            cout<<" cls = "<<tmpcls<<endl;
            return 0;
        }
        if( fabs(tmpcls - (1-CL)) > precision ) {
            if( tmpcls > 1-CL ){
                while ( tmpcls > (1-CL) ) {
                    mu_up *= 2.; 
                    vdata_global = cms->Get_v_data();
                    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
                    _inputNuisances = cms->Get_norminalPars();
                    L_data_condmin =  MinuitFit(3, tmp, tmperr, mu_up, pars, false, debug);
                    vdata_global = cms->Get_AsimovData(0);
                    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
                    if(!bConstructAsimovbFromNominal)_inputNuisances = cms->Get_fittedParsInData_b();
                    L_asimovb_condmin =  MinuitFit(3, tmp, tmperr, mu_up, pars, false, debug);

                    if(L_data_condmin - L_data_glbmin <0 ) {
                        cout<<"3 L_data_condmin = "<<L_data_condmin<<"  < "<<"L_data_glbmin="<<L_data_glbmin<<", exit"<<endl;
                        exit(1);
                    }
                    if(L_asimovb_condmin - L_asimovb_glbmin <0 ) {
                        cout<<"4 L_asimovb_condmin = "<<L_asimovb_condmin<<"  < "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
                        exit(1);
                    }
                    tmpcls = (  1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) ) / ( ROOT::Math::normal_cdf(
                                sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin)
                                ) ); 
                    nsteps_in_asymptotic +=2;

                    if(debug)cout<<"tmpcls = "<<tmpcls<<endl;
                    if(debug)cout<<"mu = "<<mu<<" mu_up="<<mu_up<<" mu_down="<<mu_down<<endl;
                }
            }else{
                while ( tmpcls < (1-CL) ) {
                    mu_down /= 2.; 
                    vdata_global = cms->Get_v_data();
                    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
                    _inputNuisances = cms->Get_norminalPars();
                    L_data_condmin =  MinuitFit(3, tmp, tmperr, mu_down, pars, false, debug);
                    vdata_global = cms->Get_AsimovData(0);
                    cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
                    if(!bConstructAsimovbFromNominal)_inputNuisances = cms->Get_fittedParsInData_b();
                    L_asimovb_condmin =  MinuitFit(3, tmp, tmperr, mu_down, pars, false, debug);

                    if(L_data_condmin - L_data_glbmin <0 ) {
                        cout<<"5 L_data_condmin = "<<L_data_condmin<<"  < "<<"L_data_glbmin="<<L_data_glbmin<<", exit"<<endl;
                        exit(1);
                    }
                    if(L_asimovb_condmin - L_asimovb_glbmin <0 ) {
                        cout<<"6 L_asimovb_condmin = "<<L_asimovb_condmin<<"  < "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
                        exit(1);
                    }
                    tmpcls = (  1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) ) / ( ROOT::Math::normal_cdf(
                                sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin)
                                ) ); 

                    nsteps_in_asymptotic +=2;
                    if(debug)cout<<"tmpcls = "<<tmpcls<<endl;
                    if(debug)cout<<"mu = "<<mu<<" mu_up="<<mu_up<<" mu_down="<<mu_down<<endl;
                }
            }
            int iteration = 0;
            while ( fabs(tmpcls - (1-CL)) > precision and iteration<=20 ) {
                mu = (mu_up + mu_down)/2.;
                if(mu==mu_up and mu==mu_down) {cout<<" WARNING **** infinite loop , mu=mu_up=mu_down="<<mu<<endl; break;}

                vdata_global = cms->Get_v_data();
                cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset());
                _inputNuisances = cms->Get_norminalPars();
                L_data_condmin =  MinuitFit(3, tmp, tmperr, mu, pars, false, debug);
                vdata_global = cms->Get_AsimovData(0);
                cms->SetTmpDataForUnbinned(cms->Get_v_pdfs_roodataset_asimovb());
                if(!bConstructAsimovbFromNominal)_inputNuisances = cms->Get_fittedParsInData_b();
                L_asimovb_condmin =  MinuitFit(3, tmp, tmperr, mu, pars, false, debug);

                if(L_data_condmin - L_data_glbmin <0 ) {
                    cout<<"5 L_data_condmin = "<<L_data_condmin<<"  < "<<"L_data_glbmin="<<L_data_glbmin<<", exit"<<endl;
                    exit(1);
                }
                if(L_asimovb_condmin - L_asimovb_glbmin <0 ) {
                    cout<<"6 L_asimovb_condmin = "<<L_asimovb_condmin<<"  < "<<"L_asimovb_glbmin="<<L_asimovb_glbmin<<", exit"<<endl;
                    exit(1);
                }
                tmpcls = (  1 - ROOT::Math::normal_cdf( sqrt(L_data_condmin - L_data_glbmin) ) ) / ( ROOT::Math::normal_cdf(
                            sqrt(L_asimovb_condmin - L_asimovb_glbmin) - sqrt(L_data_condmin - L_data_glbmin)
                            ) ); 

                if(tmpcls > 1-CL) {
                    mu_down = mu;
                }else{
                    mu_up = mu;
                }

                if(debug)cout<<"tmpcls = "<<tmpcls<<endl;
                if(debug)cout<<"mu = "<<mu<<" mu_up="<<mu_up<<" mu_down="<<mu_down<<endl;
                nsteps_in_asymptotic +=2;
                iteration ++;
            }
        }
        if(debug)cout<<"tmpcls = "<<tmpcls<<endl;
        if(debug)cout<<"mu = "<<mu<<" mu_up="<<mu_up<<" mu_down="<<mu_down<<endl;
        r95 = mu;
        if(debug) cout<<"Addtional number of fits in asymptotic =  "<<nsteps_in_asymptotic<<endl;

        preHints_obs = r95;
    }


    cout<<"------------------------------------------------------------"<<endl;
    if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
    cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<r95<<endl;
    cout<<"------------------------------------------------------------"<<endl;

    SaveResults(jobname+"_asymptoticObsLimit", HiggsMass, r95, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    if(pars) delete pars;
    return true;
}

bool runBayesian(){
    CountingModel*cms=cms_global;
    BayesianBase bys(cms, 1-CL, 1.e-3);
    bys.SetDebug(debug);
    bys.SetNumToys(toysBayesian);
    bys.SetPreToys(toysPreBayesian);
    bys.SetCrossSectionPrior(prior);
    double rtmp;
    if(bCalcObservedLimit){
        double *parsFitted, *parsErrLow, *parsErrUp;
        parsFitted = new double[cms->Get_max_uncorrelation()+1];
        parsErrUp= new double[cms->Get_max_uncorrelation()+1];
        parsErrLow = new double[cms->Get_max_uncorrelation()+1];
        rtmp = bys.Limit(1-CL, -99999., true, parsFitted, parsErrLow, parsErrUp, flatPriorCropNsigma, bRunOnlyWithBestFittedNuisances_bayesian, inputMu_bayesian);

        if(nTries>1)cout<<"try "<<1<<": R at 95\% CL = "<<rtmp<<endl;

        vector<double> rtries;
        double hint = rtmp;
        rtries.push_back(rtmp);
        double avgR=0, errR=0;
        for(int i=1; i<nTries; i++){
            rtmp = bys.Limit(1-CL, hint, false, parsFitted, parsErrLow, parsErrUp, flatPriorCropNsigma);
            cout<<"try "<<i+1<<": R at 95\% CL = "<<rtmp<<endl;
            if(rtmp<=0) cout<<"  - -- ---- r is meaningless, skipped "<<endl;
            else{
                rtries.push_back(rtmp);
                avgR=0;
                for(int j=0; j<rtries.size(); j++) avgR+=rtries[j]; avgR/=float(rtries.size());	
                hint = avgR;
                errR=0;
                for(int i=0; i<rtries.size(); i++) errR+= (rtries[i]-avgR)*(rtries[i]-avgR); errR = sqrt(errR)/(float)(rtries.size());
                cout<<"try current arg = "<<avgR<<" +/- "<<errR<<endl;
            }
        }
        if(nTries<=1) avgR=rtmp;
        errR=0;
        for(int i=0; i<rtries.size(); i++) errR+= (rtries[i]-avgR)*(rtries[i]-avgR); errR = sqrt(errR)/(float)(rtries.size());

        std::sort(rtries.begin(), rtries.end());
        int idmedian = int(rtries.size()*0.5); if(idmedian>=rtries.size()) idmedian=rtries.size()-1;
        if(idmedian>0)cout<<"  median value of all tries = "<<rtries[idmedian]<<endl;


        cout<<"------------------------------------------------------------"<<endl;
        if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
        cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<avgR;
        if(nTries>1)cout<<" +/- "<<errR<<endl;
        else cout<<endl;
        cout<<"------------------------------------------------------------"<<endl;

        preHints_obs = avgR;

        SaveResults(jobname+"_bysObsLimit", HiggsMass, avgR, errR, 0, 0, 0, 0, 0, 0, 0, 0);

        if(bPlots)	{
            bys.PosteriorPdf();
            //start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimit ppdf: "<< (cur_time - start_time)/1000. << " millisec\n";
            TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
            DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", (jobname+"_postpdf").Data(), pt);	
            pdfr.setDrawPosteriorPdf(rtmp);
            pdfr.draw();
            pdfr.getGraph()->SetMarkerSize(0.01);
            pdfr.save();
        }
    }

    if(doExpectation && sFileLimitsDistribution==""){
        difflimits.clear();
        cms->SetSignalScaleFactor(1.);
        LimitBands lb(&bys, cms);	
        //	lb.SetTossPseudoDataConvention(tossPseudoDataConvention);
        //	bayesian toys generation should follow LEP convention
        lb.SetDebug(debug);
        int noutcomes = toys;
        lb.BysLimitBands(1-CL, noutcomes, toysBayesian);
        rmean=lb.GetBysLimitMean();
        difflimits = lb.GetDifferentialLimitsBys();
        saveExpectation();
    }
}

bool saveExpectation(){
    // plot the distribution of limits of noutcomes
    // and the cummulative pdf
    TString ts=jobname; ts+=method; ts+="_limits";

    if(sFileLimitsDistribution!=""){
        TTree * tree = (TTree*)GetTObject(sFileLimitsDistribution.Data(), "T");
        difflimits = GetVectorFrom(tree, "brT");
        rmean = 0;
        for(int i=0; i<difflimits.size(); i++)rmean+=difflimits[i];
        if(difflimits.size()==0) rmean=0; else rmean/=float(difflimits.size());
    }else 	FillTree(ts, difflimits);


    vector<double> all_calculated_R95s;
    vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
    SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

    // show bands without interpolation,  basically using step function
    double GreenBandLow = (1- 0.683)/2.; //1 sigma
    double GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
    double YellowBandLow = (1- 0.955)/2.; //2 sigma
    double YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
    double rmedian2, rm1s2, rm2s2, rp1s2, rp2s2;	
    bool brmedian2=false, brm1s2=false, brm2s2=false, brp1s2=false, brp2s2=false;	
    for(int i=0; i<cummulativeProbabilities.size(); i++){
        if(cummulativeProbabilities[i]>=GreenBandLow && !brm1s2) { rm1s2 = all_calculated_R95s[i]; brm1s2 = true; } 
        if(cummulativeProbabilities[i]>=GreenBandHigh && !brp1s2) { rp1s2 = all_calculated_R95s[i]; brp1s2 = true; } 
        if(cummulativeProbabilities[i]>=YellowBandLow && !brm2s2) { rm2s2 = all_calculated_R95s[i]; brm2s2 = true; } 
        if(cummulativeProbabilities[i]>=YellowBandHigh && !brp2s2) { rp2s2 = all_calculated_R95s[i]; brp2s2 = true; } 
        if(cummulativeProbabilities[i]>=0.5 && !brmedian2) { rmedian2 = all_calculated_R95s[i]; brmedian2 = true; } 
    }
    cout<<"------------NO INTERPOLATION--------------------------------"<<endl;
    cout<<"BandsNoInterpolation R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
    if(HiggsMass>0)cout<<"MassPoint "<<HiggsMass<<" , ";
    printf("BANDS %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s2, rm1s2, rmean, rp1s2, rp2s2, rmedian2);
    cout<<"------------------------------------------------------------"<<endl;

    SaveResults(jobname+method+"_limitbands", HiggsMass, preHints_obs, 0, 0, 0, rm2s2, rm1s2, rmedian2, rmean, rp1s2, rp2s2);
    if(bPlots){
        TCanvas *c=new TCanvas("cme","cme");
        c->SetLogy(1);
        TH1F *h=new TH1F("h",";r=#frac{#sigma_{95%CL}}{#sigma_{SM}}; entries", 200, 0, 15);	
        for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
        h->Draw();
        Save(c, (jobname+method+"_differential_limits").Data());

        TPaveText *pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
        string ssave= (jobname+method+"_cummulativeR").Data();
        string stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
        PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
                rm1s2, rp1s2, rm2s2, rp2s2, ssave, stitle, pt);
        plotRvsP.draw();
        plotRvsP.getGraph()->SetMarkerSize(1);
        plotRvsP.save();
    }

    preHints_median = rmedian2; preHints_p2s = rp2s2; preHints_p1s = rp2s2; preHints_m2s = rm2s2; preHints_m1s = rm1s2;
    return true;

}
void processParameters(int argc, const char* argv[])
{
    vector<TString> allargs;
    for(int i=1; i<argc; i++){
        allargs.push_back(TString(argv[i]));
    }

    for(int i=0; i<allargs.size(); i++){
        if(allargs[i]=="--help" or allargs[i]=="-h" ) PrintHelpMessage();
        //cout<<allargs[i]<<endl;
        if(allargs[i].BeginsWith("-") and !allargs[i].IsFloat()){
            TString key = allargs[i];
            vector<TString> values;
            for(int j=i+1; j<allargs.size(); j++){
                if(allargs[j].BeginsWith("-") and !allargs[j].IsFloat()) break;
                values.push_back(allargs[j]);
            }
            options.insert(TStrPair(key, values));
        }
    }

    TStrMap::iterator p;

    cout<<"Your arguments: "<<endl;
    for(p=options.begin(); p!=options.end();++p){
        cout<<p->first<<" ";
        for(int j=0; j<p->second.size(); j++){
            cout<<p->second[j]<<" ";
        }
        cout<<endl;
    }

    vector<TString> tmpv;
    tmpv = options["--readLimitsFromFile"]; 
    if( tmpv.size()!=1 ) { }
    else sFileLimitsDistribution = tmpv[0];

    tmpv = options["--calcPValueFromFile"]; 
    if( tmpv.size()!=1 ) { }
    else sFileBonlyM2lnQ = tmpv[0];

    // 
    tmpv = options["--makeCLsBands"];
    if( tmpv.size()!=1 ) { makeCLsBands= 0; }
    else {
        makeCLsBands = tmpv[0].Atoi();
    }

    vector<TString> tmpcards = options["-d"];
    if(tmpcards.size()==0) tmpcards = options["--datacards"];
    cout<<endl<<"You are trying to combine following "<<tmpcards.size()<<" cards:"<<endl;
    if(tmpcards.size()==0) { 
        if(sFileLimitsDistribution=="" && makeCLsBands<2 && sFileBonlyM2lnQ==""){
            cout<<"*** No data card specified, please use option \"-d or --datacards\""<<endl; exit(0); 
        }
    }else{
        for(int i=0; i<tmpcards.size(); i++){
            FileStat_t buf;
            if(gSystem->GetPathInfo(tmpcards[i], buf) ==1 ) { cout<<tmpcards[i] << " not found, skipped"<<endl;}
            else {
                cout<<tmpcards[i]<<endl;
                datacards.push_back(tmpcards[i]);
            }
        }
        cout<<" .... valid data cards = "<<datacards.size()<<endl;
        if(datacards.size()<=0){ cout<< " please provide valid data cards "<<endl; exit(0); }
    }

    if(isWordInMap("--bForceSymmetryError", options)) bForceSymmetryError = true;
    if(isWordInMap("--bMultiSigProcShareSamePDF", options)) bMultiSigProcShareSamePDF = true;

    tmpv= options["-L"];if(tmpv.size()==0)tmpv=options["--LoadLibraries"];
    librariesToBeLoaded = tmpv;

    tmpv= options["-m"];if(tmpv.size()==0)tmpv=options["--mass"];
    if( tmpv.size()!=1 ) { HiggsMass = -1; }
    else HiggsMass = tmpv[0].Atof();

    if(isWordInMap("--plot", options)) bPlots = 1;

    // limit or significance
    tmpv = options["--significance"]; 
    if( tmpv.size()!=1 ) { calcsignificance = 0; }
    else calcsignificance = tmpv[0].Atoi();

    tmpv = options["--plot"]; 
    if( tmpv.size()!=1 ) { }
    else bPlots = tmpv[0].Atoi();

    tmpv = options["-M"]; if(tmpv.size()!=1) tmpv = options["--method"];
    if( tmpv.size()!=1) { 
        if(sFileLimitsDistribution=="" && makeCLsBands<2 && sFileBonlyM2lnQ==""){
            cout<<"ERROR No method specified, please use option \"-M or --method\" "<<endl; exit(0);
        }
    }else {
        method = tmpv[0];
        if( calcsignificance && 
                (method!="ProfiledLikelihood" and method!="Hybrid" and method!="ProfileLikelihood")
          ){cout<<"ERROR You are trying to use "<<method<<", which is not supported currently to calculate Significance"<<endl; exit(0);}
        if( !calcsignificance && 
                (method!="ProfiledLikelihood" and method!="Hybrid" and method!="Bayesian" and method!="FeldmanCousins" and method!="ProfileLikelihood" and method!="AsymptoticCLs" and method!="Asymptotic" and method!="ScanningMuFit" and method!="MaxLikelihoodFit" and method!="MultiDimFit")
          ){cout<<"ERROR You are trying to use "<<method<<", which is not supported currently to calculate limit "<<endl; exit(0);}
    }

    if(method=="Asymptotic") {bConstructAsimovbFromNominal = false;}
    else {bConstructAsimovbFromNominal = true;}
    if(isWordInMap("--bConstructAsimovbFromNominal", options)) bConstructAsimovbFromNominal = true;
    if(isWordInMap("--bConstructAsimovsbFromFit", options)) bConstructAsimovsbFromFit = true;

    tmpv = options["-v"]; if(tmpv.size()!=1) tmpv = options["--verbose"]; if(tmpv.size()!=1) tmpv = options["--debug"]; 
    if( tmpv.size()!=1 ) { debug = 0; }
    else debug = tmpv[0].Atoi();

    tmpv = options["-n"]; if(tmpv.size()!=1) tmpv = options["--name"];
    if( tmpv.size()!=1 ) { jobname = (datacards.size()>0?datacards[0]:"")+"_"+method; }
    else jobname = tmpv[0];

    tmpv = options["-D"]; if(tmpv.size()==0) tmpv = options["--dataset"];
    if( tmpv.size() ==0 ) { dataset = "data_obs"; }
    else { 
        dataset = tmpv[0];
        if(dataset!="data_obs" and dataset!="asimov_sb" and dataset!="asimov_b"){cout<<"ERROR: dataset option must be one of data_obs, asimov_sb and asimov_b"<<endl; exit(0);}
        if(dataset=="asimov_sb" and tmpv.size()==2) InjectingSignalStrength = tmpv[1].Atof();
    }


    tmpv = options["--InjectMu"]; 
    if( tmpv.size()!=1 ) { }
    else InjectingSignalStrength = tmpv[0].Atof();


    tmpv = options["--PhysicsModel"];
    if( tmpv.size()!=1 ) { sPhysicsModel = "StandardModelHiggs"; }
    else {
        sPhysicsModel = tmpv[0];
        if(sPhysicsModel!="StandardModelHiggs" and sPhysicsModel!="ChargedHiggs" and sPhysicsModel!="CvCfHiggs" and sPhysicsModel!="C5Higgs")  {
            cout<<"ERROR --PhysicsModel can only be StandardModelHiggs or ChargedHiggs or CvCfHiggs or C5Higgs as arg"<<endl; exit(0);}
    }

    tmpv = options["--doExpectation"]; 
    if( tmpv.size()!=1 ) { doExpectation = 0; }
    else doExpectation = tmpv[0].Atoi();

    tmpv = options["-t"]; if(tmpv.size()!=1) tmpv = options["--toys"];
    if( tmpv.size()!=1 ) { toys = 1000; }
    else toys = tmpv[0].Atoi();

    tmpv = options["-s"]; if(tmpv.size()!=1) tmpv = options["--seed"];
    if( tmpv.size()!=1 ) { seed = 1234; }
    else seed = tmpv[0].Atoi();

    tmpv = options["-S"]; if(tmpv.size()!=1) tmpv = options["--systematics"];
    if( tmpv.size()!=1 ) { systematics = 1; }
    else systematics = tmpv[0].Atoi();

    tmpv = options["-C"]; if(tmpv.size()!=1) tmpv = options["--CL"];
    if( tmpv.size()!=1 ) { CL = 0.95; }
    else {
        CL = tmpv[0].Atof();
        if(CL<=0 or CL>=1) {cout<<"ERROR  CL must be in the range (0, 1)"<<endl; exit(0);}
    }

    // Bayesian specific options
    tmpv = options["--prior"];
    if( tmpv.size()!=1 ) { prior = flat; }
    else {
        if(tmpv[0]=="flat"){prior=flat;}
        else if(tmpv[0]=="corr"){prior=corr;}
        else if(tmpv[0]=="1/sqrt(r)"){prior=prior_1overSqrtS;}
        else { cout<<"ERROR  supported baysian signal priors: flat , corr, 1/sqrt(r) "<<endl; exit(0);}
    }

    tmpv = options["-tB"]; if(tmpv.size()!=1) tmpv = options["--toysBayesian"];
    if( tmpv.size()!=1 ) { toysBayesian = 10000; }
    else {
        toysBayesian = tmpv[0].Atoi();
        if(toysBayesian<0) toysBayesian = 1;
    }

    tmpv = options["-tPB"]; if(tmpv.size()!=1) tmpv = options["--toysPreBayesian"];
    if( tmpv.size()!=1 ) { toysPreBayesian = 10; }
    else {
        toysPreBayesian = tmpv[0].Atoi();
        if(toysPreBayesian<0) toysPreBayesian = 1;
    }

    tmpv = options["--tries"]; 
    if( tmpv.size()!=1 ) { }
    else {
        nTries = tmpv[0].Atoi();
        if(nTries<0) nTries= 1;
    }

    tmpv = options["--flatPriorCropNsigma"]; 
    if( tmpv.size()!=1 ) { flatPriorCropNsigma= 3; }
    else { flatPriorCropNsigma = tmpv[0].Atof(); }

    if(isWordInMap("--bysRunFit", options)) bRunOnlyWithBestFittedNuisances_bayesian = true;

    tmpv = options["--bysInputMu"]; 
    if( tmpv.size()!=1 ) {}
    else { inputMu_bayesian = tmpv[0].Atof(); }

    // FeldmanCousins specific options
    tmpv = options["--lowerLimit"]; 
    if( tmpv.size()!=1 ) { lowerLimit = false; }
    else {
        if(tmpv[0].Atoi()==0) lowerLimit=false;
        else lowerLimit = true;
    }

    tmpv = options["--toysFactor"]; 
    if( tmpv.size()!=1 ) { toysFactor = 1; }
    else {
        toysFactor = tmpv[0].Atoi();
        if(toysFactor<1) toysFactor = 1;
    }

    tmpv = options["--adaptiveSampling"]; 
    if( tmpv.size()!=1 ) { adaptiveSampling = true; }
    else {
        if(tmpv[0].Atoi()==0) adaptiveSampling=false;
        else adaptiveSampling = true;
    }

    tmpv = options["--bQuickEstimateInitialLimit"]; 
    if( tmpv.size()!=1 ) { bQuickEstimateInitialLimit = true; }
    else {
        if(tmpv[0].Atoi()==0) bQuickEstimateInitialLimit=false;
        else bQuickEstimateInitialLimit = true;
    }

    tmpv = options["--initialRmin"]; 
    if( tmpv.size()!=1 ) { initialRmin = 1; }
    else {
        initialRmin = tmpv[0].Atof();
        //if(initialRmin<0) initialRmin = 1;
    }

    tmpv = options["--initialRmax"]; 
    if( tmpv.size()!=1 ) { initialRmax = 20; }
    else {
        initialRmax = tmpv[0].Atof();
        //if(initialRmax<=0) initialRmax = 20;
    }

    tmpv = options["-rMin"]; if(tmpv.size()!=1) tmpv = options["--rMin"];
    if( tmpv.size()!=1 ) { customRMin = 0; }
    else {
        customRMin = tmpv[0].Atof();
    }

    tmpv = options["-rMax"]; if(tmpv.size()!=1) tmpv = options["--rMax"];
    if( tmpv.size()!=1 ) { customRMax = 0; }
    else {
        customRMax = tmpv[0].Atof();
    }



    // Hybrid specific options
    tmpv = options["-tH"]; if(tmpv.size()!=1) tmpv = options["--toysHybrid"];
    if( tmpv.size()!=1 ) { toysHybrid = 10000; }
    else {
        toysHybrid = tmpv[0].Atoi();
        if(toysHybrid<0) toysHybrid = 100;
    }

    tmpv = options["--singlePoint"]; 
    if( tmpv.size()!=1 ) { singlePoint = false; }
    else { singlePoint = true; testR = tmpv[0].Atof(); }

    tmpv = options["--scanRs"]; 
    if( tmpv.size()!=1 ) { scanRs= 0; }
    else { scanRs= tmpv[0].Atoi(); nSteps = tmpv[0].Atoi(); }

    tmpv = options["--scanMs"]; 
    if( tmpv.size()!=1 ) { scanMs= 0; }
    else { scanMs= tmpv[0].Atoi(); nMsteps = tmpv[0].Atoi(); }

    tmpv = options["--testStat"]; 
    if( tmpv.size()!=1 ) { testStat = 1; }
    else {
        if(tmpv[0]=="LEP") testStat=1;
        else if(tmpv[0]=="TEV") testStat=2;
        else if(tmpv[0]=="Atlas") testStat=3;
        else if(tmpv[0]=="AtlasAllowMuHatNeg") testStat=32;
        else if(tmpv[0]=="LHC") testStat=5;
        else if(tmpv[0]=="PL") testStat=6;
        else {cout<<"ERROR Unimplemented testStat: "<<tmpv[0]<<". Supported: LEP, TEV, Atlas, AtlasAllowMuHatNeg "<<endl; exit(0); }
    }

    tmpv = options["--rule"]; 
    if( tmpv.size()!=1 ) { rule = 1; }
    else { 
        if(tmpv[0]=="CLs") rule=1;
        else if(tmpv[0]=="CLsb" or tmpv[0]=="CLsplusb") rule=2;
        else {cout<<"ERROR Unimplemented rule: "<<tmpv[0]<<". Supported: CLs,  CLsb "<<endl; exit(0); }
    }

    tmpv = options["--clsAcc"];
    if( tmpv.size()!=1 ) { clsAcc = 0.001; }
    else clsAcc = tmpv[0].Atof();

    tmpv = options["--rAbsAcc"];
    if( tmpv.size()!=1 ) { rAbsAcc = 0.01; }
    else rAbsAcc = tmpv[0].Atof();

    tmpv = options["--rRelAcc"];
    if( tmpv.size()!=1 ) { rRelAcc = 0.01; }
    else rRelAcc = tmpv[0].Atof();

    // ProfiledLikelihood specific options
    tmpv = options["--OneOrTwoSided"];
    if( tmpv.size()!=1 ) { oneside = 1; }
    else {
        oneside = tmpv[0].Atoi();
        if(oneside!=1 and oneside!=2)  {cout<<"ERROR --OneOrTwoSided can only have 1 or 2 as arg"<<endl; exit(0);}
    }
    tmpv = options["--PLalgorithm"];
    if( tmpv.size()!=1 ) { PLalgorithm = "Minos"; }
    else {
        PLalgorithm = tmpv[0];
        if(PLalgorithm !="Minos" and PLalgorithm!="Migrad")  {cout<<"ERROR --PLalgorithm can only be Minos or Migrad as arg"<<endl; exit(0);}
    }


    //  how to toss toys in building -2lnQ distribution
    tmpv = options["--tossToyConvention"];
    if( tmpv.size()!=1 ) { tossToyConvention = 0; }
    else {
        tossToyConvention = tmpv[0].Atoi();
    }
    // how to toss pseudo data for band
    tmpv = options["--tossPseudoDataConvention"];
    if( tmpv.size()!=1 ) { tossPseudoDataConvention = 0; }
    else {
        tossPseudoDataConvention = tmpv[0].Atoi();
    }

    if(isWordInMap("--freq", options)){
        tossToyConvention = 1;
        tossPseudoDataConvention = 1;
        UseBestEstimateToCalcQ = 0;
        if(calcsignificance)
            testStat = 6; // LHC  for one-sided upper limit
        else testStat=5; //PL for significance evaluation
    }
    // 
    if(isWordInMap("--UseBestEstimateToCalcQ", options))
        UseBestEstimateToCalcQ= 1; 

    if(isWordInMap("--bReadPars", options)) bReadPars = true;
    if(isWordInMap("--bWritePars", options)) bWritePars = true;
    if(isWordInMap("--bNotCalcQdata", options)) bNotCalcQdata = true;
    if(isWordInMap("--bNotCalcCLssbb", options)) bNotCalcCLssbb = true;
    if(isWordInMap("--bSaveM2lnQ", options)) bSaveM2lnQ = true;
    if(isWordInMap("--bSkipM2lnQ", options)) bSkipM2lnQ = true;
    if(isWordInMap("--bOnlyEvaluateQdata", options)) bOnlyEvaluateQdata= true;
    if(isWordInMap("--bWriteToys", options)) bWriteToys = true;

    tmpv = options["--ExpectationHints"];
    if( tmpv.size()!=1 ) { ExpectationHints= ""; }
    else {
        ExpectationHints= tmpv[0];
        //if(ExpectationHints!="Asymptotic" || ExpectationHints!="Bayesian" || ExpectationHints!="PreComputed") {
        //	cout<<"ERROR:  arg of ExpectationHints must be one of the following: "<<endl;
        //	cout<<"Asymptotic, Bayesian, PreComputed"<<endl;
        //}
    }


    tmpv = options["--scanRAroundBand"];
    if( tmpv.size()!=1 ) { scanRAroundBand= "all"; }
    else {
        scanRAroundBand= tmpv[0];
        if(scanRAroundBand!="all" && scanRAroundBand!="obs" && scanRAroundBand!="-1" && scanRAroundBand!="-2"
                &&scanRAroundBand!="0" && scanRAroundBand!="1" && scanRAroundBand!="2") {
            cout<<"ERROR:  arg of scanRAroundBand must be one of the following: "<<endl;
            cout<<"-2, -1, 0, 1, 2, obs, all"<<endl;
            exit(1);

        }
    }

    tmpv = options["--nToysForCLsb"];
    if( tmpv.size()!=1 ) { nToysForCLsb = -1; }
    else nToysForCLsb = tmpv[0].Atoi();
    tmpv = options["--nToysForCLb"];
    if( tmpv.size()!=1 ) { nToysForCLb = -1; }
    else nToysForCLb = tmpv[0].Atoi();

    tmpv = options["--fileFittedPars"];
    if( tmpv.size()!=1 ) { fileFittedPars = ""; }
    else fileFittedPars = tmpv[0];

    tmpv = options["--M2lnQGridFile"];
    if( tmpv.size()!=1 ) { bM2lnQGridPreComputed=false; }
    else {
        bM2lnQGridPreComputed = true;
        sFileM2lnQGrid= tmpv[0];
    }

    tmpv = options["--bCalcObservedLimit"];
    if( tmpv.size()!=1 ) { bCalcObservedLimit = true; }
    else {
        if(tmpv[0].Atoi()==0) bCalcObservedLimit = false;
        else bCalcObservedLimit = true;
    }

    if(isWordInMap("--bOnlyEvalCL_forVR", options)) bOnlyEvalCL_forVR = true;
    tmpv = options["-vR"];
    if( tmpv.size()==0 && bOnlyEvalCL_forVR) { cout<<"ERROR: args of -vR empty while bOnlyEvalCL_forVR=true"<<endl; exit(1); }
    else{
        if(tmpv.size()>=1)vR_toEval = GetListToEval("-vR",tmpv);
        //else if(tmpv.size()>1) vM_toEval.push_back(tmpv[1].Atof());
        if(tmpv.size() >0)sbinningMU = tmpv[0];
    }

    tmpv = options["-vM"];
    if( scanMs and tmpv.size()==0) { cout<<"ERROR: args of -vM empty "<<endl; exit(1); }
    else{
        //vM_toEval = GetListToEval("-vM",tmpv);
        if(tmpv.size()>=1)vM_toEval = GetListToEval("-vM",tmpv);
        //else if(tmpv.size()>1) vM_toEval.push_back(tmpv[1].Atof());
        if(tmpv.size() >0)sbinningMH = tmpv[0];
    }


    if(makeCLsBands>=2){
        bSkipM2lnQ = 1; doExpectation = 0;
    }

    tmpv = options["--minuitSTRATEGY"]; 
    if( tmpv.size()!=1 ) { minuitSTRATEGY= 0; }
    else { minuitSTRATEGY = tmpv[0].Atoi(); }

    tmpv = options["--minuitPrintLevel"]; 
    if( tmpv.size()!=1 ) { minuitPrintLevel= -1; }
    else { minuitPrintLevel = tmpv[0].Atoi(); }


    tmpv = options["--minuitTolerance"]; 
    if( tmpv.size()!=1 ) { minuitTolerance= 0.001; }
    else { minuitTolerance= tmpv[0].Atof(); }

    tmpv = options["--nuisancesRange"]; 
    if( tmpv.size()!=1 ) { nuisancesRange = 5.; }
    else { nuisancesRange = tmpv[0].Atof(); }

    tmpv = options["--maximumFunctionCallsInAFit"]; 
    if( tmpv.size()!=1 ) {}
    else { maximumFunctionCallsInAFit = tmpv[0].Atoi(); }

    tmpv = options["--ManualNuisanceFeeding"]; // FIXME jaehyeok 
    if( tmpv.size()!=1 ) { ManualNuisanceFeeding=false; }
    else { 
        ManualNuisanceFeeding = true;
        FileNuisanceFeeding = tmpv[0]; 
    }

    // for ScanningMuFit
    if(isWordInMap("--bAlsoExtract2SigmErrorBar", options)) bAlsoExtract2SigmErrorBar= true;
    if(isWordInMap("--bConstrainMuPositive", options)) bConstrainMuPositive= true;

    // for LHC-CLs method
    if(isWordInMap("--bFixNuisancsAtNominal", options)) bFixNuisancsAtNominal = true;
    if(isWordInMap("--bFixNuisancesAtBestFit", options)) bFixNuisancesAtBestFit= true;
    if(isWordInMap("--fastScan", options)) bFixNuisancesAtBestFit= true;
    cout<<"MakeSureYouAreHere bFixNuisancesAtBestFit ="<<bFixNuisancesAtBestFit<<endl;


    if(isWordInMap("--bDumpFitResults", options)) bDumpFitResults = true; 
    if(isWordInMap("--bRunMinuitContour", options)) bRunMinuitContour = true; 

    tmpv = options["--PrintParameter"]; 
    if( tmpv.size()==2 ) { PrintParameterFrom = tmpv[0].Atoi(); 
        PrintParameterTo = tmpv[1].Atoi(); }
    else if(tmpv.size()==1){ PrintParameterFrom = tmpv[0].Atoi(); 
        PrintParameterTo = tmpv[0].Atoi(); }
    else { PrintParameterFrom = -1; PrintParameterTo=-1; }

    tmpv = options["--maxsets_caching"]; 
    if( tmpv.size()!=1 ) { maxsets_forcaching = 0; }
    else { maxsets_forcaching = tmpv[0].Atoi(); }

    tmpv = options["--POIs"];
    if(tmpv.size()==4 or tmpv.size()==8){ // currently support only at most two additional POIs 
        double initv = tmpv[1].Atof();
        double minv = tmpv[2].Atof();
        double maxv = tmpv[3].Atof();
        if(minv>maxv or initv<minv or initv>maxv) { cout<<" POI args are wrong: name "<<tmpv[0]<<" initV "<<tmpv[1]<<" minV "<<tmpv[2]<<" maxV "<<tmpv[3]<<endl; exit(1);}
        structPOI poi(tmpv[0], initv, 0, 0, minv, maxv);	
        addtionalPOIs.push_back(poi);
        if(tmpv.size()==8){
            initv = tmpv[5].Atof();
            minv = tmpv[6].Atof();
            maxv = tmpv[7].Atof();
            if(minv>maxv or initv<minv or initv>maxv) { cout<<" POI args are wrong: name "<<tmpv[4]<<" initV "<<tmpv[5]<<" minV "<<tmpv[6]<<" maxV "<<tmpv[7]<<endl; exit(1);}
            structPOI poi2(tmpv[4], initv, 0, 0, minv, maxv);	
            addtionalPOIs.push_back(poi2);
        }	
    }else if(tmpv.size()==0){}
    else{ cout<<" --POIs args[4 or 8 args]"<<endl;  exit(1);}

    tmpv= options["--Couplings"];
    vCouplingsDef = tmpv;

    tmpv= options["--smbr"];
    if(tmpv.size()!=0)sSMBrFile = tmpv[0];

    tmpv= options["--CV"];
    for(int i=0; i<tmpv.size(); i++) CVsetting.push_back(tmpv[i]);
    tmpv= options["--CF"];
    for(int i=0; i<tmpv.size(); i++) CFsetting.push_back(tmpv[i]);

    if(isWordInMap("--scanCvCf", options)) scanCvCf= true; 

    tmpv = options["-vCv"];
    if( scanCvCf and tmpv.size()==0) { cout<<"ERROR: args of -vCv empty "<<endl; exit(1); }
    else{
        if(tmpv.size()==1)vCv_toEval = GetListToEval("-vCv",tmpv);
        else if(tmpv.size()>1) vCv_toEval.push_back(tmpv[1].Atof());
        if(tmpv.size() >0)sbinningCV = tmpv[0];
    }
    tmpv = options["-vCf"];
    if( scanCvCf and tmpv.size()==0) { cout<<"ERROR: args of -vCf empty "<<endl; exit(1); }
    else{
        if(tmpv.size()==1)vCf_toEval = GetListToEval("-vCf",tmpv);
        else if(tmpv.size()>1) vCf_toEval.push_back(tmpv[1].Atof());
        if(tmpv.size()>0)sbinningCF = tmpv[0];
    }

    tmpv= options["--Cvv"];
    for(int i=0; i<tmpv.size(); i++) CvvSetting.push_back(tmpv[i]);
    tmpv= options["--Cgg"];
    for(int i=0; i<tmpv.size(); i++) CggSetting.push_back(tmpv[i]);
    tmpv= options["--Ctt"];
    for(int i=0; i<tmpv.size(); i++) CttSetting.push_back(tmpv[i]);
    tmpv= options["--Cbb"];
    for(int i=0; i<tmpv.size(); i++) CbbSetting.push_back(tmpv[i]);
    tmpv= options["--Cglgl"];
    for(int i=0; i<tmpv.size(); i++) CglglSetting.push_back(tmpv[i]);
    if(isWordInMap("--scanCX", options)) scanCX= true; 
    tmpv = options["-vCvv"];
    if( scanCX and tmpv.size()==0) { cout<<"ERROR: args of -vCvv empty "<<endl; exit(1); }
    else{
        vCvv_toEval = GetListToEval("-vCvv",tmpv);
    }
    if( scanCX and tmpv.size()==0) { cout<<"ERROR: args of -vCgg empty "<<endl; exit(1); }
    else{
        vCgg_toEval = GetListToEval("-vCgg",tmpv);
    }
    if( scanCX and tmpv.size()==0) { cout<<"ERROR: args of -vCtt empty "<<endl; exit(1); }
    else{
        vCtt_toEval = GetListToEval("-vCtt",tmpv);
    }
    if( scanCX and tmpv.size()==0) { cout<<"ERROR: args of -vCbb empty "<<endl; exit(1); }
    else{
        vCbb_toEval = GetListToEval("-vCbb",tmpv);
    }
    if( scanCX and tmpv.size()==0) { cout<<"ERROR: args of -vCglgl empty "<<endl; exit(1); }
    else{
        vCglgl_toEval = GetListToEval("-vCglgl",tmpv);
    }


    tmpv=options["--ErrEstAlgo"]; 
    if(tmpv.size()!=0) { ErrEstAlgo=tmpv[0]; if(ErrEstAlgo!="Minos" and ErrEstAlgo!="Bisect") {cout<<" ErrEstAlgo only supports Minos and Bisect"<<endl; exit(1);}}

    tmpv=options["--mainPOI"];
    if(tmpv.size()==1){
        mainPOI = tmpv[0];	
    }
    tmpv=options["--scanPOI"]; // use with --mainPOI 
    if(tmpv.size()){
        scanPOI=true;
        if(tmpv.size()>=1)scanPOIrange = GetListToEval("POI",tmpv);
        if(tmpv.size() >0)sbinningPOI = tmpv[0];
    }
    tmpv=options["--PrintPdfEvlCycle"];
    if(tmpv.size()!=0) {PrintPdfEvlCycle=tmpv[0].Atof();}
    tmpv=options["--PrintFuncCallCycle"];
    if(tmpv.size()!=0) {_printFuncCallCycle = tmpv[0].Atof();}

    if(isWordInMap("--doMemoryCheck", options)) doMemoryCheck=true;
    if(isWordInMap("--DoNotRunGlobalFit", options)) DoNotRunGlobalFit=true;
    if(isWordInMap("--GenerateToysAtBestFitSB", options)) GenerateToysAtBestFitSB=true;
    if(isWordInMap("--NoErrorEstimate", options)) NoErrorEstimate =true;

    tmpv=options["--ndof"];
    if(tmpv.size()!=0) {NumDegreeOfFreedom=tmpv[0].Atof();}

    tmpv=options["--intRemoveBins"];
    if(tmpv.size()!=0) {intRemoveBins=tmpv[0].Atoi();}

    tmpv=options["--RebinObservables"];
    if(tmpv.size()!=0){
        if(tmpv.size() % 4 ==0){
            for(int i=0; i<tmpv.size(); i+=4){
                rebin_observables.push_back(tmpv[i]);
                rebin_observables_nbins.push_back(tmpv[i+1].Atoi());
                rebin_observables_xmin.push_back(tmpv[i+2].Atof());
                rebin_observables_xmax.push_back(tmpv[i+3].Atof());
            }
        }else{
            cout<<"ERROR.  --RebinObservables should have 4/8/12/... args,  --> obsname nbin xmin xmax "<<endl; exit(1); 
        }
    }

    if(isWordInMap("--newExpected", options)) newExpectedAsymBand= true;

    tmpv=options["--ToysAtDifferentSignalStrength"];
    if(tmpv.size()!=0) {ToysAtDifferentSignalStrength=tmpv[0].Atoi();}


    if(isWordInMap("--useHist", options)) useHist = true; 

    tmpv = options["--loadToysFromFile"]; 
    if( tmpv.size()!=1 ) { }
    else sToysFromFile = tmpv[0];

    tmpv = options["--MinimizingApproach"]; 
    if( tmpv.size()!=1 ) { }
    else sMinimizingApproach = tmpv[0];


    tmpv = options["--FixParWhenGenToys"];
    if(tmpv.size() % 2 ==0){ // currently support only at most two additional POIs 
        for(int i=0; i<tmpv.size(); i++){
            vsParToBeFixedWhenGenToys.push_back(tmpv[i]);
            vdParToBeFixedWhenGenToys.push_back(tmpv[i+1].Atof());
            i++;
        }
    }else{ cout<<" --FixParWhenGenToys parName parVal [2,4,6,8.... args]"<<endl;  exit(1);}

    if(isWordInMap("--scanNuisance", options)){
        tmpv = options["--scanNuisance"];
        if( tmpv.size()<2) { cout<<"ERROR: args of --scanNuisance should be \" NuisanceName [min,max,step] \" "<<endl; exit(1); }
        else{
            vNuis_toEval= GetListToEvalForNuis(scanNuisance, tmpv);
        }
    }



    printf("\n\n[ Summary of configuration in this job: ]\n");
    if(sPhysicsModel!="StandardModelHiggs")cout<<" PhysicsModel:  "<<sPhysicsModel<<endl;
    cout<<"  Calculating "<<(calcsignificance?"significance":"limit")<<" with "<<method<<" method "<<endl;
    if(HiggsMass>0) cout<<" higgs mass = "<<HiggsMass<<endl;
    if(!bCalcObservedLimit) cout<<" not calc observed one"<<endl;
    cout<<"  datacards: "; for(int i=0; i<datacards.size(); i++) cout<<datacards[i]<<" "; cout<<endl;
    cout<<"  "<<(systematics?"use systematics":"not use systematics")<<endl;
    cout<<"  dataset is "<<dataset<<endl;
    if(!calcsignificance) cout<<"  target confidence level = "<<CL<<endl;
    cout<<"  "<<(doExpectation?"also calc expectation bands":"do not calc expectation bands")<<endl;
    if(doExpectation){
        cout<<"  number of outcomes to build expecation bands: "<<toys<<endl;
        if(tossPseudoDataConvention){
            cout<<"  tossPseudoDataConvention = "<<tossPseudoDataConvention<<endl;
        }
    }

    if(bOnlyEvalCL_forVR){
        cout<<" bOnlyEvalCL_forVR: ";
        for(int i=0; i<vR_toEval.size(); i++) cout<< vR_toEval[i] << " ";
        cout<<endl;
    }

    if(method=="Bayesian") { 
        cout<<"  prior = ";
        if(prior==flat) cout<<"flat"<<endl;
        else if(prior==corr) cout<<"corr"<<endl;
        else if(prior==prior_1overSqrtS) cout<<"1/sqrt(r)"<<endl;
        else cout<<"Unknow"<<endl;

        cout<<"  toysPreBayesian = "<<toysPreBayesian<<endl;
        cout<<"  toysBayesian = "<<toysBayesian<<endl;

        if(flatPriorCropNsigma!=3) cout<<" flatPriorCropNsigma: "<<flatPriorCropNsigma<<endl;
    }else if(method=="ProfiledLikelihood" or method=="ProfileLikelihood"){
        if(!calcsignificance)cout<<(oneside==1?"  one sided":"  two sided")<<endl;
        cout<<"  algorithm to extract limit: "<<PLalgorithm<<endl;
    }else if(method=="FeldmanCousins"){
        cout<<"  calc "<<(lowerLimit?"lower limit":"upper limit")<<endl;
        cout<<"  "<<(adaptiveSampling?"use adaptiveSampling":"do not use adaptiveSampling")<<endl;
        cout<<"  toysFactor = "<<toysFactor<<endl;
        if(adaptiveSampling==false) cout<<"   number of toys per point = "<<toysHybrid<<endl;
        cout<<"  do "<<(bQuickEstimateInitialLimit?"":"not")<<" estimate initial limit with bayesian method"<<endl;
        if(!bQuickEstimateInitialLimit){
            cout<<"   initialRmin = "<<initialRmin<<", initialRmax = "<<initialRmax<<endl;
            if(initialRmax<initialRmin) {cout<<"  initialRmin can't be > initialRmax"<<endl; exit(0);}
        }
    }else if(method=="Hybrid"){
        if(calcsignificance==false){
            cout<<"  testStat = "; if(testStat==1) cout<<"LEP"; if(testStat==2)cout<<"TEV"; if(testStat==3)cout<<"Atals"; 
            if(testStat==32)cout<<"AtlasAllowMuHatNeg"; if(testStat==5)cout<<"LHC"; if(testStat==6) cout<<"PL";
            cout<<endl;
            cout<<"  rule     = "; if(rule==1) cout<<"CLs"; if(rule==2)cout<<"CLsb";
            cout<<endl;
        }
        cout<<"  number of toys to build Q distribution = "<<toysHybrid<<endl;
        if(tossToyConvention){
            cout<<" tossToyConvention = "<<tossToyConvention<<endl;
        }
        if(singlePoint){
            cout<<" testing single point with r = "<<testR<<endl;
        }
        cout<<"  UseBestEstimateToCalcQ: "<<(UseBestEstimateToCalcQ?"true":"false")<<endl;

        if(bM2lnQGridPreComputed){
            cout<<" Use pre-computed grid (-2lnQ) from file : "<<sFileM2lnQGrid<<endl;
        }

        if(bFixNuisancsAtNominal) cout<<" The nuisances are fixed at nominal values instead of being fitted to data "<<endl;
    }else if(method=="Asymptotic"){
        cout<<"  RUNNING Asymptotic limit"<<endl;
        cout<<"  for expected : "<<scanRAroundBand<<endl;
    }else if(method=="AsymptoticCLs"){
        cout<<"  RUNNING Asymptotic CLs"<<endl;
    }else if(method=="ScanningMuFit" or method=="MaxLikelihoodFit" or method=="MultiDimFit"){
        cout<<"  RUNNING MaxLikelihoodFit"<<endl;
    }
    if(makeCLsBands) cout<<"  makeCLsBands = "<<makeCLsBands<<endl;

    cout<<"  random number generator seed: "<<seed<<endl;
    cout<<"  debug level = "<<debug<<endl;
    cout<<"  job name: "<<jobname<<endl;
    cout<<"  plotLevel = "<<bPlots<<endl;
    if(librariesToBeLoaded.size()>0){
        cout<<"    custom libraries to be loaded: "<<endl;
        for(int i=0; i<librariesToBeLoaded.size(); i++) cout<<" "<<librariesToBeLoaded[i]<<endl;
    }
    if(isWordInMap("minuitSTRATEGY", options))cout<<"  minuitSTRATEGY="<<minuitSTRATEGY<<endl;
    if(isWordInMap("minuitTolerance", options))cout<<"  minuitTolerance="<<minuitTolerance<<endl;
    if(isWordInMap("maximumFunctionCallsInAFit", options))cout<<"  maximumFunctionCallsInAFit="<<maximumFunctionCallsInAFit<<endl;
    if(isWordInMap("nuisancesRange", options))cout<<"  nuisancesRange = +/- "<<nuisancesRange<<" sigma "<<endl;

    if(bMultiSigProcShareSamePDF) cout<<" MultiSigProcShareSamePDF "<<endl;
    if(bForceSymmetryError) cout<<" ForceSymmetryError"<<endl;
    cout<<endl<<endl;

    if(isWordInMap("--bRedefineObservableRange", options)){
        bRedefineObservableRange = true;
        tmpv = options["--ObservableRangeMax"]; 
        if( tmpv.size()!=1 ){} 
        else { ObservableRangeMax = tmpv[0].Atof(); }
        tmpv = options["--ObservableRangeMin"]; 
        if( tmpv.size()!=1 ){} 
        else { ObservableRangeMin = tmpv[0].Atof(); }
        tmpv = options["--ObservableBins"]; 
        if( tmpv.size()!=1 ){} 
        else { ObservableBins= tmpv[0].Atoi(); }
        cout<< " ObservableBins = " << ObservableBins <<endl;

    }

    fflush(stdout);	

    // check duplicate options ,  non-exist options 
}
void PrintHelpMessage(){
    printf("Usage: ./lands.exe [options] \n");                                                                                                                 
    printf("Allowed options: \n");                                                                                                                 
    printf("-h [ --help ]                         Produce help message \n"); 
    printf("-v [ --verbose ] arg (=0)             Verbosity level \n"); 
    printf("-L [ --LoadLibraries]                 custom libs to be loaded \n"); 
    printf("--plot 	                              make plots when appropriate \n");
    printf("--PhysicsModel arg (=StandardModelHiggs) could be StandardModelHiggs, ChargedHiggs, CvCfHiggs, C5Higgs  \n");
    printf("-n [ --name ] arg                     Name of the job,  default is \"datacard\"+\"method\" \n"); 
    printf("-d [ --datacards ] args               Datacard files,  can contain \"*, ?\" \n"); 
    printf("-D [ --dataset ] arg (=data_obs)      Dataset for observed limit,  data_obs, asimov_b, \"asimov_sb [mu]\"\n"); 
    printf("-M [ --method ] arg                   Method to extract upper limit. Supported methods are: Bayesian, FeldmanCousins, Hybrid, ProfiledLikelihood \n"); 
    printf("		                      ScanningMuFit, MaxLikelihoodFit, MultiDimFit\n"); 
    printf("-m [ --mass ] arg                     input higgs mass \n"); 
    printf("--doExpectation arg (=0)              i.e calc expected bands and mean/median values     \n"); 
    printf("--bOnlyEvalCL_forVR                   only evaluate CL for each pseudo data, should with -vR option \n");
    printf("-vR args                              sth like \"1.2 1.3 1.4 [1.5,3.0,x1.05] [3.0,10.0,0.5]\" \n");
    printf("--readLimitsFromFile arg (fileName)   i.e read expectations from merged files \n");
    printf("--bCalcObservedLimit arg (=1)         i.e calc observed limit values     \n"); 
    printf("-t [ --toys ] arg (=1000)             Number of Toy MC extractions for expectation \n"); 
    printf("--tossPseudoDataConvention arg (=0)   choose convention for tossing toys to build -2lnQ distribution. 0 (LEP, TEVATRON) or 1 (LHC)\n");
    printf("-s [ --seed ] arg (=1234)             Toy MC random seed \n"); 
    printf("-S [ --systematics ] arg (=1)         if 0, then will not use systematics  \n"); 
    printf("-C [ --cl ] arg (=0.95)               Confidence Level \n"); 
    printf("--significance arg (=0)               Compute significance instead of upper limit,  supported methods: ProfiledLikelihood and Hybrid (CLb) \n"); 
    printf("--calcPValueFromFile arg (fileName)   i.e calc p-values from merged files \n");
    printf(" \n"); 
    printf("Bayesian specific options: \n"); 
    printf("--prior arg (=flat)            	      Prior to use: \'flat\' (default), \'1/sqrt(r)\', \'corr\' \n"); 
    printf("-tB [ --toysBayesian ] arg (=10000)   number of toys used to average out nuisance paramereters in Bayesian method     \n"); 
    printf("-tPB [ --toysPreBayesian ] arg (=100)   number of toys used to average out nuisance paramereters in Bayesian pre-estimation \n"); 
    printf("--tries arg (=1)                      number of tries for observed limit, if more than 1, will print out the average value and standard deviation of those tries\n");
    printf("--flatPriorCropNsigma arg (=3)        fit to data and crop the flat prior range by Nsigma\n");
    printf("--bysRunFit --bysInputMu [double] \n");
    printf(" \n"); 
    printf("FeldmanCousins specific options: \n"); 
    printf("--lowerLimit arg (=0)                 Compute the lower limit instead of the upper limit \n"); 
    printf("--toysFactor arg (=1)                 Increase the toys per point by this factor w.r.t. the minimum from adaptive sampling \n"); 
    printf("--adaptiveSampling arg (=1)           currently only implemented in FeldmanCousins,  turn on (=1) by default.  =0 off \n"); 
    printf("--bQuickEstimateInitialLimit arg (=1) quickly estimate initial limit from bayesian technique, turn off by 0  \n"); 
    printf("--initialRmin arg (=1)                only effective when bQuickEstimateInitialLimit=0 \n"); 
    printf("--initialRmax arg (=20)               only effective when bQuickEstimateInitialLimit=0 \n"); 
    printf(" \n"); 
    printf("Hybrid specific options: \n"); 
    printf("-tH [ --toysHybrid ] arg (=10000)     Number of Toy MC extractions to compute CLs+b, CLb and CLs \n"); 
    printf("--clsAcc arg (=0.001)                 Absolute accuracy on CLs to reach to terminate the scan \n"); 
    printf("--rAbsAcc arg (=0.01)                 Absolute accuracy on r to reach to terminate the scan \n"); 
    printf("--rRelAcc arg (=0.01)                 Relative accuracy on r to reach to terminate the scan \n"); 
    printf("--rule arg (=CLs)                     Rule to use: CLs, CLsb \n"); 
    printf("--testStat arg (=LEP)                 Test statistics: LEP, TEV, Atlas, AtlasAllowMuHatNeg, LHC (only for limit),  PL (only for significance). \n"); 
    printf("--singlePoint arg (=float)            Just compute CLsb/CLb/CLs values for a given value of r \n"); 
    printf("--scanRs arg (=numSteps)              scanning CLs vs. r,  r from initialRmin to initialRmax with numSteps \n"); 
    printf("--tossToyConvention arg (=0)          choose convention for tossing toys to build -2lnQ distribution. 0 (LEP, TEVATRON) or 1 (LHC)\n");
    printf("--UseBestEstimateToCalcQ arg (=1)     0: randomized nuisances; 1: use best estimate of nuisances  --> to calc Q\n");
    printf("--freq                                shorcut to the configuration of LHC-type frequestist method \n");
    printf("                                      (i.e. --tossToyConvention 1 --UseBestEstimateToCalcQ 0  --tossPseudoDataConvention 1 --testStat LHC) \n");

    printf("--bWriteToys (=0)\n");
    printf("--bReadPars (=0)\n");
    printf("--bWritePars (=0)\n");
    printf("--bNotCalcQdata (=0)\n");
    printf("--bOnlyEvaluateQdata (= 0)\n");
    printf("--bNotCalcCLssbb (= 0)\n");
    printf("--bSaveM2lnQ (= 0)\n");
    printf("--bSkipM2lnQ (= 0)\n");
    printf("--nToysForCLb (= -1)\n");
    printf("--nToysForCLsb (= -1)\n");
    printf("--fileFittedPars (=\"\")\n");
    printf("--M2lnQGridFile (= \"filename\")\n");

    printf("--makeCLsBands arg (=0)               0: skip it;  1: run toys and make it;  2: read toys and make it, with --M2lnQGridFile\n");

    printf(" \n"); 
    printf("ProfiledLikelihood specific options: \n"); 
    printf("--OneOrTwoSided arg (=2)              1 sided limit -lnL = 1.345;  2 sided limit -lnL = 1.921 \n"); 
    printf("--PLalgorithm arg (=Minos)            algorithms for ProfileLikelihood approximation: Minos, Migrad \n"); 
    printf("--minuitSTRATEGY arg (=0)             0: no calculation of secondary derivative, 1: yes \n"); 
    printf("--minuitTolerance arg (=0.001)        tolerance in minuit fit \n"); 
    printf("--maximumFunctionCallsInAFit arg (=5000)  \n"); 
    printf("--nuisancesRange arg (=5.)            nuisances viriation range allowed, when the significan is >> 5, then, need set it to be 10 or larger)\n"); 

    printf("\n");
    printf("random arguments, need to be cleaned\n");
    printf("--bForceSymmetryError                 make the asymmetry uncertainties  be symmetry in a conservative way inside code\n");
    printf("--bMultiSigProcShareSamePDF           e.g. in hhiggs4l, the ggH, WH, ZH, ttH and VBF have the same signal pdf, when evaluating LL, we can avoid duplicated calculations \n");
    printf("--bFixNuisancsAtNominal               fix nuisances to nominal value instead of fitting to data in case of LHC-CLs\n");
    printf("--scanRAroundBand  (=all)            scan signal strengths around band (-2, -1, 0, 1, 2, obs, all)\n");
    printf("--bConstructAsimovbFromNominal	      construct asimov_b dataset from nominal nuisances   (default is from fitted nuisances) \n");
    printf("--bConstructAsimovsbFromFit	      construct asimov_sb dataset from global fitted nuisances   (default is from nominal nuisances) \n");
    printf("--bRedefineObservableRange, --ObservableRangeMin, --ObservableRangeMax, --ObservableBins :  for hhiggs4l pdf extraction\n");
    printf("--bDumpFitResults                     dump fit results and also the fitted shape if there is any shape input (only binned supported) \n");
    printf("--ManualNuisanceFeeding               feed post-fit nuisances from another model to get post-fit shapes for fittedShape_floatMu (definition of nuisances should be same between two models ) \n");
    printf("--maxsets_caching arg (=0)            number of sets of cached pdf values,  the larger it is, the larger memory required    \n");
    printf("--PrintParameter arg1 arg2 (= -1 -1)    print intermediate values of the parameters specified in the range [arg1, arg2]    \n");
    printf("--rMin arg (=0)                       force signal strength to be >= rMin \n");
    printf("--rMax arg (=0)                       force signal strength to be <= rMax \n");
    printf("-vM args                              sth like \"1.2 1.3 1.4 [1.5,3.0,x1.05] [3.0,10.0,0.5]\", for scanning mass  \n");
    printf("--scanMs arg (=0)                     scanning LLR vs. M \n"); 
    printf("--bAlsoExtract2SigmErrorBar           for Maximum LL fit\n"); 
    printf("--minuitPrintLevel arg(=-1)           minuit print level in fit, for debugging \n"); 
    printf("--bRunMinuitContour                   Run MNCONT for two POIs, requiring computing time \n");
    printf("--POIs args (='NAME InitVal Min Max') Specify other POIs (like MH), which should be in the input workspace\n");
    printf("--Couplings (='ggH:[0.9,0,5]:*|ggH,hzz4l_4mu|ggH,hzz2l2nu_ee|ggH  parName:[initVal,min,max]:FinalState|ProductionMode')\n");
    printf("--mainPOI arg                         remove other POIs, let them be nuisances \n");
    printf("--smbr arg (=FilePath)                for SMHiggsBuilder, path of the file containing SM Br\n");
    printf("--CV init min max --CF init min max   settings for CvCfHiggs,  i.e. fit range\n");
    printf("--scanCvCf                            to be used with -vCv, -vCf\n");
    printf("-vCv/-vCf args                        sth like \"[3.0,10.0,0.5]\", for scanning CV vs. CF, if there are second arg, then just scan that only value\n");
    printf("--ErrEstAlgo arg                      Minos, Bisect\n");
    printf("--PrintPdfEvlCycle arg                for debugging\n");
    printf("--PrintFuncCallCycle arg                for debugging\n");
    printf("--doMemoryCheck  		      only run upto model construction, to terminate the valgrind, more usage in future\n");
    printf("--DoNotRunGlobalFit                   for scanning POIs space, don't run global fit, and save the fMin instead of delta Q\n");
    printf("--RebinObservables name nbin xmin xmax .........  for rebin observables  \n");
    printf("--ndof arg (=1)                       number of degree of freedom,   for extracting errors in n dimensions \n");
    printf("--GenerateToysAtBestFitSB 	      for generating toys at the best fit for s+b (floating mu), for signal hypotheses separation\n"); 
    printf("--NoErrorEstimate 		      do not call minos to calculate error bar \n");
    printf("--bFixNuisancesAtBestFit or --fastScan     fix nuisances to global best fit value when doing 1D poi scan without systematics (not removed from cards)  \n");
    printf("--intRemoveBins arg (=0)               0: remove only bins with total bkg<=0,  1: remove bins with total sig<=0,  2: both \n"); 
    printf("--InjectMu arg(=1)                    for InjectingSignalStrength, for MLLxsBands toys,  and asimov_sb,  it will overwrite the one associated with asimov_sb \n"); 
    printf("--ToysAtDifferentSignalStrength arg(=0)  for MLLxsBands toys,  0 for b-only toys, 1 for toys at InjectMu, and 2 for the best fitted mu\n");
    printf("--useHist                              With this, take hist as a channel, no expansion, no bin-by-by stat err\n");   
    printf("--newExpected 			     new definition of Asymptotic CLs bands \n"); 
    printf("--loadToysFromFile string            use with doExpectation,  reading saved toys from file \n");
    printf("--MinimizingApproach string            currently supports Normal, Cascade\n");
    printf("--scanNuisance nuisName [min,max,step] scanning Nuisance, now only work in MultiDimFit \n");


    printf(" \n");
    printf("------------------some comand lines-----------------------------------------------\n");
    printf("            *run LHC-CLs from hints extracted from asymptotic limits*\n");
    printf("lands.exe -d cards -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n JOBNAME\n");
    printf("            *extract CLs values with bands from precomputed grid*\n");
    printf("lands.exe  -M Hybrid --makeCLsBands 2 --M2lnQGridFile fileContainsM2lnQ_R=1.root\n");
    printf("            *extract UL expectation from merged file contains limit tree*\n");
    printf("lands.exe --doExpectation 1 --readLimitsFromFile fileContains_R_distribuition.root\n");
    printf("            *calculate p-value from merged file contains -2lnQ of b-only toys*\n");
    printf("lands.exe --significance 1 -M Hybrid --testStat PL/LEP --calcPValueFromFile fileContains_-2lnQ@bonly.root\n");
    printf("            *prepare grid of -2lnQ *\n");
    printf("lands.exe -M Hybrid --freq --bNotCalcCLssbb --bSaveM2lnQ --nToysForCLsb 1500 --nToysForCLb 500 --singlePoint rValue  -n JobName -s Seed -d cards*.txt\n");
    printf("            *make limit bands with throwing pseudo data, \n");
    printf("            *for each pseudo data, evaluate testStat for all available grid Rs from grid file*\n");
    printf("lands.exe -d datacards*txt -M Hybrid --freq --bCalcObservedLimit 0  --M2lnQGridFile grid.root  --doExpectation 1 -t XXX --bSkipM2lnQ 1\n");
    printf("            *asymptotic cls bands\n");
    printf("lands.exe -d data.txt -M Asymptotic -D asimov_b\n");
    printf("            *scanning mu with fit and make plots\n");
    printf("lands.exe -d data.txt -M ScanningMuFit --scanRs 20 --initialRmin 0 --initialRmax 2  --plot\n");
    printf("            *maximum LL fit with floating mu and mu=0,  printing fitted results with -v 1 \n");
    printf("lands.exe -d data.txt -M MaxLikelihoodFit --scanRs 1 -vR 0 -v 1 --minuitSTRATEGY 1\n");

    exit(0);
}
