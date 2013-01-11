
#include "Enums.h"

void runZFakeLooper();

void doZFakeRate() {

    //
    // load libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../../Core/SmurfTree.h+");
    gSystem->Load("libSmurfFakeLooper.so");

    bool runEle = false;
    bool runMu = true;

    // Moriond
    runZFakeLooper(runEle, runMu, MORIOND);

}

void runZFakeLooper(bool runEle, bool runMu, Era era)
{

    TString eraName = "Moriond";

    //
    // configure looper
    //

    ZFakeLooper *looper = new ZFakeLooper();
    looper->setLumi(19467.0);

    std::vector<unsigned int> ptThresholds;

    // data sample
    TChain *ch_data_m = new TChain("tree");
    ch_data_m->Add("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/data_3l.root");
    TChain *ch_wz = new TChain("tree");
    ch_wz->Add("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/wz.root");
    TChain *ch_zz = new TChain("tree");
    ch_zz->Add("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/zz.root");
    TChain *ch_wglll = new TChain("tree");
    ch_wglll->Add("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/wglll.root");


    //
    // run muon fake rates
    //

    if (runMu) {
        ptThresholds.clear();
        ptThresholds.push_back(0);
        ptThresholds.push_back(5);
        ptThresholds.push_back(10);
        ptThresholds.push_back(15);
        ptThresholds.push_back(20);
        ptThresholds.push_back(25);
        ptThresholds.push_back(30);
        ptThresholds.push_back(50);
        ptThresholds.push_back(80);
        looper->Loop(true,  ch_data_m,  "Data_MuonFakeRate_M2", ptThresholds);
        looper->Loop(false,  ch_wz,  "WZ_MuonFakeRate_M2", ptThresholds);
        looper->Loop(false,  ch_zz,  "ZZ_MuonFakeRate_M2", ptThresholds);
        //looper->Loop(false,  ch_wglll,  "wglll_MuonFakeRate_M2", ptThresholds);
    }

    //
    // save and manipulate histograms
    //

    const TString outFile = Form("histos_ZFakeLooper_%s.root", eraName.Data());
    saveHist(outFile.Data());
    deleteHistos();

    //
    // tidy up
    //

    delete looper;
    delete ch_data_m;
    delete ch_wz;
    delete ch_zz;
    delete ch_wglll;
}

