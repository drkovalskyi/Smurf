#include "Enums.h"

void doPlot(TString era, bool defaultOnly = false) {

    //
    // load libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gROOT->ProcessLine(".L ../../../Core/LeptonTree.h+");
    gSystem->Load("libSmurfFakeLooper.so");
    gSystem->mkdir("plots");

    std::vector<unsigned int> ptThresholds;

    // muons
    if (!defaultOnly) {
        ptThresholds.push_back(0);
        ptThresholds.push_back(5);
        ptThresholds.push_back(10);
        ptThresholds.push_back(20);
        ptThresholds.push_back(25);
        ptThresholds.push_back(30);
        ptThresholds.push_back(50);
        ptThresholds.push_back(80);
    }
    ptThresholds.push_back(15);
    //printFakeRate("met20", "MuonFakeRate_M2", ptThresholds);
    printFakeRate("met20mt15mll_"+era, "MuonFakeRate_M2", "Muon FR with MET20, MT15 and Mll Veto", ptThresholds);

    // electrons
    ptThresholds.clear();

    if (!defaultOnly) {
        ptThresholds.push_back(15);
        ptThresholds.push_back(20);
        ptThresholds.push_back(25);
        ptThresholds.push_back(30);
        ptThresholds.push_back(35);
        ptThresholds.push_back(40);
        ptThresholds.push_back(45);
        ptThresholds.push_back(50);
        ptThresholds.push_back(80);
    }
    ptThresholds.push_back(35);

    //printFakeRate("met20", "ElectronFakeRate_V4", ptThresholds);
    //printFakeRate("met20mt15mll_"+era, "ElectronFakeRate_V4", "Electron FR with MET20, MT15 and Mll Veto", ptThresholds);

}

