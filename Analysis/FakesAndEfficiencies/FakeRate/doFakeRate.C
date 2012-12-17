
#include "Enums.h"

void runFakeLooper(Option option);

void doFakeRate() {

    //
    // load libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../../Core/LeptonTree.h+");
    gSystem->Load("libSmurfFakeLooper.so");

    bool runEle = true;
    bool runMu = true;

    // HCP
    //runFakeLooper(MET20, runEle, runMu, HCP);
    runFakeLooper(MET20MT15MLL, runEle, runMu, HCP);

    // Moriond
    runFakeLooper(MET20MT15MLL, runEle, runMu, MORIOND);

}

void runFakeLooper(Option option, bool runEle, bool runMu, Era era)
{

    //
    // configure looper
    //

    FakeLooper *looper = new FakeLooper(option);
    std::vector<unsigned int> ptThresholds;

    // luminosities in pb
    // see http://www.t2.ucsd.edu/tastwiki/bin/view/Smurf/Moriond2013FRContamination
    TString tag         = "V00-02-09";
    TString eraName     = "HCP";
    if (era == HCP) {
        // HCP good run list
        //tag = "V00-02-08";
        //looper->setMuonLumi(1.337, 14.926);     
        //looper->setElectronLumi(3.028, 15.87);
        tag = "V00-02-09";
        eraName = "HCP";
        looper->setMuonLumi(1.551, 16.592);
        looper->setElectronLumi(3.8553, 18.3181);
    } else if (era == MORIOND) {
        // Moriond good run list
        tag = "V00-02-09";
        eraName = "Moriond";
        looper->setMuonLumi(2.011, 24.787);
        looper->setElectronLumi(4.772, 23.871);
    } else {
        std::cout << "Invalid era " << era << std::endl;
        return;
    }

    // W+jets sample
    TChain *ch_wjets_m = new TChain("leptons");
    ch_wjets_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/FR_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/merged_m.root");
    TChain *ch_wjets_e = new TChain("leptons");
    ch_wjets_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/FR_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/merged_e.root");

    // Drell-Yan sample
    TChain *ch_dy_m = new TChain("leptons");
    ch_dy_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/FR_DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/merged_m.root");
    TChain *ch_dy_e = new TChain("leptons");
    ch_dy_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/FR_DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/merged_e.root");

    // data sample
    TChain *ch_data_m = new TChain("leptons");
    ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012A-13Jul2012-v1_AOD_190456_193621/merged_"+eraName+".root");
    ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/merged_"+eraName+".root");
    ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012B-13Jul2012-v4_AOD_193834_196531/merged_"+eraName+".root");
    ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012C-24Aug2012-v1_AOD_198022_198523/merged_"+eraName+".root");
    ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012C-PromptReco-v2_AOD_198934_203755/merged_"+eraName+".root");
    //ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012D-PromptReco-v1_AOD_203773_208299/merged_"+eraName+".root");
    ch_data_m->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleMu_Run2012D-PromptReco-v1_AOD_203773_208726/merged_"+eraName+".root");

    TChain *ch_data_e = new TChain("leptons");
    ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012A-13Jul2012-v1_AOD_190456_193621/merged_"+eraName+".root");
    ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/merged_"+eraName+".root");
    ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012B-13Jul2012-v1_AOD_193834_196531/merged_"+eraName+".root");
    ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012C-24Aug2012-v1_AOD_198022_198523/merged_"+eraName+".root");
    ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012C-PromptReco-v2_AOD_198934_203755/merged_"+eraName+".root");
    //ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012D-PromptReco-v1_AOD_203773_208299/merged_"+eraName+".root");
    ch_data_e->Add("/smurf/dlevans/LeptonTree/"+tag+"/DoubleElectron_Run2012D-PromptReco-v1_AOD_203773_208726/merged_"+eraName+".root");

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

        looper->Loop(false, ch_wjets_m, "WJets_MuonFakeRate_M2",  LeptonTree::QCDFakeMu,  ptThresholds);
        looper->Loop(false, ch_dy_m,    "DY_MuonFakeRate_M2",     LeptonTree::QCDFakeMu,  ptThresholds);
        looper->Loop(true,  ch_data_m,  "Data_MuonFakeRate_M2",   LeptonTree::QCDFakeMu,  ptThresholds);
    }

    //
    // run electron fake rates
    //

    if (runEle) {
        ptThresholds.clear();
        ptThresholds.push_back(15);
        ptThresholds.push_back(20);
        ptThresholds.push_back(25);
        ptThresholds.push_back(30);
        ptThresholds.push_back(35);
        ptThresholds.push_back(40);
        ptThresholds.push_back(45);
        ptThresholds.push_back(50);
        ptThresholds.push_back(80);

        looper->Loop(false, ch_wjets_e, "WJets_ElectronFakeRate_V4",  LeptonTree::QCDFakeEle,  ptThresholds);
        looper->Loop(false, ch_dy_e,    "DY_ElectronFakeRate_V4",     LeptonTree::QCDFakeEle,  ptThresholds);
        looper->Loop(true,  ch_data_e,  "Data_ElectronFakeRate_V4",   LeptonTree::QCDFakeEle,  ptThresholds);
    }

    //
    // save and manipulate histograms
    //

    TString optName = "met20";
    if (option == MET20MT15) optName = "met20mt15";
    if (option == MET20MT15MLL) optName = "met20mt15mll";

    const TString outFile = Form("histos_FakeLooper_%s_%s.root", optName.Data(), eraName.Data());
    saveHist(outFile.Data());
    deleteHistos();

    //
    // tidy up
    //

    delete looper;
    delete ch_wjets_m;
    delete ch_wjets_e;
    delete ch_dy_m;
    delete ch_dy_e;
    delete ch_data_m;
    delete ch_data_e;

}

