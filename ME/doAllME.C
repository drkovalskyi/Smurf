#include "TFile.h"
#include "TVar.hh"
#include "TString.h"

void doAllME()
{
  gROOT->ProcessLine(".L runME_test.C++");
  using namespace std;
  
  /* =====   ME arguements  
       1 - input file
       2 - random seeds
       3 - smearing level
       4 - number of steps in the integration
       5 - error threshold
       6 - Mass, if 0 - calculate for all processes and Higgs mass hypothesis
  */

  int process=TVar::HWW;
  int seed = 10; // random seed generator
  int SmearLevel = 1; 
  int ncalls = 100000;
  int maxevt = 10;
  int higgsMass = 0;
  double error = 0.0;

  // location of the smurf trees
  // TString smurfFDir = "/smurf/data/Run2011_Spring11_SmurfV3/tas-zerojet/";
  TString smurfFDir = "/smurf/yygao/data/Run2011_Spring11_SmurfV5/mitf-zerojet/";
  TString meFDir = "datame/";

  // Flags for files to run over 
  // (0 and 1 are easier to modify)

  bool runData = 0;
  bool runData_lfake = 0;
  bool runMC_lfake = ;
  bool runqqWW = 0;
  bool runggWW = 0;
  bool runWjets = 0;
  bool runWgamma = 0;
  bool runttbar = 0;
  bool runWZ = 0;
  bool runZZ = 0;
  bool runDYee = 0;
  bool runDYmm = 0;
  bool runDYtt = 0;
  bool runtW = 0;
  bool runHWW115 = 0;
  bool runHWW120 = 1;
  bool runHWW130 = 0;
  bool runHWW140 = 0;
  bool runHWW150 = 0;
  bool runHWW160 = 0;
  bool runHWW170 = 0;
  bool runHWW180 = 0;
  bool runHWW190 = 0;
  bool runHWW200 = 0;
  bool runHWW210 = 0;
  bool runHWW220 = 0;
  bool runHWW230 = 0;
  bool runHWW250 = 0;
  bool runHWW300 = 0;
  
  
  if(runData_lfake)
    runME_test(smurfFDir, "data_lfake.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runMC_lfake)
    runME_test(smurfFDir, "mc_l_fake_pu.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runData)
    runME_test(smurfFDir, "data.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);
  
  if(runqqWW)
    runME_test(smurfFDir, "qqww.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runggWW)
    runME_test(smurfFDir, "ggww.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runZZ)
    runME_test("/smurf/data/Run2011_Spring11_SmurfV3/mitf-zerojets-wwselection/", "zz.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runWZ)
    runME_test("/smurf/data/Run2011_Spring11_SmurfV3/mitf-zerojets-wwselection/", "wz.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);
  
  if(runDYee)
    runME_test(smurfFDir, "dyee.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runDYmm)
    runME_test(smurfFDir, "dymm.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runDYtt)
    runME_test(smurfFDir, "dytt.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);
  
  if(runWjets)
    //runME_test(smurfFDir, "wjets_pythia.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);
    runME_test(smurfFDir, "wjets.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runWgamma)
    runME_test(smurfFDir, "wgamma.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runttbar)
    runME_test(smurfFDir, "ttbar.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);
  
  if(runtW)
    runME_test(smurfFDir, "tw.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runHWW115)
    runME_test(smurfFDir, "hww115.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runHWW120)
    runME_test(smurfFDir, "hww120.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);

  if(runHWW130)
    runME_test(smurfFDir, "hww130.root", meFDir, seed, SmearLevel, ncalls, error, maxevt, 0, TVar::DEBUG);
  
}

