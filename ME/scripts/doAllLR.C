#include "TFile.h"
#include "../TVar.hh"
#include "TString.h"

#include "commonFunction.h"

void doAllLR()
{
  gROOT->ProcessLine(".L LR.C+");
  gROOT->ProcessLine(".L makeOverlay.C+");
  using namespace std;
  
  /* =====   ME arguements  
       1 - smurfFDir: original smurf ntuples location
       2 - meFDir: ME input directory
       3 - outputDir: the directory of the re-calculated LR and plots
       
       For each mass hypothesis, one need to specify
       1 - The higgs Process in terms of TVar::Process
       2 - massCut by hand
       3 - The cut-based/BDT-based final results in terms of nsig_CB, nbkg_CB, nsig_BDT, nbkg_BDT
  */

  int maxevt = -1; // run over all events
  
  // location of the smurf trees
  TString smurfFDir = "/Users/yanyan/CMS/SnT/WW/ME/SmurfV3/datasmurf/";
  TString meFDir = "/Users/yanyan/CMS/SnT/WW/ME/SmurfV3/datame/";
  TString outputDir = "../output/";


  // HWW 130
  // make output root files containing the LR for each process, given a certain higgs mass
  LR(smurfFDir, meFDir, "ww_ME.root", outputDir, maxevt, proc_HWW130, 80.0);
  LR(smurfFDir, meFDir, "hww130_ME.root", outputDir, maxevt, proc_HWW130, 80.0);
  LR(smurfFDir, meFDir, "ttbar_ME.root", outputDir, maxevt, proc_HWW130, 80.0);
  LR(smurfFDir, meFDir, "wjets_ME.root", outputDir, maxevt, proc_HWW130, 80.0);
  // ok... now open the resulting root files for each process and overlay the LR distributions
  // and calculate the fom based on S/B
  makeOverlay(TVar::HWW130, smurfFDir, outputDir, "_LR_hww130.root", 80.0, 15.0, 1.9+16.8+43.1, 6.14, 9.79+1.12);


  // HWW 160
  LR(smurfFDir, meFDir, "ww_ME.root", outputDir, maxevt, proc_HWW160, 100.0);
  LR(smurfFDir, meFDir, "hww160_ME.root", outputDir, maxevt, proc_HWW160, 100.0);
  LR(smurfFDir, meFDir, "ttbar_ME.root", outputDir, maxevt, proc_HWW160, 100.0);
  LR(smurfFDir, meFDir, "wjets_ME.root", outputDir, maxevt, proc_HWW160, 100.0);
  makeOverlay(TVar::HWW160, smurfFDir, outputDir, "_LR_hww160.root", 100.0, 30.9, 1.3+17.3, 11.05, 1.79+0.14);


  // HWW 200
  LR(smurfFDir, meFDir, "ww_ME.root", outputDir, maxevt, proc_HWW200, 130.0);
  LR(smurfFDir, meFDir, "hww200_ME.root", outputDir, maxevt, proc_HWW200, 130.0);
  LR(smurfFDir, meFDir, "ttbar_ME.root", outputDir, maxevt, proc_HWW200, 130.0);
  LR(smurfFDir, meFDir, "wjets_ME.root", outputDir, maxevt, proc_HWW200, 130.0);
  makeOverlay(TVar::HWW200, smurfFDir, outputDir, "_LR_hww200.root", 130.0, 0.0, 0.0, 8.05, 8.32+0.94);
  
  
  
}

