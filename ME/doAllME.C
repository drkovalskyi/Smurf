#include "TFile.h"
#include "TVar.hh"
#include "TString.h"

void doAllME()
{
  gROOT->ProcessLine(".L runME_WWHWW.C+");
  using namespace std;
  
  /* =====   ME arguements  
       1 - input file
       2 - random seed
       3 - smearing level
       4 - number of steps in the integration
       5 - error threshold
       6 - Mass, if 0 - calculate for all processes and Higgs mass hypothesis
  */

  int process=TVar::HWW;
  int seed = 10; // random seed generator
  int SmearLevel = 1; 
  int ncalls = 100000;
  int maxevt = 1000;
  int higgsMass = 0;
  double error = 0.0;

  // location of the smurf trees
  TString inputDir = "/smurf/data/Run2010_Fall10U_SmurfWW_V1/tas/";
  TString outputDir = "./";
  
  // Flags for files to run over 
  // (0 and 1 are easier to modify)

  bool runWW = 0;
  bool runHWW120 = 0;
  bool runHWW130 = 0;
  bool runHWW140 = 0;
  bool runHWW150 = 0;
  bool runHWW160 = 1;
  bool runHWW170 = 0;
  bool runHWW180 = 0;
  bool runHWW190 = 0;
  bool runHWW200 = 0;
  bool runHWW210 = 0;
  bool runHWW220 = 0;
  bool runHWW230 = 0;
  bool runHWW250 = 0;
  bool runHWW300 = 0;


  if(runWW)
    runME_WWHWW(inputDir, "ww.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
  
  if(runHWW120)
    runME_WWHWW(inputDir, "hww120.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
  
 if(runHWW130)
   runME_WWHWW(inputDir, "hww130.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 

 if(runHWW140)
   runME_WWHWW(inputDir, "hww140.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW150)
   runME_WWHWW(inputDir, "hww150.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW160)
   runME_WWHWW(inputDir, "hww160.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW170)
   runME_WWHWW(inputDir, "hww170.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
  
 if(runHWW180)
   runME_WWHWW(inputDir, "hww180.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW190)
   runME_WWHWW(inputDir, "hww190.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW200)
   runME_WWHWW(inputDir, "hww200.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 

 if(runHWW210)
   runME_WWHWW(inputDir, "hww210.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW220)
   runME_WWHWW(inputDir, "hww220.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
  
 if(runHWW230)
   runME_WWHWW(inputDir, "hww230.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW250)
   runME_WWHWW(inputDir, "hww250.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 
 if(runHWW300)
   runME_WWHWW(inputDir, "hww300.root", outputDir, seed, SmearLevel, ncalls, error, higgsMass, maxevt); 
 

}

