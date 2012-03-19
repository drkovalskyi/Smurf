#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "TRint.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCut.h"

//###################
//# main function
//###################

void writesubscript(char *smurfFDir = "/uscms_data/d2/ygao/Run2011_Summer11_SmurfV7_42X/mitf-alljets/WW", int njet = 0) 
{  
  
  // const std::string pwd = "/uscms/home/ygao/HWW2012/CMSSW_4_2_3_Head/src/";
  const std::string pwd = "/uscms/home/ygao/HWW2012/CMSSW_4_2_3_ME_HWW0j/src/";
  const std::string submitFileName = Form("submitAllHWW%ij.sh", njet);
  FILE *text = fopen(submitFileName.c_str(), "w"); 
  fputs("#!/bin/sh\n\n", text);
  
  int nEventsPerJobBkg = 250;
  int nEventsPerJobHiggs = 2000;
  
  writesubscriptsingle(pwd, smurfFDir, "data", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "data-emb-tau123", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "qqww", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "ggww", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "ttbar", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "ttbar_mg", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "tw", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "tw_ds", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "wz", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "zz_py", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "wjets", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "wjets_data", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "wjets_PassFail", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "wgamma", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "wg3l", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "dyee", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "dymm", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "dytt", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "dyee_LooseMET", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "dymm_LooseMET", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "ww_mcnlo", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "ww_mcnlo_up", njet, nEventsPerJobBkg, text);
  writesubscriptsingle(pwd, smurfFDir, "ww_mcnlo_down", njet, nEventsPerJobBkg, text);

  // signal
  writesubscriptsingle(pwd, smurfFDir, "hww115", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww120", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww130", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww140", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww150", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww160", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww170", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww180", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww190", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww200", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww250", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww300", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww350", njet, nEventsPerJobHiggs, text);  
  writesubscriptsingle(pwd, smurfFDir, "hww400", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww450", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww500", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww550", njet, nEventsPerJobHiggs, text);
  writesubscriptsingle(pwd, smurfFDir, "hww600", njet, nEventsPerJobHiggs, text);

}

void writesubscriptsingle(string pwd, string smurfFDir, string fileName, int njet, int nEventsPerJob, FILE *text) 
{

  TFile* fin = new TFile(Form("%s/%ij/%s.root", smurfFDir.c_str(), njet, fileName.c_str()), "READ");
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  fputs(Form("./batchSubmit.sh %s %s %s/%ij/ %i %i WW\n", pwd.c_str(), fileName.c_str(), smurfFDir.c_str(), njet,  ch->GetEntries()/nEventsPerJob+1, nEventsPerJob), text);
  fin->Close();
}  

