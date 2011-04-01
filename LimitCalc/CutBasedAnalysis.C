#include "../Core/SmurfTree.h"
#include <string>
#include "TH1F.h"
#include <fstream>
typedef std::vector<std::pair<std::string,SmurfTree::DataType> > Samples;
typedef std::pair<std::string,SmurfTree::DataType> Sample;

//
// If you want to change things, look at the following parts:
//  gScale - set lumi 
//  splitByChannels - 1 or 4 channels
// Event selection can be specified in PassedFinalSelection():
//  HWWCuts_Run2010Paper - reference 2010 HWW cut based selection
//  HWWCuts_SmurfV1 - similar to HWWCuts_Run2010Paper, but low pt for M<=160 is 10  

// Global Config Parameters
const double gScale = 1.0; // fb^-1
const bool splitByChannels = false; // fb^-1

struct LimitInfo{
  double observed;
  double exp_m2sig;
  double exp_m1sig;
  double exp_mean;
  double exp_p1sig;
  double exp_p2sig;
};

struct SigInfo{
  double observed;
  double exp_m2sig;
  double exp_m1sig;
  double exp_mean;
  double exp_p1sig;
  double exp_p2sig;
};

// declaration 
bool HWWCuts_Run2010Paper(SmurfTree& tree, SmurfTree::DataType sig_type);
bool HWWCuts_SmurfV1(SmurfTree& tree, SmurfTree::DataType sig_type);
bool HWWCuts_test(SmurfTree& tree, SmurfTree::DataType sig_type);
void MakeLandsCard(const Sample& sig, const Samples& bkgs, std::string file_name);

// definition
bool PassedFinalSelection(SmurfTree& tree, SmurfTree::DataType sig_type){
  return HWWCuts_SmurfV1(tree,sig_type);
  // HWWCuts_Run2010Paper(tree,sig_type)
}

void CutBasedAnalysis(std::string path = "data/")
{
  Samples sigSamples;
  Samples bkgSamples;
  bkgSamples.push_back(Sample(path+"ww.root",     SmurfTree::qqww));
  bkgSamples.push_back(Sample(path+"ggww.root",   SmurfTree::ggww));
  bkgSamples.push_back(Sample(path+"wjets.root",  SmurfTree::wjets));
  bkgSamples.push_back(Sample(path+"ttbar.root",  SmurfTree::ttbar));
  bkgSamples.push_back(Sample(path+"tw.root",     SmurfTree::tw));
  bkgSamples.push_back(Sample(path+"dyee.root",   SmurfTree::dyee));
  bkgSamples.push_back(Sample(path+"dymm.root",   SmurfTree::dymm));
  bkgSamples.push_back(Sample(path+"dytt.root",   SmurfTree::dytt));
  bkgSamples.push_back(Sample(path+"wz.root",     SmurfTree::wz));
  bkgSamples.push_back(Sample(path+"zz.root",     SmurfTree::zz));

  sigSamples.push_back(Sample(path+"hww130.root", SmurfTree::hww130));
  sigSamples.push_back(Sample(path+"hww160.root", SmurfTree::hww160));
  sigSamples.push_back(Sample(path+"hww200.root", SmurfTree::hww200));
  sigSamples.push_back(Sample(path+"hww250.root", SmurfTree::hww250));
  
  for(Samples::const_iterator sig=sigSamples.begin(); sig!=sigSamples.end(); ++sig)
    {
      MakeLandsCard(*sig, bkgSamples, SmurfTree::name(sig->second)+".card");
    }
}

void AddEvents(const Sample& sample, TH1& hist, SmurfTree::DataType sig_type){
  SmurfTree tree;
  tree.LoadTree(sample.first.c_str());
  tree.InitTree();
  Long64_t nentries = tree.tree_->GetEntries();
  for (Long64_t i = 0; i < nentries; i++){
    tree.tree_->GetEntry(i);
    if ( PassedFinalSelection(tree,sig_type) )
      hist.Fill(tree.type_+0.5,tree.scale1fb_*gScale);
  }
}

void MakeLandsCard(const Sample& sig, const Samples& bkgs, std::string file_name){
  TH1F hsig("hsig","hsig",4,0,4);
  hsig.SetDirectory(0);
  hsig.Sumw2();
  AddEvents(sig,hsig,sig.second);
  TH1F hbkg("hbkg","hbkg",4,0,4);
  hbkg.SetDirectory(0);
  hbkg.Sumw2();
  for(Samples::const_iterator bkg=bkgs.begin(); bkg!=bkgs.end(); ++bkg)
    AddEvents(*bkg,hbkg,sig.second);
  std::ofstream fout(file_name.c_str());
  fout << "Signal target: " << SmurfTree::name(sig.second) << "\n";
  fout << "jmax   1  number of backgrounds\n"; 
  fout << "kmax   1  number of nuisance parameters\n";
  if ( splitByChannels ) {
    fout << "imax   4  number of channels\n";
    fout << "Observation  0 0 0 0\n";
  } else {
    fout << "imax   1  number of channels\n";
    fout << "Observation  0\n";
  }
  if ( splitByChannels ){
    fout << "bin        1     1    2    2    3    3    4    4\n";
    fout << "process    0     1    0    1    0    1    0    1\n";
    fout << "rate " << hsig.GetBinContent(1) << " " << hbkg.GetBinContent(1) 
	 << " " << hsig.GetBinContent(2) << " " << hbkg.GetBinContent(2)
	 << " " << hsig.GetBinContent(3) << " " << hbkg.GetBinContent(3)  
	 << " " << hsig.GetBinContent(4) << " " << hbkg.GetBinContent(4) 
	 << "\n";
    fout << "1  lnN     1.1   1.1  1.1  1.1  1.1  1.1  1.1  1.1\n";
  } else {
    fout << "bin        1     1\n";
    fout << "process    0     1\n";
    fout << "rate " << hsig.Integral() << " " << hbkg.Integral() << "\n";
    fout << "1  lnN  1.1  1.1\n";
  }
  fout.close();
  std::cout << "Card " << file_name << " is done " << std::endl; 
}

bool HWWCuts_Run2010Paper(SmurfTree& tree, SmurfTree::DataType sig_type){
  switch (sig_type){
  case SmurfTree::hww120: 
    return tree.lep1_.pt()>20 && tree.lep2_.pt()>20 &&
      tree.dilep_.mass()<40 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww130: 
  case SmurfTree::hww140:
    return tree.lep1_.pt()>25 && tree.lep2_.pt()>20 &&
      tree.dilep_.mass()<45 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww150:
    return tree.lep1_.pt()>27 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww160:
    return tree.lep1_.pt()>30 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww170:
    return tree.lep1_.pt()>34 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww180:
    return tree.lep1_.pt()>36 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<60 && fabs(tree.dPhi_)<M_PI/180*70;
  case SmurfTree::hww190:
    return tree.lep1_.pt()>38 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<80 && fabs(tree.dPhi_)<M_PI/180*90;
  case SmurfTree::hww200:
    return tree.lep1_.pt()>40 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<90 && fabs(tree.dPhi_)<M_PI/180*100;
  case SmurfTree::hww210:
    return tree.lep1_.pt()>44 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<110 && fabs(tree.dPhi_)<M_PI/180*110;
  case SmurfTree::hww220:
    return tree.lep1_.pt()>48 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<120 && fabs(tree.dPhi_)<M_PI/180*120;
  case SmurfTree::hww230:
    return tree.lep1_.pt()>52 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<130 && fabs(tree.dPhi_)<M_PI/180*130;
  case SmurfTree::hww250:
    return tree.lep1_.pt()>55 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<150 && fabs(tree.dPhi_)<M_PI/180*140;
  default: return false;
  }
}

bool HWWCuts_SmurfV1(SmurfTree& tree, SmurfTree::DataType sig_type){
  switch (sig_type){
  case SmurfTree::hww120: 
    return tree.lep1_.pt()>20 && tree.lep2_.pt()>10 &&
      tree.dilep_.mass()<40 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww130: 
  case SmurfTree::hww140:
    return tree.lep1_.pt()>25 && tree.lep2_.pt()>10 &&
      tree.dilep_.mass()<45 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww150:
    return tree.lep1_.pt()>27 && tree.lep2_.pt()>10 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww160:
    return tree.lep1_.pt()>30 && tree.lep2_.pt()>10 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww170:
    return tree.lep1_.pt()>34 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww180:
    return tree.lep1_.pt()>36 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<60 && fabs(tree.dPhi_)<M_PI/180*70;
  case SmurfTree::hww190:
    return tree.lep1_.pt()>38 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<80 && fabs(tree.dPhi_)<M_PI/180*90;
  case SmurfTree::hww200:
    return tree.lep1_.pt()>40 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<90 && fabs(tree.dPhi_)<M_PI/180*100;
  case SmurfTree::hww210:
    return tree.lep1_.pt()>44 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<110 && fabs(tree.dPhi_)<M_PI/180*110;
  case SmurfTree::hww220:
    return tree.lep1_.pt()>48 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<120 && fabs(tree.dPhi_)<M_PI/180*120;
  case SmurfTree::hww230:
    return tree.lep1_.pt()>52 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<130 && fabs(tree.dPhi_)<M_PI/180*130;
  case SmurfTree::hww250:
    return tree.lep1_.pt()>55 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<150 && fabs(tree.dPhi_)<M_PI/180*140;
  default: return false;
  }
}

bool HWWCuts_test(SmurfTree& tree, SmurfTree::DataType sig_type){
  switch (sig_type){
  case SmurfTree::hww120: 
    return tree.lep1_.pt()>20 &&
      tree.dilep_.mass()<40 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww130: 
  case SmurfTree::hww140:
    // return tree.dPhi_<2 && tree.met_<20+(2-tree.dPhi_)/2*50;
    //    return tree.lep1_.pt()>25 &&
    //      tree.dilep_.mass()<45 && fabs(tree.dPhi_)<M_PI/180*60;
    return tree.dilep_.mass()<50&&tree.lep1_.pt()<50&&tree.met_<60;
  case SmurfTree::hww150:
    return tree.lep1_.pt()>27 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww160:
    return tree.lep1_.pt()>30 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww170:
    return tree.lep1_.pt()>34 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww180:
    return tree.lep1_.pt()>36 &&
      tree.dilep_.mass()<60 && fabs(tree.dPhi_)<M_PI/180*70;
  case SmurfTree::hww190:
    return tree.lep1_.pt()>38 &&
      tree.dilep_.mass()<80 && fabs(tree.dPhi_)<M_PI/180*90;
  case SmurfTree::hww200:
    return tree.lep1_.pt()>40 &&
      tree.dilep_.mass()<90 && fabs(tree.dPhi_)<M_PI/180*100;
  case SmurfTree::hww210:
    return tree.lep1_.pt()>44 &&
      tree.dilep_.mass()<110 && fabs(tree.dPhi_)<M_PI/180*110;
  case SmurfTree::hww220:
    return tree.lep1_.pt()>48 &&
      tree.dilep_.mass()<120 && fabs(tree.dPhi_)<M_PI/180*120;
  case SmurfTree::hww230:
    return tree.lep1_.pt()>52 &&
      tree.dilep_.mass()<130 && fabs(tree.dPhi_)<M_PI/180*130;
  case SmurfTree::hww250:
    return tree.lep1_.pt()>55 &&
      tree.dilep_.mass()<150 && fabs(tree.dPhi_)<M_PI/180*140;
  default: return false;
  }
}



//
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationInput
//
// First lines are solely for bookkeeping. Here, the user can put any
// information (analysis version, Higgs mass point, date, anything)

// *imax* specifies how many channels (bins) will be combined; there is
// no limit, e.g. a binned histogram can be unfolded in such a table
// (imax = the number of bins)

// *jmax* specifies how many backgrounds will be individually
// accounted for

// *kmax* specifies how many independent sources of systematics will
// be provided in the table

// *Observation* Here one lists the event counts observed in each
// channel (there should be imax numbers here)

// bin indicates which channel that the column of data below belongs to

//   processes This first process row helps one visually identify background by abbreviated name (not used in the analysis itself)

//   processes The second process row enumerates separately tracked processes from 0 (signal) to jmax (backgrounds). Overall, there are imax x (jmax +1) columns of data

//   rate specifies the expected signal and background rates for each channel (expectations derived from Monte Carlo or data-driven measurements)

// Each next line corresponds to one independent source of systematic errors that may affect signal and backgrounds event yields.

// First entry in the row enumerates a source of systematic errors
// 							    Second entry specifies the type of systematic error pdf:
// 							    lnN stands for lognormal (recommended option)
// 							    trG --- truncated Gaussian (this option currently works with LandS only; we never bothered to transport this option into RooStats)
// 							    other options will be added later (gamma distribution, flat, etc.)
// 							    Next imax x (jmax +1) entries specify the impact of this source of uncertainties on each individual process. Here, a user should list parameters describing width of the systematic error pdf. All errors in one row are assumed to be 100% correlated, but do not have to be of the same scale. E.g., muon efficiency uncertainty of 1% would imply 4% error on H→ZZ→4μ, 2% error on H→ZZ→2e2μ, and 0% error on H→ZZ→4e; (κ=1.04, 1.02, 1.00 in the context of the lognormal pdf). Entering positive and negative numbers will imply that the errors between corresponding bins have a 100% anti -correlation (this allows one to introduce systematics errors in shapes of signal and background distributions when channels represent bins of some observable, e.g. MVA output).
// 							    At the end of the line one can add any text (name the systematic error source, remind what drives its scale, etc.)
// 							    Each next line represent a new independent (i.e. uncorrelated) source of errors.
