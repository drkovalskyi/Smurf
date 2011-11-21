#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->Load("lands.so");
  gSystem->CompileMacro("PlotExpectedLimits.C","k");
  // PlotExpectedLimits("output/hww-1000pb_0-1_cut*","H #rightarrow WW #rightarrow 2l2#nu + 0/1 jets");
  // plotBand();
}
#endif

#include "TStyle.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include "TMath.h"
#include "TROOT.h"
#include "TPRegexp.h"

#ifndef __CINT__
#include "../../LandS/include/PlotUtilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
// #include "../../LandS/include/CRandom.h"
// #include "../../LandS/include/CountingModel.h"
// #include "../../LandS/include/BayesianBase.h"
// #include "../../LandS/include/LimitBands.h"
// #include "../../LandS/include/UtilsROOT.h"

struct LimitInfo{
  double exp_m2sig;
  double exp_m1sig;
  double exp_mean;
  double exp_median;
  double exp_p1sig;
  double exp_p2sig;
  double mass;
  double observed;
};
// void AddLimits(std::vector<LimitInfo>& limits, const char* card, double mass){
//   if (gSystem->AccessPathName(card)){
//     std::cout << "Card " << card << " doesn't exist, skip it" << std::endl;
//     return;
//   }
//   int seed = 1234;
//   int noutcomes = 1000;
//   int nexps = 20000;
//   CRandom *rdm = new CRandom(seed);  //initilize a random generator
//   // set up couting model
//   CountingModel *cms=new CountingModel();
//   cms->SetRdm(rdm);
//   ConfigureModel(cms, card); 
//   cms->SetUseSystematicErrors(true);
//   cms->UseAsimovData(-1);
//   // cms->Print(0);
//   cms->RemoveChannelsWithExpectedSignal0orBkg0();
//   // cms->Print(0);

//   lands::BayesianBase bys(cms, 0.05, 1.e-3);
//   bys.SetNumToys(nexps);
//   bys.SetDebug(0);
//   // if(XSprior==20)
//   // bys.SetCrossSectionPrior(corr);
//   // else 
//   bys.SetCrossSectionPrior(flat);
//   // double rtmp;
//   // rtmp = bys.Limit();
//   // cout<<"------------------------------------------------------------"<<endl;
//   // cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
//   // cout<<"------------------------------------------------------------"<<endl;

//   cms->SetSignalScaleFactor(1.);
  
//   lands::LimitBands lb(&bys, cms);	
//   lb.SetDebug(0);
//   lb.BysLimitBands(0.05, noutcomes, nexps);
  
//   LimitInfo limit;
//   limit.mass = mass;
//   limit.exp_m2sig = lb.GetBysLimit(-2);
//   limit.exp_m1sig = lb.GetBysLimit(-1);
//   limit.exp_mean  = lb.GetBysLimitMean();
//   limit.exp_p1sig = lb.GetBysLimit(+1);
//   limit.exp_p2sig = lb.GetBysLimit(+2);

//   limits.push_back(limit);
// }
double getObservedLimit(const char* file)
{
  if (! file) return -1;
  // cout << file << endl;
  string line;
  ifstream myfile (file);
  double observed(-1);
  if (myfile.is_open())
    {
      while ( myfile.good() )
	{
	  getline (myfile,line);
	  TPRegexp pattern("Observed.*?= (.+)$");
	  // ("Observed.*?=\ *(.+)$");
	  TObjArray* matches = pattern.MatchS(line.c_str());
	  if (matches && matches->GetLast()==1){
	    TString result(((TObjString*)matches->At(1))->GetString());
	    observed = atof(result.Data());
	  }
	}
      myfile.close();   
    }
  return observed;
}

void AddLimits_Bayesian(std::vector<LimitInfo>& limits, const char* dir, const char* name, double mass)
{
  TChain* chain = new TChain("T");
  chain->Add(Form("%s/output/%s-%.0f-bayesian*_limits_tree.root",dir,name,mass));
  double result(0);
  chain->SetBranchAddress("brT", &result);
  Long64_t nentries = chain->GetEntries();
  if (nentries<=0) return;
  std::vector<double> values; 
  double sum(0);
  for (Long64_t i = 0; i < nentries; i++){
    chain->GetEvent(i);
    if(fabs(result)<1E3){
      values.push_back(result);
      sum+=result;
    }
  }
  std::sort(values.begin(),values.end());
  const double trialFactor = 1;
  const double prob1S = TMath::Erfc(1/sqrt(2))/2/trialFactor;
  const double prob2S = TMath::Erfc(2/sqrt(2))/2/trialFactor;
  std::cout << mass << std::endl;
  std::cout << "\tN toys: " << nentries << std::endl;
  std::cout << "\t-2 sigma: \t" << values.at(int(prob2S*values.size())) << std::endl;
  std::cout << "\t-1 sigma: \t" << values.at(int(prob1S*values.size())) << std::endl;
  std::cout << "\tmedian: \t"   << values.at(int(values.size()/2)) << std::endl;
  std::cout << "\tmean: \t"     << sum/values.size() << std::endl;
  std::cout << "\t+1 sigma: \t" << values.at(values.size()-int(prob1S*values.size())) << std::endl;
  std::cout << "\t+2 sigma: \t" << values.at(values.size()-int(prob2S*values.size())) << std::endl;

  LimitInfo limit;
  limit.mass = mass;
  limit.exp_m2sig  = values.at(int(prob2S*values.size()));
  limit.exp_m1sig  = values.at(int(prob1S*values.size()));
  limit.exp_mean   = sum/values.size() ;
  limit.exp_median = values.at(int(values.size()/2));
  limit.exp_p1sig  = values.at(values.size()-int(prob1S*values.size()));
  limit.exp_p2sig  = values.at(values.size()-int(prob2S*values.size()));

  limit.observed = getObservedLimit(Form("%s/output/%s-%.0f-bayesian.observed",dir,name,mass));
  printf("observed: %f\n",limit.observed);
  limits.push_back(limit);
}

void AddLimits_CLs_Asymptotic(std::vector<LimitInfo>& limits, const char* dir, const char* name, double mass)
{
  LimitInfo limit;
  limit.mass = mass;
  limit.observed = getObservedLimit(Form("%s/output/%s-%.0f-cls_asym.observed",dir,name,mass));
  printf("observed: %f\n",limit.observed);

  const char* file = Form("%s/output/%s-%.0f-cls_asym.log",dir,name,mass);
  string line;
  ifstream myfile (file);
  if (myfile.is_open())
    {
      while ( myfile.good() )
        {
          getline (myfile,line);
          TPRegexp pattern("BANDS");
          if (pattern.Match(line.c_str())){
	    // std::cout << line << std::endl;
	    stringstream ss;
	    ss << line;
	    std::string s;
	    ss >> s;
	    ss >> limit.exp_m2sig >> limit.exp_m1sig >> limit.exp_mean >> limit.exp_p1sig >> limit.exp_p2sig >> limit.exp_median;
	    // std::cout << r2m << " : " << r1m << " : " << r1p << " : " << r2p << std::endl;
	    // break;
            // TString result(((TObjString*)matches->At(1))->GetString());
            // myfile.close();
            // return atof(result.Data());
          }
        }
      myfile.close();
    }
  limits.push_back(limit);
}


void PlotExpectedLimits(std::vector<LimitInfo>& limits, const char* title){
  if (limits.empty()){
    std::cout << "No limits to plot." << std::endl;
    return;
  }
  double *limits_m1s = new double[limits.size()];
  double *limits_m2s = new double[limits.size()];
  double *limits_p1s = new double[limits.size()];
  double *limits_p2s = new double[limits.size()];
  double *limits_mean = new double[limits.size()];
  double *mass_points = new double[limits.size()];
  double *observed = new double[limits.size()];
  bool showObserved = true;
  for(unsigned int i=0; i<limits.size(); ++i){
    limits_m2s[i]  = limits.at(i).exp_m2sig;
    limits_m1s[i]  = limits.at(i).exp_m1sig;
    // limits_mean[i] = limits.at(i).exp_mean;
    limits_mean[i] = limits.at(i).exp_median;
    limits_p1s[i]  = limits.at(i).exp_p1sig;
    limits_p2s[i]  = limits.at(i).exp_p2sig;
    mass_points[i] = limits.at(i).mass;
    observed[i] = limits.at(i).observed;
    if (limits.at(i).observed<0) showObserved=false;
  }
  gROOT->SetStyle("Plain");
  gStyle->SetCanvasDefH      (600);
  gStyle->SetCanvasDefW      (800);
  gStyle->SetTitleOffset( 1, "y");
  TPaveText *pt = lands::SetTPaveText(0.5, 0.95, 0.8, 0.95); //SetTPaveText(0.5, 0.95, 0.8, 0.95)
  pt->AddText(title);
  lands::PlotWithBelts* lb = 0;
  float yMax = 10; 
  if (showObserved){
    lb = new lands::PlotWithBelts(limits_m1s, limits_p1s, limits_m2s, limits_p2s,
				  limits_mean, observed, limits.size(), mass_points, 
				  "limits", pt,
				  // 100, 305, 0, yMax, false, 
				  100, 605, 0, yMax, false, 
				  // 100, 205, 0, yMax, false, 
				  ";Higgs mass, m_{H} [GeV/c^{2}]; 95% CL Limit on #sigma/#sigma_{SM}; observed ");
    lb->plot();
    lb->drawLegend("95% CL exclusion: median","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "observed");
  } else {
    lb = new lands::PlotWithBelts(limits_m1s, limits_p1s, limits_m2s, limits_p2s,
				  limits_mean, limits.size(), mass_points, 
				  "limits", pt,
				  100, 305, 0, yMax, false, 
				  // 100, 205, 0, yMax, false, 
				  ";Higgs mass, m_{H} [GeV/c^{2}]; 95% CL Limit on #sigma/#sigma_{SM} ");
    lb->plot();
    lb->drawLegend("95% CL exclusion: median","95% CL exclusion: 68% band", "95% CL exclusion: 95% band");
  }
  // lb->drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band");
  lands::MoveLegend(lb->getLegend(),0.45,0.7);
  // lands::MoveLegend(lb->getLegend(),0.18,0.7);
  lb->getMeanGraph()->SetLineWidth(2);
  lb->getMeanGraph()->SetMarkerStyle(20);
  /*
  lb->getMeanGraph()->SetLineStyle(1);
  lb->getMeanGraph()->SetLineColor(kBlue);
  lb->getMeanGraph()->SetLineWidth(2);
  // 	lb->getObsGraph()->SetLineStyle(2);
  // 	lb->getObsGraph()->SetLineWidth(1);
  lb->getLine()->SetLineColor(kRed);
  */
  lb->save();
}

void ReportLimits(std::vector<LimitInfo>& limits, const char* title, const char* dir, const char* name){
  if (limits.empty()){
    std::cout << "No limits to plot." << std::endl;
    return;
  }
//   bool showObserved = true;
//   for(unsigned int i=0; i<limits.size(); ++i){
//     if (limits.at(i).observed<0) showObserved=false;
//   }
  printf("Mass \tObserved \tMedian Expected \t68%% C.L. band \t95%% C.L. band\n");
  for(unsigned int i=0; i<limits.size(); ++i){
    printf("%.0f \t%0.1f \t%0.1f \t[%0.1f, %0.1f] \t[%0.1f, %0.1f]\n",
	   limits.at(i).mass,
	   limits.at(i).observed, limits.at(i).exp_median, 
	   limits.at(i).exp_m1sig, limits.at(i).exp_p1sig,
	   limits.at(i).exp_m2sig, limits.at(i).exp_p2sig);
  }
  // print latex format
  for(unsigned int i=0; i<limits.size(); ++i){
    printf("%.0f & %0.1f & %0.1f & [%0.1f, %0.1f] & [%0.1f, %0.1f] \\\\\n",
	   limits.at(i).mass,
	   limits.at(i).observed, limits.at(i).exp_median, 
	   limits.at(i).exp_m1sig, limits.at(i).exp_p1sig,
	   limits.at(i).exp_m2sig, limits.at(i).exp_p2sig);
  }
}

#endif


void PlotExpectedLimits(const char* file, const char* title)
{
  std::vector<LimitInfo> limits;
  std::ifstream fin(file);
  while(!fin.eof()){
    LimitInfo limit;
    fin >> limit.mass >> limit.observed >> limit.exp_m2sig >> limit.exp_m1sig
	>> limit.exp_median >> limit.exp_p1sig >> limit.exp_p2sig;
    std::cout << limit.mass << "\t" << limit.observed << "\t" << limit.exp_m2sig << "\t" 
	      << limit.exp_m1sig << "\t" << limit.exp_median << "\t" 
	      << limit.exp_p1sig << "\t" << limit.exp_p2sig << std::endl;
    limits.push_back(limit);
  }
  PlotExpectedLimits(limits,title);
}

void PlotExpectedLimits(const char* dir, const char* name, 
			const char* title, const char* type=0)
{
  std::vector<LimitInfo> limits;
  const unsigned int npoints = 19;
  unsigned int points[npoints] = 
    {110,115,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  for(unsigned i=0; i<npoints; ++i){
    if (type==0 || strcmp(type,"Bayesian")==0)
      AddLimits_Bayesian(limits, dir, name, points[i]);
    if (strcmp(type,"CLs-asymptotic")==0)
      AddLimits_CLs_Asymptotic(limits, dir, name, points[i]);
  }
  PlotExpectedLimits(limits,title);
  ReportLimits(limits,title,dir,name);
  /*

  TCanvas* c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  //   c1->SetFillColor(42);
  //   c1->SetGrid();
  //   c1->GetFrame()->SetFillColor(21);
  //   c1->GetFrame()->SetBorderSize(12);
  Int_t n = 0;
  //   Double_t x[10]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
  //   Double_t y[10]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
  //   Double_t ex[10] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
  //   Double_t ey[10] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
  if (
  TGraphErrors* gr = new TGraphErrors(n,x,y,ex,ey);
  gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("ALP");
  */
}


void PlotExpectedLimits()
{
  std::vector<LimitInfo> limits;
//   AddLimits(limits, "hww120.card", 120);
//   AddLimits(limits, "hww130.card", 130);
//   AddLimits(limits, "hww140.card", 140);
//   AddLimits(limits, "hww150.card", 150);
//   AddLimits(limits, "hww160.card", 160);
//   AddLimits(limits, "hww170.card", 170);
//   AddLimits(limits, "hww180.card", 180);
//   AddLimits(limits, "hww190.card", 190);
//   AddLimits(limits, "hww200.card", 200);
//   AddLimits(limits, "hww250.card", 250);
  
  PlotExpectedLimits(limits,"");
  /*

  TCanvas* c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  //   c1->SetFillColor(42);
  //   c1->SetGrid();
  //   c1->GetFrame()->SetFillColor(21);
  //   c1->GetFrame()->SetBorderSize(12);
  Int_t n = 0;
  //   Double_t x[10]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
  //   Double_t y[10]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
  //   Double_t ex[10] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
  //   Double_t ey[10] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
  if (
  TGraphErrors* gr = new TGraphErrors(n,x,y,ex,ey);
  gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("ALP");
  */
}

void plotBand()
{
  // gSystem->Load("PlotUtilities_standalone_cc.so");
	gStyle->SetCanvasDefH      (600);
	gStyle->SetCanvasDefW      (800);
	gStyle->SetTitleOffset( 1, "y");

	const int npoints = 5;
	double limits_obs[npoints] = {
	  41.74677 ,
	  24.58797 ,
	  24.73013 ,
	  58.37274 ,
	 150.18342 , 
	};
	double limits_1fb_obs[npoints] = {
	  3.19987 ,
	  1.86819 ,
	  1.57866 ,
	  3.04254 ,
	  6.74371 , 
	};
	double bands_1fb[npoints][5] = {
		// -2 sigma, -1 sigma,   mean,     +1 sigma,  +2 sigma
		{ 2.36292  ,  3.10099   , 4.99050   , 6.50240  ,  8.84537},	
		{ 1.34476  ,  1.81795   , 2.92033   , 3.77545  ,  5.19139},
		{ 1.10501  ,  1.33298   , 2.18627   , 2.68424  ,  3.99269},
		{ 2.47239  ,  2.47239   , 3.85914   , 4.53227  ,  6.67285},
		{ 6.21168  ,  6.21168   , 8.17544   , 9.06174  , 12.84551},
	};
	double bands[npoints][5] = {
		// -2 sigma, -1 sigma,   mean,     +1 sigma,  +2 sigma
		{ 43.13062 ,  43.13062 ,  48.80107  , 46.82769 ,  65.03096},
		{ 25.52106 ,  25.52106 ,  28.82003  , 27.57525 ,  38.26959},
		{ 26.90173 ,  26.90173 ,  28.75666  , 26.90173 ,  35.56380},
		{ 65.68063 ,  65.68063 ,  67.36003  , 65.68063 ,  72.34968},
		{170.76523 , 170.76523 , 172.93633  ,170.76523 , 170.76523},
	};
	double xpoints[npoints] = {250, 300, 400, 500, 600};

	bool projectingRLimitLogY = true;
	double projectingXmin = 190, projectingXmax = 610;
	double projectingRLimitYmin = 0.5, projectingRLimitYmax = 200;
	string projectingRLimitXYtitles = ";Higgs mass, m_{H} [GeV/c^{2}]; 95% CL Limit on #sigma/#sigma_{SM} ";
	string ssave= "./limits_35pb";

	double limits_m1s[npoints], limits_m2s[npoints], limits_p1s[npoints], limits_p2s[npoints], limits_mean[npoints];
	double limits_1fb_m1s[npoints], limits_1fb_m2s[npoints], limits_1fb_p1s[npoints], limits_1fb_p2s[npoints], limits_1fb_mean[npoints];
	for(int n=0; n<npoints; n++){
		limits_m2s[n]=bands[n][0];
		limits_m1s[n]=bands[n][1];
		limits_mean[n]=bands[n][2];
		limits_p1s[n]=bands[n][3];
		limits_p2s[n]=bands[n][4];
		limits_1fb_m2s[n]=bands_1fb[n][0];
		limits_1fb_m1s[n]=bands_1fb[n][1];
		limits_1fb_mean[n]=bands_1fb[n][2];
		limits_1fb_p1s[n]=bands_1fb[n][3];
		limits_1fb_p2s[n]=bands_1fb[n][4];
	}

	TPaveText *pt;
	pt = lands::SetTPaveText(0.3, 0.7, 0.5, 0.9);
	pt->AddText("HZZ #rightarrow 2l2#nu");
	lands::PlotWithBelts lb(
			limits_m1s, limits_p1s, limits_m2s, limits_p2s,
			limits_mean, npoints,  
			xpoints, ssave+"_1", pt,
			projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
	lb.plot();
	lb.drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band");
	//lands::MoveLegend(lb.getLegend(),0.45,0.7);
	lands::MoveLegend(lb.getLegend(),0.3,0.7);
	lb.getMeanGraph()->SetLineStyle(1);
	lb.getMeanGraph()->SetLineColor(kBlue);
	lb.getMeanGraph()->SetLineWidth(2);
// 	lb.getObsGraph()->SetLineStyle(2);
// 	lb.getObsGraph()->SetLineWidth(1);
	lb.getLine()->SetLineColor(kRed);
// 	lands::PlotWithBelts lb2(
// 			limits_1fb_m1s, limits_1fb_p1s, limits_1fb_m2s, limits_1fb_p2s,
// 			limits_1fb_mean, limits_1fb_obs, npoints,  
// 			xpoints, ssave+"_2", pt,
// 			projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
// 	lb2.plot();
// 	lb2.getObsGraph()->SetMarkerStyle(1);//remove marker
// 	lb2.drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: mean(nosys)");
// 	lb2.getMeanGraph()->SetLineStyle(1);
// 	lb2.getMeanGraph()->SetLineColor(kBlue);
// 	lb2.getMeanGraph()->SetLineWidth(2);
// 	lb2.getObsGraph()->SetLineStyle(2);
// 	lb2.getObsGraph()->SetLineWidth(2);
// 	lb2.getLine()->SetLineColor(kRed);
// 	lb.getYellowGraph()->Draw("f");
// 	lb.getGreenGraph()->Draw("f");
// 	lb.getMeanGraph()->Draw("l");
// 	lb.getObsGraph()->SetLineWidth(2);
// 	lb.getObsGraph()->Draw("l");
// 	lb2.save();	

// 	lb.getCanvas()->cd();
// 	lb2.getYellowGraph()->Draw("f");
// 	lb2.getGreenGraph()->Draw("f");
// 	lb2.getMeanGraph()->Draw("l");
// 	lb2.getObsGraph()->SetLineWidth(2);
// 	lb2.getObsGraph()->Draw("l");
//	lb.save();	


}
