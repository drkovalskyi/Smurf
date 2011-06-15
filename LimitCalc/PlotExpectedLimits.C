#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->Load("lands.so");
  gSystem->CompileMacro("PlotExpectedLimits.C","k");
  PlotExpectedLimits("output/hww-1000pb_0-1_cut*","H #rightarrow WW #rightarrow 2l2#nu + 0/1 jets");
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

#ifndef __CINT__
#include "../../LandS/include/PlotUtilities.h"
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

void AddLimits(std::vector<LimitInfo>& limits, const char* pattern1, const char* pattern2, double mass){
  TChain* chain = new TChain("T");
  chain->Add(Form("%s%s",pattern1,pattern2));
  double result;
  chain->SetBranchAddress("brT", &result);
  Long64_t nentries = chain->GetEntries();
  if (nentries<=0) return;
  std::vector<double> values; 
  double sum(0);
  for (Long64_t i = 0; i < nentries; i++){
    chain->GetEvent(i);
    values.push_back(result);
    sum+=result;
  }
  std::sort(values.begin(),values.end());
  const double prob1S = TMath::Erfc(1/sqrt(2))/2;
  const double prob2S = TMath::Erfc(2/sqrt(2))/2;
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
  for(unsigned int i=0; i<limits.size(); ++i){
    limits_m2s[i]  = limits.at(i).exp_m2sig;
    limits_m1s[i]  = limits.at(i).exp_m1sig;
    // limits_mean[i] = limits.at(i).exp_mean;
    limits_mean[i] = limits.at(i).exp_median;
    limits_p1s[i]  = limits.at(i).exp_p1sig;
    limits_p2s[i]  = limits.at(i).exp_p2sig;
    mass_points[i] = limits.at(i).mass;
  }
  gROOT->SetStyle("Plain");
  gStyle->SetCanvasDefH      (600);
  gStyle->SetCanvasDefW      (800);
  gStyle->SetTitleOffset( 1, "y");
  TPaveText *pt = lands::SetTPaveText(0.5, 0.95, 0.8, 0.95); //SetTPaveText(0.5, 0.95, 0.8, 0.95)
  pt->AddText(title);
  lands::PlotWithBelts* lb = new lands::PlotWithBelts(limits_m1s, limits_p1s, limits_m2s, limits_p2s,
						      limits_mean, limits.size(), mass_points, 
						      "limits", pt,
						      100, 300, 0, 5, false, 
						      ";Higgs mass, m_{H} [GeV/c^{2}]; 95% CL Limit on #sigma/#sigma_{SM} ");
  lb->plot();
  lb->drawLegend("95% CL exclusion: median","95% CL exclusion: 68% band", "95% CL exclusion: 95% band");
  // lb->drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band");
  lands::MoveLegend(lb->getLegend(),0.45,0.7);
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

#endif

void PlotExpectedLimits(const char* pattern, const char* title="H #rightarrow WW #rightarrow 2l2#nu + 0/1 jets")
{
  std::vector<LimitInfo> limits;
  AddLimits(limits, pattern, "115-*.root", 115);
  AddLimits(limits, pattern, "120-*.root", 120);
  AddLimits(limits, pattern, "130-*.root", 130);
  AddLimits(limits, pattern, "140-*.root", 140);
  AddLimits(limits, pattern, "150-*.root", 150);
  AddLimits(limits, pattern, "160-*.root", 160);
  AddLimits(limits, pattern, "170-*.root", 170);
  AddLimits(limits, pattern, "180-*.root", 180);
  AddLimits(limits, pattern, "190-*.root", 190);
  AddLimits(limits, pattern, "200-*.root", 200);
  AddLimits(limits, pattern, "250-*.root", 250);
  AddLimits(limits, pattern, "300-*.root", 300);
  
  PlotExpectedLimits(limits,title);
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
	pt = lands::SetTPaveText(0.2, 0.7, 0.5, 0.9);
	pt->AddText("HZZ #rightarrow 2l2#nu");
	lands::PlotWithBelts lb(
			limits_m1s, limits_p1s, limits_m2s, limits_p2s,
			limits_mean, npoints,  
			xpoints, ssave+"_1", pt,
			projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
	lb.plot();
	lb.drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band");
	lands::MoveLegend(lb.getLegend(),0.45,0.7);
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
