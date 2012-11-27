
#include "../../LandS/test/PlotUtilities_standalone.cc"
int getMedianBin(TH1F& *h);
void hypoSeparation(TString hypo0_name,  TString hypo0_LLH_forToys0, TString hypo1_LLH_forToys0,  
		    TString hypo1_name,  TString hypo0_LLH_forToys1, TString hypo1_LLH_forToys1, int nBins, double xmin, double xmax){
	bool debug = true;

	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetNdivisions(505, "XYZ");
	gStyle->SetPaintTextFormat(".1f");               	
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);	
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.01, "XYZ"); 
	gStyle->SetLabelSize(0.05, "XYZ"); 

	gStyle->SetOptTitle(0);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleFontSize(0.06);      

	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ"); 
	gStyle->SetTitleSize(0.05, "XYZ");
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYOffset(1.4);

	// Likelihoods 
	vector<double> vh0t0, vh0t1, vh1t0, vh1t1;  // toys 
	vector<double> dh0, dh1;  // data

	TTree *tvh0t0 = (TTree*) GetTObject(hypo0_LLH_forToys0.Data(), "T");
	TTree *tvh0t1 = (TTree*) GetTObject(hypo0_LLH_forToys1.Data(), "T");
	TTree *tvh1t0 = (TTree*) GetTObject(hypo1_LLH_forToys0.Data(), "T");
	TTree *tvh1t1 = (TTree*) GetTObject(hypo1_LLH_forToys1.Data(), "T");

	TTree *tdh0= (TTree*) GetTObject(hypo0_LLH_forToys0.Data(), "T_data");
	TTree *tdh1= (TTree*) GetTObject(hypo1_LLH_forToys0.Data(), "T_data");

	vh0t0 = GetVectorFrom(tvh0t0, "likelihood");
	vh0t1 = GetVectorFrom(tvh0t1, "likelihood");
	vh1t0 = GetVectorFrom(tvh1t0, "likelihood");
	vh1t1 = GetVectorFrom(tvh1t1, "likelihood");
	if(vh0t0.size() != vh1t0.size()) { cout<<"Entries not match !"<<endl; exit(1);}  
	if(vh0t1.size() != vh1t1.size()) { cout<<"Entries not match !"<<endl; exit(1);}  

	dh0= GetVectorFrom(tdh0, "limit");
	dh1= GetVectorFrom(tdh1, "limit");

	vector<double> vts0, vts1; // test statistics --> estimator 
	double dts;
	// for(int i=0; i<vh0t0.size(); i++) vts0.push_back(2*(vh0t0[i] - vh1t0[i]));
	// for(int i=0; i<vh0t1.size(); i++) vts1.push_back(2*(vh0t1[i] - vh1t1[i]));
	// dts= 2*(dh0[0] - dh1[0]);
	
	for(int i=0; i<vh0t0.size(); i++) vts0.push_back(2*(vh1t0[i] - vh0t0[i]));
	for(int i=0; i<vh0t1.size(); i++) vts1.push_back(2*(vh1t1[i] - vh0t1[i]));
	dts= 2*(dh1[0] - dh0[0]);
	

	// plotting   
	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
	TCanvas *cCanvas = new TCanvas("cM2logQ","cM2logQ", 500, 500);

	TH1F *hPdfM2logQ= new TH1F("hPdfM2logQ","",100, 0, 0);
	for(int i=0; i<(int)vts0.size(); i++){
		hPdfM2logQ->Fill(vts0[i]);
	}
	for(int i=0; i<(int)vts1.size(); i++){
		hPdfM2logQ->Fill(vts1[i]);
	}
	TH1F* hTS0 = new TH1F("hPdfM2logQ_0","",nBins, xmin, xmax);
	TH1F* hTS1 = new TH1F("hPdfM2logQ_1","",nBins, xmin, xmax);
	delete hPdfM2logQ;
	hTS0->SetStats(0);
	hTS1->SetStats(0);
	for(int i=0; i<(int)vts0.size(); i++){
		hTS0->Fill(vts0[i]);
	}
	for(int i=0; i<(int)vts1.size(); i++){
		hTS1->Fill(vts1[i]);
	}


	double *q_0 = new double[vts0.size()]; 
	for(int i=0; i<vts0.size(); i++) q_0[i]=vts0[i];
	double *q_1 = new double[vts1.size()]; 
	for(int i=0; i<vts1.size(); i++) q_1[i]=vts1[i];
	sort(q_0, q_0+vts0.size());
	sort(q_1, q_1+vts1.size());

	// calculation of observed and expected significance    // taken from Yanyan Gao (FNAL)

	TH1F *S_H0 = hTS0; 
	TH1F *S_H1 = hTS1; 

	int nBins = S_H0->GetNbinsX();
	double norm0 = S_H0->Integral( 1, nBins );
	double norm1 = S_H1->Integral( 1, nBins );
	double int0(0.), int1(0.);
	double diff = 10.;
	double coverage = 0.;

	int nbin_eq = 0; 
	for (int i = 1; i <= nBins; i++){

		int0 = S_H0->Integral(1,i)/norm0;
		int1 = S_H1->Integral(i,nBins)/norm1;

		if (fabs(int0-int1) < diff){
			diff = fabs(int0-int1);
			coverage = 0.5*(int0+int1);
			nbin_eq = i;
		}
	}  


	std::cout << "coverage : " << coverage << ", for bin " << nbin_eq << "\n";
	double sepH = 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
	std::cout << "histogram separation is: " <<  sepH << ", with sigma coverage: " << coverage << std::endl;


	// 
	// Calculate the fraction of events of the hypothesis 2 in the hyp1's plot
	//
	int bin_median_H1 = getMedianBin(S_H1);
	std::cout << "bin_median_H1 = " << bin_median_H1 << "\n";

	double frac_H0_beyondH1Median = S_H0->Integral(bin_median_H1, nBins) / norm0;
	double sepH0vsH1 = ROOT::Math::normal_quantile_c(frac_H0_beyondH1Median, 1.0);
	std::cout << "frac of H0 histogram beyond the H1 median " << frac_H0_beyondH1Median << ", correspond to " << sepH0vsH1 << " sigma\n";


	int bin_median_H0 = getMedianBin(S_H0);
	std::cout << "bin_median_H0 = " << bin_median_H0 << "\n";

	double frac_H1_beyondH0Median = S_H1->Integral(1, bin_median_H0) / norm1;
	double sepH1vsH0 = ROOT::Math::normal_quantile_c(frac_H1_beyondH0Median, 1.0);
	std::cout << "frac of H1 histogram beyond the H0 median " << frac_H1_beyondH0Median << ", correspond to " << sepH1vsH0 << " sigma\n";


	// SM plotting settins
	// to be consistent with HZZ
	// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HZZ4L_Combination/CombinationPy/SignalSeparation/extractSignificanceStats.C?revision=1.9&view=markup
	hTS0->SetLineColor(kMagenta-3);
	hTS0->SetFillColor(kMagenta-3);
	hTS0->SetLineWidth(2);
	hTS0->SetFillStyle(3605);
	hTS0->SetXTitle("-2 #times ln(L_{2_{min}^{+}}/L_{0^{+}})");
	hTS0->SetYTitle("Generated experiments");
	// hTS0->GetXaxis()->CenterTitle();
	// hTS0->GetYaxis()->CenterTitle();
	
	hTS1->SetLineColor(kBlue+1);
	hTS1->SetFillColor(kBlue+1);
	hTS1->SetLineWidth(2);
	hTS1->SetFillStyle(3695);
	hTS1->SetXTitle("-2 #times ln(L_{2_{min}^{+}}/L_{0^{+}})");
	hTS1->SetYTitle("Generated experiments");
	// hTS1->GetXaxis()->CenterTitle();
	// hTS1->GetYaxis()->CenterTitle();

	
	hTS0->Rebin(10);
	hTS1->Rebin(10);

        hTS0->SetMaximum(max(hTS0->GetMaximum(), hTS1->GetMaximum())*1.15);

	hTS0->Draw();
	hTS1->Draw("same");

	double ymin, ymax;
	ymin=1e-2;
	ymax=hTS0->GetMaximum();
	if(ymax<hTS1->GetMaximum()) ymax=hTS1->GetMaximum();
	TLine *lineOne = new TLine(dts, ymin, dts, ymax);
	lineOne->SetLineWidth(2);
	lineOne->SetLineStyle(1);
	lineOne->SetLineColor(kBlack);
	// lineOne->Draw("same");

	TLegend *legend = new TLegend(0.25,0.65,0.45,0.85);
	legend->SetFillColor(0);
	legend->SetTextSize(0.05);
	legend->SetBorderSize(0);
	legend->AddEntry(hTS0, "SM, 0^{+}", "f");
	legend->AddEntry(hTS1, "2_{min}^{+}", "f");
	// legend->AddEntry(lineOne, "data", "l");
	legend->Draw("same");
	pt->Draw();
	Save(cCanvas, "hypo_separation");

}
int getMedianBin(TH1F& *h) {
  
  if ( h == 0 ) {
    return 0;
  }
  int nBins = h->GetNbinsX();
  
  double norm = h->Integral(1, nBins);
  int bin_median = 1;
  double diff = 1;

  for ( int i = 1; i < nBins; i++) {
    double frac = h->Integral(1, i) / norm;
    double diff_bin = fabs(frac - 0.5 );
    if ( diff_bin < diff ) {
      diff = diff_bin;
      bin_median = i;
    }
    // std::cout << "Bin " << i << ", frac " << frac << ", diff_bin " << diff_bin << "\n";
  }
  
  return bin_median;
}


