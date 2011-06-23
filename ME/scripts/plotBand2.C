int plotBand2(int nsel = 0, TString resdir = "output/limits/", TString input = "bdt", double posX = 0.3, double posY = 0.7){
  gROOT->ProcessLine(".L PlotUtilities_standalone.cc++");
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");    
  gStyle->SetCanvasDefH      (600);
  gStyle->SetCanvasDefW      (800);
  gStyle->SetTitleOffset( 1, "y");
  
  const int nPoints = 12;
  double xpoints[nPoints] = {115, 120,130,140,150,160,170,180,190,200,250,300};
  double projectingXmin = 100, projectingXmax = 310;
  double projectingRLimitYmin = 0.1, projectingRLimitYmax = 5;
  
  double obs0j[nPoints],exp0j_2sm[nPoints],exp0j_1sm[nPoints],exp0j[nPoints],exp0j_1sp[nPoints],exp0j_2sp[nPoints], exp0j_median[nPoints];
  
  char sb[50];
  int i = 0;
  sprintf(sb,"%s/results_0j_%s.txt", resdir.Data(), input.Data());
  cout << sb <<endl;
  ifstream infile0(sb);
  while (infile0>>obs0j[i]
	 >>exp0j_2sm[i]
	 >>exp0j_1sm[i]
	 >>exp0j[i]
	 >>exp0j_1sp[i]
	 >>exp0j_2sp[i]
	 >>exp0j_median[i]	 
	 ){ 
    i++;
    if (i >= nPoints) break;
    
  }
  
  bool projectingRLimitLogY = true;
  string projectingRLimitXYtitles = ";Higgs mass, m_{H} [GeV/c^{2}]; 95% CL Limit on #sigma/#sigma_{SM} ";
  string ssave = "limits_";
  if     (nsel == 0) ssave+="0j_"+input;
  
  TPaveText *pt;
  pt = SetTPaveText(0.5, 0.95, 0.8, 0.95);
  double limits_m1s[nPoints], limits_m2s[nPoints], limits_p1s[nPoints], limits_p2s[nPoints], limits_mean[nPoints], limits_obs[nPoints];
  for(int n=0; n<nPoints; n++){
    if     (nsel == 0) {
      cout << obs0j[n]<<endl;
      limits_obs[n] =obs0j[n];
      limits_m2s[n] =exp0j_2sm[n];
      limits_m1s[n] =exp0j_1sm[n];
      //       limits_mean[n]=exp0j[n];
      limits_mean[n]=exp0j_median[n];
      limits_p1s[n] =exp0j_1sp[n];
      limits_p2s[n] =exp0j_2sp[n];
      pt->AddText("H #rightarrow WW #rightarrow 2l2#nu + 0 jets");
    }
  }
  
  cout << "Print latex styel tables..........\n";
  for(int i = 0; i < nPoints; i++) {
    printf("%i & %4.2f & %4.2f &  %4.2f &  %4.2f &  %4.2f \\\\\n ", xpoints[i], limits_m2s[i], limits_m1s[i], limits_mean[i],  limits_p1s[i], limits_p2s[i]);
    //printf("%i & %4.2f & %4.2f &  %4.2f &  %4.2f &  %4.2f & %4.2f \\\\\n ", xpoints[i], limits_obs[i], limits_m2s[i], limits_m1s[i], limits_mean[i],  limits_p1s[i], limits_p2s[i]);
  }
  
  cout << "Print twiki styel tables..........\n";
  // cout << "|*mH*|*Observed*|*-2sigma*|*-1sigma*|*mean*|*median*|*+1sigma*|*+2sigma*|\n";
  cout << "|*mH*|*-2sigma*|*-1sigma*|*mean*|*median*|*+1sigma*|*+2sigma*|\n";
  for(int i = 0; i < nPoints; i++) {
    printf("|%i|%4.2f|%4.2f|%4.2f|%4.2f|%4.2f|%4.2f|\n ", xpoints[i], exp0j_2sm[i], exp0j_1sm[i], exp0j[i], exp0j_median[i], exp0j_1sp[i], exp0j_2sp[i]);
  }
  
  PlotWithBelts lb(
		   limits_m1s, limits_p1s, limits_m2s, limits_p2s,
		   //limits_mean, limits_obs, nPoints,  
		   limits_mean, nPoints,  
		   xpoints, ssave+"_1", pt,
		   projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
  lb.plot(posX, posY);
  //lb.drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: Observed");
  lb.drawLegend("95% CL exclusion: median","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: Observed");
  
}
