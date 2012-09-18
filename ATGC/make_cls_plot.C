bool getParameters(const char* name, double& x, double& y){
  TPMERegexp re("x(\-?[0-9]*\.?[0-9]*e?\-?[0-9]*)_y(\-?[0-9]*\.?[0-9]*e?\-?[0-9]*)");
  re.Match(name);
  if (re.NMatches()<3){
    x=0;
    y=0;
    return false;
  }
  x = atof(re[1].Data());
  y = atof(re[2].Data());
  return true;
}
TGraph2D* make_cls_plot(const char* pattern = "cards/higgsCombine*root"){
  TChain* chain = new TChain("limit");
  Int_t nfiles = chain->Add(pattern);
  printf("Number of files found: %d\n",nfiles);
  Double_t limit(0);
  chain->SetBranchAddress("limit", &limit);
  // chain->GetCurrentFile()->GetName();
  Int_t nevent = chain->GetEntries();
  TGraph2D* graphs[6];
  for (int i=0;i<6;++i) graphs[i] = new TGraph2D();
  double x(0);
  double y(0);
  for (Int_t i=0;i<nevent;i++) {
    chain->GetEvent(i);
    Int_t index = i%6;
    if (!getParameters(chain->GetCurrentFile()->GetName(),x,y)) continue;
    if (fabs(x+0.05)<1e-3&&fabs(y)<1e-3)
      cout << i << "\t" << index << "\t" << x << "\t" << y << "\t" << limit << endl;
    graphs[index]->SetPoint(int(i/6),x,y,limit);
  }
  new TCanvas("c1","c1",500,500);
  graphs[5]->Draw("colz");
  return graphs[5];
}
