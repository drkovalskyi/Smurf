{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  h = MakePlot(MakeDataset("smurf/qqww.root","blah"),"Blah");
  new TCanvas("c1","c1",500,500);
  h->Draw();
}
