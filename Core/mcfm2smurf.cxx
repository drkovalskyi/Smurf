{
  gSystem->CompileMacro("mcfm2smurf.C","k");
  TFile* f = TFile::Open("WWqqbr_tota_mstw8nl_80__80__average_scale.root");
  TTree* tree = (TTree*)f->Get("h10");
  mcfm2smurf processor(tree);
  processor.MakeSmurfNtuple("test.root",61);
  // processor.GetEntry(0);
  // processor.Show();

}
