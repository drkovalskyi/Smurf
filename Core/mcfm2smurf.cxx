{
  gSystem->CompileMacro("mcfm2smurf.C","k");
  // TFile* f = TFile::Open("WWqqbr_tota_mstw8nl_80__80__average_scale.root");
  TFile* f = TFile::Open("WZbbar_tota_mstw8lo_85__85__id71_200k_M1.root");
  TTree* tree = (TTree*)f->Get("h10");
  mcfm2smurf processor(tree);
  processor.MakeSmurfNtuple("WZbbar_tota_mstw8lo_85__85__id71_200k_M1_smurf_isoalted.root",71);
  // processor.GetEntry(0);
  // processor.Show();

}
