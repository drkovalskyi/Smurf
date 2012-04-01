{
  gSystem->CompileMacro("mcfm2smurf.C","k");
  // TFile* f = TFile::Open("WWqqbr_tota_mstw8nl_80__80__average_scale.root");
  TFile* f = TFile::Open("/afs/cern.ch/work/d/dmytro/atgc/mcfm/WWqqbr_tota_mstw8nl_80__80__kg100_lg50_kZ100_lZ50_gZ100.root");
  TTree* tree = (TTree*)f->Get("h10");
  mcfm2smurf processor(tree);
  processor.MakeSmurfNtuple("WWqqbr_tota_mstw8nl_80__80__kg100_lg50_kZ100_lZ50_gZ100_smurf.root",61);
  // processor.GetEntry(0);
  // processor.Show();

}
