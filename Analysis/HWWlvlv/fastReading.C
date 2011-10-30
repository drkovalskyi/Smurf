void fastReading
(
 TString sigInputFile    = "data/ntuples_115train_1jets_hww115.root",
 TString bgdInputFile    = "data/ntuples_115train_1jets_backgroundA_skim2.root",
 TString datInputFile    = "data/ntuples_115train_1jets_data_2l_skim2.root",
 TString systInputFile   = "data/ntuples_115train_1jets_hww_syst_skim3.root"
 )
{
  TString sigFile1 = sigInputFile;
  TString bgdFile1 = bgdInputFile;
  TString datFile1 = datInputFile;

  TChain *chsignal = new TChain("tree");
  chsignal->Add(sigFile1);

  TChain *chbackground = new TChain("tree");
  chbackground->Add(bgdFile1);

  TChain *chdata = new TChain("tree");
  chdata->Add(datFile1);

  TTree *signal     = (TTree*) chsignal;
  TTree *background = (TTree*) chbackground;
  TTree *data       = (TTree*) chdata;

  TChain *chsystInputFile = new TChain("tree");
  chsystInputFile->Add(systInputFile);
  TTree *treeSyst = (TTree*) chsystInputFile;
  
  printf("sig: %d bck: %d data: %d syst: %d\n",chsignal->GetEntries(),chbackground->GetEntries(),chdata->GetEntries(),chsystInputFile->GetEntries());
}
