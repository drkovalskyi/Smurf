void check_bdt(int mass){
  TFile *_file0 =
  TFile::Open(Form("/data/smurf/data/LP2011/mitf_mva_1.5fb/ntuples_%dtrain_1jets_data_2l_METVeto_ZVeto.root",mass));
  TTree* tree = (TTree*)_file0->Get("tree");
  printf ("Mass: %d\n", mass);
  tree->Scan(Form("bdtg_hww%d_1jet_ww",mass),"event==582505290||event==325565599");
  tree->Scan("lep1->pt():lep2->pt():mt:dilep->mass():dPhi*180/TMath::Pi():dR","event==582505290||event==325565599");
}
