{
  TFile *_file0 = TFile::Open("/afs/cern.ch/work/d/dmytro/atgc/mcfm/WWqqbr_tota_mstw8nl_80__80__delg1_z_________delk_z_________delk_g_________lambda_z_p0.5000_lambda_g_p0.5000_smurf.root");
  tree->Draw("min(lep1.pt(),199)>>hlz(9,20,200)","scale1fb");
  hlz->SetDirectory(0);
  h1 = hlz;
  TFile *_file0 = TFile::Open("/afs/cern.ch/work/d/dmytro/atgc/mcfm/WWqqbr_tota_mstw8nl_80__80__delg1_z_p0.5000_delk_z_p0.5000_delk_g_________lambda_z_________lambda_g_________smurf.root");
  tree->Draw("min(lep1.pt(),199)>>hgz(9,20,200)","scale1fb");
  h2 = hgz;
  h2->SetDirectory(0);
  TFile *_file0 = TFile::Open("/afs/cern.ch/work/d/dmytro/atgc/mcfm/WWqqbr_tota_mstw8nl_80__80__delg1_z_p0.7500_delk_z_p0.7500_delk_g_________lambda_z_________lambda_g_________smurf.root");
  tree->Draw("min(lep1.pt(),199)>>hgz2(9,20,200)","scale1fb");
  h3 = hgz2;
  h3->SetDirectory(0);

  cc = new TCanvas("cc","",600,600);
  h1->GetXaxis()->SetTitle("P_{T}, [GeV]");
  h1->SetLineWidth(2);
  h1->SetTitle("L_{Z} vs #delta g^{Z}_{1}");
  h1->Draw();
  h2->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->Draw("same");
  h3->SetLineWidth(2);
  h3->SetLineColor(kRed);
  h3->Draw("same");
  cc->Print("cc.root");
}
