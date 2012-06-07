#include <algorithm>
const char* n1 = "/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dytt.root";
const char* n2 = "dytt_from_z.root";
void compareTree(const char* var, Int_t nbins, double xlow, double xup, const char* cut="")
{
  f1 = TFile::Open(n1);
  t1 = (TTree*)f1->Get("tree");
  std::string command1(Form("%s>>h1(%d,%f,%f)",var,nbins,xlow,xup));
  std::string command2(Form("%s>>h2(%d,%f,%f)",var,nbins,xlow,xup));
  t1->Draw(command1.c_str(),cut);
  h1 = (TH1*)gDirectory->Get("h1");
  h1->SetDirectory(0);
  h1->Scale(1/h1->Integral());
  f2 = TFile::Open(n2);
  t2 = (TTree*)f2->Get("tree");
  t2->Draw(command2.c_str(),cut);
  h2 = (TH1*)gDirectory->Get("h2");
  h2->SetDirectory(0);
  h2->Scale(1/h2->Integral());
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(2);
  double hmax = std::max(h1->GetMaximum(),h2->GetMaximum());
  hmax *= 1.2;
  h1->SetMaximum(hmax);
  h2->SetMaximum(hmax);
  h1->Draw();
  h2->Draw("same");
}

void compareTree(){
  c1 = new TCanvas("c1","c1",3*300,4*300);
  c1->Divide(3,4);
  c1->cd(1);
  const char* cut = "lep1.pt()>20&&(cuts&518)==518&&njets==0";
  compareTree("lep1.pt()",100,0,100,cut);
  c1->cd(2);
  compareTree("lep2.pt()",100,0,100,cut);
  c1->cd(3);
  compareTree("dilep.mass()",100,0,100,cut);
  c1->cd(4);
  compareTree("met",100,0,100,cut);
  c1->cd(5);
  compareTree("mt",100,0,100,cut);
  c1->cd(6);
  compareTree("dilep.pt()",100,0,100,cut);
  c1->cd(7);
  compareTree("pmet",100,0,50,cut);
  c1->cd(8);
  compareTree("pTrackMet",100,0,50,cut);
  c1->cd(9);
  compareTree("min(pmet,pTrackMet)",100,0,50,cut);
  c1->cd(10);
  compareTree("dPhi",100,0,3.2,cut);
  c1->cd(11);
  compareTree("dPhiDiLepJet1",100,0,3.2,cut);
  c1->cd(12);
  compareTree("dPhiDiLepMET",100,0,3.2,cut);

  c2 = new TCanvas("c2","c2",3*300,4*300);
  c2->Divide(3,4);
  c2->cd(1);
  const char* cut = "lep1.pt()>20&&(cuts&518)==518&&njets==1";
  compareTree("lep1.pt()",100,0,100,cut);
  c2->cd(2);
  compareTree("lep2.pt()",100,0,100,cut);
  c2->cd(3);
  compareTree("dilep.mass()",100,0,100,cut);
  c2->cd(4);
  compareTree("met",100,0,100,cut);
  c2->cd(5);
  compareTree("mt",100,0,100,cut);
  c2->cd(6);
  compareTree("dilep.pt()",100,0,100,cut);
  c2->cd(7);
  compareTree("pmet",100,0,50,cut);
  c2->cd(8);
  compareTree("pTrackMet",100,0,50,cut);
  c2->cd(9);
  compareTree("min(pmet,pTrackMet)",100,0,50,cut);
  c2->cd(10);
  compareTree("dPhi",100,0,3.2,cut);
  c2->cd(11);
  compareTree("dPhiDiLepJet1",100,0,3.2,cut);
  c2->cd(12);
  compareTree("dPhiDiLepMET",100,0,3.2,cut);


}
