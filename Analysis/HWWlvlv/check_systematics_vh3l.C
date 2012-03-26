void check_systematics_vh3l(TString inputFilevh3l = "vh3l120.input.root"){
//gInterpreter->ExecuteMacro("/home/ceballos/root/macros/MITStyle.C");
gStyle->SetOptStat(0);
TFile *finputFilevh3l = TFile::Open(inputFilevh3l,"READ");
finputFilevh3l->cd();
TCanvas* cvsvh3l[90];
for(int i=0; i<90; i++) cvsvh3l[i] = new TCanvas(Form("cvh3l_%d",i),Form("cvh3l_%d",i),100,100,700,700);
// Begin making nice histograms
atributes(histo_VHtt,"BDT Output",1,"Events",1);
atributes(histo_VHww,"BDT Output",1,"Events",1);
atributes(histo_WZ,"BDT Output",1,"Events",1);
atributes(histo_ZZ,"BDT Output",1,"Events",1);
atributes(histo_Wjets,"BDT Output",1,"Events",1);

atributes(histo_VHtt_CMS_MVAVHttStatBounding_vh3lUp  ,"BDT Output",2,"Events",2);
atributes(histo_VHtt_CMS_MVAJESBoundingUp            ,"BDT Output",2,"Events",2);
atributes(histo_VHtt_CMS_MVALepEffBoundingUp         ,"BDT Output",2,"Events",2);
atributes(histo_VHtt_CMS_MVALepResBoundingUp         ,"BDT Output",2,"Events",2);
atributes(histo_VHtt_CMS_MVAMETResBoundingUp         ,"BDT Output",2,"Events",2);
atributes(histo_WZ_CMS_MVAWZStatBounding_vh3lUp      ,"BDT Output",2,"Events",2);
atributes(histo_WZ_CMS_MVAJESBoundingUp              ,"BDT Output",2,"Events",2);
atributes(histo_WZ_CMS_MVALepEffBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_WZ_CMS_MVALepResBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_WZ_CMS_MVAMETResBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_VHww_CMS_MVAVHwwStatBounding_vh3lUp  ,"BDT Output",2,"Events",2);
atributes(histo_VHww_CMS_MVAJESBoundingUp            ,"BDT Output",2,"Events",2);
atributes(histo_VHww_CMS_MVALepEffBoundingUp         ,"BDT Output",2,"Events",2);
atributes(histo_VHww_CMS_MVALepResBoundingUp         ,"BDT Output",2,"Events",2);
atributes(histo_VHww_CMS_MVAMETResBoundingUp         ,"BDT Output",2,"Events",2);
atributes(histo_ZZ_CMS_MVAZZStatBounding_vh3lUp      ,"BDT Output",2,"Events",2);
atributes(histo_ZZ_CMS_MVAJESBoundingUp              ,"BDT Output",2,"Events",2);
atributes(histo_ZZ_CMS_MVALepEffBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_ZZ_CMS_MVALepResBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_ZZ_CMS_MVAMETResBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWjetsStatBounding_vh3lUp,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWBoundingUp             ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAMCWBoundingUp           ,"BDT Output",2,"Events",2);
atributes(histo_WZ_CMS_MVAWZNLOBoundingUp            ,"BDT Output",2,"Events",2);
atributes(histo_ZZ_CMS_MVAZZNLOBoundingUp            ,"BDT Output",2,"Events",2);

atributes(histo_VHtt_CMS_MVAVHttStatBounding_vh3lDown  ,"BDT Output",4,"Events",4);
atributes(histo_VHtt_CMS_MVAJESBoundingDown            ,"BDT Output",4,"Events",4);
atributes(histo_VHtt_CMS_MVALepEffBoundingDown         ,"BDT Output",4,"Events",4);
atributes(histo_VHtt_CMS_MVALepResBoundingDown         ,"BDT Output",4,"Events",4);
atributes(histo_VHtt_CMS_MVAMETResBoundingDown         ,"BDT Output",4,"Events",4);
atributes(histo_WZ_CMS_MVAWZStatBounding_vh3lDown      ,"BDT Output",4,"Events",4);
atributes(histo_WZ_CMS_MVAJESBoundingDown              ,"BDT Output",4,"Events",4);
atributes(histo_WZ_CMS_MVALepEffBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_WZ_CMS_MVALepResBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_WZ_CMS_MVAMETResBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_VHww_CMS_MVAVHwwStatBounding_vh3lDown  ,"BDT Output",4,"Events",4);
atributes(histo_VHww_CMS_MVAJESBoundingDown            ,"BDT Output",4,"Events",4);
atributes(histo_VHww_CMS_MVALepEffBoundingDown         ,"BDT Output",4,"Events",4);
atributes(histo_VHww_CMS_MVALepResBoundingDown         ,"BDT Output",4,"Events",4);
atributes(histo_VHww_CMS_MVAMETResBoundingDown         ,"BDT Output",4,"Events",4);
atributes(histo_ZZ_CMS_MVAZZStatBounding_vh3lDown      ,"BDT Output",4,"Events",4);
atributes(histo_ZZ_CMS_MVAJESBoundingDown              ,"BDT Output",4,"Events",4);
atributes(histo_ZZ_CMS_MVALepEffBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_ZZ_CMS_MVALepResBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_ZZ_CMS_MVAMETResBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWjetsStatBounding_vh3lDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWBoundingDown             ,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAMCWBoundingDown           ,"BDT Output",4,"Events",4);
atributes(histo_WZ_CMS_MVAWZNLOBoundingDown	       ,"BDT Output",4,"Events",4);
atributes(histo_ZZ_CMS_MVAZZNLOBoundingDown	       ,"BDT Output",4,"Events",4);

// End making nice histograms

cvsvh3l[0]->cd();
histo_VHww_CMS_MVAVHwwStatBounding_vh3lUp  ->Draw("hist");
histo_VHww_CMS_MVAVHwwStatBounding_vh3lDown->Draw("hist,same");
histo_VHww->Draw("e,same");

cvsvh3l[1]->cd();
histo_VHtt_CMS_MVAVHttStatBounding_vh3lUp  ->Draw("hist");
histo_VHtt_CMS_MVAVHttStatBounding_vh3lDown->Draw("hist,same");
histo_VHtt->Draw("e,same");

cvsvh3l[2]->cd();
histo_WZ_CMS_MVAWZStatBounding_vh3lUp  ->Draw("hist");
histo_WZ_CMS_MVAWZStatBounding_vh3lDown->Draw("hist,same");
histo_WZ->Draw("e,same");

cvsvh3l[3]->cd();
histo_ZZ_CMS_MVAZZStatBounding_vh3lUp  ->Draw("hist");
histo_ZZ_CMS_MVAZZStatBounding_vh3lDown->Draw("hist,same");
histo_ZZ->Draw("e,same");

cvsvh3l[4]->cd();
histo_Wjets_CMS_MVAWjetsStatBounding_vh3lUp  ->Draw("hist");
histo_Wjets_CMS_MVAWjetsStatBounding_vh3lDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// Lep eff syst
cvsvh3l[5]->cd();
histo_VHww_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_VHww_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_VHww->Draw("e,same");

cvsvh3l[6]->cd();
histo_VHtt_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_VHtt_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_VHtt->Draw("e,same");

cvsvh3l[7]->cd();
histo_WZ_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_WZ_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_WZ->Draw("e,same");

cvsvh3l[8]->cd();
histo_ZZ_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ZZ_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ZZ->Draw("e,same");

// Lep Res
cvsvh3l[9]->cd();
histo_VHww_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_VHww_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_VHww->Draw("e,same");

cvsvh3l[10]->cd();
histo_VHtt_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_VHtt_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_VHtt->Draw("e,same");

cvsvh3l[11]->cd();
histo_WZ_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_WZ_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_WZ->Draw("e,same");

cvsvh3l[12]->cd();
histo_ZZ_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ZZ_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ZZ->Draw("e,same");

// METRes
cvsvh3l[13]->cd();
histo_VHww_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_VHww_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_VHww->Draw("e,same");

cvsvh3l[14]->cd();
histo_VHtt_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_VHtt_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_VHtt->Draw("e,same");

cvsvh3l[15]->cd();
histo_WZ_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_WZ_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_WZ->Draw("e,same");

cvsvh3l[16]->cd();
histo_ZZ_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ZZ_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ZZ->Draw("e,same");

// JES
cvsvh3l[17]->cd();
histo_VHww_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_VHww_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_VHww->Draw("e,same");

cvsvh3l[18]->cd();
histo_VHtt_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_VHtt_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_VHtt->Draw("e,same");

cvsvh3l[19]->cd();
histo_WZ_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_WZ_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_WZ->Draw("e,same");

cvsvh3l[20]->cd();
histo_ZZ_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ZZ_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ZZ->Draw("e,same");

// Wjets Syst1
cvsvh3l[21]->cd();
histo_Wjets_CMS_MVAWBoundingUp  ->Draw("hist");
histo_Wjets_CMS_MVAWBoundingDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// Wjets Syst2
cvsvh3l[22]->cd();
histo_Wjets_CMS_MVAMCWBoundingDown->Draw("hist");
histo_Wjets->Draw("e,same");
histo_Wjets_CMS_MVAMCWBoundingUp  ->Draw("hist,same");

// WZ NLO
cvsvh3l[23]->cd();
histo_WZ_CMS_MVAWZNLOBoundingUp  ->Draw("hist");
histo_WZ_CMS_MVAWZNLOBoundingDown->Draw("hist,same");
histo_WZ->Draw("e,same");

// ZZ NLO
cvsvh3l[24]->cd();
histo_ZZ_CMS_MVAZZNLOBoundingDown->Draw("hist");
histo_ZZ_CMS_MVAZZNLOBoundingUp  ->Draw("hist,same");
histo_ZZ->Draw("e,same");

for(int i=0; i<=24; i++) cvsvh3l[i]->SaveAs(Form("plots/cvsvh3l_%d.pdf",i));
for(int i=0; i<=24; i++) cvsvh3l[i]->SaveAs(Form("plots/cvsvh3l_%d.eps",i));
}
void atributes(TH1D *histo, Char_t xtitle[]="", Int_t COLOR = 1, Char_t ytitle[]="Fraction", Int_t style = 1){
  histo->ResetAttLine();
  histo->ResetAttFill();
  histo->ResetAttMarker();
  histo->GetYaxis()->SetNdivisions(505);
  histo->GetXaxis()->SetNdivisions(505);
  histo->SetTitle("");
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.8);
  histo->SetLineStyle(style);
  histo->SetLineWidth(4);
  histo->SetLineColor(COLOR);
  histo->SetMarkerStyle(kFullDotLarge);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->SetTitleOffset(1.);
  histo->GetYaxis()->SetTitleOffset(1.3*histo->GetYaxis()->GetTitleOffset());
  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->CenterTitle(kTRUE);
}
