void check_systematics(TString inputFilewwof = "hwwof_0j.input.root",
                       TString inputFilewwsf = "hwwsf_0j.input.root" ){

TFile *fInputFilewwof = TFile::Open(inputFilewwof,"READ");
fInputFilewwof->cd();
TCanvas* cvswwof[50];
for(int i=0; i<50; i++) cvswwof[i] = new TCanvas(Form("cwwof_%d",i),Form("cwwof_%d",i),100,100,700,700);
// WW syst
cvswwof[0]->cd();
histo_qqWW_CMS_MVAWWBounding_hwwof_0jUp  ->Draw("hist");
histo_qqWW_CMS_MVAWWBounding_hwwof_0jDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[1]->cd();
histo_ggWW_CMS_MVAWWBounding_hwwof_0jUp  ->Draw("hist");
histo_ggWW_CMS_MVAWWBounding_hwwof_0jDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

// stat. syst
cvswwof[2]->cd();
histo_ttH_CMS_MVAttHStatBounding_hwwof_0jUp  ->Draw("hist");
histo_ttH_CMS_MVAttHStatBounding_hwwof_0jDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwof[3]->cd();
histo_ZH_CMS_MVAZHStatBounding_hwwof_0jUp  ->Draw("hist");
histo_ZH_CMS_MVAZHStatBounding_hwwof_0jDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwof[4]->cd();
histo_WH_CMS_MVAWHStatBounding_hwwof_0jUp  ->Draw("hist");
histo_WH_CMS_MVAWHStatBounding_hwwof_0jDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwof[5]->cd();
histo_qqH_CMS_MVAqqHStatBounding_hwwof_0jUp  ->Draw("hist");
histo_qqH_CMS_MVAqqHStatBounding_hwwof_0jDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwof[6]->cd();
histo_ggH_CMS_MVAggHStatBounding_hwwof_0jUp  ->Draw("hist");
histo_ggH_CMS_MVAggHStatBounding_hwwof_0jDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwof[7]->cd();
histo_qqWW_CMS_MVAqqWWStatBounding_hwwof_0jUp  ->Draw("hist");
histo_qqWW_CMS_MVAqqWWStatBounding_hwwof_0jDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[8]->cd();
histo_ggWW_CMS_MVAggWWStatBounding_hwwof_0jUp  ->Draw("hist");
histo_ggWW_CMS_MVAggWWStatBounding_hwwof_0jDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwof[9]->cd();
histo_VV_CMS_MVAVVStatBounding_hwwof_0jUp  ->Draw("hist");
histo_VV_CMS_MVAVVStatBounding_hwwof_0jDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwof[10]->cd();
histo_Top_CMS_MVATopStatBounding_hwwof_0jUp  ->Draw("hist");
histo_Top_CMS_MVATopStatBounding_hwwof_0jDown->Draw("hist,same");
histo_Top->Draw("e,same");

cvswwof[11]->cd();
histo_Zjets_CMS_MVAZjetsStatBounding_hwwof_0jUp  ->Draw("hist");
histo_Zjets_CMS_MVAZjetsStatBounding_hwwof_0jDown->Draw("hist,same");
histo_Zjets->Draw("e,same");

cvswwof[12]->cd();
histo_Wjets_CMS_MVAWjetsStatBounding_hwwof_0jUp  ->Draw("hist");
histo_Wjets_CMS_MVAWjetsStatBounding_hwwof_0jDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

cvswwof[13]->cd();
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwof_0jUp  ->Draw("hist");
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwof_0jDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwof[14]->cd();
histo_Ztt_CMS_MVAZttStatBounding_hwwof_0jUp  ->Draw("hist");
histo_Ztt_CMS_MVAZttStatBounding_hwwof_0jDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Lep eff syst
cvswwof[15]->cd();
histo_ttH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwof[16]->cd();
histo_ZH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwof[17]->cd();
histo_WH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_WH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwof[18]->cd();
histo_qqH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwof[19]->cd();
histo_ggH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwof[20]->cd();
histo_qqWW_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[21]->cd();
histo_ggWW_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwof[22]->cd();
histo_VV_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_VV_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwof[23]->cd();
histo_Wgamma_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwof[24]->cd();
histo_Ztt_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Lep Res
cvswwof[25]->cd();
histo_ttH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwof[26]->cd();
histo_ZH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwof[27]->cd();
histo_WH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_WH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwof[28]->cd();
histo_qqH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwof[29]->cd();
histo_ggH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwof[30]->cd();
histo_qqWW_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[31]->cd();
histo_ggWW_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwof[32]->cd();
histo_VV_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_VV_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwof[33]->cd();
histo_Top_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_Top_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_Top->Draw("e,same");

cvswwof[34]->cd();
histo_Wgamma_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwof[35]->cd();
histo_Ztt_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// METRes
cvswwof[36]->cd();
histo_ttH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwof[37]->cd();
histo_ZH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwof[38]->cd();
histo_WH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_WH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwof[39]->cd();
histo_qqH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwof[40]->cd();
histo_ggH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwof[41]->cd();
histo_qqWW_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[42]->cd();
histo_ggWW_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwof[43]->cd();
histo_VV_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_VV_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwof[44]->cd();
histo_Wgamma_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwof[45]->cd();
histo_Ztt_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Wjets syst
cvswwof[46]->cd();
histo_Wjets_CMS_MVAWBounding_hwwof_0jUp  ->Draw("hist");
histo_Wjets_CMS_MVAWBounding_hwwof_0jDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// ggH syst
cvswwof[47]->cd();
histo_ggH_CMS_MVAggHBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVAggHBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

TFile *fInputFilewwsf = TFile::Open(inputFilewwsf,"READ");
fInputFilewwsf->cd();
TCanvas* cvswwsf[50];
for(int i=0; i<50; i++) cvswwsf[i] = new TCanvas(Form("cwwsf_%d",i),Form("cwwsf_%d",i),100,100,700,700);
// WW syst
cvswwsf[0]->cd();
histo_qqWW_CMS_MVAWWBounding_hwwsf_0jUp  ->Draw("hist");
histo_qqWW_CMS_MVAWWBounding_hwwsf_0jDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[1]->cd();
histo_ggWW_CMS_MVAWWBounding_hwwsf_0jUp  ->Draw("hist");
histo_ggWW_CMS_MVAWWBounding_hwwsf_0jDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

// stat. syst
cvswwsf[2]->cd();
histo_ttH_CMS_MVAttHStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_ttH_CMS_MVAttHStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwsf[3]->cd();
histo_ZH_CMS_MVAZHStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_ZH_CMS_MVAZHStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwsf[4]->cd();
histo_WH_CMS_MVAWHStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_WH_CMS_MVAWHStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwsf[5]->cd();
histo_qqH_CMS_MVAqqHStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_qqH_CMS_MVAqqHStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwsf[6]->cd();
histo_ggH_CMS_MVAggHStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_ggH_CMS_MVAggHStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwsf[7]->cd();
histo_qqWW_CMS_MVAqqWWStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_qqWW_CMS_MVAqqWWStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[8]->cd();
histo_ggWW_CMS_MVAggWWStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_ggWW_CMS_MVAggWWStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwsf[9]->cd();
histo_VV_CMS_MVAVVStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_VV_CMS_MVAVVStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwsf[10]->cd();
histo_Top_CMS_MVATopStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_Top_CMS_MVATopStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_Top->Draw("e,same");

cvswwsf[11]->cd();
histo_Zjets_CMS_MVAZjetsStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_Zjets_CMS_MVAZjetsStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_Zjets->Draw("e,same");

cvswwsf[12]->cd();
histo_Wjets_CMS_MVAWjetsStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_Wjets_CMS_MVAWjetsStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

cvswwsf[13]->cd();
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwsf[14]->cd();
histo_Ztt_CMS_MVAZttStatBounding_hwwsf_0jUp  ->Draw("hist");
histo_Ztt_CMS_MVAZttStatBounding_hwwsf_0jDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Lep eff syst
cvswwsf[15]->cd();
histo_ttH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwsf[16]->cd();
histo_ZH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwsf[17]->cd();
histo_WH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_WH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwsf[18]->cd();
histo_qqH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwsf[19]->cd();
histo_ggH_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwsf[20]->cd();
histo_qqWW_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[21]->cd();
histo_ggWW_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwsf[22]->cd();
histo_VV_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_VV_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwsf[23]->cd();
histo_Wgamma_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwsf[24]->cd();
histo_Ztt_CMS_MVALepEffBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVALepEffBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Lep Res
cvswwsf[25]->cd();
histo_ttH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwsf[26]->cd();
histo_ZH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwsf[27]->cd();
histo_WH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_WH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwsf[28]->cd();
histo_qqH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwsf[29]->cd();
histo_ggH_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwsf[30]->cd();
histo_qqWW_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[31]->cd();
histo_ggWW_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwsf[32]->cd();
histo_VV_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_VV_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwsf[33]->cd();
histo_Top_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_Top_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_Top->Draw("e,same");

cvswwsf[34]->cd();
histo_Wgamma_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwsf[35]->cd();
histo_Ztt_CMS_MVALepResBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVALepResBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// METRes
cvswwsf[36]->cd();
histo_ttH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwsf[37]->cd();
histo_ZH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwsf[38]->cd();
histo_WH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_WH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwsf[39]->cd();
histo_qqH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwsf[40]->cd();
histo_ggH_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwsf[41]->cd();
histo_qqWW_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[42]->cd();
histo_ggWW_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwsf[43]->cd();
histo_VV_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_VV_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwsf[44]->cd();
histo_Wgamma_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwsf[45]->cd();
histo_Ztt_CMS_MVAMETResBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVAMETResBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Wjets syst
cvswwsf[46]->cd();
histo_Wjets_CMS_MVAWBounding_hwwsf_0jUp  ->Draw("hist");
histo_Wjets_CMS_MVAWBounding_hwwsf_0jDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

cvswwsf[47]->cd();
histo_ggH_CMS_MVAggHBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVAggHBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

// Zjets syst
cvswwsf[48]->cd();
histo_Zjets_CMS_MVAZBounding_hwwsf_0jUp  ->Draw("hist");
histo_Zjets_CMS_MVAZBounding_hwwsf_0jDown->Draw("hist,same");
histo_Zjets->Draw("e,same");

for(int i=0; i<49; i++) cvswwof[i]->SaveAs(Form("plots/cvswwof_%d.eps",i));
for(int i=0; i<49; i++) cvswwsf[i]->SaveAs(Form("plots/cvswwsf_%d.eps",i));
}
