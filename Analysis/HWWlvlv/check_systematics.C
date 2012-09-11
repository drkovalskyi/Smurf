void check_systematics(TString inputFilewwof = "hwwof_2j.input_8TeV.root",
                       TString inputFilewwsf = "hwwsf_2j.input_8TeV.root" ){
gInterpreter->ExecuteMacro("/home/ceballos/root/macros/MITStyle.C");
gStyle->SetOptStat(0);
TFile *fInputFilewwof = TFile::Open(inputFilewwof,"READ");
fInputFilewwof->cd();
TCanvas* cvswwof[90];
for(int i=0; i<90; i++) cvswwof[i] = new TCanvas(Form("cwwof_%d",i),Form("cwwof_%d",i),100,100,700,700);
// Begin making nice histograms
atributes(histo_ggH,"BDT Output",1,"Events",1);
atributes(histo_ggWW,"BDT Output",1,"Events",1);
atributes(histo_qqH,"BDT Output",1,"Events",1);
atributes(histo_qqWW,"BDT Output",1,"Events",1);
atributes(histo_Top,"BDT Output",1,"Events",1);
atributes(histo_ttH,"BDT Output",1,"Events",1);
atributes(histo_VV,"BDT Output",1,"Events",1);
atributes(histo_Wgamma,"BDT Output",1,"Events",1);
atributes(histo_WH,"BDT Output",1,"Events",1);
atributes(histo_Wjets,"BDT Output",1,"Events",1);
atributes(histo_ZH,"BDT Output",1,"Events",1);
atributes(histo_Zjets,"BDT Output",1,"Events",1);
atributes(histo_Ztt,"BDT Output",1,"Events",1);

//atributes(histo_ggH_CMS_MVAggHBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVAggHStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVAggWWStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVAqqHStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAqqWWStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAWWBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAWWNLOBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVATopBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVATopStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVAttHStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVAVVStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVAWgammaStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVAWHStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWjetsStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWMCBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVAZHStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Zjets_CMS_MVAZjetsStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVAZttStatBounding_hwwof_2jUp  ,"BDT Output",2,"Events",2);

//atributes(histo_ggH_CMS_MVAggHBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVAggHStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVAggWWStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVAqqHStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAqqWWStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAWWBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAWWNLOBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVATopBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVATopStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVAttHStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVAVVStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVAWgammaStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVAWHStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWjetsStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWMCBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVAZHStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Zjets_CMS_MVAZjetsStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVAZttStatBounding_hwwof_2jDown,"BDT Output",4,"Events",4);
// End making nice histograms

// WW syst1
cvswwof[0]->cd();
histo_qqWW_CMS_MVAWWBounding_hwwUp  ->Draw("hist");
histo_qqWW_CMS_MVAWWBounding_hwwDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

// Top syst
cvswwof[1]->cd();
histo_Top_CMS_MVATopBounding_hwwUp  ->Draw("hist");
histo_Top_CMS_MVATopBounding_hwwDown->Draw("hist,same");
histo_Top->Draw("e,same");

// stat. syst
cvswwof[2]->cd();
histo_ttH_CMS_MVAttHStatBounding_hwwof_2jUp  ->Draw("hist");
histo_ttH_CMS_MVAttHStatBounding_hwwof_2jDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwof[3]->cd();
histo_ZH_CMS_MVAZHStatBounding_hwwof_2jUp  ->Draw("hist");
histo_ZH_CMS_MVAZHStatBounding_hwwof_2jDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwof[4]->cd();
histo_WH_CMS_MVAWHStatBounding_hwwof_2jUp  ->Draw("hist");
histo_WH_CMS_MVAWHStatBounding_hwwof_2jDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwof[5]->cd();
histo_qqH_CMS_MVAqqHStatBounding_hwwof_2jUp  ->Draw("hist");
histo_qqH_CMS_MVAqqHStatBounding_hwwof_2jDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwof[6]->cd();
histo_ggH_CMS_MVAggHStatBounding_hwwof_2jUp  ->Draw("hist");
histo_ggH_CMS_MVAggHStatBounding_hwwof_2jDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwof[7]->cd();
histo_qqWW_CMS_MVAqqWWStatBounding_hwwof_2jUp  ->Draw("hist");
histo_qqWW_CMS_MVAqqWWStatBounding_hwwof_2jDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[8]->cd();
histo_ggWW_CMS_MVAggWWStatBounding_hwwof_2jUp  ->Draw("hist");
histo_ggWW_CMS_MVAggWWStatBounding_hwwof_2jDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwof[9]->cd();
histo_VV_CMS_MVAVVStatBounding_hwwof_2jUp  ->Draw("hist");
histo_VV_CMS_MVAVVStatBounding_hwwof_2jDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwof[10]->cd();
histo_Top_CMS_MVATopStatBounding_hwwof_2jUp  ->Draw("hist");
histo_Top_CMS_MVATopStatBounding_hwwof_2jDown->Draw("hist,same");
histo_Top->Draw("e,same");

cvswwof[11]->cd();
histo_Zjets_CMS_MVAZjetsStatBounding_hwwof_2jUp  ->Draw("hist");
histo_Zjets_CMS_MVAZjetsStatBounding_hwwof_2jDown->Draw("hist,same");
histo_Zjets->Draw("e,same");

cvswwof[12]->cd();
histo_Wjets_CMS_MVAWjetsStatBounding_hwwof_2jUp  ->Draw("hist");
histo_Wjets_CMS_MVAWjetsStatBounding_hwwof_2jDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

cvswwof[13]->cd();
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwof_2jUp  ->Draw("hist");
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwof_2jDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwof[14]->cd();
histo_Ztt_CMS_MVAZttStatBounding_hwwof_2jUp  ->Draw("hist");
histo_Ztt_CMS_MVAZttStatBounding_hwwof_2jDown->Draw("hist,same");
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

// JES
cvswwof[46]->cd();
histo_ttH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwof[47]->cd();
histo_ZH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwof[48]->cd();
histo_WH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_WH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwof[49]->cd();
histo_qqH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwof[50]->cd();
histo_ggH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwof[51]->cd();
histo_qqWW_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[52]->cd();
histo_ggWW_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwof[53]->cd();
histo_VV_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_VV_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwof[54]->cd();
histo_Wgamma_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwof[55]->cd();
histo_Ztt_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Wjets1 syst
cvswwof[56]->cd();
histo_Wjets_CMS_MVAWBounding_hwwUp  ->Draw("hist");
histo_Wjets_CMS_MVAWBounding_hwwDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// Wjets2 syst
cvswwof[57]->cd();
histo_Wjets_CMS_MVAWMCBounding_hwwUp  ->Draw("hist");
histo_Wjets_CMS_MVAWMCBounding_hwwDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// ggH syst
//cvswwof[58]->cd();
//histo_ggH_CMS_MVAggHBoundingUp  ->Draw("hist");
//histo_ggH_CMS_MVAggHBoundingDown->Draw("hist,same");
//histo_ggH->Draw("e,same");

// WW syst2
cvswwof[59]->cd();
histo_qqWW_CMS_MVAWWNLOBounding_hwwUp  ->Draw("hist");
histo_qqWW_CMS_MVAWWNLOBounding_hwwDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwof[60]->cd();
histo_Top_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_Top_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_Top->Draw("e,same");

TFile *fInputFilewwsf = TFile::Open(inputFilewwsf,"READ");
fInputFilewwsf->cd();
TCanvas* cvswwsf[90];
for(int i=0; i<90; i++) cvswwsf[i] = new TCanvas(Form("cwwsf_%d",i),Form("cwwsf_%d",i),100,100,700,700);
// Begin making nice histograms
atributes(histo_ggH,"BDT Output",1,"Events",1);
atributes(histo_ggWW,"BDT Output",1,"Events",1);
atributes(histo_qqH,"BDT Output",1,"Events",1);
atributes(histo_qqWW,"BDT Output",1,"Events",1);
atributes(histo_Top,"BDT Output",1,"Events",1);
atributes(histo_ttH,"BDT Output",1,"Events",1);
atributes(histo_VV,"BDT Output",1,"Events",1);
atributes(histo_Wgamma,"BDT Output",1,"Events",1);
atributes(histo_WH,"BDT Output",1,"Events",1);
atributes(histo_Wjets,"BDT Output",1,"Events",1);
atributes(histo_ZH,"BDT Output",1,"Events",1);
atributes(histo_Zjets,"BDT Output",1,"Events",1);
atributes(histo_Ztt,"BDT Output",1,"Events",1);

//atributes(histo_ggH_CMS_MVAggHBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVAggHStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVAggWWStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ggWW_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqH_CMS_MVAqqHStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAqqWWStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAWWBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_qqWW_CMS_MVAWWNLOBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVATopBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_Top_CMS_MVATopStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ttH_CMS_MVAttHStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_VV_CMS_MVAVVStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wgamma_CMS_MVAWgammaStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_WH_CMS_MVAWHStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWjetsStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Wjets_CMS_MVAWMCBounding_hwwUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_ZH_CMS_MVAZHStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Zjets_CMS_MVAZBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Zjets_CMS_MVAZjetsStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVAJESBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVALepEffBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVALepResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVAMETResBoundingUp  ,"BDT Output",2,"Events",2);
atributes(histo_Ztt_CMS_MVAZttStatBounding_hwwsf_2jUp  ,"BDT Output",2,"Events",2);

//atributes(histo_ggH_CMS_MVAggHBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVAggHStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVAggWWStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ggWW_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqH_CMS_MVAqqHStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAqqWWStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAWWBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_qqWW_CMS_MVAWWNLOBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVATopBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_Top_CMS_MVATopStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ttH_CMS_MVAttHStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_VV_CMS_MVAVVStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Wgamma_CMS_MVAWgammaStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_WH_CMS_MVAWHStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWjetsStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Wjets_CMS_MVAWMCBounding_hwwDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_ZH_CMS_MVAZHStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Zjets_CMS_MVAZBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Zjets_CMS_MVAZjetsStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVAJESBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVALepEffBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVALepResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVAMETResBoundingDown,"BDT Output",4,"Events",4);
atributes(histo_Ztt_CMS_MVAZttStatBounding_hwwsf_2jDown,"BDT Output",4,"Events",4);
// End making nice histograms

// WW syst
cvswwsf[0]->cd();
histo_qqWW_CMS_MVAWWBounding_hwwUp  ->Draw("hist");
histo_qqWW_CMS_MVAWWBounding_hwwDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

// Top syst
cvswwsf[1]->cd();
histo_Top_CMS_MVATopBounding_hwwUp  ->Draw("hist");
histo_Top_CMS_MVATopBounding_hwwDown->Draw("hist,same");
histo_Top->Draw("e,same");

// stat. syst
cvswwsf[2]->cd();
histo_ttH_CMS_MVAttHStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_ttH_CMS_MVAttHStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwsf[3]->cd();
histo_ZH_CMS_MVAZHStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_ZH_CMS_MVAZHStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwsf[4]->cd();
histo_WH_CMS_MVAWHStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_WH_CMS_MVAWHStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwsf[5]->cd();
histo_qqH_CMS_MVAqqHStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_qqH_CMS_MVAqqHStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwsf[6]->cd();
histo_ggH_CMS_MVAggHStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_ggH_CMS_MVAggHStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwsf[7]->cd();
histo_qqWW_CMS_MVAqqWWStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_qqWW_CMS_MVAqqWWStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[8]->cd();
histo_ggWW_CMS_MVAggWWStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_ggWW_CMS_MVAggWWStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwsf[9]->cd();
histo_VV_CMS_MVAVVStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_VV_CMS_MVAVVStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwsf[10]->cd();
histo_Top_CMS_MVATopStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_Top_CMS_MVATopStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_Top->Draw("e,same");

cvswwsf[11]->cd();
histo_Zjets_CMS_MVAZjetsStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_Zjets_CMS_MVAZjetsStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_Zjets->Draw("e,same");

cvswwsf[12]->cd();
histo_Wjets_CMS_MVAWjetsStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_Wjets_CMS_MVAWjetsStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

cvswwsf[13]->cd();
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_Wgamma_CMS_MVAWgammaStatBounding_hwwsf_2jDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwsf[14]->cd();
histo_Ztt_CMS_MVAZttStatBounding_hwwsf_2jUp  ->Draw("hist");
histo_Ztt_CMS_MVAZttStatBounding_hwwsf_2jDown->Draw("hist,same");
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

// JES
cvswwsf[46]->cd();
histo_ttH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ttH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ttH->Draw("e,same");

cvswwsf[47]->cd();
histo_ZH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ZH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ZH->Draw("e,same");

cvswwsf[48]->cd();
histo_WH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_WH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_WH->Draw("e,same");

cvswwsf[49]->cd();
histo_qqH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_qqH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_qqH->Draw("e,same");

cvswwsf[50]->cd();
histo_ggH_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ggH_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ggH->Draw("e,same");

cvswwsf[51]->cd();
histo_qqWW_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_qqWW_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[52]->cd();
histo_ggWW_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_ggWW_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_ggWW->Draw("e,same");

cvswwsf[53]->cd();
histo_VV_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_VV_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_VV->Draw("e,same");

cvswwsf[54]->cd();
histo_Wgamma_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_Wgamma_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_Wgamma->Draw("e,same");

cvswwsf[55]->cd();
histo_Ztt_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_Ztt_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_Ztt->Draw("e,same");

// Wjets1 syst
cvswwsf[56]->cd();
histo_Wjets_CMS_MVAWBounding_hwwUp  ->Draw("hist");
histo_Wjets_CMS_MVAWBounding_hwwDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// Wjets2 syst
cvswwsf[57]->cd();
histo_Wjets_CMS_MVAWMCBounding_hwwUp  ->Draw("hist");
histo_Wjets_CMS_MVAWMCBounding_hwwDown->Draw("hist,same");
histo_Wjets->Draw("e,same");

// ggH syst
//cvswwsf[58]->cd();
//histo_ggH_CMS_MVAggHBoundingUp  ->Draw("hist");
//histo_ggH_CMS_MVAggHBoundingDown->Draw("hist,same");
//histo_ggH->Draw("e,same");

// WW syst2
cvswwsf[59]->cd();
histo_qqWW_CMS_MVAWWNLOBounding_hwwUp  ->Draw("hist");
histo_qqWW_CMS_MVAWWNLOBounding_hwwDown->Draw("hist,same");
histo_qqWW->Draw("e,same");

cvswwsf[60]->cd();
histo_Top_CMS_MVAJESBoundingUp  ->Draw("hist");
histo_Top_CMS_MVAJESBoundingDown->Draw("hist,same");
histo_Top->Draw("e,same");

// Zjets syst
cvswwsf[61]->cd();
histo_Zjets_CMS_MVAZBounding_hwwsf_2jUp  ->Draw("hist");
histo_Zjets_CMS_MVAZBounding_hwwsf_2jDown->Draw("hist,same");
histo_Zjets->Draw("e,same");

for(int i=0; i<=60; i++) cvswwof[i]->SaveAs(Form("plots/cvswwof_%d.pdf",i));
for(int i=0; i<=61; i++) cvswwsf[i]->SaveAs(Form("plots/cvswwsf_%d.pdf",i));
for(int i=0; i<=60; i++) cvswwof[i]->SaveAs(Form("plots/cvswwof_%d.eps",i));
for(int i=0; i<=61; i++) cvswwsf[i]->SaveAs(Form("plots/cvswwsf_%d.eps",i));
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
