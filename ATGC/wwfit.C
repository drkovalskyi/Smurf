#include "TFile.h"
#include "RooWorkspace.h"
#include "TProfile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TSystem.h"
#include "RooNDKeysPdf.h"
#include "RooATGCPdf.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "TH2.h"
#include "TMarker.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
// #include "RooLognormal.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooMsgService.h"
#include "../Core/SmurfTree.h"
#include "TRandom3.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooCachedPdf.h"
#include "RooKeysPdf.h"
#include "TStyle.h"
#include "RooDataHist.h"

double ww_expected = 774.0+45.9;
const double ww_uncertainty = sqrt(55.8*55.8+14.1*14.1);
const double top_expected = 121.6;
const double top_uncertainty = 23.4;
const double wjets_expected = 58.3+19.8;
const double wjets_uncertainty = sqrt(21.3*21.3+5.7*5.7);
const double wz_expected = 18.5;
const double wz_uncertainty = 1.9;
const double zz_expected = 10.6;
const double zz_uncertainty = 1.0;
const double dy_expected = 11.7;
const double dy_uncertainty = 6.8;

const bool addOverflowBin = true;

// const unsigned int Nbins = 18;
const unsigned int Nbins = 9;
// const unsigned int Nbins = 5;
const double minPt = 20;
const double maxPt = 200;

const double rangeX = 0.5;
const double rangeY = 0.5;
const double lumi = 4.9;

const int smooth = 0; // don't use it unless you have to

// const char* oldSampleSelection = "selected==1&&unique==1";
const char* oldSampleSelection = "selected==1&&unique==1&&dilpt>45";

RooRealVar* var_pt1;
RooRealVar* x_par;
RooRealVar* y_par;
RooRealVar* var_dummy;

RooRealVar* n_ww;
RooRealVar* n_top;
RooRealVar* n_wjets;
RooRealVar* n_wz;
RooRealVar* n_zz;
RooRealVar* n_dy;

RooATGCPdf* atgcPdf;
RooAbsPdf*  pdf_bkg;
RooDataSet* glb_data;
RooAbsPdf*  cpdf;
RooAbsPdf*  cSigPdf;
RooAbsPdf*  cBkgPdf;

RooAbsPdf* n_ww_con;
RooAbsPdf* n_top_con;
RooAbsPdf* n_wjets_con;
RooAbsPdf* n_wz_con;
RooAbsPdf* n_zz_con;
RooAbsPdf* n_dy_con;

RooDataSet* MakeDatasetWithOverflow(const RooDataSet* ds, const char* varname){
  RooDataSet* ds_out = dynamic_cast<RooDataSet*>(ds->emptyClone());
  assert(ds_out);
  for (Int_t i=0 ; i<ds->numEntries() ; i++){
    ds->get(i);
    RooArgSet argset(*ds->get(i));
    RooRealVar* var = dynamic_cast<RooRealVar*>(argset.find(varname));
    assert(var);
    var->setVal(TMath::Min(maxPt-0.1,var->getVal()));
    var->setRange(minPt,maxPt);
    ds_out->add(argset,ds->weight());
  }
  return ds_out;
}

class Sample{
  std::string m_file_name; // root file that contains reference point dataset
  std::string m_name;      // derived dataset name
  RooAbsData* m_dataset;
  double m_refX;
  double m_refY;
  TH1*   m_hist;        // simple histogram as a pdf
  TH1*   m_hist_keys;   // smoothed pdf
public:
  TH1* hist() {return m_hist;}
  TH1* hist_keys() {return m_hist_keys;}
  double refX() {return m_refX;}
  double refY() {return m_refY;}
  Sample(RooAbsData* dataset, double refx, double refy ){
    m_refX = refx;
    m_refY = refy;
    m_dataset = dataset;
    ((TTree*)m_dataset->tree())->Draw(Form("pt1>>h(%u,%f,%f)",Nbins,minPt,maxPt),"weight","e goff");
    m_hist = (TH1*)gDirectory->Get("h");
    m_hist->SetTitle(Form("%s: %0.2f, %s: %0.2f",x_par->GetTitle(),refx,y_par->GetTitle(),refy));
    m_hist->GetXaxis()->SetTitle("Pt, [GeV]");
    m_hist->SetDirectory(0);
    if (smooth) m_hist->Smooth(smooth);
    // for (Int_t i=0; i<=hist->GetNbinsX(); ++i)
    // hist->SetBinContent(i,hist->GetBinContent(i)+1e-5);
    m_hist->Scale(1/m_hist->Integral());
    // m_hist->GetYaxis()->SetRangeUser(0.001,0.25);
    m_hist->SetStats(kFALSE);
//     RooNDKeysPdf keys_pt1("keys_pt1","keys_pt1",*var_pt1,*((RooDataSet*)m_dataset),"am");
//     m_hist_keys = (TH1F*)keys_pt1.createHistogram("hist_keys",*var_pt1);
//     m_hist_keys->Scale(m_hist_keys->Integral());
//     m_hist_keys->SetLineColor(kRed);
//     m_hist_keys->SetLineWidth(2);
//     m_hist_keys->SetDirectory(0);
//     m_hist_keys->SetStats(kFALSE);
  }
  // old samples
  Sample( const char* file, const char* name, double refx, double refy ){
    m_file_name = file;
    m_name = name;
    m_refX = refx;
    m_refY = refy;
    TFile* f = TFile::Open(file);
    assert(f);
    RooDataSet* ds = (RooDataSet*)f->Get("ww");
    ds->SetName(name);
    if (addOverflowBin) ds = MakeDatasetWithOverflow(ds,"pt1");
    m_dataset = ds->reduce(*var_pt1,oldSampleSelection);
    ((TTree*)m_dataset->tree())->Draw(Form("pt1>>h(%u,%f,%f)",Nbins,minPt,maxPt),"","e goff");
    m_hist = (TH1*)gDirectory->Get("h");
    m_hist->SetTitle(Form("%s: %0.2f, %s: %0.2f",x_par->GetTitle(),refx,y_par->GetTitle(),refy));
    m_hist->GetXaxis()->SetTitle("Pt, [GeV]");
    m_hist->SetDirectory(0);
    if (smooth) m_hist->Smooth(smooth);
    // for (Int_t i=0; i<=hist->GetNbinsX(); ++i)
    // hist->SetBinContent(i,hist->GetBinContent(i)+1e-5);
    m_hist->Scale(1/m_hist->Integral());
    // m_hist->GetYaxis()->SetRangeUser(0.001,0.25);
    m_hist->SetStats(kFALSE);
//     RooNDKeysPdf keys_pt1("keys_pt1","keys_pt1",*var_pt1,*((RooDataSet*)m_dataset),"am");
//     m_hist_keys = (TH1F*)keys_pt1.createHistogram("hist_keys",*var_pt1);
//     m_hist_keys->Scale(m_hist_keys->Integral());
//     m_hist_keys->SetLineColor(kRed);
//     m_hist_keys->SetLineWidth(2);
//     m_hist_keys->SetDirectory(0);
//     m_hist_keys->SetStats(kFALSE);
  }
  Sample():m_dataset(0),m_refX(0),m_refY(0),m_hist(0),m_hist_keys(0){}
};

RooDataSet* MakeDataset(const char* file, const char* dataset_name, bool isData=false){
  
  RooRealVar pt1("pt1","pt1",minPt,maxPt);
  RooRealVar weightVar("weight","weight",1.0);
  RooArgSet variables(pt1,weightVar);
  RooDataSet* dataset = new RooDataSet(dataset_name, dataset_name, variables,
				       RooFit::WeightVar(weightVar));
  SmurfTree tree;
  tree.LoadTree(file);
  tree.InitTree();
  Long64_t nEntries = tree.tree_->GetEntries();
  Long64_t nSelected(0); 
  int i_permille_old = 0;
  for (Long64_t i = 0; i < nEntries; i++){
    int i_permille = (int)floor(1000 * i / float(nEntries));
    if (i_permille != i_permille_old) {
      // xterm magic from L. Vacavant and A. Cerri
      printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
	     "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
      fflush(stdout);
      i_permille_old = i_permille;
    }
    tree.tree_->GetEntry(i);
    if ( (tree.cuts_ & SmurfTree::FullSelection) != SmurfTree::FullSelection ||
	 (isData && (tree.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger )||
	 tree.lep2_.pt()<20 || 
	 //	 tree.lep1_.pt()>maxPt ||
	 tree.dilep_.pt()<45 ||
	 tree.njets_>0 ) continue;
    if (tree.type_==0||tree.type_==3){
      if ( tree.jet1_.pt()>15 && tree.dPhiDiLepJet1_>2.8798 ) continue;
      if ( tree.dilep_.mass()<20 ) continue;
    }
    if (addOverflowBin) 
      pt1.setVal(min(tree.lep1_.pt(),maxPt-0.1));
    else
      pt1.setVal(tree.lep1_.pt());
    double weight = 1.0;
    if (tree.dstype_ != SmurfTree::data)
      weight = tree.scale1fb_*lumi;
    dataset->add(variables,weight);
    nSelected++;
  }
  printf("Processed file %s\n\tSelected \t%d out of %d events\n",file,int(nSelected),int(nEntries));
  return dataset;
}

RooDataSet* MakeOldDataset(const char* file){
  TFile* f = TFile::Open(file);
  assert(f);
  RooDataSet* ds = (RooDataSet*)f->Get("ww");
  ds->SetName("ds_ww");
  if (addOverflowBin) ds = MakeDatasetWithOverflow(ds,"pt1");
  RooDataSet* dataset = dynamic_cast<RooDataSet*>(ds->reduce(*var_pt1,oldSampleSelection));
  assert(dataset);
  return dataset;
}

RooAbsPdf* MakePdfFromDataset(RooDataSet* data, RooAbsReal* variable){
  // if number of events is small, make a keys pdf. Otherwise use a histogram
  
  std::string name(Form("pdf_%s",data->GetName()));
  RooAbsPdf* outPdf(0);
  // if ( data->numEntries() < 5000 ){
  if ( data->numEntries() < 0 ){
    outPdf = new RooKeysPdf(name.c_str(),name.c_str(),*variable,*data, RooKeysPdf::MirrorBoth);
  } else {
    // Create a binned dataset 
    RooDataHist* hist = new RooDataHist("hist","hist",*var_pt1,*data) ;
    // Make a pdf
    outPdf = new RooHistPdf(name.c_str(),name.c_str(),*variable, *hist);
  }
  return outPdf;
}

void DrawPdf(RooAbsPdf* ipdf, const char* iname, const char* title, 
	     EColor fillColor = kMagenta, EColor markerColor=kBlue, const char* drawOptions =""){
  std::string name(Form("h_%s",iname));
  TH1* hpdf = ipdf->createHistogram(name.c_str(),*var_pt1);
  hpdf->SetTitle(iname);
  hpdf->GetXaxis()->SetTitle(title);
  hpdf->Scale(hpdf->Integral());
  hpdf->SetFillColor(fillColor);
  hpdf->SetMarkerColor(markerColor);
  hpdf->SetMarkerStyle(20);
  hpdf->SetLineWidth(2);
  hpdf->Draw(drawOptions);
}

// ============================================================================== //

Sample samples[11]; 

// http://root.cern.ch/root/html/RooLognormal.html
// http://en.wikipedia.org/wiki/Log-normal_distribution
double meanLogNormal(double mean, double sigma){
  return 1/sqrt(mean*mean+sigma*sigma);
}
double sigmaLogNormal(double mean, double sigma){
  return exp(sqrt(log(1+sigma*sigma/mean/mean)));
}

void setDefaults()
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  var_pt1   = new RooRealVar("pt1","pt1",minPt,maxPt);
  var_dummy = new RooRealVar("var_dummy","var_dummy",0);
  var_pt1->setBins(Nbins);

  n_ww    = new RooRealVar("n_ww","n_ww",       ww_expected, ww_expected/4, ww_expected*4);
  n_top   = new RooRealVar("n_top","n_top",     top_expected, top_expected/4, top_expected*4);
  n_wjets = new RooRealVar("n_wjets","n_wjets", wjets_expected, wjets_expected/4, wjets_expected*4);
  n_wz    = new RooRealVar("n_wz","n_wz",       wz_expected, wz_expected/4, wz_expected*4);
  n_zz    = new RooRealVar("n_zz","n_zz",       zz_expected, zz_expected/4, zz_expected*4);
  n_dy    = new RooRealVar("n_dy","n_dy",       dy_expected, dy_expected/4, dy_expected*4);

  n_ww_con    = new RooGaussian("n_ww_con",    "WW uncertainty",    *n_ww,    RooFit::RooConst(ww_expected),RooFit::RooConst(ww_uncertainty));
  n_top_con   = new RooGaussian("n_top_con",   "Top uncertainty",   *n_top,   RooFit::RooConst(top_expected),RooFit::RooConst(top_uncertainty));
  n_wjets_con = new RooGaussian("n_wjets_con", "Wjets uncertainty", *n_wjets, RooFit::RooConst(wjets_expected),RooFit::RooConst(wjets_uncertainty));
  n_wz_con    = new RooGaussian("n_wz_con",    "WZ uncertainty",    *n_wz,    RooFit::RooConst(wz_expected),RooFit::RooConst(wz_uncertainty));
  n_zz_con    = new RooGaussian("n_zz_con",    "ZZ uncertainty",    *n_zz,    RooFit::RooConst(zz_expected),RooFit::RooConst(zz_uncertainty));
  n_dy_con    = new RooGaussian("n_dy_con",    "DY uncertainty",    *n_dy,    RooFit::RooConst(dy_expected),RooFit::RooConst(dy_uncertainty));

}

void setSigPdf_LZ_GZ()
{
  x_par = new RooRealVar("x_par","#lambda_{Z}",-1,1);
  y_par = new RooRealVar("y_par","#Delta g^{Z}_{1}",-1,1);
  double norm = ww_expected/3.09163;
  
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(3,3);

  atgcPdf = new RooATGCPdf("pdf", "pdf", *var_pt1, *n_ww, *x_par, *y_par);
  cSigPdf = new RooProdPdf("cSigPdf","model with constraint",RooArgSet(*atgcPdf,*n_ww_con)) ;
  
  Int_t i=0;
  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz100_lz0_gz1000000.root", "sm_sm",0,0);

  // samples[i] = Sample(MakeDataset("smurf/qqww.root","ds_ww"),0,0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*3.09163, norm*0.0547394), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz175_lz0_gz1750000.root", "sm_p",0,0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*6.12261, 0.112729*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz25_lz0_gz250000.root", "sm_m",0,-0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*6.01277, 0.112194*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root", "p_sm", 0.5, 0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*5.53635, 0.103594*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz175_lz50_gz1750000.root", "p_p",0.5,0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*9.43854, 0.178215*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz25_lz50_gz250000.root", "p_m",0.5,-0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*7.22219, 0.13752*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz1000000.root", "m_sm",-0.5,0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*5.20204, 0.0981948*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz175_lzm50_gz1750000.root", "m_p",-0.5,0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*6.89533, 0.131155*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz25_lzm50_gz25000.root", "m_m",-0.5,-0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*9.72764, 0.183288*norm), samples[i].hist());
  i++;			    

    /*
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*547.541, norm*1.947), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz175_lz0_gz1750000.root", "sm_p",0,0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1284.979, 2.875*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz25_lz0_gz250000.root", "sm_m",0,-0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1370.808, 2.905*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root", "p_sm", 0.5, 0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1223.832, 3.414*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz175_lz50_gz1750000.root", "p_p",0.5,0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*2245.131, 5.315*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz25_lz50_gz250000.root", "p_m",0.5,-0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1766.015, 4.324*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz1000000.root", "m_sm",-0.5,0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1193.115, 3.026*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz175_lzm50_gz1750000.root", "m_p",-0.5,0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1644.440, 4.276*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz25_lzm50_gz25000.root", "m_m",-0.5,-0.75);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*2295.906, 4.985*norm), samples[i].hist());
  i++;
    */

  atgcPdf->build();

  c1->cd(1);
  x_par->setVal(0); y_par->setVal(0);
  DrawPdf(atgcPdf,"h_sm_sm","",kWhite,kBlue,"same p0");

  c1->cd(2);
  x_par->setVal(0); y_par->setVal(0.75);
  DrawPdf(atgcPdf,"h_sm_p","",kWhite,kBlue,"same p0");

  c1->cd(3);
  x_par->setVal(0); y_par->setVal(-0.75);
  DrawPdf(atgcPdf,"h_sm_m","",kWhite,kBlue,"same p0");

  c1->cd(4);
  x_par->setVal(0.5); y_par->setVal(0);
  DrawPdf(atgcPdf,"h_p_sm","",kWhite,kBlue,"same p0");

  c1->cd(5);
  x_par->setVal(0.5); y_par->setVal(0.75);
  DrawPdf(atgcPdf,"h_p_p","",kWhite,kBlue,"same p0");

  c1->cd(6);
  x_par->setVal(0.5); y_par->setVal(-0.75);
  DrawPdf(atgcPdf,"h_p_m","",kWhite,kBlue,"same p0");

  c1->cd(7);
  x_par->setVal(-0.5); y_par->setVal(0);
  DrawPdf(atgcPdf,"h_m_sm","",kWhite,kBlue,"same p0");

  c1->cd(8);
  x_par->setVal(-0.5); y_par->setVal(0.75);
  DrawPdf(atgcPdf,"h_m_p","",kWhite,kBlue,"same p0");
  
  c1->cd(9);
  x_par->setVal(-0.5); y_par->setVal(-0.75);
  DrawPdf(atgcPdf,"h_m_m","",kWhite,kBlue,"same p0");
  c1->Print("pdfs.pdf");

}

void setSigPdf_LZ_KG()
{
  x_par = new RooRealVar("x_par","#lambda_{Z}",-1,1);
  y_par = new RooRealVar("y_par","#Delta#kappa_{#gamma}",-1,1);
  double norm = ww_expected/547.541;
  
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(3,3);

  atgcPdf = new RooATGCPdf("pdf", "pdf", *var_pt1, *n_ww, *x_par, *y_par);
  cSigPdf = new RooProdPdf("cSigPdf","model with constraint",RooArgSet(*atgcPdf,*n_ww_con)) ;
  
  Int_t i=0;
  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz100_lz0_gz1000000.root", "sm_sm",0,0);
  // samples[i] = Sample(MakeDataset("smurf/qqww.root","ds_ww"),0,0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*547.541, norm*1.947), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg170_lg0_gg100_kz779311_lz0_gz100.root", "sm_p",0,0.70);
  // samples[i] = Sample(MakeDataset("/afs/cern.ch/work/d/dmytro/atgc/mcfm/WWqqbr_tota_mstw8nl_80__80__kg170_lg0_kZ7795_lZ0_gZ100_smurf.root","sm_p"),0,0.70);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*617.147, 1.604*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg30_lg0_gg100_kz1220689_lz0_gz100.root", "sm_m",0,-0.70);
  // samples[i] = Sample(MakeDataset("/afs/cern.ch/work/d/dmytro/atgc/mcfm/WWqqbr_tota_mstw8nl_80__80__kg30_lg0_kZ12205_lZ0_gZ100_smurf.root","sm_p"),0,-0.70);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*633.939 , 1.671*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root", "p_sm", 0.5, 0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1223.832, 3.414*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg170_lg50_gg100_kz779311_lz50_gz100.root", "p_p",0.5,0.70);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1295.486, 3.224*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg30_lg50_gg100_kz1220689_lz50_gz100.root", "p_m",0.5,-0.70);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1319.705, 2.920*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz1000000.root", "m_sm",-0.5,0);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1193.115, 3.026*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg170_lgm50_gg100_kz779311_lzm50_gz100.root", "m_p",-0.5,0.70);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1267.807, 2.971*norm), samples[i].hist());
  i++;			    

  c1->cd(i+1);
  samples[i] = Sample("samples/processed_data_job_WW_1j_kg30_lgm50_gg100_kz1220689_lzm50_gz100.root", "m_m",-0.5,-0.70);
  samples[i].hist()->Draw();
  atgcPdf->addPoint(Measurement(samples[i].refX(), samples[i].refY(), norm*1279.312, 2.999*norm), samples[i].hist());
  i++;			    
  atgcPdf->build();
  
  c1->cd(1);
  x_par->setVal(0); y_par->setVal(0);
  DrawPdf(atgcPdf,"h_sm_sm","",kWhite,kBlue,"same p0");

  c1->cd(2);
  x_par->setVal(0); y_par->setVal(0.75);
  DrawPdf(atgcPdf,"h_sm_p","",kWhite,kBlue,"same p0");

  c1->cd(3);
  x_par->setVal(0); y_par->setVal(-0.75);
  DrawPdf(atgcPdf,"h_sm_m","",kWhite,kBlue,"same p0");

  c1->cd(4);
  x_par->setVal(0.5); y_par->setVal(0);
  DrawPdf(atgcPdf,"h_p_sm","",kWhite,kBlue,"same p0");

  c1->cd(5);
  x_par->setVal(0.5); y_par->setVal(0.75);
  DrawPdf(atgcPdf,"h_p_p","",kWhite,kBlue,"same p0");

  c1->cd(6);
  x_par->setVal(0.5); y_par->setVal(-0.75);
  DrawPdf(atgcPdf,"h_p_m","",kWhite,kBlue,"same p0");

  c1->cd(7);
  x_par->setVal(-0.5); y_par->setVal(0);
  DrawPdf(atgcPdf,"h_m_sm","",kWhite,kBlue,"same p0");

  c1->cd(8);
  x_par->setVal(-0.5); y_par->setVal(0.75);
  DrawPdf(atgcPdf,"h_m_p","",kWhite,kBlue,"same p0");
  
  c1->cd(9);
  x_par->setVal(-0.5); y_par->setVal(-0.75);
  DrawPdf(atgcPdf,"h_m_m","",kWhite,kBlue,"same p0");
  c1->Print("pdfs2.pdf");
}

void setBkgPdf()
{
  TCanvas* c9 = new TCanvas("c9","c9",400,800);
  c9->Divide(2,4);

  c9->cd(1);
  RooDataSet* ds_ww_pt1 = MakeDataset("smurf/qqww.root","ds_ww");
  // RooDataSet* ds_ww_pt1 = MakeDataset("smurf/ww_mcnlo.root","ds_ww");
  RooAbsPdf*  pdf_ww    = MakePdfFromDataset(ds_ww_pt1,var_pt1);
  DrawPdf(pdf_ww, "WW","Leading lepton pt, [GeV]");

  c9->cd(2);
  RooDataSet* ds_wjets_pt1 = MakeDataset("smurf/wjets.root","ds_wjets");
  RooAbsPdf*  pdf_wjets    = MakePdfFromDataset(ds_wjets_pt1,var_pt1);
  DrawPdf(pdf_wjets, "Wjets","Leading lepton pt, [GeV]");

  c9->cd(3);
  RooDataSet* ds_ttbar_pt1 = MakeDataset("smurf/ttbar.root","ds_ttbar");
  RooAbsPdf*  pdf_ttbar    = MakePdfFromDataset(ds_ttbar_pt1,var_pt1);
  DrawPdf(pdf_ttbar, "TTbar","Leading lepton pt, [GeV]");

  c9->cd(4);
  RooDataSet* ds_tw_pt1 = MakeDataset("smurf/tw.root","ds_tw");
  RooAbsPdf*  pdf_tw    = MakePdfFromDataset(ds_tw_pt1,var_pt1);
  DrawPdf(pdf_tw, "tW","Leading lepton pt, [GeV]");
  
  c9->cd(5);
  RooDataSet* ds_wz_pt1 = MakeDataset("smurf/wz.root","ds_wz");
  RooAbsPdf*  pdf_wz    = MakePdfFromDataset(ds_wz_pt1,var_pt1);
  DrawPdf(pdf_wz, "WZ","Leading lepton pt, [GeV]");

  c9->cd(6);
  RooDataSet* ds_zz_pt1 = MakeDataset("smurf/zz_py.root","ds_zz");
  RooAbsPdf*  pdf_zz    = MakePdfFromDataset(ds_zz_pt1,var_pt1);
  DrawPdf(pdf_zz, "ZZ","Leading lepton pt, [GeV]");

  c9->cd(7);
  RooDataSet* ds_dyee_pt1 = MakeDataset("smurf/dyee.root","ds_dyee");
  // RooAbsPdf*  pdf_dyee    = MakePdfFromDataset(ds_dyee_pt1,var_pt1);
  // DrawPdf(pdf_dyee, "DYee","Leading lepton pt, [GeV]");

  c9->cd(8);
  RooDataSet* ds_dymm_pt1 = MakeDataset("smurf/dymm.root","ds_dyee");
  // RooAbsPdf*  pdf_dymm    = MakePdfFromDataset(ds_dymm_pt1,var_pt1);
  // DrawPdf(pdf_dymm, "DYmm","Leading lepton pt, [GeV]");

  c9->cd(7);
  RooDataSet* ds_dy_pt1 = new RooDataSet("ds_dy","ds_dy",ds_dyee_pt1,*var_pt1);
  ds_dy_pt1->append(*ds_dymm_pt1);
  RooAbsPdf*  pdf_dy    = MakePdfFromDataset(ds_dy_pt1,var_pt1);
  DrawPdf(pdf_dy, "DY","Leading lepton pt, [GeV]");

  c9->cd(8);
  // extended pdfs
  RooAbsPdf* epdf_wjets = new RooExtendPdf("epdf_wjets","epdf_wjets",*pdf_wjets,*n_wjets);
  RooAbsPdf* epdf_top = new RooExtendPdf("epdf_top","epdf_top",*pdf_ttbar,*n_top);
  RooAbsPdf* epdf_wz = new RooExtendPdf("epdf_wz","epdf_wz",*pdf_wz,*n_wz);
  RooAbsPdf* epdf_zz = new RooExtendPdf("epdf_zz","epdf_zz",*pdf_zz,*n_zz);
  RooAbsPdf* epdf_dy = new RooExtendPdf("epdf_dy","epdf_dy",*pdf_dy,*n_dy);
  pdf_bkg = new RooAddPdf("pdf_bkg","pdf_bkg",RooArgList(*epdf_wjets,*epdf_top, *epdf_wz, *epdf_zz, *epdf_dy));
  cBkgPdf = new RooProdPdf("cBkgPdf","model with constraint",RooArgSet(*pdf_bkg,*n_top_con,*n_wjets_con,*n_wz_con,*n_zz_con,*n_dy_con)) ;
  DrawPdf(cBkgPdf, "Background PDF","Leading lepton pt, [GeV]");
}

void setOldBkgPdf()
{
  TCanvas* c9 = new TCanvas("c9","c9",600,900);
  c9->Divide(2,3);

  TFile* f = TFile::Open("samples/processed_data_final.root");
  assert(f);
  
  c9->cd(1);
  // ww
  RooAbsData* ds_ww = (RooAbsData*)f->Get("ww");
  ds_ww->SetName("ds_ww");
  RooAbsData* ds_ww_pt1 = ds_ww->reduce(*var_pt1,oldSampleSelection);
  RooNDKeysPdf* pdf_ww = new RooNDKeysPdf("pdf_ww","pdf_ww",*var_pt1,*((RooDataSet*)ds_ww_pt1),"am");
  TH1F* hpdf_ww = (TH1F*)pdf_ww->createHistogram("hpdf_ww",*var_pt1);
  hpdf_ww->SetTitle("WW");
  hpdf_ww->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_ww->Scale(hpdf_ww->Integral());
  hpdf_ww->Draw();

  c9->cd(2);
  // wjets
  RooAbsData* ds_wjets = (RooAbsData*)f->Get("wjets");
  ds_wjets->SetName("ds_wjets");
  RooAbsData* ds_wjets_pt1 = ds_wjets->reduce(*var_pt1,oldSampleSelection);
  RooNDKeysPdf* pdf_wjets = new RooNDKeysPdf("pdf_wjets","pdf_wjets",*var_pt1,*((RooDataSet*)ds_wjets_pt1),"am");
  TH1F* hpdf_wjets = (TH1F*)pdf_wjets->createHistogram("hpdf_wjets",*var_pt1);
  hpdf_wjets->SetTitle("WJets");
  hpdf_wjets->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_wjets->Scale(hpdf_wjets->Integral());
  hpdf_wjets->Draw();

  c9->cd(3);
  // ttbar
  RooAbsData* ds_ttbar = (RooAbsData*)f->Get("ttbar");
  ds_ttbar->SetName("ds_ttbar");
  RooAbsData* ds_ttbar_pt1 = ds_ttbar->reduce(*var_pt1,oldSampleSelection);
  RooNDKeysPdf* pdf_ttbar = new RooNDKeysPdf("pdf_ttbar","pdf_ttbar",*var_pt1,*((RooDataSet*)ds_ttbar_pt1),"am");
  TH1F* hpdf_ttbar = (TH1F*)pdf_ttbar->createHistogram("hpdf_ttbar",*var_pt1);
  hpdf_ttbar->SetTitle("TTbar");
  hpdf_ttbar->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_ttbar->Scale(hpdf_ttbar->Integral());
  hpdf_ttbar->Draw();

  c9->cd(4);
  // tW
  RooAbsData* ds_tw = (RooAbsData*)f->Get("tw");
  ds_tw->SetName("ds_tw");
  RooAbsData* ds_tw_pt1 = ds_tw->reduce(*var_pt1,oldSampleSelection);
  RooNDKeysPdf* pdf_tw = new RooNDKeysPdf("pdf_tw","pdf_tw",*var_pt1,*((RooDataSet*)ds_tw_pt1),"am");
  TH1F* hpdf_tw = (TH1F*)pdf_tw->createHistogram("hpdf_tw",*var_pt1);
  hpdf_tw->SetTitle("tW");
  hpdf_tw->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_tw->Scale(hpdf_tw->Integral());
  hpdf_tw->Draw();
  
  c9->cd(5);
  // wz
  RooAbsData* ds_wz = (RooAbsData*)f->Get("wz");
  ds_wz->SetName("ds_wz");
  RooAbsData* ds_wz_pt1 = ds_wz->reduce(*var_pt1,oldSampleSelection);
  RooNDKeysPdf* pdf_wz = new RooNDKeysPdf("pdf_wz","pdf_wz",*var_pt1,*((RooDataSet*)ds_wz_pt1),"am");
  TH1F* hpdf_wz = (TH1F*)pdf_wz->createHistogram("hpdf_wz",*var_pt1);
  hpdf_wz->SetTitle("WZ");
  hpdf_wz->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_wz->Scale(hpdf_wz->Integral());
  hpdf_wz->Draw();
  
  c9->cd(6);
  // zz
  RooAbsData* ds_zz = (RooAbsData*)f->Get("zz");
  ds_zz->SetName("ds_zz");
  RooAbsData* ds_zz_pt1 = ds_zz->reduce(*var_pt1,oldSampleSelection);
  RooNDKeysPdf* pdf_zz = new RooNDKeysPdf("pdf_zz","pdf_zz",*var_pt1,*((RooDataSet*)ds_zz_pt1),"am");
  TH1F* hpdf_zz = (TH1F*)pdf_zz->createHistogram("hpdf_zz",*var_pt1);
  hpdf_zz->SetTitle("ZZ");
  hpdf_zz->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_zz->Scale(hpdf_zz->Integral());
  hpdf_zz->Draw();

  // extended pdfs
  RooAbsPdf* epdf_wjets = new RooExtendPdf("epdf_wjets","epdf_wjets",*pdf_wjets,*n_wjets);
  RooAbsPdf* epdf_top = new RooExtendPdf("epdf_top","epdf_top",*pdf_ttbar,*n_top);
  RooAbsPdf* epdf_wz = new RooExtendPdf("epdf_wz","epdf_wz",*pdf_wz,*n_wz);
  RooAbsPdf* epdf_zz = new RooExtendPdf("epdf_zz","epdf_zz",*pdf_zz,*n_zz);
  pdf_bkg = new RooAddPdf("pdf_bkg","pdf_bkg",RooArgList(*epdf_wjets,*epdf_top, *epdf_wz, *epdf_zz));
  // pdf_bkg = new RooAddPdf("pdf_bkg","pdf_bkg",RooArgList(*epdf_wjets,*epdf_top, *epdf_wz));
  cBkgPdf = new RooProdPdf("cBkgPdf","model with constraint",RooArgSet(*pdf_bkg,*n_top_con,*n_wjets_con,*n_wz_con,*n_zz_con)) ;
  // cBkgPdf = new RooProdPdf("cBkgPdf","model with constraint",RooArgSet(*pdf_bkg,*n_top_con,*n_wjets_con,*n_wz_con)) ;
}


// not used
// need to restore
/*
void sensitivity()
{
  TFile *f = TFile::Open("atgc.root");
  assert(f);
  RooWorkspace *ww = (RooWorkspace*)f->Get("ww");
  assert(ww);
  
  const unsigned int nPoints = 21;
  std::vector<double> sum(nPoints,0);
  std::vector<double> sum2(nPoints,0);
    
  // loop over WW dataset and make small samples
  // perform Likelihood difference calculation for each
  // value of the parameter
  unsigned int n(0);
  TTree* iTree = const_cast<TTree*>(((RooDataSet*)ww->data("ds_ww_pt1"))->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  unsigned int firstEntry = 0;
  y_par->setVal(0);
  data=0;
  unsigned int nEvents(ww_expected);
  while ( firstEntry <= size - nEvents ){
    n++;
    TTree* tree = iTree->CopyTree("","",nEvents,firstEntry);
    assert(tree);
    firstEntry+=nEvents;
    if (data) delete data;
    data = new RooDataSet("data","data",tree, *ww->var("pt1"));
    // continue;
    for ( unsigned int i = 0; i<nPoints; ++i ){
      x_par->setVal(0);
      RooAbsReal* nll = cSigPdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(*n_ww));
      nll->addServer(*var_dummy);
      double nll_sm(nll->getVal());
      delete nll;
      x_par->setVal(-0.5+i*0.05);
      nll = cSigPdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(*n_ww));
      nll->addServer(*var_dummy);
      double nll_atgc(nll->getVal());
      delete nll;
      double sig = sqrt(2.*fabs(nll_sm-nll_atgc));
      sum[i]+=sig;
      sum2[i]+=sig*sig;
    }
  }
  // return;

  // p->Draw();
  if (n==0) return;
  // show results
  TCanvas* c2 = new TCanvas("c2","Significance",500,500);
  //c1->SetFillColor(42);
  c2->SetGrid();
  Float_t xval[nPoints];
  Float_t yval[nPoints];
  Float_t yerr[nPoints];
  for ( unsigned int i = 0; i<nPoints; ++i ){
    xval[i] = -0.5+i*0.05;
    yval[i] = sum[i]/n;
    yerr[i] = sqrt(sum2[i]/n-pow(sum[i]/n,2));
  }
  TGraphErrors* gr = new TGraphErrors(nPoints,xval,yval,(Float_t*)0,yerr);
  // gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");

  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  y_par->setConstant(0);
  x_par->setConstant(0);

  // fit 
//   {
//     RooAbsReal* nll = cSigPdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(*n_ww));
//     RooMinuit m(*nll);
//     m.migrad();
//     m.hesse();
//     RooPlot* p1 = m.contour(*x_par,*y_par);
//     m.migrad();
//     m.hesse();
//     TCanvas* c4 = new TCanvas("c4","c4",500,500);
//     c4->SetGrid();
//     p1->Draw();
//   }

  // fit the reference point samples
  TCanvas* c5 = new TCanvas("c5","c5",800,800);
  c5->Divide(3,3);
  TCanvas* c6 = new TCanvas("c6","c6",800,800);
  c6->Divide(3,3);
  for (unsigned int i=0; i<9; ++i){
    x_par->setVal(samples[i].refX());
    y_par->setVal(samples[i].refY());
    Int_t n = cSigPdf->expectedEvents(RooArgSet());
    cout << "i: " << i << " \tExpected number of events: " << n << endl;
    TTree* iTree = const_cast<TTree*>(((RooDataSet*)samples[i].dataset)->tree());
    assert(iTree);
    TTree* tree = iTree->CopyTree("","",n,0);
    assert(tree);
    RooDataSet*ds = new RooDataSet("ds","ds",tree, *ww->var("pt1"));
    c5->cd(i+1);
    gPad->SetGrid();
    // x_par->setVal(0);
    // y_par->setVal(0);
    RooAbsReal* nll = cSigPdf->createNLL(*ds,RooFit::Extended(),RooFit::Constrain(*n_ww));
    double minNLL = nll->getVal();
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    TMarker* markerFit = new TMarker(x_par->getVal(),y_par->getVal(),3);
    RooPlot* p = m.contour(*x_par,*y_par);
    m.migrad();
    m.hesse();
    p->Draw();
    c6->cd(i+1);
    const Int_t nBins = 200;
    RooMsgService::instance().setStreamStatus(1,0);
    TH2F* nll2d = new TH2F("nll2d",Form("#sqrt{2 #Delta logL}: %s: %0.2f, %s: %0.2f",x_par->GetTitle(),samples[i].refX(),x_par->GetTitle(),samples[i].refY()),
			   nBins,-1,1,nBins,-1,1);
    nll2d->GetXaxis()->SetTitle(x_par->GetTitle());
    nll2d->GetYaxis()->SetTitle(y_par->GetTitle());
    nll2d->SetDirectory(0);
    for ( Int_t xi = 1; xi<=nBins; ++xi ){
      cout << ".";
      cout.flush();
      for ( Int_t yi = 1; yi<=nBins; ++yi )
	{
	  x_par->setVal(-0.95+2.0*(xi-1)/nBins);
	  y_par->setVal(-0.95+2.0*(yi-1)/nBins);
	  RooAbsReal* nll = cSigPdf->createNLL(*ds,RooFit::Extended(),RooFit::Constrain(*n_ww));
	  nll->addServer(*var_dummy);
	  double sig = sqrt(2*fabs(nll->getVal()-minNLL));
	  nll2d->SetBinContent(xi,yi,sig);
	  delete nll;
	}
    }
    nll2d->Draw("colz");
    RooMsgService::instance().setStreamStatus(1,1);
    TMarker* marker = new TMarker(samples[i].refX(),samples[i].refY(),20);
    // marker->SetColor(kBlue);
    marker->Draw();
    markerFit->Draw();
  }
  
}
*/

void ww1DFits(const char* file = "smurf/qqww.root", bool smurfFormat=true, int N=-1)
{
  RooAbsData* dataset(0);
  if (!smurfFormat){
    dataset = MakeOldDataset(file);
  } else {
    dataset = MakeDataset(file,"qqww");
  }
  //// RooAbsData* dataset = MakeDataset("samples/WW1j_n0_kg94_lg0_gg100_kZ100_lZ0_gZ100_fastsim429_v5.root","qqww");

  // loop over WW dataset and make small samples
  // perform Likelihood difference calculation for each
  // value of the parameter
  // TTree* iTree = const_cast<TTree*>(((RooDataSet*)ww->data("ds_ww_pt1"))->tree());
  TTree* iTree = const_cast<TTree*>(((RooDataSet*)dataset)->tree());
  assert(iTree);
  iTree->Draw(Form("pt1>>h(%u,%f,%f)",Nbins,minPt,maxPt),"weight","goff");
  TH1* hist = (TH1*)gDirectory->Get("h");

  unsigned int size = iTree->GetEntries();
  double meanN = hist->Integral(0,Nbins+1);
  if (N>0) meanN = N;

  printf("Average number of events: %0.0f\t Total number of events: %d\n",
	 meanN, int(iTree->GetEntries()) );
  unsigned int firstEntry = 0;
  y_par->setVal(0);
  glb_data=0;
  TRandom3 generator;

  TH1F* hx = new TH1F("hx","Fit of WW SM events",40,-.2,.2);
  hx->GetXaxis()->SetTitle(x_par->getTitle());
  hx->SetFillColor(kMagenta);
  TH1F* hxpull = new TH1F("hxpull","Fit of WW SM events",40,-5,5);
  hxpull->GetXaxis()->SetTitle(Form("%s/#sigma",x_par->getTitle().Data()));
  hxpull->SetFillColor(kMagenta);
  TH1F* hy = new TH1F("hy","Fit of WW SM events",40,-0.3,0.3);
  hy->GetXaxis()->SetTitle(y_par->getTitle());
  hy->SetFillColor(kMagenta);
  TH1F* hypull = new TH1F("hypull","Fit of WW SM events",40,-5,5);
  hypull->GetXaxis()->SetTitle(Form("%s/#sigma",y_par->getTitle().Data()));
  hypull->SetFillColor(kMagenta);
  while ( firstEntry < size ){
    int n = generator.Poisson(meanN);
    if ( firstEntry+n > size ) break;
    TTree* tree = iTree->CopyTree("","",n,firstEntry);
    assert(tree);
    firstEntry+=n;
    if (glb_data) delete glb_data;
    glb_data = new RooDataSet("data","data",tree, *var_pt1);
    
    x_par->setVal(0);
    x_par->setConstant(0);
    y_par->setVal(0);
    y_par->setConstant(1);
    // assert(atgcPdf->fitTo(*glb_data,RooFit::InitialHesse(true),RooFit::Strategy(2),RooFit::Minos(),RooFit::Save())->status()==0);
    atgcPdf->fitTo(*glb_data,RooFit::Minos());
    hx->Fill(x_par->getVal());
    hxpull->Fill(x_par->getVal()/x_par->getError());
//     if(x_par->getVal()>0)
//       hxpull->Fill(x_par->getVal()/fabs(x_par->getErrorLo()));
//     else
//       hxpull->Fill(x_par->getVal()/fabs(x_par->getErrorHi()));

    x_par->setVal(0);
    x_par->setConstant(1);
    y_par->setVal(0);
    y_par->setConstant(0);
    // assert(atgcPdf->fitTo(*glb_data,RooFit::InitialHesse(true),RooFit::Strategy(2),RooFit::Minos(),RooFit::Save())->status()==0);
    atgcPdf->fitTo(*glb_data);
    hy->Fill(y_par->getVal());
    hypull->Fill(y_par->getVal()/y_par->getError());
//     if(y_par->getVal()>0)
//       hypull->Fill(y_par->getVal()/fabs(y_par->getErrorLo()));
//     else
//       hypull->Fill(y_par->getVal()/fabs(y_par->getErrorHi()));
  }

  // show results
  TCanvas* c7 = new TCanvas("c7","Fit results",800,400);
  c7->Divide(2,1);
  c7->cd(1);
  hx->Draw();
  c7->cd(2);
  hy->Draw();
  TCanvas* c71 = new TCanvas("c71","Fit results",800,400);
  c71->Divide(2,1);
  c71->cd(1);
  c71->cd(1);
  hxpull->Draw();
  c71->cd(2);
  hypull->Draw();
}

TH1F* wwATGC1DFit(const char* file, const char* name, double lz, double dkz)
{
  x_par->setRange(-1.5,1.5);
  y_par->setRange(-1.5,1.5);
  assert(lz==0||dkz==0);
  
  RooAbsData* dataset = MakeOldDataset(file);
  
  // loop over WW dataset and make small samples
  // perform Likelihood difference calculation for each
  // value of the parameter
  unsigned int n(0);
  TTree* iTree = const_cast<TTree*>(((RooDataSet*)dataset)->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  unsigned int firstEntry = 0;
  y_par->setVal(0);
  glb_data=0;
  // TRandom::Poisson

  // TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW with anomalous couplings",40,-1,1);
  TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW with anomalous couplings",60,0,1.5);
  h->SetDirectory(0);
  unsigned int nEvents(ww_expected);
  while ( firstEntry <= size - nEvents ){
    n++;
    TTree* tree = iTree->CopyTree("","",nEvents,firstEntry);
    assert(tree);
    firstEntry+=nEvents;
    if (glb_data) delete glb_data;
    glb_data = new RooDataSet("data","data",tree, *var_pt1);
    cout << "Num entries: " << glb_data->numEntries() << endl;

    x_par->setVal(lz);
    y_par->setVal(dkz);
    if ( lz==0 ){
      x_par->setConstant(1);
      y_par->setConstant(0);
    } else {
      x_par->setConstant(0);
      y_par->setConstant(1);
    }      
    atgcPdf->fitTo(*glb_data);
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;
    if ( lz==0 )
      h->Fill(fabs(y_par->getVal()));
    else
      h->Fill(fabs(x_par->getVal()));
  }
  h->SetFillColor(kMagenta);
  return h;
}

void wwATGC1DFits()
{
  TCanvas* c8 = new TCanvas("c8","Fit results",800,800);
  c8->Divide(2,2);
  
  c8->cd(1);
  // TH1F* h1 = wwATGC1DFit("samples/processed_data_WW_1j_kg100_lg0_gg100_kz175_lz0_gz175_fastsim386_v1.root",
  TH1F* h1 = wwATGC1DFit("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz175_lz0_gz1750000.root",
			 "sm_p",0,0.75);
  h1->SetTitle("#lambda_{Z}=0, #Delta g^{Z}_{1}=0.75");
  h1->GetXaxis()->SetTitle("|#Delta g^{Z}_{1}|");
  h1->Draw();

  c8->cd(2);
  // TH1F* h2 = wwATGC1DFit("samples/processed_data_WW_1j_kg100_lg0_gg100_kz25_lz0_gz25_fastsim386_v1.root",
  TH1F* h2 = wwATGC1DFit("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz25_lz0_gz250000.root",
			 "sm_m",0,-0.75);
  h2->SetTitle("#lambda_{Z}=0, #Delta g^{Z}_{1}=-0.75");
  h2->GetXaxis()->SetTitle("|#Delta g^{Z}_{1}|");
  h2->Draw();

  c8->cd(3);
  // TH1F* h3 = wwATGC1DFit("samples/processed_data_WW_1j_kg100_lg50_gg100_kz100_lz50_gz100_fastsim386_v1.root",
  TH1F* h3 = wwATGC1DFit("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root",
			 "p_sm",0.5,0);
  h3->SetTitle("#lambda_{Z}=0.5, #Delta g^{Z}_{1}=0");
  h3->GetXaxis()->SetTitle("|#lambda_{Z}|");
  h3->Draw();
    
  c8->cd(4);
  // TH1F* h4 = wwATGC1DFit("samples/processed_data_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz100_fastsim386_v1.root",
  TH1F* h4 = wwATGC1DFit("samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz1000000.root",
			 "p_sm",-0.5,0);
  h4->SetTitle("#lambda_{Z}=-0.5, #Delta g^{Z}_{1}=0");
  h4->GetXaxis()->SetTitle("|#lambda_{Z}|");
  h4->Draw();
  c8->Print("fit_wwATGC_mc_1D_abs.pdf");

}

void wwATGC1DLzKgFits()
{
  TCanvas* c8 = new TCanvas("c8","Fit results",800,400);
  c8->Divide(2,1);
  
  c8->cd(1);
  // TH1F* h1 = wwATGC1DFit("samples/processed_data_WW_1j_kg170_lg0_gg100_kz779311_lz0_gz100.root",
  TH1F* h1 = wwATGC1DFit("samples/processed_data_job_WW_1j_kg170_lg0_gg100_kz779311_lz0_gz100.root",
			 "sm_p",0,0.7);
  h1->SetTitle(Form("%s 0.7",y_par->GetTitle()));
  h1->GetXaxis()->SetTitle(y_par->GetTitle());
  h1->Draw();

  c8->cd(2);
  // TH1F* h2 = wwATGC1DFit("samples/processed_data_WW_1j_kg30_lg0_gg100_kz1220689_lz0_gz100.root",
  TH1F* h2 = wwATGC1DFit("samples/processed_data_job_WW_1j_kg30_lg0_gg100_kz1220689_lz0_gz100.root",
			 "sm_m",0,-0.7);
  h2->SetTitle(Form("%s -0.7",y_par->GetTitle()));
  h2->GetXaxis()->SetTitle(y_par->GetTitle());
  h2->Draw();
  c8->Print("fit_wwATGC_mc_1D_abs2.pdf");
}

void prepareForDataFits(){
  glb_data = MakeDataset("smurf/data.root","data",true);
  cpdf = new RooAddPdf("cpdf","combined pdf",RooArgList(*cSigPdf,*cBkgPdf));
}  

TH2* visualizeCorrelations(const TMatrixTSym<double>& matrix, double precision = 0.01 ){
  unsigned int nColumns = matrix.GetNcols();
  unsigned int nRows = matrix.GetNrows();
  TH2D* hist = new TH2D("cor","Correlation matrix",nColumns,0,nColumns,nRows,0,nRows);
  hist->SetDirectory(0);
  for (unsigned int i=0;i<nColumns;++i)
    for (unsigned int j=0;j<nRows;++j){
      hist->SetBinContent(i+1,j+1,round(matrix(j,i)/precision)*precision);
    }
  return hist;
}

void DrawCorrelationMatrix(const TMatrixTSym<double>& matrix ){
  TCanvas* cc1 = new TCanvas("cc1","cc1",500,500);
  Int_t palette[7];
  palette[3] = kWhite;
  for (unsigned int i=0;i<3;++i){
    palette[2-i] = kGray+i;
    palette[4+i] = kGray+i;
  }
  gStyle->SetPalette(7,palette);

  TH2* hist = visualizeCorrelations(matrix);
  hist->SetMaximum(1.0);
  hist->SetMinimum(-1.0);
  hist->GetXaxis()->SetBinLabel(1,"DY");
  hist->GetXaxis()->SetBinLabel(2,"Top");
  hist->GetXaxis()->SetBinLabel(3,"Wjets");
  hist->GetXaxis()->SetBinLabel(4,"WW");
  hist->GetXaxis()->SetBinLabel(5,"WZ");
  hist->GetXaxis()->SetBinLabel(6,"ZZ");
  hist->GetXaxis()->SetBinLabel(7,"#lambda_{Z}");
  hist->GetXaxis()->SetBinLabel(8,"#Delta g^{Z}_{1}");
  hist->GetYaxis()->SetBinLabel(1,"DY");
  hist->GetYaxis()->SetBinLabel(2,"Top");
  hist->GetYaxis()->SetBinLabel(3,"Wjets");
  hist->GetYaxis()->SetBinLabel(4,"WW");
  hist->GetYaxis()->SetBinLabel(5,"WZ");
  hist->GetYaxis()->SetBinLabel(6,"ZZ");
  hist->GetYaxis()->SetBinLabel(7,"#lambda_{Z}");
  hist->GetYaxis()->SetBinLabel(8,"#Delta g^{Z}_{1}");
  hist->Draw("coltext");
  cc1->Update();

//   TPaletteAxis* paletteAxis = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject( "palette" );
//   if ( paletteAxis ){
//     paletteAxis->SetLabelSize( 0.03 );
//     paletteAxis->SetX1NDC( paletteAxis->GetX1NDC() + 0.02 );
//   }    
//   hist->Draw("sametext");
}

void fitData( bool makeContourAllIn = true,
	      bool fitATGCx = true,
	      bool fitATGCy = true,
	      bool yieldFitNoATGC = false,
	      bool makeContourWWOnly = false
	     )
{
  prepareForDataFits();
  // TCanvas* c10 = new TCanvas("c10","c10",1000,1000);
  // c10->Divide(2,2);
  // c10->cd(1);
  // First plot data and pdf with the default parameters without any fit
  /*
  ((TTree*)ds_data_pt1->tree())->Draw(Form("pt1>>h(%u,%f,%f)",Nbins,minPt,maxPt));
  glb_data = dynamic_cast<RooDataSet*>(ds_data_pt1);
  RooNDKeysPdf* pdf_data = new RooNDKeysPdf("pdf_data","pdf_data",*var_pt1,*((RooDataSet*)ds_data_pt1),"am");
  TH1F* hpdf_data = (TH1F*)pdf_data->createHistogram("hpdf_data",*var_pt1);
  hpdf_data->SetTitle("DATA");
  hpdf_data->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_data->Scale(hpdf_data->Integral());
  c10->cd(2);
  hpdf_data->Draw();
  */
  
  // RooArgSet constrainedParams(*n_ww,*n_top,*n_wjets,*n_wz,*n_zz);
  RooArgSet constrainedParams(*n_ww,*n_top,*n_wjets,*n_wz,*n_zz,*n_dy);
  if (1) {
    n_top->setConstant(0);
    n_wjets->setConstant(0);
    n_wz->setConstant(0);
    n_zz->setConstant(0);
  }

  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  y_par->setConstant(0);
  x_par->setConstant(0);
  y_par->setVal(0);
  x_par->setVal(0);

  // signal only pdf
  if (makeContourWWOnly)
  {
    TCanvas* c11 = new TCanvas("c11","c11",500,500);
    c11->SetGrid();
    y_par->setVal(0);
    x_par->setVal(0);
    RooAbsReal* nll = atgcPdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    RooMinuit m(*nll);
    m.migrad();
    assert(m.save()->status()==0);
    m.hesse();
    assert(m.save()->status()==0);
    // RooPlot* p1 = m.contour(*x_par,*y_par,sqrt(2.3),sqrt(6.0));
    RooPlot* p1 = m.contour(*x_par,*y_par,sqrt(2.3),0);
    p1->SetTitle("68% and 95% C.L.");
    m.migrad();
    m.hesse();
    p1->Draw();
  }
  
  // signal + background pdf 
  if (makeContourAllIn)
  {
    TCanvas* c12 = new TCanvas("c12","c12",500,500);
    c12->SetGrid();
    y_par->setVal(0);
    x_par->setVal(0);
    RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    RooMinuit m(*nll);
    // RooMinimizer m(*nll);
    // m.setMinimizerType("Minuit2");
    m.migrad();
    assert(m.save()->status()==0);
    m.hesse();
    assert(m.save()->status()==0);
    RooPlot* p1 = m.contour(*x_par,*y_par,sqrt(2.3),sqrt(6.0));
    p1->SetTitle("68% and 95% C.L.");
    m.migrad();
    m.hesse();
    p1->Draw();
  }

  if (0){
    RooMsgService::instance().setStreamStatus(1,0);
    // fits
    x_par->setConstant(0);
    y_par->setConstant(0);
    x_par->setVal(0);
    y_par->setVal(0);
    cpdf->fitTo(*glb_data);
    double minNLL = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams))->getVal();
    TMarker* marker = new TMarker(x_par->getVal(),y_par->getVal(),20);

    TCanvas* c12_2 = new TCanvas("c12_2","c12_2",500,500);
    c12_2->SetGrid();
    const Int_t nBins = 200;
    TH2F* nll2d = new TH2F("nll2d","#sqrt{2 #Delta logL} for 35.5/pb of data", nBins,-1,1,nBins,-1,1);
    nll2d->GetXaxis()->SetTitle(x_par->GetTitle());
    nll2d->GetYaxis()->SetTitle(y_par->GetTitle());
    nll2d->SetDirectory(0);
    for ( Int_t xi = 1; xi<=nBins; ++xi ){
      cout << ".";
      cout.flush();
      for ( Int_t yi = 1; yi<=nBins; ++yi )
	{
	  x_par->setVal(-0.95+2.0*(xi-1)/nBins);
	  y_par->setVal(-0.95+2.0*(yi-1)/nBins);
	  RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams)/*,RooFit::NumCPU(2),RooFit::Verbose(0)*/);
	  nll->addServer(*var_dummy);
	  // double sig = sqrt(2*fabs(nll->getVal()-minNLL));
	  nll2d->SetBinContent(xi,yi,2*fabs(nll->getVal()-minNLL));
	  delete nll;
	}
    }
    nll2d->SetStats(kFALSE);
    nll2d->Draw("colz");
    // marker->SetMarkerColor(kBlue);
    marker->Draw();
    RooMsgService::instance().setStreamStatus(1,1);
  }

  if (yieldFitNoATGC){
    // yield fit
    x_par->setConstant(1);
    y_par->setConstant(1);
    // n_wjets->setConstant(1);
    // n_wjets->setVal(wjets_expected);
    x_par->setVal(0);
    y_par->setVal(0);
    // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());
    RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    // RooMinuit m(*nll);
    RooMinimizer m(*nll);
    m.setMinimizerType("Minuit2");
    m.setErrorLevel(0.5); //68% C.L. 
    m.migrad();
    m.hesse();
    m.minos();
  }

  if(fitATGCx)
  {
    // fits
    x_par->setConstant(0);
    y_par->setConstant(1);
    // n_wjets->setConstant(1);
    // n_wjets->setVal(wjets_expected);
    x_par->setVal(0);
    y_par->setVal(0);
    // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());
    RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    // RooMinuit m(*nll);
    RooMinimizer m(*nll);
    m.setMinimizerType("Minuit2");
    m.setErrorLevel(3.84*0.5); //95% C.L. 
    m.migrad();
    m.hesse();
    m.minos();
  }

  if (fitATGCy)
  {
    x_par->setConstant(1);
    y_par->setConstant(0);
    x_par->setVal(0);
    y_par->setVal(0);
    // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());
    RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    // RooMinuit m(*nll);
    RooMinimizer m(*nll);
    m.setMinimizerType("Minuit2");
    m.setErrorLevel(3.84*0.5); //95% C.L. 
    m.migrad();
    m.hesse();
    m.minos();
  }
}

void fitTop()
{
  TFile* f = TFile::Open("samples/processed_data_final.root");
  assert(f);
  RooAbsData* ds_ttbar = (RooAbsData*)f->Get("ttbar");
  ds_ttbar->SetName("ds_ttbar");
  RooAbsData* ds_ttbar_pt1 = ds_ttbar->reduce(*var_pt1,oldSampleSelection);

  TTree* iTree = const_cast<TTree*>(((RooDataSet*)ds_ttbar_pt1)->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  unsigned int nEvents(ww_expected);
  TTree* tree = iTree->CopyTree("","",nEvents,0);
  assert(tree);
  glb_data = new RooDataSet("data","data",tree, *var_pt1);

  RooAbsPdf* combinedPdf = new RooAddPdf("combinedPdf","combined pdf",RooArgList(*atgcPdf,*pdf_bkg));

  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  y_par->setConstant(0);
  x_par->setConstant(0);
  y_par->setVal(0);
  x_par->setVal(0);

  // signal + background pdf 
//   {
//     TCanvas* c13 = new TCanvas("c13","c13",500,500);
//     c13->SetGrid();
//     y_par->setVal(0);
//     x_par->setVal(0);
//     RooAbsReal* nll = combinedPdf->createNLL(*data,RooFit::Extended());
//     RooMinuit m(*nll);
//     m.migrad();
//     m.hesse();
//     RooPlot* p1 = m.contour(*x_par,*y_par);
//     m.migrad();
//     m.hesse();
//     p1->Draw();
//   }

  {
    // fits
    x_par->setConstant(0);
    y_par->setConstant(0);
    x_par->setVal(0);
    y_par->setVal(0);
    combinedPdf->fitTo(*glb_data);
    double minNLL = combinedPdf->createNLL(*glb_data)->getVal();
    TMarker* marker = new TMarker(x_par->getVal(),y_par->getVal(),20);

    TCanvas* c13_2 = new TCanvas("c13_2","c13_2",500,500);
    c13_2->SetGrid();
    const Int_t nBins = 200;
    TH2F* nll2d = new TH2F("nll2d","TTbar Monte Carlo as WW aTGC", nBins,-1,1,nBins,-1,1);
    nll2d->SetTitle("2#Delta#ln L");
    nll2d->GetXaxis()->SetTitle("|#lambda_{Z}|");
    nll2d->GetYaxis()->SetTitle("|#Delta#kappa_{Z}|");
    nll2d->SetDirectory(0);
    for ( Int_t xi = 1; xi<=nBins; ++xi ){
      cout << ".";
      cout.flush();
      for ( Int_t yi = 1; yi<=nBins; ++yi )
	{
	  x_par->setVal(-0.95+2.0*(xi-1)/nBins);
	  y_par->setVal(-0.95+2.0*(yi-1)/nBins);
	  RooAbsReal* nll = combinedPdf->createNLL(*glb_data,RooFit::Extended());
	  nll->addServer(*var_dummy);
	  // double sig = sqrt(2*fabs(nll->getVal()-minNLL));
	  nll2d->SetBinContent(xi,yi,2*fabs(nll->getVal()-minNLL));
	  delete nll;
	}
    }
    nll2d->Draw("colz");
    // marker->SetMarkerColor(kBlue);
    marker->Draw();
  }

  // fits
  x_par->setConstant(0);
  y_par->setConstant(1);
  x_par->setVal(0);
  y_par->setVal(0);
  combinedPdf->fitTo(*glb_data, RooFit::Minos());
  
  x_par->setConstant(1);
  y_par->setConstant(0);
  x_par->setVal(0);
  y_par->setVal(0);
  combinedPdf->fitTo(*glb_data, RooFit::Minos());

}

TH1* MakePlot(RooDataSet* ds, const char* name){
  TTree* iTree = const_cast<TTree*>(ds->tree());
  assert(iTree);
  iTree->Draw(Form("pt1>>%s(%u,%f,%f)",name,Nbins,minPt,maxPt),"weight","goff");
  TH1* hist = (TH1*)gDirectory->Get(name);
  hist->Sumw2();
  hist->Scale(1/hist->Integral());
  return hist;
}

void compareDistributions()
{
  var_pt1   = new RooRealVar("pt1","pt1",minPt,maxPt);
  var_dummy = new RooRealVar("var_dummy","var_dummy",0);
  var_pt1->setBins(Nbins);

  new TCanvas("cc","cc",500,500);
  // TH1* h1 = MakePlot(MakeOldDataset("samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz100_lz0_gz1000000.root"),"Sherpa-FastSim");
  TH1* h1 = MakePlot(MakeOldDataset("samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root"),"Sherpa-FastSim");
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2);
  h1->Draw("hist");
//   TH1* h2 = MakePlot(MakeDataset("smurf/qqww.root","qqww1"),"Madgraph-FullSim");
//   h2->SetLineColor(kBlack);
//   h2->SetLineWidth(2);
//   //  h2->Draw("same hist c");
//   h2->Draw("same hist");
  TH1* h3 = MakePlot(MakeDataset("WWqqbr_tota_mstw8nl_80__80__kg100_lg50_kZ100_lZ50_gZ100_smurf.root","qqww2"),"MCFM-GenOnly");
  h3->SetLineColor(kRed);
  h3->SetLineWidth(2);
  h3->Draw("same hist");
} 

void getDataYield(){
  MakeDataset("smurf/data.root","data",true);
}

void yieldDataFit()
{
  cpdf = new RooAddPdf("cpdf","combined pdf",RooArgList(*cSigPdf,*cBkgPdf));
  x_par->setConstant(0);
  y_par->setConstant(0);
  x_par->setVal(0);
  y_par->setVal(0);
  // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());                                                                                                                                
  RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(RooArgSet(*n_top,*n_wjets,*n_wz,*n_zz,*n_dy)));
  // RooMinuit m(*nll);                                                                                                                                                         
  RooMinimizer m(*nll);
  m.setMinimizerType("Minuit2");
  m.setErrorLevel(0.5); //68% C.L.                                                                                                                                              
  m.migrad();
  m.hesse();
  m.minos();
  
  DrawCorrelationMatrix(m.lastMinuitFit()->correlationMatrix());

  // Plot
  RooPlot* frame = var_pt1->frame(RooFit::Title("Leading lepton Pt distribution in Data")) ;
  glb_data->plotOn(frame) ;
  cpdf->plotOn(frame,RooFit::Precision(1e-5)) ;
  new TCanvas("cc2","cc2",600,600);
  frame->Draw();
}

TH1F* toy_wwATGC1DFit(const char* name, double lz, double dkz, unsigned int ntoys=1000)
{
  x_par->setRange(-1,1);
  y_par->setRange(-1,1);
  assert(lz==0||dkz==0);

  unsigned int nEvents(ww_expected);
  glb_data=0;
  // TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW with anomalous couplings",40,-1,1);
  TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW toy MC data with anomalous couplings",60,0,1.5);
  h->SetDirectory(0);

  for(unsigned int i=0; i<ntoys; ++i){
    x_par->setVal(lz);
    y_par->setVal(dkz);
    // TRandom::Poisson
    RooDataSet* dataset = atgcPdf->generate(*var_pt1,nEvents);
    if (glb_data) delete glb_data;
    glb_data = dataset;

    if ( lz==0 ){
      x_par->setConstant(1);
      y_par->setConstant(0);
    } else {
      x_par->setConstant(0);
      y_par->setConstant(1);
    }      
    atgcPdf->fitTo(*glb_data);
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;
    if ( lz==0 )
      h->Fill(fabs(y_par->getVal()));
    else
      h->Fill(fabs(x_par->getVal()));
  }
  h->SetFillColor(kMagenta);
  return h;
}

TH1F* toy_fullATGC1DFit(const char* name, double lz, double dkz, unsigned int ntoys=1000, bool fit1 = true) // when fit1 false fit the other variable
{
  x_par->setRange(-1,1);
  y_par->setRange(-1,1);
  cpdf = new RooAddPdf("cpdf","combined pdf",RooArgList(*cSigPdf,*cBkgPdf));

  unsigned int nEvents(ww_expected+top_expected+wjets_expected+wz_expected+zz_expected+dy_expected);
  glb_data=0;
  // TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW with anomalous couplings",40,-1,1);
  TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW toy MC data with anomalous couplings",100,-0.5,0.5);
  h->SetDirectory(0);
  for(unsigned int i=0; i<ntoys; ++i){
    x_par->setVal(lz);
    y_par->setVal(dkz);
    // TRandom::Poisson
    RooDataSet* dataset = cpdf->generate(*var_pt1,nEvents);
    if (glb_data) delete glb_data;
    glb_data = dataset;

    if ( fit1 ){
      x_par->setConstant(0);
      y_par->setConstant(1);
    } else {
      x_par->setConstant(1);
      y_par->setConstant(0);
    }  
    RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(RooArgSet(*n_top,*n_wjets,*n_wz,*n_zz,*n_dy)));
    // RooMinuit m(*nll);                                                                                                                                                         
    RooMinimizer m(*nll);
    m.setMinimizerType("Minuit2");
    m.setErrorLevel(0.5); //68% C.L.                                                                                                                                              
    m.migrad();
    
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;
    if ( fit1 )
      h->Fill(x_par->getVal());
    else
      h->Fill(y_par->getVal());
  }
  h->SetFillColor(kMagenta);
  return h;
}

void toy_wwATGC1DFits()
{
  TCanvas* c8 = new TCanvas("c8toy","Fit results",800,800);
  c8->Divide(2,2);
  
  c8->cd(1);
  TH1F* h1 = toy_wwATGC1DFit("sm_p",0,0.75);
  h1->SetTitle("#lambda_{Z}=0, #Delta g^{Z}_{1}=0.75");
  h1->GetXaxis()->SetTitle("|#Delta g^{Z}_{1}=0.75|");
  h1->Draw();

  c8->cd(2);
  TH1F* h2 = toy_wwATGC1DFit("sm_m",0,-0.75);
  h2->SetTitle("#lambda_{Z}=0, #Delta g^{Z}_{1}=-0.75");
  h2->GetXaxis()->SetTitle("|#Delta g^{Z}_{1}|");
  h2->Draw();

  c8->cd(3);
  TH1F* h3 = toy_wwATGC1DFit("p_sm",0.5,0);
  h3->SetTitle("#lambda_{Z}=0.5, #Delta g^{Z}_{1}=0");
  h3->GetXaxis()->SetTitle("|#lambda_{Z}|");
  h3->Draw();
    
  c8->cd(4);
  TH1F* h4 = toy_wwATGC1DFit("p_sm",-0.5,0);
  h4->SetTitle("#lambda_{Z}=-0.5, #Delta g^{Z}_{1}=0");
  h4->GetXaxis()->SetTitle("|#lambda_{Z}|");
  h4->Draw();
}

void toy_wwATGC1DLzKgFits()
{
  TCanvas* c8 = new TCanvas("c8toy2","Fit results",800,400);
  c8->Divide(2,1);
  
  c8->cd(1);
  TH1F* h1 = toy_wwATGC1DFit("sm_p",0,0.7);
  h1->SetTitle(Form("%s 0.7",y_par->GetTitle()));
  h1->GetXaxis()->SetTitle(y_par->GetTitle());
  h1->Draw();

  c8->cd(2);
  TH1F* h2 = toy_wwATGC1DFit("sm_m",0,-0.7);
  h2->SetTitle(Form("%s -0.7",y_par->GetTitle()));
  h2->GetXaxis()->SetTitle(y_par->GetTitle());
  h2->Draw();
}

void toy_fullATGC1DFits()
{
  TCanvas* c81 = new TCanvas("c8fulltoy_gz1","Fit results",800,800);
  c81->cd();
  TH1F* h1 = toy_fullATGC1DFit("toy_gz1",0,0,1000,false);
  h1->SetTitle("#lambda_{Z}=0, #Delta g^{Z}_{1}=0");
  h1->GetXaxis()->SetTitle("#Delta g^{Z}_{1}");
  h1->Draw();

  TCanvas* c82 = new TCanvas("c8fulltoy_lz","Fit results",800,800);
  c82->cd();
  TH1F* h3 = toy_fullATGC1DFit("toy_lz",0,0,1000,true);
  h3->SetTitle("#lambda_{Z}=0, #Delta g^{Z}_{1}=0");
  h3->GetXaxis()->SetTitle("#lambda_{Z}");
  h3->Draw();
}

void toy_fullATGC1DLzKgFits()
{
  TCanvas* c8 = new TCanvas("c8fulltoy_kg","Fit results",800,800);
  c8->cd(1);
  TH1F* h1 = toy_fullATGC1DFit("toy_kg",0,0,1000,false);
  h1->SetTitle(Form("%s 0.7",y_par->GetTitle()));
  h1->GetXaxis()->SetTitle(y_par->GetTitle());
  h1->Draw();
}

// just one sample (fitted like signal)
void simpleFitForATGC(const char* file, bool oldformat=false)
{
  // n_ww->setConstant(0);
  x_par->setRange(-1.5,1.5);
  y_par->setRange(-1.5,1.5);

  RooDataSet* dataset(0);
  if (!oldformat) 
    dataset = MakeDataset(file,"ds");
  else 
    dataset = MakeOldDataset(file);
  glb_data = dataset;
  cout << "Num entries: " << glb_data->numEntries() << endl;

  x_par->setVal(0);
  y_par->setVal(0);
  // first fit
  x_par->setConstant(0);
  y_par->setConstant(1);
  cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
  cout << "First fit (Y is constant)" << endl;
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;
  atgcPdf->fitTo(*glb_data);
  cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
  cout << "First fit is done" << endl;
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;

  x_par->setVal(0);
  y_par->setVal(0);
  // second fit
  x_par->setConstant(1);
  y_par->setConstant(0);
  cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
  cout << "Second fit (X is constant)" << endl;
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;
  atgcPdf->fitTo(*glb_data);
  cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
  cout << "Second fit is done" << endl;
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;
}
