#include <vector>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"

//
// ** apply phil's recoil corrections **
// 
// usage: 
//    double met=rawMetValue;
//    double metphi=rawMetPhiValue;
//    RecoilCorrector corrector;
//    corrector->Correct(met,metphi,GenZPt,GenZPhi,leptonPt,leptonPhi);
//    printf("corrected met: %10.2f%10.2f\n",met,metphi);
//
// where leptonPt, leptonPhi are dilepton kinematics for z->ll and single lepton kinematics for w->lnu
//

using namespace std;

class RecoilCorrector
{
  
public:
  RecoilCorrector(string iNameZDat, int iSeed=0xDEADBEEF);
  RecoilCorrector(string iNameZDat1, string iPrefix, int iSeed=0xDEADBEEF);
  void Correct(double &met, double &metphi, double iGenPt, double iGenPhi, double iLepPt, double iLepPhi);
  void Correct(double &pfmet, double &pfmetphi, double &trkmet, double &trkmetphi, 
	       double iGenPt, double iGenPhi, double iLepPt, double iLepPhi,double iFluc=0);
  
protected:
  void readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,
		  std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit,
		  std::string iFName,std::string iPrefix); 
  void readCorr(std::string iName,int iType=2);

  void metDistribution(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
		       double iLepPt,double iLepPhi,TRandom1 *iRand,
		       TF1 *iU1RZFit, 
		       TF1 *iU1MSZFit, 
		       TF1 *iU1S1ZFit,
		       TF1 *iU1S2ZFit,
		       TF1 *iU2MSZFit,		   
		       TF1 *iU2S1ZFit, 
		       TF1 *iU2S2ZFit, 
		       TF1 &iU1U2Corr, 
		       int iFluc=0);

  void metDistribution(double &iPFMet,double &iPFMPhi,double &iTKMet,double &iTKMPhi,
		       double iGenPt,double iGenPhi,
		       double iLepPt,double iLepPhi,TRandom1 *iRand,
		       TF1 *iU1RZPFFit,  TF1 *iU1RZTKFit, 
		       TF1 *iU1MSZPFFit, TF1 *iU1MSZTKFit, 
		       TF1 *iU1S1ZPFFit, TF1 *iU1S1ZTKFit,
		       TF1 *iU1S2ZPFFit, TF1 *iU1S2ZTKFit,
		       TF1 *iU2MSZPFFit, TF1 *iU2MSZTKFit,		   
		       TF1 *iU2S1ZPFFit, TF1 *iU2S1ZTKFit, 
		       TF1 *iU2S2ZPFFit, TF1 *iU2S2ZTKFit, 
		       TF1 &iPFU1U2Corr, TF1 &iTKU1U2Corr,
		       TF1 &iPFTKU1Corr, TF1 &iPFTKU2Corr,
		       TF1 &iPFTKU1MCorr,TF1 &iPFTKU2MCorr,
		       int iFluc=0);

  double calculate(int iMet,double iEPt,double iEPhi,double iWPhi,double iU1,double iU2);
  double getError(double iVal,TF1 *iZDatFit);
  double getError2(double iVal,TF1 *iFit);
  double getCorError2(double iVal,TF1 *iFit);
  double mag(double iV0,double iV1,double iV2,double iV3);
  double correlatedSeed(double iVal, double iCorr1,double iCorr2,double iCorr3,double iSeed0,double iSeed1,double iSeed2,double iSeed3);

  TRandom1 *fRandom; 
  vector<TF1*> fF1U1Fit; vector<TF1*> fF1U1RMSSMFit; vector<TF1*> fF1U1RMS1Fit; vector<TF1*> fF1U1RMS2Fit; 
  vector<TF1*> fF1U2Fit; vector<TF1*> fF1U2RMSSMFit; vector<TF1*> fF1U2RMS1Fit; vector<TF1*> fF1U2RMS2Fit; 
  vector<TF1*> fF2U1Fit; vector<TF1*> fF2U1RMSSMFit; vector<TF1*> fF2U1RMS1Fit; vector<TF1*> fF2U1RMS2Fit; 
  vector<TF1*> fF2U2Fit; vector<TF1*> fF2U2RMSSMFit; vector<TF1*> fF2U2RMS1Fit; vector<TF1*> fF2U2RMS2Fit; 
  TF1 fF1U1U2Corr;     TF1 fF2U1U2Corr;
  TF1 fF1F2U1Corr;     TF1 fF1F2U2Corr;
  TF1 fF1F2U1U2Corr;   TF1 fF1F2U2U1Corr;
};

//-----------------------------------------------------------------------------------------------------------------------------------------
  RecoilCorrector::RecoilCorrector(string iNameZDat,std::string iPrefix, int iSeed) {

  fRandom = new TRandom1(iSeed);

  // get fits for Z data
  readRecoil(fF1U1Fit,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,iNameZDat,iPrefix);
  if(iPrefix == "PF") readCorr  (iNameZDat,0);
  if(iPrefix == "TK") readCorr  (iNameZDat,1);  
}

RecoilCorrector::RecoilCorrector(string iNameZ, int iSeed) {

  fRandom = new TRandom1(iSeed);
  // get fits for Z data
  readRecoil(fF1U1Fit,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,iNameZ,"PF");
  readRecoil(fF2U1Fit,fF2U1RMSSMFit,fF2U1RMS1Fit,fF2U1RMS2Fit,fF2U2Fit,fF2U2RMSSMFit,fF2U2RMS1Fit,fF2U2RMS2Fit,iNameZ,"TK");
  readCorr  (iNameZ);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::Correct(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi) {
    metDistribution(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi,fRandom,
		    fF1U1Fit     [0],
		    fF1U1RMSSMFit[0],
		    fF1U1RMS1Fit [0],
		    fF1U1RMS2Fit [0],
		    fF1U2RMSSMFit[0],
		    fF1U2RMS1Fit [0],
		    fF1U2RMS2Fit [0],
		    fF1U1U2Corr);
}

void RecoilCorrector::Correct(double &pfmet, double &pfmetphi, double &trkmet, double &trkmetphi, 
                              double lGenPt, double lGenPhi, double lepPt, double lepPhi,double iFluc) {
  metDistribution(pfmet,pfmetphi,trkmet,trkmetphi,lGenPt,lGenPhi,lepPt,lepPhi,fRandom,
		  fF1U1Fit     [0],fF2U1Fit     [0],
		  fF1U1RMSSMFit[0],fF2U1RMSSMFit[0],
		  fF1U1RMS1Fit [0],fF2U1RMS1Fit [0],
		  fF1U1RMS2Fit [0],fF2U1RMS2Fit [0],
		  fF1U2RMSSMFit[0],fF2U2RMSSMFit[0],
		  fF1U2RMS1Fit [0],fF2U2RMS1Fit [0],
		  fF1U2RMS2Fit [0],fF2U2RMS2Fit [0],
		  fF1U1U2Corr     ,fF2U1U2Corr     ,
		  fF1F2U1Corr     ,fF1F2U2Corr     ,
		  fF1F2U1U2Corr   ,fF1F2U2U1Corr   ,
		  iFluc
		  );
}

//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,
		                 std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit,
		                 std::string iFName,std::string iPrefix) {
//  if(!getenv("CMSSW_BASE")) {
//    printf("error! RecoilCorrector called without input files. Define CMSSW_BASE or add by hand.\n");
//    assert(0);
//  }
  TFile *lFile  = new TFile(iFName.c_str());
  iU1Fit.push_back    ( (TF1*) lFile->FindObjectAny((iPrefix+"u1Mean_0").c_str()));
  iU1MRMSFit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1MeanRMS_0").c_str()));
  iU1RMS1Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1RMS1_0").c_str()));
  iU1RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1RMS2_0").c_str()));
  iU2Fit    .push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2Mean_0").c_str()));
  iU2MRMSFit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2MeanRMS_0").c_str()));
  iU2RMS1Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS1_0").c_str()));
  iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS2_0").c_str()));
  iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS2_0").c_str()));
  lFile->Close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::readCorr(std::string iName,int iType) {
  TFile *lFile = new TFile(iName.c_str());
  std::stringstream pSS1,pSS2,pSS3,pSS4,pSS5,pSS6;
  if(iType != 1) {pSS1  << "u1u2pfCorr" ;   fF1U1U2Corr  = *((TF1*) lFile->FindObjectAny(pSS1.str().c_str())); }
  if(iType != 0) {pSS2  << "u1u2tkCorr";    fF2U1U2Corr  = *((TF1*) lFile->FindObjectAny(pSS2.str().c_str())); }
  if(iType <  2) {lFile->Close(); return;}
  pSS3  << "pftku1Corr";    fF1F2U1Corr   = *((TF1*) lFile->FindObjectAny(pSS3.str().c_str()));
  pSS4  << "pftku2Corr";    fF1F2U2Corr   = *((TF1*) lFile->FindObjectAny(pSS4.str().c_str()));
  pSS5  << "pftkum1Corr";   fF1F2U1U2Corr = *((TF1*) lFile->FindObjectAny(pSS5.str().c_str()));
  pSS6  << "pftkum2Corr";   fF1F2U2U1Corr = *((TF1*) lFile->FindObjectAny(pSS6.str().c_str()));
  lFile->Close();
}


//-----------------------------------------------------------------------------------------------------------------------------------------
// Met Prediction from Recoil Model
//-----------------------------------------------------------------------------------------------------------------------------------------
// Total Recoil Model is defined as : Recoil(u1) = pFrac*Gaussian1(u1;pSigma1,Mean) + (1-pFrac)*Gaussian2(u1;pSigma2;Mean)
// 4 Parameters: Mean, MeanRMS, Sigma1 (width of Gaus1), Sigma2 (width of Gaus2)
// pFrac can be solved in terms of MeanRMS, Sigma1, and Sigma2: pFrac = (MeanRMS - Sigma2) / (Sigma1 - Sigma2)
//
// For U2, we explicitly assume that the Mean is 0
// The fit was performed without this constraint but it has been verified (at least for PFMet) that the mean
// returned by the fit is consistent with a straight line at 0. 
// 
void RecoilCorrector::metDistribution(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
		                      double iLepPt,double iLepPhi,TRandom1 *iRand,
		                      TF1 *iU1RZDatFit,
		                      TF1 *iU1MSZDatFit, 
		                      TF1 *iU1S1ZDatFit,
		                      TF1 *iU1S2ZDatFit,
		                      TF1 *iU2MSZDatFit, 		   
		                      TF1 *iU2S1ZDatFit, 
		                      TF1 *iU2S2ZDatFit, 		                      
		                      TF1 &iU1U2Corr, 
		                      int iFluc) {
  double lRescale = sqrt((TMath::Pi())/2.);		     
  double pU1Mean       = iU1RZDatFit->Eval(iGenPt);
  double pU2Mean       = 0; 
  double pMeanRMS1 = iU1MSZDatFit->Eval(iGenPt)*lRescale;
  double pMeanRMS2 = iU2MSZDatFit->Eval(iGenPt)*lRescale;
  double pSigma1_1 = iU1S1ZDatFit->Eval(iGenPt)*lRescale*iU1MSZDatFit->Eval(iGenPt);
  double pSigma1_2 = iU1S2ZDatFit->Eval(iGenPt)*lRescale*iU1MSZDatFit->Eval(iGenPt);
  double pSigma2_1 = iU2S1ZDatFit->Eval(iGenPt)*lRescale*iU2MSZDatFit->Eval(iGenPt);
  double pSigma2_2 = iU2S2ZDatFit->Eval(iGenPt)*lRescale*iU2MSZDatFit->Eval(iGenPt);
  
  //Uncertainty propagation
  if(iFluc != 0) { 
    double lEUR1    = getError(iGenPt,iU1RZDatFit);
    double lEUS1_1  = getError(0.,iU1S1ZDatFit);
    double lEUS1_2  = getError(0.,iU1S2ZDatFit);
    double lEU1Frac = getError(iGenPt,iU1MSZDatFit);
    double lEUS2_1  = getError(0.,iU2S1ZDatFit);
    double lEUS2_2  = getError(0.,iU2S2ZDatFit);
    double lEU2Frac = getError(iGenPt,iU2MSZDatFit);
    
    //Modify all the different parameters the choice of signs makes it maximal
    pU1Mean   = pU1Mean   - iFluc*lEUR1;                //Recoil
    pMeanRMS1 = pMeanRMS1 + iFluc*(lEU1Frac);           //Mean RMS 
    pSigma1_1 = pSigma1_1 - iFluc*lEUS1_1*pMeanRMS1;    //Sigma 1 smalles sigma
    pSigma1_2 = pSigma1_2 + iFluc*lEUS1_2*pMeanRMS1;    //Sigma 2 (Maximal when oppsite sigma 1)
    pMeanRMS2 = pMeanRMS2 + iFluc*(lEU2Frac);           //Mean RMS for U2
    pSigma2_1 = pSigma2_1 - iFluc*lEUS2_1*pMeanRMS2;    //Sigma 1 U2
    pSigma2_2 = pSigma2_2 + iFluc*(lEUS2_2)*pMeanRMS2;
  }

  //Solve for pFrac in terms of MeanRMS, Sigma1, and Sigma2
  double pFrac1 = (pMeanRMS1 - pSigma1_2)/(pSigma1_1 - pSigma1_2);
  double pFrac2 = (pMeanRMS2 - pSigma2_2)/(pSigma2_1 - pSigma2_2);

  //Now sample for the MET distribution
  double pVal0 = iRand->Uniform(0,1);
  double pVal1 = iRand->Uniform(0,1);
  double pCorr1     = iRand->Gaus(0,1);     double pCorr2     = iRand->Gaus(0,1);  
  double pCorrT1    = iRand->Gaus(0,1);     double pCorrT2    = iRand->Gaus(0,1);  

  //Defines the slope of the correlation in the space of (sigma1, sigma2)
  //U1U2Corr gives the correlation of sigma1^2 vs sigma2^2, 
  //a result of a linear fit in the space of (sigma1^2,sigma2^2)
  //We find that Sqrt( Corr(sigma1^2,sigma2^2) ) gives a better model for Corr(sigma1, sigma2)
  //in Monte Carlo simulation.
  double CorrU1U2  = sqrt(TMath::Max(iU1U2Corr.Eval(iGenPt),0.));

  //***************
  //Is this right?
  // (pCorr1,pCorr2) is the random sampling in the uncorrelated axes defined by (1,CorrU1U2) and (CorrU1U2,1) 
  // [2-vectors expressed in the (sigma1,sigma2) coordinate system]
  // Then the projection done below, is a coordinate transformation that brings (pCorr1,pCorr2) back to the rectangular
  // (sigma1,sigma2) coordinate system.
  //
  // So explicitly :
  // (pVal1,pVal2) = CoordinateTransformation[ (pCorr1,pCorr2) ] 
  // where CoordinateTransformation is a linear transformation taking 2vectors expressed in the 
  // coordinate system defined by (1,CorrU1U2) and (CorrU1U2,1) to the rectangular coordinate system 
  // defined by (0,1) and (1,0).
  // 
  //****************

  //Project the random vector (pCorr1, pCorr2) onto the axes defined by the correlation slope CorrU1U2.
  //For the projection onto the axis orthogonal to CorrU1U2, we use the trick of taking the 
  //transposed vector (pCorr2,pCorr1), since it is equivalent to transposing the axis direction CorrU1U2.
  double pVal1_1 = correlatedSeed(pSigma1_1,CorrU1U2,0.,0.,pCorr1,pCorr2,0.,0.);
  double pVal2_1 = correlatedSeed(pSigma2_1,CorrU1U2,0.,0.,pCorr2,pCorr1,0.,0.);
  double pVal1_2 = correlatedSeed(pSigma1_2,CorrU1U2,0.,0.,pCorrT1,pCorrT2,0.,0.);
  double pVal2_2 = correlatedSeed(pSigma2_2,CorrU1U2,0.,0.,pCorrT2,pCorrT1,0.,0.);

  //Compute actual recoil components:
  //Recoil_{U1/U2} = pFrac * Gaussian[pMean,pVal_1] + (1-pFrac)*Gaussian[pMean,pVal_2]
  double pU1   = (pVal0 < pFrac1)*(pVal1_1+pU1Mean)+(pVal0 > pFrac1)*(pVal1_2+pU1Mean);
  double pU2   = (pVal1 < pFrac2)*(pVal2_1+pU2Mean)+(pVal1 > pFrac2)*(pVal2_2+pU2Mean);
  //pU1   = (lVal0 < pFrac1)*iRand->Gaus(pU1,pSigma1_1)+(lVal0 > pFrac1)*iRand->Gaus(pU1,pSigma1_2);
  //pU2   = (lVal1 < pFrac2)*iRand->Gaus(pU2,pSigma2_1)+(lVal1 > pFrac2)*iRand->Gaus(pU2,pSigma2_2);

  //Compute Met given recoil components (U1, U2), and Dilepton direction
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  return;
}

void RecoilCorrector::metDistribution(double &iPFMet,double &iPFMPhi,double &iTKMet,double &iTKMPhi,
				      double iGenPt,double iGenPhi,
		                      double iLepPt,double iLepPhi,TRandom1 *iRand,
		                      TF1 *iU1RPFFit, TF1 *iU1RTKFit,
		                      TF1 *iU1MSPFFit,TF1 *iU1MSTKFit, 
		                      TF1 *iU1S1PFFit,TF1 *iU1S1TKFit,
		                      TF1 *iU1S2PFFit,TF1 *iU1S2TKFit,
		                      TF1 *iU2MSPFFit,TF1 *iU2MSTKFit, 		      
		                      TF1 *iU2S1PFFit,TF1 *iU2S1TKFit, 
		                      TF1 *iU2S2PFFit,TF1 *iU2S2TKFit, 
		                      TF1 &iPFU1U2Corr, TF1 &iTKU1U2Corr,
		                      TF1 &iPFTKU1Corr, TF1 &iPFTKU2Corr,
		                      TF1 &iPFTKU1MCorr,TF1 &iPFTKU2MCorr,
		                      int iFluc) {
  //Important constants re-scaling of sigma on left and mean wpt of W resbos on right
  double lRescale  = sqrt((TMath::Pi())/2.); //double lPtMean = 16.3; //==> tuned for W bosons
  ///
  double pPFU1       = iU1RPFFit->Eval(iGenPt);
  double pPFU2       = 0;
  double pPFSigma1_1 = iU1S1PFFit->Eval(iGenPt)*iU1MSPFFit->Eval(iGenPt)*lRescale;
  double pPFSigma1_2 = iU1S2PFFit->Eval(iGenPt)*iU1MSPFFit->Eval(iGenPt)*lRescale;
  double pPFFrac1    = iU1MSPFFit->Eval(iGenPt)                         *lRescale;
  double pPFSigma2_1 = iU2S1PFFit->Eval(iGenPt)*iU2MSPFFit->Eval(iGenPt)*lRescale;
  double pPFSigma2_2 = iU2S2PFFit->Eval(iGenPt)*iU2MSPFFit->Eval(iGenPt)*lRescale;
  double pPFFrac2    = iU2MSPFFit->Eval(iGenPt)                         *lRescale;
  if(pPFSigma1_1 > pPFSigma1_2) {double pT = pPFSigma1_2; pPFSigma1_2 = pPFSigma1_1; pPFSigma1_1 = pT;}
  if(pPFSigma2_1 > pPFSigma2_2) {double pT = pPFSigma2_2; pPFSigma2_2 = pPFSigma2_1; pPFSigma2_2 = pT;}
  
  double pTKU1       = iU1RTKFit->Eval(iGenPt);
  double pTKU2       = 0;
  double pTKSigma1_1 = iU1S1TKFit->Eval(iGenPt)*iU1MSTKFit->Eval(iGenPt)*lRescale;
  double pTKSigma1_2 = iU1S2TKFit->Eval(iGenPt)*iU1MSTKFit->Eval(iGenPt)*lRescale;
  double pTKFrac1    = iU1MSTKFit->Eval(iGenPt)                         *lRescale;
  double pTKSigma2_1 = iU2S1TKFit->Eval(iGenPt)*iU2MSTKFit->Eval(iGenPt)*lRescale;
  double pTKSigma2_2 = iU2S2TKFit->Eval(iGenPt)*iU2MSTKFit->Eval(iGenPt)*lRescale;
  double pTKFrac2    = iU2MSTKFit->Eval(iGenPt)                         *lRescale;
  if(pTKSigma1_1 > pTKSigma1_2) {double pT = pTKSigma1_2; pTKSigma1_2 = pTKSigma1_1; pTKSigma1_1 = pT;}
  if(pTKSigma2_1 > pTKSigma2_2) {double pT = pTKSigma2_2; pTKSigma2_2 = pTKSigma2_1; pTKSigma2_2 = pT;}

  //Uncertainty propagation
  if(iFluc != 0) { 
    double lEUR1    = getError(iGenPt ,iU1RPFFit);
    double lEUS1_1  = getError(iGenPt ,iU1S1PFFit);
    double lEUS1_2  = getError(iGenPt ,iU1S2PFFit);
    double lEU1Frac = getError(iGenPt ,iU1MSPFFit)*lRescale;
    double lEUS2_1  = getError(iGenPt ,iU2S1PFFit);
    double lEUS2_2  = getError(iGenPt ,iU2S2PFFit);
    double lEU2Frac = getError(iGenPt ,iU2MSPFFit)*lRescale;
    //cout << "===> " << pPFSigma1_1 << " -- " << iU1S2PFFit->GetParError(0) << " -- " << lEUS1_1 << endl;
    //Modify all the different parameters the choice of signs makes it maximal
    if(iU1S1PFFit->Eval(iGenPt) > 1) {double pPF = lEUS1_1; lEUS1_1 = lEUS1_2; lEUS1_1 = pPF;}
    if(iU2S1PFFit->Eval(iGenPt) > 1) {double pPF = lEUS2_1; lEUS2_1 = lEUS2_2; lEUS2_1 = pPF;}

    pPFU1       = pPFU1       ;// iFluc*lEUR1;              //Recoil
    pPFFrac1    = pPFFrac1    ;//+ iFluc*(lEU1Frac);        //Mean RMS 
    pPFSigma1_1 = pPFSigma1_1 + iFluc*lEUS1_1*pPFFrac1;    //Sigma 1 smalles sigma
    pPFSigma1_2 = pPFSigma1_2 - iFluc*lEUS1_2*pPFFrac1;    //Sigma 2 (Maximal when oppsite sigma 1)
    pPFFrac2    = pPFFrac2    ;//+ iFluc*(lEU2Frac);        //Mean RMS for U2
    pPFSigma2_1 = pPFSigma2_1 ;//+ iFluc*lEUS2_1*pPFFrac2  ;    //Sigma 1 U2
    pPFSigma2_2 = pPFSigma2_2 ;//- iFluc*(lEUS2_2)*pPFFrac2;
    
    lEUR1    = getError(iGenPt,iU1RTKFit);
    lEUS1_1  = getError(iGenPt,iU1S1TKFit);
    lEUS1_2  = getError(iGenPt,iU1S2TKFit);
    lEU1Frac = getError(iGenPt,iU1MSTKFit)*lRescale;
    lEUS2_1  = getError(iGenPt,iU2S1TKFit);
    lEUS2_2  = getError(iGenPt,iU2S2TKFit);
    lEU2Frac = getError(iGenPt,iU2MSTKFit)*lRescale;
    if(iU1S1TKFit->Eval(iGenPt) > 1) {double pPF = lEUS1_1; lEUS1_1 = lEUS1_2; lEUS1_1 = pPF;}
    if(iU2S1TKFit->Eval(iGenPt) > 1) {double pPF = lEUS2_1; lEUS2_1 = lEUS2_2; lEUS2_1 = pPF;}
    //Modify all the different parameters the choice of signs makes it maximal
    pTKU1       = pTKU1        ;//+ iFluc*lEUR1;              //Recoil
    pTKFrac1    = pTKFrac1     ;//+ iFluc*(lEU1Frac);        //Mean RMS 
    pTKSigma1_1 = pTKSigma1_1  +  iFluc*lEUS1_1 *pTKFrac1;    //Sigma 1 smalles sigma
    pTKSigma1_2 = pTKSigma1_2  -  iFluc*lEUS1_2 *pTKFrac1;    //Sigma 2 (Maximal when oppsite sigma 1)
    pTKFrac2    = pTKFrac2     ;//+ iFluc*(lEU2Frac);        //Mean RMS for U2
    pTKSigma2_1 = pTKSigma2_1  ;//+ iFluc*lEUS2_1  *pTKFrac2;    //Sigma 1 U2
    pTKSigma2_2 = pTKSigma2_2  ;//+ iFluc*(lEUS2_2)*pTKFrac2;
  }
  //Caculat the proper fraction
  pPFFrac1 = (pPFFrac1-pPFSigma1_2)/(pPFSigma1_1-pPFSigma1_2);
  pPFFrac2 = (pPFFrac2-pPFSigma2_2)/(pPFSigma2_1-pPFSigma2_2);

  pTKFrac1 = (pTKFrac1-pTKSigma1_2)/(pTKSigma1_1-pTKSigma1_2);
  pTKFrac2 = (pTKFrac2-pTKSigma2_2)/(pTKSigma2_1-pTKSigma2_2);
  //if(iGenPt > 60) cout << "===> " << pFrac1 << " -- " <<pFrac2 << endl;
  double pPFVal0 = iRand->Uniform(0,1);
  double pPFVal1 = iRand->Uniform(0,1);

  double pTKVal0 = iRand->Uniform(0,1);
  double pTKVal1 = iRand->Uniform(0,1);

  double pPFCorr1     = iRand->Gaus(0,1);     double pPFCorr2     = iRand->Gaus(0,1);  
  double pPFCorrT1    = iRand->Gaus(0,1);     double pPFCorrT2    = iRand->Gaus(0,1);  

  double pTKCorr1     = iRand->Gaus(0,1);     double pTKCorr2     = iRand->Gaus(0,1);  
  double pTKCorrT1    = iRand->Gaus(0,1);     double pTKCorrT2    = iRand->Gaus(0,1);  

  double lPFU1U2  = sqrt(TMath::Max(iPFU1U2Corr.Eval(iGenPt),0.));//iPFU1U2Corr->Eval(iGenPt) ,0.);
  double lTKU1U2  = sqrt(TMath::Max(iTKU1U2Corr.Eval(iGenPt),0.));//iTKU1U2Corr->Eval(iGenPt) ,0.);
  double lPFTKU1  = iPFTKU1Corr.Eval(iGenPt); //TMath::Max(iPFTKU1Corr->Eval(iGenPt) ,0.);
  double lPFTKU2  = iPFTKU2Corr.Eval(iGenPt);//TMath::Max(iPFTKU2Corr->Eval(iGenPt) ,0.);
  double lPFTKU1M = 0;//iPFTKU1MCorr.Eval(iGenPt);
  double lPFTKU2M = 0;//iPFTKU2MCorr.Eval(iGenPt);
  //if(iGenPt < 20.) lPFU1U2/=2.;

  double pPFVal1_1 = correlatedSeed(pPFSigma1_1,lPFU1U2,lPFTKU1,lPFTKU1M,pPFCorr1,pPFCorr2,pTKCorr1,pTKCorr2);
  double pPFVal2_1 = correlatedSeed(pPFSigma2_1,lPFU1U2,lPFTKU2,lPFTKU2M,pPFCorr2,pPFCorr1,pTKCorr2,pTKCorr1);
  double pTKVal1_1 = correlatedSeed(pTKSigma1_1,lTKU1U2,lPFTKU1,lPFTKU2M,pTKCorr1,pTKCorr2,pPFCorr1,pPFCorr2);
  double pTKVal2_1 = correlatedSeed(pTKSigma2_1,lTKU1U2,lPFTKU2,lPFTKU1M,pTKCorr2,pTKCorr1,pPFCorr2,pPFCorr1);

  double pPFVal1_2 = correlatedSeed(pPFSigma1_2,lPFU1U2,lPFTKU1,lPFTKU1M,pPFCorrT1,pPFCorrT2,pTKCorrT1,pTKCorrT2);
  double pPFVal2_2 = correlatedSeed(pPFSigma2_2,lPFU1U2,lPFTKU2,lPFTKU2M,pPFCorrT2,pPFCorrT1,pTKCorrT2,pTKCorrT1);
  double pTKVal1_2 = correlatedSeed(pTKSigma1_2,lTKU1U2,lPFTKU1,lPFTKU2M,pTKCorrT1,pTKCorrT2,pPFCorrT1,pPFCorrT2);
  double pTKVal2_2 = correlatedSeed(pTKSigma2_2,lTKU1U2,lPFTKU2,lPFTKU1M,pTKCorrT2,pTKCorrT1,pPFCorrT2,pPFCorrT1);

  pPFU1   = (pPFVal0 < pPFFrac1)*(pPFVal1_1+pPFU1)+(pPFVal0 > pPFFrac1)*(pPFVal1_2+pPFU1);
  pPFU2   = (pPFVal1 < pPFFrac2)*(pPFVal2_1+pPFU2)+(pPFVal1 > pPFFrac2)*(pPFVal2_2+pPFU2);
  pTKU1   = (pTKVal0 < pTKFrac1)*(pTKVal1_1+pTKU1)+(pTKVal0 > pTKFrac1)*(pTKVal1_2+pTKU1);
  pTKU2   = (pTKVal1 < pTKFrac2)*(pTKVal2_1+pTKU2)+(pTKVal1 > pTKFrac2)*(pTKVal2_2+pTKU2);

  iPFMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pPFU1,pPFU2);
  iPFMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pPFU1,pPFU2);
  iTKMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pTKU1,pTKU2);
  iTKMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pTKU1,pTKU2);
  return;
  //Not used right now
  iPFTKU2Corr.Eval(0);
  iPFTKU1MCorr.Eval(0);
  iPFTKU2MCorr.Eval(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::calculate(int iMet,double iEPt,double iEPhi,double iWPhi,double iU1,double iU2) { 
  double lMX = -iEPt*cos(iEPhi) - iU1*cos(iWPhi) + iU2*sin(iWPhi);
  double lMY = -iEPt*sin(iEPhi) - iU1*sin(iWPhi) - iU2*cos(iWPhi);
  if(iMet == 0) return sqrt(lMX*lMX + lMY*lMY);
  if(iMet == 1) {if(lMX > 0) {return atan(lMY/lMX);} return (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX); } 
  if(iMet == 2) return lMX;
  if(iMet == 3) return lMY;
  return lMY;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::getCorError2(double iVal,TF1 *iFit) { 
  double lE = sqrt(iFit->GetParError(0))  + iVal*sqrt(iFit->GetParError(2));
  if(fabs(iFit->GetParError(4)) > 0) lE += iVal*iVal*sqrt(iFit->GetParError(4));
  return lE*lE;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::getError2(double iVal,TF1 *iFit) { 
  double lE2 = iFit->GetParError(0) + iVal*iFit->GetParError(1) + iVal*iVal*iFit->GetParError(2);
  if(fabs(iFit->GetParError(3)) > 0) lE2 += iVal*iVal*iVal*     iFit->GetParError(3);
  if(fabs(iFit->GetParError(4)) > 0) lE2 += iVal*iVal*iVal*iVal*iFit->GetParError(4);
  if(fabs(iFit->GetParError(5)) > 0 && iFit->GetParameter(3) == 0) lE2 += iVal*iVal*               iFit->GetParError(5);
  if(fabs(iFit->GetParError(5)) > 0 && iFit->GetParameter(3) != 0) lE2 += iVal*iVal*iVal*iVal*iVal*iFit->GetParError(5);
  if(fabs(iFit->GetParError(6)) > 0) lE2 += iVal*iVal*iVal*iVal*iVal*iVal*iFit->GetParError(6);
  return lE2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::getError(double iVal,TF1 *iZDatFit) {
  return sqrt(getError2(iVal,iZDatFit));
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::mag(double iV0,double iV1,double iV2,double iV3) { 
  return sqrt(iV0*iV0 + iV1*iV1 + iV2*iV2 + iV3*iV3);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::correlatedSeed(double iVal, double iCorr1,double iCorr2,double iCorr3,double iSeed0,double iSeed1,double iSeed2,double iSeed3) { 
  double lMag = mag(1.,iCorr1,iCorr2,iCorr3); 
  double lVal = ((1./lMag)*iSeed0 + (iCorr1/lMag)*iSeed1 + (iCorr2/lMag)*iSeed2 + (iCorr3/lMag)*iSeed3)*iVal;
  return lVal;
}
