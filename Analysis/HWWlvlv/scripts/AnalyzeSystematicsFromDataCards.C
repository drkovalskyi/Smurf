#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TLorentzVector.h>
#endif

//====================================================================================================

void AnalyzeSystematicsFromDataCards(string inputcardfile = "/data/smurf/data/cards/HWW2l/Full2011/124/hwwof_0j_cut.txt")
{
  
  ifstream ifs;
  ifs.open(inputcardfile.c_str());
  assert(ifs.is_open()); 
  
  string line; 
  
  cout << "Processing file: " << inputcardfile << endl;

  double SignalRate;
  vector<string> processLabel;
  vector<double> processRate;
  vector<double> processRateUncertaintySqr;
  vector<double> processRateLeadingUncertainty;
  vector<string> processRateLeadingUncertaintyLabel;

  for (UInt_t i=0; i < 12; ++i) {
    processRateUncertaintySqr.push_back(0.0);
    processRateLeadingUncertainty.push_back(0.0);
    processRateLeadingUncertaintyLabel.push_back("N/A");
  }


  getline(ifs,line);
  getline(ifs,line);
  getline(ifs,line);
  getline(ifs,line);
  getline(ifs,line);

  //Read Process Labels
  getline(ifs,line);  
  stringstream ss1(line);  
  string tmp;
  ss1 >> tmp;
  for (UInt_t i=0; i < 12; ++i) {
    ss1 >> tmp;
    processLabel.push_back(tmp);
    //cout << tmp << endl;
  }

  getline(ifs,line);

  //Read Process Rates
  getline(ifs,line);  
  stringstream ss2(line);  
  ss2 >> tmp;
//   cout << line << endl;
//   cout << tmp << endl;
  for (UInt_t i=0; i < 12; ++i) {
    double tmpValue;
    ss2 >> tmpValue;
    processRate.push_back(tmpValue);
//     cout << tmpValue << endl;
  }
  assert(processRate.size() == processRateUncertaintySqr.size());

  //Read all systematics
  Bool_t DoneReadingFile = kFALSE;
  getline(ifs,line);
  while (ifs.good()) {

//     cout << line << endl;

    stringstream ss3(line);  
    string SystematicLabel;
    ss3 >> SystematicLabel;
    ss3 >> tmp;

//      cout << SystematicLabel << endl;

    for (UInt_t i=0; i < 12; ++i) {
      double tmpValue = 0.0;
//      ss3 >> tmpValue;
      ss3 >> tmp;
      if (tmp != "-") tmpValue = fabs( atof(tmp.c_str()) - 1.0 );
//       cout << tmp << " : " << tmpValue << "  ";
//       cout << tmpValue*processRate[i] << endl;

      processRateUncertaintySqr[i] += pow(tmpValue*processRate[i],2);
      if (tmpValue*processRate[i] > processRateLeadingUncertainty[i]) {
        processRateLeadingUncertainty[i] = tmpValue*processRate[i];
        processRateLeadingUncertaintyLabel[i] = SystematicLabel;
      }
    }

    getline(ifs,line);
    if (!ifs.good()) DoneReadingFile = kTRUE;
  }

  SignalRate = processRate[0] + processRate[1] + processRate[2] + processRate[3];
  
//   cout << "Signal Rate: " << SignalRate << endl;
//   cout << "--------------------------------------------------" << endl;
//   double TotalBkgUncertaintySqr = 0;
//   for (UInt_t i=4; i < 12; ++i) {
//     cout << "Process: " 
//          << setw(10) << left << processLabel[i]
//          << setw(10) << left << processRate[i]
//          << setw(10) << left << TMath::Sqrt(processRateUncertaintySqr[i])
//          << setw(10) << left << processRateLeadingUncertaintyLabel[i] 
//          << "(" << processRateLeadingUncertainty[i] << ")" 
//          << endl;
//     TotalBkgUncertaintySqr += processRateUncertaintySqr[i];
//   }
//   cout << "Total Bkg Uncertainty : " << TMath::Sqrt(TotalBkgUncertaintySqr) << endl;



  cout << "| *Signal Rate* | *" << SignalRate << "* | | | |"
       << endl;
  cout << "| *Bkg Process* | *Rate* | *Rate Uncertainty* | *Leading Systematic* | *Leading Systematic Uncertainty* |" << endl;
  double TotalBkgUncertaintySqr = 0;

  char buffer[200];
  for (UInt_t i=4; i < 12; ++i) {
    cout << "| " 
         << processLabel[i] << " | ";

    sprintf(buffer,"%.2f",processRate[i]);
    cout << buffer << " | ";
    sprintf(buffer,"%.2f",TMath::Sqrt(processRateUncertaintySqr[i]));
    cout << buffer << " | ";
    cout << processRateLeadingUncertaintyLabel[i] << " | ";
    sprintf(buffer,"%.2f",processRateLeadingUncertainty[i]);
    cout << buffer << " |";
    cout << endl;
    TotalBkgUncertaintySqr += processRateUncertaintySqr[i];
  }

  sprintf(buffer,"%.2f",TMath::Sqrt(TotalBkgUncertaintySqr));
  cout << "| *Total Bkg Uncertainty* | *" << buffer << "* | | | |"
       << endl;



  ifs.close();
    
}
