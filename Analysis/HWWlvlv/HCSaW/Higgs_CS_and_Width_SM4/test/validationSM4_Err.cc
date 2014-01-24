#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "HiggsCSandWidthSM4.cc"

using namespace std;

int main()
{

  ofstream fileOut;
  char* fileName = "HiggsSM4Systematics.h";
  fileOut.open(fileName);

  const Int_t nPoints = 145;
  double sqrts = 7;
  int ID = 10;
  double mass[nPoints] = {110.0,111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,120.0,121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,136.0,137.0,138.0,139.0,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,162.0,164.0,166.0,168.0,170.0,172.0,174.0,176.0,178.0,180.0,182.0,184.0,186.0,188.0,190.0,192.0,194.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,295.0,300.0,305.0,310.0,315.0,320.0,325.0,330.0,335.0,340.0,345.0,350.0,360.0,370.0,380.0,390.0,400.0,420.0,440.0,450.0,460.0,480.0,500.0,520.0,540.0,550.0,560.0,580.0,600.0};
  HiggsCSandWidthSM4 *myCSW = new HiggsCSandWidthSM4();
  double a[nPoints],b[nPoints],c[nPoints];

  for( int i=0; i<nPoints; i++) {
    a[i] = myCSW->HiggsBRErr_Hff    (ID, mass[i], sqrts);
    b[i] = myCSW->HiggsBRErr_HVV    (ID, mass[i], sqrts);
    c[i] = myCSW->HiggsBRErr_Hgluglu(ID, mass[i], sqrts);
  }
  fileOut << "Double_t HiggsSM4Systematics_HiggsBRErr_Hff(Double_t mH){" << endl;
  fileOut << "  Double_t mass[" << nPoints << "]={";
  for( int i=0; i<nPoints; i++) {fileOut << mass[i]; if(i!=nPoints-1) fileOut << ",";}
  fileOut << "};" << endl;
  fileOut << "  Double_t val[" << nPoints << "]={";
  for( int i=0; i<nPoints; i++) {fileOut << a[i]; if(i!=nPoints-1) fileOut << ",";}
  fileOut << "};" << endl;
  fileOut << "  Int_t massIndex = -1;" << endl;
  fileOut << "  for (UInt_t m=0; m < " << nPoints << " ; ++m) {" << endl;
  fileOut << "    if (mH == mass[m]) {massIndex = m; break;} " << endl;
  fileOut << "  }" << endl;
  fileOut << "  assert(massIndex >= 0);" << endl;
  fileOut << "  return val[massIndex];" << endl;
  fileOut << "}" << endl;

  fileOut << "Double_t HiggsSM4Systematics_HiggsBRErr_HVV(Double_t mH){" << endl;
  fileOut << "  Double_t mass[" << nPoints << "]={";
  for( int i=0; i<nPoints; i++) {fileOut << mass[i]; if(i!=nPoints-1) fileOut << ",";}
  fileOut << "};" << endl;
  fileOut << "  Double_t val[" << nPoints << "]={";
  for( int i=0; i<nPoints; i++) {fileOut << b[i]; if(i!=nPoints-1) fileOut << ",";}
  fileOut << "};" << endl;
  fileOut << "  Int_t massIndex = -1;" << endl;
  fileOut << "  for (UInt_t m=0; m < " << nPoints << " ; ++m) {" << endl;
  fileOut << "    if (mH == mass[m]) {massIndex = m; break;} " << endl;
  fileOut << "  }" << endl;
  fileOut << "  assert(massIndex >= 0);" << endl;
  fileOut << "  return val[massIndex];" << endl;
  fileOut << "}" << endl;

  fileOut << "Double_t HiggsSM4Systematics_HiggsBRErr_Hgluglu(Double_t mH){" << endl;
  fileOut << "  Double_t mass[" << nPoints << "]={";
  for( int i=0; i<nPoints; i++) {fileOut << mass[i]; if(i!=nPoints-1) fileOut << ",";}
  fileOut << "};" << endl;
  fileOut << "  Double_t val[" << nPoints << "]={";
  for( int i=0; i<nPoints; i++) {fileOut << c[i]; if(i!=nPoints-1) fileOut << ",";}
  fileOut << "};" << endl;
  fileOut << "  Int_t massIndex = -1;" << endl;
  fileOut << "  for (UInt_t m=0; m < " << nPoints << " ; ++m) {" << endl;
  fileOut << "    if (mH == mass[m]) {massIndex = m; break;} " << endl;
  fileOut << "  }" << endl;
  fileOut << "  assert(massIndex >= 0);" << endl;
  fileOut << "  return val[massIndex];" << endl;
  fileOut << "}" << endl;

  fileOut.close();
 return 0;
}
