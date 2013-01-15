#ifndef PROCESSFAKERATE_H
#define PROCESSFAKERATE_H

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

enum DataType {
    WJETS,
    DY,
    WZ,
    ZZ,
    DATA,
};

void printline(FILE* fout, TH2F* h2, const char* caption, bool doInt);
void format(TH1F *hist, DataType dataType);
void printStack(TFile *f, const char* name, const char* fname, Float_t SF_WJets, Float_t Err_WJets, Float_t SF_DY, Float_t Err_DY, bool ss=false);
void printZStack(TFile *f, const char* name, const char* fname, Float_t SF_EWK, Float_t Err_EWK);
void compareSS(TFile *f, const char* name, const char* fname);

void printValidationHistograms(const char* file, const char* name, const char* bin);
void printFakeRate(const char* file, const char *name, const char* friendlyName, const std::vector<unsigned int>& ptThresholds);


void printZFakeRate(const char* file, const char *name, const char* friendlyName);


// get scale factor for W+jets MC
void getWJetsScaleFactor(TFile *f, const char* name, float &sf, float &err);

// get the DY scale factor from DY MC
void getDYScaleFactor(TFile *f, const char* name, float &sf, float &err);


void setUncertainty(TH1F *hist, const float &err);
void setUncertainty(TH2F *hist, const float &err);


#endif

