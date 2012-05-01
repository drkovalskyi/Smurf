#include "TH2.h"
#include "TMatrixTSym.h"
#include <math.h>
#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TList.h"
#include "TStyle.h"

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
  cc->Update();

//   TPaletteAxis* paletteAxis = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject( "palette" );
//   if ( paletteAxis ){
//     paletteAxis->SetLabelSize( 0.03 );
//     paletteAxis->SetX1NDC( paletteAxis->GetX1NDC() + 0.02 );
//   }    
//   hist->Draw("sametext");
}

