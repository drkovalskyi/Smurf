// --------------------------------------------------------
// This macro strips off uninteresting columns and rows 
//  * makes N x N TH2 into n x n TH2
//	* currently works with output of combine fits   
//	* It first removes boring X slices and then Y slices 
//	* only work with cases of N x N ---> n x n  
//  * need to set number of interesting nuisance parameters ()
//  * 2D histogram name convention 
//     - h2in    : N x N 
//     - h2outX  : n x N 
//     - h2outXY : n x n 
//
//  Usage : root[] .x strip_histo.C("mlfit.root") 
//
// --------------------------------------------------------

//
// List of nuisance parameters to be drawn  
//
const int nNui = 5;
char *nuisance_out[nNui] = {"CMS_hww_1j_WW_8TeV", "CMS_hww_1j_ttbar_8TeV", "CMS_hww_Wg3l", "FakeRate", "QCDscale_WW_EXTRAP"}; 

//
// This function compare a string(nui_histo) with a list of strings (nuisance_out)
// If nui_histo exists in nui_list[nNui], it returns true.
// If not, it returns false.
//
bool comparenames(char* nui_list[nNui], char* nui_histo) { 

	bool match = false;

	for(int i=0; i<nNui; i++) { 
		TString str_nui_list = nui_list[i] ;
		TString str_nui_histo = nui_histo ;
		if( str_nui_histo.Contains(str_nui_list) ) match = true; 
	}
	return match;
}

//
// B/W color scheme for the color blind
// from Dave who got from Dima
//
void loadColorPalette() { 
	
	Int_t palette[7];
	palette[3] = kWhite;
	for (unsigned int i=0;i<3;++i){
		palette[2-i] = kGray+i;
		palette[4+i] = kGray+i;
	}
	gStyle->SetPalette(7,palette);

}

//
// Main function
//
void strip_histo(TString infile) {
	
	// load B/W color scheme from Dave who got from Dima 	
	loadColorPalette();

	// Style  
  	gStyle->SetPaintTextFormat("6.2f"); 

   	// input file
	TFile *File = TFile::Open(infile, "READ"); 
	// get the covariance matrix(TH2 histogram)
	// fit_s is a result of a fit with mu floating
	RooFitResult *fit_s = (RooFitResult*) File->Get("fit_s"); 
    TH2F *h2in   = (TH2F*) fit_s->correlationHist(); 
	
	//
	// strip off unnecessary X slices  
	//
	//  N x N --> n x N

	// Axes of N x N
	TAxis *axisX = h2in->GetXaxis(); 	 
	TAxis *axisY = h2in->GetYaxis();	
	int nbinsX = axisX->GetNbins();	
	int nbinsY = axisY->GetNbins();	 
	
	// output histogram after stripping unnecessary X slices
	TH2F *h2outX 	= new TH2F("h2outX", "h2outX", nNui,0,nNui, nbinsY,0,nbinsY);
	int h2outXbin = 1; // bin index for h1outX
	
	for (int ibinX = 1; ibinX <= nbinsX; ibinX++) { 
		if( !comparenames(nuisance_out, axisX->GetBinLabel(ibinX)) ) continue;
			
		// project an X slice
		TH1D *h1Xslice = h2in->ProjectionY("h1Xslice",ibinX,ibinX);	 
		h2outX->GetXaxis()->SetBinLabel(h2outXbin, axisX->GetBinLabel(ibinX));
		// put them into 2D 
		for (int ibinY = 1; ibinY <= nbinsY; ibinY++) {  
			h2outX->SetBinContent(h2outXbin, ibinY, h1Xslice->GetBinContent(ibinY)); 
			h2outX->GetYaxis()->SetBinLabel(ibinY, axisY->GetBinLabel(ibinY));
		}
		h2outXbin++;
	} 
	
	//
	// strip off unnecessary Y slices  
	//
	//  n x N --> n x n 
	
	// Axes of n x N
	TAxis *axisXoutX = h2outX->GetXaxis(); 	 
	TAxis *axisYoutX = h2outX->GetYaxis(); 	  
	int nbinsXoutX = axisXoutX->GetNbins();	
	int nbinsYoutX = axisYoutX->GetNbins();	 
	
	// output histogram after stripping unnecessary X slices
	TH2F *h2outXY 	= new TH2F("h2outXY", "h2outXY", nNui,0,nNui, nNui,0,nNui);
	int h2outYbin = 1; // bin index for h1outXY
	
	for (int ibinY = 1; ibinY <= nbinsYoutX; ibinY++) { 
		if( !comparenames(nuisance_out, axisYoutX->GetBinLabel(ibinY)) ) continue;
			
		// project an Y slice
		TH1D *h1Yslice = h2outX->ProjectionX("h1Yslice",ibinY,ibinY);	  
				
		h2outXY->GetYaxis()->SetBinLabel(h2outYbin, axisYoutX->GetBinLabel(ibinY));
		// put them into 2D 
		for (int ibinX = 1; ibinX <=nbinsXoutX; ibinX++) {  
			h2outXY->SetBinContent(ibinX, h2outYbin, h1Yslice->GetBinContent(ibinX)); 
			h2outXY->GetXaxis()->SetBinLabel(ibinX, axisXoutX->GetBinLabel(ibinX));
		}
		h2outYbin++;
	}

	// 
	// Draw and save histogram
	// 
	TCanvas* c = new TCanvas("c","",800,600);
    c->SetLeftMargin(0.25);
    c->SetBottomMargin(0.2);
    c->SetRightMargin(0.2);
	c->SetGridx();	
	c->SetGridy();	

	h2outXY->SetTitle("Covariance");	
	h2outXY->SetMinimum(-1.);	
	h2outXY->SetMaximum(1.);	
	h2outXY->SetMarkerSize(2);	
	h2outXY->SetStats(0);	
	h2outXY->Draw("colz text");	 
	
	c->Print("covairance.pdf");
	
}
