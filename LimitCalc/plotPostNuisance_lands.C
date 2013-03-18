//
// Script that visualizes the post-fit nuisance parameters
//
// Usage :
//  root -b plotPostNuisance_combine.C'("lands.log")' 
// 
//  Arguments : 
//    * lands.log : log of ML fit by LandS



void printaline(int index, char* nuiname, Double_t centralval, Double_t uncert){

    cout.width(5);  cout << index << " ::: ";
    cout.width(40); cout << nuiname << " : ";
    cout.width(10);  cout << Form("%1.3f", centralval);
    cout.width(5);  cout << " +/- ";
    cout.width(10);  cout << Form("%1.3f", uncert) << endl;

}

void plotPostNuisance_lands(char* log="ana_Moriond13_2D_V1/output/limits_0j_shape_of_125-fit-125-all.log") {

    //
    // Containers for post-fit nuisances
    //
    std::vector<TString> nuisance;  TString nuisance_tmp;   // nuisance name 
    std::vector<float> central;     float central_tmp;      // post-fit central value  : normalized to the input uncertainty 
    std::vector<float> uncert;      float uncert_tmp;       // post-fit uncertainty    : normalized to the input uncertainty 
    TString tmp[11]; 
    bool SBfit = true;      // This is a flag to read only the first fit(S+B fit) done by GetFitResults.pl

    cout << "-------------------------------------------------------------------------------" << endl;
    cout << "Index :::               Nuisance name              :       Pull +/- uncertainty" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    // 
    // Read each line in the input log file and store 
    // post-fit results to the containers
    // 

    string line;
    ifstream inlog(log); 

    if (inlog.is_open())
    {
        while ( inlog.good() )
        { 
            // get a line from the input file
            getline (inlog,line);

            if( line.find("Real")!=string::npos ) SBfit=false; 
            if(!SBfit) continue;
            
            // skip line if it does not contain "rate" 
            if( line.find("par")==string::npos ) continue;
            if( line.find("name")!=string::npos ) continue;
            if( line.find("signal_strength")!=string::npos ) continue; 

            stringstream stream(line); 
            for(int i=0; i<11; i++) { 
                if (i==1) stream >> nuisance_tmp;
                else if (i==2) stream >> central_tmp;
                else if (i==4) stream >> uncert_tmp;
                else stream >> tmp[i];  
            }
            
            central.push_back(central_tmp);
            uncert.push_back(uncert_tmp);
            nuisance.push_back(nuisance_tmp);
           
            printaline( central.size(), nuisance_tmp, central_tmp, uncert_tmp);
            
            if( !inlog.good() ) continue;
        }
    }
    inlog.close();  
    cout << "-------------------------------------------------------------------------------" << endl;

    //
    // Visualize the result
    //
    // * gray area corresponds to post-fit cenral value and pre-fit uncertainty
    // * solid lines correspond to post-fit central value and post-fit uncertainty
    // Purpose of gray area is to compare pre/post-fit uncertainty and
    // check how better data constrains a nuisance than prediction
    TCanvas *c = new TCanvas("c", "c", 800, 500); 
 //   c->Divide(1,2);
  
    // store nuisances with pull >=0.5 
    // value is set ~20 lines below : if(TMath::Abs(central[i])>=0.5)
    int nchangednuisance = 0;
    std::vector<TString> changednuisance;

    // histogram
    c->cd(1);
    int nbins = central.size()+2;
    TH1D* h1        =   new TH1D("h1", "h1", nbins, -0.5, nbins-0.5);           // post-fit central value and post-fit uncertainty
    TH1D* h1_ref    =   new TH1D("h1_ref", "h1_ref", nbins, -0.5, nbins-0.5);   // post-fit central value and pre-fit uncertainty

    cout << endl;
    cout << " ****************************** LARGE VARIATION ****************************** " << endl;
    for(int i=0; i<nbins-2; i++) {
        h1->SetBinContent(i+2, central[i]);
        h1->SetBinError(i+2, uncert[i]); 
        h1_ref->SetBinContent(i+2, central[i]);
        h1_ref->SetBinError(i+2, 1.); 
        if(TMath::Abs(central[i])>=0.5) { 
            nchangednuisance++; 
            changednuisance.push_back(Form("Index(%2i)  %s : %1.2f +/- %1.2f", i+1, nuisance[i].Data(), central[i], uncert[i])); 
            
            printaline( i+1, nuisance[i].Data(), central[i], uncert[i]);
        }
    }
    cout << " ****************************************************************************** " << endl;
    h1_ref->SetStats(0);
    h1_ref->SetTitle("");
    h1_ref->SetYTitle("Pull");
    h1_ref->SetFillColor(18); 
    h1_ref->SetMinimum(-2.); 
    h1_ref->SetMaximum(2.); 
    //h1_ref->SetFillStyle(3004); 
    h1->SetLineColor(kBlack); 

    TH1D* h1_zero    =   new TH1D("h1_zero", "h1_zero", 1, -0.5, nbins-0.5);
    h1_zero->SetBinContent(1,0); 
    h1_zero->SetLineStyle(2); 

    h1_ref->Draw("E2");
    h1->Draw("same E1");
    h1_zero->Draw("histo same"); 

    TString plotname = log;
    plotname.ReplaceAll(".log", "_Postnuisance.pdf");
    c->SaveAs(Form("%s", plotname.Data()));

/*   
    // nuisances with large change 
    c->cd(2);
    TLatex l;
    l.DrawLatex(0.1, 0.9, "Nuicance : pull +/- postfit uncertainty");
    l.SetTextSize(0.04);

    for(int i=0; i<nchangednuisance; i++) {
        l.DrawLatex(0.1, 0.7-0.08*i, changednuisance[i]);
    }
*/

}
