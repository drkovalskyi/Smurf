void SetBinContentError(TH1F* &h1, int index, float content=0, float error=0){
    h1->SetBinContent(index, content);
    h1->SetBinError(index, error);
}

void makecardfortoy(int TEST=0, int toynumber=1, char* TESTNAME){

    // options
    bool doNominal      = 0;
    bool doNewDefault   = 0;
    bool doCR1          = 0;
    bool doCR2          = 0;
    if(TEST==0) doNominal       = 1;
    if(TEST==1) doNewDefault    = 1;
    if(TEST==2) doCR1           = 1;
    if(TEST==3) doCR2           = 1; 
    if (doNominal+doNewDefault+doCR1+doCR2>1) {
        cout << "!!!! Aborted :: multiple tests selected !!!! "<< endl; 
        return;
    }

    char* test=""; 
    if(doNewDefault)   test = "_newdefault";
    if(doCR1)       test = "_CR1";
    if(doCR2)       test = "_CR2";
    
    char* toyrootfile = "FullqqWWSyst_PseudoData_sb_seed12345.root"; 
    char* card = Form("125/hwwof_0j_shape_8TeV%s.txt", test); 

    // Get obs from a pseudo data
    TFile   *toyRootFile = TFile::Open(toyrootfile); 
    TH1F    *h1_data = (TH1F*)toyRootFile->Get(Form("j0of_%i", toynumber));  
    h1_data->SetTitle("histo_new_Data");
    h1_data->SetName("histo_new_Data");

    //
    // chop off the unnecessary region
    //
    bool removethisbin = false; 
    for(int i=1; i<127; i++){

        // remove CR2 : mll>60 && mT<120
        if ( doCR1 &&
                ( (i>=4 && i<=9)    || (i>=13 && i<=18) || (i>=22 && i<=27)
                  || (i>=31 && i<=36)  || (i>=40 && i<=45) || (i>=49 && i<=54) )
           )  removethisbin = true;

        // remove CR1 : mT>120
        if ( doCR2 &&
                i>=55
           )  removethisbin = true;
        
        if(removethisbin) SetBinContentError(h1_data, i); 

        removethisbin = false;
    } 
    
    gSystem->Exec(Form("cp 125/hwwof_0j.input_8TeV%s.root 125/toycards_%s/hwwof_0j.input_8TeV%s_%i.root", test, TESTNAME, test, toynumber));
    TFile *inputrootfile = new TFile(Form("125/toycards_%s/hwwof_0j.input_8TeV%s_%i.root", TESTNAME, test, toynumber), "UPDATE");
    gROOT->cd();
    inputrootfile->cd();
    h1_data->SetDirectory(0); h1_data->Write();
    inputrootfile->Close();

    unsigned int obs = h1_data->Integral();
    
    //
    // print out a new card
    // 
    ofstream fout;
    fout.open(Form("125/toycards_%s/hwwof_0j_shape_8TeV%s_%i.txt", TESTNAME,test,toynumber));

    string line;
    ifstream infile (card);
    if (infile.is_open())
    {
        while ( infile.good() )
        {
            // get a line from input file
            getline (infile,line);
            
            if( line.find("Observation")!=string::npos )  {
                fout << "Observation  " << obs << endl;    
            } else if( line.find("shapes ")!=string::npos )  {
                TString tstringline = line;
                tstringline.ReplaceAll(Form("125/hwwof_0j.input_8TeV%s.root", test), 
                                        Form("125/toycards_%s/hwwof_0j.input_8TeV%s_%i.root",TESTNAME,test,toynumber)); 
                tstringline.ReplaceAll("histo_Data","histo_new_Data"); 
                fout << tstringline << endl;
            } else if( line.find("bin")!=string::npos )  {
                fout << "bin j0of j0of j0of j0of j0of j0of j0of j0of j0of j0of j0of j0of j0of j0of" << endl;
            } else {
                fout << line << endl;
            }



            if( !infile.good() ) continue;

        }
    }
    fout.close();
    infile.close();
    toyRootFile->Close();
}

void makecardfortoy_runAll(int TOTALTOY=10, char* TESTNAME) {

    for(int itoy=0; itoy<TOTALTOY; itoy++) { 
        
        cout << "... Making cards using toy "<< itoy << endl;
        
        makecardfortoy(0,itoy,TESTNAME); // nominal 
        makecardfortoy(1,itoy,TESTNAME); // new default 
        makecardfortoy(2,itoy,TESTNAME); // CR1
        makecardfortoy(3,itoy,TESTNAME); // CR2
    }
}
