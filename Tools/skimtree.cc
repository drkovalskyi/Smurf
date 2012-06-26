#include <iostream>
#include <cstdlib>
#include <string>
#include "TTree.h"
#include "../Core/SmurfTree.h"

bool passed(const SmurfTree& tree){
  return std::min(tree.pmet_,tree.pTrackMet_)>20;
}

int main(int argc, char** argv) {
  if (argc != 3){
    printf("Usage:\n\tskimtree <input_file> <output_file>\n");
    exit(1);
  }
  std::string input(argv[1]);
  std::string output(argv[2]);

  int i_permille_old = 0;
  SmurfTree tree;
  tree.LoadTree(input.c_str());
  tree.InitTree();
  TFile* f = TFile::Open(output.c_str(),"RECREATE");
  assert(f);
  TTree* newtree = tree.tree_->CloneTree(0);
  Long64_t nEntries = tree.tree_->GetEntries();
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
    if ( passed(tree) ) newtree->Fill();
  }
  newtree->Write();
  f->Close();
  return 0;
}
/*
  // int shift = 1;
  std::string tmp;
  for (int it = 1; it < argc; it++) {
    tmp = std::string(argv[it]);
    // print list of bad files                                                                                                                                                                       
    if (tmp == "-b") {
      doPrintBad = true;
      ++shift;
      continue;
    }
    // print list of good files                                                                                                                                                                      
    if (tmp == "-g") {
      doPrintGood = true;
      ++shift;
      continue;
    }
    // move bad files                                                                                                                                                                                
    if (tmp == "-x") {
      doMove = true;
      ++shift;
      continue;
    }
    // prefix rfio                             


#include "TSystem.h"
#include "Tools.h"

int main(int argc, char** argv) {
  std::vector<std::string> printbad;
  std::vector<std::string> printgood;


void SmurfAnalysisNew::addSample(SmurfTree::DataType type, const char* pattern){
  // FIXME: add proper handling of patterns
  samples[type].push_back(pattern);
}

void SmurfAnalysisNew::addProcessor(Processor* process){
  processors.push_back(process);
}

void SmurfAnalysisNew::addDefaultProcessors(){
  for ( SampleMap::const_iterator
	  sample = samples.begin(); sample != samples.end(); ++sample )
    {
      SmurfTree::DataType type = sample->first;
      bool foundProcessor(false);
      for ( std::vector<Processor*>::const_iterator itr=processors.begin(); itr!= processors.end(); ++itr )
	for ( std::vector<SmurfTree::DataType>::const_iterator t = (*itr)->m_types.begin(); t!=(*itr)->m_types.end(); ++t )
	  if ( type == *t ) foundProcessor = true;
      
      if (!foundProcessor) processors.push_back( new Processor(type) );
    }
}


const Double_t SmurfAnalysisNew::ptBins[SmurfAnalysisNew::nPtPoints] = {10., 15., 20., 25., 30., 35.};
const Double_t SmurfAnalysisNew::etaBins[SmurfAnalysisNew::nEtaPoints] = {0.0, 1.0, 1.479, 2.0, 2.5};
const char*    SmurfAnalysisNew::types[4] = {"mm","me","em","ee"};

void SmurfAnalysisNew::processSamples(){
  // add default processors for samples that don't have any
  // all samples are relevant
  for ( SampleMap::const_iterator
	  sample = samples.begin(); sample != samples.end(); ++sample )
    {
      SmurfTree::DataType type = sample->first;
      const std::vector<std::string>& files = sample->second;
      for (std::vector<std::string>::const_iterator file = files.begin();
	   file != files.end(); ++file)
	{
	  FileStat_t buf;
	  // check if file exists
	  if (gSystem->GetPathInfo(file->c_str(), buf)) {
	    printf("ERROR: file not found\n\t%s\n", file->c_str());
	    continue;
	  }
	  printf("processing file %s\n", file->c_str());
  
	  int i_permille_old = 0;
	  SmurfTree tree;
	  tree.LoadTree(file->c_str());
	  tree.InitTree();
	  Long64_t nEntries = tree.tree_->GetEntries();
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
	    if ( type == SmurfTree::data && m_json )
	      if( !goodrun(tree.run_, tree.lumi_) ) continue;
	    double weight = 1.0;
	    if ( type != SmurfTree::data ) weight = tree.scale1fb_*m_lumi;
	    for ( std::vector<Processor*>::iterator processor = processors.begin();
		  processor != processors.end(); ++processor )
	      {
		(*processor)->processEvent(tree, type, higgs_type_, m_selectionType, weight);
	      }
	  }
	}
    }
  printf("Done processing samples\n");
}

void SmurfAnalysisNew::showYields(YieldType type,double lumi)
{
  const char* pm = "+/-";
  const char* groups[3] = {"Data","Standard Model","Higgs"};
  double weight = 1;
  if (lumi>0) weight = lumi/m_lumi;
  Yield total;
  for ( unsigned int i=0; i<3; ++i ){
    printf("\n%s process event yields:\n",groups[i]);    
    printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);

    for ( std::vector<Processor*>::iterator processor = processors.begin();
	  processor != processors.end(); ++processor )
      {
	SmurfTree::DataType sample = (*processor)->m_types.front();
	if ( sample==SmurfTree::data ) {
	  if ( i!=0 ) continue;
	} else {
	  if ( isStandardModel(sample) ) {
	    if (i!=1) continue;
	  } else {
	    if (i!=2) continue;
	  }
	}
	std::vector<Yield> yields = (*processor)->getYields(type);
	for ( std::vector<Yield>::const_iterator yield = yields.begin(); yield != yields.end(); ++yield ){
	  printf("%s", yield->name.c_str());
	  if ( yield->type == SmurfTree::data ){
	    for (unsigned int j=0; j<4; ++j)
	      printf(" \t   %5.0f    ", yield->yield[j]*weight);
	    printf(" \t    %5.0f\n", yield->total_yield()*weight);
	  } else {
	    for (unsigned int j=0; j<4; ++j)
	      printf(" \t%5.2f%s%4.2f", yield->yield[j]*weight, pm, sqrt(yield->err2[j])*weight);
	    printf(" \t%5.2f%s%4.2f\n", yield->total_yield()*weight, pm, yield->total_error()*weight);
	  }      
	  if ( isStandardModel(yield->type) ) total.add(*yield, weight);
	}
      }
    if (i==1) {
      printf("Total:");
      for (unsigned int j=0; j<4; ++j)
	printf(" \t%5.2f%s%4.2f", total.yield[j], pm, sqrt(total.err2[j]));
      printf(" \t%5.2f%s%4.2f\n", total.total_yield(), pm, total.total_error());
    }
  }
  
//   printf("\nHiggs expected event yields:\n");;    
//   printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);
//   for ( ProcessorMap::iterator process = processes.begin();
// 	process != processes.end(); ++process )
//     {
//       SmurfTree::DataType sample = process->first;
//       Yield yield = process->second->getYield();
//       if (isStandardModel(sample)||sample==SmurfTree::data) continue;
//       printf("%s", SmurfTree::name(sample).c_str());
//       for (unsigned int j=0; j<4; ++j)
// 	printf(" \t%5.2f%s%4.2f", yield->yield[j]*weight, pm, sqrt(yield->err2[j])*weight);
//       printf(" \t%5.2f%s%4.2f\n", yield->total_yield()*weight, pm, yield->total_error()*weight);
//     }
}
*/
