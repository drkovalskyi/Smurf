#include "Smurf/ProcessingAndSkimming/interface/Monitor.h"
#include <fstream>
#include "TH1.h"

MonitorEventId::MonitorEventId(const SmurfTree& tree){
  run = tree.run_;
  event = tree.event_;
  lumi = tree.lumi_;
}

MonitorEventId::MonitorEventId(){
  run = 0;
  event = 0;
  lumi = 0;
}

Entry::Entry()
{
  for (unsigned int i=0; i<5; ++i){
    nhyp[i] = 0;
    nevt[i] = 0;
    seen[i] = false;
  }
}

void Monitor::count(const SmurfTree& tree, SmurfTree::Type type, const char* name, double weight)
{
  std::vector<Entry>::iterator itr = counters.begin();
  while (itr != counters.end() && itr->name != name) itr++;
  Entry* entry(0);
  if ( itr == counters.end() ){
    counters.push_back(Entry());
    entry = &counters.back();
    entry->name = name;
  } else {
    entry = &*itr;
  }
  MonitorEventId id(tree);
  entry->nhyp[type]++;
  entry->nhyp[4]++;
  if (id != entry->lastEvent){
    if (keepEventList) entry->events.push_back(id);
    for (unsigned int i=0; i<5; ++i) entry->seen[i] = false;
    entry->nevt[type]++;
    entry->nevt[4]++;
    entry->nevt_weighted[type]+=weight;
    entry->nevt_weighted[4]+=weight;
    entry->seen[type] = true;
    entry->lastEvent = id;
  } else {
    if ( !entry->seen[type] ){
      entry->nevt[type]++;
      entry->nevt_weighted[type]+=weight;
      entry->seen[type] = true;
    }
  }
}

std::string HypothesisTypeName(unsigned int i){
  switch(i){
  case SmurfTree::mm:
    return "mm";
  case SmurfTree::ee:
    return "ee";
  case SmurfTree::me:
    return "me";
  case SmurfTree::em:
    return "em";
  default:
    return "all";
  }
}

void Monitor::print() const
{
  std::cout << "Total number of processed events: \t" << nEvtProcessed << std::endl;
  unsigned int order[5] = {SmurfTree::mm,SmurfTree::me,SmurfTree::em,SmurfTree::ee,4};

  std::cout << HypothesisTypeName(order[0]) << " / " << HypothesisTypeName(order[1]) << "/" << HypothesisTypeName(order[2]) << "/" <<
    HypothesisTypeName(order[3]) << "/" << HypothesisTypeName(order[4]) << std::endl;
  for ( unsigned int i=0; i<counters.size(); ++i ){
    std::cout << Form("%-40s \thyps: %u/%u/%u/%u/%u \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		      counters[i].nhyp[order[0]],counters[i].nhyp[order[1]],counters[i].nhyp[order[2]],counters[i].nhyp[order[3]],counters[i].nhyp[order[4]],
		      counters[i].nevt[order[0]],counters[i].nevt[order[1]],counters[i].nevt[order[2]],counters[i].nevt[order[3]],counters[i].nevt[order[4]])      
	      << std::endl;
    std::ofstream cut_file(Form("cut-%d.txt",i));
    cut_file << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		     counters[i].nevt[order[0]],counters[i].nevt[order[1]],counters[i].nevt[order[2]],counters[i].nevt[order[3]],counters[i].nevt[order[4]]) << "\n";
    for ( std::vector<MonitorEventId>::const_iterator id=counters[i].events.begin();
	  id!=counters[i].events.end(); ++id ){
      cut_file << id->run << "\t" << id->lumi << "\t" << id->event <<"\n";
    }
    cut_file.close();
  }
  if (nEvtProcessed>0){
    for ( unsigned int i=0; i<counters.size(); ++i ){
      std::cout << Form("%-40s \thyps: %f/%f/%f/%f/%f \tnevts: %f/%f/%f/%f/%f", counters[i].name.c_str(),
			counters[i].nhyp[order[0]]/((float)nEvtProcessed),counters[i].nhyp[order[1]]/((float)nEvtProcessed),
			counters[i].nhyp[order[2]]/((float)nEvtProcessed),counters[i].nhyp[order[3]]/((float)nEvtProcessed),
			counters[i].nhyp[order[4]]/((float)nEvtProcessed),
			counters[i].nevt[order[0]]/((float)nEvtProcessed),counters[i].nevt[order[1]]/((float)nEvtProcessed),
			counters[i].nevt[order[2]]/((float)nEvtProcessed),counters[i].nevt[order[3]]/((float)nEvtProcessed),
			counters[i].nevt[order[4]]/((float)nEvtProcessed)) 
		<< std::endl;
    }
  }
}

void Monitor::makeHistograms(const char* prefix) const
{
  TH1F* hist[4];
  TH1F* histw[4];
  for (unsigned int i=0; i<4; i++){
    hist[i]  = new TH1F(Form("%s_hcuts_%s", prefix, HypothesisTypeName(i).c_str()), 
			"Number of events vs cuts", counters.size(), 0, counters.size() );	
    histw[i] = new TH1F(Form("%s_hcuts_weighted_%s", prefix, HypothesisTypeName(i).c_str()), 
		       "Number of weighted events vs cuts", counters.size(), 0, counters.size() );	
    for ( unsigned int j=0; j<counters.size(); ++j ){
      hist[i]->SetBinContent(j+1,counters[j].nevt[i]);
      hist[i]->GetXaxis()->SetBinLabel(j+1,counters[j].name.c_str());
      histw[i]->SetBinContent(j+1,counters[j].nevt_weighted[i]);
      histw[i]->GetXaxis()->SetBinLabel(j+1,counters[j].name.c_str());
    }
  }
}
