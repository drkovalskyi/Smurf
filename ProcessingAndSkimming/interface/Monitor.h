#ifndef SMURF_PROCESSINGANDSKIMMING_MONITOR_H
#define SMURF_PROCESSINGANDSKIMMING_MONITOR_H
#include <vector>
#include <string>
#include "Smurf/Core/SmurfTree.h"

struct MonitorEventId {
  unsigned long int run, event, lumi;
  // --------------------------------------------------------------- //
  MonitorEventId(const SmurfTree&);
  MonitorEventId();
  bool operator < (const MonitorEventId& id) const{
    if (run != id.run) return run < id.run;
    if (lumi != id.lumi) return lumi < id.lumi;
    return event < id.event;
  }
  bool operator == (const MonitorEventId& id) const{
    return (run==id.run) && (lumi==id.lumi) && (event==id.event);
  }
  bool operator != (const MonitorEventId& id) const{
    return ! operator == (id);
  }
};

struct Entry {
  unsigned int nhyp[5];
  unsigned int nevt[5];
  double nhyp_weighted[5];
  double nevt_weighted[5];
  bool seen[5];
  MonitorEventId lastEvent;
  std::vector<MonitorEventId> events;
  std::string name;
  // -------------------------------------------------------------- //
  Entry();
};

struct Monitor{
  void count(const SmurfTree&, SmurfTree::Type type, const char* name, double weight=1.0);
  void print() const;
  void makeHistograms(const char* prefix) const;
  Monitor(bool iKeepEventList=false):nEvtProcessed(0),
                                          keepEventList(iKeepEventList){}
  // -------------------------------------- //
  std::vector<Entry> counters;
  unsigned int nEvtProcessed;
  bool keepEventList;
};
#endif
