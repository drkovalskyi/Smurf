#include <iomanip>
#include <sstream> 
#include <fstream>
#include <iostream>

#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom.h"

#include "SmurfTree.h"

/*
  TreeValidator:

  class to validate SmurfTree against each other. Events that fails the validation are divided 
  into three categories: 

  + category A: events that are present in the reference SmurfTree but not in the test SmurfTree.
  + category B: events that are present in the test SmurfTree but not in the reference SmurfTree.
  + category A: events that are present in both SmurfTrees but missmatch.

  The reference SmurfTree is defined in the constructor. At the same time logfiles are opened for 
  each event category. These files contain a headline with common information on the reference & 
  test SmurfTree and for each failing event in the corresponding category the event number, run 
  number and a set of variables, which describe the event. 

  During the validation process first the list of all available branches in the test SmurfTree is 
  tested against all available branches in ref_. Fails will issue Warnings but not stop the 
  further validation process. These Warnings will not be part of the logfiles. 

  The content validation is performed on all comparable SmurfTree variables (at the moment this
  does not yet include Lorentz vector types -- sorry). Each failing event will be dumped into the 
  logfile of the corresponding error category.

  At the end of the validation process a short summary is given on the number of events, which 
  are present in both SmurfTrees and the number of events that failed on category A or B.
*/

class TreeValidator{
public:
  /// we call Envelope a map min|mean|max values of a branch in a SmurfTree mapped against the 
  /// corresponding variables name 
  //typedef std::map<std::string, std::vector<double> >  Envelope;  

public:
  /// default constructor (takes ownership of reference SmurfTree); file is the name of the file, 
  /// that contains the reference SmurfTree; opens logfiles for error categories A, B and C
  TreeValidator(std::string file);
  //TreeValidator(std::string file);
  /// default destructor; deletes the reference SmurfTree and closes logfiles
  ~TreeValidator();
  /// global validation function, file is the name of the file containing the SmurfTree to be 
  /// validated against the reference tree 
  bool validate(std::string file);

private:
  /// fill map of events versus entry; thus the reference tree will automatically be ordered by 
  /// increasing event numbers
  void initEventList();
  /// fill list of reference branches; all available branches in the test SmurfTree will be compared 
  /// against all available branches in the reference SmurfTree
  void initBranchList();
  /// find the envolope for a given set of variables in a SmurfTree (see above what envelope is)
  std::map<std::string, std::vector<double> > envelope(SmurfTree* smurf, std::vector<std::string> variables);
  /// comapre available branches in smurf against ref_ 
  bool validateBranches(SmurfTree* smurf);
  /// compare content of all comparable branches in smurf against corresponding branches in ref_
  bool validateContent(SmurfTree* smurf);
  /// issue a standardized statement (by default labeled as INFO)
  void info(std::string issue, std::string=std::string("INFO"));
  /// issue a standardized warning 
  void warning(std::string issue) {info(issue, std::string("WARNING"));};
  /// write header to file (used for output files categoryA, categoryB and categoryC)
  void fileHeader(std::ofstream& file, std::string value, SmurfTree* smurf);
  /// write line to file
  void fileLine(std::ofstream& file, unsigned int idx=0, SmurfTree* smurf=0);
  /// check whether event is in eventList_ or not
  bool isInEventList(unsigned int event) const {return eventList_.find(event)!=eventList_.end();}
  /// retrieve entry of event in SmurfTree ref_
  unsigned int getEntry(unsigned int event) const { return isInEventList(event) ? eventList_.find(event)->second : -9999; }
  /// erase event from eventList_ value off event list; this is used to check the presence of events 
  /// in ref_ at the same tiem as the presence of events in the test SmurfTree
  void eraseFromEventList(unsigned int event) { if( isInEventList(event) ){ eventList_.erase(eventList_.find(event)); } }

private:
  /// reference SmurfTree
  SmurfTree* ref_;
  /// name of the file that contains the reference SmurfTree (needed for headers of logfiles) 
  std::string refName_;
  /// logfile for events that are present in ref_ but not in the test SmurfTree
  std::ofstream categoryA_; 
  /// logfile for events that are present in the test SmurfTree but not in ref_
  std::ofstream categoryB_; 
  /// logfile for events that are present in both SmurfTrees but missmatch in at least one of 
  /// the testable varibales
  std::ofstream categoryC_; 
  /// list of all available branches in ref_
  std::vector<std::string> branchList_;
  /// list of events in ref_; the first element (index) is 'event', the second element (value)
  /// is the entry in the TTree; this the events in eventList_get automatically ordered
  std::map<unsigned int, unsigned int> eventList_;
  /// increased verbosity for debuging
  bool debug_; 
};

TreeValidator::TreeValidator(std::string file):refName_(file), categoryA_("categoryA.txt"), categoryB_("categoryB.txt"), categoryC_("categoryC.txt"), debug_(true) 
{
  // load SmurfTree for ref_
  ref_= new SmurfTree();
  ref_->LoadTree(refName_.c_str()); ref_->InitTree();
  // init branch and event list
  initBranchList(); initEventList();
};

TreeValidator::~TreeValidator() 
{
  // close logfiles
  categoryA_.close();
  categoryB_.close();
  categoryC_.close();
  // cleanup
  delete ref_;
};

bool
TreeValidator::validate(std::string file)
{
  bool validated=true;
  // load SmurfTree to be tested (i.e. value)
  SmurfTree* value = new SmurfTree();
  value->LoadTree(file.c_str()); value->InitTree();
  
  // set logfile header
  fileHeader(categoryA_, file, value);
  fileHeader(categoryB_, file, value);
  fileHeader(categoryC_, file, value);

  // check that all branch names are the same
  validated |= validateBranches(value);
  // check that content is the same
  validated |= validateContent(value);  

  // clean up abd return
  delete value;
  return validated;
}

void 
TreeValidator::info(std::string issue, std::string level)
{
  std::cout << " --- " << level << ": " << issue << std::endl;
}

std::map<std::string, std::vector<double> >
TreeValidator::envelope(SmurfTree* smurf, std::vector<std::string> vars)
{
  std::map<std::string, std::vector<double> > envolope;
  std::vector<double> mins, means, maxes;
  unsigned int nevent = smurf->tree_->GetEntriesFast();
  for(unsigned int ientry=0; ientry<nevent; ++ientry){
    smurf->tree_->GetEntry(ientry);
    // loop variables in vars
    unsigned int ivar=0;
    for(std::vector<std::string>::const_iterator var=vars.begin(); var!=vars.end(); ++var, ++ivar){
      // fill min value
      if(mins.size()<=ivar ){ mins .push_back(smurf->Get(*var)); } 
      else { if(smurf->Get(*var)<mins .at(ivar)){ mins .at(ivar)=smurf->Get(*var); } }
      // fill max value
      if(maxes.size()<=ivar){ maxes.push_back(smurf->Get(*var)); } 
      else { if(smurf->Get(*var)>maxes.at(ivar)){ maxes.at(ivar)=smurf->Get(*var); } }
      // fill mean value
      if(means.size()<=ivar){ means.push_back(smurf->Get(*var)); } 
      else { means.at(ivar)+=1./ientry*smurf->Get(*var)-means.at(ivar); }
    }
  }
  // re-key min|mean|max against variables in vars
  unsigned int ivar=0;
  for(std::vector<std::string>::const_iterator var=vars.begin(); var!=vars.end(); ++var, ++ivar){
    std::vector<double> buffer; 
    buffer.push_back(mins .at(ivar)); 
    buffer.push_back(means.at(ivar)); 
    buffer.push_back(maxes.at(ivar)); 
    envolope[*var] = buffer;
  }
  return envolope;
}

void 
TreeValidator::fileLine(std::ofstream& file, unsigned int idx, SmurfTree* smurf)
{
  if(!smurf){
    // if called with smurf=0 this is taken as the headline
    file << "   "
	 << std::setw(10) << std::right << " [idx] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [run] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [event] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [nvtx] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [sumet] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [met] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [njets] " << std::setw(15) << std::right
	 << std::setw(10) << std::right << " [evtype] " << std::setw(15) << std::right
	 << " \n";
  }
  else{ 
    file << std::setw(10) << std::right << idx 
	 << std::setw(10) << std::right << smurf->run_ 
	 << std::setw(11) << std::right << smurf->event_ 
	 << std::setw( 9) << std::right << smurf->nvtx_ 
	 << std::setw(10) << std::right << smurf->sumet_ 
	 << std::setw(12) << std::right << smurf->met_ 
	 << std::setw( 7) << std::right << smurf->njets_ 
	 << std::setw( 9) << std::right << smurf->evtype_ 
	 << " \n";
  }
  return;
}

void 
TreeValidator::fileHeader(std::ofstream& file, std::string value, SmurfTree* smurf)
{
  // prepare run and event as varibales to determine min/max for
  std::vector<std::string> vars; vars.push_back(std::string("run"  )); vars.push_back(std::string("event"));
  std::map<std::string, std::vector<double> > refEnv = envelope(ref_ , vars), valEnv = envelope(smurf, vars);

  // set up header for logfiles
  file << " SmurfTree Validation \n"
       << " \n"
       << " " << std::setw(52) << std::right <<  "N(Events)"
       << " " << std::setw(10) << std::right <<  "minRun   "
       << " " << std::setw(10) << std::right <<  "maxRun   "
       << " " << std::setw(10) << std::right <<  "minEvent "
       << " " << std::setw(10) << std::right <<  "minEvent "
       << " \n"
       << " Tree Reference:   " << std::setw(20) << std::right << refName_
       << " " << std::setw(10) << std::right <<  ref_->tree_->GetEntriesFast()
       << " " << std::setw(10) << std::right <<  refEnv.find("run"  )->second.at(0) // min
       << " " << std::setw(10) << std::right <<  refEnv.find("run"  )->second.at(2) // max
       << " " << std::setw(10) << std::right <<  refEnv.find("event")->second.at(0)
       << " " << std::setw(10) << std::right <<  refEnv.find("event")->second.at(2)
       << " \n"
       << " Tree Test     :   " << std::setw(20) << std::right << value
       << " " << std::setw(10) << std::right <<  smurf->tree_->GetEntriesFast()
       << " " << std::setw(10) << std::right <<  valEnv.find("run"  )->second.at(0) // min
       << " " << std::setw(10) << std::right <<  valEnv.find("run"  )->second.at(2) // max
       << " " << std::setw(10) << std::right <<  valEnv.find("event")->second.at(0)
       << " " << std::setw(10) << std::right <<  valEnv.find("event")->second.at(2)
       << " \n"
       << " \n"
       << " Branches ["; 
  for(unsigned int idx=0; idx<branchList_.size(); ++idx){
    file << branchList_[idx]; if(idx<branchList_.size()-1) file << " ,";
  }
  file << "] \n";
  file << "  \n";
  if(file==categoryA_){
    file << " Events in "  << refName_ << " which are not present in " << value << ":"
	 << " \n";
  }
  else if(file==categoryB_){
    file << " Events in "  << value << " which are not present in " << refName_ << ":"
	 << " \n";
  }
  else{
    file << " Events which are present in both files, but missmatch: \n"
	 << " \n";
  }
  file << " \n";
  fileLine(file, 0, 0);
  file << " \n";
  return;
}

void
TreeValidator::initEventList()
{
  if(debug_) {info(std::string("Building up event list for reference SmurfTree...")); }
  unsigned int nevent = (int)ref_->tree_->GetEntriesFast();
  for(unsigned int ientry=0; ientry<nevent; ++ientry){
    ref_->tree_->GetEntry(ientry); eventList_[ref_->event_]=ientry;
  }
}

void 
TreeValidator::initBranchList()
{
  if(debug_) {info(std::string("Inizialising list of available TBranches in reference SmurfTree...")); }
  TObjArray* branches = ref_->tree_->GetListOfBranches();
  for(int idx=0; idx<branches->GetEntriesFast(); ++idx){
    TBranch* branch = dynamic_cast<TBranch*>(branches->At(idx));
    if(!branch){ warning("TObject in list of TBranches is not recognized as TBranch."); return; }
    else{
      branchList_.push_back(branch->GetName()); 
    }
  }
  return;
}  


bool 
TreeValidator::validateBranches(SmurfTree* smurf)
{
  bool validated=true;
  if(debug_) {info(std::string("Comparing list of TBranches in test SmurfTree with reference SmurfTree...")); }
  TObjArray* branches = smurf->tree_->GetListOfBranches();
  std::vector<std::string> branchList = branchList_;
  for(int idx=0; idx<branches->GetEntriesFast(); ++idx){
    TBranch* branch = dynamic_cast<TBranch*>(branches->At(idx));
    if(!branch){ warning("TObject in list of TBranches is not recognized as TBranch."); validated|=false; }
    else{
      std::vector<std::string>::iterator branchIt = std::find(branchList.begin(), branchList.end(), std::string(branch->GetName()));
      if(branchIt == branchList.end()){
	validated|=false;
	warning(std::string("Branch ").append(branch->GetName()).append(" not in list of TBranches in reference TTree."));
      }
      else{
	branchList.erase(branchIt);
      }
    }
  }
  if(!branchList.empty()){
    std::string issue("The following branches which are present in the referecne SmurfTree are missing in the test SmurfTree: \n");
    for(std::vector<std::string>::const_iterator branch=branchList.begin(); branch!=branchList.end(); ++branch){
      issue += std::string(" + ").append(*branch).append("\n");
    }
    warning(issue);
    validated|=false;
  }
  return validated;
}

bool 
TreeValidator::validateContent(SmurfTree* smurf){
  bool validated=true;
  if(debug_) {info(std::string("Validating content of SmurfTree...")); }
  unsigned int inBoth=0, catBErr=0, catCErr=0;
  unsigned int nevent = (int)smurf->tree_->GetEntries();
  for(unsigned int ientry=0; ientry<nevent; ++ientry){
    smurf->tree_->GetEntry(ientry);
    if(!isInEventList(smurf->event_)){
      // this is an event of category B; Note: the events in value are not sorted 
      // thus this will also not either be the case in the corresponding logfile
      fileLine(categoryB_, ++catBErr, smurf); validated|=false; continue;
    } ++inBoth; // event is present in both SmurfTrees

    ref_->tree_->GetEntry(getEntry(smurf->event_));
    std::vector<std::string> fails = ref_->Compare(smurf);
    if(!fails.empty()) {
      // up to now we only print the failing varibales here; the varibales in the 
      // logfile belong to a standard selection
      if(debug_){ warning(std::string("Found inconsistentcy in event"));
	for(std::vector<std::string>::const_iterator fail=fails.begin(); fail!=fails.end(); ++fail){
	  std::cout << "fail in " << *fail << " ref[" << ref_->Get(*fail) << "]"
		    << " test[" << smurf->Get(*fail) << "]" << std::endl; 
	}
      }
      // these are errors of category C print the variables of the failing events below 
      // each other
      fileLine(categoryC_, ++catCErr, smurf); fileLine(categoryC_, catCErr, smurf); categoryC_ << " \n";
      validated|=false;
    }
    // if the event was present in both SmurfTrees erase is from the list; those events
    // which are left over are errors of category B 
    eraseFromEventList(smurf->event_);
  }

  unsigned int catAErr=0;
  if(!eventList_.empty()){
    for(std::map<unsigned int, unsigned int>::const_iterator ievent=eventList_.begin(); ievent!=eventList_.end(); ++ievent){
      ref_->tree_->GetEntry(ievent->second);
      // events in ref_ are sorted in increasing order of 'event', so is the output in 
      // the logfile
      fileLine(categoryA_, ++catAErr, ref_);
      validated|=false;
    }
  }
  
  // prepare a tail to summarize the content validation 
  ostringstream tail;
  tail << "-----------------------------------------------------------------------------"     << "\n"
       << "     Events present in both SmurfTrees                              : " << inBoth  << "\n"
       << "     Events present in test but NOT in reference (category A error) : " << catAErr << "\n"
       << "     Events present in reference but NOT in test (category B error) : " << catBErr << "\n"
       << "\n"
       << "     Check event lists in corresponding output files.";
  info(tail.str());
  info("-----------------------------------------------------------------------------");
  return validated;
}

void validateTree(){
  TreeValidator validate(std::string("hww160_dk.root"));
  validate.validate(std::string("hww160_gc.root"));
return;
}
