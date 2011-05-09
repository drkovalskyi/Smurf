#include <string>

enum process {
  proc_WW,
  proc_Wpj,
  proc_Wmj,
  proc_Wpg,
  proc_Wmg,
  proc_ZZ,
  proc_HWW120,
  proc_HWW130,
  proc_HWW140,
  proc_HWW150,
  proc_HWW160,
  proc_HWW170,
  proc_HWW180,
  proc_HWW190,
  proc_HWW200,
  proc_HWW210,
  proc_HWW220,
  proc_HWW230,
  proc_HWW250,
  proc_HWW300,
  kNProc
};

// this is needed to interface with the smurf trees
std::string name(int  dstype){
  switch (dstype){
  case proc_WW:     return "ww";
  case proc_Wpj:    return "wjets";
  case proc_Wmj:    return "wjets";
  case proc_Wpg:    return "wgamma";
  case proc_Wmg:    return "wgamma";
  case proc_ZZ:     return "zz";
  case proc_HWW120: return "hww120";
  case proc_HWW130: return "hww130";
  case proc_HWW140: return "hww140";
  case proc_HWW150: return "hww150";
  case proc_HWW160: return "hww160";
  case proc_HWW170: return "hww170";
  case proc_HWW180: return "hww180";
  case proc_HWW190: return "hww190";
  case proc_HWW200: return "hww200";
  case proc_HWW210: return "hww210";
  case proc_HWW220: return "hww220";
  case proc_HWW230: return "hww230";
  case proc_HWW250: return "hww250";
  case proc_HWW300: return "hww300";
  default:     return "uknown";
  }
};
