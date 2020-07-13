#ifndef __PROCESSES_H_INCLUDED__   // if x.h hasn't been included yet...
#define __PROCESSES_H_INCLUDED__   //   #define this so the compiler knows it has been included

#ifdef __CLING__
#include <map>
#include <vector>
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TObject.h"

#include <string>
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

#include "ExRootTreeReader.h"
#include "ExRootTreeWriter.h"
#include "ExRootTreeBranch.h"
#include "ExRootResult.h"

#include "DelphesClasses.h"
#include "TClonesArray.h"

using namespace std;

class Processes {
public:
  Processes();
  std::map <string, double> ttB;
  std::map <string, double> tt;
  std::map <string, double> tB;
  std::map <string, double> Stops;
  std::map <string, std::map <string, double> > Catalog;
};


#endif 
