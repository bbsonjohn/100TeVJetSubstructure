#ifndef ANALYSISJET_H
#define ANALYSISJET_H

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

#include "external/fastjet/internal/base.hh"
#include "external/fastjet/internal/numconsts.hh"
#include "external/fastjet/PseudoJet.hh"
#include "external/fastjet/JetDefinition.hh"
#include "external/fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace std;

int JetPushBack(vector <fastjet::PseudoJet> fjInputs, TClonesArray *br_Particle);

int SingleJetMass(vector <fastjet::PseudoJet> sortedJets, const double mJetCut, const double ptJetCut, const double etaJetCut);

int DiJetResonance(vector <fastjet::PseudoJet> sortedJets, const double mJetCut, const double ptJetCut, const double etaJetCut);

#endif 

