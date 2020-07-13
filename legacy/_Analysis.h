#ifndef ANALYSIS_H
#define ANALYSIS_H
 
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

int CohenJetPT(TClonesArray *br_Jet, const int, const double ptJetCut, const double etaJetCut);

int CohenMuonInJet(TClonesArray *br_Jet, TClonesArray *br_EFlowMuon, const double ptMuInJetCut, const int nJetMuons, const int nLeadJets, 
const double jetR, const double etaCut);

int CohenNoIsoLepton(TClonesArray *br_Electron, TClonesArray *br_Muon, const double ptLepton, const double ptIsolationRatio, const int isoLTolerance, const double etaCut);

int CohenPhiIsoMET(TClonesArray *br_MissingET, TClonesArray *br_Jet, const double phiCut, const double JetpTCut, const double etaCut);

int CohenMET(TClonesArray *br_MissingET, const double eMetCut);

int CutScalarHT(TClonesArray *br_ScalarHT, const double scalarHTCut);

int ZJetDistribution(TClonesArray *br_Jet, const double ZVarCut, const double jetPtCut, const int minJetsEvent, const int leadingJets);

#endif
