#include "Analysis.h"

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

//-------------------------------------------------------------------------------------------------------------------
int CohenJetPT(TClonesArray *br_Jet, const int njetreq = 2, const double ptJetCut = 1000., const double etaJetCut = 2.5){

    int nJetsEvent = br_Jet->GetEntriesFast();
    int jetsPassed = 0;

    if(nJetsEvent < njetreq) return 0;

    Jet* jet;
    for(int i=0; i < nJetsEvent; i++){
      jet = (Jet*) br_Jet->At(i);
      if(fabs(jet->Eta) > etaJetCut) continue;
      if(jet->PT > ptJetCut) jetsPassed += 1; 
      if(jetsPassed == njetreq) return 1;
    }
    return 0;
}

//-------------------------------------------------------------------------------------------------------------------
int CohenMuonInJet(TClonesArray *br_Jet, TClonesArray *br_EFlowMuon, const double ptMuInJetCut = 1000., const int nJetMuons = 1, const int nLeadJets = 2, const double jetR = 0.5, const double etaCut = 2.5){

    const double PI = 3.14159265359;
    int muonInJet = 0;
    int nMuonEntries = br_EFlowMuon->GetEntriesFast();
    int jetsLoopSize = nLeadJets;

    Muon *muon;
    Jet *jet;

    if(br_Jet->GetEntriesFast() < nLeadJets) jetsLoopSize = br_Jet->GetEntriesFast();

    for(int i=0; i < jetsLoopSize; i++) {
       jet = (Jet*) br_Jet->At(i);
       if (fabs(jet->Eta) > etaCut) continue;

       for (int j = 0; j < nMuonEntries; ++j) {
          muon = (Muon *) br_EFlowMuon->At(j);
	  if(muon->PT > ptMuInJetCut && fabs(muon->Eta) < etaCut) {
  	     double dphi = fabs(muon->Phi - jet->Phi);
	     if (dphi > PI) dphi = 2*PI - dphi;
	     double deta = fabs(muon->Eta - jet->Eta);
	     if (sqrt(dphi*dphi + deta*deta) < jetR) muonInJet += 1;
	     if (muonInJet == nJetMuons) return 1;
           }
	}
     }
     return 0;
}

//-------------------------------------------------------------------------------------------------------------------
int CohenNoIsoLepton(TClonesArray *br_Electron, TClonesArray *br_Muon, const double ptLepton = 35., const double ptIsolationRatio = 0.1, const int isoLTolerance = 0, const double etaCut = 2.5){

     int nElectronEntries = br_Electron->GetEntriesFast();
     int nMuonEntries = br_Muon->GetEntriesFast();
     int nIsoLeptons = 0;
     int counter = 0;
     Electron *electron;
     Muon *muon;

     for (int i = 0; i < nElectronEntries; ++i) {
     electron = (Electron *) br_Electron->At(i);
     if (electron->PT > ptLepton && fabs(electron->Eta) < etaCut) nIsoLeptons += 1 ;
     ++counter;
     }
     for (int i = 0; i < nMuonEntries; ++i) {
     muon = (Muon *) br_Muon->At(i);
     if (muon->PT > ptLepton && fabs(muon->Eta) < etaCut) nIsoLeptons += 1 ;
     ++counter;
     }

     if (nIsoLeptons > isoLTolerance) return 0;
     if (counter == nMuonEntries + nElectronEntries) return 1;
}

//-------------------------------------------------------------------------------------------------------------------
int CohenPhiIsoMET(TClonesArray *br_MissingET, TClonesArray *br_Jet, const double phiCut = 1.0, const double JetpTCut = 200., const double etaCut = 2.5){

    const double PI = 3.14159265359;
    double jetdPhi = 0.;
    int nJetsEvent = br_Jet->GetEntriesFast();

    MissingET *met = (MissingET*) br_MissingET->At(0);
    Jet *jet;

    for (int i = 0; i < nJetsEvent; ++i){
	jet = (Jet*) br_Jet->At(i);
	if (jet->PT < JetpTCut || fabs(jet->Eta) > etaCut) continue;
        jetdPhi = fabs(jet->Phi - met->Phi);
	if (jetdPhi > PI) jetdPhi = 2*PI - jetdPhi;
	if (jetdPhi < phiCut) return 0;
	if (i == nJetsEvent) return 1;
    }
}

//-------------------------------------------------------------------------------------------------------------------

int CohenMET(TClonesArray *br_MissingET, const double eMetCut = 3000.0){

    MissingET *met = (MissingET*) br_MissingET->At(0);
    if (met->MET >  eMetCut) return 1;
    return 0;
}

//-------------------------------------------------------------------------------------------------------------------

int CutScalarHT(TClonesArray *br_ScalarHT, const double scalarHTCut = 2000.0){

    ScalarHT *sht = (ScalarHT*) br_ScalarHT->At(0);
    if (sht->HT >  scalarHTCut) return 1;
    return 0;
}

//-------------------------------------------------------------------------------------------------------------------


int ZJetDistribution(TClonesArray *br_Jet, const double ZVarCut = 0.8, const double jetPtCut = 200., const int minJetsEvent = 4, const int leadingJets = 2) {

    int jetsLoopSize = minJetsEvent;
    double z = -1.;
    double x = 0. ,y = 0.;

    if(br_Jet->GetEntriesFast() < minJetsEvent) jetsLoopSize = br_Jet->GetEntriesFast();
    Jet *lastjet = (Jet*) br_Jet->At(jetsLoopSize-1);
    if(lastjet->PT < jetPtCut) return 0;

    for(int i=0; i < jetsLoopSize; i++){
	Jet *jet = (Jet*) br_Jet->At(i);
	if (jet->PT > jetPtCut && i < leadingJets) y += jet->PT;
	else if(jet->PT > jetPtCut) x += jet->PT;
    }
    z = 2.*x/y;
    if(z > ZVarCut) return 1;
    return 0;
}
	
