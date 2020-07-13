#include "AnalysisJet.h"

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

int JetPushBack(vector <fastjet::PseudoJet> fjInputs, TClonesArray *br_Particle) {
	for(int i=0; i<br_Particle->GetEntriesFast(); i++) {
		GenParticle *particle = (GenParticle *) br_Particle->At(i);
		TLorentzVector p4 = particle->P4();
        fjInputs.push_back( fastjet::PseudoJet(p4.X(), p4.Y(), p4.Z(), p4.E()) );
    }
	return 0;
}

//-------------------------------------------------------------------------------------------------------------------

int SingleJetMass(vector <fastjet::PseudoJet> sortedJets, const double mJetCut, const double ptJetCut = 200., const double etaJetCut = 2.5){
	for (int i = 0; i < sortedJets.size(); ++i) { 
		if (sortedJets[i].perp() > ptJetCut && abs(sortedJets[i].eta()) < etaJetCut) {
		if (sortedJets[i].m() > mJetCut) return 1;
	    }
	}

    return 0;
}

//-------------------------------------------------------------------------------------------------------------------

int DiJetResonance(vector <fastjet::PseudoJet> sortedJets, const double mJetCut, const double ptJetCut = 200., const double etaJetCut = 2.5){

	for (int i = 0; i < sortedJets.size()-1; ++i) { 
		if (sortedJets[i].perp() > ptJetCut && abs(sortedJets[i].eta()) < etaJetCut) {
			if (sortedJets[i+1].perp() > ptJetCut && abs(sortedJets[i+1].eta()) < etaJetCut) {
				fastjet::PseudoJet jet_vec(sortedJets[i].px()+sortedJets[i+1].px(),sortedJets[i].py()+sortedJets[i+1].py(),sortedJets[i].pz()+sortedJets[i+1].pz(),sortedJets[i].E()+sortedJets[i+1].E());
				if (abs(jet_vec.m()) < mJetCut) return 1;
			}
	    }
	}

    return 0;
}
//-------------------------------------------------------------------------------------------------------------------
