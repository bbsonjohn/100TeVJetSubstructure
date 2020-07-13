/*
Cut Functions: Return 1 when cut is passed. Return 0 when cut is failed.
*/

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

#include "external/fastjet/internal/base.hh"
#include "external/fastjet/internal/numconsts.hh"
#include "external/fastjet/PseudoJet.hh"
#include "external/fastjet/JetDefinition.hh"
#include "external/fastjet/ClusterSequence.hh"

#include "examples/Analysis.h"
#include "examples/processes.h"

using namespace fastjet;

using namespace std;


int main()
{

  TChain chain("Delphes");

  Processes smbkg;

  const int file_error_tolerence = 2;
  char inputFile[150];
  char outputFile[150];

  const double luminosity = 2000.;
  double weight = 0.;


  char sm_process[30] = "ttB";
  int mgluino = 2000;

  for(int msquark = 3000; msquark < 9000; msquark = msquark + 2500){
  sprintf(inputFile, "/media/john/EC7A174B7A1711C6/Linux/Stops100TeV/%d_%d/001.root", msquark, mgluino);
  sprintf(outputFile, "OutputGridFiles/JetMass/point_%d_mstop_%d_mgluino_analysis.root", msquark, mgluino);

  TFile* outFile = new TFile(outputFile, "RECREATE");

  TH1F *scalarHT = new TH1F("scalarHT","Scalar Sum HT", 100, .0, 50000.);
  TH1F *missingET = new TH1F("missingET","MET", 100, .0, 10000.);
  TH1F *numberOfJet = new TH1F("numberOfJet","Number Of Jet", 20, .0, 20.);
  TH1F *massjet[6];
  for(int i = 0; i < 6; i++){
	char name[20];
	char title[100];
	sprintf(name, "massjet%i", i+1);
	sprintf(title, "Jet %i Mass", i+1);
	massjet[i] = new TH1F(name,title, 100, .0, 2000.0);
  }
/*
  double Rparam = 0.5;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;


    fjInputs.resize(0);
*/

  for(std::map<string,double>::iterator it = smbkg.Stops.begin(); it != smbkg.Stops.end(); it++) {
	string msquark_ss;
	msquark_ss.append(it->first);
	int msquark_s = atoi(msquark_ss.c_str());
	if(msquark_s == msquark){
      weight = 1000.*luminosity*(it->second);
      break;
	}
  }
  	  chain.Add(inputFile);

	  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	  Long64_t numberOfEntries = treeReader->GetEntries();

	  TClonesArray *branchEvent  	  = treeReader->UseBranch("Event");
  	  TClonesArray *branchJet 	  = treeReader->UseBranch("Jet");
  	  TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
  	  TClonesArray *branchMuon	  = treeReader->UseBranch("Muon");
  	  TClonesArray *branchElectron  = treeReader->UseBranch("Electron");
  	  TClonesArray *branchParticle  = treeReader->UseBranch("Particle");
  	  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  	  TClonesArray *branchScalarHT  = treeReader->UseBranch("ScalarHT");
	  double n_entries = numberOfEntries;
	  if(n_entries < 1) {
		cout << inputFile << " contains insufficient entries!" << endl;
		continue;
	  }
      weight = weight/n_entries;
      cout << "Reading: " << inputFile << " with bin weight " << weight << endl; 

  	  const double inclusive_jet_pt = 0.0;
  	  const double PI = 3.14159265359;
  	  const double jetR = 0.5;

	  const double scalarHTcut = 3500;
	  const double missingETcut = 500;


  // Loop over all events
  		for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    		treeReader->ReadEntry(entry);

            JetDefinition *definition;
            definition = new JetDefinition(cambridge_algorithm, 0.5);
			fastjet::ClusterSequence::print_banner();

			if(CutScalarHT(branchScalarHT, scalarHTcut) < 1) continue;
			if(CohenMET(branchMissingET, missingETcut) < 1) continue;
			if(branchJet->GetEntriesFast() < 6) continue;

			MissingET *met = (MissingET*) branchMissingET->At(0);
			ScalarHT *sht = (ScalarHT*) branchScalarHT->At(0);

			numberOfJet->Fill(branchJet->GetEntriesFast(),weight);
			scalarHT->Fill(sht->HT,weight);
			
			missingET->Fill(met->MET,weight);
			

			for(int njet = 0; njet < 6; ++njet){
				Jet* jet = (Jet*) branchJet->At(njet);

				massjet[njet]->Fill(jet->Mass,weight);
			}		
   		}


  outFile->Write();
  outFile->Close();

  cout << endl << "Output File - " << outputFile << endl;

  delete outFile;
  }

}


//-------------------------------------------------------------------------------------------------------------
bool check_File (const char * file) {
    ifstream f(file);
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

    



