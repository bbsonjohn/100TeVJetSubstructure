// main71.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Richard Corke.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

/*
 * Simple example of fastjet analysis. Roughly follows analysis of:
 * T. Aaltonen et al. [CDF Collaboration],
 * Measurement of the cross section for stop quark measurement sqrt(s)=100$ TeV
 * 
 * arXiv:1406.4512v2 [hep-ex]
 *
 * Cuts:
 *   

 *   mT(W)        > 20GeV
 */
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

// This is the minimal interface needed to access FastJet.
// A more sophisticated interface is demonstrated in main72.cc.
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;



int main(int argc, char* argv[]) {

  const double PI = 3.1415926;
  int nEvent = 10000;
  int nAbort = 10;

  // Cuts　1
  const int cut1Enable = 0;
  const double jetR = 0.5;
  const double ptJetCut = 1000;
  const double etaJetCut = 2.5;
  // Cuts 2
  const int cut2Enable = 1;
  const double ptMuInJetCut = 200.;
  // Cuts 3
  const int cut3Enable = 0;
  const double ptLeptonCut = 35.;
  const double ptIsolationRatio = 0.1;
  const double etaLeptonCut = 2.5;
  // Cuts 4
  const int cut4Enable = 0;
  const double metJetpTCut = 200.;
  const double phiMetJetCut = 1.0;
  // Cuts 5
  const int cut5Enable = 0;
  const double eMetCut = 3000.;

  // Root Settings
  TApplication theApp("hist", &argc, argv);

  TFile* outFile = new TFile("boostedStop.root", "RECREATE");

  TH1F *rawMissingPt = new TH1F("missingPtUncut","missing Pt (Uncut)", 20, .0, 10000.0);
  TH1F *missingPt = new TH1F("missingPt","missing Pt", 20, .0, 10000.0);
  TH1F *nJets = new TH1F("nJets","n jets", 30, 0., 30);
  TH1F *jet0Pt = new TH1F("jet0Pt","Jet 0 Pt", 50, 500, 8000);
  TH1F *jet1Pt = new TH1F("jet1Pt","Jet 1 Pt", 50, 500, 8000);
  TH1F *jetMetSeparation = new TH1F("jetMetSeparation","Jet MET Separation", 50, -0.1, PI+0.1);
  TH1F *rawJetMetSeparation = new TH1F("jetMetSeparationUncut","Jet MET Separation (Uncut)", 50, -0.1, PI+0.1);
  TH1F *muInJet0Efficiency = new TH1F("muInJet0Efficiency","Efficiency Finding Muon in Leading Jet", 9, 0, 8000);
  TH1F *METEfficiency = new TH1F("METEfficiency","MET after Cut", 20, 0, 10000);
  TH1F *phiEfficiency = new TH1F("phiEfficiency","Minimum MET-Jet delta phi after Cut", 25, 0, PI); 
  TH1F *njetscut = new TH1F("njetscut","NJet < 2", 3, -1, 1);
  TH1F *nmucut = new TH1F("nmucut","Nmu < 1", 3, -1, 1); 
  TH1F *isolatedlveto = new TH1F("isolatedlveto","isolated lepton", 3, -1, 1); 
  TH1F *deltaphicut = new TH1F("deltaphicut","Delta Phi MET-Jet < 1.0", 3, -1, 1); 
  TH1F *metcut = new TH1F("metcut","MET < MET Cut", 3, -1, 1); 
  TH1F *totalEvents = new TH1F("totalEvents","Total Number of Events", 3, -1, 1); 
  TH1F *cutEvents = new TH1F("cutEvents","Number of Events after Cut", 3, -1, 1);
  TH1F *jetpTCut = new TH1F("jetpTCut","Leading or sub-leading JetpT < Cut", 3, -1, 1);


  // Generator
  Pythia pythia;

  pythia.readString("PhaseSpace:pTHatMin = 25.0");

  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("HadronLevel:Hadronize = on");


  // Initialisation, p pbar @ 1.96 TeV
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = boostedStop4000GeV.lhe");
  pythia.init();



  // Fastjet analysis - select algorithm and parameters
  double Rparam = jetR;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;


  bool firstEvent = true;
  int iAbort = 0;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    totalEvents->Fill(0);

    double mpx = 0.;
    double mpy = 0.;
    double mpz = 0.;
    double mee = 0.;
    double mpT = 0.;

    fjInputs.resize(0);

    // Loop over event record to decide what to pass to FastJet
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Final state only
      if (!pythia.event[i].isFinal())        continue;
      if (pythia.event[i].pT() < 0.0)        continue;

          if( !pythia.event[i].isVisible()) {
              mpx += pythia.event[i].px();
              mpy += pythia.event[i].py();
              mpz += pythia.event[i].pz();
              mee += pythia.event[i].e();
          }
 
      // No neutrinos
      if ( pythia.event[i].isVisible()) {
           fjInputs.push_back( fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()) );
      }
        // Only |eta| < 3.6
      if (fabs(pythia.event[i].eta()) > 5.0) continue;


      // Store as input to Fastjet
    }

    Vec4 met_vec(mpx,mpy,mpz,mee);

     mpT = sqrt(mpx*mpx+mpy*mpy);
     rawMissingPt->Fill(mpT) ;

    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // For the first event, print the FastJet details
    if (firstEvent) {
      pythia.process.list();
      cout << endl << "-----------------------------------" << endl << endl;
      cout << "Ran " << jetDef->description() << endl;
      cout << "Strategy adopted by FastJet was "
           << clustSeq.strategy_string() << endl << endl;
      firstEvent = false;
    }

    //BEGIN: check cut criteria 1 
    inclusiveJets = clustSeq.inclusive_jets(25.0);
    sortedJets    = sorted_by_pt(inclusiveJets);

    int jet0Pass = -1;
    int jet1Pass = -1;

    if (cut1Enable == 0) {
	jet0Pass = 0;
	jet1Pass = 1;
	}
    else {
	if (sortedJets.size() < 2) {
	continue;
	}
	for (int i = 0; i < sortedJets.size(); ++i) { 
	    if (sortedJets[i].perp() > ptJetCut && abs(sortedJets[i].eta()) < etaJetCut) {
	    jet0Pass = i;
	    break;
	    }
	}
	for (int i = jet0Pass+1; i < sortedJets.size(); ++i) { 
	    if (sortedJets[i].perp() > ptJetCut && abs(sortedJets[i].eta()) < etaJetCut) {
	    jet1Pass = i;
	    break;
	    }
	}
	if (jet0Pass >= 0 && jet1Pass >= 0) {
	   njetscut->Fill(0);
	   jetpTCut->Fill(0);
	   continue;
	}
    }


    jet0Pt->Fill(sortedJets[jet0Pass].perp()); 
    jet1Pt->Fill(sortedJets[jet1Pass].perp());
    //END: check cut criteria 1 

    bool muonInJet = false;
    bool isolatedl = false;

    for (int i = 0; i < pythia.event.size(); ++i) {

        if (!pythia.event[i].isFinal())        continue;
        if (pythia.event[i].pT() < 0.0)        continue;

        //BEGIN: check cut criteria 2 
	if (abs(pythia.event[i].id()) == 13 && pythia.event[i].pT() > ptMuInJetCut) {
	   double phiDist = abs(sortedJets[jet0Pass].phi() - pythia.event[i].phi());
	   if (phiDist > PI) phiDist = 2*PI - phiDist;
	   double etaDist = abs(sortedJets[jet0Pass].eta() - pythia.event[i].eta());
	   if (sqrt(phiDist*phiDist + etaDist*etaDist) < jetR) {
		muonInJet = true;
		muInJet0Efficiency->Fill(sortedJets[jet0Pass].perp());
		}

	   phiDist = abs(sortedJets[jet1Pass].phi() - pythia.event[i].phi());
	   if (phiDist > PI) phiDist = 2*PI - phiDist;
	   etaDist = abs(sortedJets[jet1Pass].eta() - pythia.event[i].eta());
	   if (sqrt(phiDist*phiDist + etaDist*etaDist) < jetR) muonInJet = true;
	}
        //END: check cut criteria 2 

        //BEGIN: check cut criteria 3 
	if (abs(pythia.event[i].id()) == 13 || abs(pythia.event[i].id()) == 11) {
	   if (pythia.event[i].pT() > ptLeptonCut && pythia.event[i].eta() < etaLeptonCut) {
		double isolationpT = 0.0;
		for (int j = 0; j < pythia.event.size(); ++j) {
		   if (j == i) continue;
      		   if (!pythia.event[i].isFinal())        continue;
      		   if (pythia.event[i].pT() < 0.0)        continue;

		   double phiDist = abs(pythia.event[j].phi() - pythia.event[i].phi());
		   if (phiDist > PI) phiDist = 2*PI - phiDist;
		   double etaDist = abs(pythia.event[j].eta() - pythia.event[i].eta());
		   if (sqrt(phiDist*phiDist + etaDist*etaDist) < jetR)  isolationpT += pythia.event[j].pT();
		   }

		if (isolationpT < ptIsolationRatio*pythia.event[i].pT()) isolatedl = true; 
	    }
	}
        //END: check cut criteria 3 
  
    }

    if (cut2Enable == 1) {  
	if (muonInJet == false) {
	    continue;	
	}
	else {
	    nmucut->Fill(0);	    
	}
    }

    if (cut3Enable == 1) {
	if (isolatedl == true) {
	    continue;
	    }
	else {
	    isolatedlveto->Fill(0);
	    }
    }

    //BEGIN: check cut criteria 4 

    bool isolatedMet = true;
    double jetMetPhi = 0.;
    double maxPhiDist = 10.; //any number > 2*PI

   for (int i = 0; i < sortedJets.size(); ++i){
	if (sortedJets[i].perp() > metJetpTCut && abs(sortedJets[i].eta()) < etaJetCut){
		Vec4 jet_vec(sortedJets[i].px(),sortedJets[i].py(),sortedJets[i].pz(),sortedJets[i].E());
		jetMetPhi = phi(jet_vec,met_vec);
		if (jetMetPhi > PI) jetMetPhi = 2*PI - jetMetPhi;
		rawJetMetSeparation->Fill(jetMetPhi);
		if (jetMetPhi < maxPhiDist) maxPhiDist = jetMetPhi;
		if (jetMetPhi < phiMetJetCut) isolatedMet = false;
	}
    	
   }
    if (cut4Enable == 1) {
	if (isolatedMet == false) {
	    continue;	
	}
	else {
	    deltaphicut->Fill(0);
	}
    }

    jetMetSeparation->Fill(jetMetPhi);

    //END: check cut criteria 4

    //BEGIN: check cut criteria 5 

     missingPt->Fill(mpT) ;

    if (cut5Enable == 1) {
	if (mpT < eMetCut) {
	    continue;
	}
	else {
	    metcut->Fill(0);
	}
    }

    //END: check cut criteria 5 

   phiEfficiency->Fill(maxPhiDist);
   METEfficiency->Fill(mpT);
   nJets->Fill(sortedJets.size());
   cutEvents->Fill(0);

  // Output histograms
  }


  std::cout << "\nDouble click on the histogram window to quit.\n";


  // Save histogram on file and close file.


  nJets->Write();
  missingPt->Write();
  rawMissingPt->Write();
  jetMetSeparation->Write();
  rawJetMetSeparation->Write();
  jet0Pt->Write();
  jet1Pt->Write();
  muInJet0Efficiency->Write();
  phiEfficiency->Write();
  METEfficiency->Write();
  njetscut->Write();
  nmucut->Write();
  isolatedlveto->Write();
  deltaphicut->Write();
  metcut->Write();
  jetpTCut->Write();
  totalEvents->Write();
  cutEvents->Write();

  jetMetSeparation->Draw();
  gPad->WaitPrimitive();


  delete outFile;

  // Done.

  return 0;
}
