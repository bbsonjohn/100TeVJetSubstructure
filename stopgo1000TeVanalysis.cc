// main71.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Richard Corke.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

/*
 * Simple example of fastjet analysis. Roughly follows analysis of:
 * T. Aaltonen et al. [CDF Collaboration],
 * Measurement of the cross section for W-boson production in association
 * with jets in ppbar collisions at sqrt(s)=1.96$ TeV
 * Phys. Rev. D 77 (2008) 011108
 * arXiv:0711.4044 [hep-ex]
 *
 * Cuts:
 *   ET(elec)     > 20GeV
 *   |eta(elec)|  < 1.1
 *   ET(missing)  > 30GeV
 *   ET(jet)      > 20GeV
 *   |eta(jet)|   < 2.0
 *   deltaR(elec, jet) > 0.52
 * Not used:
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

  const double Rparam = 0.4;
  const double jetDefEtaCut = 3.6;
  const double jetDefPtCut = 0.0;

  const double jetPtCut = 25.0;
  const int nJetCut = 2;


  const int nAbort = 10;
  const int nEvent = 100000;

  bool firstEvent = true;

  
  // Root Settings
  TApplication theApp("hist", &argc, argv);

  TFile* outFile = new TFile("stopgo1000TeVanalysis.root", "RECREATE");

  TH1F *missingPt = new TH1F("missingPt","missing Pt", 100, .0, 8000.0);
  TH1F *nJets = new TH1F("nJets","n jets", 35, -0.5, 34.5);
  TH1F *jet0Pt = new TH1F("jet0Pt","Jet 0 Pt", 100, -50, 5000);
  TH1F *jet1Pt = new TH1F("jet1Pt","Jet 1 Pt", 100, -50, 5000);
  TH1F *jet2Pt = new TH1F("jet2Pt","Jet 2 Pt", 100, -50, 5000);
  TH1F *jet3Pt = new TH1F("jet3Pt","Jet 3 Pt", 100, -50, 5000);
  TH1F *jet4Pt = new TH1F("jet4Pt","Jet 4 Pt", 100, -50, 5000);
  TH1F *jet5Pt = new TH1F("jet5Pt","Jet 5 Pt", 100, -10, 5000);
  TH1F *jet6Pt = new TH1F("jet6Pt","Jet 6 Pt", 100, -10, 1000);
  TH1F *jet7Pt = new TH1F("jet7Pt","Jet 7 Pt", 100, -10, 1000);
  TH1F *jet8Pt = new TH1F("jet8Pt","Jet 8 Pt", 100, -10, 1000);
  TH1F *jet9Pt = new TH1F("jet9Pt","Jet 9 Pt", 100, -10, 1000);
  TH1F *jetMetSeparation = new TH1F("jetMetSeparation","Jet MET Separation", 80, -0.01 , 3.19);
  TH1F *jetSeparation = new TH1F("jetSeparation","Jet Separation", 100, -0.1 , 7.5);
  TH1F *jetPhiSeparation = new TH1F("jetPhiSeparation","Jet Phi Separation", 80, -0.01 , 3.19);
  TH1F *nearest2JetPhi = new TH1F("nearest2JetPhi","Nearest 2 Jet Phi Separation", 80, -0.01 , 1.59);
  TH1F *nearest3JetPhi = new TH1F("nearest3JetPhi","Nearest 3 Jet Phi Separation", 80, -0.01 , 1.59);
  TH1F *nearest4JetPhi = new TH1F("nearest4JetPhi","Nearest 4 Jet Phi Separation", 80, -0.01 , 1.59);
  TH1F *nearest5JetPhi = new TH1F("nearest5JetPhi","Nearest 5 Jet Phi Separation", 80, -0.01 , 1.59);
  TH1F *nearest6JetPhi = new TH1F("nearest6JetPhi","Nearest 6 Jet Phi Separation", 80, -0.01 , 1.59);
  TH1F *jet0MetSeparation = new TH1F("jet0MetSeparation","Jet 0 Phi Met Separation", 80, -0.01 , 3.19);
  TH1F *jet1MetSeparation = new TH1F("jet1MetSeparation","Jet 1 Phi Met Separation", 80, -0.01 , 3.19);
  TH1F *jet2MetSeparation = new TH1F("jet2MetSeparation","Jet 2 Phi Met Separation", 80, -0.01 , 3.19);
  TH1F *jet3MetSeparation = new TH1F("jet3MetSeparation","Jet 3 Phi Met Separation", 80, -0.01 , 3.19);
  TH1F *metPhiDistribution = new TH1F("metPhiDistribution","MET Phi Distribution", 80, -0.01 , 6.39);

  // Generator
  Pythia pythia;

  pythia.readString("PhaseSpace:pTHatMin = 25.0");

  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("HadronLevel:Hadronize = on");

  // SLHA file
  pythia.readString("SLHA:readFrom = 2");
  pythia.readString("SLHA:file = godecay.slha");
  pythia.readString("SLHA:useDecayTable = true");

  // Initialisation, p pbar @ 1.96 TeV
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = stopGo1000TeVAnalysis_10stop_2Go_1N.lhe");
  pythia.init();





  // Fastjet analysis - select algorithm and parameters

  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;

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

    double mass = 0.0;
    double mpx = 0;
    double mpy = 0;
    double mpz = 0;
    double jetMetPhi = 0;

    fjInputs.resize(0);

    // Loop over event record to decide what to pass to FastJet
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Final state only
      if (!pythia.event[i].isFinal())        continue;

          if( !pythia.event[i].isVisible() || pythia.event[i].id() == 1000022 ) {
              mpx += pythia.event[i].px();
              mpy += pythia.event[i].py();
	      mpz += pythia.event[i].pz();
          }
 
      if ( pythia.event[i].isVisible() && pythia.event[i].id() != 1000022) {

           fjInputs.push_back( fastjet::PseudoJet( pythia.event[i].px(),
             pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e() ) );
      }

      if (fabs(pythia.event[i].eta()) > jetDefEtaCut) continue;

      if (pythia.event[i].pT() < jetDefPtCut)        continue;


      // Store as input to Fastjet
    }
     double mpt = sqrt(mpx*mpx + mpy*mpy);
     missingPt->Fill( mpt ) ;
     double metPhi = acos(mpx/mpt) ;
     if (mpy < 0) metPhi = 2*PI - metPhi;
     
     metPhiDistribution->Fill(metPhi);
     




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

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveJets = clustSeq.inclusive_jets(jetPtCut);
    sortedJets    = sorted_by_pt(inclusiveJets);
    if (sortedJets.size() < nJetCut) continue;

// Jet Momenta 

	   jet0Pt->Fill(sqrt(sortedJets[0].px()*sortedJets[0].px()+sortedJets[0].py()*sortedJets[0].py())) ; 
	   jet1Pt->Fill(sqrt(sortedJets[1].px()*sortedJets[1].px()+sortedJets[1].py()*sortedJets[1].py()));
	   if (sortedJets.size() > 2) jet2Pt->Fill(sqrt(sortedJets[2].px()*sortedJets[2].px()+sortedJets[2].py()*sortedJets[2].py()));
	   if (sortedJets.size() > 3) jet3Pt->Fill(sqrt(sortedJets[3].px()*sortedJets[3].px()+sortedJets[3].py()*sortedJets[3].py()));
	   if (sortedJets.size() > 4) jet4Pt->Fill(sqrt(sortedJets[4].px()*sortedJets[4].px()+sortedJets[4].py()*sortedJets[4].py()));
	   if (sortedJets.size() > 5) jet5Pt->Fill(sqrt(sortedJets[5].px()*sortedJets[5].px()+sortedJets[5].py()*sortedJets[5].py()));
	   if (sortedJets.size() > 6) jet6Pt->Fill(sqrt(sortedJets[6].px()*sortedJets[6].px()+sortedJets[6].py()*sortedJets[6].py()));
	   if (sortedJets.size() > 7) jet7Pt->Fill(sqrt(sortedJets[7].px()*sortedJets[7].px()+sortedJets[7].py()*sortedJets[7].py()));
	   if (sortedJets.size() > 8) jet8Pt->Fill(sqrt(sortedJets[8].px()*sortedJets[8].px()+sortedJets[8].py()*sortedJets[8].py()));
	   if (sortedJets.size() > 9) jet9Pt->Fill(sqrt(sortedJets[9].px()*sortedJets[9].px()+sortedJets[9].py()*sortedJets[9].py()));

   for (int j = 0; j < 4; ++j){


      for (int k = 0; k < j; ++k){
//         Vec4 p1(sortedJets[k].px(), sortedJets[k].py(), sortedJets[k].pz(), sortedJets[k].e());
	   double phiDist = abs(sortedJets[k].phi()-sortedJets[j].phi());
	   double etaDist = abs(sortedJets[k].eta()-sortedJets[j].eta());

 	   if (phiDist > PI) phiDist = 2*PI - phiDist;
	   jetPhiSeparation->Fill(phiDist);
	   jetSeparation->Fill(sqrt(phiDist*phiDist + etaDist*etaDist));
   }

    //average distance between Jets and MET

    jetMetPhi = abs(sortedJets[j].phi() - metPhi);
    if (jetMetPhi > PI) jetMetPhi = 2*PI - jetMetPhi;
    jetMetSeparation->Fill(jetMetPhi);

  }

// calculate difference in phi between MET and Jet

    double jetMetPhi1 = abs(sortedJets[0].phi() - metPhi);
    if (jetMetPhi1 > PI) jetMetPhi1 = 2*PI - jetMetPhi1;
    jet0MetSeparation->Fill(jetMetPhi1);

    double jetMetPhi2 = abs(sortedJets[1].phi() - metPhi);
    if (jetMetPhi2 > PI) jetMetPhi2 = 2*PI - jetMetPhi2;
    jet1MetSeparation->Fill(jetMetPhi2);

    double jetMetPhi3 = abs(sortedJets[2].phi() - metPhi);
    if (jetMetPhi3 > PI) jetMetPhi3 = 2*PI - jetMetPhi3;
    jet2MetSeparation->Fill(jetMetPhi3);

    double jetMetPhi4 = abs(sortedJets[3].phi() - metPhi);
    if (jetMetPhi4 > PI) jetMetPhi4 = 2*PI - jetMetPhi4;
    jet3MetSeparation->Fill(jetMetPhi4);

    double jetMetPhi5 = abs(sortedJets[4].phi() - metPhi);
    if (jetMetPhi5 > PI) jetMetPhi5 = 2*PI - jetMetPhi5;
    jet3MetSeparation->Fill(jetMetPhi5);

    double jetMetPhi6 = abs(sortedJets[5].phi() - metPhi);
    if (jetMetPhi6 > PI) jetMetPhi6 = 2*PI - jetMetPhi6;
    jet3MetSeparation->Fill(jetMetPhi6);

//Find the closest jet distance among nth highest momentum jets

    if (jetMetPhi1 > PI/2) jetMetPhi1 = PI - jetMetPhi1;

    if (jetMetPhi2 > PI/2) jetMetPhi2 = PI - jetMetPhi2;
    if (jetMetPhi1 < jetMetPhi2) jetMetPhi2 = jetMetPhi1;
    nearest2JetPhi->Fill(jetMetPhi2);    

    if (jetMetPhi3 > PI/2) jetMetPhi3 = PI - jetMetPhi3;
    if (jetMetPhi2 < jetMetPhi3) jetMetPhi3 = jetMetPhi2;
    nearest3JetPhi->Fill(jetMetPhi3); 
 
    if (jetMetPhi4 > PI/2) jetMetPhi4 = PI - jetMetPhi4;
    if (jetMetPhi3 < jetMetPhi4) jetMetPhi4 = jetMetPhi3;
    nearest4JetPhi->Fill(jetMetPhi4);  

    if (jetMetPhi5 > PI/2) jetMetPhi5 = PI - jetMetPhi5;
    if (jetMetPhi4 < jetMetPhi5) jetMetPhi5 = jetMetPhi4;
    nearest5JetPhi->Fill(jetMetPhi5);  

    if (jetMetPhi6 > PI/2) jetMetPhi6 = PI - jetMetPhi6;
    if (jetMetPhi5 < jetMetPhi6) jetMetPhi6 = jetMetPhi5;
    nearest6JetPhi->Fill(jetMetPhi6);  




   nJets->Fill(sortedJets.size());
  // Output histograms
  }

 std::cout << "\nDouble click on the histogram window to quit.\n";


  // Save histogram on file and close file.


  nJets->Write();
  missingPt->Write();
  jetSeparation->Write();
  jetPhiSeparation->Write();
  jetMetSeparation->Write();
  jet0MetSeparation->Write();
  jet1MetSeparation->Write();
  jet2MetSeparation->Write();
  jet3MetSeparation->Write();
  metPhiDistribution->Write();
  nearest2JetPhi->Write();
  nearest3JetPhi->Write();
  nearest4JetPhi->Write();
  nearest5JetPhi->Write();
  nearest6JetPhi->Write();
  jet0Pt->Write();
  jet1Pt->Write();
  jet2Pt->Write();
  jet3Pt->Write();
  jet4Pt->Write();
  jet5Pt->Write();
  jet6Pt->Write();
  jet7Pt->Write();
  jet8Pt->Write();
  jet9Pt->Write();


  jetMetSeparation->Draw();
  gPad->WaitPrimitive();


  delete outFile;

  // Done.

  return 0;
}
