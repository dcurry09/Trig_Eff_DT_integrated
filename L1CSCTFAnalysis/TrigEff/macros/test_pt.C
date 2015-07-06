#include<exception>
#include<TFile.h>
#include<TChain.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<TGraphErrors.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
#include <iostream>
#include <fstream>
#include <TGraph2D.h>

using namespace std;

// start of main function
void test_pt() {

  gROOT->Clear();
  gStyle->SetOptStat(111111);

  // open the root file
  TFile* file = TFile::Open("dataset/MinBias_2012AB.root");

  TTree *csctfTTree;
  TTree *recoMuons;

  recoMuons  = (TTree*) file -> Get("recoMuons");
  csctfTTree = (TTree*) file -> Get("csctfTTree");

  // Access needed variables
  int Run, Event, muonSize;
  vector<float>* ptReco = new vector<float>();

  recoMuons->SetBranchAddress("Run"  , &Run  );
  recoMuons->SetBranchAddress("Event", &Event);
  recoMuons->SetBranchAddress("muonSize", &muonSize );
  recoMuons->SetBranchAddress("gmrPt" , &ptReco );

  // Initialize histograms
  TH1F* hPt = new TH1F("hPt", "", 100, 0, 20);
  
  // Loop over the events
  for (int iEvt=0; iEvt < recoMuons->GetEntries(); iEvt++) {
    
    recoMuons ->GetEntry(iEvt);
    
    if ( ( iEvt % 10000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);
    
    // Loop over recponstructed muons
    for (int iReco=0; iReco < muonSize; iReco++) {
      
      // get the pT for each muon and fill histogram
      cout << "pT is = " << ptReco->at(iReco) << endl;
      hPt ->Fill ( ptReco->at(iReco) ); 
      
    }// end reco loop
  } // end event loop

  // Print the histogram


} //end main function
