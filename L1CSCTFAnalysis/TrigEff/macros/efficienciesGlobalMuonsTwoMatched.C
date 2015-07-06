#include<exception>
#include<vector>
#include<iostream>

#include<TSystem.h>
#include<TFile.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<TMath.h>
#include<TROOT.h>
#include<TStyle.h>
#include<TString.h>
#include<TMultiGraph.h>

// ---- CP error ----
#include<TEfficiency.h>//for Clopper Pearson
//--------

using namespace std;

#include "include/effi_constants.h"
#include "include/uf_style_macros.h"
#include "include/uf_common_tools.h"
#include "include/effi_histos.h"
#include "include/MapCounter.h"
#include "include/EffiEvent.h"

/*====================================================================================================
This script takes in an ntuple of RECO and RAW data and performes efficiency checks.  Within the
ntuple RECO muon segments have been matched to RAW LCTs.

Ntuplizer: plugins/TrigEff.C

Histograms, Trees, Branches defined in includes.

MapCounter class provided for easy counting.

Written by David Curry and GP DiGiovanni
====================================================================================================*/

TTree* recoMuons;
TTree* csctfTTree;
TLegend *tl;

#include "src/EffiEvent.cc"

void efficienciesGlobalMuonsTwoMatched(int  printLevel =     0,
				       int  graphLevel =     0,
				       bool isSave = !true) {

  gROOT->Clear();
  gStyle->SetOptStat(111111);  

  EffiEvent evt; //a class to initialize all variables and tree tasks
 
  //------- Choose which data set to analyze: A,B,C,D or all -----------------

  TFile* file = TFile::Open("merge_2012C.root");
  //TFile* file = TFile::Open("/cms/data/store/user/dcurry/trigeff/2012C_new/MinimumBiasTrigEffNtuple_new_C_allevt_782_1_oRj.root");
  //TFile* file = TFile::Open("dataset/MinBias_2012D_merge.root");
  //TFile* file = TFile::Open("dataset/MinBias_2012AB.root");
  //TFile* file = TFile::Open("dataset/MinBias_2012C_new.root");
  //TFile* file = TFile::Open("dataset/MinBias_2012D_new.root");


  evt.AttachToFile( file );

  //---------------------------------------------------------------------

  initScaleValues();

  // Fill Hist with phi window data
  for (int x=1 ; x<65 ; x++) {
    float gbl_eta = x*(1.5/64) + 0.9;
    hphiwindow -> Fill( gbl_eta , eta_dphi_ME1toME2[x]); 
  }

  // Creat mapcounter to replace all explicit global counters
  MapCounter _mapc;

  // counters
  int nMuons              = 0; // # muons in the denominator
  int nMuons_total        = 0; // # global muons total the show up in event record
  int nMuonsNoTrigger     = 0; // # muons with no trigg    er info
  //int nOverlapMuonsNoTrigger = 0;
  int num_evt_multi_muon = 0;
  int num_SizeTrk_is_Zero = 0; // # times no tracks made by CSCTF
  int me1_count           = 0; // # of LCT's per station
  int me2_count           = 0;
  int me3_count           = 0;
  int me4_count           = 0;
  int dr_eta              = 0;
  int dr_phi              = 0;
  int endcap_front_count  = 0;
  int endcap_back_count   = 0;
  int how_many_not_in_dr_eta   = 0;
  int how_many_not_in_dr_phi   = 0;
  int how_many_not_in_either   = 0;
  int num_muons_failed_deta    = 0;
  int num_muons_failed_dphi_zerotrk = 0;
  int num_muons_failed_deta_zerotrk = 0;
  int num_muons_failed_dphi    = 0;
  int num_muons_failed_deta_dphi_zerotrk = 0;
  int num_muons_failed_deta_dphi = 0;
  int num_muons_with_at_least_one_match = 0;
  int num_muons_with_at_least_one_match_zerotrk = 0;
  int num_overlap_DT_CSC_dphi_fail = 0;
  int num_multi_muon_event  = 0;
  int num_single_muon_event = 0;
  int how_many_ME1ME2 = 0;
  int how_many_ME1ME4 = 0;
  int how_many_ME1ME3 = 0;
  int how_many_ME2ME3 = 0;
  int how_many_ME2ME4 = 0;
  int how_many_ME3ME4 = 0;
  int num_muons_fail_different_sector = 0;
  int num_muons_fail_different_sector_zerotrk = 0;
  int num_noTrigger_sector_edge_zerotrk = 0;
  int num_noTrigger_sector_edge = 0;
  int failed_muons_with_event_trigger = 0;
  int num_overlap_muons = 0;
  int num_total_turn_on = 0;
  int num_overlap_sector_fail = 0;
  int num_muons_all_lcts_one_station = 0;
  int lctable = 0;
  int matched = 0;
  int segments = 0;
  int debug_matched = 0;
  int num_overlap_should_have_made = 0;

  // ------------------------------------------------------------------
  // Loop over the events
  // ------------------------------------------------------------------

  //for (int iEvt=0; iEvt < evt.GetEntries(); iEvt++) {

  for (int iEvt=0; iEvt < 100000; iEvt++) {

    evt.GetEntry(iEvt);

    if ( ( iEvt % 10000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);

    bool is_triggered   = false;
    int how_many_global = 0;
    bool event_trigger  = false;

    if (printLevel > 1) {
      cout << "=========================  Starting new event loop! ==================================" << endl;
      cout << "RUN = " <<  evt.Run   << endl;
      cout << "EVT = " <<  evt.Event << endl;
      cout << "LUMI = " << evt.Lumi << endl;
    }

    // Print out a List of all lcts in the event record for error checking
    if (printLevel > 2) {
      cout << "======= List of All Lcts in Event ======= " << endl;

      // Loop over event lcts(trigPrim)      
      for (int iLct=0; iLct < evt.SizeLCTs; iLct++) {
	
	cout << "Looping over Lct # " << iLct << endl; 
        cout << " Chamber = " << evt.lctChamber    -> at(iLct) << endl;
        cout << " Station = " << evt.lctStation    -> at(iLct) << endl << endl;
      }
    }


    // Block for testing/various histograms
    // Loop over tracks
    for (int iRaw=0; iRaw<evt.SizeTrk; iRaw++) {

     double phi1 = -999;
     double phi2 = -999;

     int nLcts = evt.NumLCTsTrk->at(iRaw);
      //loop over the LCT belonging to the track
     
      for (int iLCT=0; iLCT<nLcts; iLCT++) {
	
	if (evt.trLctStation[iRaw][iLCT] == 1) phi1 = evt.trLctglobalPhi[iRaw][iLCT];
	if (evt.trLctStation[iRaw][iLCT] == 2) phi2 = evt.trLctglobalPhi[iRaw][iLCT];

	//cout << "Phi 1 = " << phi1 << endl;
	//cout << "Phi 2 = " << phi2 << endl;
	
      }
      
      // csc dphi plots
      if (phi1 != -999 && phi2 != -999) {
	double dphi = phi1 - phi2;
	//if ( abs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - abs(dphi);	
	//double dphi = abs( abs( abs(phi1 - phi2) - 3.14 ) - 3.14 );
	hdphi_csc -> Fill(dphi);
      }
      
    }
    

    // First reco muon loop to check for single or dimuon status and how many muon triggers occurred.

    for (int iReco=0; iReco < evt.muonSize; iReco++) { 

      if ( evt.IsGlobalRecoMuon(iReco) ) {
	nMuons_total    += 1;
	how_many_global += 1;
      } else {
	if (printLevel > 2) cout << " Not a Global Muon. Go to Next Muon or event. " << endl;
	continue;
      }

      // --------------------------------------------------------------
      // routine to identify if the global muon has two segments, and at least one csc segment
      // --------------------------------------------------------------
      bool hasTwoSegsMatched = false;
      int counterSegs=0;
      segments += evt.muonNsegs->at(iReco);
   
      for (int iSeg=0; iSeg < evt.muonNsegs->at(iReco); iSeg++) {
		
	if (evt.muon_cscsegs_islctable[iReco][iSeg]) lctable ++;
	if (evt.muon_cscsegs_ismatched[iReco][iSeg]) matched ++ ;
	if (evt.muon_cscsegs_ismatched[iReco][iSeg])  counterSegs++;
	
	if (counterSegs > 1) {
	  hasTwoSegsMatched=true;      
	  continue;
	}
	
	if (is_triggered) continue;
	
	// look only at matched segments
	if (!evt.muon_cscsegs_ismatched[iReco][iSeg]) continue;
        int segLCTid = evt.muon_cscsegs_lctId[iReco][iSeg];
	if (segLCTid == -999 || segLCTid >= evt.SizeLCTs) continue;

	//loop over the CSCTF raw information
        for (int iRaw=0; iRaw<evt.SizeTrk; iRaw++) {
	  
          int nLcts = evt.NumLCTsTrk->at(iRaw);
	  //loop over the LCT belonging to the track
	  
          for (int iLCT=0; iLCT<nLcts; iLCT++) {
	    if (evt.trLctEndcap   [iRaw][iLCT] != evt.lctEndcap   ->at(segLCTid)) continue;
            if (evt.trLctSector   [iRaw][iLCT] != evt.lctSector   ->at(segLCTid)) continue;
            if (evt.trLctStation  [iRaw][iLCT] != evt.lctStation  ->at(segLCTid)) continue;
            if (evt.trLctglobalPhi[iRaw][iLCT] != evt.lctglobalPhi->at(segLCTid)) continue;
            if (evt.trLctglobalEta[iRaw][iLCT] != evt.lctglobalEta->at(segLCTid)) continue;
	    is_triggered = true;
	      
	  } // end lct loop
        } //end track loop
      } //end segment loop

      if (hasTwoSegsMatched) debug_matched++;

    }  //end reco loop

    // Does this event contain more than one global muon?
    if ( how_many_global > 1 ) {
      num_evt_multi_muon += 1;
    }
    
    // Assign this event single or Di muon status
    bool single_muon_evt = false;
    bool multi_muon_evt  = false;
    if ( how_many_global > 1 )  multi_muon_evt = true;
    if ( how_many_global == 1 ) single_muon_evt = true;
    
    // For Dimuon events make a plot of dR between muons
    if ( multi_muon_evt ) {
      double muon_eta[2];
      double muon_phi[2];
      int counter_muon = 0;
      // delta R = sqrt( (muon1_eta - muon2_eta)^2 - (muon1_phi - muon2_phi)^2 );
      
      for ( int iReco=0; iReco < evt.muonSize; iReco++) {
	
	if (!evt.IsGlobalRecoMuon(iReco)) continue;

	muon_eta[counter_muon] = evt.etaReco->at(iReco);
	muon_phi[counter_muon] = evt.phiReco->at(iReco);
	
	counter_muon +=1;
      }
      
      // calculate dR.  Fill Dimuon histograms
      double dphi = abs( abs( abs((muon_phi[0] - muon_phi[1])) - 3.14 ) - 3.14 ); 
      double dR = sqrt( (muon_eta[0] - muon_eta[1])*(muon_eta[0] - muon_eta[1]) + (dphi*dphi) );
      
      //      if ( (muon_eta[0] < 0 && muon_eta[1] > 0) ||
      //	   (muon_eta[0] > 0 && muon_eta[1] < 0) ) continue;
      hDR -> Fill(dR); 
      
      double deta = muon_eta[0] - muon_eta[1];
      h_dimuon_eta -> Fill(deta);
      
      //double dphi = muon_phi[0] - muon_phi[1];
      h_dimuon_phi -> Fill(dphi);      
    }
    
    // Assign this event trigger status - did a trigger occur on any muon in event
    if ( is_triggered ) event_trigger = true;
       
    // Now check each muon individually and fill hists 
    for (int iReco=0; iReco < evt.muonSize; iReco++) { 
      // -----------------------------------------------------------
      
      // different muon types .. checks built into EffiEvent
      bool isEndcap_plus = false;
      bool isGblMuon            = evt.IsGlobalRecoMuon (iReco);
      bool isStandAloneMuonOnly = evt.IsStandaloneMuon (iReco);
      bool isTrackerMuonOnly    = evt.IsTrackerMuon    (iReco);
      bool isStdAndTrkButNotGbl = evt.IsStdTrkNotGlb   (iReco);
    
      int counter=0;
      if (isGblMuon)            {hTypeMu->Fill(0); counter++;}
      if (isStandAloneMuonOnly) {hTypeMu->Fill(1); counter++;}
      if (isTrackerMuonOnly)    {hTypeMu->Fill(2); counter++;}
      if (isStdAndTrkButNotGbl) {hTypeMu->Fill(3); counter++;}
        
      // sanity check
      if (counter>1) {
        cout << "Error: the same muon fills two category!"
             << " Review their definitions\n";
        continue;
      }
      
      // pick your selection
      // --------------------------------------------------------------
      // here you select global, standalone only, tracker only
      if (!isGblMuon) continue;
      
      if (printLevel>1) {
        cout << "\n----------------- New Global ------------------------\n" << endl;
        cout << "How many global muons in this event =  " << how_many_global << endl ;
	cout << "RUN = "  << evt.Run   << endl;
	cout << "EVT = "  << evt.Event << endl;
	cout << "LUMI = " << evt.Lumi  << endl;
      }
      
      // --------------------------------------------------------------
      // routine to identify if the global muon has two segments matched
      // --------------------------------------------------------------
      bool hasTwoSegsMatched = false;
      bool is_overlap_muon   = false;
      bool has_csc_seg       = false;
      int counterSegs   = 0;
      int sectorTrk_phi = 0;

      // Set CSC Muon_seg values
      int  lctEtaBitSeg[MAX_SEGS_STD];
      int  lctPhiBitSeg[MAX_SEGS_STD];
      int  lctEndcapSeg[MAX_SEGS_STD];
      int lctStationSeg[MAX_SEGS_STD];
      int  lctSectorSeg[MAX_SEGS_STD];
      int      lctBxSeg[MAX_SEGS_STD];
      int    isdtlctSeg[MAX_SEGS_STD];

      // Set DT Muon_seg values
      int dt_lctStationSeg[MAX_SEGS_STD];
    
      if (printLevel>2) cout << "List of Matched Segments:\n" ;
      
      for (int iSeg=0;iSeg<evt.muonNsegs->at(iReco);iSeg++) {

	if (printLevel > 2) cout << "\n\nChecking Segment # " << iSeg << endl;
	if ( evt.muon_isdtseg[iReco][iSeg] == 0) { 
	  if (printLevel > 2) cout << "Segment is CSC" << endl; 
	  has_csc_seg = true;
	}
	
	if ( evt.muon_isdtseg[iReco][iSeg] == 1 ) { 
	  if (printLevel > 2) cout << "Segemnt is DT" << endl; 
	  if ( !(evt.muon_cscsegs_ismatched[iReco][iSeg] == -999) ) is_overlap_muon = true;
	}
	
	if (printLevel > 2) cout << "muon_cscsegs_ismatched = " << evt.muon_cscsegs_ismatched[iReco][iSeg]<< endl;
	 
        // Check is seg is matched to Lct 
	if (evt.muon_cscsegs_ismatched[iReco][iSeg]) {
	  int id = evt.muon_cscsegs_lctId[iReco][iSeg];
	  if (printLevel > 2) cout << "muon_cscsegs_lctId = " << id << endl;
	  if (id == -999) continue;
	  if (id > evt.SizeLCTs-1) continue;

	  if ( evt.muon_isdtseg[iReco][iSeg] == 0) {
	    lctEtaBitSeg[counterSegs]  = evt.lctglobalEta->at(id); 
	    lctPhiBitSeg[counterSegs]  = evt.lctglobalPhi->at(id); 
	    lctEndcapSeg[counterSegs]  = evt.lctEndcap   ->at(id); 
	    lctStationSeg[counterSegs] = evt.lctStation  ->at(id); 
	    lctSectorSeg[counterSegs]  = evt.lctSector   ->at(id); 
	    lctBxSeg[counterSegs]      = evt.lctBx       ->at(id);
	    isdtlctSeg[counterSegs]    = 0;

	    if (printLevel > 2) {
	      cout << " lctEtaBitSeg ="  << lctEtaBitSeg[counterSegs] << endl;
	      cout << " lctPhiBitSeg ="  << lctPhiBitSeg[counterSegs] << endl;
	      cout << " lctStationSeg =" << evt.lctStation->at(id) << endl;
	      cout << " lctSectorSeg ="  << evt.lctSector->at(id) << endl;
	      cout << " lctEndcapSeg ="  << evt.lctEndcap->at(id) << endl;
	      cout << " lctBxSeg ="      << lctBxSeg[counterSegs] << endl;
	    }

	         // Fill phi window hist
	    //int eta_bit = lctglobalEta->at(id) >> 1;
	    //hphiwindow -> Fill ( trig_prim_eta -> at(id) , eta_dphi_ME1toME2[eta_bit] );
	     		    


	  }
	  
	  if ( evt.muon_isdtseg[iReco][iSeg] == 1) {

	    if (printLevel > 1) cout << "---> here1" << endl;
	    dt_lctStationSeg[counterSegs] = evt.lctStation->at(id);
	    if (printLevel > 1) cout << "---> here2" << endl;
	    isdtlctSeg[counterSegs]       = 1;
	    if (printLevel > 1) cout << "---> here3" << endl;
	    if (printLevel > 2 && evt.muon_isdtseg[iReco][iSeg] == 1) {	    
	      cout << " lctEtaBitSeg ="  << lctEtaBitSeg[counterSegs] << endl;
	      cout << " lctPhiBitSeg ="  << lctPhiBitSeg[counterSegs] << endl;
	      cout << " lctStationSeg =" << evt.lctStation->at(id) << endl;
	      cout << " lctSectorSeg ="  << evt.lctSector->at(id) << endl;
	      cout << " lctWheelSeg ="   << evt.lctWheel->at(id);
 	      cout << " lctBxSeg ="      << evt.lctBx->at(id) << endl;
	    }
	  }
	  
	  counterSegs +=1;
	             	
	} // end ismatched 
      } // end loop over segs
	

      if (counterSegs > 1) hasTwoSegsMatched = true;

      if (printLevel > 2) {
	cout << "countersegs = " << counterSegs << endl;
	cout << " -> hasTwoSegsMatched = " << hasTwoSegsMatched << endl;
      }
      
      // only Global muons with at least two segments matched
      if (!hasTwoSegsMatched) continue;

      // Only Global Muons with at least one CSC segment
      if (!has_csc_seg) continue;

      // Only keep overlap Global Muons that have DT station1 and CSC station 1, 2, 3.
      bool wrong_station = false;
      if (is_overlap_muon) {
	for (int iSeg=0; iSeg<counterSegs; iSeg++) {
	  
	  if (isdtlctSeg[iSeg] == 1) 
	    if ( dt_lctStationSeg[iSeg] == 2 ||
		 dt_lctStationSeg[iSeg] == 3 ||
		 dt_lctStationSeg[iSeg] == 4) wrong_station = true;	
	  
	  if (isdtlctSeg[iSeg] == 0) 
	    if ( lctStationSeg[iSeg] == 4 ) wrong_station = true;	  
	}
      }
      
      if ( wrong_station ) continue;

      // ------- Only keep muons whos segments are matched to LCT's from the same endcap -------
      // ------- Only keep those muons whose segments have LCT's from more than one station ---------
      // Clear all station id's.... ONLY GOOD FOR CSC SEGMENTS!
      
      bool are_lcts_one_station = false;

	me1_count = 0;
	me2_count = 0;
	me3_count = 0;
	me4_count = 0;
	
	endcap_front_count = 0;
	endcap_back_count  = 0;
	
	// loop over segments
	for (int iSeg=0; iSeg < evt.muonNsegs->at(iReco); iSeg++) {
	  
	  // look only at matched segments to LCT's
	  if (!evt.muon_cscsegs_ismatched[iReco][iSeg]) continue;
	  
	  // look only at csc segs
	  if (evt.muon_isdtseg[iReco][iSeg] == 1) continue;

	  int segLCTid = evt.muon_cscsegs_lctId[iReco][iSeg];
	  if (segLCTid == -999) continue;
	  // count which station lct is in
	  switch(evt.lctStation->at(segLCTid)){
	  case 1: me1_count +=1; break;
	  case 2: me2_count +=1; break;
	  case 3: me3_count +=1; break;
	  case 4: me4_count +=1; break;
	  default: break; 
	  }
	  
	  // count which endcap lct is in
	  if (evt.lctEndcap->at(segLCTid) == 1) {
	    endcap_front_count +=1;
	    isEndcap_plus = true;
	  }

	  if (evt.lctEndcap->at(segLCTid) == 2) endcap_back_count +=1;
	  
	  
	} // end segment loop
	
	if (printLevel>1) {
	  cout << "me1_count = " << me1_count << endl
	       << "me2_count = " << me2_count << endl
	       << "me3_count = " << me3_count << endl
	       << "me4_count = " << me4_count << endl ;
	}
 	
	// keep only multi station muons
	if(me1_count == 0 && me2_count == 0 && me3_count == 0) {
	  if (printLevel>1) 
	    cout << " -> There are LCT's in one station only" << endl;
	  continue; 
	  are_lcts_one_station = true;
	}
	
	if(me1_count == 0 && me2_count == 0 && me4_count == 0) {
	  if (printLevel>1) 
	    cout << " -> There are LCT's in one station only" << endl;
	  continue; 
	  are_lcts_one_station = true;
	}
	
	if(me1_count == 0 && me3_count == 0 && me4_count == 0) {
	  if (printLevel>1) 
	    cout << " -> There are LCT's in one station only" << endl;
	  continue;
	  are_lcts_one_station = true;
	}
	
	if(me2_count == 0 && me3_count == 0 && me4_count == 0) {
	  if (printLevel>1) 
	    cout << " -> There are LCT's in one station only" << endl;
	  continue; 
	  are_lcts_one_station = true;
	}
	
	// Keep only the muons with two or more lcts in one endcap
	if (endcap_front_count < 2 ) {
	  if (endcap_back_count  < 2 ) {  
	    
	    if (printLevel > 1 ) {
	      cout << " endcap_front_count = " << endcap_front_count << endl;
	      cout << " -> There are not at least two LCT's in the same endcap" << endl;
	    }
	    continue; 
	  }
	} 

      
      // --------------  Filters to look only at certain event/muons ---------------------------------------------------------------------------
      // =======================================================================================================================================
      float eta = abs(evt.etaReco -> at(iReco));
      if ( eta < 0.9 || eta > 2.4 ) continue;
      
      // Single or Di muon event Filter
      if ( !single_muon_evt ) continue;
      //if ( !multi_muon_evt  ) continue;
      
      // Apply a fliter to only look at muons in certain eta ranges
      //if ( abs(evt.etaReco -> at(iReco)) > 1.3 ) continue ; 
      
      //if ( (abs(evt.etaReco -> at(iReco)) <= 1.3) || (abs(evt.etaReco -> at(iReco)) >= 1.7) ) continue;
      
      //if ( (abs(evt.etaReco -> at(iReco)) < 1.7) || (abs(evt.etaReco -> at(iReco)) >= 2.1) ) continue;
      
      //if ( abs(evt.etaReco -> at(iReco)) < 2.1 ) continue ;
      
      // Overlap Muon Filter
      //if ( is_overlap_muon ) continue; 
      
      // Station and Chamber Filter
      //if ( evt.lctStation -> at(iReco) != 1 && evt.lctChamber -> at(iReco) != 3) continue;
      
      // Pt
      if (evt.ptReco -> at(iReco) < 1) continue;
      
      //fill the denominator    
      if (printLevel>1) cout << " -> Denominator is filled " << endl;
      
      if (eta < 1.3 && !is_overlap_muon) _mapc << "Denominator: eta < 1.3";
      if (eta <= 1.7 && eta >= 1.3)      _mapc << "Denominator: 1.3 <= eta <= 1.7";
      if (eta <= 2.1 && eta >= 1.7)      _mapc << "Denominator: 1.7 <= eta <= 2.1";
      if (eta > 2.1)                     _mapc << "Denominator: eta > 2.1";
      if (is_overlap_muon)               _mapc << "Denominator: Overlap Muon";
      nMuons+=1;
      
     if ( is_overlap_muon ) {
       num_overlap_muons +=1;
       if (printLevel > 1) cout << "This is an overlap muon" << endl; 
     }
     else if(printLevel > 1) cout << "This is not an overlap muon" << endl; 
     
     hNMuonsVsPt  -> Fill ( evt.ptReco ->at(iReco) );
     hNMuonsVsEta -> Fill ( evt.etaReco->at(iReco) );
     hNMuonsVsPhi -> Fill ( evt.phiReco->at(iReco) );
     overlap_eta  -> Fill ( evt.etaReco->at(iReco) );

     if (!isEndcap_plus) new_hNMuonsVsPhi -> Fill ( evt.phiReco->at(iReco) );
     
     
     bool isTriggered  = false;
     float pt_measured = -999; // used for pt turn on curves
     float etaWinner   = -999;
     int modeWinner    = -999;
     
     //-------------- Numerator Time! --------------------------------------------------------------------------------
     // loop over segments
     
     for (int iSeg=0; iSeg < evt.muonNsegs->at(iReco); iSeg++) {
             
       if (isTriggered) continue;
       
       // look only at matched segments
       if (!evt.muon_cscsegs_ismatched[iReco][iSeg]) continue;
        	
	int segLCTid = evt.muon_cscsegs_lctId[iReco][iSeg];	
        if (segLCTid == -999 || segLCTid > evt.SizeLCTs) continue;	
     
	if (printLevel > 1) cout << "******************* Start Trigger Check! *********************************************\n\n";

	if ( evt.muon_isdtseg[iReco][iSeg] == 1 ) { if (printLevel > 1) cout << "Muon Segment is DT" << endl; }
	else if (printLevel > 1) cout << "Muon Segment is CSC" << endl;

	if (printLevel>1 && evt.muon_isdtseg[iReco][iSeg] == 0) {
          cout  << "******************* Segment LCT ID's *************************************************\n";
	  cout  << "muon_cscsegs_lctId[iReco][iSeg]  = " << segLCTid          << std::endl;
	  cout  << "lctEndcap->at(segLCTid)   =" << evt.lctEndcap->at(segLCTid)   << std::endl
		<< "lctSector->at(segLCTid)   =" << evt.lctSector->at(segLCTid)   << std::endl
		<< "lctStation->at(segLCTid)  =" << evt.lctStation->at(segLCTid)  << std::endl
		<< "lctglobalPhi->at(segLCTid)=" << evt.lctglobalPhi->at(segLCTid)<< std::endl
		<< "lctglobalEta->at(segLCTid)=" << evt.lctglobalEta->at(segLCTid)<< std::endl
		<< "lctBx->at(segLCTid)="        << evt.lctBx->at(segLCTid)       << std::endl;    
	}
	
	if (printLevel>1 && evt.muon_isdtseg[iReco][iSeg] == 1) {
          cout  << "******************* Segment LCT ID's *************************************************\n";
	  cout  << "muon_cscsegs_lctId[iReco][iSeg]  = " << segLCTid          << std::endl;
          cout  << "lctSector->at(segLCTid)   =" << evt.lctSector->at(segLCTid)   << std::endl
		<< "lctStation->at(segLCTid)  =" << evt.lctStation->at(segLCTid)  << std::endl
		<< "lctChamber->at(segLCTid)  =" << evt.lctChamber->at(segLCTid)  << std::endl
		<< "lctglobalPhi->at(segLCTid)=" << evt.lctglobalPhi->at(segLCTid)<< std::endl
		<< "lctglobalEta->at(segLCTid)=" << evt.lctglobalEta->at(segLCTid)<< std::endl
		<< "lctBx->at(segLCTid)="        << evt.lctBx->at(segLCTid)       << std::endl;
        }
	
	if (printLevel>1)
	  cout << " Number of CSCTF Tracks in this event = " << evt.SizeTrk << endl;
	
        //loop over the CSCTF raw track information
	
        for (int iRaw=0; iRaw < evt.SizeTrk; iRaw++) {   
	  if (printLevel > 1) cout << "\n---- Looping over Track # " << iRaw << endl;
	  
	  int nLcts = evt.NumLCTsTrk->at(iRaw);
	  if (printLevel>1) cout << " nLcts = " << nLcts << endl;

	  //loop over the LCTs belonging to the track
          for (int iLCT=0; iLCT<nLcts; iLCT++) {
	    
	    // Fill 2D plot with LCTs vs eta for the track
	    if ( iLCT == 0 ) {
	      float gbl_eta = evt.trLctglobalEta[iRaw][iLCT] * (1.5/128) + 0.9;
	      twod_lct_eta -> Fill ( nLcts , gbl_eta );  
	    }
	    

	    // First Check CSC Lcts
	    if ( evt.isdttrlct[iRaw][iLCT] == 0 ) {
	      if (printLevel>1) {
		cout << "Lct is CSC" << endl;
		cout << "----------------------CSCTF Track ID's ------------------------------------\n";
		cout       << "iLCT=" << iLCT << endl;
		cout       << "trLctEndcap[iRaw][iLCT]   =" << evt.trLctEndcap[iRaw][iLCT]    << std::endl
			   << "trLctSector[iRaw][iLCT]   =" << evt.trLctSector[iRaw][iLCT]    << std::endl
			   << "trLctStation[iRaw][iLCT]  =" << evt.trLctStation[iRaw][iLCT]   << std::endl
			   << "trLctglobalPhi[iRaw][iLCT]=" << evt.trLctglobalPhi[iRaw][iLCT] << std::endl
			   << "trLctglobalEta[iRaw][iLCT]=" << evt.trLctglobalEta[iRaw][iLCT] << std::endl
			   << "trLctBx[iRaw][iLCT]       =" << evt.trLctBx[iRaw][iLCT]        << std::endl;
	      }
	      
	      if (evt.trLctEndcap[iRaw][iLCT]    != evt.lctEndcap->at(segLCTid))    continue;
	      if (evt.trLctSector[iRaw][iLCT]    != evt.lctSector->at(segLCTid))    continue;
	      if (evt.trLctStation[iRaw][iLCT]   != evt.lctStation->at(segLCTid))   continue;
	      if (evt.trLctglobalPhi[iRaw][iLCT] != evt.lctglobalPhi->at(segLCTid)) continue;
	      if (evt.trLctglobalEta[iRaw][iLCT] != evt.lctglobalEta->at(segLCTid)) continue;
	      if (evt.trLctBx[iRaw][iLCT]        != evt.lctBx->at(segLCTid))        continue;
	      
	      if (printLevel > 1) cout << " ----> IS TRIGGERD!" << endl;      
	      
	      isTriggered = true;
	      pt_measured = evt.PtTrk   -> at(iRaw);
	      modeWinner  = evt.ModeTrk -> at(iRaw);
	      etaWinner   = evt.EtaTrk  -> at(iRaw);

	    } // end CSC Lct trigger check
	    
	    // Check Dt Lcts
	    if ( evt.isdttrlct[iRaw][iLCT] == 1) {
	      
	      if (printLevel > 1) {
		cout << "Lct is DT" << endl;
		cout << "----------------------CSCTF Track ID's ------------------------------------\n";
		cout       << "iLCT=" << iLCT << endl;
		cout       << "trLctSector[iRaw][iLCT]   =" << evt.trLctSector[iRaw][iLCT]    << std::endl
			   << "trLctStation[iRaw][iLCT]  =" << evt.trLctStation[iRaw][iLCT]   << std::endl
			   << "trLctChamber[iRaw][iLCT]  =" << evt.trLctChamber[iRaw][iLCT]   << std::endl
			   << "trLctglobalPhi[iRaw][iLCT]=" << evt.trLctglobalPhi[iRaw][iLCT] << std::endl
			   << "trLctglobalEta[iRaw][iLCT]=" << evt.trLctglobalEta[iRaw][iLCT] << std::endl
			   << "trLctBx[iRaw][iLCT]       =" << evt.trLctBx[iRaw][iLCT]        << std::endl;
	      }
	      
	      if (evt.trLctSector[iRaw][iLCT]    != evt.lctSector->at(segLCTid))    continue;
              if (evt.trLctStation[iRaw][iLCT]   != evt.lctStation->at(segLCTid))   continue;
	      if (evt.trLctChamber[iRaw][iLCT]   != evt.lctChamber->at(segLCTid))   continue;
              if (evt.trLctglobalPhi[iRaw][iLCT] != evt.lctglobalPhi->at(segLCTid)) continue;
              if (evt.trLctglobalEta[iRaw][iLCT] != evt.lctglobalEta->at(segLCTid)) continue;
              if (evt.trLctBx[iRaw][iLCT]        != evt.lctBx->at(segLCTid))        continue;

              if (printLevel > 1) cout << " ----> IS TRIGGERD!" << endl;

              isTriggered = true;
	      pt_measured = evt.PtTrk   -> at(iRaw);
	      modeWinner  = evt.ModeTrk -> at(iRaw);	      
	      etaWinner   = evt.EtaTrk  -> at(iRaw);


	    } // end DT Lct trigger check
          } // end lct loop
        } //end track loop
     } //end segment loop
       

     // ================================================================================================================================================================
     // ----------------------- If the muon was triggered create Pt Turn On Curves --------------------------------------------------------------

     // only analyze tracks with quality = 3 
     bool quality_3 = false;
     if (modeWinner == 2 && etaWinner > 1.2 && etaWinner < 2.1)     quality_3 = true;
     if (modeWinner == 4 && etaWinner < 2.1)                        quality_3 = true;
     if (modeWinner == 5)                                           quality_3 = true;
     if (modeWinner == 11 || modeWinner == 12 || modeWinner == 14 ) quality_3 = true; 
     
     if ( isTriggered && quality_3 ) {
       
       num_total_turn_on += 1;
       // Specify Pt cut off(GeV) 5, 7, 10
       // Access Pt True(tracker Pt) and Pt Measured from trigger check above
       float pt_true  = evt.ptReco -> at(iReco);
       
       if ( (pt_true > 5 && pt_measured > 5) ||
	    (pt_true < 5 && pt_measured > 5) ) pt_turn_pass_5 -> Fill(pt_true);
       
       else pt_turn_fail_5 -> Fill(pt_true);
	
       if ( (pt_true > 7 && pt_measured > 7) ||
	    (pt_true < 7 && pt_measured > 7) ) pt_turn_pass_7 -> Fill(pt_true);
       
       else pt_turn_fail_7 -> Fill(pt_true);
       
       if ( (pt_true > 10 && pt_measured > 10) ||
	    (pt_true < 10 && pt_measured > 10) ) pt_turn_pass_10 -> Fill(pt_true);
       
       else pt_turn_fail_10 -> Fill(pt_true);

     } // end isTriggered block


     // ================================================================================================================================================================
     // ----------------------- If the muon was not triggered calculate segment windows and find why it did not trigger ------------------------------------------------
     
     if ( !isTriggered ) {

       if (eta < 1.3 && !is_overlap_muon) _mapc << "No Trigger: eta < 1.3";
       if (eta <= 1.7 && eta >= 1.3)      _mapc << "No Trigger: 1.3 <= eta <= 1.7";
       if (eta <= 2.1 && eta >= 1.7)      _mapc << "No Trigger 1.7 <= eta <= 2.1";
       if (eta > 2.1)                     _mapc << "No Trigger: eta > 2.1";
       if (is_overlap_muon)               _mapc << "No Trigger: Overlap Muon";
       
       nMuonsNoTrigger +=1;
       if ( event_trigger )   failed_muons_with_event_trigger += 1;
       if ( single_muon_evt ) num_single_muon_event  += 1;
       if ( multi_muon_evt )  num_multi_muon_event   += 1;
       if ( evt.SizeTrk == 0 )    num_SizeTrk_is_Zero    +=1;
       
       if (printLevel > 1) cout << "\nHouston we have a problem" << endl;
       
       
       // Fill Mode vs. Eta plots for tracks that did not match
       int ModeTrk  = -999;
       float EtaTrk = -999;
       int num_trks_fail = 0; 
	 
       for (int iRaw=0; iRaw < evt.SizeTrk; iRaw++) {
	 
	 ModeTrk = evt.ModeTrk -> at(iRaw);
	 EtaTrk  = abs( evt.EtaTrk  -> at(iRaw) );
	 
	 //cout << "ModeTrk = " <<  ModeTrk << endl;
	 //cout << "EtaTrk = " <<   EtaTrk << endl;
	 num_trks_fail += 1;
	 
	 switch (ModeTrk) {
	 
	 case 1:  hmode1  -> Fill(EtaTrk);
	 case 2:  hmode2  -> Fill(EtaTrk);
	 case 3:  hmode3  -> Fill(EtaTrk);
	 case 4:  hmode4  -> Fill(EtaTrk);
	 case 5:  hmode5  -> Fill(EtaTrk);
	 case 6:  hmode6  -> Fill(EtaTrk);
	 case 7:  hmode7  -> Fill(EtaTrk);
	 case 8:  hmode8  -> Fill(EtaTrk);
	 case 9:  hmode9  -> Fill(EtaTrk);
	 case 10: hmode10 -> Fill(EtaTrk);
	 case 11: hmode11 -> Fill(EtaTrk);
	 case 12: hmode12 -> Fill(EtaTrk);
	 case 13: hmode13 -> Fill(EtaTrk);
	 case 14: hmode14 -> Fill(EtaTrk);

	 }
	 

       } // end loop over tracks
  
       bool seg_pair_in_window = false;
       bool seg_pair_in_window_but_too_close = false;
       how_many_not_in_dr_eta = 0;
       how_many_not_in_dr_phi = 0;
       how_many_not_in_either = 0;
       int num_segs_different_sector  = 0;
       int num_seg_pair_in_window = 0;
       int num_segment_pairs = 0;
       int num_segs_same_station = 0;
       int num_segs_not_same_endcap = 0;
       int pairMode = 666;
       bool dt_csc_in_window = false;
       
       // If muon is not an overlap muon then do the following.  Analyze overlap muons later
       if ( !is_overlap_muon ) {
	 
	 if (printLevel > 1) cout << "Now checking non overlap muon segments" << endl;
	 
	 for (int i=0; i<counterSegs; i++) {
	   
	   for (int j=i+1; j<counterSegs; j++) {
	     
	     // Dont compare a segment with itself.  Dont compare two segs in different sectors or same station. They will never match.  
	     // Also dont compare segments in different endcaps(for halo muons)
	     //if ( i == j ) continue;
	     num_segment_pairs += 1;
	     if ( lctSectorSeg[i] != lctSectorSeg[j] ) {
	       num_segs_different_sector +=1;
	       if ( printLevel > 1 ) cout << "Segment " << i << " and " << j << " are in different sector. Skip further comparisons." << endl;
	       continue;
	     }
	     if (lctStationSeg[i] == lctStationSeg[j]) {
	       if ( printLevel > 1 ) cout << "Segment " << i << " and " << j << " are in same station. Skip further comparisons." << endl;
	       num_segs_same_station +=1;
	       continue;
	     }
	     if ( lctEndcapSeg[i] != lctEndcapSeg[j] ) {
	       if ( printLevel > 1 ) cout << "Segment " << i << " and " << j << " are in different endcaps. Skip further comparisons." << endl;
	       num_segs_not_same_endcap +=1;
	       continue;
	     }
	     
	     // access the correct eta to dr phi window from arrays
	     int eta_bit = lctEtaBitSeg[j] >> 1;
	     bool in_dr_eta = false;
	     bool in_dr_phi = false;
	     bool local_seg_pair_in_window = false;
	     
	     // Count how many times segment LCT's are not in eta/phi windows.
	     // First the eta windows. Different stations have different windows.
	     if ( (lctStationSeg[i] == 1 && lctStationSeg[j] == 2) || (lctStationSeg[i] == 1 && lctStationSeg[j] == 3)) {
	       if ( abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) <= 4) in_dr_eta = true; 
	     }
	     
	     else if ( abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) <= 6) in_dr_eta = true; 
	     
	     // Now the phi windows which require a vector of phi window values depending on eta and ME station combo
	     // Get eta bit for this segment pair. Find asscoiated phi window
	     if ( lctStationSeg[i] == 1 && lctStationSeg[j] == 2 ) {
	       int dphi_window = eta_dphi_ME1toME2[eta_bit];	       
	       if ( abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] ) <= dphi_window ) in_dr_phi = true; 
	     }
	     
	     if ( (lctStationSeg[i] == 1 && lctStationSeg[j] == 3) || (lctStationSeg[i] == 1 && lctStationSeg[j] == 4) ) {
	       int dphi_window = eta_dphi_ME1toME3[eta_bit];
	       if ( abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] ) <= dphi_window ) in_dr_phi = true; 
	     }
	     
	     if ( (lctStationSeg[i] == 2 && lctStationSeg[j] == 3) || (lctStationSeg[i] == 2 && lctStationSeg[j] == 4) || (lctStationSeg[i] == 3 && lctStationSeg[j] == 4)) 
	       if (abs (lctStationSeg[i]/4 - lctStationSeg[j]/4) <= 127) in_dr_phi = true; 
	     
	     
	     // Fill hist for phi cooridiante of pair that failed phi window
	     if ( !in_dr_phi ) hNoTrigger_phicuts -> Fill ( evt.phiReco->at(iReco) );
	     	     
	     // Now tell me if these segments were in the windows and which ones.  Keep track of pairs that were in windows.
	     if ( (in_dr_phi) && (in_dr_eta) ) {
	       local_seg_pair_in_window = true;
	       num_seg_pair_in_window +=1;
	       if ( printLevel > 1 ) cout << "  Segment " << i << " and " << j << " are in eta and phi window " << endl;
	     }
	     else if ( (in_dr_phi) ) {
	       how_many_not_in_dr_eta += 1;
	       if ( printLevel > 1 ) cout << " Segment " << i << " and " << j << " are in phi window, but not eta " << endl; 
	     }
	     else if ( (in_dr_eta) ) {
	       how_many_not_in_dr_phi += 1;
	       if ( printLevel > 1 ) cout << " Segment " << i << " and " << j << " are in eta window, but not phi " << endl; 
	     }
	     if ( (!in_dr_eta) && (!in_dr_phi) ) {
	       how_many_not_in_either += 1;
	       if ( printLevel > 1 ) cout << " Segment " << i << " and " << j << " are not in eta or phi window " << endl;
	     }
	     
	     // Which stations are segment i and segment j in? Keep track of total phi failures between station pairs.
	     if ( !in_dr_phi ) {
	       if ( (lctStationSeg[i] == 1 || lctStationSeg[j] == 1) && (lctStationSeg[i] == 2 || lctStationSeg[j] == 2 ) )  how_many_ME1ME2 += 1;
	       if ( (lctStationSeg[i] == 1 || lctStationSeg[j] == 1) && (lctStationSeg[i] == 3 || lctStationSeg[j] == 3 ) )  how_many_ME1ME3 += 1;
	       if ( (lctStationSeg[i] == 1 || lctStationSeg[j] == 1) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4 ) )  how_many_ME1ME4 += 1;
	       if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 3 || lctStationSeg[j] == 3 ) )  how_many_ME2ME3 += 1;
	       if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4 ) )  how_many_ME2ME4 += 1;
	       if ( (lctStationSeg[i] == 3 || lctStationSeg[j] == 3) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4 ) )  how_many_ME3ME4 += 1;
	     }
	     
	     // ------- What is seg pair mode? Used by CSCTF to determine sector edging... -------------------------------
	     if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 3 || lctStationSeg[j] == 3) )      pairMode = 8;
	     
	     else if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4) ) pairMode = 9;
	     
	     else if ( (lctStationSeg[i] == 3 || lctStationSeg[j] == 3) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4) ) pairMode = 10;

	   	     
	     if ( printLevel > 1 ) {
	       cout << "Segment " << i << " and segment " << j << " are in stations ME" << lctStationSeg[i] << " and ME" << lctStationSeg[j]  << endl;
	       cout << "Segment " << i << " and segment " << j << " have an eta difference of  = " << abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) << endl;
	       cout << "Segment " << i << " and segment " << j << " have a phi difference of  = " << abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] )  << endl;
	     }
	     
	     // calculate dr phi/eta between all segment LCT combos ; fill histograms
	     dr_eta = abs( lctEtaBitSeg[i]-lctEtaBitSeg[j]);
	     dr_phi = abs( lctPhiBitSeg[i]-lctPhiBitSeg[j]);
	     
	     hdr_eta  -> Fill (dr_eta = abs( lctEtaBitSeg[i]-lctEtaBitSeg[j]) );
	     hdr_phi  -> Fill (dr_phi = abs( lctPhiBitSeg[i]-lctPhiBitSeg[j]) );
	     
	     // ----- pairs(or triplets) that made tracks to close to sector edge---------------------------------------------
	     // ----- track phi is taken from lct in station 2,3,4 in that priority -----------------------------
	     if ( (lctStationSeg[i] < lctStationSeg[j]) && (lctStationSeg[i] != 1) ) {
	       sectorTrk_phi = lctPhiBitSeg[i];
	     }
	     else {sectorTrk_phi = lctPhiBitSeg[j]; }
	     
	     // ------- Check to see if track failed due to sector edging ----------------------------
	     if ( local_seg_pair_in_window ) {
	       if (printLevel > 1 ) {
		 cout << "Segments could have made a track. Track phi is = " << sectorTrk_phi << ".  Checking to see if this is too close to sector edge " << endl;
		 cout << "Track mode = " << pairMode << ". If mode is 5, 8, 9, 10, and too close to sector edge, cancel the track. " << endl;
	       }
	       if ( ( sectorTrk_phi < 128 || sectorTrk_phi > (4095-128) ) &&
		    (  (pairMode == 10 || pairMode == 8 || pairMode == 9) ) ) {
		 seg_pair_in_window_but_too_close = true;
		 if (printLevel > 1) {
		   cout << "This segment pair is in windows but track was cancelled due to track being to close to sector edge(128 phi bits)" << endl;
		 }
	       } 
	       else {
		 seg_pair_in_window = true;
		 if (printLevel > 1) cout << "This segment pair was in eta/phi windows and not too close to sector edge." << endl;
	       }
	       if (printLevel > 1) cout << "Distance from sector edge is = " << abs(sectorTrk_phi - 4096) << ".  Cant be < 128 or > 3968" << endl << endl;
	     }
	     
	   } // end j
	 }  // end i
	 // end segment loop
       
	 // Report all segment comparison results
	 if ( printLevel > 1 ) {
	   cout << "This muon had " <<   num_seg_pair_in_window     << " segment pairs in eta/phi window" << endl;
	   cout << "This muon had " <<  how_many_not_in_dr_eta      << " segment pairs not in eta window " << endl;
	   cout << "This muon had " <<  how_many_not_in_dr_phi      << " segment pairs not in phi window " << endl;
	   cout << "This muon had " <<  how_many_not_in_either      << " segment pairs not in phi and eta window " << endl;
	   cout << "This muon had " <<  num_segs_different_sector   << " segment pairs in different sectors " << endl;
	 }
	 // Make sure every muon not triggered is counted for one reason only.
	 int fail_counter = 0;	 
	 //int eta = abs(evt.etaReco -> at(iReco));
	 
	 if ( are_lcts_one_station ) {
	   if (printLevel > 1) cout << " This muon failed due to all CSC lcts being in one station " << endl;
	   num_muons_all_lcts_one_station += 1;
	   
	   if (eta < 1.3)                _mapc << "Fail All LCTS one station: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Fail All LCTS one station: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Fail All LCTS one station: 1.7 <= eta <= 2.1";
           if (eta > 2.1)                _mapc << "Fail All LCTS one station: eta > 2.1";
	   if (is_overlap_muon)          _mapc << "Fail All LCTS one station: Overlap Muon";
	
	   fail_counter +=1;
	 }

	 if ( seg_pair_in_window ) {
	   if ( printLevel > 1 ) cout << " This muon had a segment pair eta/phi window and not too close to sector edge. Should have made a track." << endl;	   
	   if ( evt.SizeTrk == 0 ) num_muons_with_at_least_one_match_zerotrk += 1;
	   if ( evt.SizeTrk > 0  ) num_muons_with_at_least_one_match += 1;

	   if (eta < 1.3)                _mapc << "Should Have made a Muon: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Should Have made a Muon: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Should Have made a Muon: 1.7 <= eta <= 2.1";
           if (eta > 2.1)                _mapc << "Should Have made a Muon: eta > 2.1";
           if (is_overlap_muon)          _mapc << "Should Have made a Muon: Overlap Muon";

	   fail_counter +=1;
	 }
	 
	 // Keep track of muons that failed for just eta, just phi window, or both
	 if ( (how_many_not_in_dr_eta > 0) && (how_many_not_in_dr_phi == 0) && (how_many_not_in_either == 0) &&
	      (!seg_pair_in_window) && (!seg_pair_in_window_but_too_close) ) {
	   if ( evt.SizeTrk == 0 ) num_muons_failed_deta_zerotrk += 1;
	   if ( evt.SizeTrk > 0  ) num_muons_failed_deta += 1;
	   if ( printLevel > 1 ) cout << " This muon failed due to eta window only " << endl;
	   
	   if (eta < 1.3)                _mapc << "Fail Eta/Phi Windows: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Fail Eta/Phi Windows: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Fail Eta/Phi Windows: 1.7 <= eta <= 2.1";
	   if (eta > 2.1)                _mapc << "Fail Eta/Phi Windows: eta > 2.1";
	   if (is_overlap_muon)          _mapc << "Fail Eta/Phi Windows: Overlap Muon";	  
	   
	   fail_counter +=1;
	 }
        
	 if ( (how_many_not_in_dr_phi > 0) && (how_many_not_in_dr_eta == 0) && (how_many_not_in_either == 0) &&
	      (!seg_pair_in_window) && (!seg_pair_in_window_but_too_close) ) {
	   if ( evt.SizeTrk == 0 ) num_muons_failed_dphi_zerotrk += 1;
	   if ( evt.SizeTrk > 0  ) num_muons_failed_dphi += 1;
	   if ( printLevel > 1 ) cout << " This muon failed due to phi window only " << endl;
	   
	   if (eta < 1.3)                _mapc << "Fail Eta/Phi Windows: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Fail Eta/Phi Windows: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Fail Eta/Phi Windows: 1.7 <= eta <= 2.1";
           if (eta > 2.1)                _mapc << "Fail Eta/Phi Windows: eta > 2.1";
           if (is_overlap_muon)          _mapc << "Fail Eta/Phi Windows: Overlap Muon";
	   
	   fail_counter +=1;
	 }

	 if ( ( how_many_not_in_either > 0 && !seg_pair_in_window && !seg_pair_in_window_but_too_close ) ||
	      ( !seg_pair_in_window && how_many_not_in_dr_phi > 0 && how_many_not_in_dr_eta > 0 && how_many_not_in_either == 0 && !seg_pair_in_window_but_too_close ) ) {
	   if ( evt.SizeTrk == 0 ) num_muons_failed_deta_dphi_zerotrk += 1;
	   if ( evt.SizeTrk > 0  ) num_muons_failed_deta_dphi += 1;
	   if ( printLevel > 1 ) cout << " This muon failed due to phi and eta window " << endl;
	   
	   if (eta < 1.3)                _mapc << "Fail Eta/Phi Windows: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Fail Eta/Phi Windows: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Fail Eta/Phi Windows: 1.7 <= eta <= 2.1";
           if (eta > 2.1)                _mapc << "Fail Eta/Phi Windows: eta > 2.1";
           if (is_overlap_muon)          _mapc << "Fail Eta/Phi Windows: Overlap Muon";
	   
	   fail_counter +=1;
	 }         
	 
	 if ( (!seg_pair_in_window) && (seg_pair_in_window_but_too_close) ) {
	   if ( printLevel > 1 ) cout << " This muon failed due to segment pair in eta/phi window but is too close to sector edge." << endl;
	   if ( evt.SizeTrk == 0 ) num_noTrigger_sector_edge_zerotrk += 1;
	   if ( evt.SizeTrk > 0  ) num_noTrigger_sector_edge         += 1;

	   if (eta < 1.3)                _mapc << "Fail Sector Edge: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Fail Sector Edge: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Fail Sector Edge: 1.7 <= eta <= 2.1";
           if (eta > 2.1)                _mapc << "Fail Sector Edge: eta > 2.1";
           if (is_overlap_muon)          _mapc << "Fail Sector Edge: Overlap Muon";
	   
	   fail_counter +=1;
	 }
	 
	 if ( num_segs_different_sector >= (num_segment_pairs - num_segs_same_station - num_segs_not_same_endcap) ) {
	   if ( evt.SizeTrk == 0 )  num_muons_fail_different_sector_zerotrk += 1;
	   if ( evt.SizeTrk > 0  )  num_muons_fail_different_sector         += 1;
	   if ( printLevel > 1 ) cout << " This muon failed due to all LCTs pairs in different sectors " << endl;

	   if (eta < 1.3)                _mapc << "Fail Different Sectors: eta < 1.3";
           if (eta <= 1.7 && eta >= 1.3) _mapc << "Fail Different Sectors: 1.3 <= eta <= 1.7";
           if (eta <= 2.1 && eta >= 1.7) _mapc << "Fail Different Sectors: 1.7 <= eta <= 2.1";
           if (eta > 2.1)                _mapc << "Fail Different Sectors: eta > 2.1";
           if (is_overlap_muon)          _mapc << "Fail Different Sectors: Overlap Muon";
         
	   fail_counter +=1;
	 }

	 
	 if ( event_trigger && printLevel > 1) cout << " This muon was part of an event that had a trigger on another muon " << endl;
	 
	 // Throw a warning if failure counter is > 1.
	 if (fail_counter > 1) 
	   cout << evt.Event << "-----> Failure Counter is more than one for this muon!!!!" << endl;  
	 if (fail_counter == 0)
	   cout << evt.Event << "-----> Failure Counter is zero for this muon!!!!" << endl;
 	 
     } // end !is_muon_overlap block

     
       // ==============================================================================================================================================
       
       // -------  Now check for failure reason if muon is overlap --------------------------------
       int num_csc_muonsegs = 0;
       
       if ( is_overlap_muon ) {
	 
	 if (printLevel > 1) cout << "Now checking overlap muon CSC segments" << endl;
	 // OK, what were going to do is first check CSC seg combos.  Then DT-CSC combos
	 
	 for (int i=0; i<counterSegs; i++) {
	   if ( isdtlctSeg[i] == 0) num_csc_muonsegs += 1;
	   if ( isdtlctSeg[i] == 1) continue;
	   
           for (int j=i+1; j<counterSegs; j++) {
	     if ( isdtlctSeg[j] == 1) continue;
	     
	     // Dont compare a segment with itself.  Dont compare two segs in different sectors or same station. They will never match.
             // Also dont compare segments in different endcaps(for halo muons)
	     num_segment_pairs += 1;
	     if ( lctSectorSeg[i] != lctSectorSeg[j] ) {
	       num_segs_different_sector +=1;
               if ( printLevel > 1 ) cout << "Segment " << i << " and " << j << " are in different sector. Skip further comparisons." << endl;
	       continue;
	     }
             if (lctStationSeg[i] == lctStationSeg[j]) {
               if ( printLevel > 1 ) cout << "Segment " << i << " and " << j << " are in same station. Skip further comparisons." << endl;
               num_segs_same_station +=1;
               continue;
             }
             if ( lctEndcapSeg[i] != lctEndcapSeg[j] ) {
               if ( printLevel > 1 ) cout << "Segment " << i << " and " << j << " are in different endcaps. Skip further comparisons." << endl;   
               continue;
             }
             // access the correct eta to dr phi window from arrays
             int eta_bit = lctEtaBitSeg[j] >> 1;
             bool in_dr_eta = false;
             bool in_dr_phi = false;
             bool local_seg_pair_in_window = false;

	     // Count how many times segment LCT's are not in eta/phi windows
             // Different stations have different windows
             if ( (lctStationSeg[i] == 1 && lctStationSeg[j] == 2) || (lctStationSeg[i] == 1 && lctStationSeg[j] == 3)) {
	       if ( abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) <= 4) in_dr_eta = true; 
             }
             
	     else if ( abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) <= 6 ) in_dr_eta = true; 

             // Now the phi windows which require a vector of phi window values depending on eta and ME station combo
             // Get eta bit for this segment pair. Find asscoiated phi window
             if ( lctStationSeg[i] == 1 && lctStationSeg[j] == 2 ) {
               int dphi_window = eta_dphi_ME1toME2[eta_bit];

               if ( abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] ) <= dphi_window ) in_dr_phi = true; 
             }

	     if ( (lctStationSeg[i] == 1 && lctStationSeg[j] == 3) || (lctStationSeg[i] == 1 && lctStationSeg[j] == 4) ) {
               int dphi_window = eta_dphi_ME1toME3[eta_bit];

               if ( abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] ) <= dphi_window ) in_dr_phi = true; 
             }

             if ( (lctStationSeg[i] == 2 && lctStationSeg[j] == 3) || (lctStationSeg[i] == 2 && lctStationSeg[j] == 4) ) {

               if (abs (lctStationSeg[i]/4 - lctStationSeg[j]/4) <= 127) in_dr_phi = true; 
             }

             // Fill hist for phi cooridiante of pair that failed phi window
             if ( !in_dr_phi ) {
               hNoTrigger_phicuts -> Fill ( evt.phiReco->at(iReco) );
             }
	     
	     // Now tell me if these segments were in the windows and which ones.  Keep track of pairs that were in windows.
             // cout << "in_dr_phi = " << in_dr_phi << " and in_dr_eta = " << in_dr_eta << endl;
             if ( (in_dr_phi) && (in_dr_eta) ) {
               local_seg_pair_in_window = true;
               num_seg_pair_in_window +=1;
               if ( printLevel > 1 ) cout << "  Segment " << i << " and " << j << " are in eta and phi window " << endl;
             }
             else if ( (in_dr_phi) ) {
               how_many_not_in_dr_eta += 1;
               if ( printLevel > 1 ) cout << " Segment " << i << " and " << j << " are in phi window, but not eta " << endl;
             }
             else if ( (in_dr_eta) ) {
               how_many_not_in_dr_phi += 1;
               if ( printLevel > 1 ) cout << " Segment " << i << " and " << j << " are in eta window, but not phi " << endl;
             }
             if ( (!in_dr_eta) && (!in_dr_phi) ) {
               how_many_not_in_either += 1;
               if ( printLevel > 1 ) cout << " Segment " << i << " and " << j << " are not in eta or phi window " << endl;
             }
	     
	     // Which stations are segment i and segment j in? Keep track of total phi failures between station pairs.
             if ( !in_dr_phi ) {
               if ( (lctStationSeg[i] == 1 || lctStationSeg[j] == 1) && (lctStationSeg[i] == 2 || lctStationSeg[j] == 2 ) )  how_many_ME1ME2 += 1;
               if ( (lctStationSeg[i] == 1 || lctStationSeg[j] == 1) && (lctStationSeg[i] == 3 || lctStationSeg[j] == 3 ) )  how_many_ME1ME3 += 1;
               if ( (lctStationSeg[i] == 1 || lctStationSeg[j] == 1) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4 ) )  how_many_ME1ME4 += 1;
               if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 3 || lctStationSeg[j] == 3 ) )  how_many_ME2ME3 += 1;
               if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4 ) )  how_many_ME2ME4 += 1;
               if ( (lctStationSeg[i] == 3 || lctStationSeg[j] == 3) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4 ) )  how_many_ME3ME4 += 1;
             }
	     
	     // ------- What is seg pair mode? -------------------------------
             if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 3 || lctStationSeg[j] == 3) )      pairMode = 8;

             else if ( (lctStationSeg[i] == 2 || lctStationSeg[j] == 2) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4) ) pairMode = 9;

             else if ( (lctStationSeg[i] == 3 || lctStationSeg[j] == 3) && (lctStationSeg[i] == 4 || lctStationSeg[j] == 4) ) pairMode = 10;

             if ( printLevel > 1 ) {
               cout << "Segment " << i << " and segment " << j << " are in stations ME" << lctStationSeg[i] << " and ME" << lctStationSeg[j]  << endl;
               cout << "Segment " << i << " and segment " << j << " have an eta difference of  = " << abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) << endl;
               cout << "Segment " << i << " and segment " << j << " have a phi difference of  = " << abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] )  << endl;
             }

             // calculate dr phi/eta between all segment LCT combos ; fill histograms
             dr_eta = abs( lctEtaBitSeg[i]-lctEtaBitSeg[j]);
             dr_phi = abs( lctPhiBitSeg[i]-lctPhiBitSeg[j]);

             hdr_eta  -> Fill (dr_eta = abs( lctEtaBitSeg[i]-lctEtaBitSeg[j]) );
             hdr_phi  -> Fill (dr_phi = abs( lctPhiBitSeg[i]-lctPhiBitSeg[j]) );

	     // ----- pairs(or triplets) that made tracks to close to sector edge---------------------------------------------
             // ----- track phi is taken from lct in station 2,3,4 in that priority -----------------------------
             if ( (lctStationSeg[i] < lctStationSeg[j]) && (lctStationSeg[i] != 1) ) {
               sectorTrk_phi = lctPhiBitSeg[i];
             }
             else {sectorTrk_phi = lctPhiBitSeg[j]; }

	   
             // ------- Check to see if track failed due to sector edging ----------------------------
             if ( local_seg_pair_in_window ) {
               if (printLevel > 1 ) {
                 cout << "Segments could have made a track. Track phi is = " << sectorTrk_phi << ".  Checking to see if this is too close to sector edge " << endl;
                 cout << "Track mode = " << pairMode << ". If mode is 5, 8, 9, 10, and too close to sector edge, cancel the track. " << endl;
               }
               if ( ( sectorTrk_phi < 128 || sectorTrk_phi > (4095-128) ) &&
                    (  (pairMode == 10 || pairMode == 8 || pairMode == 9) ) ) {
                 seg_pair_in_window_but_too_close = true;
                 if (printLevel > 1) {
                   cout << "This segment pair is in windows but track was cancelled due to track being to close to sector edge(128 phi bits)" << endl;
                 }
               }
               else {
                 seg_pair_in_window = true;
                 if (printLevel > 1) cout << "This segment pair was in eta/phi windows and not too close to sector edge." << endl;
               }
               if (printLevel > 1) cout << "Distance from sector edge is = " << abs(sectorTrk_phi - 4096) << ".  Cant be < 128 or > 3968" << endl << endl;
             }
	   } // end j seg loop
	 } // end i seg loop
	 

	 // ============================================================================================================================

	 // This ends CSC seg checks. Now look at DT segs and CSC segs. Only look at CSC to DT segs.
	 bool isdt_iseg = false;
	 //bool dt_csc_in_phi_window = false;	
	 int num_dt_csc_segs_different_sector = 0;
	 int num_dt_csc_not_in_dphi  = 0;
	 int num_dt_csc_seg_pairs    = 0;
	
	 if (printLevel > 1) cout << "Now checking DT-CSC Segments " << endl;

	 for (int i=0; i<counterSegs; i++) {

	   if ( isdtlctSeg[i] == 1 ) {
	     if (printLevel >1 ) cout << " Segment " << i << " is DT " << endl;  
	     isdt_iseg = true;
	   }
	   
           for (int j=i+1; j<counterSegs; j++) {
	     if ( isdtlctSeg[j] == 1 ) {
	       if (printLevel >1 ) cout << " Segment " << j << " is DT " << endl;
	     }

	     if ( isdtlctSeg[i] == 1 && isdtlctSeg[j] == 1) continue;
	     if ( isdtlctSeg[i] == 0 && isdtlctSeg[j] == 0) continue;

	     num_dt_csc_seg_pairs += 1;
	     
	     // If youve made it this far then its time to compare CSC phi to DT phi and if both are in same sector
	     // Get phi window from LUT, which takes CSC eta bit as input
	     int eta_bit = -999;
	     int dphi    = -999;
	     
	     // Sector check
	     if ( lctSectorSeg[i] != lctSectorSeg[j] ) {
	       if (printLevel > 1) cout << " Seg " << i << " and Seg " << j << " are not in same sector." << endl;
	       num_dt_csc_segs_different_sector +=1;
	       continue;
	     }

	     if ( !isdt_iseg ) eta_bit = lctEtaBitSeg[i] >> 1;  
	     else              eta_bit = lctEtaBitSeg[j] >> 1; 
	    
	     dphi = dt_csc_dphi[eta_bit];
	     if (printLevel > 1) cout << " eta_bit = " << eta_bit << " , dphi = " << dphi << endl;
	     	    
	     if ( ( abs(lctPhiBitSeg[i] - lctPhiBitSeg[j]) ) > dphi ) {
	       if (printLevel > 1) cout << " Seg " << i << " and Seg " << j << " are not in phi window." << endl;
	       num_dt_csc_not_in_dphi += 1;	       
	     }
	     else { 
	       if(printLevel > 1) cout << " Seg " << i << " and Seg " << j << " are in phi window." << endl;	       
	       dt_csc_in_window = true;
	     }
	     	     	     
	   } // end j seg loop
	 } // end i seg loop

         // Make sure every muon not triggered is counted for one reason only.
         int fail_counter = 0;
	 
	 if ( num_csc_muonsegs < 2 && num_dt_csc_segs_different_sector > 0 && !dt_csc_in_window && num_dt_csc_not_in_dphi == 0) {
	   if (printLevel > 1) cout << "This muon failed due to DT-CSC seg pair in different sectors. " << endl;
	   num_overlap_sector_fail += 1;

	   if (is_overlap_muon) _mapc << "Fail DT-CSC different sectors: Overlap Muon";

	   fail_counter += 1;
	 }

	 if ( num_csc_muonsegs < 2 && !dt_csc_in_window && num_dt_csc_segs_different_sector == 0 ) {
	   if (printLevel > 1) cout << " This muon failed due to DT-CSC seg pair not in phi window(only one CSC seg present). " << endl;
	   num_overlap_DT_CSC_dphi_fail += 1;
	   
           if (is_overlap_muon) _mapc << "Fail DT-CSC Phi window: Overlap Muon";
	   
	   fail_counter+= 1;
	 }       	  

	 else if (dt_csc_in_window) {
	   if (printLevel > 1) cout << " This muon had DT-CSC segment pair in phi window. Should have made a track " << endl;
	   num_overlap_should_have_made += 1;

           if (is_overlap_muon) _mapc << "DT-CSC Should have made a track: Overlap Muon";

	   fail_counter+= 1;
	 }
	   
       // If DT-CSC pair was not in window, check CSC seg pairs
       if (!dt_csc_in_window && num_csc_muonsegs > 1) {
	 // Lets see if the muon was in eta/phi windows
	 if ( printLevel > 1 ) {
	   cout << "Now checking CSC seg pairs" << endl;
	   cout << "This muon had " <<   num_seg_pair_in_window   << " segment pairs in eta/phi window" << endl;
           cout << "This muon had " <<  how_many_not_in_dr_eta    << " segment pairs not in eta window " << endl;
           cout << "This muon had " <<  how_many_not_in_dr_phi    << " segment pairs not in phi window " << endl;
	   cout << "This muon had " <<  how_many_not_in_either    << " segment pairs not in phi and eta window " << endl;
	   cout << "This muon had " <<  num_segs_different_sector << " segment pairs in different sectors " << endl;
	 }

         if ( seg_pair_in_window ) {
           if ( printLevel > 1 ) cout << " This muon had a segment pair eta/phi window and not too close to sector edge. Should have made a track." << endl;
           if ( evt.SizeTrk == 0 ) num_muons_with_at_least_one_match_zerotrk += 1;
           if ( evt.SizeTrk > 0  ) num_muons_with_at_least_one_match += 1;
	   
	   if (is_overlap_muon) _mapc << "Should have made a Muon: Overlap Muon";
	   
	   fail_counter +=1;
         }
 
        // Keep track of muons that failed for just eta, just phi window, or both
	 if ( (how_many_not_in_dr_eta > 0) && (how_many_not_in_dr_phi == 0) && (how_many_not_in_either == 0) &&
	      (!seg_pair_in_window) && (!seg_pair_in_window_but_too_close) ) {
           if ( evt.SizeTrk == 0 ) num_muons_failed_deta_zerotrk += 1;
           if ( evt.SizeTrk > 0  ) num_muons_failed_deta += 1;
           if ( printLevel > 1 ) cout << " This muon failed due to eta window only " << endl;

	   if (is_overlap_muon) _mapc << "Fail Eta/Phi Windows: Overlap Muon";

           fail_counter +=1;
         }
	 
	 if ( (how_many_not_in_dr_phi > 0) && (how_many_not_in_dr_eta == 0) && (how_many_not_in_either == 0) &&
	      (!seg_pair_in_window) && (!seg_pair_in_window_but_too_close) ) {
           if ( evt.SizeTrk == 0 ) num_muons_failed_dphi_zerotrk += 1;
           if ( evt.SizeTrk > 0  ) num_muons_failed_dphi += 1;
	   if ( printLevel > 1 ) cout << " This muon failed due to phi window only " << endl;

	   if (is_overlap_muon) _mapc << "Fail Eta/Phi Windows: Overlap Muon";

           fail_counter +=1;
         }
	 
         if ( (how_many_not_in_either > 0 && !seg_pair_in_window && !seg_pair_in_window_but_too_close) ||
            (how_many_not_in_dr_phi > 0 && how_many_not_in_dr_eta > 0 && how_many_not_in_either == 0 && !seg_pair_in_window_but_too_close) ) {
	   if ( evt.SizeTrk == 0 ) num_muons_failed_deta_dphi_zerotrk += 1;
           if ( evt.SizeTrk > 0  ) num_muons_failed_deta_dphi += 1;
	   if ( printLevel > 1 ) cout << " This muon failed due to phi and eta window " << endl;

	   if (is_overlap_muon) _mapc << "Fail Eta/Phi Windows: Overlap Muon";

           fail_counter +=1;
         }
	 
	 if ( (!seg_pair_in_window) && (seg_pair_in_window_but_too_close) ) {
           if ( printLevel > 1 ) cout << " This muon failed due to segment pair in eta/phi window but is too close to sector edge." << endl;
           if ( evt.SizeTrk == 0 ) num_noTrigger_sector_edge_zerotrk += 1;
           if ( evt.SizeTrk > 0  ) num_noTrigger_sector_edge         += 1;

	   if (is_overlap_muon) _mapc << "Fail Sector Edge: Overlap Muon";

           fail_counter +=1;
         }

	 if ( num_segs_different_sector == (num_segment_pairs - num_segs_same_station) ) {
           if ( evt.SizeTrk == 0 )  num_muons_fail_different_sector_zerotrk += 1;
           if ( evt.SizeTrk > 0  )  num_muons_fail_different_sector         += 1;
           if ( printLevel > 1 ) cout << " This muon failed due to all LCTs pairs in different sectors " << endl;

	   if (is_overlap_muon) _mapc << "Fail Different Sectors: Overlap Muon";

           fail_counter +=1;
         }

         if ( event_trigger && printLevel > 1) cout << " This muon was part of an event that had a trigger on another muon " << endl;
       
	 // Throw a warning if failure counter is > 1.
	 if (fail_counter > 1)
	   cout << evt.Event << "-----> Failure Counter is more than one for this muon!!!!" << endl;
	 
	 if (fail_counter == 0)
	   cout << evt.Event << "-----> Failure Counter is zero for this muon!!!!" << endl;
	 
       } // end !dt_csc_in_window block 
       
       } // end is overlap muon block 
       
     } //end isTriggered if block
         
       // Fill the numerator for match       
     if (isTriggered) {
       if (printLevel>1) cout << " -> Numerator is filled " << endl ;
       
       //hMode->Fill(modeWinner);
       //twod_mode_eta -> Fill( modeWinner, evt.etaReco->at(iReco) ); 
       hNTriggeredMuonsVsPt  -> Fill ( evt.ptReco ->at(iReco) );
       hNTriggeredMuonsVsEta -> Fill ( evt.etaReco->at(iReco) );
       hNTriggeredMuonsVsPhi -> Fill ( evt.phiReco->at(iReco) );

       if (!isEndcap_plus) new_hNTriggeredMuonsVsPhi -> Fill ( evt.phiReco->at(iReco) );

     }
     
     // Fill the Hist for not matched
     if (!isTriggered) {
       
       hNnoTriggeredMuonsVsEta -> Fill ( evt.etaReco->at(iReco) );
       hNoTriggerPhi -> Fill ( evt.phiReco->at(iReco) );  
       twod_phi_eta -> Fill (evt.etaReco->at(iReco), evt.phiReco->at(iReco));
       
     }
     
    } // end loop on reco muons
    
  } // end Loop evts


  //------------------------------------------------------------------------- 
  cout << "----------------------------------------"
    "----------------------------------------" << endl;
  //  Set graphLevel to make hists or not --------------------------------------------------------------------------- 
  if ( graphLevel > 1 ) {
  
    // hTypeMu
    // ---------------------------------------------------------------------------
    SetStyleh1(hTypeMu, 629, 1, 2, "#mu Type");
    
    hTypeMu->GetXaxis()->SetBinLabel(1,"GL");
    hTypeMu->GetXaxis()->SetBinLabel(2,"STA");
    hTypeMu->GetXaxis()->SetBinLabel(3,"TRK");
    hTypeMu->GetXaxis()->SetBinLabel(4,"TRK+STA != GBL");
    
    TCanvas* TypeMu = new TCanvas("TypeMu", "", 0, 0, 400, 400);
    
    hTypeMu->SetMinimum(0);
    hTypeMu->Draw("histo text");
    
    if ( graphLevel > 1 ) PrintIt(TypeMu, "Reco'd Muon Type");
    
    if (isSave) {
      TypeMu -> SaveAs(png+"TypeMuCollisions-TwoSegMatched.png");
      TypeMu -> SaveAs(eps+"TypeMuCollisions-TwoSegMatched.eps");
      TypeMu -> SaveAs(rootPlot+"TypeMuCollisions-TwoSegMatched.root");
    }
    
    
    // --------------------------------------------------------------------------- 
    // hMode
    // --------------------------------------------------------------------------- 
    SetStyleh1(hMode, 629, 1, 2, "Mode");
    
    hMode->SetMinimum(0);
    
    TCanvas* Mode = new TCanvas("Mode", "", 420, 0, 400, 400);
    
    hMode -> Draw("histo text");
    
    if( graphLevel > 1)  PrintIt(Mode, "Mode");
    
    if (isSave) {
      Mode -> SaveAs(png+"GBL-ModeTriggerMu-TwoSegMatched.png");
      Mode -> SaveAs(eps+"GBL-ModeTriggerMu-TwoSegMatched.eps");
      Mode -> SaveAs(rootPlot+"GBL-ModeTriggerMu-TwoSegMatched.root");
    }
    
    
    /*  TCanvas* Pt = new TCanvas("Pt(eta = 2.2)", "", 420, 0, 400, 400);
	hPt -> Draw();
	PrintIt(Pt, "Pt(eta = 2.2)");
    */
    
    // ---------------------------------------------------------------------------
    // Efficiency Vs Pt With ClopperPearsonBinomialInterval
    // ---------------------------------------------------------------------------
    double pt[8];
    double ptErr[8];
    
    double eff[8];
    double * eefflCP = new double[8];
    double * eeffhCP = new double[8];
    
    for(int iPt = 0; iPt < 8; ++iPt) {
      
      double lowEdge = hNMuonsVsPt->GetXaxis()->GetBinLowEdge(iPt+1);
      double upEdge  = hNMuonsVsPt->GetXaxis()->GetBinUpEdge(iPt+1);
      
      pt[iPt]    = (upEdge-lowEdge)/2 + lowEdge;
      ptErr[iPt] = (upEdge-lowEdge)/2;
      
      // if (printLevel>1)
      //  cout << "pt="<<pt[iPt]<<", ptErr="<<ptErr[iPt];
      
      double num = hNTriggeredMuonsVsPt -> GetBinContent(iPt+1);
      double den = hNMuonsVsPt -> GetBinContent(iPt+1);
       
      // compute the efficiency
      if (den !=0 ) eff[iPt] = num/den;
      else          eff[iPt] = 0;
      
      // compute the error
      double  cp_lower = TEfficiency::ClopperPearson(den, num, 0.683 , 1);
      double  cp_upper = TEfficiency::ClopperPearson(den, num, 0.683 , 0);
      
      eefflCP[iPt] = eff[iPt] - cp_lower;
      eeffhCP[iPt] = cp_upper - eff[iPt];
         
    }
    
    // --- DRAW THE EFFICIENCY VS PT --- 
    gEffvsPt = new TGraphAsymmErrors(8, pt, eff, ptErr, ptErr, eefflCP, eeffhCP);
    gEffvsPt -> SetTitle("");
    gEffvsPt -> SetMinimum(0.8);
    gEffvsPt -> SetMaximum(1.02);
    
    gEffvsPt -> SetLineColor(629);
    gEffvsPt -> SetLineWidth(2);
    gEffvsPt -> SetMarkerStyle(23);
    gEffvsPt -> SetMarkerSize(0.8);
    
    gEffvsPt -> GetXaxis()-> SetTitle("p_T(#mu) [GeV/c]");
    gEffvsPt -> GetYaxis()-> SetTitle("Efficiency");
    gEffvsPt -> GetYaxis()-> SetTitleOffset(1.35);
    
    gEffvsPt -> GetXaxis()-> SetNdivisions(509);
    gEffvsPt -> GetYaxis()-> SetNdivisions(514);
    
    //gEffvsPt -> GetXaxis()->SetRangeUser(0,20);
    TCanvas* EffvsPt = new TCanvas("EffvsPt", "", 0, 500, 600, 600);
    
    EffvsPt->SetGridx();
    EffvsPt->SetGridy();
    
    gEffvsPt -> Draw("AP");
    
    if( graphLevel > 1)  PrintIt(EffvsPt, "CSCTF Efficiency Curve, "+TitleExt);
    
    if (isSave) {
      EffvsPt -> SaveAs(png+"GBL-EffvsPt-TwoSegMatched.png");
      EffvsPt -> SaveAs(eps+"GBL-EffvsPt-TwoSegMatched.eps");
      EffvsPt -> SaveAs(rootPlot+"GBL-EffvsPt-TwoSegMatched.root");
    }
    
    // ---------------------------------------------------------------------------
    // Efficiency Vs Eta
    // ---------------------------------------------------------------------------
    double eta[etaNBins];
    double etaErr[etaNBins];
    
    double effEta[etaNBins];
    double * eefflCPEta = new double[etaNBins];
    double * eeffhCPEta = new double[etaNBins];
    
    for(int iEta = 0; iEta < etaNBins; ++iEta) {
      
      double lowEdge = hNMuonsVsEta->GetXaxis()->GetBinLowEdge(iEta+1);
      double upEdge  = hNMuonsVsEta->GetXaxis()->GetBinUpEdge(iEta+1);
      
      eta[iEta]    = (upEdge-lowEdge)/2 + lowEdge;
      etaErr[iEta] = (upEdge-lowEdge)/2;
      
      double num = hNTriggeredMuonsVsEta -> GetBinContent(iEta+1);
      double den = hNMuonsVsEta -> GetBinContent(iEta+1);
      
      // compute the efficiency
      if (den !=0 ) effEta[iEta] = num/den;
      else          effEta[iEta] = 0;
      
      // compute the error
      double  cp_lower = TEfficiency::ClopperPearson(den, num, 0.683 , 1);
      double  cp_upper = TEfficiency::ClopperPearson(den, num, 0.683 , 0);
      
      eefflCPEta[iEta] = effEta[iEta] - cp_lower;
      eeffhCPEta[iEta] = cp_upper - effEta[iEta];
      
    }
    
    gEffvsEta = new TGraphAsymmErrors(etaNBins, eta, effEta, etaErr, etaErr, eefflCPEta, eeffhCPEta);
    gEffvsEta -> SetTitle("");
    gEffvsEta -> SetMinimum(0.8);
    gEffvsEta -> SetMaximum(1.02);
    
    gEffvsEta -> SetLineColor(629);
    gEffvsEta -> SetLineWidth(2);
    gEffvsEta -> SetMarkerStyle(23);
    gEffvsEta -> SetMarkerSize(0.8);
    
    gEffvsEta -> GetXaxis()-> SetTitle("|#eta(#mu)|");
    gEffvsEta -> GetYaxis()-> SetTitle("Efficiency");
    gEffvsEta -> GetYaxis()-> SetTitleOffset(1.35);
    
    gEffvsEta -> GetXaxis()-> SetNdivisions(509);
    gEffvsEta -> GetYaxis()-> SetNdivisions(514);
    
    TCanvas* EffvsEta = new TCanvas("EffvsEta", "", 420, 500, 600, 600);
    
  EffvsEta->SetGridx();
  EffvsEta->SetGridy();

  gEffvsEta -> Draw("AP");

  if( graphLevel > 1)  PrintIt(EffvsEta, "CSCTF Efficiency Curve, "+TitleExt);

  if (isSave) {
    EffvsEta -> SaveAs(png+"GBL-EffvsEta-TwoSegMatched.png");
    EffvsEta -> SaveAs(eps+"GBL-EffvsEta-TwoSegMatched.eps");
    EffvsEta -> SaveAs(rootPlot+"GBL-EffvsEta-TwoSegMatched.root");
  }

  // --------------------------------------------------------------------------- 
  // hNnoTriggeredMuonsVsEta
  // ---------------------------------------------------------------------------
  TCanvas* NoTrigger = new TCanvas("NoTriggerEta", "", 420, 500, 600, 600);
  //hNnoTriggeredMuonsVsEta -> GetXaxis() -> SetTitle("|#eta(#mu)|");
  //hNnoTriggeredMuonsVsEta -> GetYaxis() -> SetTitle("Count");
  //hNnoTriggeredMuonsVsEta -> SetStats(0);\  
  //hNnoTriggeredMuonsVsEta -> Draw();

  TCanvas* twodphieta = new TCanvas("No Trigger Eta vs Phi", "", 420, 500, 600, 600);
  twod_phi_eta -> GetXaxis() -> SetTitle(" Global Eta ");
  twod_phi_eta -> GetYaxis() -> SetTitle(" Global Phi ");
  twod_phi_eta -> GetZaxis() -> SetTitle("count");
  twod_phi_eta -> SetStats(0);

  twod_phi_eta -> Draw("cont");
  if( graphLevel > 1)  PrintIt(twodphieta, "Failed Muons Eta vs Phi");

  TCanvas* twodlcteta = new TCanvas("#LCT's vs Eta", "", 420, 500, 600, 600);
  twod_lct_eta -> GetXaxis() -> SetTitle(" LCT's ");
  twod_lct_eta -> GetYaxis() -> SetTitle(" Global Eta ");
  twod_lct_eta -> Draw("cont");
  if( graphLevel > 1)  PrintIt(twodlcteta, "#LCT's vs Eta");

  TCanvas* twodmodeeta = new TCanvas("Mode vs Eta", "", 420, 500, 600, 600);
  twod_mode_eta -> GetXaxis() -> SetTitle(" Mode ");
  twod_mode_eta -> GetYaxis() -> SetTitle(" |Eta| ");
  twod_mode_eta -> Draw("cont");
  if( graphLevel > 1)  PrintIt(twodmodeeta, "Mode vs Eta");

  TCanvas* twodsegphi = new TCanvas(" twodsegphi", "", 420, 500, 600, 600);
  twod_mode_eta -> GetXaxis() -> SetTitle(" Segment Global Phi ");
  twod_mode_eta -> GetYaxis() -> SetTitle(" Trig Prim Global Phi ");
  twod_mode_eta -> Draw("cont");
  //if( graphLevel > 1)
  PrintIt(twodsegphi, "");


  TCanvas* phiwindow = new TCanvas("Phi Window", "", 420, 500, 600, 600);
  hphiwindow -> GetYaxis() -> SetTitle(" Phi Window ");
  hphiwindow -> GetXaxis() -> SetTitle(" Eta Bit ");
  hphiwindow -> SetStats(0);
  hphiwindow -> GetYaxis() -> SetTitleOffset(1.35);
  hphiwindow -> Draw("cont");
  if( graphLevel > 1)  PrintIt(phiwindow, "Phi Window");

  TCanvas* overlapeta = new TCanvas("overlap eta", "", 420, 500, 600, 600);
  overlap_eta -> GetYaxis() -> SetTitle(" Count ");
  overlap_eta -> GetXaxis() -> SetTitle(" Global Eta ");
  overlap_eta -> Draw();
  overlap_eta -> SetStats(0);
  PrintIt(overlapeta, "Overlap Muon Eta");

  TCanvas* deltaR = new TCanvas("deltaR", "", 420, 500, 600, 600);
  hDR -> GetYaxis() -> SetTitle(" Count ");
  hDR -> GetXaxis() -> SetTitle(" Delta R ");
  hDR -> Draw();
  hDR -> SetStats(0);
  PrintIt(deltaR, "DiMuon Delta R");

  TCanvas* deltaEta = new TCanvas("deltaEta", "", 420, 500, 600, 600);
  h_dimuon_eta -> GetYaxis() -> SetTitle(" Count ");
  h_dimuon_eta -> GetXaxis() -> SetTitle(" Delta Eta ");
  h_dimuon_eta -> Draw();
  h_dimuon_eta -> SetStats(0);
  PrintIt(deltaEta, "DiMuon Delta Eta");

  TCanvas* deltaPhi = new TCanvas("deltaPhi", "", 420, 500, 600, 600);
  h_dimuon_phi -> GetYaxis() -> SetTitle(" Count ");
  h_dimuon_phi -> GetXaxis() -> SetTitle(" Delta Phi ");
  h_dimuon_phi -> Draw();
  h_dimuon_phi -> SetStats(0);
  PrintIt(deltaPhi, "DiMuon Delta Phi");

  // ---------------------------------------------------------------------------                                                                                                  
  // hNoTriggerPhi                                                                                                                                                      
  // ---------------------------------------------------------------------------                                                                                                                                                            

  TCanvas* NoTriggerPhi = new TCanvas("NoTriggerPhi", "", 420, 500, 600, 600);
  hNoTriggerPhi -> GetXaxis() -> SetTitle("#phi(#mu)");
  hNoTriggerPhi -> GetYaxis() -> SetTitle("Count");
  hNoTriggerPhi -> SetStats(0);
  hNoTriggerPhi -> Draw();

  if( graphLevel > 1)  PrintIt(NoTriggerPhi, "Global Phi of Failed Muons");

  if (isSave) {
    NoTrigger -> SaveAs(png+"NoTriggerPhi.png");
    NoTrigger -> SaveAs(eps+"NoTriggerPhi.eps");
    NoTrigger -> SaveAs(rootPlot+"NoTriggerPhi.png.root");
  }

  TCanvas* NoTrigger_phicuts = new TCanvas("NoTrigger_phicuts", "", 420, 500, 600, 600);
  hNoTrigger_phicuts -> GetXaxis() -> SetTitle("|#phi(#mu)|");
  hNoTrigger_phicuts -> GetYaxis() -> SetTitle("Count");
 
  hNoTrigger_phicuts -> Draw();

  if ( graphLevel > 1)  PrintIt(NoTrigger_phicuts, "NoTriggerPhi(not in phi window)");
 

  // ---------------------------------------------------------------------------                                                                                                 

// hdr_phi/eta
// ---------------------------------------------------------------------------                                                                                                

   TCanvas* dreta = new TCanvas("Delta Eta of Failed Events", "", 420, 500, 600, 600);
  hdr_eta -> GetXaxis() -> SetTitle("Delta Eta");
  hdr_eta -> GetYaxis() -> SetTitle("Count");

  hdr_eta -> Draw();

  if ( graphLevel > 1)  PrintIt(dreta, "Delta Eta of Failed Events(1.3<|eta|<1.7)");

  if (isSave) {
    NoTrigger -> SaveAs(png+"dr_eta.png");
    //NoTrigger -> SaveAs(eps+"dr_eta.eps");
    //NoTrigger -> SaveAs(rootPlot+"dr_eta.png.root");
  }
  
  TCanvas* drphi = new TCanvas("Delta Phi of Failed Events(1.3<|eta|<1.7)", "", 420, 500, 600, 600);
  hdr_phi -> GetXaxis() -> SetTitle("Delta Phi");
  hdr_phi -> GetYaxis() -> SetTitle("Count");

  hdr_phi -> Draw();

  if ( graphLevel > 1) PrintIt(drphi, "Delta Phi of Failed Events(1.3<|eta|<1.7)");

  if (isSave) {
    NoTrigger -> SaveAs(png+"dr_phi.png");
    //NoTrigger -> SaveAs(eps+"dr_eta.eps");
    //NoTrigger -> SaveAs(rootPlot+"dr_eta.png.root");
  }


  // ---------------------------------------------------------------------------
  // Efficiency Vs Phi
  // ---------------------------------------------------------------------------
  double phi[phiNBins];
  double phiErr[phiNBins];
   
  double effPhi[phiNBins];
  double * eefflCPPhi = new double[phiNBins];
  double * eeffhCPPhi = new double[phiNBins];
 
  for(int iPhi = 0; iPhi < phiNBins; ++iPhi) {
  
    double lowEdge = hNMuonsVsPhi->GetXaxis()->GetBinLowEdge(iPhi+1);
    double upEdge  = hNMuonsVsPhi->GetXaxis()->GetBinUpEdge(iPhi+1);

    phi[iPhi]    = (upEdge-lowEdge)/2 + lowEdge;
    phiErr[iPhi] = (upEdge-lowEdge)/2;
    
    // if (printLevel>1)
    // cout << "phi="<<phi[iPhi]<<", phiErr="<<phiErr[iPhi];
    
    double num = hNTriggeredMuonsVsPhi -> GetBinContent(iPhi+1);
    double den = hNMuonsVsPhi -> GetBinContent(iPhi+1);

    // compute the efficiency
    if (den !=0 ) effPhi[iPhi] = num/den;
    else          effPhi[iPhi] = 0;
     
    // compute the error
    double  cp_lower = TEfficiency::ClopperPearson(den, num, 0.683 , 1);
    double  cp_upper = TEfficiency::ClopperPearson(den, num, 0.683 , 0);

    eefflCPPhi[iPhi] = effPhi[iPhi] - cp_lower;
    eeffhCPPhi[iPhi] = cp_upper - effPhi[iPhi];
    
  }

  gEffvsPhi = new TGraphAsymmErrors(phiNBins, phi, effPhi, phiErr, phiErr, eefflCPPhi, eeffhCPPhi);
  gEffvsPhi -> SetTitle("");
  gEffvsPhi -> SetMinimum(0.8);
  gEffvsPhi -> SetMaximum(1.02);

  gEffvsPhi -> SetLineColor(629);
  gEffvsPhi -> SetLineWidth(2);
  gEffvsPhi -> SetMarkerStyle(23);
  gEffvsPhi -> SetMarkerSize(0.8);

  gEffvsPhi -> GetXaxis()-> SetTitle("#phi(#mu)");
  gEffvsPhi -> GetYaxis()-> SetTitle("Efficiency");
  gEffvsPhi -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsPhi -> GetXaxis()-> SetNdivisions(509);
  gEffvsPhi -> GetYaxis()-> SetNdivisions(514);

  TCanvas* EffvsPhi = new TCanvas("EffvsPhi", "", 840, 500, 600, 600);

  EffvsPhi->SetGridx();
  EffvsPhi->SetGridy();
  
  gEffvsPhi -> Draw("AP");
  
  if ( graphLevel > 1)  PrintIt(EffvsPhi, "CSCTF Efficiency Curve, "+TitleExt);

  if (isSave) {
    EffvsPhi -> SaveAs(png+"GBL-EffvsPhi-TwoSegMatched.png");
    EffvsPhi -> SaveAs(eps+"GBL-EffvsPhi-TwoSegMatched.eps");
    EffvsPhi -> SaveAs(rootPlot+"GBL-EffvsPhi-TwoSegMatched.root");
  }
 } // close graphLevel if block


  /*
  // ---------------- Create Pt Turn On Curves --------------------
  // these define the x(pt) axis location and width of x(pt) error(bin width).
  double pt_5     [pt_turn_bins];
  double pt_err_5 [pt_turn_bins];
  
  double pt_7     [pt_turn_bins];
  double pt_err_7 [pt_turn_bins];

  double pt_10    [pt_turn_bins];
  double pt_err_10[pt_turn_bins];

  // these define the y(efficiency) axis and y error as determined by counting statistics
  double pt_eff_5     [pt_turn_bins];
  double pt_eff_err_5 [pt_turn_bins];

  double pt_eff_7     [pt_turn_bins]; 
  double pt_eff_err_7 [pt_turn_bins];
  
  double pt_eff_10    [pt_turn_bins];
  double pt_eff_err_10[pt_turn_bins];

  double pass = 0;
  double fail = 0;
  double eff  = 0;

  // loop to determine turn on efficiency (y axis coordinates)
  for (int i=0; i < pt_turn_bins; i++) {
    
    // =========== 5 GeV ==================================
    pass = pt_turn_pass_5 -> GetBinContent(i);
    fail = pt_turn_fail_5 -> GetBinContent(i);
  
    if ( (pass + fail) != 0 ) eff = pass/(pass+fail);

    pt_eff_5[i] = eff;

    pt_eff_err_5[i] = sqrt( ( ( (pass+fail) / num_total_turn_on)*(1-( (pass+fail) / num_total_turn_on) ) / num_total_turn_on ) );

    //get x(pt) axis coordinates
    double lowEdge_5 = pt_turn_on_5 -> GetXaxis() -> GetBinLowEdge(i+1);
    double upEdge_5  = pt_turn_on_5 -> GetXaxis() -> GetBinUpEdge(i+1);

    pt_5[i]     = (upEdge_5-lowEdge_5)/2 + lowEdge_5;  
    pt_err_5[i] = (upEdge_5-lowEdge_5)/2;


    // =========== 7 GeV ==================================
    pass  = pt_turn_pass_7 -> GetBinContent(i);
    fail  = pt_turn_fail_7 -> GetBinContent(i);
      
    if ( (pass + fail) != 0 ) eff = pass/(pass+fail);

    pt_eff_7[i] = eff;
    
    pt_eff_err_7[i] = sqrt( ( ( (pass+fail) / num_total_turn_on)*(1-( (pass+fail) / num_total_turn_on) ) / num_total_turn_on ) );
    
    //get x(pt) axis coordinates
    double lowEdge_7 = pt_turn_on_7 -> GetXaxis() -> GetBinLowEdge(i+1);
    double upEdge_7  = pt_turn_on_7 -> GetXaxis() -> GetBinUpEdge(i+1);
    
    pt_7[i]     = (upEdge_7-lowEdge_7)/2 + lowEdge_7;
    pt_err_7[i] = (upEdge_7-lowEdge_7)/2;
  
    // =========== 10 GeV ==================================
    pass  = pt_turn_pass_10 -> GetBinContent(i);
    fail  = pt_turn_fail_10 -> GetBinContent(i);
  
    if ( (pass + fail) != 0 ) eff = pass/(pass+fail);

    pt_eff_10[i] = eff;

    pt_eff_err_10[i] = sqrt( ( ( (pass+fail) / num_total_turn_on)*(1-( (pass+fail) / num_total_turn_on) ) / num_total_turn_on ) );

    //get x(pt) axis coordinates
    double lowEdge_10 = pt_turn_on_10 -> GetXaxis() -> GetBinLowEdge(i+1);
    double upEdge_10  = pt_turn_on_10 -> GetXaxis() -> GetBinUpEdge(i+1);

    pt_10[i]     = (upEdge_10-lowEdge_10)/2 + lowEdge_10;
    pt_err_10[i] = (upEdge_10-lowEdge_10)/2;

    // Print out plot coordinates
    //cout << "5 GeV = "  << pt_5[i]  << " , " << pt_eff_5[i] << endl;
    //cout << "7 GeV = "  << pt_7[i]  << " , " << pt_eff_7[i] << endl;
    //cout << "10 GeV = " << pt_10[i] << " , " << pt_eff_10[i] << endl << endl;
  
  }
  
  // Create error function fits to turn on curves
  // The error function fit we wish to use -- courtesy of Bobby. Thanks Bobby.
  // TF1( NAME, FIT, LOWER BOUND, UPPER BOUND)
  //TF1 *fit = new TF1("fit", "(0.5*TMath::Erf((x/[0] + 1.0)/(TMath::Sqrt(2.0)*[1]))+0.5*TMath::Erf((x/[0] - 1.0)/(TMath::Sqrt(2.0)*[1])))*([2] + [3]*x)", 0, 100);

  // NEw fit courtesy of Ivan
  TF1 *fit5  = new TF1("fit5",  "[0]* (1 + TMath::Erf( (x - [1])*[2] ))/2.0", 0, 100);  
  TF1 *fit7  = new TF1("fit7",  "[0]* (1 + TMath::Erf( (x - [1])*[2] ))/2.0", 0, 100);
  TF1 *fit10 = new TF1("fit10", "[0]* (1 + TMath::Erf( (x - [1])*[2] ))/2.0", 0, 100);

  // The error function parameters.
  // [0] = x-value at 50% threshold (starting point should be 50% point)
  // [1] = the resolution( the range should be where x valye goes from 5% to 95% efficiency)
  // [2] = plateau efficiency
  // [3] = controls the dip in plateau efficiency for large x
  
  // The fit needs reasonable ranges for the parameters.

  //These are good for red
  fit7 -> SetParameter(0,1.0);
  fit7 -> SetParameter(1,5.0);
  fit7 -> SetParameter(2,0.5);

  fit7 -> SetParLimits(0,1,10);
  fit7 -> SetParLimits(1,1,10);
  fit7 -> SetParLimits(2,0,10);

  fit5 -> SetParameter(0,1.0);
  fit5 -> SetParameter(1,5.0);
  fit5 -> SetParameter(2,0.5);

  fit5 -> SetParLimits(0,1,10);
  fit5 -> SetParLimits(1,1,10);
  fit5 -> SetParLimits(2,0,10);

  fit10 -> SetParameter(0,1.0);
  fit10 -> SetParameter(1,5.0);
  fit10 -> SetParameter(2,0.5);

  fit10 -> SetParLimits(0,1,10);
  fit10 -> SetParLimits(1,1,10);
  fit10 -> SetParLimits(2,0,10);

  
  fit->SetParLimits(0, 3, 7);
  fit->SetParLimits(1, 0.02, 7);
  fit->SetParLimits(2, 0.5, 1);
  fit->SetParLimits(3, 1e-7,1e-4);
   

  /*

  TCanvas* turn_on = new TCanvas("turn_on", "", 840, 500, 600, 600);
  turn_on -> SetLogx();
  
  // TGraphErrors(# bins, x vector, y vector, x error vector , y error vector)
  TGraphErrors *turn_err_7 = new TGraphErrors(pt_turn_bins, pt_7, pt_eff_7, pt_err_7, pt_eff_err_7);
  TGraph *fit_graph_7      = new TGraph(pt_turn_bins, pt_7, pt_eff_7);
  fit_graph_7 -> Fit("fit7", "r");
  turn_err_7 -> SetMinimum(0);
  turn_err_7 -> GetXaxis()-> SetTitle("Pt(GeV)");
  turn_err_7 -> GetYaxis()-> SetTitle("Efficiency");
  //turn_err_7 -> GetYaxis()-> SetTitleOffset(1.35);
  turn_err_7 -> SetTitle("CSCTF Pt Turn-On Curve: 2.1<|eta|");
  turn_err_7 -> SetMarkerColor(kRed);
  turn_err_7 -> SetLineColor(kRed);
  turn_err_7 -> SetMarkerStyle(21);
  turn_err_7 -> SetMarkerSize(0.6);
  turn_err_7 -> SetLineWidth(2);
  turn_err_7  -> Draw("AP");
  fit_graph_7 -> Draw("Psame");
  
  
  //TCanvas* turn_on10 = new TCanvas("turn_on10", "", 840, 500, 600, 600);
  //turn_on10 -> SetLogx();
  TGraphErrors *turn_err_10 = new TGraphErrors(pt_turn_bins, pt_10, pt_eff_10, pt_err_10, pt_eff_err_10);
  TGraph *fit_graph_10      = new TGraph(pt_turn_bins, pt_10, pt_eff_10);
  fit_graph_10 -> Fit("fit10", "r");
  fit_graph_10 -> GetFunction("fit10")->SetLineColor(kBlue);
  turn_err_10 -> SetMarkerColor(kBlue);
  turn_err_10 -> SetLineColor(kBlue);
  turn_err_10 -> SetMarkerStyle(21);
  turn_err_10 -> SetMarkerSize(0.6);
  turn_err_10 -> SetLineWidth(2);
  turn_err_10 -> Draw("Psame");
  //turn_err_10 -> Draw("AP");
  fit_graph_10 -> Draw("Psame");

  
  //TCanvas* turn_on5 = new TCanvas("turn_on5", "", 840, 500, 600, 600);
  //turn_on5 -> SetLogx();
  TGraphErrors *turn_err_5 = new TGraphErrors(pt_turn_bins, pt_5, pt_eff_5, pt_err_5, pt_eff_err_5);
  TGraph *fit_graph_5      = new TGraph(pt_turn_bins, pt_5, pt_eff_5);
  fit_graph_5 -> Fit("fit5", "r");  
  fit_graph_5 -> GetFunction("fit5")->SetLineColor(kOrange+3);
  turn_err_5 -> SetMarkerColor(kOrange+3);
  turn_err_5 -> SetLineColor(kOrange+3);
  turn_err_5 -> SetMarkerStyle(21);
  turn_err_5 -> SetMarkerSize(0.6);
  turn_err_5 -> SetLineWidth(2);
  turn_err_5 -> Draw("Psame");
  //turn_err_5 -> Draw("AP");
  fit_graph_5 -> Draw("Psame");  

  TLine *line_10 = new TLine(10,0,10,1.1);
  line_10 -> SetLineWidth(2);
  line_10 -> SetLineStyle(kDashed);
  line_10 -> SetLineColor(kBlue);
  line_10 -> Draw("same");

  TLine *line_7 = new TLine(7,0,7,1.1);
  line_7 -> SetLineWidth(2);
  line_7 -> SetLineStyle(kDashed);
  line_7 -> SetLineColor(kRed);
  line_7 -> Draw("same");

  TLine *line_5 = new TLine(5,0,5,1.1);
  line_5 -> SetLineWidth(2);
  line_5 -> SetLineStyle(kDashed);
  line_5 -> SetLineColor(kOrange+3);
  line_5 -> Draw("same");

  //turn_on -> BuildLegend();
  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  leg -> SetFillStyle(0);
  leg -> AddEntry(turn_err_5, "5 GeV",  "l");
  leg -> AddEntry(turn_err_7, "7 GeV",  "l");
  leg -> AddEntry(turn_err_10,"10 GeV", "l");
  leg -> Draw("same");
  
  turn_on -> Update();
  

  // ------------------- Create Mode vs. Eta plots ---------------------

  TCanvas *mode_eta= new TCanvas("mode_eta","",700,500);
  mode_eta -> SetTitle("CSCTF Mode(eta)");
  hmode14 -> GetXaxis()-> SetTitle("|Eta|");
  hmode14 -> GetYaxis()-> SetTitle("Count");
  hmode14 -> SetStats(0);
 
  hmode1  -> SetFillColor(kBlue);
  hmode2  -> SetFillColor(kYellow);
  hmode3  -> SetFillColor(kBlack); 
  hmode4  -> SetFillColor(kTeal);
  hmode5  -> SetFillColor(kGreen);
  hmode6  -> SetFillColor(kRed);
  hmode7  -> SetFillColor(kViolet);
  hmode8  -> SetFillColor(kOrange);
  hmode9  -> SetFillColor(kGray);
  hmode10 -> SetFillColor(44);
  hmode11 -> SetFillColor(32);
  hmode12 -> SetFillColor(40);
  hmode13 -> SetFillColor(38);
  hmode14 -> SetFillColor(41);

  //Make each higher mode add the previous modes to it
  for (int i=0; i < 50; i++) {
    hmode2  -> SetBinContent(i, hmode1->GetBinContent(i)  + hmode2->GetBinContent(i) );
    hmode3  -> SetBinContent(i, hmode2->GetBinContent(i)  + hmode3->GetBinContent(i) );
    hmode4  -> SetBinContent(i, hmode3->GetBinContent(i)  + hmode4->GetBinContent(i) );  
    hmode5  -> SetBinContent(i, hmode4->GetBinContent(i)  + hmode5->GetBinContent(i) );
    hmode6  -> SetBinContent(i, hmode5->GetBinContent(i)  + hmode6->GetBinContent(i) );
    hmode7  -> SetBinContent(i, hmode6->GetBinContent(i)  + hmode7->GetBinContent(i) );
    hmode8  -> SetBinContent(i, hmode7->GetBinContent(i)  + hmode8->GetBinContent(i) );
    hmode9  -> SetBinContent(i, hmode8->GetBinContent(i)  + hmode9->GetBinContent(i) );
    hmode10 -> SetBinContent(i, hmode9->GetBinContent(i)  + hmode10->GetBinContent(i) );  
    hmode11 -> SetBinContent(i, hmode10->GetBinContent(i) + hmode11->GetBinContent(i) );
    hmode12 -> SetBinContent(i, hmode11->GetBinContent(i) + hmode12->GetBinContent(i) );
    hmode13 -> SetBinContent(i, hmode12->GetBinContent(i) + hmode13->GetBinContent(i) );
    hmode14 -> SetBinContent(i, hmode13->GetBinContent(i) + hmode14->GetBinContent(i) );
}

  hmode14 -> Draw();
  hmode13 -> Draw("same");
  hmode12 -> Draw("same");
  hmode11 -> Draw("same");
  hmode10 -> Draw("same");
  hmode9  -> Draw("same");
  hmode8  -> Draw("same");
  hmode7  -> Draw("same");
  hmode6  -> Draw("same");
  hmode5  -> Draw("same");
  hmode4  -> Draw("same");
  hmode3  -> Draw("same");
  hmode2  -> Draw("same");
  hmode1  -> Draw("same");
  
  mode_eta -> BuildLegend();
  PrintIt(mode_eta, "CSCTF Mode(eta)");
  */
  //TCanvas *mdoe_eta14 = new TCanvas("mode_eta14","",700,500);
  //hmode14 -> SetFillColor(kRed);
  //hmode14 -> Draw("same");
 
/* 
  TCanvas* NoTrigger = new TCanvas("NoTriggerEta", "", 420, 500, 600, 600);
  hNnoTriggeredMuonsVsEta -> SetStats(0);
  SetStyleh1(hNnoTriggeredMuonsVsEta,
             2,
             0,
	     1,
             "|#eta(#mu)|",
             "Count",
             1,
             3);

  hNnoTriggeredMuonsVsEta -> Draw();
 
  // Set Title
  //TString NoTriggerEta;
  PrintIt(NoTrigger, "NoTriggerEta");
  */



/*TCanvas* cscphibit = new TCanvas("cscphibit", "", 420, 500, 600, 600);
  CSCphibit -> GetXaxis() -> SetTitle(" CSC Global Phi");
  CSCphibit -> GetYaxis() -> SetTitle(" CSC Phi Bit ");
  CSCphibit -> GetYaxis()-> SetTitleOffset(1.50);
  //twod_seg_phi -> SetStats(0);

  CSCphibit -> Draw("cont");
  PrintIt(cscphibit, "CSC Global Phi vs. Phi Bit");
  */  

  /*
  TCanvas* twodsegphi = new TCanvas(" twodsegphi", "", 420, 500, 600, 600);
  twod_seg_phi -> GetXaxis() -> SetTitle(" Segment Global Phi ");
  twod_seg_phi -> GetYaxis() -> SetTitle(" Trig Prim Global Phi ");
  twod_seg_phi -> GetYaxis()-> SetTitleOffset(1.50);
  twod_seg_phi -> SetStats(0);
 
  twod_seg_phi -> Draw("cont");
  PrintIt(twodsegphi, "Segment Phi vs. Lct Phi");
  */
  /*TCanvas* twodsegeta = new TCanvas(" twodsegeta", "", 420, 500, 600, 600);
  twod_seg_eta -> GetXaxis() -> SetTitle(" Segment Global Eta ");
  twod_seg_eta-> GetYaxis() -> SetTitle(" Trig Prim Global Eta ");
  twod_seg_eta -> GetYaxis()-> SetTitleOffset(1.50);
  twod_seg_eta -> SetStats(0);

  twod_seg_eta -> Draw("cont");
  PrintIt(twodsegeta, "Segment Eta vs. Lct Eta");
  


  TCanvas* deltasegphi = new TCanvas("deltasegphi", "", 420, 500, 600, 600);
  delta_seg_phi -> GetXaxis() -> SetTitle("Global Phi");
  delta_seg_phi -> GetYaxis() -> SetTitle("Count");
  delta_seg_phi -> GetYaxis()-> SetTitleOffset(1.50);
  //delta_seg_phi -> SetStats(0);
  delta_seg_phi -> Draw();

  PrintIt(deltasegphi, "Segment Phi - TrigPrimitive Phi");
  */
  /*
  TCanvas* deltasegeta = new TCanvas("deltasegeta", "", 420, 500, 600, 600);
  delta_seg_eta -> GetXaxis() -> SetTitle("Global Eta");
  delta_seg_eta -> GetYaxis() -> SetTitle("Count");
  delta_seg_eta -> GetYaxis()-> SetTitleOffset(1.50);
  delta_seg_eta -> SetStats(0);
  delta_seg_eta -> Draw();

  PrintIt(deltasegeta, "Segment Eta - TrigPrimitive Eta");
*/   
  

  /*
  // make new phi plot
  double phi[phi_bin];
  double phiErr[phi_bin];
  
  double effPhi[phi_bin];
  double * eefflCPPhi = new double[phi_bin];
  double * eeffhCPPhi = new double[phi_bin];
  
  for(int iPhi = 0; iPhi < phi_bin; ++iPhi) {
    
    double lowEdge = new_hNMuonsVsPhi->GetXaxis()->GetBinLowEdge(iPhi+1);
    double upEdge  = new_hNMuonsVsPhi->GetXaxis()->GetBinUpEdge(iPhi+1);
    
    phi[iPhi]    = (upEdge-lowEdge)/2 + lowEdge;
    phiErr[iPhi] = (upEdge-lowEdge)/2;
    
    cout << "Bin # " << iPhi << " = " << phi[iPhi] << endl; 
    
    double num = new_hNTriggeredMuonsVsPhi -> GetBinContent(iPhi+1);
    double den = new_hNMuonsVsPhi          -> GetBinContent(iPhi+1);
    
    // compute the efficiency
    if (den !=0 ) effPhi[iPhi] = num/den;
    else          effPhi[iPhi] = 0;
    
    // compute the error
    double  cp_lower = TEfficiency::ClopperPearson(den, num, 0.683 , 1);
    double  cp_upper = TEfficiency::ClopperPearson(den, num, 0.683 , 0);
    
    eefflCPPhi[iPhi] = effPhi[iPhi] - cp_lower;
    eeffhCPPhi[iPhi] = cp_upper - effPhi[iPhi];
    
    cout << "Eff Phi bin " << iPhi << " = " << effPhi[iPhi] << endl;
    cout << " Numerator = " << num << " , Denominator = " << den << endl;
    
    
  }
  
  gEffvsPhi = new TGraphAsymmErrors(50, phi, effPhi, phiErr, phiErr, eefflCPPhi, eeffhCPPhi);
  gEffvsPhi -> SetTitle("");
  gEffvsPhi -> SetMinimum(0.8);
  gEffvsPhi -> SetMaximum(1.02);

  gEffvsPhi -> SetLineColor(629);
  gEffvsPhi -> SetLineWidth(2);
  gEffvsPhi -> SetMarkerStyle(23);
  gEffvsPhi -> SetMarkerSize(0.8);

  gEffvsPhi -> GetXaxis()-> SetTitle("#phi(#mu)");
  gEffvsPhi -> GetYaxis()-> SetTitle("Efficiency");
  gEffvsPhi -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsPhi -> GetXaxis()-> SetNdivisions(509);
  gEffvsPhi -> GetYaxis()-> SetNdivisions(514);

  TCanvas* EffvsPhi = new TCanvas("EffvsPhi", "", 840, 500, 600, 600);

  EffvsPhi->SetGridx();
  EffvsPhi->SetGridy();

  gEffvsPhi -> Draw("AP");
  PrintIt(EffvsPhi, "CSCTF Efficiency Curve, "+TitleExt);
  */


  TCanvas* c10 = new TCanvas("c10", "", 420, 500, 600, 600);
  hdphi_csc -> Draw();
  
  



  // ---------------------------------------------------------------------------
  // printout
  // ---------------------------------------------------------------------------


  cout << " ------------------------- End Analysis ----------------------------------------------------------" << endl;
  /*
  cout << " lctable count = " << lctable << endl;
  cout << " matched count  = " << matched << endl;
  cout << " segments count = " << segments << endl;
  cout << " debug muons = " << debug_muons << endl;
  cout << " debug matched = " << debug_matched << endl;
  cout << "debug eta = " << debug_eta << endl;
  */  
  //cout << "Number of Global muons in all the events        = " << nMuons_total << endl;
  //cout << "Number of events with more than one global muon = " << num_evt_multi_muon << endl;

  cout << "Number of Global muons matched to two LCTs      = " << nMuons << endl;
  cout << "   of these " << num_overlap_muons << " were overlap muons. " << endl; 
  cout << "Number of times CSCTF did not trigger on a muon  = "<< nMuonsNoTrigger << endl;
  //cout << "   Of these " << nOverlapMuonsNoTrigger << " were overlap muons. " << endl;
  cout << "Number of times CSCTF did not produce a track(failed muons) = " << num_SizeTrk_is_Zero << endl << endl;
       
/* cout << "Number of times failed muon is from an event with trigger on another muon = " << failed_muons_with_event_trigger << endl;
   cout << " Single Muon fail = " << num_single_muon_event << endl;
   cout << " Di Muon fail     = " << num_multi_muon_event << endl;
   cout << "  Number of failed single muon events when CSCTF makes a track = " << num_single_muon_event_withTrk << endl;
   cout << "  Number of failed multi muon events when CSCTF makes a track = " << num_multi_muon_event_withTrk << endl;
   cout << "Number of times a pair of segments failed phi extrapolation in ME1+ME2 = " << how_many_ME1ME2 << endl;
   cout << "                                                               ME1+ME3 = " << how_many_ME1ME3 << endl;
   cout << "                                                               ME1+ME4 = " << how_many_ME1ME4 << endl;
   cout << "                                                               ME2+ME3 = " << how_many_ME2ME3 << endl;
   cout << "                                                               ME2+ME4 = " << how_many_ME2ME4 << endl;
   cout << "                                                               ME3+ME4 = " << how_many_ME3ME4 << endl << endl;
   
   cout << "Number of failed muons due to DT-CSC being in wrong station = " << num_overlap_fail_station << endl;
   cout << "Number of failed muons due to DT-CSC Phi windows = " << num_overlap_DT_CSC_dphi_fail << endl;
   cout << "Number of failed muons due to DT-CSC being in wrong sector = " << num_overlap_sector_fail << endl;
   cout << "Number of failed overlap muons that should have made a track = " << num_overlap_should_have_made << endl;
   
   cout << "Number of failed non-overlap muons due to all CSC lcts being in one station = " << num_muons_all_lcts_one_station << endl;
   cout << "Number of failed muons(at least one track) with at least one pair of segments in eta/phi windows and not too close to sector edge is " << num_muons_with_at_least_one_match << endl;
   cout << "Number of muons(at least one track) that failed due to phi window only = " << num_muons_failed_dphi << endl;
   cout << "Number of muons(at least one track) that failed due to eta window only = " << num_muons_failed_deta << endl;
   cout << "Number of muons(at least one track) that failed due to eta and phi window = " << num_muons_failed_deta_dphi << endl;
   cout << "Number of muons(at least one track) with all lct pairs in different sectors = " << num_muons_fail_different_sector << endl;
   cout << "Number of muons(at least one track) that had LCT pair in eta/phi windows but was cancelled due to sector edge to close = " << num_noTrigger_sector_edge << endl;
   cout << "Number of failed muons(zero tracks) with at least one pair of segments in eta/phi windows and not too close to sector edge is " << num_muons_with_at_least_one_match_zerotrk << endl;
       cout << "Number of muons(zero track) that failed due to phi window only = " << num_muons_failed_dphi_zerotrk << endl;
       cout << "Number of muons(zero track) that failed due to eta window only = " << num_muons_failed_deta_zerotrk << endl;      
       cout << "Number of muons(zero track) that failed due to eta and phi window = " << num_muons_failed_deta_dphi_zerotrk << endl;
       cout << "Number of muons(zero track) that had LCT pair in eta/phi windows but was cancelled due to sector edge too close = " << num_noTrigger_sector_edge_zerotrk << endl;
       cout << "Number of muons(zero track) that failed due to all lct pairs in different sectors = " << num_muons_fail_different_sector_zerotrk << endl;
      */

       cout << "Total Failures: " << endl;
       cout << " Eta/Phi Windows = " << num_muons_failed_dphi + num_muons_failed_dphi_zerotrk + num_muons_failed_deta_zerotrk + num_muons_failed_deta_dphi_zerotrk +
	                               num_muons_failed_deta + num_muons_failed_deta_dphi + num_overlap_DT_CSC_dphi_fail << endl;
       cout << " LCTs Different Sector = " << num_muons_fail_different_sector + num_muons_fail_different_sector_zerotrk + num_overlap_sector_fail << endl;
       cout << " Sector Edging = " << num_noTrigger_sector_edge_zerotrk + num_noTrigger_sector_edge << endl;
       cout << " Unaccounted Failures = " << num_muons_with_at_least_one_match_zerotrk + num_muons_with_at_least_one_match << endl;
       cout << " Sum of failures = " << num_muons_failed_dphi + num_muons_failed_dphi_zerotrk + num_muons_failed_deta_zerotrk + num_muons_failed_deta_dphi_zerotrk +
                                        num_muons_failed_deta + num_muons_failed_deta_dphi + num_overlap_DT_CSC_dphi_fail +
                                        num_muons_fail_different_sector + num_muons_fail_different_sector_zerotrk +
                                        num_noTrigger_sector_edge_zerotrk + num_noTrigger_sector_edge + num_overlap_sector_fail +
                                	num_muons_with_at_least_one_match_zerotrk + num_muons_with_at_least_one_match + num_overlap_should_have_made << endl;
       
       
       cout << "Map Counters: Different Eta regions " << endl;
       cout << _mapc << endl << endl;
  


} // end function

double DR(double diffeta,
          double diffphi) {

  double diffetaSquare =diffeta*diffeta;
  double diffphiSquare =diffphi*diffphi;

  return sqrt(diffetaSquare+diffphiSquare);
}

double computeError(double num, double den) {

  if (den==0) return 0;
  double efficiency = num/den;

  return sqrt( ((1-efficiency) * efficiency)/den );

}

