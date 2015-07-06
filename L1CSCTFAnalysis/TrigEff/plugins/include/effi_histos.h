#ifndef __EFFI_HISTOS__
#define __EFFI_HISTOS__ (1)

#include <TH1.h>
#include "TH1F.h"
#include <TH2.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

// ------------------------------------------------------------------------------
// which kind of muon?
// ------------------------------------------------------------------------------
TH1F* hTypeMu        = new TH1F("hTypeMu","", 4,    0,   4); 
TH1F* hTriggerCounts = new TH1F("hTriggerCounts","",10,0,10);

// ------------------------------------------------------------------------------
// Phi/Eta Resolutions
// ------------------------------------------------------------------------------
TH1F* h_dimuon_eta = new TH1F("h_dimuon_eta", "", 100,  -3.5,  3.5);
TH1F* h_dimuon_phi = new TH1F("h_dimuon_phi", "", 100,  0, Pi);
TH1F* hDeltaPhi = new TH1F("hDeltaPhi", "", 100,    Pi,   Pi);
TH1F* hDeltaEta = new TH1F("hDeltaEta", "", 100,  -1.0,  1.0);
TH1F *hDR       = new TH1F("hDR"    , "",   100,     0,  10);
TH1F *hMode     = new TH1F("hMode",     "",  15,     0,   15);
TH1F* hNnoTriggeredMuonsVsEta = new TH1F("hNnoTriggeredMuonsVsEta","",100,0.5, 3.0);
TH1F* hNoTriggerPhi = new TH1F("hNoTriggerPhi", "", 100,    Pi,   Pi);
TH1F* hIntegralPhi = new TH1F("hIntegralPhi", "", 100,    Pi,   Pi);
TH1F* hPt = new TH1F("hPt", "", 100, 0, 20);

TH2F* CSCphibit = new TH2F("CSCphibit", "", 100,  -Pi,   Pi, 10000, 0, 4096);
TH2F* twod_phi_eta = new TH2F("2d_phi_eta", "", 1000, 0, 3, 1000, Pi, Pi);
TH2F* twod_lct_eta = new TH2F("2d_lct_eta", "", 100, 1, 5, 100, 0.8, 2.5);
TH2F* twod_mode_eta = new TH2F("2d_mode_eta", "", 15, 0, 15, 100, 0, 3);
TH2F* hphiwindow = new TH2F("hphiwindow", "", 100,  0.8,  2.5, 140 , 0 , 140);
TH2F* twod_seg_phi = new TH2F("2d_seg_phi", "", 100, -Pi, Pi, 100, -Pi, Pi);
TH2F* twod_seg_eta = new TH2F("2d_seg_eta", "", 100, -2, 1.5, 100, -2, 1.5);  

// ------------------------------------------------------------------------------
// Histograms for efficiencies
// ------------------------------------------------------------------------------
TH1F* hNTriggeredMuonsVsPt = new TH1F("hNTriggeredMuonsVsPt","",8,scalePt);
TH1F* hNMuonsVsPt = new TH1F("hNMuonsVsPt","",8,scalePt);
TH1F* overlap_eta = new TH1F("overlap_eta", "", 100,  -1.5,  1.5);
TH1F* hNTriggeredMuonsVsEta = new TH1F("hNTriggeredMuonsVsEta","",etaNBins,scaleEta);
TH1F* hNMuonsVsEta = new TH1F("hNMuonsVsEta","",etaNBins,scaleEta);
TH1F* delta_seg_phi = new TH1F("delta_seg_phi", "", 100, -0.5, 0.5);
TH1F* delta_seg_eta = new TH1F("delta_seg_eta", "", 100, -0.4, 0.4);
TH1F* hNTriggeredMuonsVsPhi = new TH1F("hNTriggeredMuonsVsPhi","",phiNBins,scalePhi);
TH1F* hNMuonsVsPhi = new TH1F("hNMuonsVsPhi","",phiNBins,scalePhi);
TH1F* hdr_eta = new TH1F("hdr_eta", "", 100, -1, 8);
TH1F* hdr_phi = new TH1F("hdr_phi", "", 100, -1, 1000);
TH1F* hNoTrigger_phicuts = new TH1F("hNoTrigger_phicuts", "", 100,  Pi,  Pi);

// ------------------------------------------------------------------------------
// Efficiencies Graphs
// ------------------------------------------------------------------------------
TGraphAsymmErrors* gEffvsPt;
TGraphAsymmErrors* gEffvsEta;
TGraphAsymmErrors* gEffvsPhi;

// ------------------------------------------------------------------------------
// Pt Turn On Curves
// ------------------------------------------------------------------------------
int pt_turn_bins = 17;
Double_t scale_turn[18] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 25, 35, 39.5, 50, 75, 100};

TH1F* pt_turn_pass_5  = new TH1F("pt_turn_pass_5",  "", pt_turn_bins, scale_turn);
TH1F* pt_turn_pass_7  = new TH1F("pt_turn_pass_7",  "", pt_turn_bins, scale_turn);
TH1F* pt_turn_pass_10 = new TH1F("pt_turn_pass_10", "", pt_turn_bins, scale_turn);


TH1F* pt_turn_fail_5  = new TH1F("pt_turn_fail_5",  "", pt_turn_bins, scale_turn);
TH1F* pt_turn_fail_7  = new TH1F("pt_turn_fail_7",  "", pt_turn_bins, scale_turn);
TH1F* pt_turn_fail_10 = new TH1F("pt_turn_fail_10", "", pt_turn_bins, scale_turn);

TH1F* pt_turn_on_5    = new TH1F("pt_turn_on_5",    "", pt_turn_bins, scale_turn);
TH1F* pt_turn_on_7    = new TH1F("pt_turn_on_7",    "", pt_turn_bins, scale_turn);
TH1F* pt_turn_on_10   = new TH1F("pt_turn_on_10",   "", pt_turn_bins, scale_turn);



// Mode plots
int eta_bins = 64;
TH1F* hmode_eta = new TH1F("hmode_eta", "", eta_bins, 0.8, 2.8);
TH1F* hmode1    = new TH1F("hmode1", "", eta_bins, 0.8, 2.8);
TH1F* hmode2    = new TH1F("hmode2", "", eta_bins, 0.8, 2.8);
TH1F* hmode3    = new TH1F("hmode3", "", eta_bins, 0.8, 2.8);
TH1F* hmode4    = new TH1F("hmode4", "", eta_bins, 0.8, 2.8);
TH1F* hmode5    = new TH1F("hmode5", "", eta_bins, 0.8, 2.8);
TH1F* hmode6    = new TH1F("hmode6", "", eta_bins, 0.8, 2.8);
TH1F* hmode7    = new TH1F("hmode7", "", eta_bins, 0.8, 2.8);
TH1F* hmode8    = new TH1F("hmode8", "", eta_bins, 0.8, 2.8);
TH1F* hmode9    = new TH1F("hmode9", "", eta_bins, 0.8, 2.8);
TH1F* hmode10   = new TH1F("hmode10", "", eta_bins, 0.8, 2.8);
TH1F* hmode11   = new TH1F("hmode11", "", eta_bins, 0.8, 2.8);
TH1F* hmode12   = new TH1F("hmode12", "", eta_bins, 0.8, 2.8);
TH1F* hmode13   = new TH1F("hmode13", "", eta_bins, 0.8, 2.8);
TH1F* hmode14   = new TH1F("hmode14", "", eta_bins, 0.8, 2.8);


#endif
