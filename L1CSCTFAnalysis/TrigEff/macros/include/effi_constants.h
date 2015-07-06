#ifndef __EFFI_CONSTANTS__
#define __EFFI_CONSTANTS__ (1)

#include "TString.h"
#include "TMath.h"

// ------------------------------------------------------------------------------
// variables
// ------------------------------------------------------------------------------
const double PhiStep = (62*TMath::Pi()/180)/4096;
const double EtaStep = 1.6/128; //(2.5-0.9)/(2^7)

const double Pi = TMath::Pi();

// ------------------------------------------------------------------------------
// switches
// ------------------------------------------------------------------------------

const int MAX_MUONS = 100; 
const int MAX_CSC_RECHIT = 48;
const int MAX_TRK_SEGS = 100;
const int MAX_CSCTF_TRK = 36;
const int MAX_LCTS_PER_TRK = 4;
const int MAX_SEGS_STD=16; 

// ===========================================================================
TString eps = "eps/";
TString png = "png/";
TString rootPlot = "rootPlot/";

// --------------------------------------------------------------------------- 
// All global muons
TString SaveExt("GblMu");
TString TitleExt("Gbl #mu");
// ===========================================================================

// --------------------------------------------------------------------------- 
// Efficiency Variables Phi
// --------------------------------------------------------------------------- 
// variable bin size Phi 
const int phiNBins = 8;
Double_t scalePhi[phiNBins+1];

// --------------------------------------------------------------------------
// Array of delta phi depending on eta.  v = [ eta , delta phi ]

int eta_dphi_ME1toME2[64] = {127, 127, 127, 127, 57, 45, 41,  42,  42,  31, 
			     127, 127,  29,  28, 29, 30, 35,  37,  34,  34,
			      36,  37,  37,  37, 39, 40, 52, 126, 104, 104,
			      87,  90,  93,  87, 85, 82, 80,  79,  82,  79,
			      79,  77,  75,  75, 67, 67, 69,  68,  67,  65, 
			      65,  64,  60,  58, 57, 57, 57,  55,  53,  51, 
			      49,  46,  36, 127};

int eta_dphi_ME1toME3[64] = {127, 127, 127, 127, 127, 127, 40, 80, 80, 64, 
			     127, 127,  62,  41,  41,  45, 47, 48, 47, 46, 
			      47,  50,  52,  51,  53,  54, 55, 73, 82, 91, 
			      91,  94, 100,  99,  95,  94, 95, 91, 96, 96, 
			      94,  94,  88,  88,  80,  80, 84, 84, 79, 78, 
			      80,  78,  75,  72,  70,  71, 69, 71, 71, 66, 
			      61,  60,  43, 127};

int dt_csc_dphi[64]       = {127, 127, 127, 127,  90,  78,  76,  76,  66,  65, 
			      59,  90,  50,  49,  37, 127, 127, 127, 127, 127, 
			     127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
			     127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
			     127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
			     127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 
			     127, 127, 127, 127};


  // variable bin size Eta 
const int etaNBins = 8;
Double_t scaleEta[etaNBins+1];
  
const  Double_t scalePt[9] = {0, 2, 3, 4, 8, 15, 25, 35, 39.5};

static void initScaleValues(){
  double phiMin = -TMath::Pi();

  scalePhi[0]=phiMin;

  for (int iPhi=1;iPhi<phiNBins+1;iPhi++)
    scalePhi[iPhi] = phiMin + (TMath::TwoPi()*iPhi/phiNBins); 

  double etaMin = 0.9;
  scaleEta[0]=etaMin;

  for (int iEta=1;iEta<etaNBins+1;iEta++)
    scaleEta[iEta] = etaMin + (1.6*iEta/etaNBins); 
}  


#endif
