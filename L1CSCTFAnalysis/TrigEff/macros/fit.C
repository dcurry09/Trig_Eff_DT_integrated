#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"



double x[17]  = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 20, 30, 37.25, 44.75, 62.5, 87.5};

double y5[17] = {0, 0, 0.302736, 0.47936, 0.727172, 0.878876, 0.948303, 0.972536, 0.980288, 0.983845, 0.987268, 0.987909, 0.986143, 0.982857,  0.970833, 0.97561, 0.954023}; 

double y7[17] = {0, 0, 0.275364, 0.439841, 0.673502,  0.826439, 0.896762, 0.934192, 0.95531, 0.969709, 0.976377, 0.983181, 0.981415, 0.980571, 0.970833, 0.971545, 0.954023};

TF1 *fit5 = new TF1("fit5","[0]* (1 + TMath::Erf( (x - [1])*[2] ))/2.0",0,100);
TF1 *fit7 = new TF1("fit7","[0]* (1 + TMath::Erf( (x - [1])*[2] ))/2.0",0,100);

//TH1F *hist = new TH1F("hist","",1,1,100);

void fit() {
  
  fit5 -> SetParameter(0,4.0);
  fit5 -> SetParameter(1,5.0);
  fit5 -> SetParameter(2,2.0);
  
  fit5 -> SetParLimits(0,1,10);
  fit5 -> SetParLimits(1,1,10);
  fit5 -> SetParLimits(2,0,5);
  
  fit7 -> SetParameter(0,1.0);
  fit7 -> SetParameter(1,5.0);
  fit7 -> SetParameter(2,0.5);
  
  fit7 -> SetParLimits(0,1,10);
  fit7 -> SetParLimits(1,1,10);
  fit7 -> SetParLimits(2,0,10);

  TCanvas *c5 = new TCanvas("c5","",600,600);  
  TGraph* gr5 = new TGraph(17, x, y5); 
  c5 -> SetLogx();
  gr5 -> Fit("fit5", "r");
  gr5 -> GetFunction("fit5")->SetLineColor(kBlue);
  gr5 -> SetMarkerColor(kBlue);
  gr5 -> SetLineColor(kBlue);
  gr5 -> SetMarkerStyle(21);
  gr5 -> SetMarkerSize(0.6);
  gr5 -> SetLineWidth(2);
  gr5 -> Draw("AP"); 
 
  TCanvas *c7 = new TCanvas("c7","",600,600);
  TGraph* gr7 = new TGraph(17, x, y7);
  c7 -> SetLogx();
  gr7 -> Fit("fit7", "r");
  gr7 -> GetFunction("fit7")->SetLineColor(kRed);
  gr7 -> SetMarkerColor(kRed);
  gr7 -> SetLineColor(kRed);
  gr7 -> SetMarkerStyle(21);
  gr7 -> SetMarkerSize(0.6);
  gr7 -> SetLineWidth(2);
  gr7 -> Draw("AP");



}// end fit
