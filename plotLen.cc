{
#include "TLegend.h"
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Draw histos filled by Geant4 simulation
  //

  // Open file filled by Geant4 simulation
  TFile* file = new TFile("Test.root");

  // Create a canvas and divide it into 2x2 pads
 // TCanvas* c1 = new TCanvas("c1", "",0., 0., 1000, 1000);

  // Draw Eabs histogram in the pad 1
 //// c1->cd(1);
 
 /*TH1F* hist1 = (TH1F*)file->Get("Number-Norm");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Distance , mm");
 hist1->GetYaxis()->SetTitle("Number of muons");*/
 
 TH1F* hist1 = (TH1F*)file->Get("Length");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Distance , mm");
 hist1->GetYaxis()->SetTitle("Number of muons");
 
 /*TH1F* hist1 = (TH1F*)file->Get("Start Time");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Time of Hits, ns");
 hist1->GetYaxis()->SetTitle("Number of muons");

 
 
  TH1F* hist2 = (TH1F*)file->Get("Time3");
 hist2->SetFillColor(kYellow);
 hist2->SetFillStyle(3001);
 hist2->Draw("HIST BOX SAME");
 
 
  TH1F* hist3 = (TH1F*)file->Get("Time4");
 hist3->SetFillColor(kBlue);
 hist3->SetFillStyle(3001);
 hist3->Draw("HIST BOX SAME");
 
 
  TH1F* hist4 = (TH1F*)file->Get("Time5");
 hist4->SetFillColor(kRed);
 hist4->SetFillStyle(3001);
 hist4->Draw("HIST BOX SAME");
 
 
  TH1F* hist5 = (TH1F*)file->Get("Time6");
 hist5->SetFillColor(kGreen);
 hist5->SetFillStyle(3001);
 hist5->Draw("HIST BOX SAME");
 
 
  TH1F* hist6 = (TH1F*)file->Get("Time7");
 hist6->SetFillColor(kOrange);
 hist6->SetFillStyle(3001);
 hist6->Draw("HIST BOX SAME");
 
 
  TH1F* hist7 = (TH1F*)file->Get("Time8");
 hist7->SetFillColor(kPink);
 hist7->SetFillStyle(3001);
 hist7->Draw("HIST BOX SAME");
 
 
  TH1F* hist8 = (TH1F*)file->Get("Time9");
 hist8->SetFillColor(kTeal);
 hist8->SetFillStyle(3001);
 hist8->Draw("HIST BOX SAME");

 TH1F* hist9 = (TH1F*)file->Get("Time10");
 hist9->SetFillColor(kMagenta);
 hist9->SetFillStyle(3001);
 hist9->Draw("HIST BOX SAME");
 
 TH1F* hist10 = (TH1F*)file->Get("Time11");
 hist10->SetFillColor(kYellow-20);
 hist10->SetFillStyle(3001);
 hist10->Draw("HIST BOX SAME");

 TH1F* hist11 = (TH1F*)file->Get("Time12");
 hist11->SetFillColor(kSpring);
 hist11->SetFillStyle(3001);
 hist11->Draw("HIST BOX SAME");
 
 TH1F* hist12 = (TH1F*)file->Get("Time13");
 hist12->SetFillColor(kAzure);
 hist12->SetFillStyle(3001);
 hist12->Draw("HIST BOX SAME");
 
 TH1F* hist13 = (TH1F*)file->Get("Time14");
 hist13->SetFillColor(kCyan);
 hist13->SetFillStyle(3001);
 hist13->Draw("HIST BOX SAME");
 
 TH1F* hist14 = (TH1F*)file->Get("Time15");
 hist14->SetFillColor(kRed);
 hist14->SetFillStyle(3001);
 hist14->Draw("HIST BOX SAME");*/
  
 

 }
