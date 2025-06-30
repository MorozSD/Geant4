{
#include "TLegend.h"
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Open file filled by Geant4 simulation
  TFile* file = new TFile("Test.root");
 /*auto rdf = ROOT::RDF::FromCSV("Z_T.csv");
 auto h = rdf.Histo2D("Z0","Time");
 h->Draw();*/
  
/* TH2F* hist1 = (TH2F*)rdf->Get("T-Z0");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST COL");
 hist1->GetXaxis()->SetTitle(" Z0, cm ");
 hist1->GetYaxis()->SetTitle("time, ns");*/
 
 TH2F* hist1 = (TH2F*)file->Get("Z0-theta dist pZ");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST COL");
 hist1->GetXaxis()->SetTitle(" theta, rad ");
 hist1->GetYaxis()->SetTitle("Z0-Z_nach, mm");
 
 
 
/* TH2F* hist1 = (TH2F*)file->Get("Zraz-R2Adj");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST COL");
 hist1->GetXaxis()->SetTitle(" R2Adj ");
 hist1->GetYaxis()->SetTitle("Z0-Z_nach, mm");*/
 
/* TH1F* hist1 = (TH1F*)file->Get("Znach");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Znach , mm");
 hist1->GetYaxis()->SetTitle("Number of tracks");
 
 TH1F* hist2 = (TH1F*)file->Get("Z_result");
 hist2->SetFillColor(kYellow);
 hist2->SetFillStyle(3001);
 hist2->Draw("HIST BOX SAME");*/
 
 
 /*TH1F* hist2 = (TH1F*)file->Get("Z1");
 hist2->SetFillColor(kYellow);
 hist2->SetFillStyle(3001);
 hist2->Draw("HIST BOX SAME");*/
 
 /* TH1F* hist3 = (TH1F*)file->Get("Z2");
 hist3->SetFillColor(kBlue);
 hist3->SetFillStyle(3001);
 hist3->Draw("HIST BOX SAME");
 
  TH1F* hist4 = (TH1F*)file->Get("Z3");
 hist4->SetFillColor(kRed);
 hist4->SetFillStyle(3001);
 hist4->Draw("HIST BOX SAME");
 
  TH1F* hist5 = (TH1F*)file->Get("Z4");
 hist5->SetFillColor(kGreen);
 hist5->SetFillStyle(3001);
 hist5->Draw("HIST BOX SAME");
 
  TH1F* hist6 = (TH1F*)file->Get("Z5");
 hist6->SetFillColor(kYellow);
 hist6->SetFillStyle(3001);
 hist6->Draw("HIST BOX SAME");
 
  TH1F* hist7 = (TH1F*)file->Get("Z6");
 hist7->SetFillColor(kBlue);
 hist7->SetFillStyle(3001);
 hist7->Draw("HIST BOX SAME");
 
  TH1F* hist8 = (TH1F*)file->Get("Z7");
 hist8->SetFillColor(kRed);
 hist8->SetFillStyle(3001);
 hist8->Draw("HIST BOX SAME");
 
  TH1F* hist9 = (TH1F*)file->Get("Z8");
 hist9->SetFillColor(kGreen);
 hist9->SetFillStyle(3001);
 hist9->Draw("HIST BOX SAME");
 
  TH1F* hist10 = (TH1F*)file->Get("Z9");
 hist10->SetFillColor(kYellow);
 hist10->SetFillStyle(3001);
 hist10->Draw("HIST BOX SAME");
 
  TH1F* hist11 = (TH1F*)file->Get("Z10");
 hist11->SetFillColor(kPink);
 hist11->SetFillStyle(3001);
 hist11->Draw("HIST BOX SAME");
 
  TH1F* hist12 = (TH1F*)file->Get("Z11");
 hist12->SetFillColor(kOrange);
 hist12->SetFillStyle(3001);
 hist12->Draw("HIST BOX SAME");
 
  TH1F* hist13 = (TH1F*)file->Get("Z12");
 hist13->SetFillColor(kBlue);
 hist13->SetFillStyle(3001);
 hist13->Draw("HIST BOX SAME");
 
  TH1F* hist14 = (TH1F*)file->Get("Z13");
 hist14->SetFillColor(kOrange);
 hist14->SetFillStyle(3001);
 hist14->Draw("HIST BOX SAME");
 
  TH1F* hist15 = (TH1F*)file->Get("Z14");
 hist15->SetFillColor(kGreen);
 hist15->SetFillStyle(3001);
 hist15->Draw("HIST BOX SAME");*/
 
 /*TH1F* hist1 = (TH1F*)file->Get("Length");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Distance , mm");
 hist1->GetYaxis()->SetTitle("Number of muons");*/
 
 /*TH1F* hist2 = (TH1F*)file->Get("Length_0");
 hist2->SetFillColor(kGreen);
 hist2->SetFillStyle(3001);
 hist2->Draw("HIST BOX SAME");
 
 TH1F* hist3 = (TH1F*)file->Get("Length_10");
 hist3->SetFillColor(kBlue);
 hist3->SetFillStyle(3001);
 hist3->Draw("HIST BOX SAME");
 
 TH1F* hist4 = (TH1F*)file->Get("Length_100");
 hist4->SetFillColor(kRed);
 hist4->SetFillStyle(3001);
 hist4->Draw("HIST BOX SAME");
 
 auto legend = new TLegend(0.1,0.7,0.48,0.9);
 legend->SetHeader("All Limits","C"); // option "C" allows to center the header
 legend->AddEntry(hist1,"No Limits","f");
 legend->AddEntry(hist2,"Edep > 5 keV","f");
 legend->AddEntry(hist3,"Edep > 10 keV","f");
 legend->AddEntry(hist4,"Edep > 50 keV","f");
 legend->Draw();*/
 
 /*TH1F* hist1 = (TH1F*)file->Get("Edep");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Edep , eV");
 hist1->GetYaxis()->SetTitle("Number of muons");*/
 
 
 /*TH1F* hist1 = (TH1F*)file->Get("Start Time");
 hist1->SetFillColor(kBlack);
 hist1->SetFillStyle(3001);
 hist1->Draw("HIST BOX");
 hist1->GetXaxis()->SetTitle("Time, ns");
 //hist1->GetXaxis()->SetTitle("|Z_{T} - Z_{c}|, mm");
 hist1->GetYaxis()->SetTitle("Number of intersections");

 
 
 TH1F* hist2 = (TH1F*)file->Get("Time3");
 hist2->SetFillColor(kYellow);
 hist2->SetFillStyle(3001);
 hist2->Draw("HIST BOX SAME");*/
 
 /*
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
 hist14->Draw("HIST BOX SAME");
 
 TH1F* hist15 = (TH1F*)file->Get("Time16");
 hist15->SetFillColor(kBlue);
 hist15->SetFillStyle(3001);
 hist15->Draw("HIST BOX SAME");*/
  
 

 }
