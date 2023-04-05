{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Draw histos filled by Geant4 simulation
  //

  // Open file filled by Geant4 simulation
  TFile* file = new TFile("B3Test.root");

  // Create a canvas and divide it into 2x2 pads
 // TCanvas* c1 = new TCanvas("c1", "",0., 0., 1000, 1000);

  // Draw Eabs histogram in the pad 1
 //// c1->cd(1);
  TH1F* hist = (TH1F*)file->Get("Length");
 // hist1->Draw("0lego2 PFC");
 hist->SetFillColor(kOrange);
 hist->SetFillStyle(3001);
 hist->Draw("HIST BOX");
 }
