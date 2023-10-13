// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RootAnalysisReader.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each event
 // G4RunManager::GetRunManager()->SetPrintProgress(1000);
  G4RunManager::GetRunManager()->SetPrintProgress(1);
  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  // analysisManager->CreateH1("0","Time", 1, 0., 80.*ns);
  // analysisManager->CreateH1("Start Time","Timing distribution", 100, -10.*ns, 1000.*ns);
  // analysisManager->SetH1XAxisTitle(0, "Time of Hits, ns");
   ///analysisManager->SetH1YAxisTitle(0, "Number of muons");
 //  analysisManager->CreateH1("Start Position","Position distribution", 100, -25.*cm, 25.*cm);
  // analysisManager->CreateH1("StartPos","Shortest distance", 100, -20.*cm, 20.*cm);
  // analysisManager->CreateH1("Length","Shortest distance", 100, -0.05*mm, 5.*mm);
   //analysisManager->CreateH1("Start Time","Start Timing distribution", 100, -1.*ns, 1200.*ns);
  // analysisManager->CreateH1("Time3","Timing distribution no drift", 100, -1.*ns, 1200.*ns);
   //analysisManager->CreateH1("Time4","Timing distribution with drift", 100, -1.*ns, 1100.*ns);
  // analysisManager->CreateH1("Z(x), Z(y)","Z", 100, -1.5*mm, 1.5*mm);
  // analysisManager->CreateH1("Zy","Z", 100, -1.5*mm, 1.5*mm);
  
   //analysisManager->CreateH1("StartPos","Z", 100, -21.0*cm, 21.0*cm);
  // analysisManager->CreateH1("Zx","X(z)", 100, -21.0*cm, 21.0*cm);
   //analysisManager->CreateH1("Zy","Y(z)", 100, -21.0*cm, 21.0*cm);
 //  analysisManager->CreateH1("R1","RAdj", 100, -0.01, 0.1);
 //  analysisManager->CreateH1("P-Value","P-Value", 100, -0.01, 0.1);
 //  analysisManager->CreateH1("R3","RAdj", 100, -0.01, 0.1);
   analysisManager->CreateH1("X(Z), Y(Z)","Z-Z0", 100, -10.0*cm, 10.0*cm);
   analysisManager->CreateH1("Zy","Z-Z0", 100, -10.0*cm, 10.0*cm);
   
  // analysisManager->CreateNtuple("Test", "STPOS");
   //analysisManager->CreateNtupleDColumn("Start Time");
   //analysisManager->CreateNtupleDColumn("Length");
   //analysisManager->CreateNtupleDColumn("Time3");
   //analysisManager->CreateNtupleDColumn("Time4");
//   analysisManager->CreateNtupleDColumn("StartPos");
 //  analysisManager->CreateNtupleDColumn("Zx");
  // analysisManager->CreateNtupleDColumn("Zy");
 //  analysisManager->FinishNtuple();
  
 /*  analysisManager->CreateH1("Time5","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time6","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time7","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time8","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time9","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time10","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time11","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time12","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time13","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time14","Timing distribution", 100, -1.*ns, 1100.*ns);
   analysisManager->CreateH1("Time15","Timing distribution", 100, -1.*ns, 1100.*ns);*/
   //analysisManager->SetH1XAxisTitle(0, "Time of Hits, ns");
   //analysisManager->SetH1YAxisTitle(0, "Number of muons");
   
  // analysisManager->CreateH1("Length","Distance to tube center", 100, 0.1 *mm, 9.*mm);
 //  analysisManager->SetH1XAxisTitle(0, "Distance");
 //  analysisManager->SetH1YAxisTitle(0, "Number of muons");
   
   //analysisManager->CreateH1("1","Number", 10, 0., 20.*MeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  // G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  

  // Open an output file
  //
  G4String fileName = "Test.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
//  // G4String fileName = "B4.xml";
  analysisManager->OpenFile(fileName);
  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance(); 
  analysisManager->Write();
  analysisManager->CloseFile();
  
  analysisManager->SetNtupleMerging(true);
 
// G4String fileName0 = "Z.root";
// analysisManager->OpenFile(fileName0); 
  
 /*using G4AnalysisReader = G4RootAnalysisReader;

// Create (or get) analysis reader
   auto analysisReader = G4AnalysisReader::Instance();
 ///  auto analysisReader1 = G4AnalysisReader::Instance();
  // auto analysisReader2 = G4AnalysisReader::Instance();
   
   analysisReader->SetVerboseLevel(1);
          
        // Define a base file name
        analysisReader->SetFileName("Test.root");
     //   analysisReader1->SetFileName("Test.root");
//analysisReader2->SetFileName("Test.root");

        // Read ntuple
        G4int ntupleId = analysisReader->GetNtuple("Test");
        G4int counter = 0;
        if ( ntupleId >= 0 ) {
        G4double Z0, Z_X, Z_Y;
        analysisReader->SetNtupleDColumn("StartPos", Z0);
      //  analysisReader1->SetNtupleDColumn("Zx", Z_X);
     //   analysisReader2->SetNtupleDColumn("Zy", Z_Y);
        G4cout << "Ntuple Start Position, reading selected column StartPos" << G4endl;
        while ( analysisReader->GetNtupleRow()) {
          G4cout << counter++ << "th entry: " << " Start Position Z0: " << G4BestUnit(Z0,"Length") << G4endl;
        //  G4cout << counter << "th entry: " << " Start Position Zx: " << Z_X << G4endl;
        //  G4cout << counter << "th entry: " << " Start Position Zy: " << Z_Y << G4endl;
        }


  //analysisManager->FillH1(0, Zras);
}*/

 
  
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
