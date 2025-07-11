//
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
/// \file eventgenerator/userPrimaryGenerator/userPrimaryGenerator.cc
/// \brief Main program of the eventgenerator/userPrimaryGenerator example
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"
#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
//#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
//#include "FTFP_BERT.hh"
#include "QBBC.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIcommand.hh"
#include "G4StepLimiterPhysics.hh"

#include "Randomize.hh"
#include <fstream>
#include <chrono>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {



  unsigned int seed;
    std::ifstream urandom("/dev/urandom", std::ios::binary);
    
    if (urandom) {
        urandom.read(reinterpret_cast<char*>(&seed), sizeof(seed));
        urandom.close();
    } else {
        // Fallback: используем время, если /dev/urandom недоступен
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    
    CLHEP::HepRandom::setTheSeed(seed);
    
    
    

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //choose the Random engine
 // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheEngine(new CLHEP::Ranlux64Engine);
  
  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  //construct the default run manager
  G4RunManager* runManager = new G4RunManager;
//auto* runManager =
//   G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  //set mandatory initialization classes
  //
  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);
  

  auto physicsList = new QBBC;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);


/*      G4VModularPhysicsList* physicsList = new QBBC;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);*/
//  auto actionInitialization = new ActionInitialization(detConstruction);
  auto actionInitialization = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);
  
 // runManager->SetUserInitialization(new DetectorConstruction);

  //runManager->SetUserInitialization(new PhysicsList);
  
  

  runManager->SetUserAction(new PrimaryGeneratorAction);

  //initialize visualization
  G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive;
   visManager->Initialize();
   UImanager->ApplyCommand("/control/execute vis.mac");          
   ui->SessionStart();
   delete ui;
  }
  else  {
   //batch mode  
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
  }

  //job termination 
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

