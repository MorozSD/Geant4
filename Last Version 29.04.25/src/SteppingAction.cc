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
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  
   G4double x = step->GetPreStepPoint()->GetPosition().x();
   G4double y = step->GetPreStepPoint()->GetPosition().y();
   G4double z = step->GetPreStepPoint()->GetPosition().z();
  G4int PDG = step->GetTrack()->GetDefinition()->GetPDGEncoding();
 // G4LogicalVolume* volume
  //  = step->GetPreStepPoint()->GetTouchableHandle()
      //->GetVolume()->GetLogicalVolume();
   G4int Track_ID = step->GetTrack()->GetTrackID();
  if( PDG == 13){    
 // G4cout << " step X = " << x << " step Y = " << y << " ID = " << Track_ID << G4endl;
  // collect energy deposited in this step

#ifdef G4MULTITHREADED
static G4Mutex stuffMutex = G4MUTEX_INITIALIZER;
G4AutoLock al(&stuffMutex);
#endif
static std::ofstream stuff_xy("test_xy.csv");
static std::ofstream stuff_z("test_z.csv");
static bool first = true;
if (first) {
first = false;
stuff_xy << "#,x/mm,y/mm,TrackID,PDG" << std::endl;
stuff_z << "#,z/mm,TrackID" << std::endl;
}
stuff_xy << "," << x << "," << y <<","<< Track_ID << "," << PDG << std::endl;
stuff_z << ","<< z << ","<< Track_ID << std::endl;
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

