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
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

 G4NistManager* nist = G4NistManager::Instance();
  // get volume of the current step
//  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
   
 /*  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
   G4Material* tube_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
   G4Material* anod_mat = nist->FindOrBuildMaterial("G4_W");
   G4double z, a, density, fractionmass;
  G4String name, symbol;
  G4int ncomponents, natoms;
  
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);
 
  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
 
 density = 1.977*kg/m3;
 G4Material* CO2 = new G4Material(name="Carbon dioxide",density, ncomponents=2);
 CO2->AddElement(elC, natoms=1);
 CO2->AddElement(elO, natoms=2);
 
 
 density = 1.4*g/cm3; 
 G4Material* Ar = nist->FindOrBuildMaterial("G4_Ar");
 G4Material* Gas_mat = new G4Material(name="Gas_mat",density,ncomponents=2);
 Gas_mat->AddMaterial(Ar, fractionmass=30*perCent);
 Gas_mat->AddMaterial(CO2, fractionmass=70*perCent);

  G4Material* mat = step->GetPreStepPoint()->GetMaterial();
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  if(edepStep>0) //G4cout << "MATERIAL" << mat << G4endl;
  {
    if(mat==world_mat) {G4cout << "MATERIAL: AIR" << G4endl; fEventAction->AddMaterial(1,0,0);}
    if(mat==tube_mat) {G4cout << "MATERIAL: PET" << G4endl;fEventAction->AddMaterial(0,0,1);}
    if(mat==anod_mat) {fEventAction->AddMaterial(0,0,0);}
    if(mat==Gas_mat) {G4cout << "MATERIAL: GAS" << G4endl;fEventAction->AddMaterial(0,1,0);}
    
  
  }*/
}



