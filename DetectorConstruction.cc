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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Element.hh"

namespace E1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
 // G4double env_sizeXY = 3*cm, env_sizeZ = 30*cm;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_Cu");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 25.*cm;
  G4double world_sizeZ  = 7.*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* tube_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* catod_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* catodprot_mat = nist->FindOrBuildMaterial("G4_Au");
  G4Material* anod_mat = nist->FindOrBuildMaterial("G4_W");
 
 
 
  G4Tubs* solidWorld =
  new G4Tubs("World", 0, 0.5*world_sizeXY, 0.5*world_sizeZ, 0., CLHEP::pi*2);                 //its name
     
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //Straw treker
  // PET Tube
  //
  //G4double tube_dPhi = CLHEP::pi*2/6;
  G4double within = 1*cm;
  G4double r_L = 0*cm, R_L = 5*cm;
 
 //Layer
 
  G4Tubs* solidL =
    new G4Tubs("L",                    //its name
       0*cm, 5*cm, 0.5*world_sizeZ, 0., CLHEP::pi*2); //its size

  G4LogicalVolume* logicL =
    new G4LogicalVolume(solidL,            //its solid
                        tube_mat,             //its material
                        "L");         //its name
  G4VPhysicalVolume* physL =
   new G4PVReplica("physL",                       //no rotation
                 logicL,         //at (0,0,0)
                    logicWorld,                //its logical volume
                    kRho,              //its name
                    4,              //its mother  volume
                   5*cm,
                   0);          //overlaps checking
                   
                   
  G4Tubs* solidPET =
    new G4Tubs("PET",                    //its name
        2.5*cm, 5*cm, 0.5*world_sizeZ, 0., CLHEP::pi*2); //its size

  G4LogicalVolume* logicPET =
    new G4LogicalVolume(solidPET,            //its solid
                        tube_mat,             //its material
                        "PET");         //its name
  G4VPhysicalVolume* physPET =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicPET,            //its logical volume
                      "PET",               //its name
                      logicL,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                     checkOverlaps);        //overlaps checking
                   
G4Tubs* solidAir =
    new G4Tubs("Air",                    //its name
        0, 2.5*cm, 0.5*world_sizeZ, 0., CLHEP::pi*2); //its size

  G4LogicalVolume* logicAir =
    new G4LogicalVolume(solidAir,            //its solid
                        tube_mat,             //its material
                        "Air");         //its name
                        
  G4VPhysicalVolume* physAir =
   new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicAir,            //its logical volume
                      "Air",               //its name
                      logicL,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                    
 
  

  //
  fScoringVolume = logicL;
  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  G4VisAttributes* PETTube = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  PETTube->SetVisibility(true);
 logicL->SetVisAttributes(PETTube);
 // logicAir->SetVisAttributes (G4VisAttributes::GetInvisible());
  
 
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
