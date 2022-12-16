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
/// \brief Implementation of the B3::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"


namespace B3
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  //DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void DetectorConstruction::DefineMaterials()
//{
 // G4NistManager* man = G4NistManager::Instance();

  //G4bool isotopes = false;

  //G4Element*  O = man->FindOrBuildElement("O" , isotopes);
 // G4Element* Si = man->FindOrBuildElement("Si", isotopes);
 // G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

 // G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
//  LSO->AddElement(Lu, 2);
//  LSO->AddElement(Si, 1);
//  LSO->AddElement(O , 5);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Gamma detector Parameters
  //
  G4double cryst_r = 10*mm, cryst_x = 54*cm;
  G4int nb_cryst = 100;
  G4int num_rings = 8;
  G4int nbn_cryst = 150;
  //
  //G4double dPhi = twopi/nbn_cryst, half_dPhi = 0.5*dPhi;
//  G4double cosdPhi = std::cos(half_dPhi);
 // G4double tandPhi = std::tan(half_dPhi);
  //
  G4double ring_R1 = 275*mm;
  G4double ring_R2 = (ring_R1+2.0*num_rings*cryst_r);
  G4double ring_r;
  //
  G4double detector_dZ = cryst_x;
  //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

  //
  // World
  //
  G4double world_sizeXY = 2.4*ring_R2;
  G4double world_sizeZ  = 1.2*detector_dZ;

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps

  //
  // ring
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2, cryst_x/2.0, 0., twopi);

  G4LogicalVolume* logicRing =
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name

  //
  // define crystal
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dr = cryst_r - gap;
  G4Tubs* solidCryst = new G4Tubs("crystal", 0.0, dr, cryst_x/2.0, 0., twopi);

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name

  // place crystals within a ring
  //
   for (G4int i_ring = 1; i_ring < num_rings+1; i_ring++) {
  ring_r = ring_R2-2.0*cryst_r*i_ring;
  G4double phi0 = twopi/nbn_cryst;
  G4double R0 = (ring_r + cryst_r)*std::sqrt(2*(1-std::cos(phi0)));
  
  while (R0 < 2.0*cryst_r)
  {
  nbn_cryst = nbn_cryst - 1;
  phi0 = twopi/nbn_cryst;
  R0 = (ring_r + cryst_r)*std::sqrt(2*(1-std::cos(phi0)));
  }
  G4int nprav_cryst = nbn_cryst;
  for (G4int icrys = 0; icrys < nprav_cryst ; icrys++) {
    G4double phi = icrys*twopi/nprav_cryst;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
    G4ThreeVector position = (ring_r+cryst_r)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);

    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps
  }
  }

  //
  // full detector
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);

  G4LogicalVolume* logicDetector =
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name

  //
  // place rings within detector
  //
 new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,0), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      logicDetector,         //its mother  volume
                      false,                 //no boolean operation
                      1,                 //copy number
                      fCheckOverlaps);       // checking overlaps


  //
  // place detector in world
  //
 new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps

  

 // Visualization attributes
  //
  
  logicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  //logicPatient->SetVisAttributes (G4VisAttributes::GetInvisible());
  
 
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;


  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

