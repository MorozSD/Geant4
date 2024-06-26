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
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"
//#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
//#include "GeoTube.h"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "EventAction.hh"
#include "GeoModelKernel/GeoBox.h"

#include "G4UserLimits.hh"




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}


/*DetectorConstruction::DetectorConstruction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();
 // fUserLimits = new G4UserLimits();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{

  // Material definition

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* tube_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* anod_mat = nist->FindOrBuildMaterial("G4_W");
  
  // CO2
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
  Gas_mat->AddMaterial(Ar, fractionmass=70*perCent);
  Gas_mat->AddMaterial(CO2, fractionmass=30*perCent);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 // Detector Parameters
  
  G4int num_rings = 30, nbn_cryst = 500;
  G4double cryst_r = 5*mm, cryst_x = 236*cm;
  G4double ring_R1 = 250*mm, ring_r, detector_dZ = cryst_x;
  G4double ring_R2 = (ring_R1+2.0*num_rings*cryst_r);
 
  
  // Definitions of Solids, Logical Volumes, Physical Volumes

  //
  // World
  //
  G4double world_sizeXY = 2.4*ring_R2;
  G4double world_sizeZ  = 1.2*detector_dZ;

  G4Box* solidWorld =
    new G4Box("World",                                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size

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
                      fCheckOverlaps);       // checking overlaps
                      
                      
  //
  // ring
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2, cryst_x/2.0, 0., twopi);

  G4LogicalVolume* logicRing =
    new G4LogicalVolume(solidRing,           //its solid
                        world_mat,         //its material
                        "Ring");             //its name

  //
  // define crystal
  //
              
  G4double dr = cryst_r;
  G4Tubs* solidCryst = new G4Tubs("crystal", 0.0, dr, cryst_x/2.0, 0., twopi);

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,        //its solid
                        world_mat,         //its material
                        "Crystal");        //its name
                        

  
// place shapes into crystals
G4double PETgap = 0.036*mm;
//G4double PETgap = 1.0*mm;
G4double r_PET = dr - PETgap, R_PET = dr;
G4Tubs* solidPET =
    new G4Tubs("PET",                    //its name
        r_PET, R_PET, 0.5*cryst_x, 0., CLHEP::pi*2); //its size

  G4LogicalVolume* logicPET =
    new G4LogicalVolume(solidPET,            //its solid
                        tube_mat,             //its material
                        "PET");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicPET,                //its logical volume
                    "PET",              //its name
                    logicCryst,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking
                
                    
  G4int copyN = 0;     
 //Gas                
 G4double r_Gas = 0.03*mm, R_Gas = r_PET ;
//G4double r_Gas = 0.05*mm, R_Gas = dr ;
  G4Tubs* solidGas =
    new G4Tubs("Gas",                    //its name
        r_Gas, R_Gas, 0.5*cryst_x, 0., CLHEP::pi*2); //its size

  G4LogicalVolume* logicGas =
    new G4LogicalVolume(solidGas,            //its solid
                       Gas_mat,             //its material
                        "Gas");         //its name

                    new G4PVPlacement(0,     //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicGas,                //its logical volume
                    "Gas",              //its name
                    logicCryst,              //its mother  volume
                    false,                   //no boolean operation
                    copyN,                       //copy number
                    fCheckOverlaps);          //overlaps checking 
                    
  //Wolfram wire                 
 G4double r_W = 0*cm, R_W = 0.03*mm;
  G4Tubs* solidW =
    new G4Tubs("W",                    //its name
        r_W, R_W, 0.5*cryst_x, 0., CLHEP::pi*2); //its size

  G4LogicalVolume* logicW =
    new G4LogicalVolume(solidW,            //its solid
                        anod_mat,             //its material
                        "W");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicW,                //its logical volume
                    "W",              //its name
                    logicCryst,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking                    
                    
                    
  
                                
        // place crystals within a ring
  //
  
  // G4int nbn_cryst = 120; used above
  
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
    copyN = i_ring*1000+icrys;

    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "Crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      copyN,                 //copy number
                      fCheckOverlaps);       // checking overlaps
                                    
  }
  }
                    
 //detector 
 
 G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);

  G4LogicalVolume* logicDetector =
    new G4LogicalVolume(solidDetector,       //its solid
                        world_mat,         //its material
                        "Detector");         //its name
                        
     new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    copyN,                       //copy number
                    fCheckOverlaps);         // checking overlaps                   

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
                      
   
                        
                         
   

 

 logicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  //logicGas->SetVisAttributes (G4VisAttributes::GetInvisible());

  G4VisAttributes* PETTube = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
 
  logicPET->SetVisAttributes(PETTube);
  
  G4VisAttributes* Wolfram = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  logicW->SetVisAttributes(Wolfram);
   G4VisAttributes* GasTube = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicGas->SetVisAttributes(GasTube);
 // logicChamber->SetVisAttributes(Tube);
  // Print materials
 //
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  G4double maxStep = 0.1*mm;
  fUserLimit = new G4UserLimits(maxStep);
  logicGas->SetUserLimits(fUserLimit);
  
  
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{

  EventAction* fEvent_Action = nullptr;
  
  // Sensitive detectors
   G4String trackerChamberSDname = "/TrackerChamberSD";
  TrackerSD* aTrackerSD = new TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection", fEvent_Action);
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  
  SetSensitiveDetector("Gas", aTrackerSD, true);
  
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector(0.,0.,1.*tesla);
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

/*void DetectorConstruction::SetMaxStepSize(G4double maxStep)
{
        maxStep = 0.001*mm;
	//myStepLimit = new G4UserLimits(maxStep);
	fUserLimits = new G4UserLimits();
	fUserLimits->SetMaxAllowedStep(maxStep);

	//logicWorld->SetUserLimits(myStepLimit);
}*/



void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
