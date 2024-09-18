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
{
	delete fStepLimit;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  static const std::string fileName = "/home/sveta/GeoModel/GeoModelExamples/geostraw.db";
	//static const std::string fileName = "/home/sveta/GeoModel/GeoModelExamples/cubik.db";
	
	// check if DB file exists. If not, return.
	std::ifstream infile(fileName.c_str());
	if (!infile.good()) 
	{
        	std::cout << "\n\tERROR!! A '" << fileName
                  << "' file does not exist!! Please, check the path of the "
                     "input file before running this program. Exiting...";
        	exit(EXIT_FAILURE);
        }
        
        infile.close();
        
        // open the DB
        GMDBManager* db = new GMDBManager(fileName);
        
        /* Open database */
        if (db->checkIsDBOpen()) 
        {
		std::cout << "OK! Database is open!\n";
	} 
	else
	{
		std::cout << "Database is not open!\n";
		exit(EXIT_FAILURE);
	}
	
	// -- testing the input database
	std::cout << "\n=== Printing the list of all tables in the input DB:" << std::endl;
	db->printAllDBTables();
	
	std::cout << "\n=== Printing the list of all GeoMaterial nodes:" << std::endl;
	db->printAllMaterials();
	
	std::cout << "\n=== Printing the list of all GeoLogVol nodes:" << std::endl;
	db->printAllLogVols();
	
	/* setup the GeoModel reader */
	GeoModelIO::ReadGeoModel readInGeo = GeoModelIO::ReadGeoModel(db);
	std::cout << "ReadGeoModel set.\n";
	
	/* build the GeoModel geometry */
	const GeoVPhysVol* world =
        readInGeo.buildGeoModel();  // builds the whole GeoModel tree in memory
                                    // and get an handle to the 'world' volume
	std::cout << "ReadGeoModel::buildGeoModel() done.\n";
	
	// --- testing the imported geometry

	    // get the GeoLogVol used for the 'world' volume
	    std::cout << "Getting the GeoLogVol used by the 'world' volume"
		      << std::endl;
	    const GeoLogVol* logVol = world->getLogVol();
	    std::cout << "'world' GeoLogVol name: " << logVol->getName() << std::endl;
	    std::cout << "'world' GeoMaterial name: "
		      << logVol->getMaterial()->getName() << std::endl;

	    // get number of children volumes
	    unsigned int nChil = world->getNChildVols();
	    std::cout << "'world' number of children: " << nChil << std::endl;

	    // loop over all children nodes
	    std::cout << "Looping over all 'volume' children (i.e., GeoPhysVol and "
		         "GeoFullPhysVol)..."
		      << std::endl;
	    for (unsigned int idx = 0; idx < nChil; ++idx) 
	    {
		PVConstLink nodeLink = world->getChildVol(idx);

		if (dynamic_cast<const GeoVPhysVol*>(&(*(nodeLink)))) {
		    std::cout << "\t"
		              << "the child n. " << idx << " ";
		    const GeoVPhysVol* childVolV = &(*(nodeLink));
		    if (dynamic_cast<const GeoPhysVol*>(childVolV))

		    {
		        const GeoPhysVol* childVol =
		            dynamic_cast<const GeoPhysVol*>(childVolV);
		        std::cout << "is a GeoPhysVol, whose GeoLogVol name is: "
		                  << childVol->getLogVol()->getName() << std::endl;
		        std::cout << " and it has  " << childVol->getNChildVols()
		                  << " child volumes\n";

		    } else if (dynamic_cast<const GeoFullPhysVol*>(childVolV)) {
		        const GeoFullPhysVol* childVol =
		            dynamic_cast<const GeoFullPhysVol*>(childVolV);
		        std::cout << "is a GeoFullPhysVol, whose GeoLogVol name is: "
		                  << childVol->getLogVol()->getName() << std::endl;
		        std::cout << " and it has  " << childVol->getNChildVols()
		                  << " child volumes\n";
		    }
		} else if (dynamic_cast<const GeoNameTag*>(&(*(nodeLink)))) {
		    std::cout << "\t"
		              << "the child n. " << idx << " is a GeoNameTag\n";
		    const GeoNameTag* childVol =
		        dynamic_cast<const GeoNameTag*>(&(*(nodeLink)));
		    std::cout << "\t\tGeoNameTag's name: " << childVol->getName()
		              << std::endl;
		    // std::cout<< " and it has  "<<childVol->getNChildVols()<<" child
		    // volumes\n";
		}
	    }

	// build the Geant4 geometry and get an hanlde to the world' volume
	ExtParameterisedVolumeBuilder* builder = new ExtParameterisedVolumeBuilder("cubik");
	std::cout << "Building G4 geometry." << std::endl;
	
	//G4LogicalVolume* g4World = builder->Build(world);
	g4World = builder->Build(world);
	
	std::cout << "This is the newly-created Geant4 G4LogicalVolume, ready to be used: " << g4World << std::endl;
	
	G4VPhysicalVolume* physWorld =
		new G4PVPlacement(0,                     //no rotation
		G4ThreeVector(),       //at (0,0,0)
		g4World,            //its logical volume
		"World",               //its name
		0,                     //its mother  volume
		false,                 //no boolean operation
		0,                     //copy number
		true);        //overlaps checking

	std::cout << "This is the newly-created Geant4 G4VPhysicalVolume, ready to be used: " << physWorld << std::endl;

	G4cout << "*********** List of registered logical volumes *************" << G4endl;
	
	std::vector <G4String> SensLVnames;
	
	std::vector<G4LogicalVolume*>* lvStore = G4LogicalVolumeStore::GetInstance();
        std::size_t nlv=lvStore->size();
        
        std::cout<<"G4LogicalVolume store size: "<<nlv<<std::endl;
        
        for (std::size_t i=0;  i<nlv; i++) 
        {
        
        	G4LogicalVolume* lv = (*lvStore)[i];
        	
        	G4String name = lv->GetName();
        	G4cout << name << G4endl;
        	
        	if (name.contains("_sens"))
        	{
        	
        		SensLVnames.push_back(name);
        		std::cout<< "Volume with "<<name << " is sensetive" << std::endl;
        	}
        
        }
        
	std::cout << "*************************CREATION SENSITIVE VOLUMES ***************" << std::endl;
	/*for(int i=0; i<SensLVnames.size(); i++)
	{
		G4String names = SensLVnames[i];
		
		std::cout << i <<" LV names "  << names << std::endl;
		cubik_box = G4LogicalVolumeStore::GetInstance()->GetVolume(names);
	
	}*/
	G4double maxStep = 0.1*mm;
  	fStepLimit = new G4UserLimits(maxStep);
 	
        for(int i=0; i<SensLVnames.size(); i++)
	{
		G4String names = SensLVnames[i];
		cubik_box = G4LogicalVolumeStore::GetInstance()->GetVolume(names);
		cubik_box->SetUserLimits(fStepLimit);
		std::cout << i <<" LV names "  << names << std::endl;
		std::cout <<" cubik "  << cubik_box << std::endl;
	}		
	
	g4World->SetVisAttributes (G4VisAttributes::GetInvisible());
	
	
	
	
return physWorld;
	
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	G4cout << " Construct !!!!!!!" << G4endl;
	//G4String names = SensLVnames[0];
	G4String BoxName = cubik_box->GetName();
	/*SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");
	cubik_box->SetSensitiveDetector(sensDet);
	//cubik_box2->SetSensitiveDetector(sensDet);*/

	G4ThreeVector fieldValue = G4ThreeVector(0.,0.,1.*tesla);
  	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  	fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  	G4AutoDelete::Register(fMagFieldMessenger);
  	
	EventAction* fEvent_Action = nullptr;
  
  // Sensitive detectors
  G4String trackerChamberSDname = "/TrackerChamberSD";
  TrackerSD* aTrackerSD = new TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection", fEvent_Action);
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  
  SetSensitiveDetector(BoxName, aTrackerSD, true);
  G4cout << " NAME " << BoxName << G4endl;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}



void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
