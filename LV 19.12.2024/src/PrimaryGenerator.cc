#include "PrimaryGenerator.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "EventAction.hh"
///#include "RandPoisson.h"
#include "G4Poisson.hh"
#include "G4AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4Cache.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{


  G4double steptime = 76*ns;
  G4ThreeVector position;
         
  for (icollnum = 0; icollnum < 131; icollnum++){
       G4double ptime = steptime*icollnum;
  
       G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime);
       G4ParticleDefinition* particleDefinition
              = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
           
       G4int num_Ev = G4Poisson(0.3), num_Mu = G4Poisson(7);
       //G4int num_Ev = 1, num_Mu = 1;
           
    //   G4cout << " time " << ptime << " num_Ev " << num_Ev << " num_Mu " << num_Mu << G4endl;
           
       for (G4int i = 0; i < num_Ev; i++){
  
           G4double Z0 = G4RandGauss::shoot(0,30*cm);
           position = G4ThreeVector(0.0,0.0, Z0);
   
           for(G4int j = 0; j < num_Mu; j++){
               G4PrimaryParticle* p = new G4PrimaryParticle(particleDefinition);
 
               G4double tetta = G4UniformRand()*2*3.14159;
               G4double phi = G4UniformRand()*2*3.14159;
               G4double px = std::sin(tetta)*std::cos(phi);
               G4double py = std::sin(tetta)*std::sin(phi);
               G4double pz = std::cos(tetta);
               G4ThreeVector direction = G4ThreeVector(px,py,pz);
  
               p->SetMomentumDirection(direction);
               p->SetKineticEnergy(1*GeV);
               vertex->SetPrimary(p);
           }
        } 
        
        //  G4cout << " just tetta = " << tetta*180/3.14159  << G4endl;
 /* if(tetta>3.14159/2 && tetta<3*3.14159/2) G4cout << " tetta1 = " << std::fabs(tetta*180/3.14159-180) << G4endl;
  else if(tetta>3*3.14159/2 && tetta<2*3.14159) G4cout << " tetta2 = " << 2*180-tetta*180/3.14159 << G4endl;
  else G4cout << " tetta3 = " << tetta*180/3.14159 << G4endl;//G4cout << "Number = " << num << G4endl;*/
  //analysisManager->FillH1(7,position.z());
        event->AddPrimaryVertex(vertex);
   }
}  
  
  
 
