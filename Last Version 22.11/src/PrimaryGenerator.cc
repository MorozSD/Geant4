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

  G4int Num_of_Part = 0;
  G4double steptime = 76*ns;
  G4ThreeVector position;
  
  G4MapCache <G4double, G4int> Position_Map;
  
  auto analysisManager = G4AnalysisManager::Instance(); 
 // G4String fileName = "TestPos.root";
 // analysisManager->OpenFile(fileName);
      
  for (icollnum = 0; icollnum < 14; icollnum++){
 // for (icollnum = 0; icollnum < 2; icollnum++){
  G4double ptime = steptime*icollnum;
 // G4double Z0 =  G4UniformRand()*20*cm-10*cm;
 
  
 // G4cout << " POSITION = " << G4BestUnit(position.z(),"Length") << G4endl;
  
  
 // analysisManager->FillNtupleDColumn(0, position.z());
 // analysisManager->AddNtupleRow();
  //analysisManager->FillH1(4,position.z());
  
  

  //EventNumber->SetF_Number(icollnum);

  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime);
  
  G4ParticleDefinition* particleDefinition
           = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
           
 G4int num_Ev = G4Poisson(0.3);
// G4int num_Ev = 0;
 G4int num_Mu = G4Poisson(7);
// G4int num = 0;
 //G4int num_Mu = 1;
 // I need 0.3
//  if(num_Ev!=0) G4cout << " POISSON EV: " << G4endl;
 //G4cout << "POISSON MU: " << num_Mu << " POISSON EV: " << num_Ev<< G4endl;
 
           
  for (G4int i = 0; i < num_Ev; i++){
  
  G4double Z0 = G4RandGauss::shoot(0,30*cm);
 
  position = G4ThreeVector(0.0,0.0, Z0);
  
  if(num_Mu !=0) G4cout << " POISSON EV: " << G4endl;
  
  for(G4int j = 0; j < num_Mu; j++){
  
 
 
  
  G4PrimaryParticle* p = new G4PrimaryParticle(particleDefinition);
 
 // analysisManager->FillH1(0, ptime);
  //analysisManager->FillNtupleDColumn(0, ptime);
  G4double tetta = G4UniformRand()*2*3.14159;
 
  G4double phi = G4UniformRand()*2*3.14159;
 // G4double pt = G4UniformRand() + 0.1;
 // G4double px = pt*cos(phi);
  G4double px = std::sin(tetta)*std::cos(phi);
  G4double py = std::sin(tetta)*std::sin(phi);
  G4double pz = std::cos(tetta);
  G4ThreeVector direction = G4ThreeVector(px,py,pz);
 // G4cout << " Primary_Vertex = " << G4endl;
 // G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  p->SetMomentumDirection(direction);
  p->SetKineticEnergy(1*GeV);
 // num++;
  vertex->SetPrimary(p);
 //  G4cout << " just tetta = " << tetta*180/3.14159  << G4endl;
 /* if(tetta>3.14159/2 && tetta<3*3.14159/2) G4cout << " tetta1 = " << std::fabs(tetta*180/3.14159-180) << G4endl;
  else if(tetta>3*3.14159/2 && tetta<2*3.14159) G4cout << " tetta2 = " << 2*180-tetta*180/3.14159 << G4endl;
  else G4cout << " tetta3 = " << tetta*180/3.14159 << G4endl;//G4cout << "Number = " << num << G4endl;*/
  //analysisManager->FillH1(3,position.z());
  }
} 
  event->AddPrimaryVertex(vertex);
 
// }  
 }
  
 /* G4MapCache <G4double, G4int >::iterator it = Position_Map.Begin();
 
 for(G4int i = 0; it != Position_Map.End(); it++, i++){
 
     G4cout << i << ") Key - Z0 truly " << G4BestUnit(it->first,"Length") << ", number " << it->second << G4endl;
 } */
}  
  
  
 // }
  
  
