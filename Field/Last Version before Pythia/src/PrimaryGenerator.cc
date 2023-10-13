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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//{  
PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
  EventAction* EventNumber = new EventAction;
 // G4double Z0 =  G4UniformRand()*200*cm;
 // G4ThreeVector position = G4ThreeVector(0.0,0.0,-1.0*cm);
  
  G4double steptime = 76*ns;
  G4ThreeVector position;
  
  auto analysisManager = G4AnalysisManager::Instance(); 
 // G4String fileName = "TestPos.root";
 // analysisManager->OpenFile(fileName);
      
  for (icollnum = 0; icollnum < 13; icollnum++){
 // for (icollnum = 0; icollnum < 1; icollnum++){
  G4double ptime = steptime*icollnum;
//  G4double Z0 =  G4UniformRand()*20*cm-10*cm;
  position = G4ThreeVector(0.0,0.0, 0.0);
  
 // G4cout << " POSITION = " << G4BestUnit(position.z(),"Length") << G4endl;
 // analysisManager->FillH1(0, position.z());
  
 // analysisManager->FillNtupleDColumn(0, position.z());
 // analysisManager->AddNtupleRow();
  //analysisManager->FillH1(4,position.z());
  
  

  EventNumber->SetF_Number(icollnum);

  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime);
  
  G4ParticleDefinition* particleDefinition
           = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
           
 G4int iptclnum = G4Poisson(0.3);
 // I need 0.3
 //G4cout << "POISSON: " << iptclnum << G4endl;
 
           
  for (G4int iptcl = 0; iptcl < iptclnum; iptcl++){
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
  
  p->SetMomentumDirection(direction);
  p->SetKineticEnergy(1*GeV);
  
  vertex->SetPrimary(p);
 }
 event->AddPrimaryVertex(vertex);  
 }
  
  
  }
