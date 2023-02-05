#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// namespace B3 
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
  //vertex A uniform on a cylinder
  //
 /// const G4double r = 2*mm;
//  const G4double zmax = 8*mm;
  //
  G4ThreeVector position = G4ThreeVector(0.0,0.0,0.0);
  G4double ptime = 0*s;
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime);
  
  G4ParticleDefinition* particleDefinition
           = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  G4PrimaryParticle* p = new G4PrimaryParticle(particleDefinition);
  
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
  event->AddPrimaryVertex(vertex);
  
  G4cout<< "000000000000000000000000000000000.\n" << "px = "<< px << "\n py = "<< py << "\n pz = "<< pz <<"  000000000000000000000000000000000.\n"<<G4endl;
  
  }
