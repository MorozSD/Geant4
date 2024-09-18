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
#include "G4Poisson.hh"
#include "G4AnalysisManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//{  
using namespace Pythia8;
//{  
PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator(),
 fPy8Gun(nullptr)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ 
 delete fPy8Gun;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
    
  G4double steptime = 76*ns;
  G4ThreeVector position;
  
  auto analysisManager = G4AnalysisManager::Instance(); 
  
  
      
  
  position = G4ThreeVector(0.0,0.0, 0.0);

          
 
 // I need 0.3

  
  fPy8Gun = new Pythia( "../share/Pythia8/xmldoc", false );
 // fPy8Gun->readString("ProcessLevel:all = off");
  fPy8Gun->readString("Beams:eCM = 27.");
  fPy8Gun->readString("SoftQCD:inelastic = on");
 // fPy8Gun->readString("PhaseSpace:pTHatMin = 20.");
  
  
  fPy8Gun->readString("Beams:idA = 2212");
  fPy8Gun->readString("Beams:idB = 2212");
  // specify how many Py8 events to print out, at either level
   // in this particular case print out a maximum of 10 events
   //
  fPy8Gun->readString("Next:numberShowProcess = 0" );
  fPy8Gun->readString("Next:numberShowEvent = 10");
  fPy8Gun->init();
  
  for (icollnum = 0; icollnum < 14; icollnum++){ 
 
 
  	G4int iptclnum = G4Poisson(0.3);
 
  	G4double Z0 = G4RandGauss::shoot(0,30*cm);
 	position = G4ThreeVector(0.0,0.0, Z0);

 	//G4cout << "POISSON: " << iptclnum << G4endl;
 	G4double ptime = steptime*icollnum;
 	G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime);
 
  	G4int nCharged = 0;
  	for (G4int iEvent = 0; iEvent < iptclnum; ++iEvent) {
  		if (!fPy8Gun->next()) continue;
       
          	for (G4int i = 0; i < fPy8Gun->event.size(); ++i){
        
        		if (fPy8Gun->event[i].isFinal()){
           
         		nCharged++;
             		//G4cout << " Number of particles = " << nCharged << G4endl;
             		G4int ID_P = fPy8Gun->event[i].id();
             		//G4cout << " ID of particle = " << ID_P << G4endl;
             		G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(ID_P);
             
  
             		G4PrimaryParticle* p = new G4PrimaryParticle(particleDefinition);
 
 		        G4double px = fPy8Gun->event[i].px();
            		G4double py = fPy8Gun->event[i].py();
            		G4double pz = fPy8Gun->event[i].pz();
            
            		G4double norm = std::sqrt(px*px+py*py+pz*pz);
            
             		G4ThreeVector direction = G4ThreeVector(px/norm,py/norm,pz/norm);
             		G4double En = fPy8Gun->event[i].e()-fPy8Gun->event[i].m();
             		p->SetMomentumDirection(direction);
  	     		p->SetKineticEnergy(En*GeV);
 	 
   			vertex->SetPrimary(p);        
}
             }
            }
            event->AddPrimaryVertex(vertex);  
            } 
         

  }
