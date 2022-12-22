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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B3::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PrimaryParticle.hh"

namespace B3
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  
 // G4int n_particle = 1;
 G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
                    = particleTable->FindParticle("proton");
 G4PrimaryParticle* p =  new G4PrimaryParticle(particle);
 
 PrimaryGeneratorAction::~PrimaryGeneratorAction()
 {
 delete p;
 delete vertex;
 }
 void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 
 G4double energy = G4UniformRand()*1*GeV + 0.1*GeV; 
 p->SetKineticEnergy( energy );
 
 G4double phi = G4UniformRand()*2*CLHEP::pi;
 G4double teta = G4UniformRand()*2*CLHEP::pi;
 
 G4double px = -std::sin(teta)*std::cos(phi);
 G4double py = -std::sin(teta)*std::sin(phi);
 G4double pz = -std::cos(teta);
 
 p->SetMomentumDirection( G4ThreeVector(px,py,pz) );

 p->SetPolarization(G4ThreeVector(px,py,pz));
 G4double time = 10*ns;
 //p->SetTime( time );
 G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(0.,0.,0.), time);
 //void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
// {
 vertex->SetPrimary( p );
 anEvent->AddPrimaryVertex( vertex );
 }
