#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "EventAction.hh"
#include <fstream>
#include <iostream>
#include <string>
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

//TrackerSD::TrackerSD(EventAction* eventAction)
//{
//  this->eventAction = eventAction;
//}

//TrackerSD tSD (evAc);
//tSD.eventAction->AddAbs(123);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection
    = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
  
  G4cout << "Инициализация" << G4endl;


 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory* )
{

  //EventAction* fEventAction = acum;
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false; 

  TrackerHit* newHit = new TrackerHit();
 // EventAction* Event = new EventAction;
  
  //static EventAction* acum = new EventAction();
 // TrackerHit* newHit_Edep = new TrackerHit(acum);
 // newHit_Edep->acum->AddAbs(edep);
 
 // acum->AddAbs(edep);
  
 // G4cout << "Energy= " << G4BestUnit(acum->getEdep(),"Energy") << G4endl;

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
//  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
//                                               ->GetCopyNumber());
//  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()->GetVolume()
//                                               ->GetCopyNo());
//Number of tube is a number of mother volume to PET, so GetVolume(1)
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                               ->GetCopyNo());
  G4ThreeVector point1 = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                               ->GetObjectTranslation();
  newHit->SetEdep(edep);
  
  G4ThreeVector point0 = aStep->GetPostStepPoint()->GetPosition();
  newHit->SetPos (point0);
  newHit->SetPosChamber (point1);
  newHit->SetTime (aStep->GetPreStepPoint()->GetGlobalTime());
 // Event->SetF_Time(aStep->GetPreStepPoint()->GetGlobalTime());

  fHitsCollection->insert( newHit );

  newHit->Print();
  
  
  // G4cout << "AAAAAAAAAAAAAAAA" << Event->GetF_Time() << G4endl;
  
  //point0 - hit, point1 - tube
  G4double L = std::sqrt(std::pow((point0.x()-point1.x()),2)+std::pow((point0.y()-point1.y()),2));
//  Event->SetF_Length(L);
  newHit->SetLength(L);
  G4cout << "OOOOOOOOOOOOOOOOOOOOOOOOOOO" << newHit->GetLength() << G4endl;

  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......'

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
   
 G4int nofHit = 1;
 G4cout
   << "-------->Hits Collection: in this event they are " << nofHit
   << " hits in the tracker chambers: " <<  G4endl;
  G4int nofHits = fHitsCollection->entries();
  
//  std::ifstream infile;
 // infile.open("/home/sveta/Project/B3Test-build/Print.txt");
  
  
 // for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();

   G4double A, B, L ;
   G4ThreeVector pointC, pointH1, pointH2;
   for (G4int i=0;i<nofHits;i++)
     {
     A = 0., B = 0., L = 0.;
      G4double time = (*fHitsCollection)[i]->GetTime();
      G4int nextnumber = 0;
      
      G4int prevnumber = (*fHitsCollection)[i]->GetChamberNb();
      if (i+1 != nofHits)
       {
        nextnumber = (*fHitsCollection)[i+1]->GetChamberNb();
        pointC = (*fHitsCollection)[i]->GetPosChamber();
        pointH1 = (*fHitsCollection)[i]->GetPos();
        B = std::sqrt(std::pow((pointC.x()-pointH1.x()),2)+std::pow((pointC.y()-pointH1.y()),2)+std::pow((pointC.z()-pointH1.z()),2));
        if(prevnumber == nextnumber)
         {
           pointH2 = (*fHitsCollection)[i+1]->GetPos();
           A = 0.5*std::sqrt(std::pow((pointH2.x()-pointH1.x()),2)+std::pow((pointH2.y()-pointH1.y()),2)+std::pow((pointH2.z()-pointH1.z()),2));
           L = std::sqrt(std::pow(B,2) - std::pow(A,2));
         }
       }
     // G4double len = (*fHitsCollection)[i]->GetLength();
     if(A==0) G4cout << "LENGTH of 1 HIT = " << G4BestUnit(B,"Length") << G4endl;
     else G4cout << "LENGTH of 2 HITS = " << G4BestUnit(L,"Length") << G4endl;
      G4cout << "OOOOOOOOOOOOOOOOOOOOOOOOOOO" <<"Time = "<< G4BestUnit(time,"Time") << "CHAMBER NUMBER" << prevnumber << G4endl;
      
      auto analysisManager = G4AnalysisManager::Instance();  
      analysisManager->FillH1(0, time);
     // analysisManager->FillH1(0, len);
      
    //  infile >> time;
      
      
 }
// infile.close();
 }
  

