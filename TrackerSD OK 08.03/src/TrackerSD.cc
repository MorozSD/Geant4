#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "EventAction.hh"

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
  EventAction* Event = new EventAction;
  
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
  newHit->SetTime (aStep->GetPreStepPoint()->GetGlobalTime());
  Event->SetF_Time(aStep->GetPreStepPoint()->GetGlobalTime());

  fHitsCollection->insert( newHit );

  newHit->Print();
  
  
   G4cout << "AAAAAAAAAAAAAAAA" << Event->GetF_Time() << G4endl;
  
  //point0 - hit, point1 - tube
  G4double L = std::sqrt(std::pow((point0.x()-point1.x()),2)+std::pow((point0.y()-point1.y()),2));
  Event->SetF_Length(L);
  G4cout << "OOOOOOOOOOOOOOOOOOOOOOOOOOO" << Event->GetF_Length() << G4endl;

  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
   
 G4int nofHit = 1;
 G4cout
   << "-------->Hits Collection: in this event they are " << nofHit
   << " hits in the tracker chambers: " <<  G4endl;
  G4int nofHits = fHitsCollection->entries();
  
  for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();

  auto analysisManager = G4AnalysisManager::Instance();  
   for (G4int i=0;i<nofHits;i++)
     {
      G4double time = (*fHitsCollection)[i]->GetTime();
      G4double number = (*fHitsCollection)[i]->GetChamberNb();
      G4cout << "OOOOOOOOOOOOOOOOOOOOOOOOOOO" <<"Time = "<< G4BestUnit(time,"Time")<< "Number = " << number << G4endl;
      analysisManager->FillH1(0, time);
 }
 }
  

