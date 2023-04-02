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
#include "G4SystemOfUnits.hh"

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

 // newHit->Print();
  
  
  // G4cout << "AAAAAAAAAAAAAAAA" << Event->GetF_Time() << G4endl;
  
  //point0 - hit, point1 - tube
 
//  Event->SetF_Length(L);
 // newHit->SetLength(L);
//  G4cout << " X = " <<point1.x()<<" Y = "<< point1.y() <<  " Z = "<<point1.z() << " EDEP = "<< G4BestUnit(edep,"Energy")  << " TRACKID = "<< newHit->GetTrackID()<< " NUMBER = " << newHit->GetChamberNb() << " POS X = " << point0.x()<<" Y = "<< point0.y() <<  " Z = "<<point0.z() << " TIME = " <<newHit->GetTime() <<G4endl;

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

   
 
   G4ThreeVector pointC, pointH1, pointH2;
   for (G4int i=0;i<nofHits;i++)
     {
    
   
      G4double time = (*fHitsCollection)[i]->GetTime();
      G4double N = 0.;
      G4int j = 0;
      while (N>=0)
      {
        N = time - j*76*ns;
        j++;
      }
      G4int potok = j-1;
   G4cout << "OOOOOOOOOOOOOOOOOOOOOOOOOOO" <<"Time = "<< G4BestUnit(time,"Time") << "POTOK NUMBER = " << potok << G4endl;   
   auto analysisManager = G4AnalysisManager::Instance();  
 //  analysisManager->FillH1(potok-1, time);
 //пыталась через switch case
  if (potok == 1) analysisManager->FillH1(0, time);
  if (potok == 2) analysisManager->FillH1(1, time);
  if (potok == 3) analysisManager->FillH1(2, time);
  if (potok == 4) analysisManager->FillH1(3, time);
  if (potok == 5) analysisManager->FillH1(4, time);
  if (potok == 6) analysisManager->FillH1(5, time);
  if (potok == 7) analysisManager->FillH1(6, time);
  if (potok == 8) analysisManager->FillH1(7, time);
  if (potok == 9) analysisManager->FillH1(8, time);
  if (potok == 10) analysisManager->FillH1(9, time);
  if (potok == 11) analysisManager->FillH1(10, time);
  if (potok == 12) analysisManager->FillH1(11, time);
  if (potok == 13) analysisManager->FillH1(12, time);
  }
  
  
 G4int ii = 0, k;
  G4int nextnumber, prevnumber;
  	
   while (ii < nofHits)
   {
     nextnumber = 0;
   
      G4double L1 = 0., L2 = 0.;
      prevnumber = (*fHitsCollection)[ii]->GetChamberNb();
      pointC = (*fHitsCollection)[ii]->GetPosChamber();
      pointH1 = (*fHitsCollection)[ii]->GetPos();
      
     for (k = ii; k<nofHits-1; k++)
     {
       nextnumber = (*fHitsCollection)[k]->GetChamberNb();
       
       if (prevnumber!=nextnumber) 
       {
         pointH2 = (*fHitsCollection)[k-1]->GetPos();
         break; 
       }
       G4cout << " K = " << k << G4endl;
     }
      if(pointH1==pointH2)
      {
    //  l1 = std::sqrt(pow(pointH1.x()-pointC.x(),2)+pow(pointH1.y()-pointC.y(),2)+pow(pointH1.z()-pointC.z(),2));
     // l2 = std::sqrt(pow(pointH2.x()-pointC.x(),2)+pow(pointH2.y()-pointC.y(),2)+pow(pointH2.z()-pointC.z(),2));
       L1 = std::sqrt(std::pow(pointH1.x()-pointC.x(),2)+std::pow(pointH1.y()-pointC.y(),2)); //формула для плоскости
       
       G4cout << " LH1 = " << L1 << G4endl;
      }
      else
      {
      L2 = ((pointH2.x()-pointH1.x())*(pointH2.y()-pointC.y())-(pointH2.y()-pointH1.y())*(pointH2.x()-pointC.x()))/std::sqrt(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2));          //фформула скрещивающихся
      G4cout << " LH2 = " << L2 << G4endl;
      }
      ii = k;
   }
      
 }
  

