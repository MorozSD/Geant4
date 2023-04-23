#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "EventAction.hh"
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

   
   G4int step = 0, k;
   G4int nextnumber, prevnumber;
   G4ThreeVector pointC, pointH1, pointH2;
   for (G4int i=0;i<nofHits;i++)
     {
     G4double time = (*fHitsCollection)[i]->GetTime();
     auto analysisManager = G4AnalysisManager::Instance();
     analysisManager->FillH1(1, time);
    // }
   /*
      G4double time = (*fHitsCollection)[i]->GetTime();
      G4double N = 0.;
      G4int j = 0;
      while (N>=0)
      {
        N = time - j*76*ns;
        j++;
      }
      G4int potok = j-1;*/
 // if (potok == 1) analysisManager->FillH1(w2, time);
 // if (potok == 2) analysisManager->FillH1(1, time);
 // if (potok == 3) analysisManager->FillH1(2, time);
 // if (potok == 4) analysisManager->FillH1(3, time);
//  if (potok == 5) analysisManager->FillH1(4, time);
//  if (potok == 6) analysisManager->FillH1(5, time);
 // if (potok == 7) analysisManager->FillH1(6, time);
//  if (potok == 8) analysisManager->FillH1(7, time);
//  if (potok == 9) analysisManager->FillH1(8, time);
 // if (potok == 10) analysisManager->FillH1(9, time);
//  if (potok == 11) analysisManager->FillH1(10, time);
//  if (potok == 12) analysisManager->FillH1(11, time);
//  if (potok == 13) analysisManager->FillH1(12, time);
 
  if (i==step)
  {
   nextnumber = 0;
   
      G4double L1 = 0., L2 = 0.;
      prevnumber = (*fHitsCollection)[i]->GetChamberNb();
      pointC = (*fHitsCollection)[i]->GetPosChamber();
      pointH1 = (*fHitsCollection)[i]->GetPos();
      G4double TIME = 0.;
   for (k = i+1; k<nofHits; k++)
     {
       nextnumber = (*fHitsCollection)[k]->GetChamberNb();
      // G4cout << " PREV = " << prevnumber << " NEXT = " << nextnumber<< G4endl;
       if (prevnumber!=nextnumber) 
       {
         pointH2 = (*fHitsCollection)[k-1]->GetPos();
         step = k;
     //    G4cout << " STEP = " << step << G4endl;
         break; 
       }
     //  G4cout << " K = " << k << G4endl;
     }
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
     
      if(pointH1==pointH2)
      {
       L1 = std::sqrt(std::pow(pointH1.x()-pointC.x(),2)+std::pow(pointH1.y()-pointC.y(),2)); //формула для плоскости
       
    //  G4cout << " LH1 = " << G4BestUnit(L1,"Length") << G4endl;
       if (L1!=0) 
       {
       //Field 1, non field 2
     // TIME=5.03009 - 0.0405048*L1 + 0.432331*L1*L1+time;
      TIME=5.01783 + 0.00976921*L1 + 0.33728*L1*L1+time; 
      analysisManager->FillH1(3, TIME);
      analysisManager->FillH1(2, L1);
       }
      }
      else
      {
    //  G4cout <<" POS1 X = " <<pointH1.x()<<" Y = "<< pointH1.y() <<  " Z = "<<pointH1.z()<< " POS2 X = " << pointH2.x()<<" Y = "<< pointH2.y() <<  " Z = "<<pointH2.z() <<  " POSC X = " <<pointC.x()<<" Y = "<< pointC.y() <<  " Z = "<<pointC.z()<< G4endl;
      L2 = std::abs(((pointH2.x()-pointH1.x())*(pointH2.y()-pointC.y())-(pointH2.y()-pointH1.y())*(pointH2.x()-pointC.x())))/std::sqrt(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2));          //фформула скрещивающихся
     // G4cout << " LH2 = " << G4BestUnit(L2,"Length") << G4endl;
      if (L2!=0) 
      {
       //Field 1, non field 2
      //TIME=5.03009 - 0.0405048*L2 + 0.432331*L2*L2+time;
      TIME=5.01783 + 0.00976921*L1 + 0.33728*L1*L1+time; 
      analysisManager->FillH1(3, TIME);
      analysisManager->FillH1(2, L2);
      }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)(y2-y)-(y2-y1)(x2-x)|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      }
    }  
  }
 }
  


