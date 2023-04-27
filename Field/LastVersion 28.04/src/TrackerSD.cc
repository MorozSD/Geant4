#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}




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


  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false; 

  TrackerHit* newHit = new TrackerHit();


  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());

  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                               ->GetCopyNo());
  G4ThreeVector point1 = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                               ->GetObjectTranslation();
  newHit->SetEdep(edep);
  G4ThreeVector worldPosition = aStep->GetPreStepPoint()->GetPosition();

  G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();


  G4ThreeVector point0 = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  
  //G4ThreeVector localPosition = 

  newHit->SetPos (point0);
  newHit->SetPosChamber (point1);
  newHit->SetTime (aStep->GetPreStepPoint()->GetGlobalTime());
 

  fHitsCollection->insert( newHit );

 // newHit->Print();


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
   
   G4int step = 0, k;
   G4int nextnumber, prevnumber;
   G4ThreeVector pointC, pointH1, pointH2;
   for (G4int i=0;i<nofHits;i++)
     {
        
  if (i==step)
  {
   nextnumber = 0;
   
      G4double L1 = 0., L2 = 0.;
      prevnumber = (*fHitsCollection)[i]->GetChamberNb();
 
      pointH1 = (*fHitsCollection)[i]->GetPos();
      G4double TIME = 0.;
   for (k = i+1; k<nofHits; k++)
     {
       nextnumber = (*fHitsCollection)[k]->GetChamberNb();
    //   G4cout << " PREV = " << prevnumber << " NEXT = " << nextnumber<< G4endl;
       if (prevnumber!=nextnumber) 
       {
         pointH2 = (*fHitsCollection)[k-1]->GetPos();
         step = k;
         G4cout << " STEP = " << step << G4endl;
         break; 
       }
     //  G4cout << " K = " << k << G4endl;
     }
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
     
  
      
   //   G4cout <<" POS1 X = " <<pointH1.x()<<" Y = "<< pointH1.y() << " POS2 X = " << pointH2.x()<<" Y = "<< pointH2.y()  << G4endl;
    //  L2=std::sqrt(std::pow(pointH2.y(),2)+std::pow(pointH2.x(),2)); 
 //   G4cout << " SQRT = " << G4BestUnit(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2),"Length") << G4endl;
    if(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2)>0)
    {
    //фформула скрещивающихся
      L2 = std::abs(((pointH2.x()-pointH1.x())*(pointH2.y())-(pointH2.y()-pointH1.y())*(pointH2.x())))/std::sqrt(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2));   
      }
     /* else
      {
      L2 = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2));  
     // Now L2 = 0, if we have only one hit in volume. So we don't write this L2
      }*/
             
  //    G4cout << " LH2 = " << G4BestUnit(L2,"Length") << G4endl;
      if (L2>0) 
      {
    
     G4double time = (*fHitsCollection)[i]->GetTime();
     auto analysisManager = G4AnalysisManager::Instance();
 analysisManager->FillH1(0, time);
 // My strange attempt to solve this problem
      G4double N = 0.;
      G4int j = 0;
      while (N>=0)
      {
        N = time - j*76*ns;
        j++;
      }
      G4int potok = j-1;
       G4double times = (*fHitsCollection)[step]->GetTime();
      TIME=times+2.71012+1.21564*L2+6.82868*L2*L2;
     
     analysisManager->FillH1(1, L2);
  //   analysisManager->FillH1(2, TIME);
    
  if (potok == 1) analysisManager->FillH1(2, TIME);
  if (potok == 2) analysisManager->FillH1(3, TIME);
  if (potok == 3) analysisManager->FillH1(4, TIME);
  if (potok == 4) analysisManager->FillH1(5, TIME);
  if (potok == 5) analysisManager->FillH1(6, TIME);
  if (potok == 6) analysisManager->FillH1(7, TIME);
  if (potok == 7) analysisManager->FillH1(8, TIME);
  if (potok == 8) analysisManager->FillH1(9, TIME);
  if (potok == 9) analysisManager->FillH1(10, TIME);
  if (potok == 10) analysisManager->FillH1(11, TIME);
  if (potok == 11) analysisManager->FillH1(12, TIME);
  if (potok == 12) analysisManager->FillH1(13, TIME);
  if (potok == 13) analysisManager->FillH1(14, TIME);
      
      
      }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)y2-(y2-y1)x2|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      
    }
    }  
  }
 
  


