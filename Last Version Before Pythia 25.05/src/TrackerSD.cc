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
  G4int Track_ID = aStep->GetTrack()->GetTrackID();
  G4int PDGE = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

  newHit->SetTrackID(Track_ID);

  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                               ->GetCopyNo());
 // G4ThreeVector point1 = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                              // ->GetObjectTranslation();
  newHit->SetEdep(edep);
  G4ThreeVector worldPosition = aStep->GetPreStepPoint()->GetPosition();

//  G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();


  G4ThreeVector point0 = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  

  newHit->SetLocalPos (point0);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetTime (aStep->GetPreStepPoint()->GetGlobalTime());
 

 /*if (PDGE == 13) 
 {*/
 
  fHitsCollection->insert( newHit );
//  G4cout << "TRACK_ID = " << Track_ID << " TYPE_OF_PARTICLE " << PDGE << G4endl;
  
 // }
 // fHitsCollection->insert( newHit );
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
   
   G4int step = 0, k, number_warnings=0, number_norm=0;
   G4int nextnumber, prevnumber;
   G4ThreeVector pointH1, pointH2;
   G4double L;
   auto analysisManager = G4AnalysisManager::Instance();
   for (G4int i=0;i<nofHits;i++)
     {
     if(i==nofHits-1) 
     {
     //Ratio = (number_warnings/(number_warnings+number_norm))*100;
     G4cout << "NUMBER OF WARNINGS = " << number_warnings << " NUMBER OF NORM " << number_norm << G4endl;//" RATIO =  " << Ratio  << G4endl;
     }
      
     
 //    G4cout << "i = " << i << G4endl;
        
  if (i==step)
  {
//  G4cout << "i after if = " << i << G4endl;
   nextnumber = 0;
   
      L = 0.;
      prevnumber = (*fHitsCollection)[i]->GetChamberNb();
 
      pointH1 = (*fHitsCollection)[i]->GetLocalPos();
      G4double TIME = 0.;
   for (k = i+1; k<nofHits; k++)
     {
       nextnumber = (*fHitsCollection)[k]->GetChamberNb();
    //   G4cout << " PREV = " << prevnumber << " NEXT = " << nextnumber<< G4endl;
       if (prevnumber!=nextnumber) 
       {
         pointH2 = (*fHitsCollection)[k-1]->GetLocalPos();
         step = k;
       //  G4cout << " STEP = " << step << G4endl;
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
      L = std::abs(((pointH2.x()-pointH1.x())*(pointH2.y())-(pointH2.y()-pointH1.y())*(pointH2.x())))/std::sqrt(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2)); 
      number_norm = 1*ns;
      analysisManager->FillH1(16, number_norm);
     
    //   G4cout << "GOOD" << G4endl;
      }
      else
      {
      L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
      number_warnings = 1*ns ;
      analysisManager->FillH1(15, number_warnings);
  //    G4cout << "WARNING" << G4endl; 
      }
             
  //    G4cout << " LH = " << G4BestUnit(L,"Length") << G4endl;
     
      if (L>0) 
      {
    
     G4double time = (*fHitsCollection)[i]->GetTime();
     
 analysisManager->FillH1(0, time);
 
// if(time<0) G4cout << "WARNING!!!!!" << G4endl;
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
      TIME=times+2.71012+1.21564*L+6.82868*L*L;
     
     analysisManager->FillH1(1, L);
    // analysisManager->FillH1(2, TIME);
    
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
 
  


