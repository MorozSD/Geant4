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
#include "PolyFit.hh"
#include "G4RootAnalysisReader.hh" 


//#include <iostream>

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
  //G4int PDGE = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4ThreeVector Pos = aStep->GetTrack()->GetVertexPosition();
  newHit->SetTrackID(Track_ID);
  
  newHit->SetTrackVertex(Pos);

  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                               ->GetCopyNo());
 // G4ThreeVector point1 = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1)
                                              // ->GetObjectTranslation();
  newHit->SetEdep(edep);
  G4ThreeVector worldPosition = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector point0 = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  
  newHit->SetLocalPos (point0);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
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
   
 G4int step = 0, k, ID_Number = 0, i_step = 0;
 G4int nextnumber, prevnumber, prevTrackID, nextTrackID, Number;
 G4ThreeVector pointH1, pointH2, pointH;
 G4double L, time_for_Track, time_for_Track_Global;
   
 G4int Array[50] {};
           
           
 for (G4int i=0;i<nofHits;i++)
    {
            
    if (i==step){
    
        nextnumber = 0;
        L = 0.;
       prevnumber = (*fHitsCollection)[i]->GetChamberNb();
       prevTrackID = (*fHitsCollection)[i]->GetTrackID();
       pointH1 = (*fHitsCollection)[i]->GetLocalPos();
      
       for (k = i+1; k<nofHits; k++){
        
          nextnumber = (*fHitsCollection)[k]->GetChamberNb();
          nextTrackID = (*fHitsCollection)[k]->GetTrackID();
           
/*          if (prevnumber!=nextnumber){ 
       
             pointH2 = (*fHitsCollection)[k-1]->GetLocalPos();
             step = k;
             Array[i_step] = step;
             i_step++;
             
             if (prevTrackID!=nextTrackID || k+1 > nofHits-1){
             
                G4double **function_X_T = new G4double* [2];
                G4double **function_Y_T = new G4double* [2];
                G4double **function_Z_T = new G4double* [2];        
                                 
                for (G4int ii = 0; ii<2; ii++)
                {
                  function_X_T[ii] = new G4double [ID_Number];    //Create array, i = 2, j = ID
           	}
           
                for (G4int ii = 0; ii<2; ii++)
             	{
                  function_Y_T[ii] = new G4double [ID_Number];    //Create array, i = 2, j = ID
           	}
           
               for (G4int ii = 0; ii<2; ii++)
               {
                  function_Z_T[ii] = new G4double [ID_Number];    //Create array, i = 2, j = ID
               }
         
              for (G4int jj = 0; jj<ID_Number; jj++){
               
                 pointH = (*fHitsCollection)[Array[jj]]->GetPos();
                 time_for_Track_Global = (*fHitsCollection)[Array[jj]]->GetTime();
                 G4double N = 0.;		// My strange attempt to solve this problem
                 G4int j = 0;
                 while (N>=0)
                    {
                    N = time_for_Track_Global - j*76*ns;
                    j++;
                    }
                   
 
       ID_Number = 0;  
       i_step = 0;   
        }
        else 
        {
           ID_Number++;
        }
        
           break; 
        }
     
}*/
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
     
  
      
   //   G4cout <<" POS1 X = " <<pointH1.x()<<" Y = "<< pointH1.y() << " POS2 X = " << pointH2.x()<<" Y = "<< pointH2.y()  << G4endl;
    //  L2=std::sqrt(std::pow(pointH2.y(),2)+std::pow(pointH2.x(),2)); 
 //   G4cout << " SQRT = " << G4BestUnit(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2),"Length") << G4endl;
     if(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2)>0){
    //фформула скрещивающихся
         L = std::abs(((pointH2.x()-pointH1.x())*(pointH2.y())-(pointH2.y()-pointH1.y())*(pointH2.x())))/std::sqrt(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2));  
    //   G4cout << "GOOD" << G4endl;
     }
     else {
         L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
     }
             
  //    G4cout << " LH = " << G4BestUnit(L,"Length") << G4endl;
     
      if (L>0) 
      {
    
  //   G4double time = (*fHitsCollection)[i]->GetTime();
     
     
 

       auto analysisManager = G4AnalysisManager::Instance();
       G4double TIME = 0.;
       G4double time = (*fHitsCollection)[step]->GetTime();
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
      //j= number of interaction - 1, we start with hist 2 => j-1+2=j+1
       G4int potok = j+1;      
       
       
       TIME=time+2.71012+1.21564*L+6.82868*L*L;
       analysisManager->FillH1(1, L);
           
       analysisManager->FillH1(potok, TIME);   // number of interaction - 1
        }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)y2-(y2-y1)x2|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      
    }
    }  
    
   
  }
} 
  


