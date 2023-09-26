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
 

 if (PDGE == 13) 
 {
 
  fHitsCollection->insert( newHit );
  G4cout << "TRACK_ID = " << Track_ID << " TYPE_OF_PARTICLE " << PDGE << G4endl;
  
  }
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
   
   G4int step = 0, k, ID_Number = 0, i_step = 0;
   G4int nextnumber, prevnumber, prevTrackID, nextTrackID, Number;
   G4ThreeVector pointH1, pointH2, pointH;
   G4double L, time_for_Track;
   
   G4int Array[50] {};
           
           
   for (G4int i=0;i<nofHits;i++)
     {
     
 //    G4cout << "i = " << i << G4endl;
        
  if (i==step)
  {
//  G4cout << "i after if = " << i << G4endl;
   nextnumber = 0;
   
      L = 0.;
      prevnumber = (*fHitsCollection)[i]->GetChamberNb();
      prevTrackID = (*fHitsCollection)[i]->GetTrackID();
 
      pointH1 = (*fHitsCollection)[i]->GetLocalPos();
      
   for (k = i+1; k<nofHits; k++)
     {
       nextnumber = (*fHitsCollection)[k]->GetChamberNb();
       nextTrackID = (*fHitsCollection)[k]->GetTrackID();
    //   G4cout << " PREV_ID = " << prevTrackID << " NEXT_ID = " << nextTrackID << G4endl;
        
       if (prevnumber!=nextnumber) 
       {
         pointH2 = (*fHitsCollection)[k-1]->GetLocalPos();
         step = k;
    //     G4cout << " STEP = " << step << G4endl;
         Array[i_step] = step;
         i_step++;
         G4int justTrackID = (*fHitsCollection)[step]->GetTrackID();
     //    G4cout << " STEP_TRACK_ID = " << justTrackID << G4endl;
         if (prevTrackID!=nextTrackID || k == nofHits-1)
         {
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
         
           for (G4int jj = 0; jj<ID_Number; jj++)
               {
            //     G4cout << "STEP_Number" << Array[jj]<< G4endl;
                 pointH = (*fHitsCollection)[Array[jj]]->GetPos();
                 time_for_Track = (*fHitsCollection)[Array[jj]]->GetTime();
                  G4double N = 0.;		// My strange attempt to solve this problem
                  G4int j = 0;
                  while (N>=0)
                   {
                     N = time_for_Track - j*76*ns;
                     j++;
                   }
                   
                 time_for_Track = N+76*ns;  
                 function_X_T[0][jj] = pointH.x();
                 G4cout << " X(T) = " << function_X_T[0][jj] << G4endl; 
                 function_Y_T[0][jj] = pointH.y();
                 G4cout << " Y(T) = " << function_Y_T[0][jj] << G4endl;
                 function_Z_T[0][jj] = pointH.z();
                 G4cout << " Z(T) = " << function_Z_T[0][jj] << G4endl;
                 function_X_T[1][jj] = time_for_Track;
                 G4cout << " T = " << function_X_T[1][jj] << G4endl;
                 function_Y_T[1][jj] = time_for_Track;
                 function_Z_T[1][jj] = time_for_Track;
               }      
           
          G4cout << "ID_Number = " << ID_Number << G4endl;
          

    G4cout << "Polynomial fit!" << G4endl;
    
    // Input values
    // **************************************************************
    size_t k = 3;                                    // Polynomial order
    G4bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    G4int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
    G4double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    G4double alphaval = 0.05;                          // Critical apha value

    /*G4double x[] = {0., 0.4, 1.0, 2.0, 4.0, 6.0};
    G4double y[] = {0., 0.16, 1.0, 4.0, 16.0, 36.};
    G4double erry[] = {0.1, 0.3, 0.2, 0.4, 0.1, 0.3}; */      // Data points (err on y) (if applicable)
    G4double x[ID_Number];
    G4double y[ID_Number];
    G4double erry[ID_Number]; 
    
    for (G4int ii = 0; ii<ID_Number; ii++){
     x[ii] = function_X_T[1][ii]; //arg t
     G4cout << "x = "<< x[ii] << G4endl;
     y[ii] = std::fabs(function_X_T[0][ii]); // function x(t)
     G4cout << "y = "  << y[ii] << G4endl;
     erry[ii] = 0.;
     }
    
    
    
    
    
    // Definition of other variables
    // **************************************************************
    size_t n = 0;                                    // Number of data points (adjusted later)
    size_t nstar = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    G4double coefbeta[k+1];                            // Coefficients of the polynomial
    G4double serbeta[k+1];                             // Standard error on coefficients
    G4double tstudentval = 0.;                         // Student t value
    G4double SE = 0.;                                  // Standard error
    
    G4double **XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
    G4double **Weights;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
    n = sizeof(x)/sizeof(G4double);
    nstar = n-1;
    if (fixedinter) nstar = n;
           
 /*   cout << "Number of points: " << n << endl;
    cout << "Polynomial order: " << k << endl;
    if (fixedinter) {
        cout << "A0 is fixed!" << endl;
    } else {
        cout << "A0 is adjustable!" << endl;
    }

    if (k>nstar) {
        cout << "The polynomial order is too high. Max should be " << n << " for adjustable A0 ";
        cout << "and " << n-1 << " for fixed A0. ";  
        cout << "Program stopped" << endl;
        return -1;
    }

    if (k==nstar) {
        cout << "The degree of freedom is equal to the number of points. ";
        cout << "The fit will be exact." << endl;  
    }*/

    XTWXInv = Make2DArray(k+1,k+1);
    Weights = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
    CalculateWeights(erry, Weights, n, wtype);
    
  //  cout << "Weights" << endl;
   // displayMat(Weights,n,n);

 /*   if (determinant(Weights,n)==0.) {
        G4cout << "One or more points have 0 error. Review the errors on points or use no weighting. ";
        G4cout << "Program stopped" << G4endl;
        return -1;
    } */

    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit(x,y,n,k,fixedinter,fixedinterval,coefbeta,Weights,XTWXInv);


    // Calculate related values
    // **************************************************************
    G4double RSS = CalculateRSS(x,y,coefbeta,Weights,std::fixed,n,k+1);
    G4double TSS = CalculateTSS(x,y,coefbeta,Weights,fixedinter,n,k+1);
    G4double R2 = CalculateR2COD(x,y,coefbeta,Weights,fixedinter,n,k+1);
    G4double R2Adj = CalculateR2Adj(x,y,coefbeta,Weights,fixedinter,n,k+1);

    if ((nstar-k)>0) {
        SE = std::sqrt(RSS/(nstar-k)); 
        tstudentval = std::fabs(CalculateTValueStudent(nstar-k, 1.-0.5*alphaval)); 
    }
    G4cout << "t-student value: " << tstudentval << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE,k,serbeta,XTWXInv);

    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k, nstar, tstudentval, coefbeta, serbeta);

    // Display statistics
    // **************************************************************
//    DisplayStatistics(n,nstar,k,RSS,R2,R2Adj,SE);
  
    // Display ANOVA table
    // **************************************************************
 //   DisplayANOVA(nstar, k, TSS, RSS);

    // Write the prediction and confidence intervals
    // **************************************************************
  //  WriteCIBands("CIBands2.dat",x,coefbeta,XTWXInv,tstudentval,SE,n,k);

    // Display the covariance and correlation matrix
    // **************************************************************
   // DisplayCovCorrMatrix(k, SE, fixedinter, XTWXInv);




       
          
           for (G4int ii = 0; ii<2; ii++)
     {
       delete [] function_X_T[ii];
     } 
     
     for (G4int ii = 0; ii<2; ii++)
     {
       delete [] function_Y_T[ii];
     }
     
      for (G4int ii = 0; ii<2; ii++)
     {
       delete [] function_Z_T[ii];
     }
    
          
          
           ID_Number = 0;  
           i_step = 0;   
         }
        else 
         {
         ID_Number++;
        // G4cout << "Id_number = " << ID_Number << G4endl;
         }
        
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
    //   G4cout << "GOOD" << G4endl;
      }
      else
      {
      L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
  //    G4cout << "WARNING" << G4endl; 
      }
             
  //    G4cout << " LH = " << G4BestUnit(L,"Length") << G4endl;
     
      if (L>0) 
      {
    
  //   G4double time = (*fHitsCollection)[i]->GetTime();
     auto analysisManager = G4AnalysisManager::Instance();
     
 
// if(time<0) G4cout << "WARNING!!!!!" << G4endl;
 // My strange attempt to solve this problem
     /* G4double N = 0.;
      G4int j = 0;
      while (N>=0)
      {
        N = time - j*76*ns;
        j++;
      }*/
  //    G4int potok = j-1;
       G4double TIME = 0.;
       G4double time = (*fHitsCollection)[step]->GetTime();
       analysisManager->FillH1(0, time);
      TIME=time+2.71012+1.21564*L+6.82868*L*L;
     analysisManager->FillH1(1, L);
     
  //   anaysislysisManager->FillH1(j, TIME);   // number of interaction - 1
  
     analysisManager->FillH1(2,TIME);
      
      }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)y2-(y2-y1)x2|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      
    }
    }  
    
   
  }
 
  


