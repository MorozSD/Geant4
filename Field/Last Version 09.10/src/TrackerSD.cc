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
//  G4cout << "TRACK_ID = " << Track_ID << " TYPE_OF_PARTICLE " << PDGE << G4endl;
  
  }
 // fHitsCollection->insert( newHit );
 // newHit->Print();
 


  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......'

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{

 using G4AnalysisReader = G4RootAnalysisReader;   
 
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
                 time_for_Track_Global = (*fHitsCollection)[Array[jj]]->GetTime();
                  G4double N = 0.;		// My strange attempt to solve this problem
                  G4int j = 0;
                  while (N>=0)
                   {
                     N = time_for_Track_Global - j*76*ns;
                     j++;
                   }
                   
                 time_for_Track = N+76*ns;  
                 function_X_T[0][jj] = pointH.x();
               //  G4cout << " X(T) = " << function_X_T[0][jj] << G4endl; 
                 function_Y_T[0][jj] = pointH.y();
               //  G4cout << " Y(T) = " << function_Y_T[0][jj] << G4endl;
                 function_Z_T[0][jj] = pointH.z();
                //function_Z_T[0][jj] = jj*jj+25;
                // G4cout << " Z(T) = " << function_Z_T[0][jj] << G4endl;
                 function_X_T[1][jj] = time_for_Track_Global;
              //   G4cout << " T = " << function_X_T[1][jj] << G4endl;
                 function_Y_T[1][jj] = time_for_Track;
                 function_Z_T[1][jj] = time_for_Track;
              // function_Z_T[1][jj] = jj;
               }      
           
          G4cout << "ID_Number = " << ID_Number << G4endl;
          size_t k = 2;                              // Polynomial order
          
    if(ID_Number > k){      

    G4cout << "Polynomial fit X(z)!" << G4endl;
    
    bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
    double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    double alphaval = 0.05;                          // Critical apha value

    // Input values
    // **************************************************************
    G4double xX[ID_Number];
    G4double yX[ID_Number];
    G4double erryX[ID_Number]; 
    
    for (G4int ii = 0; ii<ID_Number; ii++){
    // xX[ii] = std::fabs(function_Z_T[0][ii]); //arg t
    //xX[ii] = ii+0.2;
    xX[ii] = function_Z_T[0][ii];
     G4cout << "z = "<< xX[ii] << G4endl;
     yX[ii] = std::fabs(function_X_T[0][ii]); // function x(t)
  //  yX[ii] = ii*ii + 0.5;
    // yX[ii] = function_Y_T[0][ii];
     G4cout << "x = "  << yX[ii] << G4endl;
     erryX[ii] = 0.;
     }
        
    
    
    // Definition of other variables
    // **************************************************************
    size_t n = 0;                                    // Number of data points (adjusted later)
    size_t nstar = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    G4double coefbeta_X[k+1];                            // Coefficients of the polynomial
    G4double serbeta_X[k+1];                             // Standard error on coefficients
    G4double tstudentval_X = 0.;                         // Student t value
    G4double SE_X = 0.;                                  // Standard error
    
    G4double **XTWXInv_X;                                // Matrix XTWX Inverse [k+1,k+1]
    G4double **Weights_X;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
    n = sizeof(xX)/sizeof(G4double);
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

    XTWXInv_X = Make2DArray(k+1,k+1);
    Weights_X = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
    CalculateWeights(erryX, Weights_X, n, wtype);
    
  //  cout << "Weights" << endl;
   // displayMat(Weights,n,n);

 /*   if (determinant(Weights,n)==0.) {
        G4cout << "One or more points have 0 error. Review the errors on points or use no weighting. ";
        G4cout << "Program stopped" << G4endl;
        return -1;
    } */

    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit(xX,yX,n,k,fixedinter,fixedinterval,coefbeta_X,Weights_X,XTWXInv_X);


    // Calculate related values
    // **************************************************************
    G4double RSS_X = CalculateRSS(xX,yX,coefbeta_X,Weights_X,std::fixed,n,k+1);
    G4double TSS_X = CalculateTSS(xX,yX,coefbeta_X,Weights_X,fixedinter,n,k+1);
    G4double R2_X = CalculateR2COD(xX,yX,coefbeta_X,Weights_X,fixedinter,n,k+1);
    G4double R2Adj_X = CalculateR2Adj(xX,yX,coefbeta_X,Weights_X,fixedinter,n,k+1);

    if ((nstar-k)>0) {
        SE_X = std::sqrt(RSS_X/(nstar-k)); 
        tstudentval_X = std::fabs(CalculateTValueStudent(nstar-k, 1.-0.5*alphaval)); 
    }
    G4cout << "t-student value: " << tstudentval_X << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE_X,k,serbeta_X,XTWXInv_X);

    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k, nstar, tstudentval_X, coefbeta_X, serbeta_X);

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
		
    G4cout << "Polynomial fit Y(z)!" << G4endl;
	
    G4double xY[ID_Number];
    G4double yY[ID_Number];
    G4double erryY[ID_Number]; 
    G4cout << " Data: " << G4endl;
    for (G4int ii = 0; ii<ID_Number; ii++){
    xY[ii] = function_Z_T[0][ii]; //arg t
     
     yY[ii] = std::fabs(function_Y_T[0][ii]); // function x(t)
     G4cout << "z = "  << xY[ii] << G4endl;
     G4cout << "y = " << function_Y_T[0][ii] << G4endl;
     erryY[ii] = 0.;
     
  /*   G4cout << "{" << function_Z_T[1][ii] << ","<< "\t";
     //G4cout << "{" << function_Z_T[0][ii]<<"," << "\t";
     G4cout << yY[ii] <<"}"<< "\t";
     G4cout << G4endl;*/
     }
    
    G4double coefbeta_Y[k+1];                            // Coefficients of the polynomial
    G4double serbeta_Y[k+1];                             // Standard error on coefficients
    G4double tstudentval_Y = 0.;                         // Student t value
    G4double SE_Y = 0.;                                  // Standard error
    
    G4double **XTWXInv_Y;                                // Matrix XTWX Inverse [k+1,k+1]
    G4double **Weights_Y;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
 /*   n = sizeof(x)/sizeof(G4double);
    nstar = n-1;
    if (fixedinter) nstar = n;*/
           

    XTWXInv_Y = Make2DArray(k+1,k+1);
    Weights_Y = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
    CalculateWeights(erryY, Weights_Y, n, wtype);
    
    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit(xY,yY,n,k,fixedinter,fixedinterval,coefbeta_Y,Weights_Y,XTWXInv_Y);


    // Calculate related values
    // **************************************************************
    G4double RSS_Y = CalculateRSS(xY,yY,coefbeta_Y,Weights_Y,std::fixed,n,k+1);
    G4double TSS_Y = CalculateTSS(xY,yY,coefbeta_Y,Weights_Y,fixedinter,n,k+1);
    G4double R2_Y = CalculateR2COD(xY,yY,coefbeta_Y,Weights_Y,fixedinter,n,k+1);
    G4double R2Adj_Y = CalculateR2Adj(xY,yY,coefbeta_Y,Weights_Y,fixedinter,n,k+1);

    if ((nstar-k)>0) {
        SE_Y = std::sqrt(RSS_Y/(nstar-k)); 
        tstudentval_Y = std::fabs(CalculateTValueStudent(nstar-k, 1.-0.5*alphaval)); 
    }
    G4cout << "t-student value: " << tstudentval_Y << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE_Y,k,serbeta_Y,XTWXInv_Y);

    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k, nstar, tstudentval_Y, coefbeta_Y, serbeta_Y);

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
	auto analysisManager = G4AnalysisManager::Instance();	
		
	G4double Zx = SetStartingPoint(k, coefbeta_X, ID_Number, xX);
	
	//G4cout << " Z0 by X(z) = " << Zx << G4endl;
     //   analysisManager->FillH1(0, Zx);

	G4double Zy = SetStartingPoint(k, coefbeta_Y, ID_Number, xY);
	
	//G4cout << " Z0 by Y(z) = " << Zy << G4endl;
	//analysisManager->FillH1(1, Zy);
 }	

	auto analysisReader = G4AnalysisReader::Instance();
	analysisReader->SetVerboseLevel(1);
          
        // Define a base file name
     /*   analysisReader->SetFileName("Test.root");

        // Read ntuple
        G4int ntupleId = analysisReader->GetNtuple("Test");
        G4int counter = 0;
        if ( ntupleId >= 0 ) {
        G4double Z0;
        analysisReader->SetNtupleDColumn("StartPos", Z0);
        G4cout << "Ntuple Start Position, reading selected column StartPos" << G4endl;
        while ( analysisReader->GetNtupleRow() ) {
          G4cout << counter++ << "th entry: " << " Start Position Z0: " << Z0 << G4endl;
        }
}  */
          
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
       auto analysisManager = G4AnalysisManager::Instance();
       G4double TIME = 0.;
       G4double time = (*fHitsCollection)[step]->GetTime();
    //   analysisManager->FillH1(2, time);
       //analysisManager->FillNtupleDColumn(2, time);
       TIME=time+2.71012+1.21564*L+6.82868*L*L;
   //    analysisManager->FillH1(1, L);
     //analysisManager->FillNtupleDColumn(1, L);
     
  //   anaysislysisManager->FillH1(j, TIME);   // number of interaction - 1
  
    // analysisManager->FillH1(3,TIME);
     //analysisManager->FillNtupleDColumn(3, TIME); 
   //  analysisManager->AddNtupleRow();
      }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)y2-(y2-y1)x2|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      
    }
    }  
    
   
  }
 
  


