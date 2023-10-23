#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "PolyFit.hh"
#include "G4RootAnalysisReader.hh" 
//#include "G4RandGauss.hh"


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
  G4int status = aStep->GetTrack()->GetParentID();
  
 // if (edep==0.||status != 0) return false; 
 // G4cout << " Edep = " << edep << G4endl;
 if (status != 0) return false;
  TrackerHit* newHit = new TrackerHit();
  G4int Track_ID = aStep->GetTrack()->GetTrackID();
  G4int PDGE = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  
//  G4cout << " PDGE = " << PDGE << "My TrackID = " << Track_ID << " Parent TrackID = " << status<< G4endl;
  
 // newHit->SetEdep(edep);
  
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
  G4double time0 = aStep->GetPreStepPoint()->GetGlobalTime();
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetTime (time0);
  
  newHit->SetTrStat(aStep->GetTrack()->GetParentID());
  
  G4double Ek = aStep->GetTrack()->GetKineticEnergy();
  
  
  auto analysisManager = G4AnalysisManager::Instance();
 // analysisManager->FillH1(0, Ek);
          
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
   
 G4int step = 0, k, ID_Number = 0, i_step = 0, k_min = 0, k_pred = 0;
 G4int nextnumber, prevnumber, prevTrackID, nextTrackID, Number, status;
 G4ThreeVector pointH0, pointH1, pointH;
 G4double L, L_min, L_pred, time_for_Track_Global, E_dep, Popal=0., Nepopal = 0.;
   
 G4int Array[100] {};
       
G4cout << nofHits << G4endl;           
 for (G4int i=0;i<nofHits;i++)
    {
        
    if (i==step){
 
//  G4cout << " step " << step << G4endl;
       nextnumber = 0;
            
       prevnumber = (*fHitsCollection)[i]->GetChamberNb();
       
       prevTrackID = (*fHitsCollection)[i]->GetTrackID();
      // G4cout << " PREV = " << prevTrackID<<G4endl;
       pointH0 = (*fHitsCollection)[i]->GetLocalPos();
     //  L = std::sqrt(std::pow(pointH0.y(),2)+std::pow(pointH0.x(),2)); 
     /*  
       k_min = i;
       L_min = L;
       L_pred = L;
       k_pred = i;
     */
   //  G4cout << " L_nach =" << L_min << " k_nach = " << k_min << G4endl; 
       for (k = i+1; k<nofHits; k++){
         
    
           
    //   if(status)
        //  G4cout << "Parent ID" << status << G4endl;
         // G4int IDTr = (*fHitsCollection)[step]->GetTrackID();
        //  G4cout << "MY ID" << IDTr << G4endl;     
         
          
          nextnumber = (*fHitsCollection)[k]->GetChamberNb();
          nextTrackID = (*fHitsCollection)[k]->GetTrackID();
         // G4cout << " NEXT = " << nextTrackID<<G4endl;
     //     pointH1 = (*fHitsCollection)[k]->GetLocalPos();
      //   L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
      //  G4cout << " PREV = " << prevnumber << " NEXT = " << nextnumber <</* " L = " << L << " k = " << k <<*/ G4endl;
      //  G4cout << " PR ID = " << prevTrackID << " NX ID = " << nextTrackID << G4endl;  
      /*   
          if(L_min > L) { 
          
            L_min = L;
            k_min = k;
          }  
      */      
           
       
          if (prevnumber!=nextnumber){ 
      
             pointH1 = (*fHitsCollection)[k-1]->GetLocalPos();
             
            
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
     
  
      
   //   G4cout <<" POS1 X = " <<pointH1.x()<<" Y = "<< pointH1.y() << " POS2 X = " << pointH2.x()<<" Y = "<< pointH2.y()  << G4endl;
    //  L2=std::sqrt(std::pow(pointH2.y(),2)+std::pow(pointH2.x(),2)); 
 //   G4cout << " SQRT = " << G4BestUnit(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2),"Length") << G4endl;
     if(std::pow(pointH1.y()-pointH0.y(),2)+std::pow(pointH1.x()-pointH0.x(),2)>0){
    //фформула скрещивающихся
         L = std::fabs(((pointH1.x()-pointH0.x())*(pointH1.y())-(pointH1.y()-pointH0.y())*(pointH1.x())))/std::sqrt(std::pow(pointH1.y()-pointH0.y(),2)+std::pow(pointH1.x()-pointH0.x(),2));  
        }
     else {
         L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
           }
          
             
          Array[i_step] = k-1; // it is wrong
      //    G4cout << " Arr = "  <<Array[i_step]<<G4endl;
          i_step++;
             
          ID_Number++;
        
       
          G4double TIME = 0.;
    //   G4int status = (*fHitsCollection)[step]->GetTrStat();
    //   if(status)
     //  G4cout << "Parent ID" << status << G4endl;
      // G4int IDTr = (*fHitsCollection)[step]->GetTrackID();
      // G4cout << "MY ID" << IDTr << G4endl;
    //  if (status == 1){
       //G4double time = (*fHitsCollection)[k_pred]->GetTime();
       G4double time2 = (*fHitsCollection)[k-1]->GetTime();
       G4double time1 = (*fHitsCollection)[step]->GetTime();
       G4double time = (time1 + time2)/2;
       
       auto analysisManager = G4AnalysisManager::Instance();
       
       analysisManager->FillH1(0, time);
       //analysisManager->FillH1(1, status);
  
       // if(time<0) G4cout << "WARNING!!!!!" << G4endl;
       // My strange attempt to solve this problem
   /*   
       G4double N = 0.;
       G4int j = 0;
       while (N>=0)
       {
         N = time - j*76*ns;
         j++;
        }
      //j= number of interaction - 1, we start with hist 2 => j-1+2=j+1
       G4int potok = j+1;      
    */  
       
    /*   TIME=time+2.71012+1.21564*L_pred+6.82868*L_pred*L_pred;
       analysisManager->FillH1(1, L_pred);  
    //  G4cout << " potok = " << potok << G4endl;     
    */
      // TIME=time+2.71012+1.21564*L+6.82868*L*L;
       
       TIME = G4RandGauss::shoot(time+2.71012+1.21564*L+6.82868*L*L,20*ns);
       analysisManager->FillH1(1, L);
       analysisManager->FillH1(2, TIME);
          // number of interaction - 1
     // status   }
    // if L>0  }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)y2-(y2-y1)x2|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      
        step = k;
                     
             
        // block with fit     
             if (prevTrackID!=nextTrackID || k+1 == nofHits){
             G4cout << " ID_Number = " << ID_Number << G4endl;
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
                  time_for_Track_Global = (*fHitsCollection)[Array[jj]]->GetTime();
                  pointH = (*fHitsCollection)[Array[jj]]->GetPos();
                  
                  function_X_T[0][jj] = pointH.x();
                  function_Y_T[0][jj] = pointH.y();
                  function_Z_T[0][jj] = pointH.z();
                  function_X_T[1][jj] = time_for_Track_Global;
                  function_Y_T[1][jj] = time_for_Track_Global;
                  function_Z_T[1][jj] = time_for_Track_Global;
                  
                  G4cout << " x = " << function_X_T[0][jj] << G4endl;
                }  
                 
               size_t k = 2;                            // Polynomial order
               if(ID_Number > k){      

                  G4cout << "Polynomial fit Y(x)!" << G4endl;
    
                  G4bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
                  G4int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
                  G4double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
                  G4double alphaval = 0.05;                          // Critical apha value

   // Input values
 // **************************************************************
                  
                  G4int Array_Size = ID_Number + 1;
                  G4double x[Array_Size] {};
                  G4double y[Array_Size] {};
                  G4double erry[Array_Size];   
                  
                  
                
                  x[0] = 0.;
                  y[0] = 0.;
                  
         //         G4cout << " x = " << G4BestUnit(x[0],"Length")<< G4endl;
          //        G4cout << " y = " << G4BestUnit(y[0],"Length")<< G4endl;
                  
                  for (G4int ii = 1; ii<Array_Size; ii++){
                      x[ii] = function_X_T[0][ii-1];
            //          G4cout << " x = " << G4BestUnit(x[ii],"Length")<< G4endl;
                      y[ii] = std::fabs(function_Y_T[0][ii-1]); // function x(t)
            //          G4cout << " y = " << G4BestUnit(y[ii],"Length")<< G4endl;
                     
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
           
 
   		 XTWXInv = Make2DArray(k+1,k+1);
   		 Weights = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
   		 CalculateWeights(erry, Weights, n, wtype);
    

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
 	//   DisplayStatistics(n,nstar,k,RSS,R2,R2Adj,SE);
  
 
    // Display ANOVA table
    // **************************************************************
 	 //  DisplayANOVA(nstar, k, TSS, RSS, J);
    

    // Write the prediction and confidence intervals
    // **************************************************************
  	//  WriteCIBands("CIBands2.dat",x,coefbeta,XTWXInv,tstudentval,SE,n,k);

    // Display the covariance and correlation matrix
    // **************************************************************
  	 // DisplayCovCorrMatrix(k, SE, fixedinter, XTWXInv);
           
               G4double l[ID_Number];
               G4double z[ID_Number];
               G4double erry_Z[ID_Number];
               
               for (G4int ii = 0; ii<ID_Number; ii++){
                     l[ii] = std::fabs(SetLength(k, coefbeta, x[ii+1]) - SetLength(k, coefbeta, 0.)); //I'm hz why fabs
                      
                     z[ii] = function_Z_T[0][ii]; // function x(t)
                     /// erryZ[ii] = 0.;
                     G4cout << " t = " << G4BestUnit(SetLength(k, coefbeta, x[ii+1]),"Length") << " Co = " << G4BestUnit(SetLength(k, coefbeta, 0.),"Length") << G4endl;
                     G4cout << " l = " << G4BestUnit(l[ii],"Length") << " z = " << G4BestUnit(z[ii],"Length") << G4endl;
               }
               
                G4cout << "Polynomial fit Z(l)!" << G4endl;
                
                size_t k_Z = 1;
                size_t n_Z = 0;                           // Number of data points (adjusted later)
    		size_t nstar_Z = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    		G4double coefbeta_Z[k_Z+1];                            // Coefficients of the polynomial
    		G4double serbeta_Z[k_Z+1];                             // Standard error on coefficients
    		G4double tstudentval_Z = 0.;                         // Student t value
    		G4double SE_Z = 0.;                                  // Standard error
    
   		G4double **XTWXInv_Z;                                // Matrix XTWX Inverse [k+1,k+1]
    		G4double **Weights_Z;                                // Matrix Weights [n,n]
                 
                  // Initialize values
    // **************************************************************
                  n_Z = sizeof(l)/sizeof(G4double);
    		  nstar_Z = n_Z-1;
    		  if (fixedinter) nstar_Z = n_Z;
           
 
    		  XTWXInv_Z = Make2DArray(k_Z+1,k_Z+1);
    		  Weights_Z = Make2DArray(n_Z,n_Z);

    // Build the weight matrix
    // **************************************************************
    		  CalculateWeights(erry_Z, Weights_Z, n_Z, wtype);
    

    // Calculate the coefficients of the fit
    // **************************************************************
    		  PolyFit(l,z,n_Z,k_Z,fixedinter,fixedinterval,coefbeta_Z,Weights_Z,XTWXInv_Z);


    // Calculate related values
    // **************************************************************
       		 G4double RSS_Z = CalculateRSS(l,z,coefbeta_Z,Weights_Z,std::fixed,n_Z,k_Z+1);
    		 G4double TSS_Z = CalculateTSS(l,z,coefbeta_Z,Weights_Z,fixedinter,n_Z,k_Z+1);
    		 G4double R2_Z = CalculateR2COD(l,z,coefbeta_Z,Weights_Z,fixedinter,n_Z,k_Z+1);
    		 G4double R2Adj_Z = CalculateR2Adj(l,z,coefbeta_Z,Weights_Z,fixedinter,n_Z,k_Z+1);

    		 if ((nstar_Z-k_Z)>0) {
        	   SE_Z = std::sqrt(RSS_Z/(nstar_Z-k_Z)); 
        	   tstudentval_Z = std::fabs(CalculateTValueStudent(nstar_Z-k_Z, 1.-0.5*alphaval)); 
     		 }
    		 G4cout << "t-student value: " << tstudentval_Z << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
  		 CalculateSERRBeta(fixedinter,SE_Z,k_Z,serbeta_Z,XTWXInv_Z);

    // Display polynomial
    // **************************************************************
   		DisplayPolynomial(k_Z);

    // Display polynomial coefficients
    // **************************************************************
    		DisplayCoefs(k_Z, nstar_Z, tstudentval_Z, coefbeta_Z, serbeta_Z);

                Nepopal ++;
                G4double Z0 = coefbeta_Z[0];
                G4double Z_nach = (*fHitsCollection)[step-1]->GetTrackVertex().z();
                G4double ZAb = std::fabs(Z0-Z_nach);
                if(ZAb < 5*mm) Popal++;
                G4double result = Popal/Nepopal;
                G4cout << " Z0 = " << G4BestUnit(ZAb,"Length") << " Result = " << result << G4endl;
               
                analysisManager->FillH1(3, ZAb);
               
               
               }      
             // }
               
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
          
        } // end of fit
      break;    
               }  //end of prev!=next chamber number
 /*
    L_pred = L_min;
    k_pred = k_min;
 */ 
  //  G4cout << " L_pred = " << L_pred << " k_pred = " << k_pred << G4endl;
    }  
    
   
  }
} 
}   
           
   
