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
#include "G4Cache.hh"
#include "EventAction.hh"
//#include "G4RandGauss.hh"
#include <vector>
#include <math.h>


//#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName,
                     EventAction* eventAction)
 : G4VSensitiveDetector(name), fEventAction(eventAction)
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
  G4double Energy = aStep->GetTrack()->GetKineticEnergy();
  
  
/*  G4double pz = Mom.z();
            //    G4cout << " pz = " << pz << G4endl;
  G4double Dphi = 0.;
  if(pz<0) Dphi = 3.14159-std::acos(Mom.z());
  if(pz>=0) Dphi = std::acos(Mom.z());*/
  
  
 
 // if (status != 0) return false;
  
  //G4cout << " Success " << " Pt = " << G4BestUnit(Pt,"Energy") <<" P = " << G4BestUnit(P,"Energy") << " Energy = " << G4BestUnit(Energy,"Energy") << G4endl;
  
  TrackerHit* newHit = new TrackerHit();
  G4int Track_ID = aStep->GetTrack()->GetTrackID();
  G4int PDGE = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  
//  G4cout << " PDGE = " << PDGE << "My TrackID = " << Track_ID << " Parent TrackID = " << status<< G4endl;
  
 // newHit->SetEdep(edep);
  if (PDGE != 13) return false;
  G4ThreeVector Pos = aStep->GetTrack()->GetVertexPosition();
  newHit->SetTrackID(Track_ID);
 // G4cout << " px = " << Mom.x() << " py = " << Mom.y() << " pz = " << Mom.z() << G4endl;
  G4ThreeVector Mom = aStep->GetTrack()->GetVertexMomentumDirection();
  newHit->SetMomDir(Mom);
  
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
   
 G4int step = 0, k, ID_Number = 0, i_step = 0, k_min = 0, k_pred = 0, Nb;
 G4int nextnumber, prevnumber, prevTrackID, nextTrackID, Number, status, No_parabol = 0, All_Track = 0, All_Track_no_fit = 0, All_Track_dphi = 0, All_Track_coef = 0, All_Track_versh = 0;
 G4ThreeVector pointH0, pointH1, pointH;
 G4double L, L_min, L_pred, time_for_Track_Global, E_dep, Popal=0., All = 0., L_fit_Z, L_fit_YX, func, Dtime;
 G4bool No_fit = true;
 
 size_t k_fit = 2;
 
 //G4cout << " St " << No_fit << G4endl;
   
 G4int Array[100] {};
 
 G4MapCache <G4int, G4double> Z_Map;
 G4MapCache <G4int, G4double> Z0_Map;
 G4MapCache <G4double, G4int> Z0up_Map;
 
 std::vector <G4double> Z_nach_vector;
 
        
            
 for (G4int i=0;i<nofHits;i++)
    {
        
    if (i==step){
 
//  G4cout << " step " << step << G4endl;
       nextnumber = 0;
            
       prevnumber = (*fHitsCollection)[i]->GetChamberNb();
       
       prevTrackID = (*fHitsCollection)[i]->GetTrackID();
      // G4cout << " PREV = " << prevTrackID<<G4endl;
       pointH0 = (*fHitsCollection)[i]->GetLocalPos();
       L = std::sqrt(std::pow(pointH0.y(),2)+std::pow(pointH0.x(),2)); 
       
       k_min = i;
       L_min = L;
       L_pred = L;
       k_pred = i;
     
  //   G4cout << " L_nach =" << L_min << " k_nach = " << k_min << G4endl; 
       for (k = i+1; k<nofHits; k++){
         
    
           
    //   if(status)
        //  G4cout << "Parent ID" << status << G4endl;
         // G4int IDTr = (*fHitsCollection)[step]->GetTrackID();
        //  G4cout << "MY ID" << IDTr << G4endl;     
         
          
          nextnumber = (*fHitsCollection)[k]->GetChamberNb();
          nextTrackID = (*fHitsCollection)[k]->GetTrackID();
         // G4cout << " NEXT = " << nextTrackID<<G4endl;
          pointH1 = (*fHitsCollection)[k]->GetLocalPos();
         L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
      ///  G4cout << " PREV = " << prevnumber << " NEXT = " << nextnumber << " L = " << L << " k = " << k << G4endl;
       // G4cout << " PR ID = " << prevTrackID << " NX ID = " << nextTrackID << G4endl;  
         
          if(L_min > L) { 
          
            L_min = L;
            k_min = k;
          }  
            
           
       
          if (prevnumber!=nextnumber){ 
      
           //  pointH1 = (*fHitsCollection)[k-1]->GetLocalPos();
             
            
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
     
  
      
   //   G4cout <<" POS1 X = " <<pointH1.x()<<" Y = "<< pointH1.y() << " POS2 X = " << pointH2.x()<<" Y = "<< pointH2.y()  << G4endl;
    //  L2=std::sqrt(std::pow(pointH2.y(),2)+std::pow(pointH2.x(),2)); 
 //   G4cout << " SQRT = " << G4BestUnit(std::pow(pointH2.y()-pointH1.y(),2)+std::pow(pointH2.x()-pointH1.x(),2),"Length") << G4endl;
    /* if(std::pow(pointH1.y()-pointH0.y(),2)+std::pow(pointH1.x()-pointH0.x(),2)>0){
    //фформула скрещивающихся
         L = std::fabs(((pointH1.x()-pointH0.x())*(pointH1.y())-(pointH1.y()-pointH0.y())*(pointH1.x())))/std::sqrt(std::pow(pointH1.y()-pointH0.y(),2)+std::pow(pointH1.x()-pointH0.x(),2));  
        }
     else {
         L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
           }
     */    
    // G4cout << " L_result " << L_pred << " k_result " << k_pred << G4endl;      
             
          Array[i_step] = k_pred;
       //   pointH1 = (*fHitsCollection)[k_pred]->GetPos();
      ////    G4cout <<" x = " << G4BestUnit(pointH1.x(),"Length")<< " y = " << G4BestUnit(pointH1.y(),"Length") << " z = " << G4BestUnit(pointH1.z(),"Length")<< G4endl;
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
       G4double time = (*fHitsCollection)[k_pred]->GetTime();
     //  G4double time1 = (*fHitsCollection)[step]->GetTime();
      // G4double time = (time1 + time2)/2;
       
       auto analysisManager = G4AnalysisManager::Instance();
       
       analysisManager->FillH1(0, time);
       //analysisManager->FillH1(1, status);
  
       // if(time<0) G4cout << "WARNING!!!!!" << G4endl;
       // My strange attempt to solve this problem
     
     /*  G4double N = 0.;
       G4int j = 0;
       while (N>=0)
       {
         N = time - j*76*ns;
         j++;
        }
      //j= number of interaction - 1, we start with hist 2 => j-1+2=j+1
       G4int potok = j+1;  */    
      
       
    /*   TIME=time+2.71012+1.21564*L_pred+6.82868*L_pred*L_pred;
       analysisManager->FillH1(1, L_pred);  
    //  G4cout << " potok = " << potok << G4endl;     
    */
      // TIME=time+2.71012+1.21564*L+6.82868*L*L;
       
      // TIME = G4RandGauss::shoot(time+2.71012+1.21564*L+6.82868*L*L,20*ns);
  
       TIME = time+2.71012+1.21564*L_pred+6.82868*L_pred*L_pred;
       analysisManager->FillH1(1, L_pred);
       analysisManager->FillH1(2, TIME);
          // number of interaction - 1
     // status   }
    // if L>0  }
     
      
      // формула скрещивающихся, редуцированная до моего случая R = |(x2-x1)y2-(y2-y1)x2|/sqrt((y2-y1)^2+(x2-x1)^2) x1,y1 - H1, x2, y2 - H2, x,y - C
      //формула для 1 точки - Пифагор, ибо z один => делаем прямая лежит в плоскость z=z2, тогда расстояние находим между точками x2,y2,z2, x,y,z, где z2=z,  а знаяит, уходит
      
        step = k;
                     
             
        // block with fit     
             if (prevTrackID!=nextTrackID && ID_Number > k_fit){
               
                All_Track++;
                 
                  G4cout << " All_Track = " << All_Track << G4endl;
               G4int Array_Size = ID_Number+1;
               
                G4double x[Array_Size] {};
                G4double y[Array_Size] {};
                G4double erry[Array_Size];  
                   
                x[0] = 0.;
                y[0] = 0.;
                
             /*   G4cout <<" x  " <<" y  " << G4endl;
                
                 for (G4int ii = 1; ii<Array_Size; ii++){
                         pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                         x[ii] = pointH.x(); 
                         y[ii] = pointH.y(); // function x(t)
                        G4cout << x[ii]<< " " << y[ii]<< G4endl;
                   //  G4cout <<" x = " << x[ii]<< " y = " << y[ii]<< G4endl;
                         
                     }*/
                                
                
                
              //  if((*fHitsCollection)[Array[1]]->GetPos().x()*(*fHitsCollection)[Array[1]]->GetPos().y()>0){  
            //    G4cout <<" x = " << G4BestUnit(x[0],"Length")<< " y = " << G4BestUnit(y[0],"Length")<< G4endl;
                  
               
                    for (G4int ii = 1; ii<Array_Size; ii++){
                         pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                         x[ii] = pointH.x();
                         y[ii] = std::fabs(pointH.y()); // function x(t)
                        if(ii != 0){              
                       /* if(std::fabs(x[ii]) < std::fabs(x[ii-1]) || std::fabs(y[ii]) < std::fabs(y[ii-1])){ 
                          No_fit = false; 
                   //    G4cout << " WARNING " << std::fabs(x[ii]) - std::fabs(x[ii-1]) << G4endl; 
                          No_parabol++;
                          break;
                        }*/
                        }   
                       
                    //    G4cout << x[ii]<< " " << y[ii]<< G4endl;
                     G4cout <<" x = " << x[ii]<< " y = " << y[ii]<< G4endl;
                        erry[ii] = 0.;  
                     }
               //   }
               
                         G4cout << "Polynomial fit Y(x)!" << G4endl;
    
                  G4bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
                  G4int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
                  G4double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
                  G4double alphaval = 0.05;                          // Critical apha value

   // Input values
 // **************************************************************
                            
             
                  
                // G4cout << "Status " << No_fit << G4endl;
                  
             //    if(No_fit){
                All_Track_no_fit++;
                G4cout << " All_Track_no_fit = " << All_Track_no_fit << G4endl;
                 // Definition of other variables
    // **************************************************************
    		 size_t n = 0;                                    // Number of data points (adjusted later)
    		 size_t nstar = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    	  	 G4double coefbeta[k_fit+1];                            // Coefficients of the polynomial
   		 G4double serbeta[k_fit+1];                             // Standard error on coefficients
   		 G4double tstudentval = 0.;                         // Student t value
   		 G4double SE = 0.;                                  // Standard error
    
    		 G4double **XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
   		 G4double **Weights;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
    		 n = sizeof(x)/sizeof(G4double);
    		 nstar = n-1;
    		 if (fixedinter) nstar = n;
           
 
   		 XTWXInv = Make2DArray(k_fit+1,k_fit+1);
   		 Weights = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
   		 CalculateWeights(erry, Weights, n, wtype);
    

    // Calculate the coefficients of the fit
    // **************************************************************
   		 PolyFit(x,y,n,k_fit,fixedinter,fixedinterval,coefbeta,Weights,XTWXInv);


    // Calculate related values
    // **************************************************************
   		 G4double RSS = CalculateRSS(x,y,coefbeta,Weights,std::fixed,n,k_fit+1);
   		 G4double TSS = CalculateTSS(x,y,coefbeta,Weights,fixedinter,n,k_fit+1);
   		 G4double R2 = CalculateR2COD(x,y,coefbeta,Weights,fixedinter,n,k_fit+1);
   		 G4double R2Adj = CalculateR2Adj(x,y,coefbeta,Weights,fixedinter,n,k_fit+1);

    if ((nstar-k_fit)>0) {
        SE = std::sqrt(RSS/(nstar-k_fit)); 
        tstudentval = std::fabs(CalculateTValueStudent(nstar-k_fit, 1.-0.5*alphaval)); 
    }
 //   G4cout << "t-student value: " << tstudentval << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE,k_fit,serbeta,XTWXInv);

    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k_fit);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k_fit, nstar, tstudentval, coefbeta, serbeta);
                 
    
          //  G4double Versh = - coefbeta[1]/(2*coefbeta[2]);
         //   G4cout << " Versh = " <<Versh << G4endl;
            
          // G4double KOEF[3];
         //  KOEF[0] = coefbeta[0]; KOEF[1] = coefbeta[1]; KOEF[2] = coefbeta[2];   
           
                
           if(std::fabs(coefbeta[0])<1){
                     All_Track_coef++;
                     G4cout << " All_Track_coef = " << All_Track_coef << G4endl; 
                    //Test for Versh
              
                     
               G4double l[ID_Number];
               G4double z[ID_Number];
               G4double erry_Z[ID_Number];
             //  All_Track_versh++;
           //     G4cout << " All_Track_versh = " << All_Track_versh << " z = "<< G4endl;
              //// G4cout << " KOEF0" << coefbeta[0] <<  " KOEF1" << coefbeta[1] <<  " KOEF2" << coefbeta[2] << G4endl;
          //     G4cout << " coefbeta " << KOEF[0] << " " << KOEF[1] << " " << KOEF[2] << G4endl;  
               for (G4int ii = 0; ii<ID_Number; ii++){
                     l[ii] = std::fabs(SetLength(k_fit, coefbeta, x[ii+1])); 
                    
                  //   G4cout << " x = " << x[ii+1] << " l = " << l[ii]  << G4endl;
                     pointH = (*fHitsCollection)[Array[ii]]->GetPos();
                  //   G4cout << " k from Z " << Array[ii]<<G4endl;
                     z[ii] = pointH.z(); // function x(t)
                 //    G4cout << z[ii]<< " " << l[ii]<< G4endl;
                     erry_Z[ii] = 0.;
                   //   G4cout << G4BestUnit(x[ii+1],"Length")<< " " << G4BestUnit(y[ii+1],"Length") << " " << G4BestUnit(z[ii],"Length")<< G4endl;
                   
                  
                   
            //       G4cout << x[ii+1]<< " " << y[ii+1] << G4endl;
                   //  G4cout << " t = " << G4BestUnit(SetLength(k, coefbeta, x[ii+1]),"Length") << " Co = " << G4BestUnit(SetLength(k, coefbeta, 0.),"Length") << G4endl;
         //            G4cout << " l = " << G4BestUnit(l[ii],"Length") << " z = " << G4BestUnit(z[ii],"Length") << G4endl;
                /*     Nb = (*fHitsCollection)[Array[ii]]->GetChamberNb()/1000;
                     func = std::fabs(y[ii+1]-(coefbeta[2]*x[ii+1]*x[ii+1]+coefbeta[1]*x[ii+1]+coefbeta[0]));
                     G4cout << " Nb = " << Nb << " funcYX = " << G4BestUnit(func,"Length") << G4endl;
                     
                     analysisManager->FillH2(0, Nb, func);*/
                
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
    	//	 G4cout << "t-student value: " << tstudentval_Z << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
  		 CalculateSERRBeta(fixedinter,SE_Z,k_Z,serbeta_Z,XTWXInv_Z);

    // Display polynomial
    // **************************************************************
   	//	DisplayPolynomial(k_Z);

    // Display polynomial coefficients
    // **************************************************************
    //		DisplayCoefs(k_Z, nstar_Z, tstudentval_Z, coefbeta_Z, serbeta_Z);

      
                G4double Z0 = coefbeta_Z[0];
                G4double Z_nach = (*fHitsCollection)[k_pred]->GetTrackVertex().z();
             //   G4cout << " k from Z_nach and ID " << step-1 << " Z0 = " << Z0 << G4endl;
                G4int TrackID = (*fHitsCollection)[k_pred]->GetTrackID();
              //  analysisManager->FillH2(1, prevNumber, L_fit_Z);
                G4double Zraz = Z0-Z_nach;
                
               /* G4double ZAb = std::fabs(Z0-Z_nach);
                if(ZAb < 10*mm) Popal++;
                G4double result = Popal/All;*/
                
                analysisManager->FillH1(3, Z_nach);
                analysisManager->FillH1(4, Z0);
               
               
                   
	/////	}
		
                
               
               
       //        G4cout << " KOEF0_Z " << coefbeta_Z[0] <<  " KOEF1_Z " << coefbeta_Z[1] << G4endl;
   //      G4double LD, DTIME;
  //       G4ThreeVector pointHD;
               
  //             G4cout << " x0 " << " y0 " << G4endl;
               
                for (G4int ii = 0; ii<ID_Number; ii++){ 
         //       Nb = (*fHitsCollection)[Array[ii]]->GetChamberNb()/1000;
          //      pointHD = (*fHitsCollection)[Array[ii]]->GetLocalPos();
         //       LD = std::sqrt(std::pow(pointHD.y(),2)+std::pow(pointHD.x(),2)); 
         //       Dtime = (*fHitsCollection)[Array[ii]]->GetTime();
          //      DTIME = Dtime+2.71012+1.21564*LD+6.82868*LD*LD;         
     //           G4cout << " Time =  " << Dtime << " Zraz = " << Zraz << G4endl;
             //   func = std::fabs(z[ii]-(coefbeta_Z[1]*l[ii]+coefbeta_Z[0]));
              //  G4cout<< x[ii+1]<<" " <<  y[ii+1] << /*" l = " << l[ii]<< " z0 = " << z[ii]<< */G4endl;
                //     G4cout << " Nb = " << Nb << " funcZ = " << G4BestUnit(func,"Length") << G4endl;
                     
                //     analysisManager->FillH2(1, Nb, func);
               // analysisManager->FillH2(0, Z0, DTIME);
                }
               // G4double Dphi = std::atan2(z[ID_Number-1],y[ID_Number]);
                
       //         G4cout << " y = " << y[1] << " z = " << z[0] << " y/z = "<< y[1]/z[0] << G4endl;
                G4ThreeVector Vertex_Mom = (*fHitsCollection)[k_pred]->GetMomDir();
                G4double pz = Vertex_Mom.z();
            //    G4cout << " pz = " << pz << G4endl;
                G4double Dphi2 = 0.;
                if(pz<0) Dphi2 = 3.14159-std::acos(Vertex_Mom.z());
                if(pz>=0) Dphi2 = std::acos(Vertex_Mom.z());
                G4cout << " Dphi2 = " << Dphi2*180/3.14159 << G4endl;
                
           //     G4cout << " Z0 = " << G4BestUnit(Z0,"Length") << " Z_nach = " << G4BestUnit(Z_nach,"Length") <<  " ZRAZ = " << G4BestUnit(Zraz,"Length") << G4endl;
                //Z0 - reconstructed, Z_nach - truly Z
                analysisManager->FillH2(1, Dphi2, Zraz);
                
                if(Dphi2>=0.5){
                All_Track_dphi++;
                 // if(Zraz>2*cm){ 
                 if(x[1]*Versh>0 && std::fabs(Versh)<std::fabs(x[Array_Size-1])){
                     
                #ifdef G4MULTITHREADED
		static G4Mutex stuffMutex = G4MUTEX_INITIALIZER;	
		G4AutoLock al(&stuffMutex);
		#endif
		static std::ofstream stuff_xy("Wrong_xy.csv");
		static std::ofstream stuff_zl("Wrong_zl.csv");
		static bool first = true;
		if (first) {
		first = false;
		stuff_xy << "#,x/mm,y/mm,TrackID" << std::endl;
		stuff_zl << "#,z/mm,l/mm,TrackID" << std::endl;
			}
			stuff_xy << "," << 0 << "," << 0 << ","<< prevTrackID << std::endl;
			for (G4int ii = 1; ii<Array_Size; ii++){
                        stuff_xy << "," << x[ii] << "," << y[ii] <<"," << prevTrackID << std::endl;
                        }   
                        stuff_xy << "," << coefbeta[0] << "," << coefbeta[1]  << "," << coefbeta[2] << std::endl;
                        
             		for (G4int ii = 1; ii<Array_Size; ii++){
                        stuff_zl << "," << z[ii-1] << "," << l[ii-1] <<"," << prevTrackID << std::endl;
                        }   
                        stuff_zl << "," << coefbeta_Z[0] << "," << coefbeta_Z[1]  << "," << Z_nach << std::endl;        
                        }
               //  }  
                G4cout << " All_Track_dphi = " << All_Track_dphi << G4endl;
        //        G4cout <<  " TrackID = " << TrackID <<" Z0 = " << G4BestUnit(Z0,"Length") << " Z_nach = " << G4BestUnit(Z_nach,"Length") <<   " ZRAZ = " << G4BestUnit(Zraz,"Length")<< G4endl;
                if(std::fabs(Zraz)>3*cm) G4cout << " WARNING!!!! " << G4endl;
                
                
                Z0_Map[TrackID] = Z_nach; // ID, Znach
                Z_Map[TrackID] = Z0;  // ID, Z0
                Z0up_Map[Z0] = TrackID; // Z0, ID
                
                }
            /*    G4double Z_0, Z_01;
             //   G4cout << " z_N = " << z[ID_Number-1] <<" l_N = " << l[ID_Number-1]  << " z_0 = " << z[0]<< " l_0 = " << l[0] << G4endl; 
                Z_0 = z[ID_Number-1] - (z[ID_Number-1]-z[0])*l[ID_Number-1]/(l[ID_Number-1]-l[0]);
                Z_01 = (z[ID_Number-1]-0.3*mm) - ((z[ID_Number-1]-0.3*mm)-(z[0]+0.3*mm))*l[ID_Number-1]/(l[ID_Number-1]-l[0]);
                
           //     G4cout << " Z_0 = " << Z_0 << " Z_01 = " << Z_01 <<" Delta Z = " << Z_0-Z_01 << G4endl;
                if(std::fabs(Z_0-Z_nach)>1*cm){    
             //   G4cout << " WARNING " << G4endl;
                    for (G4int ii = 0; ii<ID_Number; ii++){
                   //     G4cout <<" x = " << G4BestUnit(x[ii+1],"Length")<< " y = " << G4BestUnit(y[ii+1],"Length") << " z = " << G4BestUnit(z[ii],"Length")<< G4endl;
                    }
                } 
                analysisManager->FillH1(6,Z_01-Z_0);*/
              //}
              
             
              
      //        }
            //}                  
               }//end of fit
             //  No_fit = true;
               ID_Number = 0;  
               i_step = 0; 
          
        } //end of prev!=next TrackID
     
      break;    
               }  //end of prev!=next chamber number
 
    L_pred = L_min;
    k_pred = k_min;
  
   // G4cout << " L_pred = " << L_pred << " k_pred = " << k_pred << G4endl;
    }  
    
   
  }
 } 


// G4cout << " B " << G4endl; 
 G4double H = 4.0*cm, k0, k1, Z_nach_step = 0, Z0_step, Z_end; //set our step   
 G4int Num_success = 0, Num_all = 0, Num_wrong = 0, step_N = 1;
 G4bool flag = true, test = true;
              
 G4MapCache <G4int, G4double >::iterator it = Z_Map.Begin(), it_ID, itp = Z_Map.Begin();
 G4MapCache <G4int, G4double >::iterator it0 = Z0_Map.Begin(), it0_ID, it1 = Z0_Map.Begin(), it0p = Z0_Map.Begin();
 G4MapCache <G4double, G4int >::iterator itZ0 = Z0up_Map.Begin(), it_step = Z0up_Map.Begin(), it0upp = Z0up_Map.Begin();
 
 std::vector <G4double> Z_vector;
 std::vector <G4double> Z0_vector;
 std::vector <G4int> IDZ0;
 std::vector <G4int> IDZnach;
 
 G4double Z_step, Z_last = 0.;
 auto analysisManager = G4AnalysisManager::Instance();                                        
 
 for(G4int i = 0; itp != Z_Map.End(); itp++, i++){
 
     G4cout << i << ") Key - TrackID " << itp->first <<" Z0 = "<< G4BestUnit(itp->second,"Length") << G4endl;
 }   
 
  for(G4int i = 0; it0p != Z0_Map.End(); it0p++, i++){
 
     G4cout << i << ") Key - TrackID " << it0p->first <<" Znach = "<< G4BestUnit(it0p->second,"Length") << G4endl;
 }  
 
  for(G4int i = 0; it0upp != Z0up_Map.End(); it0upp++, i++){
 
     G4cout << i << ") Key - Z0 " << G4BestUnit(it0upp->first,"Length") <<" TrackID = "<< it0upp->second << G4endl;
 } 
 

 
 

Z0_step = itZ0->first;
G4cout << " Z0_step = " << Z0_step << G4endl;
 while (it_step != Z0up_Map.End() ){
 //for(G4int j = 0; j<5; j++){
 test = true;
 it = Z_Map.Begin(); it0 = Z0_Map.Begin();
 //G4cout << " it1 -> " << it1->second << G4endl;
 Z_step = it_step->first-0.3*cm;
// if(Z_step+0.5*cm == Z_last) G4cout << " WARNING " << G4endl;
//if(it!=Z_Map.Begin()){ it--; Z_last = it->second; it++;}
 G4cout << " Z_step = " << G4BestUnit(Z_step,"Length")  << G4endl;
 G4cout << H + Z_step << " flag = " << flag << G4endl;
 for(G4int i = 0; it != Z_Map.End(); it++, i++, it0++, itZ0++){
 
     G4cout << i << ") Key - TrackID " << it->first <<" Z0 = "<< G4BestUnit(it->second,"Length") << " Z_nach =  " << G4BestUnit(it0->second,"Length") << G4endl;
     
     if(it->second <= H + Z_step && it->second > Z_step){
     G4cout << " Z0 in H " << G4endl;
        
     //   if(it->second == Z_last) G4cout << " WARNING " << G4endl;
        
        Z_vector.push_back(it->second);  //fill vector of needed Z0
        IDZ0.push_back(it->first);     //fill vector of Z0's ID
        if(it->second >= Z0_step)  Z0_step = it->second;
        G4cout << " Z0 = " << it->second << " ID_Z0 = " << it->first  << " Z0_step = " << Z0_step << G4endl;
        
        
     } 
     
     if(it0->second <= H + Z_step && it0->second > Z_step){
        IDZnach.push_back(it0->first); 
        G4cout << " Znach in H = " << it0->second << " ID_Znach = " << it0->first << G4endl;
        Z_nach_vector.push_back(it0->second);
       
     }
 }
 
it_step = Z0up_Map.Find(Z0_step);
it_step++;
G4cout << " Z_step finally " << it_step->first << G4endl;
 //work with vector
 G4double Sum = 0., Norm = Z_vector.size(), Z_norm, Z_result; //Z_norm = averange Z, Z_result = |Z_nach(truly Z0)-Z_norm|
 G4int ID;
 
 for (G4int i = 0, j = Z_vector.size(); i < j; i++){
    Sum+=Z_vector[i];
    
 
    
   //  G4cout << " Z_nach in H = " << it0_ID->second << " ID_nach = " << it0_ID->first<<  G4endl;            //fill vector of all Z_nach's ID in H
    // Z0_inH_vector.push_back(it0_ID->second);      //fill vector of Z_nach in H    
    //} 
 }
 Z_vector.clear();
 Z_norm = Sum/Norm; 
// G4cout << " Norm = " << Norm << " Sred znach " << Z_norm << G4endl;
 

 Z0_vector.clear();
 
 Z_nach_step = 0;
// G4cout << " Size = " << Z_nach_vector.size() << G4endl;
 

 for(G4int i  = 0, j = Z_nach_vector.size(); i<j; i++){
    
  //  G4cout << " Z_nach = " << Z_nach_vector[i] << G4endl; 
    
     if(Z_nach_vector[i] > Z_step && Z_nach_vector[i] < Z_step+H){
       G4cout << " Z_nach in H0 = " << Z_nach_vector[i] << " Z_step = " << Z_nach_step << G4endl;  
     if(Z_nach_step == 0){ Z_nach_step = Z_nach_vector[i]; 
        G4cout << " Z_step = " << Z_nach_step << G4endl;
     }
     else{
     if(Z_nach_step != Z_nach_vector[i]) {Num_wrong++; G4cout << " WRONG " << " Wrong_Num = "<<Num_wrong << G4endl; flag = false; Z_nach_step = Z_nach_vector[i];}
     }
     
     }
 
 
 } 

 
 
 if(flag){
    
    G4cout << " Array_Test " << G4endl;
    
    for (G4int i = 0, j = IDZnach.size(); i<j; i++){
    G4cout << " IDZnach = " << IDZnach[i] << G4endl;
    }
     for (G4int i = 0, j = IDZ0.size(); i<j; i++){
    G4cout << " IDZ0 = " << IDZ0[i] << G4endl;
    }
    
  
    
 }
 
 if(IDZnach!=IDZ0 || !flag) {Num_wrong++; G4cout << " Wrong_Num = " << Num_wrong << G4endl;}; 
 G4cout << " ID size = " << IDZ0.size() <<" ID_nach size = " << IDZnach.size() << " flag = " << flag <<" Wrongs = " << Num_wrong<< G4endl;
 IDZnach.clear();
 IDZ0.clear();
 
 
 flag = true;
 
 
 Z_nach_vector.clear();
 }// end of while
 
 
 
 it1 = Z0_Map.Begin();
 it0 = Z0_Map.Begin();

 for(G4int i = 1; it0 != Z0_Map.End(); it0++, i++){
     if(i==1) {Num_all++; G4cout << " All_Vertex = " << Num_all << G4endl;}
     G4cout << i << ") Key - TrackID " << it0->first << " Z_nach " << G4BestUnit(it0->second,"Length") << " All_Num "<< G4endl;
     if(it0!=Z0_Map.Begin()){
        if(it1->second!=it0->second) {Num_all++; G4cout << " All_Vertex = " << Num_all << G4endl;}
     it1++;
     }
 } 
 G4cout << " Vsego = " << Num_all << " " << Num_wrong << G4endl;
 G4int Cuts = All_Track -All_Track_dphi; 
 G4cout << " Cuts = " << Cuts << G4endl;
 
 #ifdef G4MULTITHREADED
static G4Mutex stuffMutex = G4MUTEX_INITIALIZER;
G4AutoLock al(&stuffMutex);
#endif
static std::ofstream stuff_v("Vertex.csv");
static std::ofstream stuff_inf("Inf_Cut.csv");
static bool first = true;
if (first) {
first = false;
stuff_v << "#,wrong,Vsego" << std::endl;
stuff_inf << "#,cuts,All_Track" << std::endl;
}
//stuff_v << "," << Num_wrong << "," << Num_all << std::endl;
stuff_inf << "," << Cuts << "," << All_Track << std::endl;

// auto analysisManager = G4AnalysisManager::Instance();
 //analysisManager->CreateH1("P","pop/all", 240, 0*cm, 60.*cm); 
}   
           
   
