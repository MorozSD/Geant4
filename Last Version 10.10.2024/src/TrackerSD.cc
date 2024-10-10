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
#include <vector>
#include <math.h>
#include "kmeans.hh"
//#include "IDManager.h"




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName,
                     EventAction* eventAction)
 : G4VSensitiveDetector(name), fEventAction(eventAction)
{
  collectionName.insert(hitsCollectionName);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD()
{}*/

SensitiveDetector::SensitiveDetector(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDetector::~SensitiveDetector()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
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

G4bool SensitiveDetector::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory* )
{


  G4double edep = aStep->GetTotalEnergyDeposit();
  G4int status = aStep->GetTrack()->GetParentID();
  G4double Energy = aStep->GetTrack()->GetKineticEnergy();
  G4int PDGE = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4double Etot = aStep->GetPreStepPoint()->GetTotalEnergy();
  G4String Name_Maybe = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
  G4int ID_Maybe = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
  
  //G4cout << " Name " << Name_Maybe << " Vol " << aStep->GetPreStepPoint()->GetTouchable()->GetVolume() << " ID " << ID_Maybe << G4endl;
  if (status != 0) return false;
  if(PDGE == 11 || PDGE == -11) return false;
  if(Etot < 100*MeV) return false;
  
 
  TrackerHit* newHit = new TrackerHit();
  
  G4ThreeVector Pos = aStep->GetTrack()->GetVertexPosition();
  newHit->SetTrackVertex(Pos);
  
  G4int Track_ID = aStep->GetTrack()->GetTrackID();
  newHit->SetTrackID(Track_ID);
 
  G4ThreeVector Mom = aStep->GetTrack()->GetVertexMomentumDirection();
  newHit->SetMomDir(Mom);
  
  //G4String Name_Maybe = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
 // G4int ID_Maybe = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
  
  //IDManager MyIDManager;
  
  newHit->SetChamberNb(ID_Maybe);
  //newHit->SetChamberLayer(MyST_ID.get_LayerNum(ID_Maybe));
  //MyIDManager.printVolume_info(ID_Maybe);
  //MyIDManager.printST_info(ID_Maybe);
 
  
  newHit->SetChamberLayer(ID_Maybe);
//  G4cout << " Name " << Name_Maybe << " ID " << ID_Maybe << G4endl;
  //G4cout << " Layer " <<  MyST_ID.get_LayerNum(ID_Maybe) << " ID " << ID_Maybe << G4endl;
 
  G4ThreeVector worldPosition = aStep->GetPreStepPoint()->GetPosition();

  G4ThreeVector pointGlob = aStep->GetPostStepPoint()->GetPosition();
  newHit->SetPos (pointGlob);
  
  G4ThreeVector point0 = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  newHit->SetLocalPos (point0);
  
  G4double time0 = aStep->GetPreStepPoint()->GetGlobalTime();
  newHit->SetTime (time0);
  
  G4int prim_charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();  
  newHit->SetCharge(prim_charge);
  
  newHit->SetTrStat(aStep->GetTrack()->GetParentID());
     
  fHitsCollection->insert( newHit );
  
 // newHit->Print();


  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......'

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{

 /*G4int nofHit = 1;
 G4cout
   << "-------->Hits Collection: in this event they are " << nofHit
   << " hits in the tracker chambers: " <<  G4endl;*/
 G4int nofHits = fHitsCollection->entries();
   
 G4int step = 0, k, ID_Number = 0, i_step = 0, k_min = 0, k_pred = 0, Nb;
 G4int nextnumber, prevnumber, prevTrackID, nextTrackID, Ignore_Track = 0, Number, status, No_parabol = 0, All_Track = 0, All_Track_no_fit = 0, All_Track_dphi = 0, All_Track_z = 0;
 G4ThreeVector pointH0, pointH1, pointH, pointGlob;
 G4double L, L_min, L_pred, time_for_Track_Global, E_dep, Popal=0., All = 0., L_fit_Z, L_fit_YX, func, Dtime;
 G4bool No_fit = true;
 
 G4double RSS_z_max = 1.0;
 
 G4int Layer_Number, BAD_Tracks = 0, ALL_FIT_Tracks = 0, ALL_Tracks = 0, RSS_z_Tracks = 0, RSS_Tracks = 0, RSS_Tracks_Bad = 0, RSS_z_Tracks_Bad = 0, Dang = 0;
 
 
 size_t k_fit = 2;
   
 G4int Array[100] {};
 
 G4MapCache <G4int, G4double> Z0_Map;
 G4MapCache <G4int, G4double> Z_Map;
 G4MapCache <G4double, G4int> Z0up_Map;
 
 std::vector <G4double> Z_nach_vector;

 
 for (G4int i=0;i<nofHits;i++)
    {
    if (i==step && (*fHitsCollection)[i]->GetTrackID() == Ignore_Track) step++;
    
    if (i==step){
 
       nextnumber = 0;    
       prevnumber = (*fHitsCollection)[i]->GetChamberNb();
       
       prevTrackID = (*fHitsCollection)[i]->GetTrackID();
       pointH0 = (*fHitsCollection)[i]->GetLocalPos();
       L = std::sqrt(std::pow(pointH0.y(),2)+std::pow(pointH0.x(),2)); 
       
       k_min = i;
       L_min = L;
       L_pred = L;
       k_pred = i;
      
       G4int itnumber = 0;
       
       for (k = i+1; k<nofHits; k++){
         
          
       nextnumber = (*fHitsCollection)[k]->GetChamberNb();
       nextTrackID = (*fHitsCollection)[k]->GetTrackID();
        
       pointH1 = (*fHitsCollection)[k]->GetLocalPos();
       L = std::sqrt(std::pow(pointH1.y(),2)+std::pow(pointH1.x(),2)); 
       pointGlob = (*fHitsCollection)[k]->GetPos();
       itnumber++;
       
      // G4cout << " PREV " << prevnumber << " NEXT " << nextnumber << " L " << L << " L_min " << L_min << G4endl;
       
       // We found the nearest point to the wire, not crossed lines
       
       if(L_min > L) { 
          L_min = L;
          k_min = k;
       } 
        
        //IF ONLY ONE TRACK ALGORITHM WON'T FIT THE TRACK
        if (prevnumber!=nextnumber){        
            
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
        Layer_Number = (*fHitsCollection)[k-1]->GetChamberLayer();
        
        Array[i_step] = k_pred;
      
        i_step++;
             
        ID_Number++;
               
        G4double TIME = 0.;
   
        G4double time = (*fHitsCollection)[k_pred]->GetTime();
     
        auto analysisManager = G4AnalysisManager::Instance();
       
        analysisManager->FillH1(0, time);
  
        if(time<0) G4cout << "WARNING!!!!!" << G4endl;
       
        TIME = time+2.71012+1.21564*L_pred+6.82868*L_pred*L_pred;
        analysisManager->FillH1(1, L_pred);
        //analysisManager->FillH1(2, TIME); 
     
        step = k;
        
        // My strange attempt to solve this problem
      G4double N = 0.;
      G4int j = 0;
      while (N>=0)
      {
        N = time - j*76*ns;
        j++;
      }
      G4int potok = j-1;
           
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
                     
       // G4cout << " TIME " << TIME << " time " << time << G4endl;       
        // block with fit  
             
   
        if (prevTrackID!=nextTrackID){
        
       
        	Ignore_Track = prevTrackID;
        	ALL_Tracks++;              
        	if(ID_Number <= k_fit){ 
        		ID_Number = 0;  
            		i_step = 0; 
            		break;
        	}
       
         	G4int Array_Size = ID_Number+1;
         	      
         	G4double x[Array_Size] {}, x_0[Array_Size] {};
         	G4double y[Array_Size] {}, y_0[Array_Size] {};
         	G4double erry[Array_Size];  
         	G4int chamb, chamb_prev = 0;
         	G4double x_prev, y_prev, X, Y;
         	G4bool Y_of_x = true, X_of_y = true, Another = true;
                
         	#ifdef G4MULTITHREADED
	 	static G4Mutex stuffMutex = G4MUTEX_INITIALIZER;
	 	G4AutoLock al(&stuffMutex);
	 	#endif
	 	static std::ofstream stuff_XY("XY_all.csv");
	 	static std::ofstream stuff_Zl("Zl_all.csv");
	 	static std::ofstream stuff_ZT("Z0.txt");
	 
	 	static bool first = true;
	 	if (first) {
	 		first = false;
			stuff_XY << "#,x,y,Track_ID,Zraz" << std::endl;
			stuff_Zl << "#,l,z,Track_ID,Zraz" << std::endl;
		}

		
                   
                x[0] = 0.;
                y[0] = 0.;
                x_prev = 0.;
                y_prev = 0.;
                chamb_prev = 0;	
                
                //Definition of polar angle
                
                G4ThreeVector Vertex_Mom = (*fHitsCollection)[k_pred]->GetMomDir();
                G4double pz = Vertex_Mom.z();
                G4double Dphi2 = 0.;
                if(pz<0) Dphi2 = 3.14159-std::acos(Vertex_Mom.z());
                if(pz>=0) Dphi2 = std::acos(Vertex_Mom.z());
                
                
               /* if(Dphi2>=0.5){ 
                	All_Track_dphi++;
                            
                // IF PARTICLE IS NEGATIVE
                	G4int charge = (*fHitsCollection)[Array[1]]->GetCharge();
                              
                	pointH = (*fHitsCollection)[Array[1]]->GetPos();
                	X = pointH.x();
                	Y = pointH.y();
               
               		if(charge <= 0){
               			for (G4int ii = 1; ii<Array_Size; ii++){
            
            				chamb = (*fHitsCollection)[Array[ii-1]]->GetChamberLayer();
                			pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                			X = std::fabs(pointH.x());
                			Y = std::fabs(pointH.y());
                	
                        	if(chamb >= chamb_prev && X > x_prev){
                        		if(chamb >= chamb_prev && Y > y_prev){
                        			chamb_prev = chamb;
                        			x_prev = X;
                        			y_prev = Y;
                        		}
                        		else{
                        			Y_of_x = false;
                        			break;
                        		}
                        	
                        	}
                        	else{	
                        		X_of_y = false;
                        		break;     
                        	}               
              
                		}
                	}
                	else{
                		for (G4int ii = 1; ii<Array_Size; ii++){
            
            				chamb = (*fHitsCollection)[Array[ii-1]]->GetChamberLayer();
                			pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                			X = std::fabs(pointH.x());
                			Y = std::fabs(pointH.y());
                		
                        		if(chamb >= chamb_prev && X > x_prev){
                        			if(chamb >= chamb_prev && Y > y_prev){
                        				chamb_prev = chamb;
                        				x_prev = X;
                        				y_prev = Y;
                        		        }
                        			else{
                        				X_of_y = false;
                        		                break;
                        			}
                        	
                        		}
                        		else{	
                        			Y_of_x = false;
                        	                break;     
                        		}               
                                                  
                		}
                
                	}
                
                
                
                	if(!Y_of_x){
               			G4cout << " Y(X) " << G4endl;
                		for (G4int ii = 1; ii<Array_Size; ii++){
           
                        		pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                        		chamb = (*fHitsCollection)[Array[ii-1]]->GetChamberLayer();
                			x[ii] = pointH.x();
                        		y[ii] = pointH.y(); // function x(t)  
                        		erry[ii] = 0.;  
                        		G4cout << " x " << x[ii] << " y " << y[ii] << G4endl;                  
                        	}	
                	}	
                	else if(!X_of_y){
                		G4cout << " X(Y) " << G4endl;
                		for (G4int ii = 1; ii<Array_Size; ii++){
           
                        		pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                        		chamb = (*fHitsCollection)[Array[ii-1]]->GetChamberLayer();
                			x[ii] = pointH.y();
                         		y[ii] = pointH.x(); // function x(t) 
                         		erry[ii] = 0.;   
                         		G4cout << " x " << x[ii] << " y " << y[ii] << G4endl;                    
                        	}
                
                	}
                	else{
                		Another = false;
                		G4cout << " Another " << G4endl;
                		for (G4int ii = 1; ii<Array_Size; ii++){
            		 
                         		pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                         		x[ii] = pointH.x();
                         		y[ii] = pointH.y(); // function x(t)
                        		erry[ii] = 0.; 
                        		G4cout << " x " << x[ii] << " y " << y[ii] << G4endl;   
                    		}
                    	}           
                   // G4cout << "Polynomial fit Y(x)!" << G4endl;
              
                    	G4double coefbeta[k_fit+1], coefbeta_inv[k_fit+1], term = 0.;                            // Coefficients of the polynomial
                 
                    	G4double R2Adj_x = 0.,R2Adj_y = 0., R2Adj_z = 0., RSS = 0., RSS_z = 0.; 
                    	FitData(x,y,coefbeta,k_fit,erry,Array_Size,prevTrackID,R2Adj_x,RSS);  
                               
         	    	if( !Another){ 
         	    		for (G4int ii = 1; ii<Array_Size; ii++){
                 			term = y[ii];
                 			y[ii] = x[ii];
                 			x[ii] = term;
                 		}         	    
         	     	FitData(x,y,coefbeta_inv,k_fit,erry,Array_Size,prevTrackID,R2Adj_y,RSS);  
         	    
         	   
         	     		if(R2Adj_x<R2Adj_y){ 
         	    		
         	    			coefbeta[0] = coefbeta_inv[0]; 
         	    			coefbeta[1] = coefbeta_inv[1]; 
         	    			coefbeta[2] = coefbeta_inv[2];      	    		
         	      		}     
         	      		else{
         	      			for (G4int ii = 1; ii<Array_Size; ii++){
                 				term = y[ii];
                 				y[ii] = x[ii];
                 				x[ii] = term;
                 			}	
         	      		}
         	     	} 	
        
               		G4double l[ID_Number];
               		G4double z[ID_Number];
               		G4double erry_Z[ID_Number];
             
               		for (G4int ii = 0; ii<ID_Number; ii++){
               
                     		l[ii] = std::fabs(SetLength(k_fit, coefbeta, x[ii+1]));             
                     		pointH = (*fHitsCollection)[Array[ii]]->GetPos();
                     		z[ii] = pointH.z(); // function x(t)
                     		erry_Z[ii] = 0.;
            
               		}
               
               // G4cout << "Polynomial fit Z(l)!" << G4endl;
                
                	size_t k_Z = 1;
                	G4double coefbeta_Z[k_Z+1];                            // Coefficients of the polynomial
                
                	FitData(l,z,coefbeta_Z,k_Z,erry_Z,Array_Size-1,prevTrackID,R2Adj_z,RSS_z); 
                
                	G4double Z0 = coefbeta_Z[0];
                        G4double Z_nach = (*fHitsCollection)[k_pred]->GetTrackVertex().z();
                	G4int TrackID = (*fHitsCollection)[k_pred]->GetTrackID();
                        G4double Zraz = Z0-Z_nach;
                G4cout << " ZRAZ " << Zraz << " RSS " << RSS <<" RSS_z " << RSS_z << G4endl;
                
              
                 	if(std::fabs(Zraz)>1*cm){
                			BAD_Tracks++;
                	}
                	if(RSS_z > RSS_z_max){
                		RSS_z_Tracks++;
                		if(std::fabs(Zraz)>1*cm) RSS_z_Tracks_Bad++; 
                	}
                	else{
                	
                		ALL_FIT_Tracks++;
                        		    
                		if(RSS_z < RSS_z_max){
                			analysisManager->FillH1(3, Z_nach);
                			analysisManager->FillH1(4, Z0);
    					analysisManager->FillH2(0, Dphi2, Zraz);
                			Z0_Map[TrackID] = Z_nach; // ID, Znach
                			Z_Map[TrackID] = Z0;  // ID, Z0
                			Z0up_Map[Z0] = TrackID; // Z0, ID 
                		//G4cout << " ZRAZ " << Zraz << " RSS " << RSS <<" RSS_z " << RSS_z << G4endl;
                		}
                
                		//Z0 - reconstructed, Z_nach - truly Z
               
                	}
                  
               }//end of fit    */
             
               		ID_Number = 0;  
               		i_step = 0; 
       
        	} //end of prev!=next TrackID
     
      		break;    
              }  //end of prev!=next chamber number
 
    L_pred = L_min;
    k_pred = k_min;
  
    }  
    
   
  }
 } 
/* 
static std::ofstream stuff_out("Output.txt"); 
static std::ofstream stuff_clust("Clust.txt");
stuff_out << BAD_Tracks <<" "<< RSS_Tracks << " " << RSS_Tracks_Bad << " " << RSS_z_Tracks << " " << RSS_z_Tracks_Bad  <<" " << ALL_FIT_Tracks << " " << ALL_Tracks << std::endl; 
//G4cout << " BAD " << BAD_Tracks << " RSS " << RSS_Tracks << " RSS_BAD " << RSS_Tracks_Bad << " RSS_z " << RSS_z_Tracks << " RSS_z_BAD " << RSS_z_Tracks_Bad  <<" ALL_FIT " << ALL_FIT_Tracks << " ALL " << ALL_Tracks << G4endl; 

G4cout  << " RSS_z " << RSS_z_Tracks << " RSS_z_BAD " << RSS_z_Tracks_Bad  << " RSS = " << RSS_z_max << G4endl; 


//Vertex reconstruction from interval H
G4double H = 11*cm, k0, k1, Z_nach_step = 0, Z0_step, Z_end; //set our step   
 G4int Num_success = 0, Num_all = 0, Num_wrong = 0, step_N = 1, Vert_wrong = 0;
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
                                    
 //Just Printing Maps
 for(G4int i = 0; itp != Z_Map.End(); itp++, i++){
 
     G4cout << i << ") Key - TrackID " << itp->first <<" Z0 = "<< G4BestUnit(itp->second,"Length") << G4endl;
 }   
 
  for(G4int i = 0; it0p != Z0_Map.End(); it0p++, i++){
 
     G4cout << i << ") Key - TrackID " << it0p->first <<" Znach = "<< G4BestUnit(it0p->second,"Length") << G4endl;
 }  
 
  for(G4int i = 0; it0upp != Z0up_Map.End(); it0upp++, i++){
 
     G4cout << i << ") Key - Z0 " << G4BestUnit(it0upp->first,"Length") <<" TrackID = "<< it0upp->second << G4endl;
 }  

 

 //Clustering
Z0_step = itZ0->first;
//G4cout << " Z0_step = " << Z0_step << G4endl;
 while (it_step != Z0up_Map.End() ){

 	test = true;
 	it = Z_Map.Begin(); it0 = Z0_Map.Begin();
 	Z_step = it_step->first-0.3*cm;
// if(Z_step+0.5*cm == Z_last) G4cout << " WARNING " << G4endl;

 	for(G4int i = 0; it != Z_Map.End(); it++, i++, it0++, itZ0++){
 
     		//G4cout << i << ") Key - TrackID " << it->first <<" Z0 = "<< G4BestUnit(it->second,"Length") << " Z_nach =  " << G4BestUnit(it0->second,"Length") << G4endl;
     
     		if(it->second <= H + Z_step && it->second > Z_step){
     			//G4cout << " Z0 in H " << G4endl;
        
     			if(it->second == Z_last) G4cout << " WARNING " << G4endl;
        
        			Z_vector.push_back(it->second);  //fill vector of needed Z0
        			IDZ0.push_back(it->first);       //fill vector of Z0's ID
        
        		if(it->second >= Z0_step)  Z0_step = it->second;
        		//G4cout << " Z0 = " << it->second << " ID_Z0 = " << it->first  << " Z0_step = " << Z0_step << G4endl;
          	} 
     
     		if(it0->second <= H + Z_step && it0->second > Z_step){
     
        		Z_nach_vector.push_back(it0->second);	 //fill vector of needed Z_nach
        		IDZnach.push_back(it0->first); 	         //fill vector of Z_nach's ID
        
        	//G4cout << " Znach in H = " << it0->second << " ID_Znach = " << it0->first << G4endl;
       		}
     	}
 

	it_step = Z0up_Map.Find(Z0_step);
	it_step++;
	//G4cout << " Z_step finally " << it_step->first << G4endl;
 	
 	//work with vector
 	G4double Sum = 0., Norm = Z_vector.size(), Z_norm, Z_result; //Z_norm = averange Z, Z_result = |Z_nach(truly Z0)-Z_norm|
 	G4int ID;
 
 	for (G4int i = 0, j = Z_vector.size(); i < j; i++){
    		Sum+=Z_vector[i];
 	} 
 
 	Z_vector.clear();
 	Z_norm = Sum/Norm;  
 	Z_nach_step = 0;

 	for(G4int i  = 0, j = Z_nach_vector.size(); i<j; i++){
    
     		if(Z_nach_vector[i] > Z_step && Z_nach_vector[i] < Z_step+H){
     			if(test){ Z_nach_step = Z_nach_vector[i]; 
        		test = false;
        
     			}
     		else{
     
     			if(Z_nach_step != Z_nach_vector[i]) {Vert_wrong++;   flag = false; Z_nach_step = Z_nach_vector[i];}
     		}
     
     		}
 
 
 	} 

 
 	if(!flag) Vert_wrong++;
 	else{
 		if(IDZnach!=IDZ0) Num_wrong++; 
 	};

 	IDZnach.clear();
 	IDZ0.clear();
 
 	flag = true;
 
 	Z_nach_vector.clear();
 }// end of while
 
 
 
 it1 = Z0_Map.Begin();
 it0 = Z0_Map.Begin();

 for(G4int i = 1; it0 != Z0_Map.End(); it0++, i++){
     if(i==1) Num_all++;
    // G4cout << i << ") Key - TrackID " << it0->first << " Z_nach " << G4BestUnit(it0->second,"Length") << " All_Num "<< G4endl;
     if(it0!=Z0_Map.Begin()){
        if(it1->second!=it0->second) {Num_all++; }
     it1++;
     }
 } 
 stuff_clust << Num_all << " " << Num_wrong << " " << Vert_wrong << std::endl;
 G4cout << " Vsego = " << Num_all << " Wrong = " << Num_wrong << " Vertex_wrong = " << Vert_wrong << G4endl; */
}   
           
