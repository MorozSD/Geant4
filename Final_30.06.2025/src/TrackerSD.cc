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
#include "ST_ID.h"
#include "IDManager.cpp"
#include <iostream>
#include <chrono>




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
   auto analysisManager = G4AnalysisManager::Instance();
       
        

  G4double edep = aStep->GetTotalEnergyDeposit();
 // analysisManager->FillH1(7, edep);
  
  G4int status = aStep->GetTrack()->GetParentID();
  G4double Energy = aStep->GetTrack()->GetKineticEnergy();
  G4int PDGE = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4double Etot = aStep->GetPreStepPoint()->GetTotalEnergy();
  
  if (status != 0) return false;
  if(PDGE == 11 || PDGE == -11) return false;
  if(Etot < 100*MeV) return false;
  //if(edep <= 2*keV) return false;
  
 
  TrackerHit* newHit = new TrackerHit();
  
  G4ThreeVector Pos = aStep->GetTrack()->GetVertexPosition();
  newHit->SetTrackVertex(Pos);
  
  G4int Track_ID = aStep->GetTrack()->GetTrackID();
  newHit->SetTrackID(Track_ID);
 
  G4ThreeVector Mom = aStep->GetTrack()->GetVertexMomentumDirection();
  newHit->SetMomDir(Mom);
  
  G4int ID_Maybe = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetCopyNo();
  
  //std::string filename = "/home/sveta/Geant4/Projects/LV_19.12/spd_straw.db";
 // IDManager*  MyIDManager = new IDManager(filename);
//  G4cout << " HIT " << G4endl;
 // G4cout << " Layer " << MyIDManager->getST_bar_layer_num(ID_Maybe) << " ID " << ID_Maybe << G4endl;
  
  newHit->SetChamberNb(ID_Maybe);
  //newHit->SetChamberLayer(MyST_ID.get_LayerNum(ID_Maybe));
  
 // G4ThreeVector Centre = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetPhysicalVolume();
  //G4cout << " BEGIN OF HIT " <<" ID " << ID_Maybe << G4endl;
 // G4cout << " ID " << ID_Maybe << G4endl;
 // G4cout << " Volume " << aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName() << G4endl;
 // G4cout << " StepLength " << G4BestUnit(aStep->GetTrack()->GetStepLength(),"Length")   << G4endl;
 // G4cout << " Displacement  " << G4BestUnit(aStep->GetDeltaPosition(),"Length") << G4endl;
  
/*  if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary && aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary){ 
  	G4cout << " WARNING! " << G4endl;
  	G4cout << " Volume " << aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName() << G4endl;
  	//G4cout << " Layer " << MyIDManager->getST_bar_layer_num(ID_Maybe) << " ID " << ID_Maybe << G4endl;
  }	*/
 // else G4cout <<  " Edep " << edep << G4endl;
 // if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary) G4cout << " First step! " << G4endl;
 // if(aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) G4cout << " Last step! " << G4endl;
 // G4cout << " Etot " << G4BestUnit(Etot,"Energy") <<" Ekin " << G4BestUnit(aStep->GetPreStepPoint()->GetKineticEnergy(),"Energy")<< " Edep " << G4BestUnit(edep,"Energy") << G4endl;
 
  G4ThreeVector worldPosition = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector pointGlob = aStep->GetPostStepPoint()->GetPosition();
  newHit->SetPos (pointGlob);
  
  G4ThreeVector point0 = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  newHit->SetLocalPos(point0);
  
  //G4cout << "LocalPos x = " <<  point0.x() << " y = " << point0.y() << " z = " << point0.z() << G4endl;
 // G4cout << "Pos x = " <<  pointGlob.x() << " y = " << pointGlob.y() << " z = " << pointGlob.z() << G4endl;
  
  G4double L = std::sqrt(std::pow(point0.y(),2)+std::pow(point0.x(),2)); 
  
  
  //G4cout << " x " << pointGlob.x() << " y " << pointGlob.y() << " Volume " << ID_Maybe << " TrackID "  << Track_ID << G4endl;
  
  
 // G4cout <<" ID " << ID_Maybe << " L " << G4BestUnit(L,"Length") << G4endl;
//  
  
  G4double time0 = aStep->GetPreStepPoint()->GetGlobalTime();
  newHit->SetTime(time0);
  
  G4int prim_charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();  
  newHit->SetCharge(prim_charge);
  
  newHit->SetTrStat(aStep->GetTrack()->GetParentID());
  
  //G4cout << " Mat " << aStep->GetPreStepPoint()->GetMaterial() << G4endl;
  
  newHit->SetEdep(edep);
          
  fHitsCollection->insert( newHit );
  
 // newHit->Print();


  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......'

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
#ifdef G4MULTITHREADED
	 	static G4Mutex stuffMutex = G4MUTEX_INITIALIZER;
	 	G4AutoLock al(&stuffMutex);
	 	#endif
	 	static std::ofstream stuff_XY("XY_all.csv");
	 	static std::ofstream stuff_Zl("Zl_all.csv");
	 	static std::ofstream stuff_ZT("Z0.txt");
	 	static std::ofstream stuff_TettaZraz("TettaZraz.txt");
	 	static std::ofstream stuff_Z0Sigma("Z0Sigma.txt");
	 	
	 	stuff_Z0Sigma << " Z0; Z_nach; Sigma " << std::endl;
	 	
	 	static bool first = true;
	 	if (first) {
	 		first = false;
			stuff_XY << "#,x,y,Track_ID,Zraz" << std::endl;
			stuff_Zl << "#,l,z,Track_ID,Zraz" << std::endl;
		}
 G4int nofHit = 1;
 G4cout
   << "-------->Hits Collection: in this event they are " << nofHit
   << " hits in the tracker chambers: " <<  G4endl;
   
   
 G4int nofHits = fHitsCollection->entries();
   
 G4int step = 0, k, ID_Number = 0, i_step = 0, k_min = 0, k_pred = 0, Nb;
 G4int nextnumber, prevnumber, prevTrackID, nextTrackID, Ignore_Track = 0, Number, status, No_parabol = 0, All_Track = 0, All_Track_no_fit = 0, All_Track_dphi = 0, All_Track_z = 0;
 G4ThreeVector pointH0, pointH1, pointH, pointGlob;
 G4double L, L_min, L_pred, time_for_Track_Global, E_dep, Popal=0., All = 0., L_fit_Z, L_fit_YX, func, Dtime;
 G4bool No_fit = true;
 
 G4double RSS_z_max = 3.5;
 
 G4int Layer_Number, BAD_Tracks = 0, ALL_FIT_Tracks = 0, ALL_Tracks = 0, RSS_z_Tracks = 0, RSS_Tracks = 0, RSS_Tracks_Bad = 0, RSS_z_Tracks_Bad = 0, Dang = 0, GOOD_Tracks = 0;
 
 
 size_t k_fit = 2;
   
 G4int Array[200] {};
 
 G4MapCache <G4int, G4double> Z0_Map;
 G4MapCache <G4int, G4double> Z_Map;
 G4MapCache <G4double, G4int> Z0up_Map; 
 G4MapCache <G4double, G4double> ZSidma_Map; 
 
 auto analysisManager = G4AnalysisManager::Instance();
 G4double NofH = 1;
/* for(G4int i=1; i<nofHits; i++){
 	if((*fHitsCollection)[i]->GetChamberNb() == (*fHitsCollection)[i-1]->GetChamberNb()) NofH++;
 	else{ 
 		analysisManager->FillH1(0, NofH); 
 		NofH = 1;
 		
 	}
 }*/
 
 
 
 
 for (G4int i=0;i<nofHits;i++)
    {
    if (i==step && (*fHitsCollection)[i]->GetTrackID() == Ignore_Track) step++;
    
    if (i==step){
 
       nextnumber = 0;    
       prevnumber = (*fHitsCollection)[i]->GetChamberNb();
       
       prevTrackID = (*fHitsCollection)[i]->GetTrackID();
       pointH0 = (*fHitsCollection)[i]->GetLocalPos();
       L = std::sqrt(std::pow(pointH0.y(),2)+std::pow(pointH0.x(),2)); 
     //  G4cout << " PREV L " << L << " Volume " << prevnumber << " ID " << prevTrackID << " x " << (*fHitsCollection)[i]->GetPos().x() << " y " << (*fHitsCollection)[i]->GetPos().y()<<  G4endl;
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
    //   pointGlob = (*fHitsCollection)[k]->GetPos();
       itnumber++;
     //  G4cout << " NEXT L " << L <<" Volume " << nextnumber << " ID " << nextTrackID<<" x " << (*fHitsCollection)[k]->GetPos().x() << " y " << (*fHitsCollection)[k]->GetPos().y() << G4endl; 
       // We found the nearest point to the wire, not crossed lines
       
       if(L_min > L) { 
          L_min = L;
          k_min = k;
       } 
        
        //IF ONLY ONE TRACK ALGORITHM WON'T FIT THE TRACK
        if (prevnumber!=nextnumber){        
            
     //в последнем объёме ошибка, ибо нет никакого nextnumber, с которым можно было бы сравнивать. Пока не знаю, как это исправить
        //Layer_Number = (*fHitsCollection)[k-1]->GetChamberLayer();
        
        Array[i_step] = k_pred;
      
        i_step++;
             
        ID_Number++;
       // G4cout << " ID_Number " << ID_Number << " L_min " << L_pred << " x " << (*fHitsCollection)[k_pred]->GetPos().x() << " y " << (*fHitsCollection)[k_pred]->GetPos().y() << G4endl;
       // if(L_pred>4.78*mm) G4cout << " WARNING " << G4endl; 
               
        G4double TIME = 0.;
   
        G4double time = (*fHitsCollection)[k_pred]->GetTime();
     
        auto analysisManager = G4AnalysisManager::Instance();
       
        
  
        if(time<0) G4cout << "WARNING!!!!!" << G4endl;
       
        TIME = time+2.71012+1.21564*L_pred+6.82868*L_pred*L_pred;
      //  G4cout << " L " << L_pred << G4endl;
        
        //analysisManager->FillH1(0, time);
        analysisManager->FillH1(1, L_pred);
  /*       if(L_pred<4.78){
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
  if (potok == 13) analysisManager->FillH1(14, TIME);   */ 
      //  analysisManager->FillH1(2, TIME);
  }      
      
     
     
        step = k;
                     
       if (prevTrackID!=nextTrackID){
       		//G4cout << " ID_Number PREV!=NEXT " << ID_Number << G4endl;
       		
       		Ignore_Track = prevTrackID;
        	ALL_Tracks++;              
        	if(ID_Number <= k_fit){ 
        		ID_Number = 0;  
            		i_step = 0; 
            		break;
        	}
       
         	G4int Array_Size = ID_Number+1;
         	
         	      
         	
         	//G4int chamb, chamb_prev = 0;
         	G4double x_prev, y_prev, X, Y;
         	G4bool Y_of_x = true, X_of_y = true, Another = false;
                
         	#ifdef G4MULTITHREADED
	 	static G4Mutex stuffMutex = G4MUTEX_INITIALIZER;
	 	G4AutoLock al(&stuffMutex);
	 	#endif
	 	static std::ofstream stuff_XYZ("XYZ_all.csv");
	 	/*static std::ofstream stuff_Zl("Zl_all.csv");
	 	static std::ofstream stuff_ZT("Z0.txt");
	 	static std::ofstream stuff_TettaZraz("TettaZraz.txt");
	 	
	 	static bool first = true;
	 	if (first) {
	 		first = false;
			stuff_XY << "#,x,y,Track_ID,Zraz" << std::endl;
			stuff_Zl << "#,l,z,Track_ID,Zraz" << std::endl;
		}*/

		//Block of tabular data
		
		G4int kolvo = 43;
		//Dphi2	
		G4double Mu[kolvo] {0.461599936974952, 0.524539908104167, 0.646253755154867, 0.654391916728016, 0.603819953775322, 0.523558543649374, 0.416596006235955, 0.411131499585062, 0.421448245224905, 0.350975874817229, 0.300058297354443, 0.261484499098923, 0.304096597473404, 0.283202723849372, 0.227766356958947, 0.239329548208845, 0.15326219074772, 0.139027584005739, 0.114987730092669, 0.0951005533289732, 0.0559333710913417, 0.00251816939492753, -0.0539454560559627, -0.0808774197646667, 	0.1264085624708, -0.151485208029522, -0.152903203937455, -0.239881559445364, -0.225858070491956, -0.278772087860684, -0.277880298732313, -0.288707407936538, -0.33261483877855, -0.267082154043321,  -0.370636366470588, -0.448969375264205, -0.429251899001585, -0.432289466361974, -0.443853127120623, -0.498211120083333, -0.466231105324948, -0.46896604559633, -0.471287105021552};
		
		/*G4double Sigma[kolvo] {1.96677460928533, 1.64902585411518, 1.7239328785328, 1.70953481115742, 1.3219108615932, 1.18387353544116, 1.17056805284436, 1.12246790080727, 0.979746433098751, 0.944082158553266, 0.889913711796883, 0.830863998462082, 0.796929235322458, 0.7562245693281, 0.711982624948978, 0.666176298623544, 0.689014418138944, 0.613509808339133, 0.597174344352142, 0.594529486084598, 0.617106453951416, 0.601533787072462, 0.612341408607582, 0.594875720896579, 0.607730176575683, 0.666001161807356, 0.657381318786117, 0.731113932957568, 0.685375079633508, 0.749845144993402, 0.79533530797906, 0.845680785083557, 0.872102064732839, 0.892316819904204, 1.04877332417197, 1.03805250213789, 1.31856466549383, 1.23030260565147, 1.23906682842292, 1.45505452224755, 1.58628635125836, 1.68360121552457, 1.85728290898924};*/
		
		G4double Sigma[kolvo] {1.96677460928533, 1.64902585411518, 1.16125689294877, 1.39474984572059, 1.4705933441111,	0.99887979801592, 0.897260536382459, 0.796730104800771, 0.784639407295803, 0.802894441325801, 0.779293119217916, 0.727206389207282, 0.674671655336421, 0.683767265357199, 0.643649011924103, 0.586311532602805, 0.583250369277866, 0.558099189922226, 0.591171410789075, 0.618871897705848, 0.583892395294992, 0.525057644857675, 0.528902007027798, 0.52485471211376, 0.53651838993862, 0.585742382933096, 0.559168011972711, 0.627122718534413, 0.638911546394693, 0.745846776866031, 0.773203126194302, 0.696053373294705, 0.716453369932465, 0.758500831188466, 0.747592944250591, 1.02557963937608, 1.08734884120276, 1.06170862838038, 1.18315682557754, 1.17191032700903, 1.44433684417546, 1.56732651089181, 1.80819139189277};

		
		G4double Tetta[kolvo] {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6,  CLHEP::pi-0.5};

//Dphi1	
	/*G4double Mu[kolvo] {0.446960273799058, 0.506624292724505, 0.417702725937008, 0.43586353338957, 0.484328831609002, 0.502795627202703, 0.367288555605721, 0.450434817831271, 0.444670361760358, 0.288600676875664, 0.3311611519, 0.276002171317874, 0.287384339537453, 0.246422952605684, 0.185237710961573, 0.180968808903993, 0.129129465931217, 0.0949223306982893, 0.114913366493238, 0.0689700725898058, 0.0279181015715864, -0.00164680013687359, -0.018975959265437, -0.0814850296745843, -0.0855712035065682, -0.132810587091415, -0.162239161569547, -0.203176873881057, -0.178740617312278, -0.246383808343462, -0.255677238525261, -0.270420412134021, -0.286409752183714, -0.360720979784394, -0.336492811873614, -0.438049587156129, -0.414182021822985, -0.427235731699847, -0.397847350395954, -0.528592395547112, -0.331876921232449, -0.385907052728707, -0.381467194873188};
	
		G4double Sigma[kolvo] {1.48476427951446, 1.46599533259675, 1.21327919045905, 1.05896133786136, 1.1860678068091, 1.04806456416195, 1.00692448540734, 0.993598978592564, 0.86302244951148, 0.82446112862165, 0.799209081627177, 0.809945306383072, 0.726132875154105, 0.69661207772245, 0.687498744824917, 0.65144442768418, 0.634374977475298, 0.626203040013838, 0.638228972874061, 0.599233152610835, 0.576206401673821, 0.577412540040871, 0.582045454401953, 0.591330986317751, 0.598308021179905, 0.634421897025579, 0.65701606061697, 0.653002176758027, 0.689373896113274, 0.687472193827243, 0.744099458499859, 0.79456048065324, 0.822889267031301, 0.85975760362293, 0.909570586820723, 0.926307760412399, 1.07167268077484, 1.30259756103903, 1.10711591642693, 1.15945100159221, 1.22543105925621, 1.27570341407025, 1.34503457682597};
		
		G4double Tetta[kolvo] {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6,  CLHEP::pi-0.5};*/
	
               
                //Definition of polar angle
                
                G4ThreeVector Vertex_Mom = (*fHitsCollection)[k_pred]->GetMomDir();
                G4double pz = Vertex_Mom.z();
                G4double Dphi2 = 0., Dphi1 = 0.;
                //if(pz<0) Dphi2 = 3.14159-std::acos(Vertex_Mom.z());
                //if(pz>=0) Dphi2 = std::acos(Vertex_Mom.z());
                pointGlob = (*fHitsCollection)[k_pred]->GetPos();
                G4ThreeVector pointNach = (*fHitsCollection)[k_pred]->GetTrackVertex();
                Dphi2 = std::acos(Vertex_Mom.z());
                Dphi1 = std::acos((pointGlob.z()-pointNach.z())/std::sqrt(std::pow(pointGlob.y(),2)+std::pow(pointGlob.z()-pointNach.z(),2)));
                //Dphi1 = 2*std::atan(pointGlob.y()/(pointGlob.z()-(*fHitsCollection)[k_pred]->GetTrackVertex().z()));
               // G4cout << "Dphi pz " << Dphi2 << " Dphi x,y " << Dphi1 << G4endl;
               // analysisManager->FillH1(7, Dphi1);
              //  analysisManager->FillH1(8, Dphi2);
       		      if(Dphi2>=0.5 && Dphi2<=(CLHEP::pi-0.5)){
                	All_Track_dphi++;
                            
                // IF PARTICLE IS NEGATIVE
                	G4int charge = (*fHitsCollection)[Array[1]]->GetCharge();
                         
                        x_prev = 0.;
                	y_prev = 0.;
              
                      
                	pointH = (*fHitsCollection)[Array[1]]->GetPos();
                	X = pointH.x();
                	Y = pointH.y();
               		
               		G4int Right_Size = 2;
               	
               		for (G4int ii = 1; ii<Array_Size; ii++){
               			pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
               			//G4cout << " Array " << G4endl; 
               			if(ii>1 && std::fabs(pointH.x()-(*fHitsCollection)[Array[ii-2]]->GetPos().x())>0.00001*mm){ Right_Size++;}
               			//if(ii!=1) G4cout << " X " << pointH.x() << " PREV X " << (*fHitsCollection)[Array[ii-2]]->GetPos().x() << " Right_Size " << Right_Size <<G4endl;
               		}
               	
               		if(charge <= 0){
               		
               		
               			for (G4int ii = 1; ii<Array_Size; ii++){
            
            				
                			pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                			X = std::fabs(pointH.x());
                			Y = std::fabs(pointH.y());
                			
          
                		//G4cout << " Cross " << CrossVolume << G4endl;
                        	if(X >= x_prev){
                        	
                        		if( Y >= y_prev){
                        		
                        		
                        			
                        			x_prev = X;
                        			y_prev = Y;
                        		}
                        		else{
                        			Y_of_x = false;
                        			Another = true;
                        			break;
                        		}
                        	
                        	}
                        	else{	
                        		X_of_y = false;
                        		Another = true;
                        		break;     
                        	}               
                		}
                	}
                	else{
                		for (G4int ii = 1; ii<Array_Size; ii++){
            
            				
                			pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                			X = std::fabs(pointH.x());
                			Y = std::fabs(pointH.y());
                		
                        		if(X > x_prev){
                        			if(Y > y_prev){
                        				
                        				x_prev = X;
                        				y_prev = Y;
                        		        }
                        			else{
                        				X_of_y = false;
                        				Another = true;
                        		                break;
                        			}
                        	
                        		}
                        		else{	
                        		
                        			Y_of_x = false;
                        			Another = true;
                        	                break;     
                        		}               
                                                  
                		}
                
                	}
                
            //    G4cout <<  " TrackID " << (*fHitsCollection)[Array[1]]->GetTrackID() <<" Array_Size " << Array_Size << " Right_Size " << Right_Size << " X_of_y " << X_of_y << " Y_of_x " << Y_of_x << " Another " << Another << G4endl;
                //if((*fHitsCollection)[Array[1]]->GetTrackID()==122) G4cout << " x " << X << " x_prev " << x_prev << " y " << Y << " y_prev " << y_prev<< " Ras " << std::fabs(X-x_prev) << " " <<std::fabs(Y-y_prev) <<  G4endl;
                	G4double x[Right_Size] {}, x_0[Right_Size] {};
         		G4double y[Right_Size] {}, y_0[Right_Size] {};
         		G4int LZ_Size = Right_Size-1;
               		G4double z[LZ_Size];
         		G4double erry[Right_Size];  
         		x[0] = 0; y[0] = 0;	
                	G4bool CrossVolume = false;
                	
                	/*G4cout << " Before Cross " << G4endl;
                	for (G4int ii = 1; ii<Array_Size; ii++){
           
                        		G4cout << " x " << (*fHitsCollection)[Array[ii-1]]->GetPos().x() << " y " << (*fHitsCollection)[Array[ii-1]]->GetPos().y() << G4endl; ;
                        }*/		
                	G4int i = 0;
                	for (G4int ii = 1; ii<Array_Size; ii++){
           
                        		pointH = (*fHitsCollection)[Array[ii-1]]->GetPos();
                        		stuff_XYZ << pointH.x() << " " << pointH.y() << " " << pointH.z() << std::endl;
                        		
                        		CrossVolume = false;
                			if(ii!=1 && std::fabs(pointH.x()-(*fHitsCollection)[Array[ii-2]]->GetPos().x())<0.00001*mm){ CrossVolume = true;}
                			//G4cout << " Cross " <<  CrossVolume << G4endl;
                			if(!CrossVolume){ 
                			i++;
                        			if(!Y_of_x || !Another){
                        				Another = false;
                		                       // G4cout << " Y (X) " << G4endl;
               						x[i] = pointH.x();
                        				y[i] = pointH.y(); // function y(x)  
                        				z[i-1] = pointH.z();
                        				erry[i] = 0.; 
                					
                				}
                				else if(!X_of_y){
                	
                					x[i] = pointH.y();
                         				y[i] = pointH.x(); // function x(t) 
                         				z[i-1] = pointH.z();
                         				erry[i] = 0.;   
                
                				}
                			
                			}		
                		}
                      /*  if(Array_Size!=Right_Size){
                        	G4cout << " After Cross " << G4endl;
                		for (G4int ii = 1; ii<Right_Size; ii++){
                	//G4cout << " x " << (*fHitsCollection)[Array[ii-1]]->GetPos().x() << " y " << (*fHitsCollection)[Array[ii-1]]->GetPos().y() << G4endl; ;
           
                        		G4cout << " x " << x[ii] << "; " << (*fHitsCollection)[Array[ii-1]]->GetPos().x() << " y " << y[ii] <<"; " <<(*fHitsCollection)[Array[ii-1]]->GetPos().y() << G4endl; ;
                        	}
                        
                        
                        }*/
                       			
                	
                	
                	     
                  //  G4cout << "Polynomial fit Y(x)!" << G4endl;
              
                    	G4double coefbeta[k_fit+1], coefbeta_inv[k_fit+1], term = 0.;                            // Coefficients of the polynomial
                 
                    	G4double R2Adj_x = 0.,R2Adj_y = 0., R2Adj_z = 0., RSS = 0., RSS_z = 0.;
                    	
                    //	auto start = std::chrono::high_resolution_clock::now(); 
                    	
                    	FitData(x,y,coefbeta,k_fit,erry,Right_Size,prevTrackID,R2Adj_x,RSS);  
                    	
                    	
                               
         	    	if( !Another){ 
         	    	
         	    		//G4cout << " Another " << G4endl;
         	    		for (G4int ii = 1; ii<Right_Size; ii++){
                 			term = y[ii];
                 			y[ii] = x[ii];
                 			x[ii] = term;
                 		}         	    
         	     	FitData(x,y,coefbeta_inv,k_fit,erry,Right_Size,prevTrackID,R2Adj_y,RSS);  
         	    
         	   
         	     		if(R2Adj_x<R2Adj_y){ 
         	    		
         	    			coefbeta[0] = coefbeta_inv[0]; 
         	    			coefbeta[1] = coefbeta_inv[1]; 
         	    			coefbeta[2] = coefbeta_inv[2];      	    		
         	      		}     
         	      		else{
         	      			for (G4int ii = 1; ii<Right_Size; ii++){
                 				term = y[ii];
                 				y[ii] = x[ii];
                 				x[ii] = term;
                 			}	
         	      		}
         	     	} 	
        
        	
               		G4double l[LZ_Size];
               		G4double erry_Z[LZ_Size];
             
               		for (G4int ii = 0; ii<LZ_Size; ii++){
               
               
                     		l[ii] = std::fabs(SetLength(k_fit, coefbeta, x[ii+1]));             
                     		//pointH = (*fHitsCollection)[Array[ii]]->GetPos();
                     		//z[ii] = pointH.z(); // function x(t)
                     		erry_Z[ii] = 0.;
            			//G4cout << "x " << x[ii+1] << " l " << l[ii] << " z " << z[ii] << G4endl;
               		}
               
         //       G4cout << "Polynomial fit Z(l)!" << G4endl;
                
                	size_t k_Z = 1;
                	G4double coefbeta_Z[k_Z+1];                            // Coefficients of the polynomial
                
                	FitData(l,z,coefbeta_Z,k_Z,erry_Z,LZ_Size,prevTrackID,R2Adj_z,RSS_z); 
                	G4double Z0 = coefbeta_Z[0];
                	
                //	auto end = std::chrono::high_resolution_clock::now();
                //    	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    		//	std::cout << "Время выполнения: " << duration.count() << " микросекунд" << std::endl;
                	G4double Z_nach = (*fHitsCollection)[k_pred]->GetTrackVertex().z();
                	G4int TrackID = (*fHitsCollection)[k_pred]->GetTrackID();
                        
                        
                      //  G4cout << " ZRAZ " << Z0-Z_nach << " RSS " << RSS_z << " RSS_max " << RSS_z_max << G4endl;	    
                        
                	if(RSS_z < RSS_z_max){
                	
                	ALL_FIT_Tracks++;
                	
                	for(G4int i = 0; i<kolvo; i++){
                	
                		if( Dphi2>=Tetta[i] && Dphi2<=Tetta[i+1]){
                			//G4cout << "  Write " << Dphi2 << " Tetta [" << Tetta[i] << ", " << Tetta[i+1] << "]" << G4endl; 
                			Z0 = Z0-Mu[i];
                			stuff_Z0Sigma << Z0 << " " << Z_nach << " " << Sigma[i] << std::endl;
                			ZSidma_Map[Z0] = Sigma[i];
                			break;
                		}	
                		//G4cout << i << " Dphi2 " << Dphi2 << " Tetta " << Tetta[i+1] << " Mu " << Mu[i+1]<<" Z0 " << Z0 << G4endl;
                	}
                	//G4cout << " Z0 " << " Dphi2 " << Dphi2 << G4endl;
                	G4double Zraz = Z0-Z_nach;
                	//analysisManager->FillH2(0, Dphi2, Zraz);
    			//stuff_TettaZraz << Dphi1 << " " << Zraz << std::endl;
    			//analysisManager->FillH2(1, Dphi1, Zraz);
                	Z0_Map[TrackID] = Z_nach; // ID, Znach
                	Z_Map[TrackID] = Z0;  // ID, Z0
                	Z0up_Map[Z0] = TrackID; // Z0, ID 
                	
                	//G4cout << " Zraz " << Zraz << G4endl;
                	/*if(coefbeta[0]<1*mm && RSS < 1.0 && Zraz > 3*mm ){
                		G4cout << " Zraz " << Zraz << G4endl;
                		analysisManager->FillH2(0, Dphi2, Zraz);
                	 
                	 }*/ // std::sizeof(x)/sizeof(x[0])
                	 
                	
                	 
                      /*  if(RSS < 1.0 && coefbeta[0] < 1*mm && Right_Size>10){
                        	G4cout << " WARNING !!!!!!! " << " Zraz = " << Zraz << " Z " << Z0 << " Znach " << Z_nach<< G4endl;
                        	for (G4int ii = 0; ii<Right_Size; ii++){
                 			G4cout << " x " << x[ii] << " y " << y[ii] << G4endl; 
                 			
                 		}   
                 		G4cout << " Size " <<  Right_Size << G4endl;  
                 		for (G4int ii = 0; ii<LZ_Size; ii++){
               				G4cout << " l " << l[ii] << " z " << z[ii] << G4endl; 
               			}
               			GOOD_Tracks++;
                 		analysisManager->FillH2(0, Dphi2, Zraz);
                 		stuff_TettaZraz << Dphi2 << " " << Zraz << std::endl;  	    
                        	if(Zraz>3*mm) BAD_Tracks++;
                        }*/
               		
                
              
                 	/*if(std::fabs(Zraz)>1*cm){
                			BAD_Tracks++;
                	}
                	if(RSS_z > RSS_z_max){
                		RSS_z_Tracks++;
                		if(std::fabs(Zraz)>1*cm) RSS_z_Tracks_Bad++; 
                	}*/
                //	else{
                	
                		//ALL_FIT_Tracks++;
                        	
                	/*	if(RSS_z < RSS_z_max){
                			 if (Zraz>5*mm){
                //			 	G4cout << " ZRAZ " << Zraz << " Znach " << Z_nach << " Z0 " << Z0 << G4endl;
                			 	G4bool Cross = false;*/
               /* if(!X_of_y){
                 G4cout << " X(Y) " << G4endl; 
               //  G4cout << " TrackID " << (*fHitsCollection)[Array[0]]->GetTrackID() << G4endl;
                 for (G4int ii = 1; ii<Right_Size; ii++){
                 	if(x[ii]==x[ii-1] && y[ii]==y[ii-1]){ G4cout << " Cross " << G4endl; Cross = true;}
                 	G4cout << " x " << x[ii] << " y " << y[ii] << G4endl;
                      }	
                 
                }  
                if(!Y_of_x){
                
                G4cout << " Y(X) " << G4endl; 
                for (G4int ii = 1; ii<Right_Size; ii++){
                	if(x[ii]==x[ii-1] && y[ii]==y[ii-1]){ G4cout << " Cross " << G4endl; Cross = true;}
                 	
                      }	
                } 
                if(!Another){
                	G4cout << " Another " << G4endl;  
                	for (G4int ii = 1; ii<Right_Size; ii++){
                	if(x[ii]==x[ii-1] && y[ii]==y[ii-1]){ G4cout << " Cross " << G4endl; Cross = true;}
                 	
                      }	
                }*/
            /*    if(!Cross){
                	G4cout << " need to check " << G4endl;
                	for (G4int ii = 1; ii<Right_Size; ii++){
                 		if(ii==Right_Size-1) G4cout << " x " << x[ii] << " y " << y[ii] << G4endl;
                 		else G4cout << " x " << x[ii] << " y " << y[ii] << " z " << z[ii] << " l " << l[ii] << G4endl;
                 	
                      }	
                }*/	
                
              /*  for (G4int ii = 0; ii<ID_Number; ii++){
                	G4cout << " l " << l[ii] << " z " << z[ii] << G4endl; 
            
               		}*/
                     
                
                	
         //       } 
                
                
                
                			//analysisManager->FillH1(3, Z_nach);
                			//analysisManager->FillH1(4, Z0);
    				//	analysisManager->FillH2(0, Dphi2, Zraz);
    					//stuff_TettaZraz << Dphi2 << " " << Zraz << std::endl;
    				//	analysisManager->FillH2(1, Dphi1, Zraz);
    				
    				
    					
                 			analysisManager->FillH2(0, Dphi2, Zraz);
                 			stuff_TettaZraz << Dphi2 << " " << Zraz << std::endl;  	    
                        		if(Zraz>3*mm) BAD_Tracks++;
                			Z0_Map[TrackID] = Z_nach; // ID, Znach
                			Z_Map[TrackID] = Z0;  // ID, Z0
                			Z0up_Map[Z0] = TrackID; // Z0, ID 
                		//G4cout << " ZRAZ " << Zraz << " RSS " << RSS <<" RSS_z " << RSS_z << G4endl;
                	//	}
                
                		//Z0 - reconstructed, Z_nach - truly Z
               
                	//}
                  }
               }//end of fit
       		
       		
       		
       		
       		
                ID_Number = 0;  
    		i_step = 0; 
               } //end of prev!=next TrackID
     
      		break; 
        	   
              }  //end of prev!=next chamber number
              
         		
               		
 
    L_pred = L_min;
    k_pred = k_min;
  
    }  
    
    
   
  }
  static std::ofstream stuff_clust("Clust.txt");
  //stuff_clust << ALL_FIT_Tracks << " " << GOOD_Tracks << " " << BAD_Tracks << std::endl;
//  G4cout << " All_Tracks " << ALL_FIT_Tracks << " GOOD " << GOOD_Tracks << G4endl;

 //G4cout << " CLUSTERISATION " << G4endl;
 // CLUSTERISATION SIGMA 

 //static std::ofstream stuff_clust("Clust.txt");
 
 G4MapCache <G4double, G4double >::iterator it = ZSidma_Map.Begin(),it1 = ZSidma_Map.Begin(), it_step = ZSidma_Map.Begin(), itp = ZSidma_Map.Begin();
 G4MapCache <G4int, G4double >::iterator itpp = Z_Map.Begin();
 G4MapCache <G4int, G4double >::iterator it_nach = Z0_Map.Begin(), it0, itZ0;
 //G4MapCache <G4double, G4int >::iterator itpp = Z_Map.Begin();
 
 std::vector <G4double> Clust_Z0;
 std::vector <G4double> Z_nach_vector;
 std::vector <G4int> IDZ0;
 std::vector <G4int> IDZnach;
 
 G4double Z_nach_step;
 G4int Vert_wrong = 0, Num_wrong = 0, Num_all = 0;
 G4bool flag = true;
 
 //Just Printing Maps
/*for(G4int i = 0; itp != ZSidma_Map.End(); itp++, i++){
 
     G4cout << i << ") Key - Z0 " << G4BestUnit(itp->first,"Length") <<" sigma = "<< itp->second << G4endl;
 }   
 G4cout << " " << G4endl;
 for(G4int i = 0; itpp != Z_Map.End(); itpp++, i++){
 
     G4cout << i << ") Key - TrackID " << itpp->first <<" Z0 = "<< itpp->second << G4endl;
 }  
 G4cout << " " << G4endl;
 for(G4int i = 0; it_nach != Z0_Map.End(); it_nach++, i++){
 
     G4cout << i << ") Key - TrackID " << it_nach->first <<" Z_nach = "<< it_nach->second << G4endl;
 }*/   
 

 while (it_step != ZSidma_Map.End()){

 	Clust_Z0.push_back(it_step->first);
 	
 	IDZ0.push_back(Z0up_Map.Find(it_step->first)->second); 
 	
 	it = ZSidma_Map.Find(it_step->first);
 	it_step++;
 	it1 = ZSidma_Map.Find(it_step->first);
 	itp = ZSidma_Map.Find(it_step->first);
 	itp++;
 	
 	
 	//G4cout << it->first << " " << it1->first  << " " << itp->first << G4endl;
 	
 	for(it1 = ZSidma_Map.Find(it_step->first); it1 != ZSidma_Map.End(); it1++){
 	//G4cout << it->first << " " << it1->first  << " " << itp->first << G4endl;
 	//G4cout << " Z1 + 3 sig " << it->first+3*it->second << " Z2 - 3 sig " << it1->first-3*it1->second << G4endl;
 		if(it->first+3*it->second > it1->first-3*it1->second){
 			Clust_Z0.push_back(it1->first);
 			IDZ0.push_back(Z0up_Map.Find(it1->first)->second); 
 			
 			
 			//G4cout << " itp " << itp->first << G4endl;
 			if(itp == ZSidma_Map.End()){
 				it_step = ZSidma_Map.End();
 				//G4cout << " Special " << G4endl;
 				//G4cout << it_step->first << G4endl;
 				break;
 			}
 			it++; itp++;
 		}
 		else{
 		it_step = ZSidma_Map.Find(it1->first);
 		//G4cout << it_step->first << G4endl;
 		break;
 			
 		}
 		
 		
 	}
 	std::sort(IDZ0.begin(), IDZ0.end());
 	G4int Len = Clust_Z0.size();
 	G4double SumZ0 = 0, NormZ0; 
 	
 	//Z norm
 	for(G4int i = 0; i<Len; i++){
 		SumZ0+=Clust_Z0[i];		 	
 	}
 	
 	NormZ0 = SumZ0/Len;
 	
 	for(it_nach = Z0_Map.Begin(); it_nach != Z0_Map.End(); it_nach++){
 
     	if(it_nach->second >=Clust_Z0[0]-2*ZSidma_Map.Find(Clust_Z0[0])->second && it_nach->second <=Clust_Z0[Len-1]+2*ZSidma_Map.Find(Clust_Z0[Len-1])->second){
     		Z_nach_vector.push_back(it_nach->second);
     		IDZnach.push_back(it_nach->first); 
     		
     	
     	}
 }   
 	//G4cout << " Clust len " << Clust_Z0.size()<< " it_step " << it_step->first << G4endl;
 	
 	
 	
 	if(Z_nach_vector.size()!=0){
 	Z_nach_step = Z_nach_vector[0];
 	
 	//G4cout << " Begin " << " STEP " << Z_nach_step << G4endl;
 	
 	for(G4int i  = 0, j = Z_nach_vector.size(); i<j; i++){
 		//G4cout << " STEP " << Z_nach_step << " Z nach " << Z_nach_vector[i] << G4endl;
    		if(Z_nach_step != Z_nach_vector[i] && std::fabs(Z_nach_step-Z_nach_vector[i])>6*mm) {
    		//G4cout << " Wrong Vert " << std::fabs(Z_nach_step-Z_nach_vector[i]) <<  G4endl;
    		Vert_wrong++;   flag = false; Z_nach_step = Z_nach_vector[i];  stuff_Z0Sigma << " WRONG = " << Z_nach_step << " AND " << Z_nach_vector[i-1]<< std::endl;
    		}
     	} 


 	if(!flag) Vert_wrong++;
 	else{
 		if(IDZnach!=IDZ0) {Num_wrong++; }; 
 	};

	Clust_Z0.clear();
 	Z_nach_vector.clear();
 	IDZnach.clear();
 	IDZ0.clear();
 
 	flag = true;
 
 	}
 	else{
 		Num_wrong++;
 		Clust_Z0.clear();
 		IDZ0.clear();
 	}
 	
 
 } // end of while
 
 itZ0 = Z0_Map.Begin();
 it0 = Z0_Map.Begin();

 for(G4int i = 1; it0 != Z0_Map.End(); it0++, i++){
     if(i==1) Num_all++;
   //  G4cout << i << ") Key - TrackID " << it0->first << " Z_nach " << G4BestUnit(it0->second,"Length") << " All_Num "<< G4endl;
     if(it0!=Z0_Map.Begin()){
        if(itZ0->second!=it0->second) {Num_all++; }
     itZ0++;
     }
 } 
 stuff_clust << ALL_Tracks << " " << All_Track_dphi << " " << ALL_FIT_Tracks << " " << BAD_Tracks << " " << Num_all << " " <<  Num_wrong << " " << Vert_wrong << std::endl;
	G4cout << " All_Tracks " << ALL_Tracks << " Dphi " << All_Track_dphi << " RSS_Z " << ALL_FIT_Tracks << " BAD " << BAD_Tracks << G4endl;
	
// stuff_clust << Num_all << " " << Num_wrong << " " << Vert_wrong << std::endl;
 //G4cout << " Vsego = " << Num_all << " Wrong = " << Num_wrong << " Vertex_wrong = " << Vert_wrong << G4endl;
// stuff_Z0Sigma << " Vsego = " << Num_all << " Wrong = " << Num_wrong << " Vertex_wrong = " << Vert_wrong << std::endl;
 
/* 
//Interval clusterisation
static std::ofstream stuff_out("Output.txt"); 
static std::ofstream stuff_clust("Clust.txt");
stuff_out << BAD_Tracks <<" "<< RSS_Tracks << " " << RSS_Tracks_Bad << " " << RSS_z_Tracks << " " << RSS_z_Tracks_Bad  <<" " << ALL_FIT_Tracks << " " << ALL_Tracks << std::endl; 
//G4cout << " BAD " << BAD_Tracks << " RSS " << RSS_Tracks << " RSS_BAD " << RSS_Tracks_Bad << " RSS_z " << RSS_z_Tracks << " RSS_z_BAD " << RSS_z_Tracks_Bad  <<" ALL_FIT " << ALL_FIT_Tracks << " ALL " << ALL_Tracks << G4endl; 

G4cout  << " RSS_z " << RSS_z_Tracks << " RSS_z_BAD " << RSS_z_Tracks_Bad  << " RSS = " << RSS_z_max << G4endl; 


//Vertex reconstruction from interval H
G4double H = 0.7*cm, k0, k1, Z_nach_step = 0, Z0_step, Z_end; //set our step   
 G4int Num_success = 0, Num_all = 0, Num_wrong = 0, step_N = 1, Vert_wrong = 0;
 G4bool flag = true, test = true;
              
 G4MapCache <G4int, G4double >::iterator it = Z_Map.Begin(), it_ID, itp = Z_Map.Begin();
 G4MapCache <G4int, G4double >::iterator it0 = Z0_Map.Begin(), it0_ID, it1 = Z0_Map.Begin(), it0p = Z0_Map.Begin();
 G4MapCache <G4double, G4int >::iterator itZ0 = Z0up_Map.Begin(), it_step = Z0up_Map.Begin(), it0upp = Z0up_Map.Begin();
 
 std::vector <G4double> Z_vector;
 std::vector <G4double> Z0_vector;
 std::vector <G4double> Z_nach_vector;
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
 G4cout << " Vsego = " << Num_all << " Wrong = " << Num_wrong << " Vertex_wrong = " << Vert_wrong << G4endl;*/ 
//}
 }

