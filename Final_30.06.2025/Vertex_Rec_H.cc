//Vertex reconstruction from interval H
G4double H = 3.0*cm, k0, k1, Z_nach_step = 0, Z0_step, Z_end; //set our step   
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
 
/* for(G4int i = 0; itp != Z_Map.End(); itp++, i++){
 
     G4cout << i << ") Key - TrackID " << itp->first <<" Z0 = "<< G4BestUnit(itp->second,"Length") << G4endl;
 }   
 
  for(G4int i = 0; it0p != Z0_Map.End(); it0p++, i++){
 
     G4cout << i << ") Key - TrackID " << it0p->first <<" Znach = "<< G4BestUnit(it0p->second,"Length") << G4endl;
 }  
 
  for(G4int i = 0; it0upp != Z0up_Map.End(); it0upp++, i++){
 
     G4cout << i << ") Key - Z0 " << G4BestUnit(it0upp->first,"Length") <<" TrackID = "<< it0upp->second << G4endl;
 }*/ 
 

 
 
/*
Z0_step = itZ0->first;
//G4cout << " Z0_step = " << Z0_step << G4endl;
 while (it_step != Z0up_Map.End() ){
 //for(G4int j = 0; j<5; j++){
 test = true;
 it = Z_Map.Begin(); it0 = Z0_Map.Begin();
 //G4cout << " it1 -> " << it1->second << G4endl;
 Z_step = it_step->first-0.3*cm;
// if(Z_step+0.5*cm == Z_last) G4cout << " WARNING " << G4endl;
//if(it!=Z_Map.Begin()){ it--; Z_last = it->second; it++;}
// G4cout << " Z_step = " << G4BestUnit(Z_step,"Length")  << G4endl;
// G4cout << H + Z_step << " flag = " << flag << G4endl;
 for(G4int i = 0; it != Z_Map.End(); it++, i++, it0++, itZ0++){
 
     //G4cout << i << ") Key - TrackID " << it->first <<" Z0 = "<< G4BestUnit(it->second,"Length") << " Z_nach =  " << G4BestUnit(it0->second,"Length") << G4endl;
     
     if(it->second <= H + Z_step && it->second > Z_step){
   //  G4cout << " Z0 in H " << G4endl;
        
     //   if(it->second == Z_last) G4cout << " WARNING " << G4endl;
        
        Z_vector.push_back(it->second);  //fill vector of needed Z0
        IDZ0.push_back(it->first);     //fill vector of Z0's ID
        if(it->second >= Z0_step)  Z0_step = it->second;
      //  G4cout << " Z0 = " << it->second << " ID_Z0 = " << it->first  << " Z0_step = " << Z0_step << G4endl;
        
        
     } 
     
     if(it0->second <= H + Z_step && it0->second > Z_step){
        IDZnach.push_back(it0->first); 
     //   G4cout << " Znach in H = " << it0->second << " ID_Znach = " << it0->first << G4endl;
        Z_nach_vector.push_back(it0->second);
       
     }
 }
*/ 
it_step = Z0up_Map.Find(Z0_step);
it_step++;
///G4cout << " Z_step finally " << it_step->first << G4endl;
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
 G4bool test = true;

 for(G4int i  = 0, j = Z_nach_vector.size(); i<j; i++){
    
  //  G4cout << " Z_nach = " << Z_nach_vector[i] << G4endl; 
    
     if(Z_nach_vector[i] > Z_step && Z_nach_vector[i] < Z_step+H){
   //    G4cout << " Z_nach in H0 = " << Z_nach_vector[i] << " Z_step = " << Z_nach_step << G4endl;  
     if(test){ Z_nach_step = Z_nach_vector[i]; 
    //    G4cout << " Z_step = " << Z_nach_step << G4endl;
        test = false;
        
     }
     else{
     
   //  G4cout << " ELSE " << G4endl; 
     if(Z_nach_step != Z_nach_vector[i]) {Num_wrong++; /*G4cout << " WRONG " << " Wrong_Num = "<<Num_wrong << G4endl; */ flag = false; Z_nach_step = Z_nach_vector[i];}
     }
     
     }
 
 
 } 

 
 
 if(flag){
    
    G4cout << " Array_Test " << G4endl;
    
    for (G4int i = 0, j = IDZnach.size(); i<j; i++){
 //   G4cout << " IDZnach = " << IDZnach[i] << G4endl;
    }
     for (G4int i = 0, j = IDZ0.size(); i<j; i++){
   // G4cout << " IDZ0 = " << IDZ0[i] << G4endl;
    }
    
  
    
 }
 
 if(IDZnach!=IDZ0 || !flag) {Num_wrong++; /*G4cout << " Wrong_Num = " << Num_wrong << G4endl;*/}; 
// G4cout << " ID size = " << IDZ0.size() <<" ID_nach size = " << IDZnach.size() << " flag = " << flag <<" Wrongs = " << Num_wrong<< G4endl;
 IDZnach.clear();
 IDZ0.clear();
 
 
 flag = true;
 
 
 Z_nach_vector.clear();
 }// end of while
 
 */
 
 it1 = Z0_Map.Begin();
 it0 = Z0_Map.Begin();

 for(G4int i = 1; it0 != Z0_Map.End(); it0++, i++){
     if(i==1) {Num_all++; G4cout << " All_Vertex = " << Num_all << G4endl;}
    // G4cout << i << ") Key - TrackID " << it0->first << " Z_nach " << G4BestUnit(it0->second,"Length") << " All_Num "<< G4endl;
     if(it0!=Z0_Map.Begin()){
        if(it1->second!=it0->second) {Num_all++; G4cout << " All_Vertex = " << Num_all << G4endl;}
     it1++;
     }
 } 
 stuff_v << "," << Num_wrong << "," << Num_all << std::endl;
 G4cout << " Vsego = " << Num_all << /*" Wrong = " << Num_wrong <<*/ G4endl;
