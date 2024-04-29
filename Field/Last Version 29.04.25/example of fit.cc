  G4bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
                  G4int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
                  G4double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
                  G4double alphaval = 0.05;                          // Critical apha value

   // Input values
 // **************************************************************          
             
                
                               
                 // Definition of other variables
    // **************************************************************
    		 size_t n = 0;                                    // Number of data points (adjusted later)
    		 size_t nstar = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    	  	// G4double coefbeta[k_fit+1];                            // Coefficients of the polynomial
   		 G4double serbeta[k_fit+1];                             // Standard error on coefficients
   		 G4double tstudentval = 0.;                         // Student t value
   		 G4double SE = 0.;                                  // Standard error
    
    		 G4double **XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
   		 G4double **Weights;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
    		 n = sizeof(x)/sizeof(G4double);
    		 G4cout << " n = " << n << G4endl;
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     // if(Zraz>2*cm){ 
             /*    if(x[1]*Versh>0 && std::fabs(Versh)<std::fabs(x[Array_Size-1])){
                     
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
                        }*/
               //  }  
               
               
               
               
               
                       
                  
                   
            //       G4cout << x[ii+1]<< " " << y[ii+1] << G4endl;
                   //  G4cout << " t = " << G4BestUnit(SetLength(k, coefbeta, x[ii+1]),"Length") << " Co = " << G4BestUnit(SetLength(k, coefbeta, 0.),"Length") << G4endl;
         //            G4cout << " l = " << G4BestUnit(l[ii],"Length") << " z = " << G4BestUnit(z[ii],"Length") << G4endl;
                /*     Nb = (*fHitsCollection)[Array[ii]]->GetChamberNb()/1000;
                     func = std::fabs(y[ii+1]-(coefbeta[2]*x[ii+1]*x[ii+1]+coefbeta[1]*x[ii+1]+coefbeta[0]));
                     G4cout << " Nb = " << Nb << " funcYX = " << G4BestUnit(func,"Length") << G4endl;
                     
                     analysisManager->FillH2(0, Nb, func);*/
                     
                     
                     
                           G4bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
                G4int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
                G4double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
                G4double alphaval = 0.05;                          // Critical apha value
               
               
                size_t n_Z = 0;                           // Number of data points (adjusted later)
    		size_t nstar_Z = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    		
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
   		DisplayPolynomial(k_Z);

    // Display polynomial coefficients
    // **************************************************************
    		DisplayCoefs(k_Z, nstar_Z, tstudentval_Z, coefbeta_Z, serbeta_Z);
