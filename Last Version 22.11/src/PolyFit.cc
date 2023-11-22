/*#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <iomanip>*/
#include "PolyFit.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define STOP 1.0e-8
#define TINY 1.0e-30

G4double incbeta(G4double a, G4double b, G4double x) {
    if (x < 0.0 || x > 1.0) return 1.0/0.0;

    if (a<=0.) {
        G4cout << "Warning: a should be >0";
        return 0.;
    }

    if (b<=0.) {
        G4cout << "Warning: b should be >0";
        return 0.;
    }


    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a+1.0)/(a+b+2.0)) {
        return (1.0-incbeta(b,a,1.0-x)); /*Use the fact that beta is symmetrical.*/
    }

    /*Find the first part before the continued fraction.*/
    const G4double lbeta_ab = std::lgamma(a)+std::lgamma(b)-std::lgamma(a+b);
    const G4double front = std::exp(std::log(x)*a+std::log(1.0-x)*b-lbeta_ab) / a;

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    G4double f = 1.0, c = 1.0, d = 0.0;

    G4int i, m;
    for (i = 0; i <= 200; ++i) {
        m = i/2;

        G4double numerator;
        if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        } else if (i % 2 == 0) {
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
        } else {
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
        }

        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (std::fabs(d) < TINY) d = TINY;
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (std::fabs(c) < TINY) c = TINY;

        const G4double cd = c*d;
        f *= cd;

        /*Check for stop.*/
        if (std::fabs(1.0-cd) < STOP) {
            return front * (f-1.0);
        }
    }

    return 1.0/0.0; /*Needed more loops, did not converge.*/
}

G4double invincbeta(G4double y,G4double alpha, G4double beta) {

    if (y <= 0.) return 0.;
    else if (y >= 1.) return 1.;
    if (alpha<=0.) {
        G4cout << "Warning: alpha should be >0";
        return 0.;
    }

    if (beta<=0.) {
        G4cout << "Warning: beta should be >0";
        return 0.;
    }


    G4double x = 0.5;
    G4double a = 0;
    G4double b = 1;
    G4double precision = 1.e-8; 
    G4double binit = y; 
    G4double bcur = incbeta(alpha,beta,x);
    
    while (std::fabs(bcur-binit)>precision) {
        
        if ((bcur-binit)<0) {
            a = x;
        }
        else {
            b = x;
        }
        x = (a+b)*0.5;
        bcur = incbeta(alpha,beta,x);

        //std::cout << x << "\t" << bcur << "\n";
   

    }

    return x;


}



// Calculate the t value for a Student distribution
// **************************************************************
G4double CalculateTValueStudent(const G4double nu, const G4double alpha) {

    G4double precision = 1.e-5;

    if (alpha<=0. || alpha >= 1.) return 0.;

    G4double x = invincbeta(2.*std::min(alpha,1.-alpha), 0.5*nu, 0.5);
    x = std::sqrt(nu*(1.-x)/x);
    return (alpha >= 0.5? x : -x);
    

}

// Cumulative distribution for Student-t
// **************************************************************
G4double cdfStudent(const G4double nu, const G4double t)
{
    G4double x = nu/(t*t+nu); 
   
    return 1.-incbeta(0.5*nu,0.5,x);
}

// Cumulative distribution for Fisher F
// **************************************************************
G4double cdfFisher(const G4double df1, const G4double df2, const G4double x) {
    G4double y = df1*x/(df1*x+df2);
    return incbeta(0.5*df1,0.5*df2,y);
}

// Initialize a 2D array
// **************************************************************
G4double **Make2DArray(const size_t rows, const size_t cols) {

    G4double **array;

    array = new G4double*[rows];
    for(size_t i = 0; i < rows; i++) {
        array[i] = new G4double[cols];
    }

    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            array[i][j] = 0.;
        }
    }
    
    return array;

}

// Transpose a 2D array
// **************************************************************
G4double **MatTrans(G4double **array, const size_t rows, const size_t cols) {

    G4double **arrayT = Make2DArray(cols,rows);

    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            arrayT[j][i] = array[i][j];
        }
    }
    
    return arrayT;

}

// Perform the multiplication of matrix A[m1,m2] by B[m2,m3]
// **************************************************************
G4double **MatMul(const size_t m1, const size_t m2, const size_t m3, G4double **A, G4double **B) {

    G4double **array = Make2DArray(m1,m3);

    for (size_t i=0; i<m1; i++) {          
        for (size_t j=0; j<m3; j++) {      
            array[i][j]=0.; 
            for (size_t m=0; m<m2; m++) {
                array[i][j]+=A[i][m]*B[m][j];
            } 
        }       
    }
    return array;

}

// Perform the multiplication of matrix A[m1,m2] by vector v[m2,1]
// **************************************************************
void MatVectMul(const size_t m1, const size_t m2, G4double **A, G4double *v, G4double *Av) {

    
    for (size_t i=0; i<m1; i++) {   
        Av[i]=0.;
        for (size_t j=0; j<m2; j++) {
            Av[i]+=A[i][j]*v[j];    
        } 
    }
   

}


// Calculates the determinant of a matrix 
// **************************************************************
G4double determinant(G4double **a, const size_t k) {

    G4double s = 1;
    G4double det = 0.;
    G4double **b = Make2DArray(k,k);
    size_t m;
    size_t n;

    if (k == 1) return (a[0][0]);

    for (size_t c=0; c<k; c++) {

        m = 0;
        n = 0;

        for (size_t i = 0; i < k; i++) {

            for (size_t j = 0; j < k; j++) {

                b[i][j] = 0;

                if (i != 0 && j != c) {

                    b[m][n] = a[i][j];
                    if (n < (k - 2)) {
                        n++;
                    }
                    else
                    {
                        n = 0;
                        m++;
                    }
                }
            }
        }

        det = det + s * (a[0][c] * determinant(b, k - 1));
        s = -1 * s;

    }
   
    return (det);

}


// Perform the 
// **************************************************************
void transpose(G4double **num, G4double **fac, G4double **inverse, const size_t r) {

    G4double **b = Make2DArray(r,r);
    G4double deter;

    for (size_t i=0; i<r; i++) {
        for (size_t j=0; j<r; j++) {
            b[i][j] = fac[j][i];
        }
    }

    deter = determinant(num, r);

    for (size_t i=0; i<r; i++) {
        for (size_t j=0; j<r; j++) {
            inverse[i][j] = b[i][j] / deter;
        }
    }

}

// Calculates the cofactors 
// **************************************************************
void cofactor(G4double **num, G4double **inverse, const size_t f)
{

    G4double **b = Make2DArray(f,f);
    G4double **fac = Make2DArray(f,f);
   
    size_t m;
    size_t n;

    for (size_t q=0; q<f; q++) {

        for (size_t p=0; p<f; p++) {

            m = 0;
            n = 0;

            for (size_t i=0; i<f; i++) {

                for (size_t j=0; j<f; j++) {

                    if (i != q && j != p) {

                        b[m][n] = num[i][j];

                        if (n < (f - 2)) {
                            n++;
                        }
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = std::pow(-1, q + p) * determinant(b, f - 1);
        }
    }

    transpose(num, fac, inverse, f);

}


// Display a matrix 
// **************************************************************
void displayMat(G4double **A, const size_t n, const size_t m)  {
 
    G4cout << "Matrix " << n << " x " << m << G4endl;
    for (size_t i=0; i<n; i++) { 
        for (size_t j=0; j<m; j++) 
            G4cout << A[i][j] << "\t"; 
        G4cout << G4endl; 
    } 
    G4cout << G4endl;

} 

// Calculate the residual sum of squares (RSS)
// **************************************************************
G4double CalculateRSS(const G4double *x, const G4double *y, const G4double *a, G4double **Weights,
const G4bool fixed, const size_t N, const size_t n) {

    G4double r2 = 0.;
    G4double ri = 0.;
    for (size_t i=0; i<N; i++) {
        ri = y[i];
        for (size_t j=0; j<n; j++) {
            ri -= a[j]*std::pow(x[i],j);
        }
        r2 += ri*ri*Weights[i][i];
    }

    return r2;

}

// Calculate the total sum of squares (TSS) 
// **************************************************************
G4double CalculateTSS(const G4double *x, const G4double *y, const G4double *a, G4double **Weights, 
const G4bool fixed, const size_t N, const size_t n) {

    G4double r2 = 0.;
    G4double ri = 0.;
    G4double sumwy = 0.;
    G4double sumweights = 0.;
    size_t begin = 0;
    if (fixed) {
        for (size_t i=begin; i<N; i++) {    
            r2+= y[i]*y[i]*Weights[i][i];
        }
    } else {


        for (size_t i=begin; i<N; i++) {    
            sumwy += y[i]*Weights[i][i];
            sumweights += Weights[i][i];
        }
 
        for (size_t i=begin; i<N; i++) {  
            ri = y[i]-sumwy/sumweights;
            r2 += ri*ri*Weights[i][i];
        }
    }

    return r2;

}

// Calculate coefficient R2 - COD
// **************************************************************
G4double CalculateR2COD(const G4double *x, const G4double *y, const G4double *a, G4double **Weights,
const G4bool fixed, const size_t N, const size_t n) {

    G4double RSS = CalculateRSS(x,y,a,Weights,fixed,N,n);
    G4double TSS = CalculateTSS(x,y,a,Weights,fixed,N,n);
    G4double R2 = 1.-RSS/TSS;

    return R2;

}

// Calculate the coefficient R2 - adjusted
// **************************************************************
G4double CalculateR2Adj(const G4double *x, const G4double *y, const G4double *a, G4double **Weights,
const G4bool fixed,const size_t N, const size_t n) {

    G4double RSS = CalculateRSS(x,y,a,Weights,fixed,N,n);
    G4double TSS = CalculateTSS(x,y,a,Weights,fixed,N,n);

    G4double dferr = N-n;
    G4double dftot = N-1;

    if (fixed) {
        dferr += 1.;
        dftot += 1.;
    }
    
    G4double R2Adj = 1.-(dftot)/(dferr)*RSS/TSS;

    return R2Adj;

}

// Perform the fit of data n data points (x,y) with a polynomial of order k
// **************************************************************
void PolyFit(const G4double *x, G4double *y, const size_t n, const size_t k, const G4bool fixedinter,
const G4double fixedinterval, G4double *beta, G4double **Weights, G4double **XTWXInv) { 
  
    // Definition of variables
    // **************************************************************
    G4double **X = Make2DArray(n,k+1);           // [n,k+1]
    G4double **XT;                               // [k+1,n]
    G4double **XTW;                              // [k+1,n]
    G4double **XTWX;                             // [k+1,k+1]

    G4double *XTWY = new G4double[k+1];
    G4double *Y = new G4double[n];

    size_t begin = 0;
    if (fixedinter) begin = 1;

    // Initialize X
    // **************************************************************
    for (size_t i=0; i<n; i++) { 
        for (size_t j=begin; j<(k+1); j++) {  // begin
          X[i][j]=std::pow(x[i],j);  
        }       
    } 

    // Matrix calculations
    // **************************************************************
    XT = MatTrans(X, n, k+1);                 // Calculate XT
    XTW = MatMul(k+1,n,n,XT,Weights);         // Calculate XT*W
    XTWX = MatMul(k+1,n,k+1,XTW,X);           // Calculate (XTW)*X

    if (fixedinter) XTWX[0][0] = 1.;  
    
    cofactor(XTWX, XTWXInv, k+1);             // Calculate (XTWX)^-1

    for (size_t m=0; m<n; m++) {
        if (fixedinter) {
            Y[m]= y[m]-fixedinterval;
        } 
        else {
            Y[m] = y[m];
        }
    } 
    MatVectMul(k+1,n,XTW,Y,XTWY);             // Calculate (XTW)*Y
    MatVectMul(k+1,k+1,XTWXInv,XTWY,beta);    // Calculate beta = (XTWXInv)*XTWY

    if (fixedinter) beta[0] = fixedinterval;

   /* G4cout << "Matrix X" << G4endl; 
    displayMat(X,n,k+1);

    G4cout << "Matrix XT" << G4endl;
    displayMat(XT,k+1,n);

    G4cout << "Matrix XTW" << G4endl;
    displayMat(XTW,k+1,n);

    G4cout << "Matrix XTWXInv" << G4endl;
    displayMat(XTWXInv,k+1,k+1);*/


}

// Calculate the polynomial at a given x value
// **************************************************************
G4double calculatePoly(const G4double x, const G4double *a, const size_t n) {

    G4double poly = 0.;

    for (size_t i=0; i<n; i++) {
        poly += a[i]*std::pow(x,i);
    }

    return poly;

}

// Calculate and write the confidence bands in a file
// **************************************************************
/*void WriteCIBands(std::string filename, const double *x, const double *coefbeta, double **XTXInv, 
const double tstudentval, const double SE, const size_t n, const size_t k) {

   
    double interval = (x[n-1]-x[0]);
    double x1,y0,y1,y2,y3,y4;
    double xstar[k+1];
    double xprod = 0.;

    ofstream output;
    output.open(filename.c_str());
    
    for (int i=0; i<101; i++) {
        x1 = x[0]+interval/100.*i;
        for (size_t j=0; j<k+1; j++) {
            xstar[j] = pow(x1,j); 
        }
        
        xprod = 0.;
        for (size_t j=0; j<(k+1); j++) {
            for (size_t m=0; m<(k+1); m++) {
                xprod += xstar[m]*xstar[j]*XTXInv[j][m];   
            }   
        }

        y0 = calculatePoly(x1, coefbeta,k+1);
        y1 = y0 - tstudentval*SE*sqrt(xprod);
        y2 = y0 + tstudentval*SE*sqrt(xprod);
        y3 = y0 - tstudentval*SE*sqrt(1+xprod);
        y4 = y0 + tstudentval*SE*sqrt(1+xprod);

        output << x1 << "\t" << y0 << "\t" << y1 << "\t" << y2 << "\t";
        output << y3 << "\t" << y4 << endl;

    }


    output.close();


}*/

// Calculate the weights matrix
// **************************************************************
void CalculateWeights(const G4double *erry, G4double **Weights, const size_t n,
const G4int type) {


    for(size_t i = 0; i < n; i++) {

        switch (type) {
            case 0:
                Weights[i][i] = 1.;
            break;
            case 1:
                Weights[i][i] = erry[i];
            break;
            case 2:
                if (erry[i]>0.) {
                    Weights[i][i] = 1./(erry[i]*erry[i]);
                } 
                else {
                    Weights[i][i] = 0.;
                }
            break;
        }
                
    }

}

// Calculate the standard error on the beta coefficients
// **************************************************************
void CalculateSERRBeta(const G4bool fixedinter, const G4double SE, size_t k, G4double *serbeta, G4double **XTWXInv) {

    size_t begin = 0;
    if (fixedinter) begin = 1;

    serbeta[0] = 0.;  
    for (size_t i=begin; i<(k+1); i++) {   
        serbeta[i] = SE*std::sqrt(XTWXInv[i][i]);
    }

}

// Display the polynomial 
// **************************************************************
void DisplayPolynomial(const size_t k) {

    G4cout << "y = ";
    for (size_t i=0; i<(k+1); i++) { 
        G4cout << "A" << i;
        if (i>0) G4cout << "X";
        if (i>1) G4cout << "^" << i;
        if (i<k) G4cout << " + ";
    }
    G4cout << G4endl << G4endl;

}

// Display the ANOVA test result
// **************************************************************
void DisplayANOVA(const size_t nstar, const size_t k, const G4double TSS, const G4double RSS, G4int J) {

    G4double MSReg = (TSS-RSS)/(k);
    G4double MSE = RSS/(nstar-k); 
    G4double FVal = MSReg/MSE; 
    G4double pFVal = 1.-cdfFisher(k,nstar-k, FVal);
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(J, pFVal);

    G4cout << "ANOVA" << G4endl;
    G4cout << "\tDF\tSum squares\tMean square\tF value\tProb>F" << G4endl;
    G4cout << "Model\t" << k << "\t" << (TSS-RSS) << "\t" << (TSS-RSS)/k << "\t" << FVal << "\t" << pFVal << G4endl;
    G4cout << "Error\t" << nstar-k << "\t" << RSS << "\t" << RSS/(nstar-k) << G4endl;  
    G4cout << "Total\t" << nstar << "\t" << TSS << G4endl << G4endl;  

}


// Display the coefficients of the polynomial
// **************************************************************
void DisplayCoefs(const size_t k, const size_t nstar, const G4double tstudentval, const G4double *coefbeta, const G4double *serbeta) {

    G4double lcibeta;    // Low confidence interval coefficients
    G4double hcibeta;    // High confidence interval coefficients

    G4cout << "Polynomial coefficients" << G4endl;
    G4cout << "Coeff\tValue\tStdErr\tLowCI\tHighCI\tStudent-t\tProb>|t|" << G4endl;

    for (size_t i=0; i<(k+1); i++) { 
        lcibeta = coefbeta[i]-tstudentval*serbeta[i];
        hcibeta = coefbeta[i]+tstudentval*serbeta[i];
        G4cout << "A" << i << "\t"; 
        G4cout << coefbeta[i] << "\t";
        G4cout << serbeta[i] << "\t";
        G4cout << lcibeta << "\t";
        G4cout << hcibeta << "\t";

        if (serbeta[i]>0) {
            G4cout << coefbeta[i]/serbeta[i] << "\t";
            G4cout << 1.-cdfStudent(nstar-k, coefbeta[i]/serbeta[i]);  
        } else {
            G4cout << "-\t-";
        }

        G4cout << G4endl;
    }

}

/*G4double GetStartingPoint()     { return f_Start_P; }; //for z
G4double GetStartingTime()     { return f_Start_T; }; // for x, y*/

G4double SetLength(const size_t k,const G4double *coefbeta, G4double x){ //for x, y
  G4double t, t1, t2;
  
 // There is some problem with complex numbers, so I should taking into account norm C
  G4double C1 = coefbeta[1]+2*coefbeta[2]*x;
 // G4cout << " C = " << C << G4endl; 
  G4double C = std::sqrt(1+C1*C1);
  G4double A = coefbeta[2];
  G4double B = coefbeta[1]; 
  G4double C0 = std::sqrt(1+B*B);
  G4double Aabs = std::abs(coefbeta[2]);
 // G4cout << " C " << C << " A " << A << " B " << B << G4endl;

 // t = (B*C)/(4*A)+(x*C)/2-(std::log(-B-2*Aabs*x+C)+std::log(B-2*Aabs*x+C))/(8*Aabs)+(std::log(1+(-2*Aabs+C)/B)-std::log(1-(-2*Aabs+C)/B))/(8*A);
  t = B/(4*A)*(C-C0)+(x*C)/2 - 1/(8*Aabs)*(std::log(std::fabs((-B-2*Aabs*x+C)*(B-2*Aabs*x+C)/(C0*C0-B*B)))) + 1/(8*A)*(std::log(std::fabs(((B-2*Aabs*x+C)*(B-C0))/((B+2*Aabs*x-C)*(B+C0)))));
  return t;
  }
/*G4double SetStartingPoint(const size_t k, const G4double *coefbeta){ 
  const G4int m = k; 
  G4double KOEF[m];
  for (G4int i = 0; i<m; i++){
    KOEF[i] = coefbeta[i];
   }
   
  return KOEF[0]; }*/	



// Display some statistics values
// **************************************************************
void DisplayStatistics(const size_t n, const size_t nstar, const size_t k, const G4double RSS, const G4double R2,
const G4double R2Adj, const G4double SE){

  
    G4cout << G4endl;
    G4cout << "Statistics" << G4endl;
    G4cout << "Number of points: " << n << G4endl;
    G4cout << "Degrees of freedom: " << nstar-k << G4endl; 
    G4cout << "Residual sum of squares: " << RSS << G4endl;
    G4cout << "R-square (COD): " << R2 << G4endl;
    G4cout << "Adj R-square: " << R2Adj << G4endl;
    G4cout << "RMSE: " << SE << G4endl << G4endl;


}

void FitData(G4double *x,G4double *y,G4double *coefbeta, const size_t k, G4double *erry){
    
             G4cout << "Polynomial fit Y(x)!" << G4endl;
    
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
    
G4cout << " coef = " << coefbeta[0] << " " << coefbeta[1] << " " << coefbeta[2] << G4endl; 
G4cout << " yx = " << y[0] << " " << x[0] << G4endl;                    
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
    
G4cout << " Correct0 " << G4endl;
    // Calculate the coefficients of the fit
    // **************************************************************
   		 PolyFit(x,y,n,k,fixedinter,fixedinterval,coefbeta,Weights,XTWXInv);
G4cout << " Correct1 " << G4endl;
G4cout << " coef = " << coefbeta[0] << " " << coefbeta[1] << " " << coefbeta[2] << G4endl;

    // Calculate related values
    // **************************************************************
   		 G4double RSS = CalculateRSS(x,y,coefbeta,Weights,std::fixed,n,k+1);
   		 G4cout << " Correct2 " << G4endl;
   		 G4double TSS = CalculateTSS(x,y,coefbeta,Weights,fixedinter,n,k+1);
   		 G4cout << " Correct3 " << G4endl;
   		 G4double R2 = CalculateR2COD(x,y,coefbeta,Weights,fixedinter,n,k+1);
   		 G4cout << " Correct4 " << G4endl;
   		 G4double R2Adj = CalculateR2Adj(x,y,coefbeta,Weights,fixedinter,n,k+1);
   		 G4cout << " Correct5 " << G4endl;

    if ((nstar-k)>0) {
        SE = std::sqrt(RSS/(nstar-k)); 
        G4cout << " Correct51 " << G4endl;
      //  tstudentval = std::fabs(CalculateTValueStudent(nstar-k, 1.-0.5*alphaval));   //doesn't work...
        tstudentval = 0;
        G4cout << " Correct6 " << G4endl;
    }
 //   G4cout << "t-student value: " << tstudentval << G4endl << G4endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE,k,serbeta,XTWXInv);
G4cout << " Correct7 " << G4endl;
    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k, nstar, tstudentval, coefbeta, serbeta);   // hz
                 
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
    
    
    
    
    
}


