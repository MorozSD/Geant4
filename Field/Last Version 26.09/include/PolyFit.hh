#ifndef POLYFIT_H
#define POLYFIT_H


//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <math.h>
//#include <cmath>
//#include <iomanip>
#include "globals.hh"
//using namespace std;




G4double incbeta(G4double a, G4double b, G4double x);
G4double invincbeta(G4double y,G4double alpha, G4double beta);
G4double CalculateTValueStudent(const G4double nu, const G4double alpha);
G4double cdfStudent(const G4double nu, const G4double t);
G4double cdfFisher(const G4double df1, const G4double df2, const G4double x);
G4double **Make2DArray(const size_t rows, const size_t cols);
G4double **MatTrans(G4double **array, const size_t rows, const size_t cols);
G4double **MatMul(const size_t m1, const size_t m2, const size_t m3, G4double **A, G4double **B);
void MatVectMul(const size_t m1, const size_t m2, G4double **A, G4double *v, G4double *Av);
G4double determinant(G4double **a, const size_t k);
void transpose(G4double **num, G4double **fac, G4double **inverse, const size_t r);
void cofactor(G4double **num, G4double **inverse, const size_t f);
void displayMat(G4double **A, const size_t n, const size_t m);
G4double CalculateRSS(const G4double *x, const G4double *y, const G4double *a, G4double **Weights,
const G4bool fixed, const size_t N, const size_t n);
G4double CalculateTSS(const G4double *x, const G4double *y, const G4double *a, G4double **Weights, 
const G4bool fixed, const size_t N, const size_t n);
G4double CalculateR2COD(const G4double *x, const G4double *y, const G4double *a, G4double **Weights,
const G4bool fixed, const size_t N, const size_t n);
G4double CalculateR2Adj(const G4double *x, const G4double *y, const G4double *a, G4double **Weights,
const G4bool fixed,const size_t N, const size_t n);
void PolyFit(const G4double *x, G4double *y, const size_t n, const size_t k, const G4bool fixedinter,
const G4double fixedinterval, G4double *beta, G4double **Weights, G4double **XTWXInv);
G4double calculatePoly(const G4double x, const G4double *a, const size_t n);
void WriteCIBands(std::string filename, const G4double *x, const G4double *coefbeta, G4double **XTXInv, 
const G4double tstudentval, const G4double SE, const size_t n, const size_t k);
void CalculateWeights(const G4double *erry, G4double **Weights, const size_t n,
const G4int type);
void CalculateSERRBeta(const G4bool fixedinter, const G4double SE, size_t k, G4double *serbeta, G4double **XTWXInv);
void DisplayPolynomial(const size_t k);
void DisplayANOVA(const size_t nstar, const size_t k, const G4double TSS, const G4double RSS);
void DisplayCoefs(const size_t k, const size_t nstar, const G4double tstudentval, const G4double *coefbeta, const G4double *serbeta);
void DisplayStatistics(const size_t n, const size_t nstar, const size_t k, const G4double RSS, const G4double R2,
const G4double R2Adj, const G4double SE);
void DisplayCovCorrMatrix(const size_t k, const G4double sigma, const G4bool fixed, G4double **XTWXInv);
void SetStartingTime  (G4double **coefbeta);
void SetStartingPoint  (G4double **coefbeta);

#endif
