// ver. 1.1
#include<iostream>
#include "matrix.hpp"
#include<stdlib.h>

#define TEST_MODE 0

#ifdef _LOCAL_MKL_
void Double_SquareMatrix::Inverse(){
    int m=row, n=column, lda=row;
    int info;
    std::vector<int> ipiv(n);
    int lwork = n;
    std::vector<double>  work(n);
    dgetrf_( &m, &n, AccessArray(), &lda, &ipiv.front(), &info);
    dgetri_( &n, AccessArray(), &lda, &ipiv.front(), &work.front(), &lwork, &info);
}
void Double_Matrix::Product(Double_Matrix& A, Double_Matrix& B){
    char transa='N',transb='N';
	int m=A.row ,n=B.column, k=A.column;
    int lda=m, ldb=k, ldc=m;
	double alpha=1.0,beta=0.0;
	dgemm_(&transa,&transb,&m,&n,&k,&alpha,A.AccessArray(),&lda,B.AccessArray(),&ldb,&beta,AccessArray(),&ldc);
}
void Double_SquareMatrix::Diagonalize(std::vector<double> &EV, std::vector<double> &ES){
    char jobz='V',uplo='U';
    int n=row;
    int lda=n, lwork=3*n-1, info;
    std::vector<double>  work(3*n);
    ES = matrix;
    dsyev_(&jobz,&uplo,&n, &ES.front(), &lda, &EV.front(), &work.front(), &lwork, &info);
    if (info!=0){
        std::cout << "Error at Double_SquareMatrix::Diagonalize in matrix.cpp";
        exit(1);
    }
}

void Complex_SquareMatrix::Inverse(){
    int m=row, n=column, lda=row;
    int info;
    std::vector<int> ipiv(n);
    int lwork = n;
    std::vector<std::complex<double> >  work(n);
    zgetrf_( &m, &n, AccessArray(), &lda, &ipiv.front(), &info);
    zgetri_( &n, AccessArray(), &lda, &ipiv.front(), &work.front(), &lwork, &info);
}
void Complex_Matrix::Product(Complex_Matrix& A, Complex_Matrix& B){
    char transa='N',transb='N';
	int m=A.row ,n=B.column, k=A.column;
    int lda=m, ldb=k, ldc=m;
	std::complex<double> alpha=1.0,beta=0.0;
	zgemm_(&transa,&transb,&m,&n,&k,&alpha,A.AccessArray(),&lda,B.AccessArray(),&ldb,&beta,AccessArray(),&ldc);
}
void Complex_SquareMatrix::Diagonalize(std::vector<double> &EV, std::vector<std::complex<double> > &ES){
    char jobz='V',uplo='U';
    int n=row;
    int lda=n, lwork=2*n-1, info;
    std::vector<std::complex<double> >  work(2*n);
    std::vector<double> rwork(3*n-2);
    ES = matrix;
    zheev_(&jobz,&uplo,&n, &ES.front(), &lda, &EV.front(), &work.front(), &lwork, &rwork.front(), &info);
    if (info!=0){
        std::cout << "Error at Complex_SquareMatrix::Diagonalize in matrix.cpp";
        exit(1);
    }
}
#else
void Double_SquareMatrix::Inverse(){
    int m=row, n=column, lda=row;
    int info;
    std::vector<int> ipiv(n);
    int lwork = n;
    std::vector<double>  work(n);
    dgetrf( &m, &n, AccessArray(), &lda, &ipiv.front(), &info);
    dgetri( &n, AccessArray(), &lda, &ipiv.front(), &work.front(), &lwork, &info);
}
void Double_Matrix::Product(Double_Matrix& A, Double_Matrix& B){
    char transa='N',transb='N';
	int m=A.row ,n=B.column, k=A.column;
    int lda=m, ldb=k, ldc=m;
	double alpha=1.0,beta=0.0;
	dgemm(&transa,&transb,&m,&n,&k,&alpha,A.AccessArray(),&lda,B.AccessArray(),&ldb,&beta,AccessArray(),&ldc);
}
void Double_SquareMatrix::Diagonalize(std::vector<double> &EV, std::vector<double> &ES){
    char jobz='V',uplo='U';
    int n=row;
    int lda=n, lwork=3*n-1, info;
    std::vector<double>  work(3*n);
    ES = matrix;
    dsyev(&jobz,&uplo,&n, &ES.front(), &lda, &EV.front(), &work.front(), &lwork, &info);
    if (info!=0){
        std::cout << "Error at Double_SquareMatrix::Diagonalize in matrix.cpp";
        exit(1);
    }
}

void Complex_SquareMatrix::Inverse(){
    int m=row, n=column, lda=row;
    int info;
    std::vector<int> ipiv(n);
    int lwork = n;
    std::vector<std::complex<double> >  work(n);
    zgetrf( &m, &n, AccessArray(), &lda, &ipiv.front(), &info);
    zgetri( &n, AccessArray(), &lda, &ipiv.front(), &work.front(), &lwork, &info);
}
void Complex_Matrix::Product(Complex_Matrix& A, Complex_Matrix& B){
    char transa='N',transb='N';
	int m=A.row ,n=B.column, k=A.column;
    int lda=m, ldb=k, ldc=m;
	std::complex<double> alpha=1.0,beta=0.0;
	zgemm(&transa,&transb,&m,&n,&k,&alpha,A.AccessArray(),&lda,B.AccessArray(),&ldb,&beta,AccessArray(),&ldc);
}
void Complex_SquareMatrix::Diagonalize(std::vector<double> &EV, std::vector<std::complex<double> > &ES){
    char jobz='V',uplo='U';
    int n=row;
    int lda=n, lwork=2*n-1, info;
    std::vector<std::complex<double> >  work(2*n);
    std::vector<double> rwork(3*n-2);
    ES = matrix;
    zheev(&jobz,&uplo,&n, &ES.front(), &lda, &EV.front(), &work.front(), &lwork, &rwork.front(), &info);
    if (info!=0){
        std::cout << "Error at Complex_SquareMatrix::Diagonalize in matrix.cpp";
        exit(1);
    }
}
#endif


// ------DOUBLE------

// General Matrix
// Constructor
Double_Matrix::Double_Matrix(){ row = column = 0; }
Double_Matrix::Double_Matrix(int n){
    row = column = n;
    matrix.resize(n*n);
}
Double_Matrix::Double_Matrix(int n, int m){
    row = n;
    column = m;
    matrix.resize(n*m);
}
Double_Matrix::Double_Matrix(const Double_Matrix& A){
    row = A.row;
    column = A.column;
    matrix.resize(row*column);
    for (int i=0; i<row*column; i++) matrix[i]=A.matrix[i];
}

// Operator
double Double_Matrix::operator () (int i){ return matrix[i]; }
double Double_Matrix::operator () (int i, int j){ return matrix[i+j*row]; }


// Interface
void Double_Matrix::SetMatrix(int n){
    row = column = n;
    matrix.resize(n*n);
}
void Double_Matrix::SetMatrix(int n, int m){
    row = n;
    column = m;
    matrix.resize(n*m);
}
double* Double_Matrix::AccessArray(){ return &matrix.front(); }
double& Double_Matrix::Element(int i){ return matrix[i]; }
double& Double_Matrix::Element(int i, int j){ return matrix[i + j*row]; }
int Double_Matrix::Row(){ return row; }
int Double_Matrix::Column(){ return column; }

// Functions
void Double_Matrix::PrintMatrix(){
    for(int i=0; i<row; i++){
        for(int j=0; j<column; j++)
            std::cout << matrix[i + j*row] << ",";
        std::cout << std::endl;
    }
}



// Square Matrix
Double_SquareMatrix::Double_SquareMatrix(){ 
    Double_Matrix();
}
Double_SquareMatrix::Double_SquareMatrix(int n){ 
    row = column = n;
    matrix.resize(n*n);
}
Double_SquareMatrix::Double_SquareMatrix(const Double_SquareMatrix& A){ 
    matrix = A.matrix;
    column = A.column;
    row = A.row;
}
int Double_SquareMatrix::Dim(){ return column; }
double Double_SquareMatrix::Trace(){
    double sum=0;
    for(int i=0; i<row; i++)sum += matrix[i + row*i];
    return sum;
}
void Double_SquareMatrix::DiagonalShift(double x){ for(int i=0; i<row; i++)matrix[i + row*i] += x; }








// ---------COMPLEX----------
Complex_Matrix::Complex_Matrix(){ row = column = 0; }
Complex_Matrix::Complex_Matrix(int n){
    row = column = n;
    matrix.resize(n*n);
}
Complex_Matrix::Complex_Matrix(int n, int m){
    row = n;
    column = m;
    matrix.resize(n*m);
}
Complex_Matrix::Complex_Matrix(const Complex_Matrix& A){
    row = A.row;
    column = A.column;
    matrix.resize(row*column);
    for (int i=0; i<row*column; i++) matrix[i]=A.matrix[i];
}

//~ // Operator
//~ inline std::complex<double>& Complex_Matrix::operator [] (int i){ return matrix[i]; }
//~ inline std::complex<double>& Complex_Matrix::operator () (int i, int j){ return matrix[i+j*row]; }
//~ inline Complex_Matrix& Complex_Matrix::operator=(const Complex_Matrix& A){ 
    //~ row = A.row;
    //~ column = A.column;
    //~ matrix = A.matrix;
    //~ return *this;
//~ }
//~ inline Complex_Matrix& Complex_Matrix::operator+=(const Complex_Matrix& A){ 
    //~ if ((column != A.column) || (row != A.row)) throw;
    //~ else {
        //~ matrix += A.matrix;
        //~ return *this;
    //~ }
//~ }
//~ inline Complex_Matrix& Complex_Matrix::operator-=(const Complex_Matrix& A){ 
    //~ if ((column != A.column) || (row != A.row)) throw;
    //~ else {
        //~ matrix -= A.matrix;
        //~ return *this;
    //~ }
//~ }
// Interface
void Complex_Matrix::SetMatrix(int n){
    row = column = n;
    matrix.resize(n*n);
}
void Complex_Matrix::SetMatrix(int n, int m){
    row = n;
    column = m;
    matrix.resize(n*m);
}
std::complex<double>* Complex_Matrix::AccessArray(){ return &matrix.front(); }
std::complex<double>& Complex_Matrix::Element(int i){ return matrix[i]; }
std::complex<double>& Complex_Matrix::Element(int i, int j){ return matrix[i + j*row]; }
int Complex_Matrix::Row(){ return row; }
int Complex_Matrix::Column(){ return column; }


// Functions
void Complex_Matrix::PrintMatrix(){
    for(int i=0; i<row; i++){
        for(int j=0; j<column; j++)
            std::cout << matrix[i + j*row] << ",";
        std::cout << std::endl;
    }
}
void Complex_Matrix::Sub(Complex_Matrix& A){
    for(int i=0; i<row; i++)
        for(int j=0; j<column; j++)
            matrix[i+j*row] -= A.matrix[i+j*row];
}
void Complex_Matrix::Sub(Complex_Matrix& A, Complex_Matrix& B){
    for(int i=0; i<row; i++)
        for(int j=0; j<column; j++)
            matrix[i+j*row] = A.matrix[i+j*row] - B.matrix[i+j*row];
}



// Square Matrix
Complex_SquareMatrix::Complex_SquareMatrix(){ 
    Complex_Matrix();
}
Complex_SquareMatrix::Complex_SquareMatrix(int n){ 
    row = column = n;
    matrix.resize(n*n);
}
Complex_SquareMatrix::Complex_SquareMatrix(const Complex_SquareMatrix& A){ 
    matrix = A.matrix;
    column = A.column;
    row = A.row;
}

int Complex_SquareMatrix::Dim(){ return column; }

std::complex<double> Complex_SquareMatrix::Trace(){
    std::complex<double> sum=0;
    for(int i=0; i<row; i++)sum += matrix[i + row*i];
    return sum;
}
void Complex_SquareMatrix::DiagonalShift(std::complex<double> x){ for(int i=0; i<row; i++)matrix[i + row*i] += x; }


#if TEST_MODE == 1
#define DIM 2
using namespace std;
const std::complex<double> I(0.0,1.0);

int main(void){
    Complex_SquareMatrix a(DIM),b(DIM),c(DIM);
    
    for(int i=0;i<DIM*DIM;i++){
        a.Element(i)=i; 
        b.Element(i)=(double)i-I*(double)i;
    }
    a.Element(1,0) = complex<double>(0.0,-2.0);
    
    cout << "a:" << endl;
    a.PrintMatrix();
    cout << "b:" << endl;
    b.PrintMatrix();
    cout << "c=a*b:" << endl;
    c.Product(a,b);
    c.PrintMatrix();
}
#endif
