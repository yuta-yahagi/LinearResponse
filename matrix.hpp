// matrix.hpp ver.1.1
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<complex>
#include<vector>

#define MKL_Complex16 std::complex<double> 
#include<mkl.h>


#ifdef _LOCAL_MKL_
extern "C" {
    void dgemm_(char* transa,char* transb,int* m,int* n,int* k,double* alpha,double* a,int* lda,double* b,int* ldb,double* beta,double* c,int* ldc);
    void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
    void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
    void dsyev_(char* jobz,char* uplo,int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
    
    void zgemm_(char* transa,char* transb,int* m,int* n,int* k,std::complex<double>* alpha,std::complex<double>* a,int* lda,std::complex<double>* b
        ,int* ldb,std::complex<double>* beta,std::complex<double>* c,int* ldc);
    void zgetrf_(int* m, int* n, std::complex<double>* a, int* lda, int* ipiv,int* info);
    void zgetri_(int* n, std::complex<double>* a, int* lda, int* ipiv, std::complex<double>* work, int* lwork, int* info);
    void zheev_(char* jobz,char* uplo,int* n, std::complex<double>* a, int* lda, double* w, std::complex<double>* work
        , int* lwork, double* rwork, int* info);
}
#endif

template <typename T>
std::vector<T>& operator+=(std::vector<T> &self, const std::vector<T> &other){
  for(int i = 0; i < (int)self.size(); i++) 
    self[i] += other[i];
    return self;
}
template <typename T>
std::vector<T>& operator-=(std::vector<T> &self, const std::vector<T> &other){
  for(int i = 0; i < (int)self.size(); i++) 
    self[i] += other[i];
    return self;
}

// DOUBLE
class Double_Matrix{
    public:
        Double_Matrix();
        Double_Matrix(int n);
        Double_Matrix(int n, int m);
        Double_Matrix(const Double_Matrix& A);
 
        double operator () (int i);
        double operator () (int i, int j);
        //~ Double_Matrix& operator=(const Double_Matrix& A);
        //~ inline Double_Matrix operator+(const Double_Matrix& A);
        //~ inline Double_Matrix operator-(const Double_Matrix& A);

        void SetMatrix(int n);
        void SetMatrix(int n, int m);
        double* AccessArray();
        double& Element(int i);
        double& Element(int i, int j);
        int Row();
        int Column();

        void PrintMatrix();
        void Product(Double_Matrix& A, Double_Matrix& B);
        
    protected:
        std::vector<double> matrix;
        int row;
        int column;
};
class Double_SquareMatrix : public Double_Matrix {
    public:
        Double_SquareMatrix();
        Double_SquareMatrix(int n);
        Double_SquareMatrix(const Double_SquareMatrix& A);
        int Dim();
        void Inverse();
        double Trace();
        void Diagonalize(std::vector<double> &EV, std::vector<double> &ES);
        void DiagonalShift(double x);
};


// COMPLEX
class Complex_Matrix{
    public:
        Complex_Matrix();
        Complex_Matrix(int n);
        Complex_Matrix(int n, int m);
        Complex_Matrix(const Complex_Matrix& A);
        
        //~ inline std::complex<double>& operator [] (int i);
        //~ inline std::complex<double>& operator () (int i, int j);
        //~ inline Complex_Matrix& operator=(const Complex_Matrix& A);
        //~ inline Complex_Matrix& operator+=(const Complex_Matrix& A);
        //~ inline Complex_Matrix& operator-=(const Complex_Matrix& A);

        void SetMatrix(int n);
        void SetMatrix(int n, int m);
        std::complex<double>* AccessArray();
        std::complex<double>& Element(int i);
        std::complex<double>& Element(int i, int j);
        int Row();
        int Column();
        void PrintMatrix();
        //~ void Add(Complex_Matrix& A);
        //~ void Add(Complex_Matrix& A, Complex_Matrix& B);
        void Sub(Complex_Matrix& A);
        void Sub(Complex_Matrix& A, Complex_Matrix& B);
        void Product(Complex_Matrix& A, Complex_Matrix& B);
    protected:
        std::vector<std::complex<double> > matrix;
        int row;
        int column;
};

class Complex_SquareMatrix : public Complex_Matrix {
    public:
        Complex_SquareMatrix();
        Complex_SquareMatrix(int n);
        Complex_SquareMatrix(const Complex_SquareMatrix& A);
        
        int Dim();
        void Inverse();
        std::complex<double> Trace();
        void Diagonalize(std::vector<double> &EV, std::vector<std::complex<double> > &ES);
        void DiagonalShift(std::complex<double> x);
};

const Complex_Matrix operator + (const Complex_Matrix& A, const Complex_Matrix& B);
const Complex_Matrix operator - (const Complex_Matrix& A, const Complex_Matrix& B);


#endif
