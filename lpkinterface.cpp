#include <iostream>
#include "cokriging.h"
#include "lpkinterface.h"
using namespace std;
void inverse(double A[], int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}
double* matrixLeftDivision(double A[], double B[], int N,int NRHS)
{
    int *IPIV = new int[N+1];
    int INFO;
    //Create temporary variables so the fortran code doesn't change
    //input variables
    double A_tmp[N*N];for(int ii =0;ii<N*N;ii++){A_tmp[ii]=A[ii];}
    double* B_rtn= new double[N];for(int ii =0;ii<N;ii++){B_rtn[ii]=B[ii];}

    dgesv_(&N,&NRHS,A_tmp,&N,IPIV,B_rtn,&N,&INFO);
    delete [] IPIV;
    return B_rtn;
}

double* matrixMultiply(double A[], double B[],int M, int N,int K)
{   
    
    //---------------------------------------------//
    // Multiple A by B
    // Inputs: 
    //    A: MxK matrix
    //    B: N+K matrix
    //    M: size of A
    //    N: size of B
    //    K: size shared by A and B
    //  Returns:
    //    C: MxN matrix
    //---------------------------------------------//

    //Create temporary variables so the fortran code doesn't change
    //input variables
    double A_tmp[K*M];for(int ii =0;ii<K*M;ii++){A_tmp[ii]=A[ii];}
    double B_tmp[K*N];for(int ii =0;ii<K*N;ii++){B_tmp[ii]=B[ii];}
    double* C_rtn = new double[M*N];for(int ii =0;ii<<M*N;ii++){C_rtn[ii]=0;}
    char transa = 'n';
    char transb = 't';
    double alpha =1;double beta = 0;
    //cout << "\nIn lpk a:\n"; 
    //Write1Darray(A_tmp,K,M);
    //cout << "In lpk b:\n"; 
    //Write1Darray(B_tmp,K,N);
    //cout << "M " << M<< " N " <<  N<< " K " << K<< endl;
    
    dgemm_(&transa, &transb, &M, &N, &K,&alpha,A_tmp,&M,B_tmp,&N,&beta,C_rtn, &K );
    return C_rtn;


}
