extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    // matrix multiplication
    void dgemm_(char* transa, char*transb,int* M,  int* N, int* K,double* alpha, double* A, int* lda, double *B, int* ldb, double* beta, double* C, int* ldc );
    // transa should be 'N', transb should be 'N'
    // M is the number of rows N is the number of columns
    
    // generate inverse of a matrix given its LU decomposition
    void dgesv_(int* N, int* NRHS, double* A, int* lda, int* IPIV, double *B, int* ldb, int* INFO );
/*  Arguments
  =========

  N       (input) INTEGER
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0.

  NRHS    (input) INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0.

  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the N-by-N coefficient matrix A.
          On exit, the factors L and U from the factorization
          A = P*L*U; the unit diagonal elements of L are not stored.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).

  IPIV    (output) INTEGER array, dimension (N)
          The pivot indices that define the permutation matrix P;
          row i of the matrix was interchanged with row IPIV(i).

  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the N-by-NRHS matrix of right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X.

  LDB     (input) INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).

  INFO    (output) INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
                has been completed, but the factor U is exactly
                singular, so the solution could not be computed.

*/


}


void inverse(double* A, int N)
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
double* matrixLeftDivision(double* A, double* B, int N,int NRHS)
{
    int *IPIV = new int[N+1];
    int INFO;
    //Create temporary variables so the fortran code doesn't change
    //input variables
    double A_tmp[N*N];for(int ii =0;ii<N*N;ii++){A_tmp[ii]=A[ii];}
    double* B_rtn = new double[N];for(int ii =0;ii<N;ii++){B_rtn[ii]=B[ii];}

    dgesv_(&N,&NRHS,A_tmp,&N,IPIV,B_rtn,&N,&INFO);
    delete IPIV;
    return B_rtn;
}

double* matrixMultiply(double* A, double* B,int M, int N,int K)
{
    //Create temporary variables so the fortran code doesn't change
    //input variables
    double A_tmp[K*M];for(int ii =0;ii<K*M;ii++){A_tmp[ii]=A[ii];}
    double B_tmp[K*N];for(int ii =0;ii<K*N;ii++){B_tmp[ii]=B[ii];}
    double* C_rtn = new double[M*N];for(int ii =0;ii<<M*N;ii++){C_rtn[ii]=0;}
    char transa = 'n';
    char transb = 't';
    double alpha =1;double beta = 0;

    dgemm_(&transa, &transb, &M, &N, &K,&alpha,A_tmp,&M,B_tmp,&N,&beta,C_rtn, &K );
    return C_rtn;


}
