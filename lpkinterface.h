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
//c calls
void inverse(double A[], int N);
double* matrixLeftDivision(double A[], double B[], int N,int NRHS);
double* matrixMultiply(double A[], double B[],int M, int N,int K);

