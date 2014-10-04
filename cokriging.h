//********************************************
// Specification file for cokriging class
//********************************************
#include "array.h"
class cokriging
{
public:
    //constructors
    cokriging();
    //Fill variables
    cokriging(double Initxe[],double Initye[],double Initxc[],double Inityc[],
             double InitthetaD[],double InitthetaC[],double Initrho,int Initnc,int Initne );
    //Print out results
    void write();
    //Generate cokriging model
    void buildModel();
    void predictor(double* x,int n);
    //destructor
    ~cokriging();
private:
    //resize variables 
    void resize();
    Arr Xe_a;
    Arr Ye_a;
    Arr Xc_a;
    Arr Yc_a;
    Arr Y_a;
    double* Xe;//expensive independent var    
    double* Ye;//expensive dependent var    
    double* Xc;//cheap independent var    
    double* Yc;//cheap dependent var   
    double* Y;
    int ne;//expensive size 
    int nc;//cheap size
    double* thetaD;//
    double* thetaC;//
    double rho;
    //Kriging resulting variables using array class
    Arr CKPsiXc_a;
    Arr UPsiXc_a;
    Arr CKPsiXe_a;
    Arr UPsiXe_a;
    Arr CKPsiXcXe_a;
    Arr CKPsiXeXc_a;
    Arr UPsiXcXe_a;
    Arr CKPsidXe_a;
    Arr UPsidXe_a;
    //Kriging resulting variables
    double* d;
    double SigmaSqrc;
    double SigmaSqrd;
    double mu;
    double* UC;
    Arr d_a;
    Arr UC_a;
    Arr muc_a;
    Arr mud_a;
    Arr mu_a;
    Arr SigmaSqrc_a;
    Arr SigmaSqrd_a;
};
// Other used functions kept seperate for information hiding from cokriging class
double* ArraybuildPsi(int n,double* x,double* theta );
double sum(double x1[],double x2[],double theta[],int p,int ii,int jj);
double* transpose(double arr[],int n);
double* transposeNoneSquare(double arr[],int nc, int nr);
double* mu_num_den(double* UPsiX,double* Y,int n,double* oneN);
//Calculates the Cholesky of a matrix S of size dxd ex: (S[d*d])  and returns D
void Cholesky(int d,double*S,double*D);
double sum_pred(double x1[],double x2[],double theta[],int p,int ii,int n);
double* c_pred(double sigma, double rho,double x1[], int n1, double x2[],int n2, double theta[]);
void Write1Darray(double A[],int m,int n);

//Arr based routines
Arr mu_num_den(Arr& UPsiX,Arr& Y,Arr& oneN);
void buildPsi(Arr& x, double* theta, Arr& CKPsixRtn);
double sum_pred(Arr x1,Arr x2,double theta[],int p,int ii);
Arr c_pred(double sigma, double rho,Arr x1, Arr x2, double theta[]);
