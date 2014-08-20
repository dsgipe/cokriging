//********************************************
// Specification file for cokriging class
//********************************************
#include "lpkinterface.h"//Used for class Arr
using namespace std;
class Arr{
public:
    //constructors and destructors
    Arr();
    Arr(const Arr& obj);//copy constructor
    Arr(double* valInit, int m,int n);
    //Initializers
    void Init(double* valInit, int m,int n);
    void Init(int m,int n);
    //Operators
    Arr& operator=(const Arr& obj);
    Arr operator/(const Arr& obj);
    Arr operator*(const Arr& obj);
    ~Arr();
    void print();
    void print(const char * message);
    double * val;
    int M;
    int N;
    Arr transpose();
    Arr cholesky();
};
struct Array{
    double * val;
    int size;
};

Arr mu_num_den(Arr& UPsiX,Arr& Y,Arr& oneN);
void buildPsi(Arr& x, double* theta, Arr& CKPsixRtn);
void Cholesky_arr(const Arr& S,Arr& D);

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
    double* CKPsiXc;
    double* UPsiXc;
    double* CKPsiXe;
    double* UPsiXe;
    double* CKPsiXcXe;
    double* CKPsiXeXc;
    double* UPsiXcXe;
    double* CKPsidXe;
    double* UPsidXe;
    double muc;
    double mud;
    double* d;
    double SigmaSqrc;
    double SigmaSqrd;
    double mu;
    double* UC;
    Arr d_a;
    Arr UC_a;
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
void Print(struct Array arr);
//Arr based routines
