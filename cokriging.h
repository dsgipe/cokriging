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
    //Generate cokriging model
    void buildModel();
    void predictor(double* x,int n);
    //destructor
    ~cokriging();
private:
    //resize variables 
    void resize();
    Arr Xe_a, Ye_a, Xc_a, Yc_a;
    Arr Y_a;
    int ne;//expensive size 
    int nc;//cheap size
    double* thetaD;//
    double* thetaC;//
    double rho;
    //Kriging resulting variables using array class
    Arr CKPsiXc_a, CKPsiXe_a, CKPsiXcXe_a, CKPsiXeXc_a, CKPsidXe_a;
    Arr UPsiXc_a, UPsiXe_a, UPsiXcXe_a, UPsidXe_a;
    //Kriging resulting variables
    Arr d_a;
    Arr UC_a;
    Arr muc_a, mud_a, mu_a;
    Arr SigmaSqrc_a, SigmaSqrd_a;
};
// Other used functions kept seperate for information hiding from cokriging class
double* ArraybuildPsi(int n,double* x,double* theta );
double sum(double x1[],double x2[],double theta[],int p,int ii,int jj);

//Arr based routines
Arr mu_num_den(Arr& UPsiX,Arr& Y,Arr& oneN);
void buildPsi(Arr& x, double* theta, Arr& CKPsixRtn);
double sum_pred(Arr x1,Arr x2,double theta[],int p,int ii);
Arr c_pred(double sigma, double rho,Arr x1, Arr x2, double theta[]);
