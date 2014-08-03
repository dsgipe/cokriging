//********************************************
// Specification file for cokriging class
//********************************************
using namespace std;
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
private:
    //resize variables 
    void resize();
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
