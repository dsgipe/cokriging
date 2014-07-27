//********************************************
// Specification file for cokriging class
//********************************************

#include <vector>

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
    void buildModel();
    //vector <vector <double> > buildPsi(int n,vector <double> x,vector <double> theta,vector <vector <double> > &UPsiX );
    //double sum(vector <double> x1,vector <double> x2,vector <double> theta,int p,int ii,int jj);
    vector <vector<double> > VecInverse(vector<vector<double> > Vec);
private:
    double* Xe;//expensive independent var    
    double* Ye;//expensive dependent var    
    double* Xc;//cheap independent var    
    double* Yc;//cheap dependent var   
    int ne;//expensive size 
    int nc;//cheap size
    double* thetaD;//
    double* thetaC;//
    double rho;
    vector <vector <double> >CKPsiXc;
    vector <vector <double> >UPsiXc;
    vector <vector <double> >CKPsiXe;
    vector <vector <double> >UPsiXe;
    vector <vector <double> >CKPsiXcXe;
    vector <vector <double> >CKPsiXeXc;
    vector <vector <double> >UPsiXcXe;
    vector <vector <double> >CKPsidXe;
    vector <vector <double> >UPsidXe;
    double muc;
    double mud;
    double* d;
    double SigmaSqrc;
    double SigmaSqrd;
    double mu;
};
// Other used functions
vector <vector <double> > buildPsi(int n,double x[],double theta[] );
double sum(double x1[],double x2[],double theta[],int p,int ii,int jj);
vector<vector<double > > chol(vector<vector<double > > PsiC);
double* transpose(double arr[],int n);
double* mu_num_den(double* UPsiX,double* Y,int n,double* oneN);

//Calculates the Cholesky of a matrix S of size dxd ex: (S[d*d])  and returns D
void Cholesky(int d,double*S,double*D);

//Funtion Used to convert a vector to an array
void vec2array(vector<vector <double> > Vec,double Array[]);
void vec2array(vector <double>  Vec,double Array[]);//overloaded function
double* vec2array(vector <double>  Vec);
