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
    cokriging(vector <double> Initxe,vector <double> Initye,vector <double> Initxc,vector <double> Inityc,
             vector<double>InitthetaD,vector<double> InitthetaC,double initrho);
    //Print out results
    void write();
    void buildModel();
    //vector <vector <double> > buildPsi(int n,vector <double> x,vector <double> theta,vector <vector <double> > &UPsiX );
    //double sum(vector <double> x1,vector <double> x2,vector <double> theta,int p,int ii,int jj);
    vector <vector<double> > VecInverse(vector<vector<double> > Vec);
private:
    vector <double> Xe;//expensive independent var    
    vector <double> Ye;//expensive dependent var    
    vector <double> Xc;//cheap independent var    
    vector <double> Yc;//cheap dependent var   
    int ne;//expensive size 
    int nc;//cheap size
    vector <double> thetaD;//
    vector <double> thetaC;//
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
vector <vector <double> > buildPsi(int n,vector <double> x,vector <double> theta );
double sum(vector <double> x1,vector <double> x2,vector <double> theta,int p,int ii,int jj);
vector<vector<double > > chol(vector<vector<double > > PsiC);
double* transpose(double arr[],int n);
double* mu_num_den(double* UPsiX,double* Y,int n,double* oneN);
void Cholesky(int d,double*S,double*D);
