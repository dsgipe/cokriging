
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;


vector <vector <double> > buildPsi(int n,vector <double> x,vector <double> theta,vector <vector <double> > &UPsiX );
double sum(vector <double> x1,vector <double> x2,vector <double> theta,int p,int ii,int jj);
vector <vector<double> > VecInverse(vector<vector<double> > Vec);
vector <vector <double> > buildPsi(int n,double x[],double theta[] );
void vec2array(vector<vector <double> > Vec,double Array[]);
void vec2array(vector <double>  Vec,double Array[]);//overloaded function
double* vec2array(vector <double>  Vec);
void CholeskyVec(int d,double*S,double*D);



