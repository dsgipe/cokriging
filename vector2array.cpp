#include "vector2array.h"

//************************************************
void vec2array(vector<vector <double> > Vec,double Array[]){
    //converts a vector to a 1d array
    int counter  = 0;
    for(int ii = 0; ii < Vec.size();ii++){
        for(int jj=0;jj<Vec[ii].size();jj++){
            Array[counter] = Vec[jj][ii];
            counter++;
        }
     }

}
//************************************************
void vec2arrayNonSquare(vector<vector <double> > Vec,double Array[]){
    //converts a vector to a 1d array
    int counter  = 0;
    for(int ii = 0; ii < Vec[0].size();ii++){
        for(int jj=0;jj<Vec.size();jj++){
            Array[counter] = Vec[jj][ii];
            counter++;
        }
     }

}
//************************************************
void vec2array(vector<double > Vec,double Array[]){
    //overloaded function takes a one -d vector to an array
    for(int ii = 0; ii < Vec.size();ii++){
        Array[ii] = Vec[ii];
     }

}
//************************************************
double* vec2array(vector<double > Vec){
    //overloaded function takes a one -d vector to an array
    double* Array = new double[Vec.size()];
    for(int ii = 0; ii < Vec.size();ii++){
        Array[ii] = Vec[ii];
     }
    return Array;
}
//************************************************
vector <vector<double> > VecInverse(vector<vector<double> > Vec){//currently unused
   //invert a square vector
   int m = Vec.size();
   if(Vec[0].size()!=m)
       cout << "Error inverse(vector): Vector is not square\n";
   double A[m*m]; 
   vec2array(Vec,A);
   return Vec;
}
//************************************************
//Vector to array cholesky wrapper
//************************************************
vector<vector<double > > chol(vector<vector<double> > PsiC){
    //wrapper of cholesky used for vectors
    int n = PsiC.size();
    double PsiCarry[n*n];
    double UPsiX[n*n];
    //convert vector to array
    int inc = 0;
    for(int ii = 0;ii<n;ii++){for(int jj=0;jj<n;jj++){PsiCarry[inc] = PsiC[ii][jj];inc++;}}
    vector<vector<double > > UPsiV = PsiC;//Initilize 
    CholeskyVec(n,PsiCarry,UPsiX); 
    inc = 0;
    for(int ii = 0;ii<n;ii++){for(int jj=0;jj<n;jj++){UPsiV[jj][ii] = UPsiX[inc];inc++;}}
    return UPsiV; 
}
//************************************************
void WriteVec(vector<vector<double> > V){
    //Print out variables, currently used for debugging
    for(int ii=0;ii<V.size();ii++){
        cout << endl;
        for(int jj = 0;jj<V[ii].size();jj++){
            cout << setiosflags(ios::fixed) << setprecision(4) << V[ii][jj] << "\t";
         }
    }
}
//************************************************
void CholeskyVec(int d,double*S,double*D){
    //solve cholesky of an array
    //inputs: 
    //    d: size
    //    S: input array
    //    D: output array. Cholesky of input array. 
    for(int k=0;k<d;++k){
       double sum=0.;
       for(int p=0;p<k;++p)sum+=D[k*d+p]*D[k*d+p];
       D[k*d+k]=sqrt(S[k*d+k]-sum);
       for(int i=k+1;i<d;++i){
          double sum=0.;
          for(int p=0;p<k;++p)sum+=D[i*d+p]*D[k*d+p];
          D[i*d+k]=(S[i*d+k]-sum)/D[k*d+k];
       }
    }
}
//************************************************
