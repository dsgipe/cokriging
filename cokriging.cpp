//************************************************
//Implementation file cokriging
//based on MATLAB code written by A I J Forrester
//Currently used for 1-D cokriging
//************************************************
#include "cokriging.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <cmath>
#include "lpkinterface.h"
using namespace std;
//************************************************
cokriging::cokriging(){
    //constructor
    Xe.push_back(0);
    Ye.push_back(0);
    Xc.push_back(0);
    Yc.push_back(0);
}
//************************************************
cokriging::cokriging(vector <double> Initxe,vector <double> Initye,vector <double> Initxc,
    vector <double> Inityc,vector<double>InitthetaD,vector<double> InitthetaC,double initrho){
    //intialize input variables
    Xe = Initxe;
    Ye = Initye;
    Xc = Initxc;
    Yc = Inityc;
    thetaD = InitthetaD;
    thetaC = InitthetaC;
    rho = initrho;
    for(int ii = 0;ii < thetaD.size();ii++){thetaD[ii] = pow(10.0,InitthetaD[ii]); }   
    for(int ii = 0;ii < thetaC.size();ii++){thetaC[ii] = pow(10.0,InitthetaC[ii]); }   
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
void Write1Darray(double A[],int m,int n){
    //prints out what is a matrix in row major format
    int counter = 0;
    int b = 0;
    for(int ii = 0; ii < m;ii++){ 
        if(ii ==0){  cout << setw(9);}
        for(int jj = 0; jj < n;jj++){ 
            cout  << setiosflags(ios::fixed) << setprecision(2)<<  A[b+counter*n] ;
            cout << setw(9);
            //cout << setiosflags(ios::fixed) << setprecision(4) << A[b+counter*n]<< "\t";
            counter ++;
            if(counter == m){
               counter =0;
               b++;
               cout << endl;
            }
        }
    }
}
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
vector <vector<double> > cokriging::VecInverse(vector<vector<double> > Vec){//currently unused
   //invert a square vector
   int m = Vec.size();
   if(Vec[0].size()!=m)
       cout << "Error inverse(vector): Vector is not square\n";
   double A[m*m]; 
   vec2array(Vec,A);
   cout << "\nA Not Inverted \n";
   Write1Darray(A,m,m);
   inverse(A,m);
   cout << "\nA inverse - i Hope \n";
   Write1Darray(A,m,m);
   return Vec;
}
//************************************************
double sum(vector <double> x1,vector <double> x2,vector <double> theta,int p,int ii,int jj){
     //sum a vector in a unique way used for constructing cokriging model
     double sumVal = 0;
         //following line is for 2d
     //for(int kk = 0; kk < thetaC.size();kk++){
         //sumVal+= pow(thetaC[kk]*abs(Xc[ii][kk]-Xc[jj][kk]),p);
     //}
     //this is for 1d
    int kk = 0;
    sumVal+= theta[kk]*pow(abs(x1[ii]-x2[jj]),p);
    return sumVal;
}
//************************************************
void cokriging::buildModel(){
    //Main function for developing the cokriging model
    ne = Ye.size();
    nc = Yc.size();
    //initialize cokriging variables
    double UPsiXc_a[nc*nc]; 
    double Yc_a[nc]; 
    double UPsidXe_a[ne*ne]; 
    double Ye_a[ne]; 
    double Y[nc+ne];
    double CKPsiXc_a[nc*nc]; 
    double CKPsiXcXe_a[nc*ne]; 
    double CKPsiXeXc_a[ne*nc]; 
    double CKPsiXe_a[ne*ne]; 
    double CKPsidXe_a[ne*ne];
    double C[(nc+ne)*(nc+ne)];
    double UC[(ne+nc)*(ne+nc)];
    vector<int> one(ne+nc,1);
    vector<double> y = Yc;y.insert(y.end(),Ye.begin(),Ye.end());//concatinate array
    int p = 2;//Curremtly a constant, but could be varied to change kriging differentiation
    double* num; //Temporary variable
    double* den; //Temporary variable
    double* oneNe = new double[ne];for(int ii=0;ii<ne;ii++){oneNe[ii] =1;} //array of ones
    double* oneNc = new double[nc];for(int ii=0;ii<nc;ii++){oneNc[ii] =1;} //array of ones
    double oneNeNc[ne+nc];for(int ii=0;ii<ne+nc;ii++){oneNeNc[ii] =1;} //array of ones
    double* dif = new double[nc];
    double* difd = new double[ne];
    //Build all the various Psi
    CKPsiXc = buildPsi(nc,Xc,thetaC);
    UPsiXc = chol(CKPsiXc);
    CKPsiXe = buildPsi(ne,Xe,thetaC);
    UPsiXe = chol(CKPsiXe);
    CKPsidXe = buildPsi(ne,Xe,thetaD);
    UPsidXe = chol(CKPsidXe);
    CKPsiXcXe.resize(Xc.size(),vector<double>(Xe.size(),0)); 
    for (int ii = 0;ii<Xc.size();ii++){
         for(int jj=0;jj<Xe.size();jj++){
             CKPsiXcXe[ii][jj] = exp(-sum(Xc,Xe,thetaC,p,ii,jj));
         }
    }

    CKPsiXeXc.resize(Xe.size(),vector<double>(Xc.size(),0)); 
    for (int ii = 0;ii<Xc.size();ii++){
         for(int jj=0;jj<Xe.size();jj++){
             CKPsiXeXc[jj][ii] = CKPsiXcXe[ii][jj];
         }
    }
   
    //Compute kriging parameters.
    //convert vector to array for lapack
    // may actually convert all vectors to array from the beginning since I don't think I need them
    
    vec2array(UPsiXc,UPsiXc_a);
    vec2array(Yc,Yc_a);
    cout << "a:\n ";
    Write1Darray(UPsiXc_a,nc,nc);
    
    // solve the rest of the kriging model
    //left it as an array since multi-dimensional may need an array; and it is more convenient
    num = mu_num_den(UPsiXc_a,Yc_a ,nc,oneNc);
    den = mu_num_den(UPsiXc_a,oneNc,nc,oneNc);
    muc = num[0]/den[0];
    d = new double[nc];
    for(int ii = 0;ii<ne;ii++){
        d[ii] = Ye[ii]-rho*Yc[nc-ne+ii];
    } 

    vec2array(Ye,Ye_a);
    vec2array(UPsidXe,UPsidXe_a);
    cout << "\n: UPsixe_a \n";
    Write1Darray(UPsidXe_a,ne,ne);
    
    num = mu_num_den(UPsidXe_a,d,ne,oneNe);
    den = mu_num_den(UPsidXe_a,oneNe,ne,oneNe);
    cout << "\nnum: " << num[0];
    cout << "\nden: " << den[0];
    mud = num[0]/den[0]; 
    for(int ii=0;ii<nc;ii++){
        dif[ii] = Yc[ii]-muc; 
    }
    num = mu_num_den(UPsiXc_a,dif,nc,dif);
    SigmaSqrc = num[0]/nc;
    for(int ii=0;ii<ne;ii++){
        difd[ii] = d[ii]-mud; 
    }
    num = mu_num_den(UPsidXe_a,difd,ne,difd);
    SigmaSqrd = num[0]/ne;
    //construct C
    double C1[nc*nc];
    double C2[nc*ne];
    double C3[ne*nc];
    double C4[ne*ne];
    vec2array(CKPsiXc,CKPsiXc_a);
    vec2arrayNonSquare(CKPsiXcXe,CKPsiXcXe_a);
    vec2arrayNonSquare(CKPsiXeXc,CKPsiXeXc_a);
    vec2array(CKPsidXe,CKPsidXe_a);
    vec2array(CKPsiXe,CKPsiXe_a);
    for(int ii=0;ii<nc*nc;ii++){C1[ii]=SigmaSqrc*CKPsiXc_a[ii];}
    for(int ii=0;ii<ne*nc;ii++){C2[ii]=rho*SigmaSqrc*CKPsiXcXe_a[ii];}
    for(int ii=0;ii<ne*nc;ii++){C3[ii]=rho*SigmaSqrc*CKPsiXeXc_a[ii];}
    for(int ii=0;ii<ne*ne;ii++){C4[ii]=rho*rho*SigmaSqrc*CKPsiXe_a[ii]+SigmaSqrd*CKPsidXe_a[ii];}
    for(int ii=0;ii<(nc+ne)*(nc+ne);ii++){ C[ii]=0; }//Initialize to 0
    int counter = 0;//The 1st quadrant upper left corner
    int rowcounter = 0;
    //Fill C with C1;
    for(int ii=0;ii<nc*nc;ii++){
        C[counter]=C1[ii];
        counter++;
        rowcounter++;
        if (rowcounter ==nc){
            counter += ne;
            rowcounter=0;
        }
    }
    //Fill C with C2
    counter = (nc+ne)*nc;//The 2nd quadrant upper left corner
    rowcounter = 0;
    for(int ii=0;ii<ne*nc;ii++){
        C[counter]=C2[ii];
        counter++;
        rowcounter++;
        if(rowcounter==nc){
            counter+=ne;
            rowcounter=0;
        }
    }
    //Fill C with C3
    counter = nc;//The 3rd quadrant upper left corner
    rowcounter = 0;
    for(int ii=0;ii<ne*nc;ii++){
        C[counter]=C3[ii];
        counter++;
        rowcounter++;
        if(rowcounter==ne){
            counter+=nc;
            rowcounter=0;
        }
    }
    //Fill C with C4
    counter = (nc+ne)*nc+nc;//The 4th quadrant upper left corner
    rowcounter = 0;
    for(int ii=0;ii<ne*ne;ii++){
        C[counter]=C4[ii];
        counter++;
        rowcounter++;
        if(rowcounter==ne){
            counter+=nc;
            rowcounter=0;
        }
    }
    //for(int ii=0;ii<ne*ne;ii++){C4[ii]=SigmaSqrd*CKPsidXe_a[ii];}
    cout <<"\nC1\n"; Write1Darray(C1,nc,nc);
    cout <<"\nC2\n"; Write1Darray(C2,ne,nc);
    cout <<"\nC3\n"; Write1Darray(C3,nc,ne);
    cout <<"\nC4\n"; Write1Darray(C4,ne,ne);
    cout <<"\nC\n"; Write1Darray(C,ne+nc,ne+nc);
    for(int ii = 0;ii< (ne+nc)*(ne+nc);ii++){UC[ii] = 0;}
    Cholesky(ne+nc,C,UC); 
    cout <<"\nUC\n"; Write1Darray(UC,ne+nc,ne+nc);
    for(int ii = 0;ii< nc;ii++){Y[ii] = Yc_a[ii];}
    for(int ii = 0;ii< ne;ii++){Y[ii+nc] = Ye_a[ii];}
    cout <<"\nY\n"; Write1Darray(Y,ne+nc,1);
    num = mu_num_den(UC,Y,ne+nc,oneNeNc);
    //Begin testing here
    den = mu_num_den(UC,oneNeNc,ne+nc,oneNeNc);
    mu = num[0]/den[0];
    cout << "\nmu\n" << mu;
}
//************************************************
double* mu_num_den(double* UPsiX,double* Y,int n,double* oneN){
    // Used to simplify some of the matrix math used in solving for kriging 
    // Inputs:
    //     UPsiX is a matrix of size nxn, 
    //     Y is a array of size nx1 or 1xn; since the array are column arrays the size doesn't matter for a 1 d vector. 
    //     oneN is also an array of size nx1 or 1xn
    // Returns:
    //    array of 1x1 solving the matrix equation One*(UPsiX\(UPsiX\Y))
    double* Solved;
    
    Solved = matrixLeftDivision(UPsiX,
    matrixLeftDivision(transpose(UPsiX,n),Y,n,1),n,1);
    
    return  matrixMultiply(oneN,Solved,1,1,n);
}
//************************************************
double* transpose(double arr[],int n){
    //transpose a row ordered array
    int counter  = 0;
    int b = 0;
    double* tran_arr = new double[n*n];
    for(int ii = 0; ii < n*n;ii++){
        tran_arr[b+counter*n] = arr[ii];
        counter++;
        if(counter == n){
           counter =0;
           b++;
        }
    }
    return tran_arr;

}
//************************************************
void Cholesky(int d,double*S,double*D){
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
    Cholesky(n,PsiCarry,UPsiX); 
    inc = 0;
    for(int ii = 0;ii<n;ii++){for(int jj=0;jj<n;jj++){UPsiV[jj][ii] = UPsiX[inc];inc++;}}
    return UPsiV; 
}
//************************************************
vector <vector <double> > buildPsi(int n,const vector <double> x,const vector <double> theta ){
    //solve for psi, not apart of the cokriging class due to information hiding, I don't want 
    // this function using any private variables.
    // However I should be able to make this more opject oriented, by converting the class to a kriging class
    // instead of a cokriging, that cokriging calls when it creates a cokriging model. Since cokriging is mainly just a 
    // series of kriging functions
    vector<vector<double> > PsiX(n,vector<double>(n,0));//initilize to zeros
    vector<vector<double> > CKPsiX(n,vector<double>(n,0));//initilize to zeros
    int p =2;
    //solve for Psi Cheap
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = ii+1; jj <n;jj++){
            PsiX[ii][jj] = exp(-sum(x,x,theta,p,ii,jj));
        }
    }
    vector<vector<int> > EyeN(n,vector<int>(n,0)); for(int ii = 0; ii<n;ii++){EyeN[ii][ii] = 1;}
    float eps = 2.2204*pow(10,-16);
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = 0; jj <n;jj++){
            CKPsiX[ii][jj] = PsiX[ii][jj] + PsiX[jj][ii]+EyeN[ii][jj]+EyeN[ii][jj]*eps;
        }
    }
    for(int ii=0;ii<x.size();ii++){cout << x[ii] << " ";} cout<<endl;
   
   return CKPsiX;

}
//************************************************
void cokriging::write(){
    // used for debugging
    //Print private variables to screen
    //Print input variables to screen   
    //cout << "\nXc: " ;for(int ii=0; ii<Xc.size();ii++){ cout << Xc[ii]<<" ";}
    //cout << "\nXe: " ;for(int ii=0; ii<Xe.size();ii++){ cout << Xe[ii]<<" ";}
    cout << "\nYc: " ;for(int ii=0; ii<Yc.size();ii++){ cout << Yc[ii]<<" ";}
    //cout << "\nYe: " ;for(int ii=0; ii<Ye.size();ii++){ cout << Ye[ii]<<" " ;}
    //
    //Print calculated variables to screen
    cout << "\nThetaD: " ;for(int ii=0; ii<thetaD.size();ii++){ cout << thetaD[ii]<<" " ;}
    cout << "\nThetaC: " ;for(int ii=0; ii<thetaC.size();ii++){ cout << thetaC[ii]<<" " ;}
    cout << "\nRho: " << rho;//for(int ii=0; ii<rho.size();ii++){ cout << rho[ii]<<" " ;}
    cout << "\nPsiXc: ";
    WriteVec(CKPsiXc); 
    cout << "\nUPsiXc: ";
    WriteVec(UPsiXc); 
    //cout << "\nPsiE: ";
    //WriteVec(CKPsiXe); 
    //cout << "\nUPsiE: ";
    //WriteVec(UPsiXe); 
    //cout << "\nPsiDE: ";
    //WriteVec(CKPsidXe); 
    //cout << "\nUPsiDE: ";
    //WriteVec(UPsidXe);
    //cout << "\nCKPsiXcXe: ";
    //WriteVec(CKPsiXcXe);
    //cout << "\nCKPsiXeXc: ";
    //WriteVec(CKPsiXeXc);
   cout << "\nmuc: "<< muc;
   cout << "\nnd: " ;for(int ii=0; ii<ne;ii++){ cout << d[ii]<<" " ;}
   cout << "\nmud: "<< mud;
   cout << "\nSigmaSqrc: "<< SigmaSqrc;
   cout << "\nSigmaSqrd: "<< SigmaSqrd;
}
//************************************************

