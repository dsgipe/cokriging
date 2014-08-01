//************************************************
//Implementation file cokriging
//based on MATLAB code written by A I J Forrester
//Currently used for 1-D cokriging
//************************************************
#include "cokriging.h"
#include <iomanip>
#include <iostream>
#include <cmath>
#include "lpkinterface.h"
using namespace std;
//************************************************
cokriging::cokriging(){
    //---------------------------------------------//
    //            constructor
    //---------------------------------------------//
}
//************************************************
cokriging::cokriging(double Initxe[],double Initye[],double Initxc[],double Inityc[],
    double InitthetaD[],double InitthetaC[],double Initrho,
    int Initnc, int Initne){
    
    int numofdim = 1;//will be passed in once I Move to 2D
    //---------------------------------------------//
    //          intialize input variables
    //---------------------------------------------//
    // Higher fidelity - expensive parameters
    ne = Initne;
    Xe = new double[ne]; for(int ii =0;ii<ne;ii++){Xe[ii] = Initxe[ii];}
    Ye = new double[ne]; for(int ii =0;ii<ne;ii++){Ye[ii] = Initye[ii];}
    // Lower fidelity - cheap parameters
    nc = Initnc;
    Xc = new double[nc]; for(int ii =0;ii<nc;ii++){Xc[ii] = Initxc[ii];}
    Yc = new double[nc]; for(int ii =0;ii<nc;ii++){Yc[ii] = Inityc[ii];}
    //linearality control variables
    rho = Initrho;
    thetaD = new double[numofdim];  
    thetaC = new double[numofdim];
    //fill theta array
    for(int ii = 0;ii < numofdim;ii++){
        thetaD[ii] = pow(10.0,InitthetaD[ii]); 
        thetaC[ii] = pow(10.0,InitthetaC[ii]);
    }   
}
//************************************************
void cokriging::resize(){

    //---------------------------------------------//
    //     Change arrays to the correct size 
    //     based on number of data points
    //---------------------------------------------//
    UPsiXc = new double[nc*nc];
    CKPsiXcXe = new double[nc*ne]; for(int ii=0;ii<ne*nc;ii++){CKPsiXcXe[ii] =0;} //array of ones
    CKPsiXeXc = new double[ne*nc]; 
    UPsidXe = new double[ne*ne];
    UPsiXe = new double[ne*ne];
    d = new double[nc];
    UC = new double[(ne+nc)*(ne+nc)];
}
//************************************************
void cokriging::buildModel(){
    //---------------------------------------------//
    //Main function for developing the cokriging model
    //---------------------------------------------//

    //---------------------------------------------//
    // initialize  and declare cokriging variables
    //---------------------------------------------//
    //resize all arrays to fit the data set
    resize();
    int p = 2;//Curremtly a constant, but could be varied to change kriging differentiation
    //---------------------------------------------//
    // initialize  and declare local variables
    //---------------------------------------------//
    double Y[nc+ne];
    double C[(nc+ne)*(nc+ne)];
    double* oneNe = new double[ne];for(int ii=0;ii<ne;ii++){oneNe[ii] =1;} //array of ones
    double* oneNc = new double[nc];for(int ii=0;ii<nc;ii++){oneNc[ii] =1;} //array of ones
    double oneNeNc[ne+nc];for(int ii=0;ii<ne+nc;ii++){oneNeNc[ii] =1;} //array of ones
    double* dif = new double[nc];
    double* difd = new double[ne];
    // used to contruct C
    double C1[nc*nc];//quadrant 1
    double C2[nc*ne];//quadrant 2
    double C3[ne*nc];//quadrant 3
    double C4[ne*ne];//quadrant 4
    //Temporary variables
    double* num;
    double* den; 
    //counters
    int rowcounter = 0;
    int counter = 0;
    int b = 0;
    //---------------------------------------------//
    // Build all the various Psi 
    // variables for cheap and expensive models
    //---------------------------------------------//
    CKPsiXc  = ArraybuildPsi(nc,Xc,thetaC);
    CKPsiXe  = ArraybuildPsi(ne,Xe,thetaC);
    CKPsidXe = ArraybuildPsi(ne,Xe,thetaD);
    //calculate cholesky, 
    //Cholesky(input square size,input square matrix,output square matrix); 
    Cholesky(ne,CKPsiXe,UPsiXe); 
    Cholesky(nc,CKPsiXc,UPsiXc); 
    Cholesky(ne,CKPsidXe,UPsidXe); 
    // fill yet another Psi variable but not square and different results
    counter = 0;b = 0;
    for (int ii = 0;ii<nc;ii++){
         for(int jj=0;jj<ne;jj++){
             CKPsiXcXe[b+counter*nc] = exp(-sum(Xc,Xe,thetaC,p,ii,jj));
             counter++; 
             if(counter == ne){
                counter =0;
                b++;
             }
         }
    }
    CKPsiXeXc = transposeNoneSquare(CKPsiXcXe,nc,ne);
    //Compute kriging parameters.
    
    // solve the rest of the kriging model
    //left it as an array since multi-dimensional may need an array; and it is more convenient
    num = mu_num_den(UPsiXc,Yc ,nc,oneNc);
    den = mu_num_den(UPsiXc,oneNc,nc,oneNc);
    muc = num[0]/den[0];
    //---------------------------------------------//
    //    calculate difference between points, 
    //---------------------------------------------//
    //    *note, it seems strange that Yc is only 
    //    evaluated at some seemingly random points.
    //    This needs to be confirmed with more tests
    //---------------------------------------------//
    //expensive
    for(int ii = 0;ii<ne;ii++){
        d[ii] = Ye[ii]-rho*Yc[nc-ne+ii];
    }  
    //cheap
    num = mu_num_den(UPsidXe,d,ne,oneNe);
    den = mu_num_den(UPsidXe,oneNe,ne,oneNe);
    mud = num[0]/den[0]; 
    for(int ii=0;ii<nc;ii++){
        dif[ii] = Yc[ii]-muc; 
    }
    //difference
    num = mu_num_den(UPsiXc,dif,nc,dif);
    SigmaSqrc = num[0]/nc;
    for(int ii=0;ii<ne;ii++){
        difd[ii] = d[ii]-mud; 
    }
    num = mu_num_den(UPsidXe,difd,ne,difd);
    SigmaSqrd = num[0]/ne;
    //---------------------------------------------//
    //                construct C
    //---------------------------------------------//
    //  C is the culmination of the kriging model,
    //  since it is cokriging, C is formulated of the 
    //  cheap, expensive, and difference models and is the 
    //  main mathematical difference between kriging
    //  and cokriging
    //---------------------------------------------//
    for(int ii=0;ii<nc*nc;ii++){C1[ii]=SigmaSqrc*CKPsiXc[ii];}
    for(int ii=0;ii<ne*nc;ii++){C2[ii]=rho*SigmaSqrc*CKPsiXcXe[ii];}
    for(int ii=0;ii<ne*nc;ii++){C3[ii]=rho*SigmaSqrc*CKPsiXeXc[ii];}
    for(int ii=0;ii<ne*ne;ii++){C4[ii]=rho*rho*SigmaSqrc*CKPsiXe[ii]+SigmaSqrd*CKPsidXe[ii];}
    for(int ii=0;ii<(nc+ne)*(nc+ne);ii++){ C[ii]=0; }//Initialize to 0
    //The 1st quadrant upper left corner
    counter = 0;
    rowcounter = 0;
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
    //The 2nd quadrant upper left corner
    counter = (nc+ne)*nc;
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
    //The 3rd quadrant upper left corner
    counter = nc;
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
    //The 4th quadrant upper left corner
    counter = (nc+ne)*nc+nc;
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
    //---------------------------------------------//
    //         combine Ye and Yc into Y
    //---------------------------------------------//
    for(int ii = 0;ii< (ne+nc)*(ne+nc);ii++){UC[ii] = 0;}
    Cholesky(ne+nc,C,UC); 
    for(int ii = 0;ii< nc;ii++){Y[ii] = Yc[ii];}
    for(int ii = 0;ii< ne;ii++){Y[ii+nc] = Ye[ii];}
    num = mu_num_den(UC,Y,ne+nc,oneNeNc);
    den = mu_num_den(UC,oneNeNc,ne+nc,oneNeNc);
    mu = num[0]/den[0];
    cout << "\nmu\n" << mu;
}
//************************************************
void Write1Darray(double A[],int m,int n){
    //---------------------------------------------//
    //prints out what is a matrix in row major format
    //inputs:
    //     Write A, in columns of m, and rows of n
    //     outputs A to the  screen 
    //---------------------------------------------//
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
double sum(double x1[],double x2[],double theta[],int p,int ii,int jj){
    //---------------------------------------------//
    //    sum a vector in a unique way used for
    //    constructing cokriging model
    //---------------------------------------------//
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
double* mu_num_den(double* UPsiX,double* Y,int n,double* oneN){
    //---------------------------------------------//
    // Used to simplify some of the matrix math used in solving for kriging 
    // Inputs:
    //     UPsiX is a matrix of size nxn, 
    //     Y is a array of size nx1 or 1xn; since the array are column arrays the size doesn't matter for a 1 d vector. 
    //     oneN is also an array of size nx1 or 1xn
    // Returns:
    //    array of 1x1 solving the matrix equation One*(UPsiX\(UPsiX\Y))
    //---------------------------------------------//
    double* Solved;
    
    Solved = matrixLeftDivision(UPsiX,
    matrixLeftDivision(transpose(UPsiX,n),Y,n,1),n,1);
    
    return  matrixMultiply(oneN,Solved,1,1,n);
}
//************************************************
double* transposeNoneSquare(double arr[],int nc,int nr){
    //---------------------------------------------//
    // Transpose a row ordered array, arr,
    // Inputs:
    //     arr: 1d array in column ordered format
    //     nc:  number of columns
    //     nr:  number of rows
    // Returns:
    //     transpose of arr
    //---------------------------------------------//
    int counter  = 0;
    int b = 0;
    double* tran_arr = new double[nc*nr];
    for(int ii = 0; ii < nc*nr;ii++){
        tran_arr[b+counter*nr] = arr[ii];
        counter++;
        if(counter == nc){
           counter =0;
           b++;
        }
    }
    return tran_arr;

}
//************************************************
double* transpose(double arr[],int n){
    //---------------------------------------------//
    //transpose a column ordered array
    //---------------------------------------------//
    // Plan on removing function for the more
    // general none square version
    //---------------------------------------------//
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
    //---------------------------------------------//
    // Solve cholesky of an array
    // Inputs: 
    //    d: size
    //    S: input array
    //    D: output array. Cholesky of input array. 
    //---------------------------------------------//
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
double* ArraybuildPsi(int n,double* x,double* theta ){
    //---------------------------------------------//
    //solve for psi, not apart of the cokriging class due to information hiding, I don't want 
    // this function using any private variables.
    // However I should be able to make this more opject oriented, by converting the class to a kriging class
    // instead of a cokriging, that cokriging calls when it creates a cokriging model. Since cokriging is mainly just a 
    // series of kriging functions
    //---------------------------------------------//
    double PsiX[n][n];//initilize to zeros
    double CKPsiX[n][n];//initilize to zeros
    int EyeN[n][n] ;
    int counter=0;int b; //increment varables
    double* CKPsixRtn = new double[n*n];
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = 0; jj<n;jj++){
            PsiX[ii][jj] = 0;
            CKPsiX[ii][jj] = 0;
            CKPsixRtn[counter]=0;
            EyeN[ii][jj]=0;
            counter++;
        } 
    } 
    // set diagonal to 1;
    for(int ii = 0; ii<n;ii++){EyeN[ii][ii] = 1;}
    int p =2;
    //solve for Psi Cheap
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = ii+1; jj <n;jj++){
            PsiX[ii][jj] = exp(-sum(x,x,theta,p,ii,jj));
        }
    }
    float eps = 2.2204*pow(10,-16);
    counter = 0;b = 0;
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = 0; jj <n;jj++){
            CKPsiX[ii][jj] = PsiX[ii][jj] + PsiX[jj][ii]+EyeN[ii][jj]+EyeN[ii][jj]*eps;
            CKPsixRtn[b+counter*n] = CKPsiX[ii][jj];
            counter++;
            if(counter == n){
               counter =0;
               b++;
            }
        }
    }
   return CKPsixRtn;
}
//************************************************
void cokriging::write(){
    //---------------------------------------------//
    // used for debugging
    // Print private variables to screen
    // Print input variables to screen   
    //---------------------------------------------//
    //cout << "\nXc: " ;for(int ii=0; ii<nc;ii++){ cout << Xc[ii]<<" ";}
    //cout << "\nXe: " ;for(int ii=0; ii<Xe.size();ii++){ cout << Xe[ii]<<" ";}
    //cout << "\nYc: " ;for(int ii=0; ii<Yc.size();ii++){ cout << Yc[ii]<<" ";}
    //cout << "\nYe: " ;for(int ii=0; ii<Ye.size();ii++){ cout << Ye[ii]<<" " ;}
    //
    //Print calculated variables to screen
    //cout << "\nThetaD: " ;for(int ii=0; ii<thetaD.size();ii++){ cout << thetaD[ii]<<" " ;}
    //cout << "\nThetaC: " ;for(int ii=0; ii<thetaC.size();ii++){ cout << thetaC[ii]<<" " ;}
    //cout << "\nRho: " << rho;//for(int ii=0; ii<rho.size();ii++){ cout << rho[ii]<<" " ;}
    //cout << "\nPsiXc: ";
   // cout << "\nUPsiXc: ";
    //cout << "\nPsiE: ";
    //cout << "\nUPsiE: ";
    //cout << "\nPsiDE: ";
    //cout << "\nUPsiDE: ";
    //cout << "\nCKPsiXcXe: ";
    //cout << "\nCKPsiXeXc: ";
   cout << "\nmuc: "<< muc;
   cout << "\nnd: " ;for(int ii=0; ii<ne;ii++){ cout << d[ii]<<" " ;}
   cout << "\nmud: "<< mud;
   cout << "\nSigmaSqrc: "<< SigmaSqrc;
   cout << "\nSigmaSqrd: "<< SigmaSqrd;
}
//************************************************

