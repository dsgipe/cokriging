//************************************************
//Implementation file cokriging
//based on MATLAB code written by A I J Forrester
//Currently used for 1-D cokriging
//************************************************
#include "cokriging.h"
#include "cstddef"//for null
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
cokriging::~cokriging(){
    //---------------------------------------------//
    //            destructor
    //---------------------------------------------//
    delete [] Xe;//expensive independent var    
    delete [] Ye;//expensive dependent var    
    delete [] Xc;//cheap independent var    
    delete [] Yc;//cheap dependent var   
    delete [] Y;
    delete [] thetaD;//
    delete [] thetaC;//
    //delete [] UPsiXcXe;
    delete [] d;
    delete [] UC;
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
    Xe_a.Init(Initxe,ne,1);
    Ye_a.Init(Initye,ne,1);
    // Lower fidelity - cheap parameters
    nc = Initnc;
    Xc_a.Init(Initxc,nc,1);
    Yc_a.Init(Inityc,nc,1);
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
    //Array class initialization
    double *ncnc_zero= new double[nc*nc];
    for(int ii=0;ii<nc*nc;ii++){ncnc_zero[ii] =0;} //initialize
    double *ncne_zero = new double[nc*ne];  
    for(int ii=0;ii<nc*ne;ii++){ncne_zero[ii] =0;} //initialize
    double *nene_zero = new double[ne*ne];  
    for(int ii=0;ii<ne*ne;ii++){nene_zero[ii] =0;} //initialize
    UPsiXc_a.Init(ncnc_zero,nc,nc);
    CKPsiXcXe_a.Init(ncne_zero,nc,ne);
    CKPsiXeXc_a.Init(ncne_zero,ne,nc);
    UPsidXe_a.Init(nene_zero,ne,ne); 
    UPsiXe_a.Init(nene_zero,ne,ne);
    CKPsiXc_a.Init(ncnc_zero,nc,nc) ;
    CKPsiXe_a.Init(nene_zero,ne,ne) ;
    CKPsidXe_a.Init(nene_zero,ne,ne);
    delete [] ncnc_zero;
    delete [] ncne_zero;  
    delete [] nene_zero;  
    // old array stle initialization
    d = new double[nc]; for(int ii=0;ii<nc;ii++){d[ii] =0;} //initialize
    UC = new double[(ne+nc)*(ne+nc)]; for(int ii=0;ii<(ne+nc)*(ne+nc);ii++){UC[ii] =0;} //initialize
    Y = new double[nc+ne]; for(int ii=0;ii<ne+nc;ii++){Y[ii] =0;} //initialize
    UC_a.Init(UC,ne+nc,ne+nc);
    Y_a.Init(Y,ne+nc,1);
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
    double C[(nc+ne)*(nc+ne)];
    double oneNe[ne];for(int ii=0;ii<ne;ii++){oneNe[ii] =1;} //array of ones
    double oneNc[nc];for(int ii=0;ii<nc;ii++){oneNc[ii] =1;} //array of ones
    double oneNeNc[ne+nc];for(int ii=0;ii<ne+nc;ii++){oneNeNc[ii] =1;} //array of ones
    Arr oneNe_a(oneNe,ne,1);
    Arr oneNc_a(oneNc,nc,1);
    double dif[nc];
    double difd[ne];
    // used to contruct C
    double C1[nc*nc];//quadrant 1
    double C2[nc*ne];//quadrant 2
    double C3[ne*nc];//quadrant 3
    double C4[ne*ne];//quadrant 4
    Arr C1_a;C1_a.Init(nc,nc);//quadrant 1
    Arr C2_a;C2_a.Init(nc,ne);//quadrant 2
    Arr C3_a;C3_a.Init(ne,nc);//quadrant 3
    Arr C4_a;C4_a.Init(ne,ne);//quadrant 4
    Arr C_a; C_a.Init(ne+nc,ne+nc);
    //counters
    int rowcounter = 0;
    int counter = 0;
    int b = 0;
    //---------------------------------------------//
    // Build all the various Psi 
    // variables for cheap and expensive models
    //---------------------------------------------//
    buildPsi(Xc_a,thetaC,CKPsiXc_a);
    buildPsi(Xe_a,thetaC,CKPsiXe_a);
    buildPsi(Xe_a,thetaD,CKPsidXe_a);
    
    //calculate cholesky, 
    UPsiXe_a  = CKPsiXe_a.cholesky(); 
    UPsiXc_a  = CKPsiXc_a.cholesky(); 
    UPsidXe_a = CKPsidXe_a.cholesky(); 
    //fill yet another Psi variable but not square and different results
    counter = 0;b = 0;
    CKPsiXcXe_a.M = ne;
    CKPsiXcXe_a.N = nc;
    for (int ii = 0;ii<nc;ii++){
         for(int jj=0;jj<ne;jj++){
             CKPsiXcXe_a.val[b+counter*nc] = exp(-sum(Xc,Xe,thetaC,p,ii,jj));
             counter++; 
             if(counter == ne){
                counter =0;
                b++;
             }
         }
    }
    CKPsiXeXc_a = CKPsiXcXe_a.transpose();
    //Compute kriging parameters.
    
    // solve the rest of the kriging model
    //left it as an array since multi-dimensional may need an array; and it is more convenient
    Arr num_a = mu_num_den(UPsiXc_a,Yc_a,oneNc_a);
    //Arr tmp = UPsiXc_a/Yc_a;
    Arr den_a = mu_num_den(UPsiXc_a,oneNc_a,oneNc_a);
    muc_a = num_a%den_a;
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
    d_a.Init(d,ne,1);
    num_a = mu_num_den(UPsidXe_a,d_a,oneNe_a);
    den_a = mu_num_den(UPsidXe_a,oneNe_a,oneNe_a);
    mud_a.Init(1,1);
    mud_a = num_a%den_a;
    for(int ii=0;ii<nc;ii++){
        dif[ii] = Yc[ii]-muc_a.val[0]; 
    }
    //difference
    for(int ii=0;ii<ne;ii++){
        difd[ii] = d[ii]-mud_a.val[0]; 
    }
    Arr dif_a(dif,nc,1);
    num_a = mu_num_den(UPsiXc_a,dif_a,dif_a);
    SigmaSqrc_a.Init(1,1); 
    SigmaSqrc_a = num_a%nc;
    Arr difd_a(difd,ne,1);
    num_a = mu_num_den(UPsidXe_a,difd_a,difd_a);
    SigmaSqrd_a.Init(1,1); 
    SigmaSqrd_a = num_a%ne;

    //---------------------------------------------//
    //                construct C
    //---------------------------------------------//
    //  C is the culmination of the kriging model,
    //  since it is cokriging, C is formulated of the 
    //  cheap, expensive, and difference models and is the 
    //  main mathematical difference between kriging
    //  and cokriging
    //---------------------------------------------//
    C1_a = times(SigmaSqrc_a.val[0],CKPsiXc_a);//
    C2_a = times(rho*SigmaSqrc_a.val[0],CKPsiXcXe_a);//
    C3_a = times(rho*SigmaSqrc_a.val[0],CKPsiXeXc_a);//
    C4_a = times(rho*rho*SigmaSqrc_a.val[0],CKPsiXe_a)+times(SigmaSqrd_a.val[0],CKPsidXe_a);//
    for(int ii=0;ii<nc*nc;ii++){C1[ii]=SigmaSqrc_a.val[0]*CKPsiXc_a.val[ii];}
    for(int ii=0;ii<ne*nc;ii++){C2[ii]=rho*SigmaSqrc_a.val[0]*CKPsiXcXe_a.val[ii];}
    for(int ii=0;ii<ne*nc;ii++){C3[ii]=rho*SigmaSqrc_a.val[0]*CKPsiXeXc_a.val[ii];}
    for(int ii=0;ii<ne*ne;ii++){C4[ii]=rho*rho*SigmaSqrc_a.val[0]*CKPsiXe_a.val[ii]+SigmaSqrd_a.val[0]*CKPsidXe_a.val[ii];}
    for(int ii=0;ii<(nc+ne)*(nc+ne);ii++){ C[ii]=0; }//Initialize to 0
    Arr tmp = C2_a.transpose();
    Arr C_12 = concatinate(C1_a,C2_a,2);
    Arr C_34 = concatinate(C3_a,C4_a,2);
    Arr C_adbg =  concatinate(C_12,C_34,1);
    cout << "============= \ndebugging \n=============\n";
    C_adbg.print("C");
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
    cout << endl;
    cout << endl;
    Write1Darray(C,ne+nc,ne+nc);
    //---------------------------------------------//
    //         combine Ye and Yc into Y
    //---------------------------------------------//
    for(int ii = 0;ii< (ne+nc)*(ne+nc);ii++){UC[ii] = 0;}
    Cholesky(ne+nc,C,UC); 
    for(int ii = 0;ii< nc;ii++){Y[ii] = Yc[ii];}
    for(int ii = 0;ii< ne;ii++){Y[ii+nc] = Ye[ii];}
    cout<<"\n\n\n\n\n:Y: ";Write1Darray(Y,ne+nc,1);
    
    double *num = mu_num_den(UC,Y,ne+nc,oneNeNc);
    double *den = mu_num_den(UC,oneNeNc,ne+nc,oneNeNc);
    mu = num[0]/den[0];
    cout << "\nmu\n" << mu;
    delete [] num;
    delete [] den;
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
    double * Trans;
    Trans = transpose(UPsiX,n);
    double * MLD1;
    double * MLD2;
    double * MLD3;
    MLD1 = matrixLeftDivision(Trans,/*\*/Y /*size info*/ ,n,1);
    MLD2 = matrixLeftDivision(UPsiX,/*\*/ MLD1,n,1);
    MLD3 = matrixMultiply(oneN, /* X */MLD2 ,1,1,n);
    cout << "MLD3_old\n";
    Write1Darray(MLD3,1,1);
    delete [] Trans;
    delete [] MLD1;    
    delete [] MLD2;
    //---------------------------------------------//
    return MLD3;
}
//************************************************
Arr mu_num_den(Arr& UPsiX,Arr& Y,Arr& oneN){
    //---------------------------------------------//
    // Used to simplify some of the matrix math used in solving for kriging 
    // Inputs:
    //     UPsiX is a matrix of size nxn, 
    //     Y is a array of size nx1 or 1xn; since the array are column arrays the size doesn't matter for a 1 d vector. 
    //     oneN is also an array of size nx1 or 1xn
    // Returns:
    //    array of 1x1 solving the matrix equation One*(UPsiX\(UPsiX\Y))
    //Arr MLD1_a;
    //MLD1_a.Init(Y.M,Y.N);
    //Solve vector algebra
    Arr MLD2_a= (UPsiX/(UPsiX.transpose()/Y));
    Arr MLD1_a = oneN*MLD2_a;
    oneN.print("one_n");
    MLD2_a.print("MLD2_a");
    MLD1_a.print("MLD3_a");
    //MLD1 = matrixLeftDivision(Trans,/*\*/Y /*size info*/ ,n,1);
    //MLD2 = matrixLeftDivision(UPsiX,/*\*/ MLD1,n,1);
    //MLD3 = matrixMultiply(oneN, /* X */MLD2 ,1,1,n);
    //---------------------------------------------//
    return MLD1_a;
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
    double* tran_arr = new double[n*n];
    int counter  = 0;
    int b = 0;
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
void Cholesky_arr(const Arr& S,Arr& D){
    //---------------------------------------------//
    // Solve cholesky of a square array
    // Inputs: 
    //    S: input array
    //    D: output array. Cholesky of input array. 
    //---------------------------------------------//
    int d = S.M;
    D.M = S.M;D.N=S.N;
    for(int k=0;k<d;++k){
        double sum=0.;
        for(int p=0;p<k;++p)sum+=D.val[k*d+p]*D.val[k*d+p];
        D.val[k*d+k]=sqrt(S.val[k*d+k]-sum);
        for(int i=k+1;i<d;++i){
           double sum=0.;
           for(int p=0;p<k;++p)sum+=D.val[i*d+p]*D.val[k*d+p];
           D.val[i*d+k]=(S.val[i*d+k]-sum)/D.val[k*d+k];
        }
    }
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
void buildPsi(Arr& x, double* theta, Arr& CKPsixRtn ){
    //---------------------------------------------//
    // Inputs: 
    //    x: expected to be a one dimensional 
    //        square row array
    //    theta: only tested for one dimension, but 
    //        hopefully be changed in the future
    // Outputs: 
    //    CKPsixRTN: an square row based 1 d array 
    //         size (array.M*array.M)
    //---------------------------------------------//
    // NOTES:
    // solve for psi, not apart of the cokriging 
    // class due to information hiding, I don't want 
    // this function modifying any private variables.
    // However I should be able to make this more 
    // opject oriented, by converting the class to a
    // kriging class instead of a cokriging, that 
    // cokriging calls when it creates a cokriging model.
    // Since cokriging is mainly just a 
    // series of kriging functions
    //---------------------------------------------//
    int n = x.M;
    CKPsixRtn.M = x.M;
    CKPsixRtn.N = x.M;
    double PsiX[n][n];//initilize to zeros
    double CKPsiX[n][n];//initilize to zeros
    int EyeN[n][n] ;
    int counter=0;int b; //increment varables
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = 0; jj<n;jj++){
            PsiX[ii][jj] = 0;
            CKPsiX[ii][jj] = 0;
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
            PsiX[ii][jj] = exp(-sum(x.val,x.val,theta,p,ii,jj));
        }
    }
    float eps = 2.2204*pow(10,-16);
    counter = 0;b = 0;
    for(int ii = 0; ii<n;ii++){ 
        for(int jj = 0; jj <n;jj++){
            CKPsiX[ii][jj] = PsiX[ii][jj] + PsiX[jj][ii]+EyeN[ii][jj]+EyeN[ii][jj]*eps;
            CKPsixRtn.val[b+counter*n] = CKPsiX[ii][jj];
            counter++;
            if(counter == n){
               counter =0;
               b++;
            }
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
}
//************************************************
void cokriging::predictor(double* x,int n){
    //---------------------------------------------//
    // Evaluates the surface for a given independent
    // variable.
    // Inputs:
    //    x - array, currently 1-d plan to develop n-d
    //    n - size or array
    // Outputs:
    //    y - currently not returned, but will be the surrogate result
    //---------------------------------------------//
    double* cc ;
    double* cd  = new double[ne];for(int ii=0;ii<ne;ii++){cd[ii]=0;}
    double* cd1;
    double* cd2;
    double* c   = new double[nc+ne];for(int ii=0;ii<nc+ne;ii++){c[ii]=0;}
    double* diffY_mu = new double[ne+nc];for(int ii=0;ii<nc+ne;ii++){diffY_mu[ii]=0;}
    double* Yout;
    cc = c_pred(SigmaSqrc,rho,Xc,nc,x,n,thetaC);
    cd1 = c_pred(SigmaSqrc,rho*rho,Xe,ne,x,n,thetaC);
    cd2 = c_pred(SigmaSqrd,1,Xe,ne,x,n,thetaD);
    //combine cd1 and cd2
    for(int ii=0;ii<ne;ii++){cd[ii]=cd1[ii]+cd2[ii];}
    //concatinate cc and cd
    for(int ii=0;ii<nc;ii++){c[ii]=cc[ii];}
    for(int ii=0;ii<ne;ii++){c[ii+nc]=cd[ii];}
    //solve the surrogate
    for(int ii=0;ii<ne+nc;ii++){diffY_mu[ii]=Y[ii]-mu;}
    double * Trans;
    Trans =  transpose(UC,ne+nc);
    double * MLD1;
    MLD1 =  matrixLeftDivision(Trans,diffY_mu,ne+nc,1);
    double * MLD2;
    MLD2 = matrixLeftDivision(UC,MLD1 ,ne+nc,1);
    Yout =  matrixMultiply(c,MLD2,1,1,nc+ne);
    delete [] MLD1;
    delete [] MLD2;
    delete [] Trans;
    cout<<"\n:Yout: "<< Yout[0]+mu;
    delete [] cc;
    delete [] cd;
    delete [] c;
    delete [] cd1;
    delete [] cd2;
    delete [] diffY_mu;
    delete [] Yout;
    //delete (cd ); 
    //delete (cd1); 
    //delete (cd2); 
    //delete (c  ); 
    //delete (diffY_mu);
    //delete (UC_t);
}
//************************************************
double* c_pred(double sigma, double rho,double x1[], int n1, double x2[],int n2, double theta[]){
    int p =2;
    double* c = new double[n1];
    for(int ii=0;ii<n1;ii++){
        //get sum value, different from previous usage.
        c[ii]=rho*sigma*exp(-sum_pred(x1,x2,theta,p,ii,n2));
    }
    return c;
}
//************************************************
double sum_pred(double x1[],double x2[],double theta[],int p,int ii,int n){
    double sumVal = 0;
    for(int jj = 0;jj<n;jj++){
        sumVal+=theta[0]*pow(abs(x1[ii]-x2[jj]),p);
    }
    return sumVal;
}
void Print(struct Array arr){

    for(int ii = 0; ii < arr.size;ii++){
        cout << " " << arr.val[ii]; 
    }
    cout << endl;
}

