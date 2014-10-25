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
    // Only require to delete pointers, 
    // Arr class handles it's own pointer destruction
    //    and does not need to be implicitly called
    delete [] thetaD;//
    delete [] thetaC;//
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
    Xe_a.Init(Initxe,ne,1);
    Ye_a.Init(Initye,ne,1);

    // Lower fidelity - cheap parameters
    nc = Initnc;
    Xc_a.Init(Initxc,nc,1);
    Yc_a.Init(Inityc,nc,1);

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
    UPsiXc_a   .Init(0.0,nc,nc);
    CKPsiXcXe_a.Init(0.0,nc,ne);
    CKPsiXeXc_a.Init(0.0,ne,nc);
    UPsidXe_a  .Init(0.0,ne,ne); 
    UPsiXe_a   .Init(0.0,ne,ne);
    CKPsiXc_a  .Init(0.0,nc,nc);
    CKPsiXe_a  .Init(0.0,ne,ne);
    CKPsidXe_a .Init(0.0,ne,ne);
    UC_a       .Init(0.0,ne+nc,ne+nc);
    Y_a        .Init(0.0,ne+nc,1);
}
//************************************************
void cokriging::buildModel(){
    //---------------------------------------------//
    //Main function for developing the cokriging model
    //---------------------------------------------//

    //resize all arrays to fit the data set
    resize();
    int p = 2;//Currently a constant, but could be varied to change kriging differentiation
    //---------------------------------------------//
    // initialize  and declare local variables
    //---------------------------------------------//
    Arr oneNe_a(1.0,ne,1);
    Arr oneNc_a(1.0,nc,1);
    Arr oneNeNc_a(1.0,ne+nc,1);
    // used to contruct C
    Arr C1_a(0.0,nc,nc);//quadrant 1
    Arr C2_a(0.0,nc,ne);//quadrant 2
    Arr C3_a(0.0,ne,nc);//quadrant 3
    Arr C4_a(0.0,ne,ne);//quadrant 4
    Arr C_a(ne+nc,ne+nc);//just create the space without initializing
    //counters
    int counter = 0;
    int b = 0;
    //---------------------------------------------//
    // Build all the various Psi 
    // variables for cheap and expensive models
    //---------------------------------------------//
    buildPsi(Xc_a,thetaC,/*&*/CKPsiXc_a);
    buildPsi(Xe_a,thetaC,/*&*/CKPsiXe_a);
    buildPsi(Xe_a,thetaD,/*&*/CKPsidXe_a);
    
    //calculate cholesky, 
    UPsiXe_a  = CKPsiXe_a.cholesky(); 
    UPsiXc_a  = CKPsiXc_a.cholesky(); 
    UPsidXe_a = CKPsidXe_a.cholesky(); 
    //fill yet another Psi variable but not square and different results
    CKPsiXcXe_a.M = ne;
    CKPsiXcXe_a.N = nc;
    for (int ii = 0;ii<nc;ii++){
         for(int jj=0;jj<ne;jj++){
             CKPsiXcXe_a.val[b+counter*nc] = exp(-sum(Xc_a.val,Xe_a.val,thetaC,p,ii,jj));
             counter++; 
             if(counter == ne){
                counter =0;
                b++;
             }
         }
    }
    CKPsiXeXc_a = CKPsiXcXe_a.transpose();
    //Compute kriging parameters.
    
    //left following as an array since multi-dimensional may need an array; and it is more convenient
    Arr num_a = mu_num_den(UPsiXc_a,Yc_a,oneNc_a);
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
    d_a.Init(0.0,ne,1);
    for(int ii = 0;ii<ne;ii++){
        d_a.push(Ye_a.element(ii,0)-rho*Yc_a.element(nc-ne+ii,0),ii,0);
    }  
    //cheap
    num_a = mu_num_den(UPsidXe_a,d_a,oneNe_a);
    den_a = mu_num_den(UPsidXe_a,oneNe_a,oneNe_a);
    mud_a.Init(1,1);
    mud_a = num_a%den_a;
    //difference
    Arr dif_a = Yc_a+(-muc_a.element(0,0));
    num_a = mu_num_den(UPsiXc_a,dif_a,dif_a);
    SigmaSqrc_a.Init(1,1); 
    SigmaSqrc_a = num_a%nc;
    Arr difd_a = d_a+(-mud_a.element(0,0));
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
    
    Arr C_12 =    concatinate(C1_a,C2_a,2);
    Arr C_34 =    concatinate(C3_a,C4_a,2);
    Arr C_adbg =  concatinate(C_12,C_34,1);
   
    //---------------------------------------------//
    //    combine Ye and Yc into Y solve for mu
    //---------------------------------------------//
    UC_a  = C_adbg.cholesky();
    Y_a   = concatinate(Yc_a,Ye_a,2);
    num_a = mu_num_den(UC_a,Y_a,oneNeNc_a);
    den_a = mu_num_den(UC_a,oneNeNc_a,oneNeNc_a);
    mu_a  = num_a%den_a;
    
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
Arr mu_num_den(Arr& UPsiX,Arr& Y,Arr& oneN){
    //---------------------------------------------//
    // Used to simplify some of the matrix math used in solving for kriging 
    // Inputs:
    //     UPsiX is a matrix of size nxn, 
    //     Y is a array of size nx1 or 1xn; since the array are column arrays the size doesn't matter for a 1 d vector. 
    //     oneN is also an array of size nx1 or 1xn
    // Returns:
    //    array of 1x1 solving the matrix equation One*(UPsiX\(UPsiX\Y))
    //---------------------------------------------//
    
    //Solve vector algebra
    Arr MLD1_a = oneN*(UPsiX/(UPsiX.transpose()/Y));
    //---------------------------------------------//
    return MLD1_a;
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
    // Zero initialization
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
    Arr x_a(x,n,1);
    Arr cc_a =  c_pred(SigmaSqrc_a.val[0],rho,Xc_a,x_a,thetaC);
    Arr cd1_a = c_pred(SigmaSqrc_a.val[0],rho*rho,Xe_a, x_a,thetaC);
    Arr cd2_a = c_pred(SigmaSqrd_a.val[0],1,Xe_a, x_a,thetaD);
    //combine cd1 and cd2
    Arr cd_a = cd1_a+cd2_a;
    //concatinate cc and cd
    Arr c_a = concatinate(cc_a,cd_a,0);
    //solve the surrogate
    Arr diffy_mu_a = Y_a + (-mu_a.element(0,0));
     
    //solve the surrogate
    Arr Yout_a = c_a*(UC_a/(UC_a.transpose()/diffy_mu_a));
    
    //output to screen
    cout<<"Xin " <<x_a.element(0,0) << " Yout: "<< Yout_a.element(0,0)<<endl;
}
//************************************************
Arr c_pred(double sigma, double rho,Arr x1, Arr x2, double theta[]){
    int p =2;
    int n1 = x1.M;
    Arr c(0.0,n1,1);
    for(int ii=0;ii<n1;ii++){
        c.push(rho*sigma*exp(-sum_pred(x1,x2,theta,p,ii)),ii,0);
    }
    return c;
}
//************************************************
double sum_pred(Arr x1,Arr x2,double theta[],int p,int ii){
    double sumVal = 0;
    for(int jj = 0;jj<x2.M;jj++){
        sumVal+=theta[0]*pow(abs(x1.element(ii,0)-x2.element(jj,0)),p);
    }
    return sumVal;
}

