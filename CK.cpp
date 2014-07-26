#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include "cokriging.h"
using namespace std;
//Reading method is temporary and will liekly be changed once 2d tests are run
//Split string used in the other split
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

//another split method - main split
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

//removed the bounding [] in the files
vector<string> FilterData(const string &s){
    vector<string> Xc = split(s,',');
    vector<string> ctmp;
    int n = Xc.size()-1;
    ctmp = split(Xc[0],'[');  Xc[0] = ctmp[1];
    ctmp = split(Xc[n],']');  Xc[n] = ctmp[0];
    return Xc;
}
vector<double> Str2Double(vector <string> StringVec){
     vector<double> DoubleVec(StringVec.size());
     for (int ii=0; ii<StringVec.size();ii++){
          DoubleVec[ii]=atof(StringVec[ii].c_str());
     }
     return DoubleVec;
}
//reads the input file
void ReadData(){
    ifstream inData;
    inData.open("input.in");
    string Lines;
    getline(inData,Lines);
    vector<string> Xe = FilterData(Lines);
    getline(inData,Lines);
    vector<string> Xc = FilterData(Lines);
    getline(inData,Lines);
    vector<string> Ye = FilterData(Lines);
    getline(inData,Lines);
    vector<string> Yc = FilterData(Lines);
    vector<double> InitXe = Str2Double(Xe);
    vector<double> InitYe = Str2Double(Ye);
    vector<double> InitXc = Str2Double(Xc);
    vector<double> InitYc = Str2Double(Yc);
    //Fill arrays need to switch with a for loop and some case statements
    //cout << "\nXc: " ;for(int ii=0; ii<InitXc.size();ii++){ cout << InitXc[ii]<<" ";}
    //cout << "\nXe: " ;for(int ii=0; ii<InitXe.size();ii++){ cout << InitXe[ii]<<" ";}
    //cout << "\nYc: " ;for(int ii=0; ii<InitYc.size();ii++){ cout << InitYc[ii]<<" ";}
    //cout << "\nYe: " ;for(int ii=0; ii<InitYe.size();ii++){ cout << InitYe[ii]<<" " ;}
    vector<double> InitthetaD(1,-1.1026);
    vector<double> InitthetaC(1,1);
    double Initrho = 1.9454;
    cokriging CK(InitXe,InitYe,InitXc,InitYc,InitthetaD,InitthetaC,Initrho);
    CK.buildModel();
//    CK.write(); 
    cout << endl;
} 

int main(){
    printf("Reading input.in\n");
    ReadData();
    return 0;
}
