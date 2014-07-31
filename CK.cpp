#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include "cokriging.h"
#include "vector2array.h"
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
    double InitthetaD[] = {-1.1026};
    double InitthetaC[] = {1};
    double Initrho = 1.9454;
    int nc = InitXc.size();
    int ne = InitXe.size();
    double* Xe_a = vec2array(InitXe);
    double* Ye_a = vec2array(InitYe);
    double* Xc_a = vec2array(InitXc);
    double* Yc_a = vec2array(InitYc);
    cokriging CK(Xe_a,Ye_a,Xc_a,Yc_a,InitthetaD,InitthetaC,Initrho,nc,ne);
    cout <<"Begin building cokriging model" << endl;
    CK.buildModel();
//    CK.write(); 
    cout << endl;
} 

int main(){
    printf("Reading input.in\n");
    ReadData();
    return 0;
}
