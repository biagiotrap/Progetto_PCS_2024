#include "Utils.hpp"
#include<iostream>
#include<fstream>
#include<sstream>
#include"Eigen/Eigen"
#include<vector>
using namespace std;
using namespace Eigen;
namespace fractureLibrary{
bool ImportData(const string &filepath,Fractures& fracture){
    ifstream file;
    file.open(filepath);
    if(file.fail()){
        return false;
    }
    string line;
    getline(file,line);
    getline(file,line);
    istringstream convertN(line);
    convertN>>fracture.FractureNumber;
    for(unsigned int i=0; i<fracture.FractureNumber;i++){
        getline(file,line);
        getline(file,line,';');
        istringstream convertId(line);
        unsigned int id;
        convertId>>id;
        fracture.Id.push_back(id);
        getline(file,line);
        istringstream convertNumV(line);
        unsigned int numV;
        convertNumV>>numV;
        fracture.VerticeNumber.push_back(numV);
        getline(file,line);
        vector<double> data;
        data.resize(fracture.VerticeNumber[i]*3);
        for(unsigned int c=0;c<fracture.VerticeNumber[i]*3;c++){
            if((c-3)%(fracture.VerticeNumber[i])==0){
                getline(file,line);
            }
            else{
                getline(file,line,';');
            }
            istringstream convertD(line);
            convertD>>data[c];
        }
        unsigned int p=0;
        for(unsigned int d=0;d<fracture.VerticeNumber[i];d++){
            Vector3d coord;
            for(unsigned int m=0;m<fracture.VerticeNumber[i]*3;m=m+4){
                coord(p)=data[m+d];
                p++;
            }
            fracture.Coordinates.push_back(coord);
            p=0;
            cout<<coord<<endl;
        }
    }
    return true;
}
}

