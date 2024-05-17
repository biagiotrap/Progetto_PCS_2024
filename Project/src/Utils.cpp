#include "Utils.hpp"
#include<iostream>
#include<fstream>
#include<sstream>
#include"Eigen/Eigen"
#include<vector>
using namespace std;
using namespace Eigen;
namespace FractureLibrary{
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
    convertN>>FractureNumber;
    for(unsigned int i=0; i<FractureNumber;i++){
        getline(file,line);
        getline(file,line,';');
        istringstream convertId(line);
        convertId>>Id.push_back(line);
        getline(file,line);
        istringstream convertNumV(line);
        convertNumV>>NumberVertices.push_back(line);
        getline(file,line);
        vector<unsigned int> data;
        data.reserve(NumberVertices*3-1);
        for(unsigned int c=0;c<NumberVertices*3-1;c++){
            istringstream convertD(line);
            if((c-3)%(NumberVertices)==0){
                getline(file,line);
                convertD>>data[c];
            }
            getline(file,line,';');
            convertD>>data[c];
        }
        unsigned int p=0;
        for(unsigned int d=0;d<NumberVertices;d++){
            Vector3d coord;
            for(unsigned int c=0;c<NumberVertices*3;c=c+4){
                c+=d;
                coord(p)=data[c];
                p++;
            }
            Coordinates.push_back(coord);
            p=0;
        }
    }


}

