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
    //apertura del file per importare i dati
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
        //Per memorizzare le coordinate, le inseriamo prima tutte in un unico vettore
        //e poi lo scomponiamo in tanti vettori quanti sono i numeri dei vertici,
        //in cui ogni vettore contiene 3 elementi che sono stati presi muovendosi nel vettore iniziale di
        //"numero di vertici" posizioni
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
        }
    }
    return true;
}

bool DefineTraces(const string &fileOutput, Fractures& fracture){
    vector<Vector3d> n;
    vector<double> d;
    Vector3d t;
    vector<Vector3d> Points;
    n.resize(fracture.FractureNumber);
    d.resize(fracture.FractureNumber);
    Points.reserve(fracture.FractureNumber*(fracture.FractureNumber-1)/2);
    int sommaParziale=0;
    for(unsigned int i=0;i<fracture.FractureNumber;i++){
        //Calcolo generatoti del piano in cui è situata la frattura
        Vector3d u;
        Vector3d v;
        u=fracture.Coordinates[sommaParziale+1]-fracture.Coordinates[sommaParziale];
        v=fracture.Coordinates[sommaParziale+2]-fracture.Coordinates[sommaParziale];
        sommaParziale+=fracture.VerticeNumber[i];
        n[i]=(u.cross(v))/u.norm()*v.norm(); //vettore perpendicolare a u e v
        d[i]=n[i].dot(fracture.Coordinates[sommaParziale]);
    }
    for(unsigned int i=0; i<n.size()-1;i++){
        for(unsigned int j=1; j<n.size();j++){
            //calcolo del sistema lineare
            j+=i;
            Vector3d t;
            t=n[i].cross(n[j]);
            MatrixXd A=MatrixXd::Ones(3, 3);
            A<<n[i], n[i+j], t;
            if(A.determinant()!=0){ //Se il determinante della matrice A contenente
  //le normali dei piani e il vettore perpenidcolare ad esse è diverso da zero vuol dire che i piani si intersecano
                Vector3d u;
                Vector3d b;
                b<<d[i],d[i+j],0;
                Vector3d x=A.lu().solve(b);
                //Una volta controllato se i due piani si intersecano bisogna
 //controllare se le fratture presenti nei due piani hanno due o meno intersezioni
 //2 intersezioni:traccia passante
 //altrimenti:traccia interna



            }
        }
    }
    return true;
}


}


