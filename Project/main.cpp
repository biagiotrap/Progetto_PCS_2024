#include <Utils.hpp>
#include"Fractures.hpp"
#include<iostream>
#include "UCDUtilities.hpp"
using namespace std;
using namespace fractureLibrary;
int main(){

    Fractures fracture;
    Traces trace;
    string filepath = "DFN/FR82_data.txt";
    string fileOutput="./Tracce.txt";
    vector<vector<unsigned int>> triangles;
    VectorXi materials;
    Gedim::UCDUtilities exporter;
    if (ImportData(filepath, fracture)) {
        cout << "Dati importati correttamente da " << filepath << endl;


        ComputeSegments(fracture);
        DefineTraces(fileOutput, fracture, trace);
        GedimInterface(fracture, triangles, materials);
        exporter.ExportPolygons("./polygons_82.inp",fracture.VerticesCoordinates , triangles,{},{}, materials);

    } else {
        cout << "Errore nell'importazione dei dati." << endl;
    }

    return 0;
}
