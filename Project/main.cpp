#include <Utils.hpp>
#include"Fractures.hpp"
#include<iostream>
using namespace std;
using namespace fractureLibrary;
int main(){

    Fractures fracture;
    string filepath = "DFN/FR3_data.txt";
    string fileOutput="./Tracce.txt";
    if (ImportData(filepath, fracture)) {
        cout << "Dati importati correttamente da " << filepath << endl;

        vector<vector<Vector3d>> segments;
        ComputeSegments(fracture, segments);

    } else {
        cout << "Errore nell'importazione dei dati." << endl;
    }
    string fileOutput="Tracce.txt";
    ImportData(filepath, fracture);
    DefineTraces(fileOutput, fracture);

    return 0;
}
