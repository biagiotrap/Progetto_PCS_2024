#include <Utils.hpp>
#include"Fractures.hpp"
#include<iostream>
using namespace std;
using namespace fractureLibrary;
int main(){
    Fractures fracture;
    string filepath = "DFN/FR3_data.txt";
    string fileOutput="Tracce.txt";
    ImportData(filepath, fracture);
    DefineTraces(fileOutput, fracture);
    return 0;
}
