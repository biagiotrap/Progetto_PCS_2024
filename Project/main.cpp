#include <Utils.hpp>
#include"Fractures.hpp"
#include<iostream>
using namespace std;
using namespace fractureLibrary;
int main(){
    Fractures fracture;
    string filepath = "DFN/FR3_data.txt";
    ImportData(filepath, fracture);
    return 0;
}
