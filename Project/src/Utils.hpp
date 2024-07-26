#ifndef __UTILS__H_
#define __UTILS_H_

#include "Fractures.hpp"

using namespace std;

namespace fractureLibrary{

// restituisce esito lettura dati, Vero se con successo, Falso altrimenti
bool ImportData(const string &filePath, Fractures& fracture);


// restituisce esito calcolo Traccia, e scrive risultato in fileOutput
bool DefineTraces(const string &fileOutput, Fractures& fracture, Traces& trace);


// restituisce la lunghezza di un segmento
double ComputeLengths(Vector3d& a, Vector3d& b);

// Ordinare usando l'algoritmo MergeSort--> da vedere quali parametri metterci dentro
void MergeSort(Traces& trace);

void Merge(vector<double>& lengths, const unsigned int& sx, const unsigned int cx, const unsigned int& dx);

void ComputeSegments(Fractures& fracture);

void Sorting(vector<double>& vec);

void GedimInterface(Fractures& fracture, vector<vector<unsigned int>>& triangles, VectorXi& materials);

}







#endif
