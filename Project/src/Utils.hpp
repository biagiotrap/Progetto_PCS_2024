#ifndef __UTILS__H_
#define __UTILS_H_

#include "Fractures.hpp"

using namespace std;

namespace fractureLibrary{

// restituisce esito lettura dati, Vero se con successo, Falso altrimenti
bool ImportData(const string &filePath, Fractures& fracture);


// restituisce esito calcolo Traccia, e scrive risultato in fileOutput
bool DefineTraces(const string &fileOutput, Fractures& fracture);

// restituisce esito verifica (passante/non passante) mediante Tips e scrive risultato in file fileOutput
bool CheckTrace(const string &fileOutput, Fractures& fracture) ;


// restituisce un vettore che contiene le lughezze delle tracce
bool ComputeLenghts(Fractures& fracture);

// Ordinare usando l'algoritmo MergeSort--> da vedere quali parametri metterci dentro
bool MergeSort();

void ComputeSegments(const Fractures& fracture, vector<vector<Vector3d>>& segments);




}







#endif
