#ifndef __FRACTURES_H__
#define __FRACTURES_H__

#include "Eigen/Eigen"
#include <map>

using namespace std;
using namespace Eigen;

// Definiamo una struttura "fractures" dove importiamo le informazioni delle fratture all'interno dei file

namespace fractureLibrary{

struct Fractures
{
    unsigned int FractureNumber ;                     // intero positivo contiene il numero totale di fratture nel file
    vector<unsigned int> VerticeNumber ;             //vettore di interi positivi contiene per ogni frattura il numero di vertici
    vector<unsigned int> Id ;                       // Vettore di interi positivi contiene l'identificativo di ogni frattura
    vector<Vector3d> Coordinates ;                   // vettore 1xFractureNumber di vettori 3x1 di double. Contiene le coordinate dei vertici
    vector<vector<Vector3d>> Segments;
    vector<vector<unsigned int>> ListVertices;
    MatrixXd VerticesCoordinates;

};
struct Traces{
    unsigned int TracesNumber=0 ;                     // Numero di traccia trovate
    map<unsigned int, bool> CoupleIdTips ;    // mappa che contiene la coppia Id traccia e esito verifica della traccia(TIPS);
    vector<Vector3d> TracesPoints;                  //vettore di vettori che contengono le coordinate degli estremi delle tracce
    vector<Vector2d> IdFractures;
    map<unsigned int,vector<unsigned int>> CoupleFracturesTraces; //Id frattura a cui si associa un vettore contenente gli id delle tracce formate dalla frattura
    vector<double> LengthsTrace;
};
}

#endif
