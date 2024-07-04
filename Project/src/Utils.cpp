#include "Utils.hpp"
#include<iostream>
#include<fstream>
#include<sstream>
#include"Eigen/Eigen"
#include<vector>
#include<iomanip>



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
    getline(file,line); //Riga "#Number of Fractures" (Scarta)
    getline(file,line); //Numero di fratture
    istringstream convertN(line);
    convertN>>fracture.FractureNumber;
    for(unsigned int i=0; i<fracture.FractureNumber;i++){ //Itera sul numero di fratture
        getline(file,line); //Riga "#FractureId; NumVertices" (Scarta)
        getline(file,line,';'); //Legge fino al ';', prende id frattura
        istringstream convertId(line);
        unsigned int id;
        convertId>>id;
        fracture.Id.push_back(id); //Aggiunto al vettore
        getline(file,line); //Prende numero vertici i-esima frattura
        istringstream convertNumV(line);
        unsigned int numV;
        convertNumV>>numV;
        fracture.VerticeNumber.push_back(numV); //Aggiunto al vettore
        getline(file,line); //Riga "#Vertices" (Scarta)
        vector<double> data;
        data.resize(fracture.VerticeNumber[i]*3); //3 Coordinate per numero di vertici
        //Per memorizzare le coordinate, le inseriamo prima tutte in un unico vettore di lunghezza 3*(numero di vertici),
        //poi lo scomponiamo in tanti vettori di tre elementi quanti sono i numeri dei vertici.
        //Ogni vettore contiene le tre coordinate prese muovendosi nel vettore iniziale di
        //(numero di vertici) posizioni
        for(unsigned int c=0;c<fracture.VerticeNumber[i]*3;c++){ //Iteriamo su 3*(numero di vertici) coordinate
            if((c-(fracture.VerticeNumber[i]-1))%(fracture.VerticeNumber[i])==0){ //L'ultima coordinata della riga non termina con un ";
                getline(file,line);
            }
            else{
                getline(file,line,';'); //Le prime ((numero di vertici)-1) sì
            }
            istringstream convertD(line);
            convertD>>data[c]; //Aggiunge al vettore "data" la coordianta c-esima
        }
        //Le coordinate sono ora da dividere tra i vettori. Ogni terna ha indici (numV)*m + d
        //con m che varia tra 0 e 3 e d che varia tra 0 e (numV)
        unsigned int p=0;
        for(unsigned int d=0;d<fracture.VerticeNumber[i];d++){
            Vector3d coord;
            for(unsigned int m=0;m<fracture.VerticeNumber[i]*3;m=m+fracture.VerticeNumber[i]){
                coord(p)=data[m+d];
                p++;
            }
            fracture.Coordinates.push_back(coord);
            p=0;
        }

    }
    file.close();

    return true;

}

void ComputeSegments(const Fractures& fracture, vector<vector<Vector3d>>& segments) { // scelto di il vettore di vettori in modo da tenere separato i segmenti di ogni frattura

    for (unsigned int i = 0; i < fracture.FractureNumber; ++i) { // Iteriamo su tutte le fratture
        vector<Vector3d> fractureSegments; // Vettore dei segmenti per la frattura corrente

        // Calcola startIndex e endIndex per la frattura corrente
        unsigned int startIndex = 0;
        for (unsigned int j = 0; j < i; ++j) { // Somma il numero di vertici delle fratture precedenti per ottenere l'indice di inizio
            startIndex += fracture.VerticeNumber[j];
        }
        unsigned int endIndex = startIndex + fracture.VerticeNumber[i]; // Calcola l'indice di fine

        // Calcola i segmenti per la frattura corrente
        for (unsigned int j = startIndex; j < endIndex; ++j) { // Itera su tutti i vertici della frattura corrente
            Vector3d segment = fracture.Coordinates[(j + 1) % endIndex] - fracture.Coordinates[j]; // Calcola il segmento come differenza tra i punti dei vertici

            // Aggiunge il segmento al vettore dei segmenti
            fractureSegments.push_back(segment);
        }

        // Aggiunge i segmenti della frattura corrente al vettore principale dei segmenti
        segments.push_back(fractureSegments);

        // Stampa i segmenti calcolati per la frattura corrente
        cout << "Frattura " << i + 1 << ":" << endl;
        for (unsigned int j = 0; j < fractureSegments.size(); ++j) { // Stampa tutti i segmenti della frattura corrente
            cout << "    Segmento " << j + 1 << ": (" << fractureSegments[j].transpose() << ")" << endl;
        }
    }
}






bool DefineTraces(const string &fileOutput, Fractures& fracture, Traces& trace, vector<vector<Vector3d>>& segments){

    double tol=1e-10;
    vector<Vector3d> n; //Contiene terne che sono i versori normali alle fratture
    vector<double> d;
    vector<Vector3d> Points;
    n.resize(fracture.FractureNumber); //Dimensione di n=1*(numero di fratture)
    d.resize(fracture.FractureNumber);//Dimensione di d=1*(numero di fratture)
    Points.reserve(fracture.FractureNumber*(fracture.FractureNumber-1)/2);
    int sommaParziale=0; //Le coordinate non sono divise per fratture. Usiamo quindi "Somma Parziale" per considerare solo i primi tre vertici di ogni frattura
    for(unsigned int i=0;i<fracture.FractureNumber;i++){
        //Calcolo generatori del piano in cui è situata la frattura
        Vector3d u;
        Vector3d v;
        u=fracture.Coordinates[sommaParziale+2]-fracture.Coordinates[sommaParziale];
        v=fracture.Coordinates[sommaParziale+1]-fracture.Coordinates[sommaParziale];
            //Passo alla prossima frattura
        cout<<setprecision(10);
        sommaParziale+=fracture.VerticeNumber[i];
        n[i]=(u.cross(v))/(u.norm()*v.norm()); //Versore perpendicolare a u e v
        d[i]=n[i].dot(fracture.Coordinates[sommaParziale]);
    }
    for(unsigned int i=0; i<n.size()-1;i++){
        for(unsigned int j=1; j<n.size();j++){
            //calcolo del sistema lineare
            j+=i;
            Vector3d t; //Vettore tangente alla traccia
            t=n[i].cross(n[j]);
            MatrixXd A=MatrixXd::Ones(3, 3);
            A<<n[i], n[j], t;
            cout<<A.determinant()<<endl;
            if(A.determinant()>tol){ //Se il determinante della matrice A contenente


                //le normali dei piani e il vettore perpenidcolare ad esse è diverso da zero vuol dire che i piani si intersecano
                Vector3d b;
                b<<d[i],d[i+j],0;
                Vector3d x=A.lu().solve(b);
                //Una volta controllato se i due piani si intersecano bisogna
                //controllare le intersezioni tra t e i segmenti delle fratture

                vector<Vector3d> intersT;
                unsigned int a=0;
                for(unsigned int s=0;s<fracture.VerticeNumber[i];i++){
                    Vector3d w=x-fracture.Coordinates[sommaParziale];
                    MatrixXd A3=MatrixXd::Ones(3, 3);
                    A3<<w, segments[i][s], t;
                    if(A3.determinant()<tol){  //verifica che le due retti siano incidenti
                        MatrixXd A2=MatrixXd::Ones(2, 3);
                        A2<<fracture.Coordinates[sommaParziale], t;
                        Vector2d coeff=(A2.transpose()).lu().solve(w);
                        if(coeff[0]>tol && coeff[0]<1-tol && coeff[1]>tol && coeff[1]<1-tol){ //verifica per vedere se il punto appartiene al segmento della frattura
                            Vector3d P=fracture.Coordinates[sommaParziale]+(coeff[0]*(segments[i][s]));
                            intersT.push_back(P);
                            a+=1;
                        }
                    }
                    if(a==2){
                        break;
                    }

                }
                for( unsigned int v=0;v<fracture.VerticeNumber[v+1];v++){
                    Vector3d w=x-fracture.Coordinates[sommaParziale];
                    MatrixXd A3=MatrixXd::Ones(3, 3);
                    A3<<w, segments[i][v], t;
                    if(A3.determinant()<tol){  //verifica che le due retti siano incidenti
                        MatrixXd A2=MatrixXd::Ones(2, 3);
                        A2<<fracture.Coordinates[sommaParziale], t;
                        Vector2d coeff=(A2.transpose()).lu().solve(w);
                        if(coeff[0]>tol && coeff[0]<1-tol && coeff[1]>tol && coeff[1]<1-tol){ //verifica per vedere se il punto appartiene al segmento della frattura
                            Vector3d P=fracture.Coordinates[sommaParziale]+(coeff[0]*(segments[i][v]));
                            intersT.push_back(P);
                            a+=1;
                        }
                    }
                    if(a==4){
                        break;
                    }
                }

            }
        }
    }
    ofstream outFile(fileOutput);
    outFile << "# Number of traces" << endl;

    outFile << trace.TracesNumber << endl;

    outFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for (unsigned int t=0;t<trace.TracesNumber;t++){
        //itera sulle tracce ecc
    }

    outFile << "# FractureId; NumTraces" << endl;
    for (unsigned int f=0;f<fracture.FractureNumber;f++){

    outFile << trace.TracesNumber << endl;

    outFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for (unsigned int t=0;t<trace.TracesNumber;t++){
        //itera sulle tracce ecc
    }

    outFile << "# FractureId; NumTraces" << endl;
    for (f=0;f<fracture.FractureNumber;f++){

        //Per ogni frattura il numero di tracce
    }
    return true;
}
// bool ComputeLenghts(Fractures& fracture){
//     outfile << "# TraceId; Tips; Length " << endl;
//     for (t=0;t<fracture.TracesNumber;t++){
//         //values
//     }
//     return true;
// }
}
}

