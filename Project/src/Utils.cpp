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
    vector<unsigned int> sumP;
    sumP.push_back(0);
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
        unsigned int q=0;
        vector<unsigned int> v1;
        v1.reserve(fracture.VerticeNumber[i]);
        for(unsigned int w=0; w<fracture.VerticeNumber[i];w++){
            q=sumP[i]+w;
            v1.push_back(q);
        }
        sumP.push_back(q+1);
        fracture.ListVertices.push_back(v1);
        getline(file,line); //Riga "#Verti;ces" (Scarta)
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

void ComputeSegments(Fractures& fracture) { // scelto di il vettore di vettori in modo da tenere separato i segmenti di ogni frattura
    vector<unsigned int> startIndex;
    Vector3d segment;
    startIndex.push_back(0);
    unsigned int numb=0;
    for (unsigned int i = 0; i < fracture.FractureNumber; ++i) { // Iteriamo su tutte le fratture
        vector<Vector3d> fractureSegments;
        // Calcola startIndex  per la frattura corrente
        numb += fracture.VerticeNumber[i];
        startIndex.push_back(numb);
        for (unsigned int j = 0; j < fracture.VerticeNumber[i]-1; ++j) {
            segment=fracture.Coordinates[startIndex[i]+j+1]- fracture.Coordinates[startIndex[i]+j];
            fractureSegments.push_back(segment);
        }
        fractureSegments.push_back(fracture.Coordinates[startIndex[i]+fracture.VerticeNumber[i]-1]- fracture.Coordinates[startIndex[i]]);
        fracture.Segments.push_back(fractureSegments);



        // Stampa i segmenti calcolati per la frattura corrente
        //cout << "Frattura " << i + 1 << ":" << endl;
        //for (unsigned int j = 0; j < fractureSegments.size(); ++j) { // Stampa tutti i segmenti della frattura corrente
            //cout << "    Segmento " << j + 1 << ": (" << fractureSegments[j].transpose() << ")" << endl;
        //}
    }
}


double ComputeLengths(Vector3d& a, Vector3d& b){
    Vector3d c=a-b;
    return c.norm();
}

bool DefineTraces(const string &fileOutput, Fractures& fracture, Traces& trace){
    double tol=1e-10;
    vector<Vector3d> n; //Contiene terne che sono i versori normali alle fratture
    vector<double> d;
    vector<Vector3d> Points;
    n.resize(fracture.FractureNumber); //Dimensione di n=1*(numero di fratture)
    d.resize(fracture.FractureNumber);//Dimensione di d=1*(numero di fratture)
    Points.reserve(fracture.FractureNumber*(fracture.FractureNumber-1)/2);
    int sommaParziale=0; //Le coordinate non sono divise per fratture. Usiamo quindi "Somma Parziale" per considerare solo i primi tre vertici di ogni frattura
    vector<unsigned int> sommeParziali;
    sommeParziali.push_back(0);
    for(unsigned int i=0;i<fracture.FractureNumber;i++){
        //Calcolo generatori del piano in cui è situata la frattura
        Vector3d u;
        Vector3d v;
        u=fracture.Segments[i][fracture.VerticeNumber[i]-1];
        v=fracture.Segments[i][0];
        //Passo alla prossima frattura
        n[i]=(u.cross(v))/(u.norm()*v.norm());        //Versore perpendicolare a u e v
        d[i]=n[i].dot(fracture.Coordinates[sommaParziale]);
        sommaParziale+=fracture.VerticeNumber[i];
        sommeParziali.push_back(sommaParziale);
    }
    for(unsigned int i=0; i<n.size()-1;i++){
        for(unsigned int j=1; j<n.size();j++){
            //calcolo del sistema lineare
            j+=i;
            Vector3d t; //Vettore tangente alla traccia
            t=n[i].cross(n[j]);
            MatrixXd A=MatrixXd::Ones(3, 3);
            A<<n[i], n[j], t;
            MatrixXd A1=A.transpose();
            if(A1.determinant()>tol || A1.determinant()<-tol){ //Se il determinante della matrice A contenente                                            //le normali dei piani e il vettore perpenidcolare ad esse è diverso da zero vuol dire che i piani si intersecano
                Vector3d b;
                b<<d[i],d[j],0;
                Vector3d x=A1.lu().solve(b);
                Vector3d r=x+3*t;
                //Una volta controllato se i due piani si intersecano bisogna
                //controllare le intersezioni tra t e i segmenti delle fratture
                vector<Vector3d> intersT;
                unsigned int a=0;
                unsigned int c=0;
                for(unsigned int s=0;s<fracture.VerticeNumber[i];s++){
                    Vector3d w;
                    if(c==fracture.VerticeNumber[i]-1){
                        c=0;
                    }
                    w=x-fracture.Coordinates[sommeParziali[i]+c];
                    MatrixXd A3=MatrixXd::Ones(3, 3);
                    A3<<w, fracture.Segments[i][s], r-x;
                    Vector3d prodVett=fracture.Segments[i][s].cross(r-x);
                    MatrixXd A3t=A3.transpose();
                    if(A3t.determinant()<tol && A3t.determinant()>-tol &&
                        (prodVett.norm()<-tol || prodVett.norm()>tol)){  //verifica che le due rette siano incidenti
                        MatrixXd A2=MatrixXd::Ones(3, 2);
                        A2<<fracture.Segments[i][s], r-x;
                        Vector2d coeff=A2.fullPivLu().solve(w);
                        if(coeff[0]>tol && coeff[0]<1+tol){ //verifica per vedere se il punto appartiene al segmento della frattura
                            Vector3d P=fracture.Coordinates[sommeParziali[i]+c]+(coeff[0]*(fracture.Segments[i][s]));
                            cout<<i<<endl;
                            intersT.push_back(P);
                            a+=1;
                        }
                    }
                    if(a==2){
                        break;
                    }
                    c++;
                }
                c=0;
                for( unsigned int v=0;v<fracture.VerticeNumber[j];v++){
                    Vector3d w;
                    if(c==fracture.VerticeNumber[j]-1){
                        c=0;
                    }
                    w=x-fracture.Coordinates[sommeParziali[j]+c];
                    MatrixXd A3=MatrixXd::Ones(3, 3);
                    A3<<w, fracture.Segments[j][v], r-x;
                    MatrixXd A3t=A3.transpose();
                    Vector3d prodVett=fracture.Segments[j][v].cross(r-x);
                    if(A3t.determinant()<tol && A3t.determinant()>-tol &&
                        (prodVett.norm()<-tol || prodVett.norm()>tol)){  //verifica che le due rette siano incidenti
                        MatrixXd A2=MatrixXd::Ones(3, 2);
                        A2<<fracture.Segments[j][v], r-x;
                        Vector2d coeff=A2.fullPivLu().solve(w);
                        if(coeff[0]>tol && coeff[0]<1+tol){ //verifica per vedere se il punto appartiene al segmento della frattura
                            Vector3d P=fracture.Coordinates[sommeParziali[j]+c]+(coeff[0]*(fracture.Segments[j][v]));
                            cout<<j<<endl;
                            intersT.push_back(P);
                            a+=1;
                        }
                    }
                    if(a==4){
                        break;
                    }
                    c++;
                }
                if(size(intersT)==4){
                    //le due fratture formano una potenziale traccia
                    // Verifica di tutte le casistiche
                    double l01=ComputeLengths(intersT[0],intersT[1]);
                    double l02=ComputeLengths(intersT[0],intersT[2]);
                    double l03=ComputeLengths(intersT[0],intersT[3]);
                    double l12=ComputeLengths(intersT[1],intersT[2]);
                    double l13=ComputeLengths(intersT[1],intersT[3]);
                    double l23=ComputeLengths(intersT[2],intersT[3]);
                    if(l01<=tol || l23<=tol){
                        continue;
                    }
                    else{
                        Vector2d Id= Vector2d::Ones();
                        Id<< i, j;
                        trace.IdFractures.push_back(Id);    //vettore contenente gli id delle fratture che coustituiscono una potenziale traccia
                        if(l02<=tol){
                            if (l13<=tol){  //caso passante
                                trace.TracesPoints.push_back(intersT[0]);
                                trace.TracesPoints.push_back(intersT[1]); 
                                trace.CoupleIdTips[trace.TracesNumber]= false;
                                trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                trace.LengthsTrace.push_back(l01);
                                trace.TracesNumber=trace.TracesNumber+1;

                            }
                            else{
                                Vector3d u=intersT[3]-intersT[1];
                                Vector3d v1=intersT[2]-intersT[1];
                                Vector3d v2=intersT[2]-intersT[3];
                                double P1=u.dot(v1);
                                double P2=u.dot(v2);
                                if(P1*P2>tol){
                                    if(l01>l03 ){
                                        trace.TracesPoints.push_back(intersT[0]);
                                        trace.TracesPoints.push_back(intersT[3]);
                                        trace.CoupleIdTips[trace.TracesNumber]= true;
                                        trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                        trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                        trace.LengthsTrace.push_back(l03);
                                        trace.TracesNumber=trace.TracesNumber+1;

                                    }
                                    else{
                                        trace.TracesPoints.push_back(intersT[0]);
                                        trace.TracesPoints.push_back(intersT[1]);
                                        trace.CoupleIdTips[trace.TracesNumber]= true;
                                        trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                        trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                        trace.LengthsTrace.push_back(l01);
                                        trace.TracesNumber=trace.TracesNumber+1;
                                    }
                                }
                            }
                        }
                        if(l03<=tol){
                            if(l12<=tol){  //caso passante
                                trace.TracesPoints.push_back(intersT[0]);
                                trace.TracesPoints.push_back(intersT[1]);
                                trace.CoupleIdTips[trace.TracesNumber]= false;
                                trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                trace.LengthsTrace.push_back(l01);
                                trace.TracesNumber=trace.TracesNumber+1;

                            }
                            else{
                                Vector3d u=intersT[2]-intersT[1];
                                Vector3d v1=intersT[0]-intersT[2];
                                Vector3d v2=intersT[0]-intersT[1];
                                double P1=u.dot(v1);
                                double P2=u.dot(v2);
                                if(P1*P2>tol){
                                    if(l01>l02){
                                        trace.TracesPoints.push_back(intersT[0]);
                                        trace.TracesPoints.push_back(intersT[2]);
                                        trace.CoupleIdTips[trace.TracesNumber]= true;
                                        trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                        trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                        trace.LengthsTrace.push_back(l02);
                                        trace.TracesNumber=trace.TracesNumber+1;


                                    }
                                    else{
                                        trace.TracesPoints.push_back(intersT[0]);
                                        trace.TracesPoints.push_back(intersT[1]);
                                        trace.CoupleIdTips[trace.TracesNumber]= true;
                                        trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                        trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                        trace.LengthsTrace.push_back(l01);
                                        trace.TracesNumber=trace.TracesNumber+1;
                                    }
                                }
                            }
                        }
                        if(l12<=tol && l03>tol){
                            Vector3d u=intersT[3]-intersT[0];
                            Vector3d v1=intersT[3]-intersT[2];
                            Vector3d v2=intersT[0]-intersT[2];
                            double P1=u.dot(v1);
                            double P2=u.dot(v2);
                            if(P1*P2>tol){
                                if(l01>l13){
                                    trace.TracesPoints.push_back(intersT[1]);
                                    trace.TracesPoints.push_back(intersT[3]);
                                    trace.CoupleIdTips[trace.TracesNumber]= true;
                                    trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                    trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                    trace.LengthsTrace.push_back(l13);
                                    trace.TracesNumber=trace.TracesNumber+1;

                                }
                                else{
                                    trace.TracesPoints.push_back(intersT[0]);
                                    trace.TracesPoints.push_back(intersT[1]);
                                    trace.CoupleIdTips[trace.TracesNumber]= true;
                                    trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                    trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                    trace.LengthsTrace.push_back(l01);
                                    trace.TracesNumber=trace.TracesNumber+1;
                                }
                            }
                        }
                        if(l13<=tol && l02>tol){
                            Vector3d u=intersT[2]-intersT[0];
                            Vector3d v1=intersT[2]-intersT[1];
                            Vector3d v2=intersT[0]-intersT[1];
                            double P1=u.dot(v1);
                            double P2=u.dot(v2);
                            if(P1*P2>tol){
                                if(l01>l12){
                                    trace.TracesPoints.push_back(intersT[1]);
                                    trace.TracesPoints.push_back(intersT[2]);
                                    trace.CoupleIdTips[trace.TracesNumber]= true;
                                    trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                    trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                    trace.LengthsTrace.push_back(l12);
                                    trace.TracesNumber=trace.TracesNumber+1;

                                }
                                else{
                                    trace.TracesPoints.push_back(intersT[0]);
                                    trace.TracesPoints.push_back(intersT[1]);
                                    trace.CoupleIdTips[trace.TracesNumber]= true;
                                    trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                    trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                    trace.LengthsTrace.push_back(l01);
                                    trace.TracesNumber=trace.TracesNumber+1;
                                }
                            }
                        }
                        if(l02>tol && l03>tol && l12>tol && l13>tol){
                            if(l02>l03 && l02>l12 && l02>l13){
                                trace.TracesPoints.push_back(intersT[1]);
                                trace.TracesPoints.push_back(intersT[3]);
                                trace.CoupleIdTips[trace.TracesNumber]= true;
                                trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                trace.LengthsTrace.push_back(l13);
                                trace.TracesNumber=trace.TracesNumber+1;

                            }
                            else if(l03>l02 && l03>l12 && l03>l13){
                                trace.TracesPoints.push_back(intersT[1]);
                                trace.TracesPoints.push_back(intersT[2]);
                                trace.CoupleIdTips[trace.TracesNumber]= true;
                                trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                trace.LengthsTrace.push_back(l12);
                                trace.TracesNumber=trace.TracesNumber+1;

                            }
                            else if(l12>l02 && l12>l03 && l12>l13){
                                trace.TracesPoints.push_back(intersT[0]);
                                trace.TracesPoints.push_back(intersT[3]);
                                trace.CoupleIdTips[trace.TracesNumber]= true;
                                trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                trace.LengthsTrace.push_back(l03);
                                trace.TracesNumber=trace.TracesNumber+1;

                            }
                            else if(l13>l02 && l13>l03 && l13>l12){
                                trace.TracesPoints.push_back(intersT[2]);
                                trace.TracesPoints.push_back(intersT[0]);
                                trace.CoupleIdTips[trace.TracesNumber]= true;
                                trace.CoupleFracturesTraces[i].push_back(trace.TracesNumber);
                                trace.CoupleFracturesTraces[j].push_back(trace.TracesNumber);
                                trace.LengthsTrace.push_back(l02);
                                trace.TracesNumber=trace.TracesNumber+1;

                            }
                        }
                    }
                }
            }
        }
    }
    ofstream outFile(fileOutput);
    outFile << "# Number of traces" << endl;

    outFile << trace.TracesNumber << endl;

    unsigned int t=0;
    for (unsigned int w=0;w<trace.TracesNumber;w++){
        outFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        outFile <<"  "<< w <<"; "<< trace.IdFractures[w][0]<<"; "<< trace.IdFractures[w][1]
        <<"; "<< trace.TracesPoints[t][0]<<"; "<<trace.TracesPoints[t][1]<<"; "
        <<trace.TracesPoints[t][2]<< "; "<< trace.TracesPoints[t+1][0]<<"; "
        <<trace.TracesPoints[t+1][1]<<"; "<<trace.TracesPoints[t+1][2]<<endl;
        t=t+2;
    }


    for(auto a=trace.CoupleFracturesTraces.begin(); a!=trace.CoupleFracturesTraces.end();a++){
        outFile << "# FractureId; NumTraces" << endl;
        vector<double> vec;
            auto v=a->second;
        outFile<< a->first<<" "<< v.size()<<endl;
        for(unsigned int b=0; b<v.size();b++){
            vec.push_back(trace.LengthsTrace[v[b]]);
        }
        Sorting(vec);
        for(unsigned int c=0; c<vec.size();c++){
            unsigned int m=0;
            for(unsigned int d=0; vec[c]!=trace.LengthsTrace[d]; d++){
                m++;
            }
            outFile<<"TraceId; Tips; Length"<<endl;
            outFile<<m<<"; "<< trace.CoupleIdTips[m]<<"; "<< vec[c]<<endl;
        }
    }


    return true;
}

void MergeSort(vector<double>& lengths, const unsigned int& sx, const unsigned int& dx){
    if(sx<dx){
        unsigned int cx= (sx+dx)/2;
        MergeSort(lengths, sx, cx);
        MergeSort(lengths, cx+1, dx);
        Merge(lengths, sx, cx, dx);
    }
}

void Merge(vector<double>& lengths, const unsigned int& sx, const unsigned int cx, const unsigned int& dx){
    unsigned int i=sx;
    unsigned int j=cx+1;
    unsigned int k=0;
    vector<double> b;
    b.reserve(dx-sx+1);
    while((i<=cx) && (j<=dx)){
        if(lengths[i]>=lengths[j]){
            b[k]=lengths[i];
            i++;
        }
        else{
            b[k]=lengths[j];
            j++;
        }
        k++;
    }
    for(; i<=cx; i++, k++){
        b[k]=lengths[i];
    }
    for(; j<=dx; j++, k++){
        b[k]=lengths[j];
    }
    for(i=sx; i<=dx; i++){
        lengths[i]=b[i-sx];
    }
}

void Sorting(vector<double>& vec){
    MergeSort(vec, 0, vec.size()-1);
}

// triangolazione dei poligoni per poterli esportare con Paraview

void GedimInterface(Fractures& fracture, vector<vector<unsigned int>>& triangles, VectorXi& materials){

    unsigned int numPoints=0;
    vector<vector<vector<unsigned int>>> triangleList(fracture.FractureNumber);
    for(unsigned int p=0; p<fracture.ListVertices.size(); p++){
        const unsigned int numPolygonVertices=fracture.VerticeNumber[p];
        numPoints+=numPolygonVertices;
        for(unsigned int v=0; v<numPolygonVertices;v++){
            const unsigned int nextVertex = fracture.ListVertices[p][(v + 1) % numPolygonVertices];
            const unsigned int nextNextVertex = fracture.ListVertices[p][(v + 2) % numPolygonVertices];
            if ((v + 2) % numPolygonVertices == 0){
                break;
            }
            vector<unsigned int> triangle_vertices = {fracture.ListVertices[p][0], nextVertex, nextNextVertex};
            triangleList[p].push_back(triangle_vertices);
        }
    }
    fracture.VerticesCoordinates=MatrixXd::Zero(3,numPoints);
    for(unsigned int d=0; d<numPoints;d++){
        for(unsigned int e=0; e<3;e++){
            fracture.VerticesCoordinates(e,d)=fracture.Coordinates[d][e];
        }
    }
    unsigned int numTotalTriangles= 0;
    for(unsigned int r = 0; r < fracture.FractureNumber; r++){
        numTotalTriangles+= triangleList[r].size();
    }
    triangles.reserve(numTotalTriangles);
    materials= VectorXi::Zero(numTotalTriangles);
    unsigned int count= 0;
    for(unsigned int m= 0; m< fracture.FractureNumber; m++){
        for(unsigned int t = 0; t < triangleList[m].size(); t++){
            triangles.push_back(triangleList[m][t]);
            materials(count)= m;
            count++;
        }
    }


}


}


