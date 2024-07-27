#ifndef __TESTFRACTURES_H
#define __TESTFRACTURES_H

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "UCDUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace fractureLibrary {

//********************************

TEST(SegmentLengthTest, TestComputeLengths) {


    Vector3d a(0.0, 0.0, 0.0);
    Vector3d b(1.0, 0.0, 0.0);

    double length = ComputeLengths(a, b);
    EXPECT_DOUBLE_EQ(length, 1.0);


    Vector3d c(0.0, 0.0, 0.0);
    Vector3d d(0.0, 1.0, 0.0);

    length = ComputeLengths(c, d);
    EXPECT_DOUBLE_EQ(length, 1.0);


    Vector3d e(0.0, 0.0, 0.0);
    Vector3d f(0.0, 0.0, 1.0);



    length = ComputeLengths(e, f);
    EXPECT_DOUBLE_EQ(length, 1.0);


    Vector3d g(1.0, 2.0, 3.0);
    Vector3d h(4.0, 5.0, 6.0);

    length = ComputeLengths(g, h);
    EXPECT_DOUBLE_EQ(length, sqrt(27.0));
}

//********************************

TEST(Sort, TestSorting) {

    vector<double> C = {1, 4, 3, 2, 6 , 5};

    Sorting(C);


    vector<double> vet = {6,5,4,3,2,1};

    EXPECT_EQ(C, vet);

}

//********************************


TEST(DefineTracesTest, BasicTest) {
    Fractures fracture;
    Traces trace;
    trace.TracesNumber = 0;


    fracture.FractureNumber = 4;
    fracture.VerticeNumber = {4, 4, 4, 4};
    fracture.Id = {1, 2, 3, 4};
    fracture.Coordinates = {
        Vector3d(0, 0, 0), Vector3d(0, 0.5, 0), Vector3d(-0.5, 0.5, 1), Vector3d(-0.5, 0, 1), // F1
        Vector3d(-0.5, 0, 0.5), Vector3d(1, 0, 0.5), Vector3d(1, 0.5, 0.5), Vector3d(-0.5, 0.5, 0.5), // F2
        Vector3d(-1, -0.5, 1.5), Vector3d(1.5, -0.5, 1.5), Vector3d(1.5, 1, 1.5), Vector3d(-1, 1, 1.5), // F3
        Vector3d(1, 0, 1), Vector3d(1, 1.5, 1), Vector3d(1, 1.5, 1.5), Vector3d(1, 0, 1.5) // F4
    };

    ComputeSegments(fracture);

    DefineTraces("output.txt", fracture, trace);


    EXPECT_EQ(trace.TracesNumber, 2);


    EXPECT_EQ(trace.LengthsTrace[0], 0.5);
    EXPECT_EQ(trace.LengthsTrace[1], 1.0);
}

//********************************

TEST(GedimTest, TestG){

    Fractures fracture;
    fracture.ListVertices = {{0,1,2,3}, {4,5,6,7}, {8,9,10,11} , {12,13,14,15}};
    vector<vector<unsigned int>> triangles;
    VectorXi materials;

    fracture.FractureNumber = 4;
    fracture.VerticeNumber = {4, 4, 4, 4};
    fracture.Coordinates = {
        Vector3d(0, 0, 0), Vector3d(0, 0.5, 0), Vector3d(-0.5, 0.5, 1), Vector3d(-0.5, 0, 1), // F1
        Vector3d(-0.5, 0, 0.5), Vector3d(1, 0, 0.5), Vector3d(1, 0.5, 0.5), Vector3d(-0.5, 0.5, 0.5), // F2
        Vector3d(-1, -0.5, 1.5), Vector3d(1.5, -0.5, 1.5), Vector3d(1.5, 1, 1.5), Vector3d(-1, 1, 1.5), // F3
        Vector3d(1, 0, 1), Vector3d(1, 1.5, 1), Vector3d(1, 1.5, 1.5), Vector3d(1, 0, 1.5) // F4
    };

    Gedim::UCDUtilities exporter;
    GedimInterface(fracture, triangles, materials);
    exporter.ExportPolygons("./polygons_4TEST.inp",fracture.VerticesCoordinates , triangles,{},{}, materials);







}

}

#endif
