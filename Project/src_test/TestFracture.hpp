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
    EXPECT_DOUBLE_EQ(length, std::sqrt(27.0));
}

//********************************

TEST(SegmentSort, TestSorting) {

    vector<double> C = {1, 4, 3, 2, 6 , 5};

    Sorting(C);

    vector<double> vet={6,5,4,3,2,1};

    EXPECT_EQ(C, vet);


}
}

#endif
