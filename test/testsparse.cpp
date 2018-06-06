//
// This file is part of the course materials for AMATH483/583 at the University
// of Washington,
// Spring 2017
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0
// International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#include "AOSMatrix.hpp"
#include "COOMatrix.hpp"
#include "CSRMatrix.hpp"
#include "Matrix.hpp"
#include "Timer.hpp"
#include "Vector.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  size_t nPoints = 32;
  size_t matDim = nPoints * nPoints;

  Matrix A(matDim, matDim);
  piscetize(A, nPoints, nPoints);
  COOMatrix ACOO(matDim, matDim);
  piscetize(ACOO, nPoints, nPoints);
  CSRMatrix ACSR(matDim, matDim);
  piscetize(ACSR, nPoints, nPoints);
  AOSMatrix AAOS(matDim, matDim);
  piscetize(AAOS, nPoints, nPoints);
  Vector x = Vector(matDim);
  randomize(x);

  Vector y0 = A * x;
  Vector y1 = ACOO * x;
  Vector y2 = ACSR * x;
  Vector y3 = AAOS * x;

  cout << "|| y0 - y1 || = " << twoNorm(y0 - y1) << endl;
  cout << "|| y0 - y2 || = " << twoNorm(y0 - y2) << endl;
  cout << "|| y0 - y3 || = " << twoNorm(y0 - y3) << endl;

  return 0;
}
