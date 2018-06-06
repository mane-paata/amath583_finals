//
// This file is part of the course materials for AMATH483/583 at the University
// of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0
// International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#include "AOSMatrix.hpp"
#include "Vector.hpp"
#include <cassert>
#include <cstddef>

Vector operator*(const AOSMatrix& A, const Vector& x) {
  assert(A.num_cols() == x.num_rows());

  Vector y(A.num_rows(), 0.0);
  matvec(A, x, y);

  return y;
}

void matvec(const AOSMatrix& A, const Vector& x, Vector& y) { A.matvec(x, y); }

void piscetize(AOSMatrix& A, size_t xpoints, size_t ypoints) {
  assert(A.num_rows() == A.num_cols());
  assert(xpoints * ypoints == A.num_rows());

  A.clear();

  for (size_t j = 0; j < xpoints; j++) {
    for (size_t k = 0; k < ypoints; k++) {
      size_t jrow = j * ypoints + k;

      if (j != 0) {
        size_t jcol = (j - 1) * ypoints + k;
        A.push_back(jrow, jcol, -1.0);
      }
      if (k != 0) {
        size_t jcol = j * ypoints + (k - 1);
        A.push_back(jrow, jcol, -1.0);
      }

      A.push_back(jrow, jrow, 4.0);

      if (k != ypoints - 1) {
        size_t jcol = j * ypoints + (k + 1);
        A.push_back(jrow, jcol, -1.0);
      }
      if (j != xpoints - 1) {
        size_t jcol = (j + 1) * ypoints + k;
        A.push_back(jrow, jcol, -1.0);
      }
    }
  }
}

void writeMatrix(const AOSMatrix& A, const std::string& filename) {
  std::ofstream outputFile(filename);
  streamMatrix(A, outputFile);
  outputFile.close();
}

void streamMatrix(const AOSMatrix& A) { A.streamMatrix(std::cout); }

void streamMatrix(const AOSMatrix& A, std::ostream& outputFile) { A.streamMatrix(outputFile); }
