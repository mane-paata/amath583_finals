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

#include "CSRMatrix.hpp"
#include "Vector.hpp"

int main() {

  CSRMatrix A(4, 4);
  piscetize(A, 2, 2);
  streamMatrix(A);

  return 0;
}
