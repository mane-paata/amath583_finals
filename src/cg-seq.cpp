//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#include <iostream>
#include "Grid.hpp"
#include "cg.hpp"

int main(int argc, char* argv[]) {
  size_t xsize              = 128, ysize = 128;
  size_t max_iters          = xsize;
  double tol                = 1.E-4;

  if (argc >= 2) xsize = ysize = std::stol(argv[1]);
  if (argc >= 3) max_iters     = std::stol(argv[2]);
  if (argc >= 4) tol           = std::stol(argv[3]);


  Grid X0(xsize, ysize), X1(xsize, ysize);

  for (size_t i = 0; i < xsize+2; ++i) {
    X1(i, 0) = X0(i, 0) = 1.0;
  }

  Stencil A;
  cg(A, X1, X0, max_iters, 1.e-6);

  return 0;
}
