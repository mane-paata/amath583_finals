//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//


#include <mpi.h>
#include "Final.hpp"
#include "mpiStencil.hpp"
#include "Stencil.hpp"
#include "Grid.hpp"


double mpi_dot(const Grid& X, const Grid& Y) {
  double sum = 0.0;

  /* write me */

  return sum;
}


size_t ir(const mpiStencil& A, Grid& x, const Grid& b, size_t max_iter, double tol) {

  /* Write me */

  return max_iter;
}


size_t cg(const mpiStencil& A, Grid& x, const Grid& b, size_t max_iter, double tol) {
  
  /* Write me */

  return max_iter;
}
