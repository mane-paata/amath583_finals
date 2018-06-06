//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#include "Timer.hpp"
#include "Vector.hpp"
#include <iostream>
#include <mpi.h>
#include <cmath>


double mpiTwoNorm(const Vector& lx) {
  double rho = dot(lx, lx), sigma = 0;
  MPI::COMM_WORLD.Allreduce(&rho, &sigma, 1, MPI::DOUBLE, MPI::SUM);
  return std::sqrt(sigma);
}


int main(int argc, char* argv[]) {
  MPI::Init();

  size_t myrank = MPI::COMM_WORLD.Get_rank();
  size_t mysize = MPI::COMM_WORLD.Get_size();
  size_t losize = 1024;
  size_t loops  = 4;
  if (argc >= 2) losize = std::stol(argv[1]);
  if (argc >= 3) loops = std::stol(argv[2]);
  size_t glsize = losize * mysize;

  double seq = 0, sigma = 0;

  if (0 != myrank) {
    glsize = 0;
  }    
  Vector gx = Vector(glsize);
  randomize(gx);

  std::cout << "vector size " << glsize << " on " << mysize << " nodes" << std::endl;

  Timer t2;
  if (0 == myrank) {
    t2.start();
    seq = std::sqrt(dot(gx, gx));
    t2.stop();
  }

  Vector lx(losize);
  MPI::COMM_WORLD.Scatter(&gx(0), losize, MPI::DOUBLE, 
			  &lx(0), losize, MPI::DOUBLE, 0);
  
  Timer t; t.start();
  for (size_t i = 0; i < loops; ++i) {
    sigma = mpiTwoNorm(lx);
  }
  t.stop();

  if (0 == myrank) {
    std::cout << "#\tglobal\tmpi\tdiff" << std::endl;
    std::cout << mysize << "\t" << seq << "\t" << sigma << "\t" << std::abs(sigma-seq) << "\t" << t2.elapsed() << "\t" << t.elapsed()/loops << std::endl;
  }

  MPI::Finalize();

  return 0;
}
