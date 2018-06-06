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

#include "Vector.hpp"
#include <iostream>
using namespace std;

int main() {

  Vector x(16);
  randomize(x);
  for (size_t i = 0; i < x.num_rows(); ++i) {
    cout << x(i) << endl;
  }
  return 0;
}
