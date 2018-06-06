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

#include "Matrix.hpp"
#include <functional>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

std::vector<
    std::pair<std::function<void(const Matrix &A, const Matrix &B, Matrix &C)>,
              const char *>>
    functionList{
        {basicMultiply, "basicMultiply"},
        {basicThreadedMultiply, "basicThreadedMultiply"},
        {hoistedMultiply, "hoistedMultiply"},
        {tiledMultiply2x2, "tiledMultiply2x2"},
        {hoistedTiledMultiply2x2, "hoistedTiledMultiply2x2"},
        {blockedTiledMultiply2x2, "blockedTiledMultiply2x2"},
        {tiledMultiply2x4, "tiledMultiply2x4"},
        {tiledMultiply4x2, "tiledMultiply4x2"},
        {tiledMultiply4x4, "tiledMultiply4x4"},
        {copyBlockedTiledMultiply2x2, "copyBlockedTiledMultiply2x2"},
        {hoistedBlockedTiledMultiply2x2, "hoistedBlockedTiledMultiply2x2"},
        {hoistedCopyBlockedTiledMultiply2x2,
         "hoistedCopyBlockedTiledMultiply2x2"},
        {hoistedCopyBlockedTiledMultiply4x4,
         "hoistedCopyBlockedTiledMultiply4x4"},
#ifdef __AVX__
        {hoistedCopyBlockedTiledMultiply2x2AVX,
         "hoistedCopyBlockedTiledMultiply2x2AVX"},
        {hoistedCopyBlockedTiledMultiply4x4AVX,
         "hoistedCopyBlockedTiledMultiply4x4AVX"}
#endif // __AVX__
    };

int main() {

  const size_t N = 8;

  Matrix A(N, N), B(N, N), C(N, N), D(N, N), E(N, N);

  randomize(A);
  randomize(B);
  C = A * B;

  for (auto &pf : functionList) {
    cout << pf.second << " " << frobeniusNorm(A) << " " << frobeniusNorm(B)
         << " " << frobeniusNorm(C) << " ";
    zeroize(D);
    pf.first(A, B, D);
    double nn = frobeniusNorm(C - D);
    cout << nn << endl;
  }

  return 0;
}
