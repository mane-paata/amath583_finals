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
#include "Timer.hpp"
#include <fstream>
#include <functional>
#include <iostream>

using namespace std;

double benchmark(size_t M, size_t N, size_t K, size_t numruns,
                 function<void(const Matrix &, const Matrix &, Matrix &)>);
void runBenchmark(function<void(const Matrix &, const Matrix &, Matrix &)> f,
                  size_t maxsize);

int main(int argc, char *argv[]) {
  size_t maxsize = (argc == 3) ? stod(argv[2]) : 4;

  if (string(argv[1]) == "mult")
    runBenchmark(basicMultiply, maxsize);
  else if (string(argv[1]) == "basicThreadedMultiply")
    runBenchmark(basicThreadedMultiply, maxsize);
  else if (string(argv[1]) == "hoistedmult")
    runBenchmark(hoistedMultiply, maxsize);
  else if (string(argv[1]) == "2x2")
    runBenchmark(tiledMultiply2x2, maxsize);
  else if (string(argv[1]) == "2x4")
    runBenchmark(tiledMultiply2x4, maxsize);
  else if (string(argv[1]) == "4x2")
    runBenchmark(tiledMultiply4x2, maxsize);
  else if (string(argv[1]) == "4x4")
    runBenchmark(tiledMultiply4x4, maxsize);
  else if (string(argv[1]) == "blocked")
    runBenchmark(blockedTiledMultiply2x2, maxsize);
  else if (string(argv[1]) == "copyblocked")
    runBenchmark(copyBlockedTiledMultiply2x2, maxsize);
  else if (string(argv[1]) == "hoisted")
    runBenchmark(hoistedTiledMultiply2x2, maxsize);
  else if (string(argv[1]) == "blockhoisted")
    runBenchmark(hoistedBlockedTiledMultiply2x2, maxsize);
  else if (string(argv[1]) == "copyblockhoisted")
    runBenchmark(hoistedCopyBlockedTiledMultiply2x2, maxsize);
  else if (string(argv[1]) == "copyblockhoisted4x4")
    runBenchmark(hoistedCopyBlockedTiledMultiply4x4, maxsize);
#ifdef __AVX__
  else if (string(argv[1]) == "copyblockhoisted4x4AVX")
    runBenchmark(hoistedCopyBlockedTiledMultiply4x4AVX, maxsize);
#endif // __AVX__
  else
    return -1;

  return 0;
}

void runBenchmark(function<void(const Matrix &, const Matrix &, Matrix &)> f,
                  size_t maxsize) {
  cout << "N\tN*N\tTime\tFlops" << endl;
  for (size_t i = 8; i <= maxsize; i *= 2) {
    size_t numruns = 8L * 1048L * 1048L * 1048L / (i * i * i) + 2;
    double t = benchmark(i, i, i, numruns, f);
    cout << i << "\t" << i * i << "\t" << t << "\t"
         << 2.0 * 1.e3 * numruns * i * i * i / t << endl;
  }
}

double benchmark(size_t M, size_t N, size_t K, size_t numruns,
                 function<void(const Matrix &, const Matrix &, Matrix &)> f) {
  Matrix A(M, K), B(K, N), C(M, N);
  randomize(A);
  randomize(B);
  randomize(C);

  Timer T;
  T.start();
  for (size_t i = 0; i < numruns; ++i) {
    f(A, B, C);
  }
  T.stop();

  zeroize(C);

  return T.elapsed();
}
