//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2017
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//
#if !defined(__MATRIX_HPP)
#define __MATRIX_HPP

#include "Vector.hpp"
#include <iostream>
#include <string>
#include <vector>

#if !defined(SLOWMATRIX)
class Matrix {
public:
  explicit Matrix(size_t M, size_t N) : num_rows_(M), num_cols_(N), storage_(num_rows_ * num_cols_) {}
  explicit Matrix(size_t M, size_t N, double init) : num_rows_(M), num_cols_(N), storage_(num_rows_ * num_cols_, init) {}

  double&       operator()(size_t i, size_t j) { return storage_[i * num_cols_ + j]; }
  const double& operator()(size_t i, size_t j) const { return storage_[i * num_cols_ + j]; }

  size_t num_rows() const { return num_rows_; }
  size_t num_cols() const { return num_cols_; }

private:
  size_t              num_rows_, num_cols_;
  std::vector<double> storage_;
};
#else
class Matrix {
public:
  Matrix(size_t rows, size_t cols) : num_rows_(rows), storage_(cols), storage_(num_rows_, std::vector<double>(storage_)) {}
  Matrix(size_t rows, size_t cols, double init)
      : num_rows_(rows), storage_(cols), storage_(num_rows_, std::vector<double>(storage_), init) {}

  double&       operator()(size_t i, size_t j) { return storage_[i][j]; }
  const double& operator()(size_t i, size_t j) const { return storage_[i][j]; }

  size_t num_rows() const { return num_rows_; }
  size_t num_cols() const { return storage_; }

private:
  size_t                           num_rows_, storage_;
  std::vector<std::vector<double>> storage_;
};
#endif

Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
double frobeniusNorm(const Matrix& A);
void   basicMultiply(const Matrix& A, const Matrix& B, Matrix& C);
void   basicThreadedMultiply(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedMultiply(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   blockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply2x4(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply4x2(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C);
void   copyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply2x2AVX(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply4x4AVX(const Matrix& A, const Matrix& B, Matrix& C);
double oneNorm(const Matrix& A);
double infinityNorm(const Matrix& A);
double frobeniusNorm(const Matrix& A);
void   zeroize(Matrix& C);
void   randomize(Matrix& A);
void   piscetize(Matrix& A, size_t xpoints, size_t ypoints);
void   writeMatrix(const Matrix& A, const std::string& filename);
void   streamMatrix(const Matrix& A);
void   streamMatrix(const Matrix& A, std::ostream& outputFile);

Vector operator*(const Matrix& A, const Vector& x);
void   basicMultiplyMV(const Matrix& A, const Vector& x, Vector& y);

#endif    // __MATRIX_HPP
