//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//
#ifndef __COO_HPP
#define __COO_HPP

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "Vector.hpp"

class COOMatrix {
public:
  COOMatrix(size_t M, size_t N) : num_rows_(M), num_cols_(N) {}

  void push_back(size_t i, size_t j, double val) {
    assert(i < num_rows_ && i >= 0);
    assert(j < num_cols_ && j >= 0);

    row_indices_.push_back(i);
    col_indices_.push_back(j);
    storage_.push_back(val);
  }

  void clear() {
    row_indices_.clear();
    col_indices_.clear();
    storage_.clear();
  }

  void reserve(size_t n) {
    assert(n >= 0);

    row_indices_.reserve(n);
    col_indices_.reserve(n);
    storage_.reserve(n);
  }

  size_t num_rows() const { return num_rows_; }
  size_t num_cols() const { return num_cols_; }
  size_t num_nonzeros() const { return storage_.size(); }

  void matvec(const Vector& x, Vector& y) const {
    for (size_t k = 0; k < storage_.size(); ++k) {
      y(row_indices_[k]) += storage_[k] * x(col_indices_[k]);
    }
  }

  void streamMatrix(std::ostream& outputFile) const {
    assert(storage_.size() == row_indices_.size() && storage_.size() == col_indices_.size());

    outputFile << "AMATH 583 COOMATRIX" << std::endl;
    outputFile << num_rows_ << " " << num_cols_ << std::endl;

    // Write data
    for (size_t i = 0; i < storage_.size(); ++i) {
      outputFile << row_indices_[i] << " ";
      outputFile << col_indices_[i] << " ";
      outputFile << storage_[i] << " ";
      outputFile << std::endl;
    }

    // Write tailer
    outputFile << "THIS IS THE END" << std::endl;
  }

private:
  size_t              num_rows_, num_cols_;
  std::vector<size_t> row_indices_, col_indices_;
  std::vector<double> storage_;
};

Vector operator*(const COOMatrix& A, const Vector& x);
void   matvec(const COOMatrix& A, const Vector& x, Vector& y);
void   piscetize(COOMatrix& A, size_t xpoints, size_t ypoints);
void   writeMatrix(const COOMatrix& A, const std::string& filename);
void   streamMatrix(const COOMatrix& A);
void   streamMatrix(const COOMatrix& A, std::ostream& outputFile);

#endif    // __COO_HPP
