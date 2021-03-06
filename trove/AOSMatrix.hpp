//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//
#ifndef __AOS_HPP
#define __AOS_HPP

#include <cassert>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "Vector.hpp"

class AOSMatrix {
private:
  typedef std::tuple<size_t, size_t, double> element;

public:
  AOSMatrix(size_t M, size_t N) : num_rows_(M), num_cols_(N) {}

  void push_back(size_t i, size_t j, double val) {
    assert(i < num_rows_ && i >= 0);
    assert(j < num_cols_ && j >= 0);

    storage_.push_back(element(i, j, val));
  }

  void clear() { storage_.clear(); }

  void reserve(size_t n) {
    assert(n >= 0);

    storage_.reserve(n);
  }

  size_t num_rows() const { return num_rows_; }
  size_t num_cols() const { return num_cols_; }
  size_t num_nonzeros() const { return storage_.size(); }

  void matvec(const Vector& x, Vector& y) const {
    for (size_t k = 0; k < storage_.size(); ++k) {
      y(std::get<1>(storage_[k])) += std::get<2>(storage_[k]) * x(std::get<0>(storage_[k]));
    }
  }

  void streamMatrix(std::ostream& outputFile) const {

    outputFile << "AMATH 583 AOSMATRIX" << std::endl;
    outputFile << num_rows_ << " " << num_cols_ << std::endl;

    // Write data
    for (size_t i = 0; i < storage_.size(); ++i) {
      outputFile << std::get<0>(storage_[i]) << " ";
      outputFile << std::get<1>(storage_[i]) << " ";
      outputFile << std::get<2>(storage_[i]) << " ";
      outputFile << std::endl;
    }

    // Write tailer
    outputFile << "THIS IS THE END" << std::endl;
  }

private:
  size_t               num_rows_, num_cols_;
  std::vector<element> storage_;
};

Vector operator*(const AOSMatrix& A, const Vector& x);
void   matvec(const AOSMatrix& A, const Vector& x, Vector& y);
void   piscetize(AOSMatrix& A, size_t xpoints, size_t ypoints);
void   writeMatrix(const AOSMatrix& A, const std::string& filename);
void   streamMatrix(const AOSMatrix& A);
void   streamMatrix(const AOSMatrix& A, std::ostream& outputFile);

#endif    // __AOS_HPP
