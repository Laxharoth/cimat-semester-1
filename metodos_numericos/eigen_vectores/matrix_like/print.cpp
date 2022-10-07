#include "matrix.hpp"
#include <iostream>

void printm(const mymtx::matrix &m, std::ostream &out) {
  for (size_t i = 0; i < m.shape_y; ++i) {
    for (auto j = m.begin(i); j < m.end(i); ++j) {
      out << *j << " ";
    }
    out << std::endl;
  }
}
void printm(const mymtx::matrix &m) { printm(m, std::cout); }
void printm(const mymtx::MatrixTraspose &m, std::ostream &out) {
  for (size_t i = 0; i < m.shape_y; ++i) {
    for (size_t j = 0; j < m.shape_x; ++j) {
      out << m(i, j) << " ";
    }
    out << std::endl;
  }
}
void printm(const mymtx::MatrixTraspose &m) { printm(m, std::cout); }
void printv(const mymtx::vector &m, std::ostream &out) {
  for (auto j = m.begin(); j < m.end(); ++j) {
    out << *j << " ";
  }
  out << std::endl;
}
void printv(const mymtx::vector &m) { printv(m, std::cout); }