#include "real_matrix.hpp"
namespace mymtx{
vector_iterator::vector_iterator(double *data, const size_t row, const size_t col):data(data),row(row),col(col){}
vector_iterator::vector_iterator(const vector_iterator &other):data(other.data),row(other.row),col(other.col){}
size_t vector_iterator::get_row() const { return row; }
size_t vector_iterator::get_col() const { return col; }
vector_iterator& vector_iterator::operator++()     {++col;return *this;}
vector_iterator vector_iterator::operator++(int)   {vector_iterator tmp(*this); operator++(); return tmp;}
vector_iterator& vector_iterator::operator--()     {--col;return *this;}
vector_iterator vector_iterator::operator--(int)   {vector_iterator tmp(*this); operator--(); return tmp;}
vector_iterator& vector_iterator::operator+=(int c){col+=c;return *this;}
vector_iterator vector_iterator::operator+(int c) const {vector_iterator tmp(*this); tmp+=c; return tmp;}
vector_iterator& vector_iterator::operator-=(int c){col-=c;return *this; }
vector_iterator vector_iterator::operator-(int c) const {vector_iterator tmp(*this); tmp-=c; return tmp;}
bool vector_iterator::operator==(const vector_iterator& rhs) const { return col==rhs.col; }
bool vector_iterator::operator!=(const vector_iterator& rhs) const { return col!=rhs.col; }
bool vector_iterator::operator>(const vector_iterator& rhs ) const { return col>rhs.col;  }
bool vector_iterator::operator<(const vector_iterator& rhs ) const { return col<rhs.col;  }
bool vector_iterator::operator>=(const vector_iterator& rhs) const { return col>=rhs.col; }
bool vector_iterator::operator<=(const vector_iterator& rhs) const { return col<=rhs.col; }
double& vector_iterator::operator[](int n){ return data[col+n]; }
double& vector_iterator::operator*() {return data[col];}
const_vector_iterator::const_vector_iterator(const double *data, const size_t row, const size_t col):data(data),row(row),col(col){}
const_vector_iterator::const_vector_iterator(const const_vector_iterator &other):data(other.data),row(other.row),col(other.col){}
size_t const_vector_iterator::get_row() const { return row; }
size_t const_vector_iterator::get_col() const { return col; }
const_vector_iterator& const_vector_iterator::operator++()     {++col;return *this;}
const_vector_iterator const_vector_iterator::operator++(int)   {const_vector_iterator tmp(*this); operator++(); return tmp;}
const_vector_iterator& const_vector_iterator::operator--()     {--col;return *this;}
const_vector_iterator const_vector_iterator::operator--(int)   {const_vector_iterator tmp(*this); operator--(); return tmp;}
const_vector_iterator& const_vector_iterator::operator+=(int c){col+=c;return *this;}
const_vector_iterator const_vector_iterator::operator+(int c) const {const_vector_iterator tmp(*this); tmp+=c; return tmp;}
const_vector_iterator& const_vector_iterator::operator-=(int c){col-=c;return *this; }
const_vector_iterator const_vector_iterator::operator-(int c) const {const_vector_iterator tmp(*this); tmp-=c; return tmp;}
bool const_vector_iterator::operator==(const const_vector_iterator& rhs) const { return col==rhs.col; }
bool const_vector_iterator::operator!=(const const_vector_iterator& rhs) const { return col!=rhs.col; }
bool const_vector_iterator::operator>(const const_vector_iterator& rhs ) const { return col>rhs.col;  }
bool const_vector_iterator::operator<(const const_vector_iterator& rhs ) const { return col<rhs.col;  }
bool const_vector_iterator::operator>=(const const_vector_iterator& rhs) const { return col>=rhs.col; }
bool const_vector_iterator::operator<=(const const_vector_iterator& rhs) const { return col<=rhs.col; }
const double& const_vector_iterator::operator[](int n) const{ return data[col+n]; }
const double& const_vector_iterator::operator*() const {return data[col];}
}