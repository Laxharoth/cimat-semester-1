#ifndef MATRIX_ITERATOR_TCC
#define MATRIX_ITERATOR_TCC
#include "matrix_like.tcc"
#include <iterator>

namespace mymtx{
template<class T>
class array_like_iterator : public std::iterator<std::forward_iterator_tag, T>{
    array_like_iterator(){}
protected:
    size_t col;
    array_like<T> *data;
    unsigned int *ref_counter;
public:
    array_like_iterator(size_t col, array_like<T> *data) 
        :col(col), data(data) { ref_counter = new unsigned int; *ref_counter = 1; }
    array_like_iterator(const array_like_iterator<T>& mit)
        :col(mit.col), data(mit.data),ref_counter(mit.ref_counter) { (*(ref_counter))++; }
    ~array_like_iterator(){
        (*ref_counter)--;
        if(*ref_counter == 0){
            delete ref_counter;
            delete data;
        }
    }
    size_t get_row() const { return data->get_row(); }
    size_t get_col() const { return col; }
    array_like_iterator<T>& operator++() {++col;return *this;}
    array_like_iterator<T> operator++(int)
    {   array_like_iterator<T> tmp(*this); operator++(); return tmp;    }
    array_like_iterator<T>& operator--() {--col;return *this;}
    array_like_iterator<T> operator--(int)
    {   array_like_iterator<T> tmp(*this); operator--(); return tmp;    }
    array_like_iterator<T>& operator+=(int c) {col+=c;return *this;}
    array_like_iterator<T> operator+(int c)
    {   array_like_iterator<T> tmp(*this); tmp+=c; return tmp;    }
    array_like_iterator<T>& operator-=(int c) {col-=c;return *this;}
    array_like_iterator<T> operator-(int c)
    {   array_like_iterator<T> tmp(*this); tmp-=c; return tmp;    }
    bool operator==(const array_like_iterator<T>& rhs) const 
    {   return col==rhs.col;  }
    bool operator!=(const array_like_iterator<T>& rhs) const 
    {   return col!=rhs.col;}
    bool operator>(const array_like_iterator<T>& rhs) const 
    {   return col>rhs.col;  }
    bool operator<(const array_like_iterator<T>& rhs) const 
    {   return col<rhs.col;  }
    bool operator>=(const array_like_iterator<T>& rhs) const 
    {   return col>=rhs.col;  }
    bool operator<=(const array_like_iterator<T>& rhs) const 
    {   return col<=rhs.col;  }
    T& operator*() {return (*data)[col];}
};
#include <cstdio>
template<class T>
class array_like_vert_iterator : public std::iterator<std::forward_iterator_tag, T>{
    array_like_vert_iterator(){}
protected:
    size_t row;
    size_t col;
    array_like<T> *data;
    unsigned int *ref_counter;
public:
    array_like_vert_iterator(size_t col, array_like<T> *data)
        :col(col), row(0), data(data){ref_counter = new unsigned int; *ref_counter = 1;}
    array_like_vert_iterator(size_t col, size_t row, array_like<T> *data) 
        :col(col), row(row), data(data){ref_counter = new unsigned int; *ref_counter = 1;}
    array_like_vert_iterator(const array_like_vert_iterator<T>& mit)
        :col(mit.col), row(mit.row), data(mit.data),ref_counter(mit.ref_counter) {(*(ref_counter))++;}
    ~array_like_vert_iterator(){
        (*ref_counter)--;
        if(*ref_counter == 0){
            delete ref_counter;
            delete data;
        }
    }
    size_t get_row() const { return data->get_row(); }
    size_t get_col() const { return col; }
    array_like_vert_iterator<T>& operator++() {data->set_row(++row);return *this;}
    array_like_vert_iterator<T> operator++(int)
    {   array_like_vert_iterator<T> tmp(*this); operator++(); return tmp;    }
    array_like_vert_iterator<T>& operator--() {data->set_row(--row);return *this;}
    array_like_vert_iterator<T> operator--(int)
    {   array_like_vert_iterator<T> tmp(*this); operator--(); return tmp;    }
    array_like_vert_iterator<T>& operator+=(int r){ data->set_row((row+=r));return *this;}
    array_like_vert_iterator<T> operator+(int r)
    {   array_like_vert_iterator<T> tmp(*this); return (tmp+=r); }
    array_like_vert_iterator<T>& operator-=(int r) {data->set_row((row-=r));return *this;}
    array_like_vert_iterator<T> operator-(int r)
    {   array_like_vert_iterator<T> tmp(*this); return (tmp-=r);    }
    bool operator==(const array_like_vert_iterator<T>& rhs) const 
    {   return row==rhs.row;  }
    bool operator!=(const array_like_vert_iterator<T>& rhs) const 
    {   return row!=rhs.row;}
    bool operator>(const array_like_vert_iterator<T>& rhs) const 
    {   return row>rhs.row;  }
    bool operator<(const array_like_vert_iterator<T>& rhs) const 
    {   return row<rhs.row;  }
    bool operator>=(const array_like_vert_iterator<T>& rhs) const 
    {   return row>=rhs.row;  }
    bool operator<=(const array_like_vert_iterator<T>& rhs) const 
    {   return row<=rhs.row;  }
    T& operator*() {return (*data)[col];}
};
}
#endif /* MATRIX_ITERATOR_TCC */
