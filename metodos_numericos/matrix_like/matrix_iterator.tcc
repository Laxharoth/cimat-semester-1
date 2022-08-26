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
public:
    array_like_iterator(size_t col, array_like<T> *data) 
        :col(col), data(data) {}
    array_like_iterator(const array_like_iterator<T>& mit)
        :col(mit.col), data(mit.data) {}
    size_t get_row() const { return data->get_row(); }
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
    {   return data==rhs.data && col==rhs.col;  }
    bool operator!=(const array_like_iterator<T>& rhs) const 
    {   return data==rhs.data && col!=rhs.col;}
    bool operator>(const array_like_iterator<T>& rhs) const 
    {   return data==rhs.data && col>rhs.col;  }
    bool operator<(const array_like_iterator<T>& rhs) const 
    {   return data==rhs.data && col<rhs.col;  }
    bool operator>=(const array_like_iterator<T>& rhs) const 
    {   return data==rhs.data && col>=rhs.col;  }
    bool operator<=(const array_like_iterator<T>& rhs) const 
    {   return data==rhs.data && col<=rhs.col;  }
    T& operator*() {return (*data)[col];}
};
}
#endif /* MATRIX_ITERATOR_TCC */
