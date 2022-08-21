#ifndef MATRIX_EXCEPTION_TCC
#define MATRIX_EXCEPTION_TCC
#include <exception>
#include <cstdlib>
#include <cstdio>

class wrong_shape_exception : public std::exception{
    char * _what;
    public:
    wrong_shape_exception(size_t m1_shape_y, size_t m1_shape_x, size_t m2_shape_y, size_t m2_shape_x);
    ~wrong_shape_exception();
    const char *what() const throw();
};

#endif /* MATRIX_EXCEPTION_TCC */
