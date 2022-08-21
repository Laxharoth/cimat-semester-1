#ifndef MATRIX_EXCEPTION_TCC
#define MATRIX_EXCEPTION_TCC
#include <exception>
#include <cstdlib>
class wrong_shape_exception : public std::exception{
    size_t m1_shape_y;
    size_t m1_shape_x;
    size_t m2_shape_y;
    size_t m2_shape_x;
    public:
    wrong_shape_exception(size_t m1_shape_y, size_t m1_shape_x, size_t m2_shape_y, size_t m2_shape_x):
        m1_shape_y(m1_shape_y), m2_shape_y(m2_shape_y), 
        m1_shape_x(m1_shape_x), m2_shape_x(m2_shape_x) {};
};

#endif /* MATRIX_EXCEPTION_TCC */
