#include "matrix_exception.hpp"

wrong_shape_exception::wrong_shape_exception(size_t m1_shape_y, size_t m1_shape_x, size_t m2_shape_y, size_t m2_shape_x){
    _what = new char[100];
    snprintf(_what, 100, "wrong shape: shape1 = ( %d, %d ), shape2 = ( %d, %d )"
        , m1_shape_y, m1_shape_x, m2_shape_y, m2_shape_x);
};
wrong_shape_exception::~wrong_shape_exception(){
    delete[] _what;
}
const char *wrong_shape_exception::what() const throw(){
    return _what;
}
