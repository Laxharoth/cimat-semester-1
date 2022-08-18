#ifndef MATRIZ_BANDA_HPP
#define MATRIZ_BANDA_HPP

class MatrizBanda{
    double** matriz;
    int left{}, right{}, size{};
class row_wrapper{
    static double default_value;
    double** &matriz;
    int &left, &right;
    int row;
public:
    row_wrapper(double** &matriz, int &left, int &right, int row);
    double &get(const int &col);
    double & operator[](const int &col);
};
public:
    double &get(const int &row, const int &col);
    row_wrapper operator[](const int &row);
    MatrizBanda(const int &left, const int &right, const int &size);
    ~MatrizBanda();
};

#endif /* MATRIZ_BANDA_HPP */
