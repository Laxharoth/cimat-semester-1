#ifndef MATRIZ_BANDA_HPP
#define MATRIZ_BANDA_HPP

class MatrizBanda{
    class row_wrapper;
    double** matriz;
    row_wrapper *wrapper;
    int left{}, right{}, size{};
class row_wrapper{
    static double default_value;
    double** &matriz;
    int &left, &right, &size;
    int row;
public:
    row_wrapper(double** &matriz, int &left, int &right, int &size, int row);
    double &get(const int &col);
    double & operator[](const int &col);
    friend class MatrizBanda;
};
public:
    row_wrapper &get(const int &row);
    row_wrapper &operator[](const int &row);
    MatrizBanda(const int &left, const int &right, const int &size);
    ~MatrizBanda();
};

#endif /* MATRIZ_BANDA_HPP */
