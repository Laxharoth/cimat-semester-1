#include "real_matrix.hpp"
#include <iostream>

void printm(const mymtx::RealMatrix &m){
    for (size_t i = 0; i < m.shape_y; ++i)
    {
        for( auto j = m.begin(i); j < m.end(i); ++j){
            std::cout << *j << " ";
        }
        std::cout << std::endl;
    }
    
}
void printm(const mymtx::MatrixTraspose &m){
    for (size_t i = 0; i < m.shape_y; ++i)
    {
        for(size_t j = 0; j < m.shape_x; ++j){
            std::cout << m(i,j) << " ";
        }
        std::cout << std::endl;
    }
    
}
void printv(const mymtx::RealVector &m){
    for( auto j = m.begin(); j < m.end(); ++j){
            std::cout << *j << " ";
    }
    std::cout << std::endl;
}