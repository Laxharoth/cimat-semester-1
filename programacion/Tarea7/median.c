#include "stdlib.h"
#include "stdio.h"
#include "quicksort.h"

typedef char (*filter)( char *mtx, //img
             int y, int x, //position to find median
             int window_size, //window size
             int dim_y, int dim_x //img dimensions
             );

int cmpchr(void *a, void *b){ return (*((char*)(a)))>(*((char*)(b)))?1:-1; }
void swapchr(void *a, void *b){ char c = (*((char*)(a))); (*((char*)(a))) = (*((char*)(b))); (*((char*)(b)))-c; }
char median( char *mtx, //img
             int y, int x, //position to find median
             int window_size, //window size
             int dim_y, int dim_x //img dimensions
             ){
    char ptr[window_size*window_size];
    int window_size_x = window_size;
    int window_size_y = window_size;
    int x_begin,x_end, y_begin, y_end;
    if( 2*x < window_size ){ window_size_x = x * 2 + 1; x_end = x; }
    if( 2*y < window_size ){ window_size_y = y*2 + 1; }
    if( 2*(x -dim_x ) > window_size ){ window_size_x =  (dim_x - x)*2 + 1; }
    if( 2*(y -dim_y ) > window_size ){ window_size_y =  (dim_y - y)*2 + 1; }
    x_begin = x - window_size_x / 2;
    x_end   = x + window_size_x / 2;
    y_begin = y - window_size_y / 2;
    y_end   = y + window_size_y / 2;
    char a = 0;
    for (size_t y_pos = y_begin, y_ptr = 0; y_pos <= y_end; ++y_pos, ++y_ptr){
        for(size_t x_pos = x_begin, x_ptr = 0; x_pos <= x_end; ++x_pos,++x_ptr,++a){
            ptr[y_ptr*window_size_x+x_ptr] = mtx[y_pos*dim_x+x_pos];
        }
    }
    quicksort(ptr, 0, a, 1,cmpchr, swapchr);
    return ptr[window_size_x*window_size_y/2] ;
}
void replace( char *mtx,int dim_y, int dim_x, filter fnc, int window_size ){
    for (size_t i = 0; i < dim_y; i++){
        for (size_t j = 0; j < dim_x; j++){
            mtx[i*dim_x+j] = fnc(mtx,i,j,window_size,dim_y,dim_x);
        }
    }   
}
