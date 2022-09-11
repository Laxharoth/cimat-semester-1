#include "quicksort.c"
#include "pgm1.c"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define char unsigned char

int cmpchr(void *a, void *b);
void swapchr(void *a, void *b);
char median( char *mtx,int y, int x, int dim_y, int dim_x );
char mean( char *mtx,int y, int x, int dim_y, int dim_x );
char gradient( char *mtx,int y, int x, int dim_y, int dim_x );
void replace( char *mtx,int dim_y, int dim_x, char (*fnc)( char *mtx,int y, int x, int dim_y, int dim_x ) );

int main(int argc, char const *argv[]){
    // char c = 3;
    // char *mtx = malloc( 6*6 );
    // for( int i = 0; i < 6; i++ )
    // for( int j = 0; j < 6; j++ ) mtx[i*6 + j] = (c*=3);
    // for( int i = 0; i < 6; i++ ){   
    // for( int j = 0; j < 6; j++ ){  printf("%d ", mtx[i*6 + j]); }printf("\n"); }
    

    // replace( mtx, 6,6,gradient );
    // for( int i = 0; i < 6; i++ ){   
    // for( int j = 0; j < 6; j++ ){  printf("%d ", mtx[i*6 + j]); }printf("\n"); }

    int cols, rows;
    char *img = pgmRead("fractal_treebarbarabarbara.ascii.pgm",&rows, &cols);
    replace( img, rows,cols,median );
    pgmWrite("fractal_treebarbarabarbara.2.ascii.pgm",rows,cols,img,"");

    return 0;
}


int cmpchr(void *a, void *b){ return (*((char*)(a)))>(*((char*)(b)))?1:-1; }
void swapchr(void *a, void *b){ char c = (*((char*)(a))); (*((char*)(a))) = (*((char*)(b))); (*((char*)(b)))-c; }
char median( char *mtx,int y, int x, int dim_y, int dim_x ){
    if( (x==0 && y==0) ||
        (x==0 && y==dim_y-1)|| 
        (x==dim_x-1 && y==0) ||
        (x==dim_x-1 && y==dim_y-1)
    ) return mtx[y*dim_x+x];
    if( x == 0 || x == dim_x-1){
        char ptr[3];
        for (size_t i = y-1,ii=0; i <= y+1; i++,ii++){
            ptr[ii] = mtx[i*dim_x+x];
        }
        quicksort(ptr, 0, 3, 1,cmpchr, swapchr);
        return ptr[1] ;    
    }
    if( y == 0 || y == dim_y-1){
        char ptr[3];
        for (size_t i = x-1,ii=0; i <= x+1; i++,ii++){
            ptr[ii] = mtx[y*dim_x+i];
        }
        quicksort(ptr, 0, 3, 1,cmpchr, swapchr);
        return ptr[1] ;    
    }
    char ptr[9];
    for (size_t i = y-1,ii=0; i <= y+1; i++,ii++){
        for(size_t j = x-1,jj=0; j <= x+1; j++,jj++){
            ptr[ii*3+jj] = mtx[i*dim_x+j];
        }
    }
    quicksort(ptr, 0, 9, 1,cmpchr, swapchr);
    return ptr[4] ;
}
char mean( char *mtx,int y, int x, int dim_y, int dim_x ){
    int sum = 0;
    if( (x==0 && y==0) ||
        (x==0 && y==dim_y-1)|| 
        (x==dim_x-1 && y==0) ||
        (x==dim_x-1 && y==dim_y-1)
    ) return mtx[y*dim_x+x];
    if( x == 0 || x == dim_x-1){
        for (size_t i = y-1; i <= y+1; i++){
            sum+= mtx[i*dim_x+x];
        }
        return (char)(sum/3);
    }
    if( y == 0 || y == dim_y-1){
        for (size_t i = x-1; i <= x+1; i++){
            sum += mtx[y*dim_x+i];
        }
        return (char)(sum/3);
    }
    for (size_t i = y-1; i <= y+1; i++){
        for(size_t j = x-1; j <= x+1; j++){
            sum += mtx[i*dim_x+j];
        }
    }
    return (char)(sum/9);
}
char gradient( char *mtx,int y, int x, int dim_y, int dim_x ){
    int sum = 0;
    if( (x==0 && y==0) ||
        (x==0 && y==dim_y-1)|| 
        (x==dim_x-1 && y==0) ||
        (x==dim_x-1 && y==dim_y-1)
    ) return mtx[y*dim_x+x];
    if( x == 0 || x == dim_x-1){
        double gy = mtx[(y+1)*dim_x+x]-mtx[y*dim_x+x];
        return (char)(abs(gy));
    }
    if( y == 0 || y == dim_y-1){
        double gx = mtx[y*dim_x+x+1]-mtx[y*dim_x+x];
        return (char)(abs(gx));
    }
    double gy = mtx[(y+1)*dim_x+x]-mtx[y*dim_x+x];
    double gx = mtx[y*dim_x+x+1]-mtx[y*dim_x+x];
    return (char)( sqrt( gy*gy + gx*gx ) );
}
void replace( char *mtx,int dim_y, int dim_x, char (*fnc)( char *mtx,int y, int x, int dim_y, int dim_x ) ){
    for (size_t i = 0; i < dim_y; i++){
        for (size_t j = 0; j < dim_x; j++){
            mtx[i*dim_x+j] = fnc(mtx,i,j,dim_y,dim_x);
        }
    }   
}