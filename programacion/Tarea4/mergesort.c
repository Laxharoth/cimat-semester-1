#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_ARRAYS 9

void mergesort( int *data, const int start, const int end );
void _merge(int *arr_1,int *arr_2, int size_1, int size_2);
int *merge(int **arr, int N, int *Dim);

int main(int argc, char const *argv[]){
    int **arr = malloc( sizeof( int** ) * N_ARRAYS );
    int *Dim = malloc( sizeof( int* ) * N_ARRAYS );
    int total = 0;
    srand(time(NULL));

    for (size_t i = 0; i < N_ARRAYS; i++){
        const int dim = 5 + rand()%14;
        arr[i] = malloc(sizeof(int)*dim);
        Dim[i] = dim;
        total+=dim;
        for(int j = 0; j < dim; j++){
            arr[i][j] = rand();
        }
        mergesort( arr[i], 0, dim );
    }
    
    int *res = merge(arr, N_ARRAYS, Dim);
    for (size_t i = 0; i < total; i++){
        printf("%d\n", res[i]);
    }
    for (size_t i = 0; i < N_ARRAYS; i++){
        free(arr[i]);
    }
    free(arr);
    free(Dim);
    free(res);

    return 0;
}
void mergesort( int *data, const int start, const int end ){
    const int mid = (start + end) / 2;
    if(end - start > 2){
        mergesort(data, start, mid);
        mergesort(data, mid, end);
    }
    _merge(data + start, data + mid, mid - start, end - mid);
}
void _merge(int *arr_1,int *arr_2, int size_1, int size_2){
    int holder[size_1 + size_2];
    int index_1 = 0, index_2 = 0;
    while( index_1 < size_1 && index_2 < size_2 ){
        if( arr_1[index_1] < arr_2[index_2]){
            holder[ index_1 + index_2 ] = arr_1[index_1];
            index_1++;
        }
        else{
            holder[ index_1 + index_2 ] = arr_2[index_2];
            index_2++;
        }
    }
    while( index_1 < size_1 ){
        holder[ index_1 + index_2 ] = arr_1[index_1];
        index_1++;
    }
    while( index_2 < size_2 ){
        holder[ index_1 + index_2 ] = arr_2[index_2];
        index_2++;
    }
    memcpy( arr_1, holder, sizeof(int) * (size_1+size_2));
}
int *merge(int **arr, int N, int *Dim){
    int n = 1;
    int *arr_1,*arr_2;
    int size_1, size_2;
    int index_1, index_2;
    int *temp;
    arr_1 = arr[0]; size_1 = Dim[0];
    while(n < N){
        arr_2 = arr[n];
        size_2 = Dim[n];
        index_1 = index_2 = 0;
        temp = malloc( sizeof(int) * (size_1 + size_2));
        int index_1 = 0, index_2 = 0;
        while( index_1 < size_1 && index_2 < size_2 ){
            if( arr_1[index_1] < arr_2[index_2]){
                temp[ index_1 + index_2 ] = arr_1[index_1];
                index_1++;
            }
            else{
                temp[ index_1 + index_2 ] = arr_2[index_2];
                index_2++;
            }
        }
        while( index_1 < size_1 ){
            temp[ index_1 + index_2 ] = arr_1[index_1];
            index_1++;
        }
        while( index_2 < size_2 ){
            temp[ index_1 + index_2 ] = arr_2[index_2];
            index_2++;
        }
        if(arr_1 != arr[0])free(arr_1);
        arr_1 = temp; size_1 = size_1 + size_2;
        n++;
    }
    return temp;
}
