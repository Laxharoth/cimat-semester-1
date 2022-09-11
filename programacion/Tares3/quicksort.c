#define get_address(start, position, size) (start + position*size)
int cmpint(void *a, void *b){
    return (*((int*)a))-(*((int*)b));
}
void swapint(void *a, void *b){
    int *aa = (int*)a;
    int *bb = (int*)b;
    int c = *aa;
    *aa = *bb;
    *bb = c;
}
void quicksort( void *data, int start, int end, unsigned int data_size, int (*cmp)(void *a, void *b), void (*swap)(void *a, void *b)){
    if(end - start <= 1) return;
    int left = start;
    int right = end-1;
    int mid = right;
    while(1){
        while( cmp( data +  left * data_size , data + mid * data_size) < 0) left++;
        while( right > 0 && cmp( data +  right * data_size , data + mid * data_size) > 0 ) right--;
        if( left >= right)break;
        swap(data +  left * data_size,data +  right * data_size);
        if( right == mid) mid = left;
        if( left == mid) mid = right;
    }
    quicksort(data, start, left,data_size,cmp,swap);
    quicksort(data, left, end,data_size,cmp,swap);
}
