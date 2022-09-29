#ifndef QUICKSORT_H
#define QUICKSORT_H
#define get_address(start, position, size) (start + position*size)
int cmpint(void *a, void *b);
void swapint(void *a, void *b);
void quicksort( void *data, int start, int end, unsigned int data_size, int (*cmp)(void *a, void *b), void (*swap)(void *a, void *b));

#endif /* QUICKSORT_H */
