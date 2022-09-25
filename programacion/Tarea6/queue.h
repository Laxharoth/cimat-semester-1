/*********************************
 * 
 * Queue implementation heavily based on tutorialspoint:
 * https://www.tutorialspoint.com/data_structures_algorithms/queue_program_in_c.htm
 * 
*********************************/

#ifndef QUEUE_H
#define QUEUE_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

struct _queue{
    void **queue_array;
    unsigned int max;
    int front;
    int rear;
    int itemCount;
};
typedef struct _queue queue;
queue *new_queue(unsigned int max);
void delete_queue(queue *q);
bool isEmpty(queue *q);
bool isFull(queue *q);
int size(queue *q);
void push(queue *q,void *data);
void *pop(queue *q);
void clear_queue(queue *q);

#endif /* QUEUE_H */
