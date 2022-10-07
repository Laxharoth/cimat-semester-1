#include "queue.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

queue *new_queue(unsigned int max) {
  queue *q = malloc(sizeof(queue));
  q->queue_array = malloc(sizeof(void **) * max);
  q->max = max;
  q->front = 0;
  q->rear = -1;
  q->itemCount = 0;
  return q;
}
void delete_queue(queue *q) {
  free(q->queue_array);
  free(q);
}
bool isEmpty(queue *q) { return q->itemCount == 0; }
bool isFull(queue *q) { return q->itemCount == q->max; }
int size(queue *q) { return q->itemCount; }
void push(queue *q, void *data) {
  if (!isFull(q)) {
    if (q->rear == q->max - 1)
      q->rear = -1;
    q->queue_array[++(q->rear)] = data;
    (q->itemCount)++;
  }
}
void *pop(queue *q) {
  if (q->itemCount == 0)
    return NULL;
  void *data = (q->queue_array)[(q->front)++];
  if (q->front == q->max)
    q->front = 0;
  (q->itemCount)--;
  return data;
}
void clear_queue(queue *q) {
  q->front = 0;
  q->rear = -1;
  q->itemCount = 0;
}
