#ifndef HASHTABLE_H
#define HASHTABLE_H
#include "tree.h"
#include <stdbool.h>

struct _hashtable{
    unsigned int size;
    unsigned int max_collisions;
    Tree **buckets;
};
typedef struct _hashtable hashtable;
hashtable *new_hashtable(unsigned int buckets, unsigned int max_collisions);
void delete_hashtable(hashtable *t);
unsigned int hash(hashtable* t,int key);
bool hash_insert(hashtable *t,int value);
bool hash_remove(hashtable *t,int key);
void rehash(hashtable *t);

#endif /* HASHTABLE_H */
