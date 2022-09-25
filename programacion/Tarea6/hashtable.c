#include "hashtable.h"
#include "queue.h"

hashtable *new_hashtable(unsigned int buckets, unsigned int max_collisions){
    hashtable *t = malloc(sizeof(hashtable));
    t->size = buckets;
    t->max_collisions = max_collisions;
    t->buckets = malloc(sizeof(Tree*)*t->size);
    for (size_t i = 0; i < t->size; i++){
        t->buckets[i] = malloc(sizeof(Tree));
    }
    return t;
}
void delete_hashtable(hashtable *t){
    for (size_t i = 0; i < t->size; i++){
        delete_Tree(t->buckets[i]);
    }
    free(t->buckets);
    free(t);
}
unsigned int hash(hashtable* t,int key){
    return key%t->size;
}
bool hash_insert(hashtable *t,int value){
    unsigned int index = hash(t,value);
    if(t->buckets[index]->nodes < t->max_collisions)
        return tree_insert(t->buckets[index],value);
    rehash(t);
    hash_insert(t,value);
}
bool hash_remove(hashtable *t,int key){
    unsigned int index = hash(t,key);
    return tree_remove(t->buckets[index],key);
}
void hash_reinsert(TreeNode *node, void *arg){
    hashtable *t = (hashtable*)arg;
    hash_insert(t,node->data);
};
void store_val(TreeNode *node, void *arg){
    unsigned int **arr = arg;
    **arr = node->data;
    ++(*arr);
}
void rehash(hashtable *t){
    //store values in queue
    unsigned int total_n_val=0;
    for (size_t i = 0; i < t->size; i++){
        total_n_val += t->buckets[i]->nodes;
    }
    unsigned int *arr = malloc(sizeof(unsigned int)*total_n_val);
    unsigned int *cur_pos = arr;
    for(size_t i = 0; i < t->size; i++){
        preorder(store_val,t->buckets[i]->root,&cur_pos);
        clearTree(t->buckets[i]);
    }    
    //resize
    t->size*=2;
    t->buckets = realloc(t->buckets,t->size*sizeof(Tree*));
    for (size_t i = t->size/2; i < t->size; i++){
        t->buckets[i] = new_tree();
    }
    //reinsert
    for (size_t i = 0; i < total_n_val; i++){
        hash_insert(t, arr[i]);
    }
    free(arr);
}