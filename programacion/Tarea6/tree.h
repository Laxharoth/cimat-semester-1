#ifndef TREE_H
#define TREE_H
#include "queue.h"
#include <stdlib.h>

struct _TreeNode{
    int data;
    struct _TreeNode *left;
    struct _TreeNode *right;
};
struct _Tree{
    struct _TreeNode *root;
    unsigned int nodes;
};
typedef struct _TreeNode TreeNode;
typedef struct _Tree Tree;

TreeNode *new_TreeNode(int data);
void delete_TreeNode(TreeNode *node);
Tree *new_tree();
void delete_Tree(Tree *tree);
void preorder( void (*callback)(TreeNode *node, void *arg) , TreeNode *curr,void *arg);
void inorder( void (*callback)(TreeNode *node, void *arg) , TreeNode *curr,void *arg);
void postorder( void (*callback)(TreeNode *node, void *arg) , TreeNode *curr,void *arg);
void levelorder( void (*callback)(TreeNode *node,void *arg), void (*lvlcallback)(void *arg) , Tree *tree, void *arg, void *lvlarg);
void clearTree(Tree *t);
void clear(TreeNode *node);
bool tree_insert(Tree* tree,int value);
bool tree_remove(Tree* tree,const int value);

#endif