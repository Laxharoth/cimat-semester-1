#include "tree.h"

TreeNode *new_TreeNode(int data) {
  TreeNode *node = malloc(sizeof(TreeNode));
  node->data = data;
  node->left = NULL;
  node->right = NULL;
}
void delete_TreeNode(TreeNode *node) { free(node); }
Tree *new_tree() {
  Tree *tree = malloc(sizeof(Tree));
  tree->root = NULL;
  tree->nodes = 0;
  return tree;
}
void delete_Tree(Tree *tree) {
  clearTree(tree);
  free(tree);
}
void preorder(void (*callback)(TreeNode *node, void *arg), TreeNode *curr,
              void *arg) {
  if (curr == NULL)
    return;
  callback(curr, arg);
  preorder(callback, curr->left, arg);
  preorder(callback, curr->right, arg);
}
void inorder(void (*callback)(TreeNode *node, void *arg), TreeNode *curr,
             void *arg) {
  if (curr == NULL)
    return;
  inorder(callback, curr->left, arg);
  callback(curr, arg);
  inorder(callback, curr->right, arg);
}
void postorder(void (*callback)(TreeNode *node, void *arg), TreeNode *curr,
               void *arg) {
  if (curr == NULL)
    return;
  postorder(callback, curr->left, arg);
  postorder(callback, curr->right, arg);
  callback(curr, arg);
}
void levelorder(void (*callback)(TreeNode *node, void *arg),
                void (*lvlcallback)(void *arg), Tree *tree, void *arg,
                void *lvlarg) {
  queue *nodes = new_queue(tree->nodes);
  queue *next = new_queue(tree->nodes);
  queue *swap;
  push(nodes, (void *)tree->root);
  while (!isEmpty(nodes)) {
    while (!isEmpty(nodes)) {
      TreeNode *curr = (TreeNode *)pop(nodes);
      if (curr == NULL)
        continue;
      callback(curr, arg);
      push(next, (void *)curr->left);
      push(next, (void *)curr->right);
    }
    if (lvlcallback != NULL)
      lvlcallback(lvlarg);
    swap = nodes;
    nodes = next;
    next = swap;
    clear_queue(next);
  }
  delete_queue(nodes);
  delete_queue(next);
}
void clearTree(Tree *t) {
  clear(t->root);
  t->root = NULL;
  t->nodes = 0;
}
void clear(TreeNode *node) {
  if (node == NULL)
    return;
  clear(node->left);
  clear(node->right);
  delete_TreeNode(node);
}
bool tree_insert(Tree *tree, int value) {
  if (tree->root == NULL) {
    tree->root = new_TreeNode(value);
    (tree->nodes)++;
    return true;
  }
  TreeNode *node = tree->root;
  while ((value < node->data && node->left != NULL) ||
         (value >= node->data && node->right != NULL)) {
    if (value < node->data)
      node = node->left;
    else
      node = node->right;
  }
  if (value == node->data)
    return false;
  if (value < node->data)
    node->left = new_TreeNode(value);
  else
    node->right = new_TreeNode(value);
  (tree->nodes)++;
  return true;
}
bool tree_remove(Tree *tree, const int value) {
  if (tree->root == NULL)
    return false;
  TreeNode *node = tree->root;
  TreeNode *prev = NULL;
  while (node != NULL && node->data != value) {
    prev = node;
    if (value < node->data)
      node = node->left;
    else
      node = node->right;
  }
  if (node == NULL)
    return false;
  // si no tiene hijos solo desacopla el nodo
  if (node->left == NULL && node->right == NULL) {
    if (node == tree->root) {
      tree->root = NULL;
      delete_TreeNode(node);
      (tree->nodes)++;
      return true;
    }
    if (value < prev->data)
      prev->left = NULL;
    else
      prev->right = NULL;
    delete_TreeNode(node);
    (tree->nodes)++;
    return true;
  }
  TreeNode *leaf = NULL;
  if (node->left != NULL) {
    leaf = node->left;
    if (leaf->right == NULL) {
      if (value < prev->data)
        prev->left = leaf;
      else
        prev->right = leaf;
      leaf->right = node->right;
      delete_TreeNode(node);
      (tree->nodes)++;
      return true;
    }
    TreeNode *prevLeaf = NULL;
    while (leaf->right != NULL) {
      prevLeaf = leaf;
      leaf = leaf->right;
    }
    prevLeaf->right = NULL;
  }
  if (node->right != NULL) {
    leaf = node->right;
    if (leaf->left == NULL) {
      if (value < prev->data)
        prev->left = leaf;
      else
        prev->right = leaf;
      leaf->left = node->left;
      delete_TreeNode(node);
      (tree->nodes)++;
      return true;
    }
    TreeNode *prevLeaf = NULL;
    while (leaf->left != NULL) {
      prevLeaf = leaf;
      leaf = leaf->left;
    }
    prevLeaf->left = NULL;
  }
  if (value < prev->data)
    prev->left = leaf;
  else
    prev->right = leaf;
  leaf->left = node->left;
  leaf->right = node->right;
  delete_TreeNode(node);
  (tree->nodes)++;
  return true;
}
