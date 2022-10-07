#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "hashtable.h"
#include "tree.h"

#if defined(_WIN32) || defined(WIN32)
#define CLEAR system("cls");
#else
#define CLEAR system("clear");
#endif

// max_depth es 15 (numero de colisiones)
void main_loop(hashtable *t);
void get_depth(TreeNode *node, unsigned int depth, unsigned int *maxdepth);
void print_node(TreeNode *node, void *arg);
void println(void *arg);
unsigned int max_depth(hashtable *t);
int main(int argc, char const *argv[]) {
  srand(time(NULL));
  int r = rand();
  hashtable *t = new_hashtable(10, 15);
  for (size_t i = 0; i < 150; i++) {
    while (!hash_insert(t, r))
      r = rand();
  }
  main_loop(t);
  delete_hashtable(t);
  return 0;
}
void main_loop(hashtable *t) {
  unsigned char choose;
  int insert;
  bool inserted;
  unsigned int index, i;
  while (1) {
    CLEAR;
    printf("1.-Desplegar todos los árboles (first breath).\n");
    printf("2.-Desplegar todos los árboles (inorden).\n");
    printf("3.-Desplegar el índice y árbol con mayor profundidad.\n");
    printf("4.-Eliminar nodo.\n");
    printf("5.-Insertar Nodo.\n");
    printf("6.-Salir.\n");
    printf("choose:\n");
    choose = getchar();
    getchar();
    switch (choose) {
    case '1':
      for (i = 0; i < t->size; i++) {
        printf("indice %d:\n", i);
        levelorder(print_node, println, t->buckets[i], NULL, NULL);
      }
      break;
    case '2':
      for (i = 0; i < t->size; i++) {
        printf("indice %d:\n", i);
        inorder(print_node, t->buckets[i]->root, NULL);
        println(NULL);
      }
      break;
    case '3':
      index = max_depth(t);
      printf("indice %d:\n", index);
      inorder(print_node, t->buckets[index]->root, NULL);
      println(NULL);
      break;
    case '4':
      printf("Ingresar valor:");
      scanf("%d", &insert);
      getchar(); // clear buffer
      inserted = hash_remove(t, insert);
      if (inserted)
        printf("%d Eliminado\n", insert);
      else
        printf("%d No se encontro\n", insert);
      break;
    case '5':
      printf("Ingresar valor:");
      scanf("%d", &insert);
      getchar(); // clear buffer
      inserted = hash_insert(t, insert);
      if (inserted)
        printf("%d Ingresado\n", insert);
      else
        printf("%d Ya ingresado\n", insert);
      break;
    default:
      break;
    }
    if (choose == '6')
      break;
    printf("press enter:\n");
    getchar();
  }
}
void get_depth(TreeNode *node, unsigned int depth, unsigned int *maxdepth) {
  if (node == NULL) {
    if (depth > *maxdepth)
      *maxdepth = depth;
    return;
  }
  get_depth(node->left, depth + 1, maxdepth);
  get_depth(node->right, depth + 1, maxdepth);
}
unsigned int max_depth(hashtable *t) {
  unsigned int index = 0;
  unsigned int max_depth = 0;
  unsigned int depth = 0;
  for (unsigned int i = 0; i < t->size; i++) {
    get_depth(t->buckets[i]->root, 0, &depth);
    if (depth > max_depth) {
      index = i;
      max_depth = depth;
    }
  }
  return index;
}
void print_node(TreeNode *node, void *arg) { printf("%d, ", node->data); }
void println(void *arg) { printf("\n"); }