#include <stdio.h>
#include <stdlib.h>

#include "quicksort.c"

void resize(void **holder, unsigned int *current_size, unsigned int data_size);
// como mstrcmp y strswap se utilizan como apuntadores a funciones es mejor
// definirlas antes de usarlas
int mstrcmp(void *a, void *b) {
  const char **aa = (const char **)a;
  const char **bb = (const char **)b;
  unsigned int i = 0, j = 0;
  while (1) {
    if ((*aa)[i] != (*bb)[j] || (*aa)[i] == 0 || (*bb)[j] == 0)
      break;
    i++;
    j++;
  }
  return ((int)((*aa)[i])) - ((int)((*bb)[j]));
}
void strswap(void *a, void *b) {
  char **aa = (char **)a;
  char **bb = (char **)b;
  char *c = (*aa);
  (*aa) = (*bb);
  (*bb) = c;
}
void mstrcpy(const char *a, char *b);
int main(int argc, char const *argv[]) {
  unsigned int n_max_names = 0;
  unsigned int n_names = 0;
  char **names = NULL;
  FILE *fp = fopen("names.txt", "r");
  char *buffer = malloc(100);
  size_t len;
  ssize_t llen;
  while ((llen = getline(&buffer, &len, fp)) != EOF) {
    // remove new line character
    buffer[llen - 1] = 0;
    char *name = malloc(llen);
    mstrcpy(buffer, name);
    if (n_max_names == n_names)
      resize((void **)(&names), &n_max_names, sizeof(char *));

    names[n_names++] = name;
  }
  // release last allocation
  quicksort(names, 0, n_names, sizeof(char *), mstrcmp, strswap);
  for (unsigned int i = 0; i < n_names; i++) {
    printf("%s\n", names[i]);
    free(names[i]);
  }
  free(buffer);
  free(names);
  return 0;
}
void resize(void **holder, unsigned int *current_size, unsigned int data_size) {
  if (holder == NULL) {
    holder = malloc(sizeof(char *) * (++(*current_size)));
  }
  (++(*current_size));
  *holder = realloc(*holder, *current_size * data_size);
}
void mstrcpy(const char *a, char *b) {
  int i = 0;
  while (a[i] != 0)
    b[i++] = a[i];
}