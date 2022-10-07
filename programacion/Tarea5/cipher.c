#include "pgm1.c"
#include <stdlib.h>
#include <string.h>

void cipher(char *img, int max_size);
char *concat(const char *s1, const char *s2);

int main(int argc, char const *argv[]) {
  int rows, cols;
  int max_size;
  char *name = "barbara.ascii.pgm";
  char *endname;
  unsigned char *img = pgmRead(name, &rows, &cols);
  max_size = rows * cols;
  while (max_size % 4 != 0)
    max_size++;
  img = realloc(img, max_size);
  cipher(img, max_size);
  endname = concat("cipher.", name);
  pgmWrite(endname, rows, cols, img, "");
  free(endname);
  return 0;
}

void cipher(char *img, int max_size) {
  int *iimg = (int *)img;
  int i = 0, j = 0;
  unsigned int *maps = malloc(3 * sizeof(int));
  unsigned int *tmp_maps = malloc(3 * sizeof(int));
  unsigned int *swap;
  unsigned int c = 16, m = 8, e = 65535;
  maps[0] = 0xF00F00;
  maps[1] = 0x00F0F0;
  maps[2] = 0xF0000F;

  for (i = 0; i < max_size / 4;) {
    int H = e & (maps[0] ^ maps[1] ^ maps[2]);
    for (j = 0; j < 3 && i < max_size / 4; j++, i++) {
      iimg[i] = iimg[i] ^ maps[j];
      tmp_maps[j] = ((c * maps[j]) + (maps[j] >> m)) + H;
    }
    swap = maps;
    maps = tmp_maps;
    tmp_maps = swap;
  }
  free(maps);
  free(tmp_maps);
}
char *concat(const char *s1, const char *s2) {
  char *result = malloc(strlen(s1) + strlen(s2) + 1);
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}
