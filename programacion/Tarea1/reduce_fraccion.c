#include <stdio.h>

int reduce_fraccion(int *_a, int *_b);

int main(int argc, char const *argv[]) {
  int a, b;
  printf("numerador:");
  scanf("%d", &a);
  printf("denominador:");
  scanf("%d", &b);
  reduce_fraccion(&a, &b);
  printf("%d/%d\n", a, b);
}

int reduce_fraccion(int *_a, int *_b) {
  int a = *_a;
  int b = *_b;
  if (a == 0)
    b = 1;
  if (b == 0)
    return -1;
  int min = a < b ? a : b;
  for (int i = min; i > 1; i--) {
    if ((a / i) * i == a && (b / i) * i == b) {
      a /= i;
      b /= i;
      i = a < b ? a : b;
    }
  }
  *_a = a;
  *_b = b;
  return 0;
}
