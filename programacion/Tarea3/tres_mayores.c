#include <stdio.h>

typedef char byte;

int tres_mayores(void *data, unsigned int *const primero,
                 unsigned int *const segundo, unsigned int *const tercero,
                 const unsigned int data_size, const unsigned int array_size,
                 int (*cmp)(void *a, void *b));
int cmp_u_int(void *a, void *b);
int cmp_byte(void *a, void *b);

int main(int argc, char const *argv[]) {
  unsigned int data[] = {
      7679, 9632, 2992, 7523, 6159, 9540, 1724, 822,  9164, 4945, 3100, 5977,
      3069, 9832, 7816, 955,  7694, 4352, 590,  6417, 843,  8977, 983,  6778,
      9157, 6423, 9413, 2777, 6702, 7109, 2775, 3081, 7085, 6060, 3330, 5110,
      994,  2051, 2584, 3949, 7420, 3586, 8461, 8061, 7369, 7466, 3672, 4813,
      6337, 3507, 3963, 5051, 3545, 2938, 1592, 6712, 4140, 8707, 4247, 2563,
      7199, 5149, 529,  1249, 3210, 1199, 2148, 2985, 133,  7617, 66,   48,
      3607, 8740, 4393, 5951, 9738, 9256, 7826, 1007, 31,   9922, 8681, 6506,
      8558, 1515, 604,  1665, 7916, 6300, 6016, 2747, 8147, 9395, 187,  9081,
      5904, 5782, 8839, 6105};
  int primero, segundo, tercero;
  tres_mayores(data, &primero, &segundo, &tercero, sizeof(unsigned int), 100,
               cmp_u_int);
  printf("unsigned int\n");
  printf("primero pos: %d segundo pos: %d tercero pos: %d\n", primero, segundo,
         tercero);
  printf("primero: %d segundo: %d tercero: %d\n", data[primero], data[segundo],
         data[tercero]);
  printf("byte\n");
  tres_mayores(data, &primero, &segundo, &tercero, sizeof(byte),
               100 * sizeof(unsigned int), cmp_byte);
  printf("primero pos: %d segundo pos: %d tercero pos: %d\n", primero, segundo,
         tercero);
  printf("primero: %d segundo: %d tercero: %d\n", ((byte *)data)[primero],
         ((byte *)data)[segundo], ((byte *)data)[tercero]);
  return 0;
}
int tres_mayores(void *data, unsigned int *const primero,
                 unsigned int *const segundo, unsigned int *const tercero,
                 const unsigned int data_size, const unsigned int array_size,
                 int (*cmp)(void *a, void *b)) {
  *primero = *segundo = *tercero = 0;
  // ordena los tres primeros
  if (cmp(data, data + data_size)) {
    *segundo = *primero;
    *primero = 1;
  } else
    *segundo = 1;
  if (cmp(data + data_size * (*primero), data + data_size * 2)) {
    *tercero = *segundo;
    *segundo = *primero;
    *primero = 2;
  } else if (cmp(data + data_size * (*segundo), data + data_size * 2)) {
    *tercero = *segundo;
    *segundo = 2;
  } else
    *tercero = 2;
  for (unsigned int i = 3 * data_size; i < array_size * data_size;
       i += data_size) {
    if (cmp(data + ((*primero) * data_size), data + i)) {
      *tercero = *segundo;
      *segundo = *primero;
      *primero = i / data_size;
    }
  }
}
int cmp_u_int(void *a, void *b) {
  unsigned int *na = (unsigned int *)a;
  unsigned int *nb = (unsigned int *)b;
  return *na < *nb;
}
int cmp_byte(void *a, void *b) {
  byte *na = (byte *)a;
  byte *nb = (byte *)b;
  return *na < *nb;
}
