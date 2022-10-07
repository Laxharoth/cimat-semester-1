#include <stdio.h>

typedef char byte;
int cantidad_a_pagar(int pesos, int *veinte, int *diez, int *cinco, int *uno);

int main(int argc, char const *argv[]) {
  int pesos = 0;
  int veinte = 0, diez = 0, cinco = 0, uno = 0;
  printf("Enter a dollar amount: ");
  scanf("%d", &pesos);
  cantidad_a_pagar(pesos, &veinte, &diez, &cinco, &uno);
  printf("$20 bills: %d\n", veinte);
  printf("$10 bills: %d\n", diez);
  printf("$5 bills: %d\n", cinco);
  printf("$1 bills: %d\n", uno);
  return 0;
}

int cantidad_a_pagar(int pesos, int *veinte, int *diez, int *cinco, int *uno) {
  int bills[] = {20, 10, 5, 1};
  int *counter[] = {veinte, diez, cinco, uno};
  (*veinte) = (*diez) = (*cinco) = (*uno) = 0;
  for (byte i = 0; i < 4; i++) {
    while (pesos - bills[i] >= 0) {
      pesos -= bills[i];
      (*(counter[i]))++;
    }
  }
  return 0;
}