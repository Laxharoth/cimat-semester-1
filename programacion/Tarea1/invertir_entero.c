#include <stdio.h>

#define scan_int_prompt(prompt, variable)                                      \
  printf(prompt);                                                              \
  scanf("%d", &variable);

int reverse_int(int a) {
  int reversed = 0;
  while (a) {
    reversed = reversed * 10 + a - ((a / 10) * 10);
    a /= 10;
  }
  return reversed;
}
int main(int argc, char const *argv[]) {
  int num;
  scan_int_prompt("Introduce Numero: ", num);
  printf("Numero invertido:%d\n", reverse_int(num));
}
