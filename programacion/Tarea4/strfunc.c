#include <stdlib.h>

typedef char *string;

int start_with(const string str1, const string str2);
int *strfind(const string str, const string sub_str, int *frec);
string *split(const string str, char chr);
string mconcat(const string str1, const string str2);
int mstrlen(const string str);

int main(int argc, char const *argv[]) {
  /* code */
  return 0;
}

/**
 * @brief Encuntra el las coincidencias del substring en el string
 *
 * @param str El string original
 * @param sub_str El string a buscar
 * @param frec El numero de coincidencias
 * @return int*
 */
int *strfind(const string str, const string sub_str, int *frec) {
  int *res = NULL;
  *frec = 0;
  int index = 0;
  while (str[index] != '\0') {
    if (!start_with(str + index, sub_str)) {
      index++;
      continue;
    }
    ++(*frec);
    if (res == NULL)
      res = malloc(sizeof(int));
    else
      res = realloc(res, sizeof(int) * (*frec));
    res[(*frec) - 1] = index;
    index++;
  }
  return res;
}

/**
 * @brief Divide un string en los substring entre el caracter dado
 *
 * @param str
 * @param chr
 * @return string*
 */
string *split(const string str, char chr) {
  string *res = malloc(sizeof(string *));
  int start = 0;
  int found = 1;
  int index = 0;
  int cpy_index;
  while (str[index] != '\0') {
    if (str[index] != chr) {
      index++;
      continue;
    }
    cpy_index = start;
    res[found - 1] = malloc(index - start + 1);
    while (cpy_index < index) {
      res[found - 1][cpy_index] = str[cpy_index];
      ++cpy_index;
    }
    res[found - 1][cpy_index] = '\0';
    found++;
    index++;
  }
  cpy_index = start;
  res[found - 1] = malloc(index - start + 1);
  while (cpy_index < index) {
    res[found - 1][cpy_index] = str[cpy_index];
    ++cpy_index;
  }
  res[found - 1][cpy_index] = '\0';
  return res;
}

/**
 * @brief Une dos string en un tercero
 *
 * @param str1
 * @param str2
 * @return string
 */
string mconcat(string str1, const string str2) {
  int len1 = mstrlen(str1);
  int len2 = mstrlen(str2);
  str1 = realloc(str1, len1 + len2 + 1);
  str1[len1 + len2] = '\0';
  for (size_t i = 0; i < len2; i++) {
    str1[len1 + i] = str2[i];
  }
  return str1;
}

int start_with(const string str1, const string str2) {
  int index = -1;
  while (str2[++index] != '\0')
    if (str1[index] != str2[index])
      return 0;
  return 1;
}
int mstrlen(const string str) {
  int len = 0;
  while (str[len] != '\0')
    ++len;
  return len;
}
