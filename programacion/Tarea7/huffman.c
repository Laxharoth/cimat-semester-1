#include "huffman.h"
#include <stdio.h>
#include <string.h>

byte_position *get_sorted_positions(unsigned int *abcedary);
Tree *gen_tree(byte_position *positions, unsigned int *leafs);
void fill_table(TreeNode *node, byte_code **current, char *code,
                unsigned int bit);
unsigned int search_in_table_by_value(const byte_code *table, const byte value);
unsigned int search_in_table_by_code(const byte_code *table, const char *code,
                                     unsigned int table_size);
int compress(const byte *file_buffer, unsigned long size,
             const char *output_file) {
  unsigned int leafs, index;
  unsigned long bit = 0, write_bytes;
  unsigned int *abcedary = handle_my_file(file_buffer, size);
  byte_position *positions = get_sorted_positions(abcedary);
  Tree *tree = gen_tree(positions, &leafs);

  unsigned int current_size = 100;
  unsigned current_byte;
  unsigned char next_bit_position;
  byte pushing_bits_size;
  char *current_code;
  byte *output_buffer = malloc(100);
  byte_code table[leafs];
  byte_code *current = table;
  char code[leafs];
  for (size_t i = 0; i < 0; i++) {
    memset(table[i].code, 0, leafs + 1);
  }
  fill_table(tree->root, &current, code, 0);
  for (size_t i = 0; i < size; i++) {
    index = search_in_table_by_value(table, file_buffer[i]);
    current_code = table[index].code;
    pushing_bits_size = strlen(current_code);
    for (unsigned int j = 0; j < pushing_bits_size; ++j) {
      current_byte = bit / 8;
      if (current_byte >= current_size) {
        output_buffer = realloc(output_buffer, current_size <<= 1);
      }
      next_bit_position = 7 - (bit % 8);
      // make sure ther is a 1 in that position
      if (current_code[j] == '1')
        output_buffer[current_byte] |= ((unsigned char)0b00000001)
                                       << next_bit_position;
      // make sure ther is a 0 in that position
      if (current_code[j] == '0')
        output_buffer[current_byte] &=
            (~(((unsigned char)0b00000001) << next_bit_position));
      bit++;
    }
  }
  FILE *f = fopen(output_file, "wb");
  write_bytes = bit / 8;
  if (bit % 8 != 0)
    ++write_bytes;
  /* HEADER */
  // original size
  fwrite((byte *)(&size), sizeof(unsigned long), 1, f);
  // write number of bits
  fwrite((byte *)(&bit), sizeof(unsigned long), 1, f);
  // write number of entries
  fwrite((byte *)(&leafs), sizeof(unsigned int), 1, f);
  /*  BODY  */
  // Write compresed file
  fwrite(output_buffer, write_bytes, 1, f);
  // write table
  for (unsigned int i = 0; i < leafs; ++i) {
    pushing_bits_size = strlen(table[i].code);
    // write size of code
    fwrite(&pushing_bits_size, sizeof(byte), 1, f);
    // write code
    fwrite(table[i].code, 1, pushing_bits_size, f);
    // write correspondent byte
    fwrite(&(table[i].value), 1, 1, f);
  }
#ifdef PRINT_TABLE
  printf("byte,codigo,frecuencia,tamano esperado\n");
  for (unsigned int i = 0; i < 255; i++) {
    if (abcedary[i] == 0)
      continue;
    index = search_in_table_by_value(table, i);
    printf("%d,\"%s\",%d,%ld\n", i, table[index].code, abcedary[i],
           strlen(table[index].code) * abcedary[i]);
  }
#endif // PRINT_TABLE

  // cleanup
  fclose(f);
  free(output_buffer);
  for (size_t i = 0; i < 0; i++)
    free(table[i].code);
  delete_Tree(tree);
  free(positions);
  free(abcedary);
  return 0;
}
int decompress(const byte *file_buffer, const char *output_file) {
  unsigned int table_entries;
  unsigned long size, total_bits;
  const byte *walkptr = file_buffer, *table_start;
  // Read header
  memcpy(&size, walkptr, sizeof(unsigned long));
  walkptr += sizeof(unsigned long);
  memcpy(&total_bits, walkptr, sizeof(unsigned long));
  walkptr += sizeof(unsigned long);
  memcpy(&table_entries, walkptr, sizeof(unsigned int));
  walkptr += sizeof(unsigned int);
  // Get table
  unsigned int bytes_in_body = ceil(total_bits / 8.0);
  table_start = walkptr + bytes_in_body;
  byte_code table[table_entries];
  for (unsigned int i = 0; i < table_entries; ++i, ++table_start) {
    // read representation size
    byte entrie_size = *table_start;
    ++table_start;
    // read representation
    for (byte j = 0; j < entrie_size; ++j, ++table_start) {
      table[i].code[j] = *table_start;
    }
    // end of representation string
    table[i].code[entrie_size] = '\0';
    // get original value
    table[i].value = *table_start;
  }
  char current_code[255];
  byte current_code_size = 0;
  // reconstruct file
  byte *original_file = malloc(sizeof(byte) * size);
  unsigned long current_byte_compresed_file;
  unsigned long current_byte_original_file = 0;
  byte bit_position;
  int table_index = -1;
  for (unsigned long current_bit = 0; current_bit < total_bits; ++current_bit) {
    current_byte_compresed_file = current_bit / 8;
    bit_position = 0b00000001 << (7 - (current_bit % 8));
    if (bit_position & walkptr[current_byte_compresed_file])
      current_code[current_code_size++] = '1';
    else
      current_code[current_code_size++] = '0';
    current_code[current_code_size] = '\0';
    if ((table_index = search_in_table_by_code(table, current_code,
                                               table_entries)) >= 0) {
      original_file[current_byte_original_file++] = table[table_index].value;
      current_code_size = 0;
    }
  }
  FILE *f = fopen(output_file, "wb");
  fwrite(original_file, size, 1, f);
  fclose(f);
  free(original_file);
  return 0;
}
unsigned int *handle_my_file(const byte *buffer, unsigned long size) {
  unsigned int *abcedary = calloc(BITS_COMB, sizeof(unsigned int));
  for (unsigned int i = 0; i < size; i++) {
    abcedary[buffer[i]]++;
  }
  return abcedary;
}
int cmp_byte_count(void *a, void *b) {
  byte_position *aa = (byte_position *)a;
  byte_position *bb = (byte_position *)b;
  return aa->count - bb->count;
}
void swap_byte_count(void *a, void *b) {
  byte_position c;
  byte_position *aa = (byte_position *)a;
  byte_position *bb = (byte_position *)b;
  c = *aa;
  *aa = *bb;
  *bb = c;
}
byte_position *get_sorted_positions(unsigned int *abcedary) {
  byte_position *positions = malloc(sizeof(byte_position) * BITS_COMB);
  for (unsigned int i = 0; i < BITS_COMB; i++) {
    positions[i].count = abcedary[i];
    positions[i].position = i;
  }
  quicksort(positions, 0, BITS_COMB, sizeof(byte_position), cmp_byte_count,
            swap_byte_count);
  return positions;
}
Tree *gen_tree(byte_position *positions, unsigned int *leafs) {
  *leafs = BITS_COMB;
  Tree *tree = new_tree();
  unsigned int start = 0;
  // skip bytes with no repetitions
  while (positions[start].count == 0) {
    start++;
    (*leafs)--;
  };
  TreeNode *lhs = new_TreeNode(positions[start]);
  byte_position parent_data;
  TreeNode *rhs;
  TreeNode *parent;
  for (unsigned int i = start + 1; i < BITS_COMB; ++i) {
    rhs = new_TreeNode(positions[i]);
    parent = new_TreeNode(parent_data);
    parent_data.count = rhs->data.count + lhs->data.count;
    if (lhs->data.count < rhs->data.count) {
      parent->left = lhs;
      parent->right = rhs;
    } else {
      parent->left = rhs;
      parent->right = lhs;
    }
    lhs = parent;
  }
  tree->root = parent;
  return tree;
}

void fill_table(TreeNode *node, byte_code **current, char *code,
                unsigned int bit) {
  if (node == NULL)
    return;
  if (node->left == NULL && node->right == NULL) {
    (**current).value = node->data.position;
    memcpy((**current).code, code, bit);
    (**current).code[bit] = 0;
    ++(*current);
    return;
  }
  code[bit] = '0';
  fill_table(node->left, current, code, bit + 1);
  code[bit] = '1';
  fill_table(node->right, current, code, bit + 1);
}
unsigned int search_in_table_by_value(const byte_code *table, const byte val) {
  unsigned int index = 0;
  while (val != table[index].value)
    ++index;
  return index;
}
unsigned int search_in_table_by_code(const byte_code *table, const char *code,
                                     unsigned int table_size) {
  for (size_t i = 0; i < table_size; i++) {
    if (strcmp(table[i].code, code) == 0)
      return i;
  }
  return -1;
}
