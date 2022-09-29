#ifndef HUFFMAN_H
#define HUFFMAN_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mystruct.h"
#include "quicksort.h"
#include "tree.h"

#define BITS_COMB 256

int cmp_byte_count( void* a, void* b);
void swap_byte_count( void* a, void* b );
/**
 * @brief Returns the occurence of each byte in the file
 * 
 * @param filename 
 * @return unsigned int* 
 */
unsigned int *handle_my_file( const byte* filename, unsigned long size);
int compress(const byte* file_buffer, unsigned long size,const char *output_file);
int decompress(const byte* file_buffer,const char *output_file);
#endif /* HUFFMAN_H */

/*********************************
 *  .compressed file description
 *********************************
 *-         Header              -
 *------------------------------- 
 *  8 Bytes  | original file size ( unsigned long )
 * - - - - - - - - - - - - - - - - 
 *  8 Bytes  | number of bits of the compressed file ( unsigned long )
 * - - - - - - - - - - - - - - - -
 *  4 Bytes  | number of entries in table ( unsigned int )
 *********************************
 *-          Body               -
 *-------------------------------
 * ceil(n/8) | Compressed file
 *     bytes |
 * - - - - - - - - - - - - - - - -
 * variable  | Codes table
 * \                             /
 *  \                           /
 *   1 Byte  | Number of characters
 *           |  used in code ('0','1') combination
 *   - - - - - - - - - - - -
 *   n Byte  | Characters in code
 *   - - - - - - - - - - -
 *   1 Byte  | Original representing Byte
 ********************************/