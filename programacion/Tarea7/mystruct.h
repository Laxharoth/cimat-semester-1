#ifndef STRUCT_H
#define STRUCT_H
typedef unsigned char byte;
struct byte_position_cout{
    unsigned int count;
    byte position;
};
struct byte_code{
    char code[255];
    byte value;
};
typedef struct byte_position_cout byte_position;
typedef struct byte_code byte_code;

#endif /* STRUCT_H */
