#include <stdio.h>
#include <stdlib.h>

#define scan_int_prompt(prompt, variable) printf(prompt);scanf("%d",&variable);

int change_base(const int num,const int base ,char *buffer);
int _change_base(const int num,const int base ,char *buffer, int *index);
int main(int argc, char const *argv[]){
    int num, base;
    scan_int_prompt("Introduce numero: ", num);
    scan_int_prompt("Introduce base: ", base);

    char *buffer = calloc(1,15);
    change_base(num,base,buffer);
    printf("%s\n",buffer);
    return 0;
}

int _change_base(const int num,const int base ,char *buffer, int *index){
    if(base<2 || base>36)return -1;
    int division,residuo = num;
    if(num >= base){
        division = num/base;
        residuo = num - division * base;
        _change_base(division,base,buffer,index);
    }
    switch (residuo){
    case 0: buffer[(*index)++] = '0';break;
    case 1: buffer[(*index)++] = '1';break;
    case 2: buffer[(*index)++] = '2';break;
    case 3: buffer[(*index)++] = '3';break;
    case 4: buffer[(*index)++] = '4';break;
    case 5: buffer[(*index)++] = '5';break;
    case 6: buffer[(*index)++] = '6';break;
    case 7: buffer[(*index)++] = '7';break;
    case 8: buffer[(*index)++] = '8';break;
    case 9: buffer[(*index)++] = '9';break;
    case 10: buffer[(*index)++] = 'A';break;
    case 11: buffer[(*index)++] = 'B';break;
    case 12: buffer[(*index)++] = 'C';break;
    case 13: buffer[(*index)++] = 'D';break;
    case 14: buffer[(*index)++] = 'E';break;
    case 15: buffer[(*index)++] = 'F';break;
    case 16: buffer[(*index)++] = 'G';break;
    case 17: buffer[(*index)++] = 'H';break;
    case 18: buffer[(*index)++] = 'I';break;
    case 19: buffer[(*index)++] = 'J';break;
    case 20: buffer[(*index)++] = 'K';break;
    case 21: buffer[(*index)++] = 'L';break;
    case 22: buffer[(*index)++] = 'M';break;
    case 23: buffer[(*index)++] = 'N';break;
    case 24: buffer[(*index)++] = 'O';break;
    case 25: buffer[(*index)++] = 'P';break;
    case 26: buffer[(*index)++] = 'Q';break;
    case 27: buffer[(*index)++] = 'R';break;
    case 28: buffer[(*index)++] = 'S';break;
    case 29: buffer[(*index)++] = 'T';break;
    case 30: buffer[(*index)++] = 'U';break;
    case 31: buffer[(*index)++] = 'V';break;
    case 32: buffer[(*index)++] = 'W';break;
    case 33: buffer[(*index)++] = 'X';break;
    case 34: buffer[(*index)++] = 'Y';break;
    case 35: buffer[(*index)++] = 'Z';break;
    }
}
int change_base(int num,int base ,char *buffer){
    int index = 0;
    _change_base(num,base ,buffer, &index);
}