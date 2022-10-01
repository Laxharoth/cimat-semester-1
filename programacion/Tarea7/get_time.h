#include <bits/types/struct_timeval.h>
#include <sys/time.h>
#include <stdlib.h>

#define GetTimeStamp(process)({\
    long elapsed;\
    struct timeval begin,end;\
    gettimeofday(&begin,NULL);\
    process;\
    gettimeofday(&end,NULL);\
    elapsed = (end.tv_sec-begin.tv_sec)*1000000 + end.tv_usec-begin.tv_usec;\
    elapsed;\
})