#ifndef FILE_MANIPS_H
#define FILE_MANIPS_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"



void saveVecToFileDetailedLongWays(fvec_t* vec, char* filename){
    FILE* file = fopen(filename, "wb");
    if(file == NULL){
        fprintf(stderr, "Error oppening file %s in code file %s  line %d\n", filename, __FILE__, __LINE__);
        return;
    }
    fprintf(file, "%d\n",vec->len);
    for(int i = 0; i < vec->len; i++){
        fprintf(file, "%lf\n", vec->vals[i]); // TODO easily swithcable prints 
    }
    fclose(file);
    return;
}

void saveVecToFile(fvec_t* vec, char* filename){
    FILE* file = fopen(filename, "wb");
    if(file == NULL){
        fprintf(stderr, "Error oppening file %s in code file %s  line %d\n", filename, __FILE__, __LINE__);
        return;
    }
    // fprintf(file, "%d\n",vec->len);
    for(int i = 0; i < vec->len; i++){
        fprintf(file, "%lf ", vec->vals[i]); // TODO easily swithcable prints 
    }
    fclose(file);
    return;
}

#endif // FILE_MANIPS_H