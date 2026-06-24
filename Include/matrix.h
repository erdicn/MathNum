#ifndef MATRIX_H
#define MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "myfloat.h"

typedef struct fmatCOO{
    int nb_cols, nb_rows;
    int nb_non_zero;
    int *row_i, *col_i;
    myfloat *vals;
} fmatCOO_t;


// CSR format in each ro it doesnt matter if it is sorter or not as long as it is in the correct row
// EXAM might be asked for MPI and openMP but not exact code but psuedocode 
typedef struct fmatCSR{
    int nb_cols, nb_rows, nb_non_zero;
    int *row_i, *col_i;// row_i is the indice in which the vals start for that row
    myfloat *vals;
} fmatCSR_t;

typedef struct fmatD{
    int nb_rows, nb_cols;
    myfloat *vals; 
} fmatD_t;

myfloat* accessFMatD(fmatD_t mat, int i, int j);

typedef struct fvec{
    int len;
    myfloat *vals;
} fvec_t;

typedef struct ivec{
    int len;
    int *vals;
} ivec_t;

fvec_t* allocateFVec(fvec_t* vec, int len);
fvec_t* reallocateFVec(fvec_t* vec, int new_len);
void freeFVec(fvec_t* vec);
ivec_t* allocateIVec(ivec_t* vec, int len);
ivec_t* reallocateIVec(ivec_t* vec, int new_len);
void freeIVec(ivec_t* vec);
fmatCSR_t* allocateFMatCSR(fmatCSR_t* mat, int nb_non_zero, int nb_rows, int nb_cols);
void freeFMatCSR(fmatCSR_t* mat);
void freeFMatCOO(fmatCOO_t* mat);
void freeFMatD(fmatD_t* mat);
fmatCOO_t* allocateFMatCOO(fmatCOO_t* mat, int nb_non_zero, int nb_rows, int nb_cols);
fmatD_t* allocateFMatD(fmatD_t* mat, int nb_rows, int nb_cols);
fmatCSR_t* initialiseFMatCSRFromCOO(fmatCSR_t* mat_writing, fmatCOO_t* mat_read);
int getRandomInt(int min, int max);
myfloat getRandomFloat(myfloat min, myfloat max);
// TODO not sure if it is the optimal comparator
int comp(const void* a, const void* b);
void printIArrValues(int* vals, int len);
void printFArrValues(myfloat* vals, int len);

// fills a random matrix with density uses rand() so initialize it accordingly
fmatCSR_t* initializeRandomFMatCSR(int nb_rows, int nb_cols, myfloat density, int rand_seed, myfloat rand_min, myfloat rand_max);
fmatD_t* initaliseRandomFMatD(fmatD_t* mat, int rand_seed, myfloat rand_min, myfloat rand_max);
fvec_t* initializeRandomFVec(int len, myfloat rand_min, myfloat rand_max);
fvec_t* linspaceFvec(fvec_t* vec, myfloat min, myfloat max, int len);
fvec_t* linspaceFvecCentered(fvec_t* vec, myfloat min, myfloat max, int len);
// also includes min too 
fvec_t* linspaceFvecIncludeMax(fvec_t* vec, myfloat min, myfloat max, int len);

fvec_t* initializeFVec(fvec_t* vec, myfloat init_val);
fvec_t* newCopyFVec(fvec_t* new_vec, fvec_t* vec);
myfloat l2NormFVec(fvec_t vec);
myfloat lnNormFVec(fvec_t vec, myfloat norm);
fmatCSR_t* initializeFMatCSRFromDenseMat(fmatD_t* D, fmatCSR_t* CSR, myfloat tolerance);
void printFMatCSRDetails(fmatCSR_t mat);
void printFArrValues(myfloat* vals, int len);
void printFVec(fvec_t* vec);
void printIArrValues(int* vals, int len);
void printFMatCSR(fmatCSR_t mat);
fvec_t* matVecProdCSRf(fmatCSR_t* mat, fvec_t* vec, fvec_t* sol_vec);
// fmatCSR_t* matMatProdFCSR(fmatCSR_t* mat1, fmatCSR_t* mat2, fmatCSR_t* mat_sol){
    
//     return mat_sol;
// }

// TODO
fvec_t* matTransposeVecProdCSRf(fmatCSR_t* mat, fvec_t* vec, fvec_t* sol_vec);
fmatCSR_t* matScalarProdCSRf(fmatCSR_t* mat, myfloat f);
fvec_t* vecScalarProdVecf(fvec_t* vec, myfloat f);
fvec_t* vecAddVecf(fvec_t* vec1, fvec_t* vec2, fvec_t* vec_sol);
fvec_t* vecAddS(fvec_t* vec1, myfloat s, fvec_t* vec_sol);
void freeCRSf(fmatCSR_t* mat);
fmatCSR_t* newCopyCSRf(fmatCSR_t* new_mat, fmatCSR_t* mat);
fmatCSR_t* newMatScalarProdCSRf(fmatCSR_t* new_mat, fmatCSR_t* mat, myfloat f);
myfloat getMaxOfVec(fvec_t* vec);
myfloat getMinOfVec(fvec_t* vec);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H */
