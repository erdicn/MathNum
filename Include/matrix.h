#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <math.h>

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

myfloat* accessFMatD(fmatD_t mat, int i, int j){
    return mat.vals + i * mat.nb_cols + j;
}

typedef struct fvec{
    int len;
    myfloat *vals;
} fvec_t;

fvec_t* allocateFVec(fvec_t* vec, int len){
    if (vec == NULL)
        vec = malloc(sizeof(fvec_t));
    vec->len  = len; 
    vec->vals = calloc(len, sizeof(myfloat));
    return vec;
}

void freeFVec(fvec_t* vec){
    if(vec) free(vec->vals);
    free(vec);
}

fmatCSR_t* allocateFMatCSR(fmatCSR_t* mat, int nb_non_zero, int nb_rows, int nb_cols){
    if( mat == NULL){
        mat = malloc(sizeof(fmatCSR_t));
    }  
    mat->nb_cols     = nb_cols;
    mat->nb_rows     = nb_rows;
    mat->nb_non_zero = nb_non_zero;

    mat->vals  = calloc(nb_non_zero       , sizeof(myfloat));
    mat->col_i = calloc(nb_non_zero       , sizeof(int)    );
    mat->row_i = calloc((mat->nb_rows + 1), sizeof(int)    );
    return mat;
}

void freeFMatCSR(fmatCSR_t* mat){
    if( mat == NULL){
        return ;
    }  

    free(mat->vals );
    free(mat->col_i);
    free(mat->row_i);
    free(mat);
    return;
}

void freeFMatCOO(fmatCOO_t* mat){
    if( mat == NULL){
        return ;
    }  

    free(mat->vals );
    free(mat->col_i);
    free(mat->row_i);
    free(mat);
    return;
}

void freeFMatD(fmatD_t* mat){
    if( mat == NULL){
        return ;
    }  

    free(mat->vals);
    free(mat);
    return;
}

fmatCOO_t* allocateFMatCOO(fmatCOO_t* mat, int nb_non_zero, int nb_rows, int nb_cols){
    if( mat == NULL){
        mat = malloc(sizeof(fmatCOO_t));
    }  
    mat->nb_cols = nb_cols;
    mat->nb_rows = nb_rows;
    mat->nb_non_zero = nb_non_zero;

    mat->vals  = calloc(nb_non_zero, sizeof(myfloat));
    mat->col_i = calloc(nb_non_zero, sizeof(int)  );
    mat->row_i = calloc(nb_non_zero, sizeof(int)  );
    return mat;
}

fmatD_t* allocateFMatD(fmatD_t* mat, int nb_rows, int nb_cols){
    if( mat == NULL){
        mat = malloc(sizeof(fmatD_t));
    }  
    mat->nb_cols = nb_cols;
    mat->nb_rows = nb_rows;

    mat->vals  = calloc(nb_cols*nb_rows, sizeof(myfloat));
    return mat;
}

fmatCSR_t* initialiseFMatCSRFromCOO(fmatCSR_t* mat_writing, fmatCOO_t* mat_read){
    if (mat_writing == NULL){
        mat_writing = allocateFMatCSR(mat_writing, mat_read->nb_non_zero, mat_read->nb_rows, mat_read->nb_cols);
    }
    // TODO if mat is not null but not correct dimension is not checked
    return NULL; // TODO not finished yet
}

int getRandomInt(int min, int max){
    return min + rand() / (RAND_MAX / (max - min + 1) + 1);
}

myfloat getRandomFloat(myfloat min, myfloat max){
    myfloat random_val_scale = rand()/(myfloat)RAND_MAX;
    return min + random_val_scale * (max - min);
}

// TODO not sure if it is the optimal comparator
int comp(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

void printIArrValues(int* vals, int len);
void printFArrValues(myfloat* vals, int len);

// fills a random matrix with density uses rand() so initialize it accordingly
fmatCSR_t* initializeRandomFMatCSR(int nb_rows, int nb_cols, myfloat density, int rand_seed, myfloat rand_min, myfloat rand_max){
    int nb_non_zero = (int) (nb_cols*nb_rows*density);
    fmatCSR_t* mat = allocateFMatCSR(NULL, nb_non_zero, nb_rows, nb_cols);
    srand(rand_seed); // TODO not sure if this is the best practice
    // bool unique = false;
    // bool zero_found = false; // since the callocated array is zeros we give it chance to put someting there
    int i;// j;
    int total_elements = nb_cols*nb_rows;
    char* used_indices = calloc(total_elements,  sizeof(char)); // bitmap goes fast like brrrrr
    for (i = 0; i < mat->nb_non_zero; i++){
        mat->vals[i]  = getRandomFloat(rand_min, rand_max);
        int candidate;
        do {
            candidate = getRandomInt(0, total_elements);
        } while(used_indices[candidate] == 1);
        
        used_indices[candidate] = 1; // Mark as used
        mat->col_i[i] = candidate;
        // // // // verry slow 
        // get random indice till unique TODO could take a lot of tif really dense
        // // // while(!unique){
        // // //     mat->col_i[i] = getRandomInt(0, nb_rows*nb_cols);
        // // //     for(j = 0; j < i; j++){
        // // //         if(mat->col_i[j] == mat->col_i[i]) break;
        // // //     } if(j == i) unique = true;
        // // // }
        // // // unique = false;
    }
    free(used_indices);

    // the faith of the indices and values gets separated
    qsort(mat->col_i, mat->nb_non_zero, sizeof(int), comp);
    // now that we have our indices sorted we need to convest them to our way of storage
    int col_i, row_i;
    mat->row_i[0] = 0;
    for(i = 0; i < mat->nb_non_zero; i++){
        row_i = mat->col_i[i] / mat->nb_cols; 
        col_i = mat->col_i[i] % mat->nb_cols;
        mat->col_i[i] = col_i;
        mat->row_i[row_i+1]++;
    }

    for(i = 0; i < mat->nb_rows; i++){
        mat->row_i[i + 1] += mat->row_i[i];
    }

    return mat;
}

fmatD_t* initaliseRandomFMatD(fmatD_t* mat, int rand_seed, myfloat rand_min, myfloat rand_max){
    assert(mat);
    assert(mat->vals);
    for(size_t i = 0; i < mat->nb_cols*mat->nb_rows; i++){
        mat->vals[i]  = getRandomFloat(rand_min, rand_max);
    }
    return mat;
}

fvec_t* initializeRandomFVec(int len, myfloat rand_min, myfloat rand_max){
    fvec_t* vec = allocateFVec(NULL, len);
     for (int i = 0; i < vec->len; i++)
        vec->vals[i]  = getRandomFloat(rand_min, rand_max);
    return vec;
}

fvec_t* linspaceFvec(fvec_t* vec, myfloat min, myfloat max, int len){
    if (vec == NULL){
        vec = allocateFVec(vec, len);
    }
    myfloat dx = (max-min)/len;
    for(int i = 0; i < len; i++){
        vec->vals[i] = min + i*dx;
    }
    return vec;
}

fvec_t* linspaceFvecCentered(fvec_t* vec, myfloat min, myfloat max, int len){
    if (vec == NULL){
        vec = allocateFVec(vec, len);
    }
    myfloat dx = (max-min)/len;
    for(int i = 0; i < len; i++){
        vec->vals[i] = min + dx/2 + i*dx;
    }
    return vec;
}

fvec_t* linspaceFvecIncludeMax(fvec_t* vec, myfloat min, myfloat max, int len){
    if (vec == NULL){
        vec = allocateFVec(vec, len);
    }
    myfloat dx = (max-min)/(len-1);
    for(int i = 0; i < len; i++){
        vec->vals[i] = min + i*dx ;
    }
    return vec;
}


fvec_t* initializeFVec(fvec_t* vec, myfloat init_val){
    if (vec == NULL){
        printf("Error vec is not allocated so cant initialise.");
        return NULL;
    }
    for(int i = 0; i < vec->len; i++){
        vec->vals[i] = init_val;
    }
    return vec;
}

fvec_t* newCopyFVec(fvec_t* new_vec, fvec_t* vec){
    new_vec = allocateFVec(new_vec, vec->len); // TODO memoy leak doesnt frees the new_vec data
    memcpy(new_vec->vals, vec->vals, sizeof(myfloat)*vec->len);
    return new_vec;
}

myfloat l2NormFVec(fvec_t vec){
    myfloat tot = 0;
    for(int i = 0; i < vec.len; i++){
        tot += vec.vals[i]*vec.vals[i];
    }
    return pow(tot, 0.5);
}

myfloat lnNormFVec(fvec_t vec, myfloat norm){
    myfloat tot = 0;
    for(int i = 0; i < vec.len; i++){
        tot += vec.vals[i]*vec.vals[i];
    }
    return pow(tot, 1./norm);
}

fmatCSR_t* initializeFMatCSRFromDenseMat(fmatD_t* D, fmatCSR_t* CSR, myfloat tolerance){
    if (tolerance < 0) tolerance = 1e-5;

    if( CSR == NULL){
        CSR = malloc(sizeof(fmatCSR_t));
    }

    CSR->nb_cols = D->nb_cols;
    CSR->nb_rows = D->nb_rows;
    
    int i, j, nb_non_zero = 0; 
    for( i = 0; i < D->nb_rows; i++){
        for( j = 0; j < D->nb_cols; j++){
            if (fabs(*accessFMatD(*D, i, j)) >= tolerance )
                nb_non_zero++;
        }   
    }

    printf("nb non 0=%d\n", nb_non_zero);
    allocateFMatCSR(CSR, nb_non_zero, D->nb_rows, D->nb_cols);
    CSR->row_i[0] = 0;

    int counter = 0, non_zero_in_row = 0;
    for( i = 0; i < D->nb_rows; i++){
        non_zero_in_row = 0;
        for( j = 0; j < D->nb_cols; j++){
            if (fabs(*accessFMatD(*D, i, j)) >= tolerance ){
                CSR->col_i[counter] = j;
                CSR->vals [counter] = *accessFMatD(*D, i, j);
                non_zero_in_row++;
                counter++;
            }
        }
        CSR->row_i[i+1] = non_zero_in_row + CSR->row_i[i] ;
    }
    return CSR;
}

void printFMatCSRDetails(fmatCSR_t mat){
    printf("Matrix is (%d %d) with %d nonzeros\n", mat.nb_rows, mat.nb_cols, mat.nb_non_zero);
}

void printFArrValues(myfloat* vals, int len){
    for(int i = 0; i < len; i++)
        printf("%.3f ", vals[i]);
    printf("\n");
}

void printFVec(fvec_t* vec){
    printFArrValues(vec->vals, vec->len);
}

void printIArrValues(int* vals, int len){
    for(int i = 0; i < len; i++)
        printf("%d ", vals[i]);
    printf("\n");
}

void printFMatCSR(fmatCSR_t mat){
    printFMatCSRDetails(mat);
    for (int i = 0; i < mat.nb_rows; i++) {
        int row_start = mat.row_i[i];
        int row_end   = mat.row_i[i+1];

        for (int j = 0; j < mat.nb_cols; j++) {
            myfloat val_to_print = 0.0f;
            for (int k = row_start; k < row_end; k++) {
                if (mat.col_i[k] == j) {
                    val_to_print = mat.vals[k];
                    break; 
                }
            }

            printf("%.3f ", val_to_print);
        }
        printf("\n");
    }
}

fvec_t* matVecProdCSRf(fmatCSR_t* mat, fvec_t* vec, fvec_t* sol_vec){
    if (mat->nb_cols != vec->len){
        printf("Mat (%d, %d), vec, %d dim not compatible\n", mat->nb_cols, mat->nb_rows, vec->len);
        return NULL;
    }
    if(sol_vec == NULL) sol_vec = allocateFVec(sol_vec, mat->nb_cols);
    if (sol_vec->len != mat->nb_cols){
        printf("Reallocating sol vec\n");
        sol_vec->len = mat->nb_cols;
        free(sol_vec->vals);
        sol_vec->vals = NULL;
    }
    if (sol_vec->vals == NULL){
        printf("Allocating sol vec\n");
        sol_vec->len = mat->nb_cols;
        sol_vec->vals = calloc( mat->nb_cols, sizeof(myfloat));
    }

    for (int i = 0; i < mat->nb_cols; i++){
        sol_vec->vals[i] = 0;
        for (int j = mat->row_i[i]; j < mat->row_i[i+1]; j++){
            sol_vec->vals[i] += mat->vals[j] * vec->vals[mat->col_i[j]];
        }
    }

    return sol_vec;
}

// fmatCSR_t* matMatProdFCSR(fmatCSR_t* mat1, fmatCSR_t* mat2, fmatCSR_t* mat_sol){
    
//     return mat_sol;
// }

// TODO
fvec_t* matTransposeVecProdCSRf(fmatCSR_t* mat, fvec_t* vec, fvec_t* sol_vec){
    if (mat->nb_rows != vec->len){
        printf("Mat transposed (%d, %d), vec, %d dim not compatible\n", mat->nb_rows, mat->nb_cols, vec->len);
        return NULL;
    }
    if(sol_vec == NULL) sol_vec = allocateFVec(sol_vec, mat->nb_rows);
    if (sol_vec->len != mat->nb_rows){
        printf("Reallocating sol vec\n");
        sol_vec->len = mat->nb_rows;
        free(sol_vec->vals);
        sol_vec->vals = NULL;
    }
    if (sol_vec->vals == NULL){
        printf("Allocating sol vec\n");
        sol_vec->len = mat->nb_rows;
        sol_vec->vals = calloc( mat->nb_rows, sizeof(myfloat));
    }
 
    // TODO change indexing for transpose
    for (int i = 0; i < mat->nb_cols; i++){
        sol_vec->vals[i] = 0;
        for (int j = mat->row_i[i]; j < mat->row_i[i+1]; j++){
            sol_vec->vals[i] += mat->vals[j] * vec->vals[mat->col_i[j]];
        }
    }

    return sol_vec;
}

fmatCSR_t* matScalarProdCSRf(fmatCSR_t* mat, myfloat f){
    for(int i = 0; i < mat->nb_non_zero; i++){
        mat->vals[i] *= f;
    }
    return mat;
}

fvec_t* vecScalarProdVecf(fvec_t* vec, myfloat f){
    for(int i = 0; i < vec->len; i++){
        vec->vals[i] *= f;
    }
    return vec;
}

fvec_t* vecAddVecf(fvec_t* vec1, fvec_t* vec2, fvec_t* vec_sol){
    if(vec_sol == NULL) vec_sol = allocateFVec(vec_sol, vec1->len);
    // TODO add more catches 
    for(int i = 0; i < vec1->len; i++){
        vec_sol->vals[i] = vec1->vals[i] + vec2->vals[i];
    }
    return vec_sol;
}

fvec_t* vecAddS(fvec_t* vec1, myfloat s, fvec_t* vec_sol){
    if(vec_sol == NULL) vec_sol = allocateFVec(vec_sol, vec1->len);
    // TODO add more catches 
    for(int i = 0; i < vec1->len; i++){
        vec_sol->vals[i] = vec1->vals[i] + s;
    }
    return vec_sol;
}

void freeCRSf(fmatCSR_t* mat){
    if(mat == NULL) return;
    free(mat->vals);
    free(mat->vals);
    free(mat->row_i);
    return;
}

fmatCSR_t* newCopyCSRf(fmatCSR_t* new_mat, fmatCSR_t* mat){
    freeCRSf(new_mat);
    new_mat = allocateFMatCSR(new_mat, mat->nb_non_zero, mat->nb_rows, mat->nb_cols);    
    memcpy(new_mat->vals , mat->vals , (mat->nb_non_zero)*sizeof(myfloat));
    memcpy(new_mat->row_i, mat->row_i, (mat->nb_rows  +1)*sizeof(int)  );
    memcpy(new_mat->col_i, mat->col_i, (mat->nb_non_zero)*sizeof(int)  );
    return new_mat;
}

fmatCSR_t* newMatScalarProdCSRf(fmatCSR_t* new_mat, fmatCSR_t* mat, myfloat f){
    new_mat = newCopyCSRf(new_mat, mat);
    for(int i = 0; i < new_mat->nb_non_zero; i++){
        new_mat->vals[i] *= f;
    }
    return new_mat;
}

myfloat getMaxOfVec(fvec_t* vec){
    myfloat max = 0;
    for(int i = 0; i  < vec->len; i++){
        max = max < vec->vals[i] ? vec->vals[i] : max;
    }
    return max;
}

myfloat getMinOfVec(fvec_t* vec){
    myfloat min = 0;
    for(int i = 0; i  < vec->len; i++){
        min = min > vec->vals[i] ? vec->vals[i] : min;
    }
    return min;
}

#endif /* MATRIX_H */
