#include <stdio.h>
#include <stdlib.h>
#include "vectors.h"

vec2D_t VecAdd(vec2D_t v1, vec2D_t v2){
    vec2D_t new_vec = {v1.x+v2.x, v1.y+v2.y};
    return new_vec;
}

vec2D_t VecxScal(vec2D_t v1, double s){
    vec2D_t new_vec = {v1.x*s, v1.y*s};
    return new_vec;
}

giant_vec_func_vals_t* InitGvecFuncVals(int len){
    giant_vec_func_vals_t* new_gvec_func = malloc(sizeof(giant_vec_func_vals_t));
    new_gvec_func->len = len;
    new_gvec_func->t = malloc(sizeof(double)*new_gvec_func->len);
    new_gvec_func->gvecs = malloc(sizeof(giant_vec_t*)*new_gvec_func->len);
    return new_gvec_func;
}

giant_vec_func_vals_t* InitGvecFuncValsInterval(int nb_intervals, double min, double max){
    int len = nb_intervals + 1;
    giant_vec_func_vals_t* new_gvec_f = InitGvecFuncVals(len);
    double h = (max - min) / nb_intervals;
    for(int i = 0; i < new_gvec_f->len; i++){
        new_gvec_f->t[i] = min + h*i;
    }
    return new_gvec_f;
}

giant_vec_t* InitGVec(int len){
    giant_vec_t* new_vec = malloc(sizeof(giant_vec_t));
    new_vec->len = len;
    new_vec->comp = malloc(sizeof(double)*new_vec->len);
    return new_vec;
}

void PrintGVec(giant_vec_t vec){
    for(int i = 0; i < vec.len; i++){
        printf("%.16lf ", vec.comp[i]);
    } puts("");
}

void PrintGiantVecFunc(giant_vec_func_vals_t gvfv){
    for(int i = 0; i < gvfv.len; i++){
        printf("%lf ", gvfv.t[i]);
        PrintGVec(*(gvfv.gvecs[i]));
    }
}

void FreeGVec(giant_vec_t* vec){
    free(vec->comp);
    free(vec);
}

void FreeGvecFuncVals(giant_vec_func_vals_t* gvec_func){
    for(int i = 0; i < gvec_func->len; i++){
        FreeGVec(gvec_func->gvecs[i]);  
    }
    free(gvec_func->t);
}

giant_vec_t* GVecAdd(giant_vec_t vec1, giant_vec_t vec2){
    if(vec1.len != vec2.len){
        PrintGVec(vec1);
        PrintGVec(vec2);
        printf("Vectors not the same length\n");
        exit(-11);
    }
    giant_vec_t* new_vec = InitGVec(vec1.len);
    for(int i = 0; i < vec1.len; i++){
        new_vec->comp[i] = vec1.comp[i] + vec2.comp[i];
    }
    return new_vec;
}

giant_vec_t* GVecxScal(giant_vec_t vec, double s){
    giant_vec_t* new_vec = InitGVec(vec.len);
    for(int i = 0; i < vec.len; i++){
        new_vec->comp[i] = vec.comp[i]*s;
    }
    return new_vec;
}

// Makes so vec1 is equal to vec2
void GVecEq(giant_vec_t* vec1, giant_vec_t vec2){
    if(vec1->len != vec2.len){
        printf("Error the two vectors are not the same length %d %d\n", vec1->len, vec2.len);
        exit(-33);
    }
    for(int i = 0; i < vec1->len; i++){
        vec1->comp[i] = vec2.comp[i];
    }
}

void EcritureGiantVecFuncVals(giant_vec_func_vals_t func_to_print, char* nom_fichier_de_print){
    FILE* fichier = fopen(nom_fichier_de_print, "w");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier_de_print);
    }
    for(int i = 0; i < func_to_print.len; i++){
        fprintf(fichier, "%.16lf ", func_to_print.t[i]);
        for(int j = 0; j < func_to_print.gvecs[i]->len; j++){
            fprintf(fichier, "%.16lf ", func_to_print.gvecs[i]->comp[j]);
        }
        fprintf(fichier, "\n");
    }
    fclose(fichier);
}