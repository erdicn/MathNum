#ifndef VECTORS_H
#define VECTORS_H

// 2D vectors
typedef struct Vec2D{
    double x;
    double y;
} vec2D_t;

vec2D_t VecAdd(vec2D_t v1, vec2D_t v2);
vec2D_t VecxScal(vec2D_t v1, double s);

// Giant vectors 
typedef struct GiantVec{
    int len;
    double* comp;
}giant_vec_t;

typedef struct GVec_Func_Vals{
    int len;
    double* t;
    giant_vec_t** gvecs;
}giant_vec_func_vals_t;

giant_vec_func_vals_t* InitGvecFuncVals(int len);
giant_vec_func_vals_t* InitGvecFuncValsInterval(int nb_intervals, double min, double max);
giant_vec_t* InitGVec(int len);
void PrintGVec(giant_vec_t vec);
void PrintGiantVecFunc(giant_vec_func_vals_t gvfv);
void FreeGVec(giant_vec_t* vec);
void FreeGvecFuncVals(giant_vec_func_vals_t* gvec_func);
giant_vec_t* GVecAdd(giant_vec_t vec1, giant_vec_t vec2);
giant_vec_t* GVecxScal(giant_vec_t vec, double s);
void GVecEq(giant_vec_t* vec1, giant_vec_t vec2);
void EcritureGiantVecFuncVals(giant_vec_func_vals_t func_to_print, char* nom_fichier_de_print);

#endif //VECTORS_H