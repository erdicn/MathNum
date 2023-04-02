#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "vectors.h"

#define EPS 3e-12


typedef double (*functype)();
typedef double (*functype1)(double);
typedef double (*functype2)(double, double);


typedef struct Tableau{
    int len;
    double* vals;
} tab_t;

typedef struct Func_Vals{
    int len;
    double* x;
    double* f;
}func_vals_t;

typedef struct Func_Vals2{
    int lenx;
    int leny;
    double* x;
    double* y;
    double** f;
}func_vals2_t;

typedef struct Vec_Func_Vals{
    int len;
    double* t;
    vec2D_t* vec;
}vec_func_vals_t;

typedef struct Value_Error{
    double val;
    double err;
}val_err;

void FreeFuncVals(func_vals_t* func);
void FreeFuncVals2(func_vals2_t* func);
void FreeVecFuncVals(vec_func_vals_t* func);

void FreeTab(tab_t* tab);
void PrintTab(tab_t tab);
tab_t* InitTab(int len_tab);
func_vals_t* InitFuncVals(int len);
func_vals_t* InitFuncValsWithInterval(int nb_intervals, double min, double max);
func_vals2_t* InitFuncVals2(int lenx, int leny);
vec_func_vals_t* InitVecFuncVals(int len);
vec_func_vals_t* InitVecFuncValsWithInterval(int nb_intervals, double min, double max);
tab_t* ReadTab(int len_tab, char* nom_fichier);
func_vals_t* ReadFuncAndValsXandF_InRows(char* nom_fichier);
func_vals_t* ReadFuncAndValsXandF_InColumns(char* nom_fichier);
void PrintFunc(func_vals_t func);
void PrintVecFuncVals(vec_func_vals_t func);
void EcritureFunc_Vals(func_vals_t func_to_print, char* nom_fichier_de_print);
void EcritureVecFuncVals(vec_func_vals_t func_to_print, char* nom_fichier_de_print);

#endif // FUNCTIONS_H