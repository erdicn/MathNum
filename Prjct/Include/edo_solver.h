#ifndef EDO_SOLVER_H
#define EDO_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"
#include "vectors.h"

typedef double (*functypeEDOsolver)(double, double, functype2, double);

double FuncTP2(double x, double y);
double SolAnalytiqueFuncTP2(double x);
double EulersMethod(double xi, double yi, functype2 func_to_solve, double h);
double RK2HeunMethod(double xi, double yi, functype2 func_to_solve, double h);
double RK2EulerMethod(double xi, double yi, functype2 func_to_solve, double h);
double RK4Method(double xi, double yi, functype2 func_to_solve, double h);
void SolFuncVal(func_vals_t* fv, functype2 func_to_solve, double cond_init, char* solver_type);
void FillFuncValWithFunction(func_vals_t* f_a_replir, functype1 func_analythique);
double CalcMoyErrFuncVals(func_vals_t f1, func_vals_t f2);

typedef vec2D_t (*functypeBIGvec)(vec2D_t, tab_t*);
typedef vec2D_t (*functypeEDOsolverBIGvec)(vec2D_t, tab_t*, functypeBIGvec, double);


vec2D_t VecEulersMethod(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h);
vec2D_t VecRK2HeunMethod(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h);
vec2D_t VecRK2EulersMethod(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h);
vec2D_t VecRK4Method(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h);
void SolVecFuncVal(vec_func_vals_t* vf, functypeBIGvec func_to_solve,
                     tab_t* vars_needed, vec2D_t cond_init, functypeEDOsolverBIGvec solver);

typedef giant_vec_t* (*functypeGIANTvec)(giant_vec_t*, tab_t*);
typedef giant_vec_t* (*functypeEDOsolverGIANTvec)(giant_vec_t*, tab_t*, functypeGIANTvec, double);

giant_vec_t* GiantEulersMethod(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h);
giant_vec_t* GiantRK2EulerMethod(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h);
giant_vec_t* GiantRK2HeunMethod(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h);
giant_vec_t* GiantRK4Method(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h);
void SolGiantVecFuncVal(giant_vec_func_vals_t* vf, functypeGIANTvec func_to_solve,
                     tab_t* vars_needed, giant_vec_t cond_init, functypeEDOsolverGIANTvec solver);
func_vals_t* GvecFuncTuFuncVals(int gvec_index, giant_vec_func_vals_t gvec_func_to_equalise);

#endif // EDO_SOLVER_H