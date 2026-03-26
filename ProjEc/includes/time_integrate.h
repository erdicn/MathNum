#ifndef TIME_INTEGRATE_H
#define TIME_INTEGRATE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"

typedef double (*functypeEDOsolver)(double, double, functype2, double);

double EulersMethodGeneral(double yi, double dt, functypegeneral func_to_solve, void* params);
double EulersMethod(double xi, double yi, functype2 func_to_solve, double h);
double RK2HeunMethod(double xi, double yi, functype2 func_to_solve, double h);
double RK2EulerMethod(double xi, double yi, functype2 func_to_solve, double h);
double RK4Method(double xi, double yi, functype2 func_to_solve, double h);
void SolFuncVal(func_vals_t* fv, functype2 func_to_solve, double cond_init, char* solver_type);
void FillFuncValWithFunction(func_vals_t* f_a_replir, functype1 func_analythique);
double CalcMoyErr(func_vals_t f1, func_vals_t f2);

#endif // TIME_INTEGRATE_H