#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "functions.h"

double Calcul_Li(int i, func_vals_t func, double x);
double Calcul_Pn(double x, func_vals_t func);
func_vals_t* Calcul_P_PointParPointAvecM(func_vals_t func, int m);
func_vals_t* Calcul_P_PointParPointAvecH(func_vals_t func, double h);


#endif // INTERPOLATION_H