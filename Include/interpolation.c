#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"

double Calcul_Li(int i, func_vals_t func, double x){
    double tot = 1;
    for(int j = 0; j < func.len; j++){
        if(j == i) continue;
        tot *= (x - func.x[j]) / (func.x[i] - func.x[j]);
    }
    return tot;
}

double Calcul_Pn(double x, func_vals_t func){
    double tot = 0;
    for(int i = 0; i < func.len; i++){
        tot += func.f[i]*Calcul_Li(i, func, x);
    }
    return tot;
}

// calculs les valeurs du polynome de langrange point par point das un intervale [lower_lim, upper_lim] avec m sous intervales avec lower_lim la 1er valeur de la fonctionet upper lim la plus grande valeur (en x)
func_vals_t* Calcul_P_PointParPointAvecM(func_vals_t func, int m){
    double h = (func.x[func.len-1] - func.x[0])/(double)m;
    //printf("%lf \n", h);
    double pas_despace;
    func_vals_t* vecP = InitFuncVals(m+1);
    for(int k = 0; k <= m; k++){
        pas_despace = func.x[0] + k * h;
        vecP->x[k] = pas_despace;
        vecP->f[k] = Calcul_Pn(vecP->x[k], func);
    } 
    return vecP;
}

func_vals_t* Calcul_P_PointParPointAvecH(func_vals_t func, double h){
    int m = (int)((double)(func.x[func.len-1] - func.x[0])/h);
    printf("%d \n", m);
    double pas_despace;
    func_vals_t* vecP = InitFuncVals(m+1);
    for(int k = 0; k <= m; k++){
        pas_despace = func.x[0] + k * h;
        vecP->x[k] = pas_despace;
        vecP->f[k] = Calcul_Pn(vecP->x[k], func);
    } 
    return vecP;
}