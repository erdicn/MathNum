#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"

// precondition les x de func sont ecarte equalement x1-x0 = x2-x1
val_err CompositeTrapezes(double lower_bound, double upper_bound, func_vals_t func){
    if (func.len < 2){
        printf("La fonction rentre a moins que 2 points besoin au moins de 2 points\n");
    }
    double h = func.x[1] - func.x[0];
    int counter = 0;
    double first = 0, last = 0, tot = 0;
    for(int i = 0; i < func.len; i ++){
        if (func.x[i] >= lower_bound - EPS && func.x[i] <= upper_bound + EPS){ // pour pouvoir integrer dans des intervals
            if (counter == 0) {first = func.f[i]; counter++;}
            tot += func.f[i];
            last = func.f[i];
        } else if (func.x[i] > upper_bound){break;}
    }
    val_err result;
    result.val = (tot-(first/2.0)-(last/2.0))*h; // on enleve le premier et le dernier car oa besoinde leur compteur une seule fois et dans la fonction on les multiplie d=par deux pour faciliter les autres 
    result.err = h*h;
    return result;
}

val_err CompositeSimposon(double lower_bound, double upper_bound, func_vals_t func){
    if (func.len < 2 && func.len % 3 != 0){
        printf("La fonction rentre a moins que 3 ou le nb de points nest pas divisable par 3\n");
    }
    double h = func.x[1] - func.x[0];
    //printf("%.16lf %.16lf\n",h, (func.x[func.len-1] - func.x[0])/(double)m);
    double last = 0, tot = 0;
    double last_multiplier = 1;
    // printf("\n\n%d\n", func.len);
    // double Ih = 0.;
    for(int i = 0; i < func.len; i++){
        if (func.x[i] >= lower_bound - EPS && func.x[i] <= upper_bound + EPS){ //-EPS et +EPS
            if (i == 0) {tot += func.f[i];}
            else {
                if ( i % 2 == 0){
                    tot += 2*func.f[i];
                    last_multiplier = 2;  
                } else {
                    tot += 4*func.f[i]; 
                    last_multiplier = 4; 
                }
                last = func.f[i];
            }
        } else if (func.x[i] > upper_bound){break;}
        // if(i %2 == 0 && i < func.len -1){
        //     Ih += (func.f[i] + 4*func.f[i+1] + func.f[i+2]) * (h/3.);
        // } 
        // printf("%.16lf %.16lf\n", Ih, (h/3.0)*(tot - (last_multiplier-1)*last));
    }

    val_err result;
    result.val = (h/3.0)*(tot - (last_multiplier-1)*last);
    result.err = h*h*h*h;

    return result;
}