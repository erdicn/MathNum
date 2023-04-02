#include <stdio.h>
#include <stdlib.h>
#include "../Include/edo_solver.h"
#include "../Include/functions.h"


int main(){
    double lower_bound = 0, upper_bound = 1; 
    int n_max = 10;
    double cond_init = 0;
    func_vals_t* fv = InitFuncValsWithInterval(n_max, lower_bound, upper_bound);
    func_vals_t* sol_analythique =  InitFuncValsWithInterval(n_max, lower_bound, upper_bound);
    
    FillFuncValWithFunction(sol_analythique, &SolAnalytiqueFuncTP2);
    EcritureFunc_Vals(*sol_analythique, "sol_analythique");
    SolFuncVal(fv, &FuncTP2, cond_init, "Euler");
    EcritureFunc_Vals(*fv, "xy_euler");
    SolFuncVal(fv, &FuncTP2, cond_init, "RK4");
    EcritureFunc_Vals(*fv, "xy_rk4");

    FreeFuncVals(fv);
    FreeFuncVals(sol_analythique);
    
    int len_h = 8;
    func_vals_t* h_Eh = InitFuncVals(len_h);
    char* names[] = {"err_euler", "err_rk2heun", "err_rk2euler", "err_rk4" };
    char sol_type[] = "solver type";
    for(int i = 0; i < 4; i++){
        int c = 0;
        for(n_max = 10; n_max <= pow(10,(len_h))+1; n_max *= 10){
            fv = InitFuncValsWithInterval(n_max, lower_bound, upper_bound);
            sol_type[0] = names[i][4];
            sol_type[2] = names[i][6];
            sol_type[3] = names[i][7];
            SolFuncVal(fv, &FuncTP2, cond_init, sol_type);
            // sol_analythique =  InitFuncValsWithInterval(n_max, lower_bound, upper_bound);
            // FillFuncValWithFunction(sol_analythique, &SolAnalytiqueFuncTP2);
            h_Eh->x[c] = fv->x[1] - fv->x[0];                          
            h_Eh->f[c] = fabs(fv->f[fv->len-1] - SolAnalytiqueFuncTP2(1.));//CalcMoyErr(*fv, *sol_analythique);
            c++;
            FreeFuncVals(fv);
            // FreeFuncVals(sol_analythique);
            printf("Finished %s %d\n", names[i], n_max);
        }
        EcritureFunc_Vals(*h_Eh, names[i]);
    }
    return 0;
}