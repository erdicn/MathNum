#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../Include/functions.h"

typedef double (*functypeEDOsolver)(double, double, functype2, double);


double FuncTP2(double x, double y){
    return x*x*cos(y);
}


double EulersMethod(double xi, double yi, functype2 func_to_solve, double h){
    double k1 = func_to_solve(xi, yi);
    double yi_plus_1 = yi + h*k1;
    return yi_plus_1;
}

double RK2HeunMethod(double xi, double yi, functype2 func_to_solve, double h){
    double k1 = func_to_solve(xi, yi);
    double k2 = func_to_solve(xi+h, yi+h*k1);
    double yi_plus_1 = yi + (h/2.)*(k1+k2);
    return yi_plus_1;
}

double RK2EulerMethod(double xi, double yi, functype2 func_to_solve, double h){
    double k1 = func_to_solve(xi, yi);
    double k2 = func_to_solve(xi+h/2., yi+(h/2.)*k1);
    double yi_plus_1 = yi + h*k2;
    return yi_plus_1;
}

double RK4Method(double xi, double yi, functype2 func_to_solve, double h){
    double k1 = func_to_solve(xi, yi);
    double k2 = func_to_solve(xi+h/2., yi+(h/2.)*k1);
    double k3 = func_to_solve(xi+h/2., yi+(h/2.)*k2);
    double k4 = func_to_solve(xi+h, yi+h*k3);
    double yi_plus_1 = yi + (h/6.)*(k1 + 2*k2 + 2*k3 + k4);
    return yi_plus_1;
}

void SolFuncVal(func_vals_t* fv, functype2 func_to_solve, double cond_init, char* solver_type){
    if (fv->len < 2){
        printf("Function has less than 2 values allocated to it that doesnt work\n");
        exit(-1);
    }
    // First we choose our solver that we entered into the function
    functypeEDOsolver solver = NULL;
    if (solver_type[0] == 'e' || solver_type[0] == 'E'){
        solver = &EulersMethod;
    } else if (solver_type[0] == 'r' || solver_type[0] == 'R'){
        if (solver_type[2] == '2'){
            if (solver_type[3] == 'h' || solver_type[3] == 'H'){
                solver = &RK2HeunMethod;
            } else if (solver_type[3] == 'e' || solver_type[3] == 'E'){
                solver = &RK2EulerMethod;
            }
        } else if(solver_type[2] == '4'){
            solver = &RK4Method; // TODO here we need a & but i am testing 
        }
    } 
    
    if (solver == NULL){
        printf("Error the function doesnt found defaulted to RK4 since it is the most accurate (str entered: %s)\n", solver_type);
        solver = &RK4Method;
    }

    fv->f[0] = cond_init; // comme tout les methodes de solutions sont pdes methodes qui resolue le futur a tout les cas cest le 0eme peut etre ameliorer ca dans le futur  
    double h = fv->x[1] - fv->x[0]; // since our h doesnt change and we chose it acoringly to n_max in main this is faster to do than writing a h each time we call this function
    for(int i = 1; i < fv->len; i++){
        fv->f[i] = solver(fv->x[i-1], fv->f[i-1], func_to_solve, h);
    }
}

void FillFuncValWithFunction(func_vals_t* f_a_replir, functype1 func_analythique){
    for(int i = 0; i < f_a_replir->len; i++){
        f_a_replir->f[i] = func_analythique(f_a_replir->x[i]);
    }
}

double CalcMoyErr(func_vals_t f1, func_vals_t f2){
    double tot = 0;
    if (f1.len != f2.len){
        printf("Les longueurs des valeurs de fonctions ne sont pas egales !\n");
        exit(1);
    }
    for(int i = 0; i < f1.len; i++){
        tot += fabs(f1.f[i] - f2.f[i]);
    }
    return tot/ f1.len; //fabs(f1.f[f1.len-1] - f2.f[f1.len-1]);
}


int main(){
    printf("Result: %lf %lf\n", RK4Method(0, 0, &FuncOldExam, 1./2.), 57./16.);
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

//fun stuff to visualise function adresses ant ot compare them how var adresses are visualised 
//int* a = malloc(sizeof(int));
//printf("%p\n%p\n%p\n%p\n", RK4Method, &RK4Method, a, &a);