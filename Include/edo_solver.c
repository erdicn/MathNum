#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"
#include "vectors.h"
#include "edo_solver.h"


double FuncTP2(double x, double y){
    return x*x*cos(y);
}

double SolAnalytiqueFuncTP2(double x){
    return 2*atan(exp(x*x*x/3)) - M_PI_2;
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

double CalcMoyErrFuncVals(func_vals_t f1, func_vals_t f2){
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


// Vec 2D Solvers

vec2D_t VecEulersMethod(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h){
    vec2D_t k1 = func_to_solve(vec, vars_needed);
    vec2D_t vec_plus_1 = VecAdd(vec, VecxScal(k1, h));
    return vec_plus_1;
}

vec2D_t VecRK2HeunMethod(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h){
    vec2D_t k1 = func_to_solve(vec, vars_needed);
    vec2D_t k2 = func_to_solve(VecAdd(vec, VecxScal(k1, h)), vars_needed);
    vec2D_t vec_plus_1 = VecAdd(vec, VecxScal(VecAdd(k1, k2), h/2));
    return vec_plus_1;
}

vec2D_t VecRK2EulersMethod(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h){
    vec2D_t k1 = func_to_solve(vec, vars_needed);
    vec2D_t k2 = func_to_solve(VecAdd(vec, VecxScal(k1, h/2)), vars_needed);
    vec2D_t vec_plus_1 = VecAdd(vec, VecxScal(k2, h));
    return vec_plus_1;
}

vec2D_t VecRK4Method(vec2D_t vec, 
                        tab_t* vars_needed, functypeBIGvec func_to_solve, double h){
    vec2D_t k1 = func_to_solve(vec, vars_needed);
    vec2D_t k2 = func_to_solve(VecAdd(vec, VecxScal(k1, h/2)), vars_needed);
    vec2D_t k3 = func_to_solve(VecAdd(vec, VecxScal(k2, h/2)), vars_needed);
    vec2D_t k4 = func_to_solve(VecAdd(vec, VecxScal(k3, h)), vars_needed);
    vec2D_t vec_plus_1 = VecAdd(vec, VecxScal( VecAdd( VecAdd(k1, k4), VecxScal(VecAdd(k2, k3), 2)), h/6));
    return vec_plus_1;
}

void SolVecFuncVal(vec_func_vals_t* vf, functypeBIGvec func_to_solve,
                     tab_t* vars_needed, vec2D_t cond_init, functypeEDOsolverBIGvec solver){
    if (vf->len < 2){
        printf("Function has less than 2 values allocated to it that doesnt work\n");
        exit(-1);
    }
    vf->vec[0].x = cond_init.x; // comme tout les methodes de solutions sont pdes methodes qui resolue le futur a tout les cas cest le 0eme peut etre ameliorer ca dans le futur  
    vf->vec[0].y = cond_init.y;
    double h = vf->t[1] - vf->t[0]; // since our h doesnt change and we chose it acoringly to n_max in main this is faster to do than writing a h each time we call this function
    for(int i = 1; i < vf->len; i++){
        vf->vec[i] = solver(vf->vec[i-1], vars_needed, func_to_solve, h);
    }
}

// Giant Vec solvers 

giant_vec_t* GiantEulersMethod(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h){
    giant_vec_t* k1 = func_to_solve(vec, vars_needed);
    giant_vec_t* k1_times_h = GVecxScal(*k1, h);          // On fait ceci car pour ajouter les vectors de taille changable on a besoin de leur passer par adresse et si on le mets pas dans une variable on peux avoir des fuites de memoire
    giant_vec_t* vec_plus_1 = GVecAdd(*vec, *k1_times_h);  
    FreeGVec(k1_times_h);
    return vec_plus_1;
}

giant_vec_t* GiantRK2EulerMethod(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h){
    giant_vec_t* k1 = func_to_solve(vec, vars_needed);
    giant_vec_t* k1_times_h_2 = GVecxScal(*k1, h/2);                   // Pour pas avoir fuite de memoire   
    giant_vec_t* vec_plus_k1_times_h_2 = GVecAdd(*vec, *k1_times_h_2); // Pour pas avoir fuite de memoire
    giant_vec_t* k2 = func_to_solve(vec_plus_k1_times_h_2, vars_needed);
    giant_vec_t* k2_times_h = GVecxScal( *k2, h);                      // Pour pas avoir fuite de memoire
    giant_vec_t* vec_plus_1 = GVecAdd(*vec, *k2_times_h); 
    FreeGVec(k1_times_h_2);                                            // Pour pas avoir fuite de memoire
    FreeGVec(vec_plus_k1_times_h_2);                                   // Pour pas avoir fuite de memoire            
    FreeGVec(k2_times_h);                                              // Pour pas avoir fuite de memoire
    return vec_plus_1;
}

giant_vec_t* GiantRK2HeunMethod(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h){
    giant_vec_t* k1 = func_to_solve(vec, vars_needed);
    giant_vec_t* k1_times_h = GVecxScal(*k1, h);                       // Pour pas avoir fuite de memoire
    giant_vec_t* vec_plus_k1_times_h = GVecAdd(*vec, *k1_times_h);     // Pour pas avoir fuite de memoire
    giant_vec_t* k2 = func_to_solve(vec_plus_k1_times_h, vars_needed);
    giant_vec_t* k1_plus_k2 = GVecAdd(*k1, *k2);                       // Pour pas avoir fuite de memoire
    giant_vec_t* k1_plus_k2_times_h_2 = GVecxScal( *k1_plus_k2, h/2);  // Pour pas avoir fuite de memoire
    giant_vec_t* vec_plus_1 = GVecAdd(*vec, *k1_plus_k2_times_h_2); 
    FreeGVec(k1_times_h);                                              // Pour pas avoir fuite de memoire
    FreeGVec(vec_plus_k1_times_h);                                     // Pour pas avoir fuite de memoire
    FreeGVec(k1_plus_k2);                                              // Pour pas avoir fuite de memoire
    FreeGVec(k1_plus_k2_times_h_2);                                    // Pour pas avoir fuite de memoire
    return vec_plus_1;
}

giant_vec_t* GiantRK4Method(giant_vec_t* vec, 
                        tab_t* vars_needed, functypeGIANTvec func_to_solve, double h){
    giant_vec_t* k1 = func_to_solve(vec, vars_needed);
    giant_vec_t* k1_times_h_2 = GVecxScal(*k1, h/2);                                                       // Pour pas avoir fuite de memoire 
    giant_vec_t* vec_plus_k1_times_h_2 = GVecAdd(*vec, *k1_times_h_2);                                     // Pour pas avoir fuite de memoire
    giant_vec_t* k2 = func_to_solve(vec_plus_k1_times_h_2, vars_needed);                                   
    giant_vec_t* k2_times_h_2 = GVecxScal(*k2, h/2);                                                       // Pour pas avoir fuite de memoire
    giant_vec_t* vec_plus_k2_times_h_2 = GVecAdd(*vec, *k2_times_h_2);                                     // Pour pas avoir fuite de memoire
    giant_vec_t* k3 = func_to_solve(vec_plus_k2_times_h_2, vars_needed);                               
    giant_vec_t* k3_times_h = GVecxScal(*k3, h);                                                           // Pour pas avoir fuite de memoire
    giant_vec_t* vec_plus_k3_times_h = GVecAdd(*vec, *k3_times_h);                                         // Pour pas avoir fuite de memoire
    giant_vec_t* k4 = func_to_solve(vec_plus_k3_times_h, vars_needed);                                    
    giant_vec_t* k1_plus_k4 = GVecAdd(*k1, *k4);                                                           // Pour pas avoir fuite de memoire     
    giant_vec_t* k2_plus_k3 = GVecAdd(*k2,*k3);                                                            // Pour pas avoir fuite de memoire 
    giant_vec_t* k2_plus_k3_times_2 = GVecxScal(*k2_plus_k3, 2);                                           // Pour pas avoir fuite de memoire
    giant_vec_t* k1_plus_2k2_plus_2k3_plus_k4 = GVecAdd(*k1_plus_k4, *k2_plus_k3_times_2);                 // Pour pas avoir fuite de memoire
    giant_vec_t* k1_plus_2k2_plus_2k3_plus_k4__times_h_6 = GVecxScal(*k1_plus_2k2_plus_2k3_plus_k4, h/6);  // Pour pas avoir fuite de memoire
    giant_vec_t* vec_plus_1 = GVecAdd(*vec, *k1_plus_2k2_plus_2k3_plus_k4__times_h_6); 
    // Je pensait de faire ceci des variables globales que on peux free a la fin de notre programe pour aleger le programme mais alors le programme ne peux pas etre thread safe car si on le laisse comme ca ca peux etre mieux pour le futur
    // Une autre option est de faire une fonction qui fait la partie (k1+2k2+2k3+k4)h/6 et ca nous alegerait de 4 creation et de free TODO peux etre pour le futur
    FreeGVec(k1_times_h_2);
    FreeGVec(vec_plus_k1_times_h_2);
    FreeGVec(k2_times_h_2);
    FreeGVec(vec_plus_k2_times_h_2);
    FreeGVec(k3_times_h);
    FreeGVec(vec_plus_k3_times_h);
    FreeGVec(k1_plus_k4);
    FreeGVec(k2_plus_k3);
    FreeGVec(k2_plus_k3_times_2);
    FreeGVec(k1_plus_2k2_plus_2k3_plus_k4);
    FreeGVec(k1_plus_2k2_plus_2k3_plus_k4__times_h_6);
    return vec_plus_1;
}

void SolGiantVecFuncVal(giant_vec_func_vals_t* vf, functypeGIANTvec func_to_solve,
                     tab_t* vars_needed, giant_vec_t cond_init, functypeEDOsolverGIANTvec solver){
    if (vf->len < 2){
        printf("Function has less than 2 values allocated to it that doesnt work\n");
        exit(-1);
    }
    vf->gvecs[0] = InitGVec(cond_init.len);
    if(vf->gvecs[0]->len != cond_init.len){
        printf("There is a problem with the lengths\n");
    }
    GVecEq(vf->gvecs[0], cond_init);
    double h = vf->t[1] - vf->t[0]; // since our h doesnt change and we chose it acoringly to n_max in main this is faster to do than writing a h each time we call this function
    for(int i = 1; i < vf->len; i++){
        vf->gvecs[i] = solver(vf->gvecs[i-1], vars_needed, func_to_solve, h);
    }
}



