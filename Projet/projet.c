#include "../Include/all.h"
#include <string.h>

double Calcul_X(double t){
    return exp(-(3./8.)*t) * (-2.*               cos((sqrt(311.) / 8.) * t) 
                              -(6./sqrt(311.)) * sin((sqrt(311.) / 8.) * t)) + 2.;
}

double Calcul_V(double t){
    return (3. / 4.) * exp(-(3. / 8.) * t) * cos((sqrt(311.) / 8.) * t) 
            + (sqrt(311.) / 4) * exp(-(3. / 8.) * t) * sin((sqrt(311.) / 8.) * t) 
            + (9. / (4. * sqrt(311.))) * sin((sqrt(311.) / 8.) * t) * exp(-(3. / 8.) * t) 
            - (3. / 4.) * exp(-(3. / 8.) * t) * cos((sqrt(311.) / 8.) * t);
}

// Ressort masse simple
double SolAnalythiqueSimpleRessort(double t, tab_t* var_tab){
    double eps = 1e-13;
    double m = var_tab->vals[0], k = var_tab->vals[1], c = var_tab->vals[2], g_s = var_tab->vals[3];
    double dis = (c*c - 4*k*m)/(m*m);
    if (dis < 0 - eps){
        double K1 = -g_s/k;
        double K2 = (-c*g_s)/(m*k*sqrt(-dis));
        return exp(-c/(2.*m) * t) * (K1*cos(sqrt(-dis)/2. * t) + K2*sin(sqrt(-dis)/2. * t)) + g_s/k;
    } else if (dis > 0 + eps){
        double C1 = (-g_s/k - (c*g_s)/(m*k*sqrt(dis))) / 2.;
        double C2 = -g_s/k - C1;
        return exp(-c/(2.*m) * t) * (C1*exp(sqrt(dis)/2. * t) + C2*exp(-sqrt(dis)/2. * t)) + g_s/k;
    } else {
        double D1 = -g_s/k;
        double D2 = -(c*g_s)/(2.*m*k);
        return (D1 + D2*t) * exp(-c/(2*m) * t) + g_s/k;
    }
}

double f1_tab(vec2D_t U, tab_t* var_tab){
    return U.y;
}

double f2_tab(vec2D_t U, tab_t* var_tab){
    if(var_tab->len != 4){
        printf("Error the tab is not made for this function\n");
        exit(-3);
    }
    double m = var_tab->vals[0], k = var_tab->vals[1], c = var_tab->vals[2], g_s = var_tab->vals[3];
    return -(k/m)*U.x - (c/m)*U.y + (g_s/m);
}

vec2D_t F_tab(vec2D_t U, tab_t* var_tab){
    vec2D_t new_vec = {.x = f1_tab(U, var_tab), .y = f2_tab(U, var_tab)};
    return new_vec;
}

giant_vec_t* gVecRessort1Func(giant_vec_t* U, tab_t* var_tab){
    giant_vec_t* new_vec = InitGVec(2);
    vec2D_t nU = {.x = U->comp[0], .y = U->comp[1]};
    new_vec->comp[0] = f1_tab(nU, var_tab);
    new_vec->comp[1] = f2_tab(nU, var_tab);
    return new_vec;
}

// TAIPEI
// giant vec is going to be x1 x2 v1 v2
//////////////////////////////////////////////////////////////////////////////////
double taipei0(giant_vec_t* vec,tab_t* vars){
    return vec->comp[2];
}

double taipei1(giant_vec_t* vec,tab_t* vars){
    return vec->comp[3];
}
// vars are m1, k1, m2, k2, c
double taipei2(giant_vec_t* vec,tab_t* vars){
    if(vars->len != 5){
        printf("Vars are not at the good length %d\n", vars->len);
        exit(-99);
    }
    if (vec->len != 4){
        printf("Vec is not at the good length %d\n", vec->len);
        exit(-66);
    }
    double x1 = vec->comp[0], x2 = vec->comp[1], v1 = vec->comp[2], v2 = vec->comp[3];
    double m1 = vars->vals[0], k1 = vars->vals[1], m2 = vars->vals[2], k2 = vars->vals[3], c = vars->vals[4];
    return -((k1+k2)*x1)/m1 + (k2*x2)/m1;
}

double taipei3(giant_vec_t* vec,tab_t* vars){
    if(vars->len != 5){
        printf("Vars are not at the good length %d\n", vars->len);
        exit(-99);
    }
    if (vec->len != 4){
        printf("Vec is not at the good length %d\n", vec->len);
        exit(-66);
    }
    double x1 = vec->comp[0], x2 = vec->comp[1], v1 = vec->comp[2], v2 = vec->comp[3];
    double m1 = vars->vals[0], k1 = vars->vals[1], m2 = vars->vals[2], k2 = vars->vals[3], c = vars->vals[4];
    return (k2*x1)/m2 - (k2*x2)/m2 + c*v1/m2 - c*v2/m2;
}

giant_vec_t* TAIPEI(giant_vec_t* vec, tab_t* vars){
    giant_vec_t* new_vec = InitGVec(vec->len);
    new_vec->comp[0] = taipei0(vec, vars);
    new_vec->comp[1] = taipei1(vec, vars);
    new_vec->comp[2] = taipei2(vec, vars);
    new_vec->comp[3] = taipei3(vec, vars);
    return new_vec;
}



// TODO cette fonction peux etre mieux otimiser mais comme on va pas lutiliser avec des str long loptimisation est negligable par rapport au reste du programme
char* MakeFilename(char* name1, char* name2){
    char* new_str = malloc(sizeof(char)*(strlen(name1)+strlen(name2)+1+4)); // +1 for \0 et +4 pour .dat
    int i1 = 0, i2 = 0;
    while(name1[i1]){ new_str[i1]    = name1[i1]; i1++; }
    while(name2[i2]){ new_str[i1+i2] = name2[i2]; i2++; }
    char* dat = ".dat";
    for(int i = 0; i < 5; i++){ new_str[i1+i2+i] = dat[i]; } // jusq a 5 pour le \0
    return new_str;
}

functypeEDOsolverBIGvec   vec_solvers[]   = {&VecEulersMethod,   &VecRK2EulersMethod,  &VecRK2HeunMethod,   &VecRK4Method  };
functypeEDOsolverGIANTvec giant_solvers[] = {&GiantEulersMethod, &GiantRK2EulerMethod, &GiantRK2HeunMethod, &GiantRK4Method};
char  s_int[20]; // 20 car ax int est 2147483647 et unpeut de marge 
char* write_file = NULL;
char* simple_filenames[4] = {"Donnees/Simple/euler", "Donnees/Simple/RK2euler", "Donnees/Simple/RK2heun", "Donnees/Simple/RK4"};

void SolveSimple();
void SolveTaipei();

int main(){

    SolveSimple();
    // SolveTaipei();
    // printf("Integral %.16lf\n", CompositeSimposon(0, gU->gvecs[gU->len - 1]->comp[0],(*GvecFuncTuFuncVals(0, *gU))).val);

    return 0;
}

// #define GIANT_VEC_SOLVE

void SolveSimple(){
    // Variables du debut
    //       kg    N/m     kg/s      N           s  
    double m = 2., k = 10., c = 1.5, g_s = 20., t_max = 10.;
    tab_t* simple_vars = InitTab(4);
    simple_vars->vals[0] = m; 
    simple_vars->vals[1] = k; 
    simple_vars->vals[2] = c;
    simple_vars->vals[3] = g_s;

    int nb_interval_ana = 1000;
    func_vals_t* sol_ana = InitFuncValsWithInterval(nb_interval_ana, 0, t_max);
    for(int i = 0; i < sol_ana->len; i++){
        sol_ana->f[i] = SolAnalythiqueSimpleRessort(sol_ana->x[i], simple_vars);
    }
    EcritureFunc_Vals(*sol_ana, "sol_ana.dat");
    FreeFuncVals(sol_ana);

    int n_max = 50; 
    vec2D_t vec_simple_cond_init = {.x = 0, .y = 0};

    vec_func_vals_t* U = InitVecFuncValsWithInterval(n_max, 0, t_max); // double dt = t_max/n_max;
    SolVecFuncVal(U, &F_tab, simple_vars, vec_simple_cond_init, &VecRK4Method);
    EcritureVecFuncVals(*U, "vec_RK4.dat");
    FreeVecFuncVals(U);

    giant_vec_t* gvec_cond_init = InitGVec(2);
    gvec_cond_init->comp[0] = vec_simple_cond_init.x;
    gvec_cond_init->comp[1] = vec_simple_cond_init.y;

    int   simple_n_maxs[7]    = {100, 500, 1000, 5000, 10000, 50000, 100000};
    // for(int i = 0; i < 7; i++) simple_n_maxs[i] *= 10;
    giant_vec_func_vals_t* gU;

    FILE* error_file = fopen("errors_simple.dat", "w");
    if(error_file == NULL){
        printf("Could not open error file\n");
        exit(-99);
    }
#ifndef GIANT_VEC_SOLVE
    for(int i = 0; i < 7; i++){
        for(int j = 0; j < 4; j++){
            // gU = InitGvecFuncValsInterval(simple_n_maxs[i], 0, t_max);
            // SolGiantVecFuncVal(gU, &gVecRessort1Func, simple_vars, *gvec_cond_init, giant_solvers[j]);
            U = InitVecFuncValsWithInterval(simple_n_maxs[i], 0, t_max);
            SolVecFuncVal(U, &F_tab, simple_vars, vec_simple_cond_init, vec_solvers[j]);
            sprintf(s_int,"_%d", simple_n_maxs[i]);
            write_file = MakeFilename(simple_filenames[j], s_int);
            EcritureVecFuncVals(*U, write_file);
            // EcritureGiantVecFuncVals(*gU, write_file);
            // printf("%.16lf %d\n", fabs(gU->gvecs[gU->len-1]->comp[0] - SolAnalythiqueSimpleRessort(gU->t[gU->len-1], simple_vars)), simple_n_maxs[i]);
            fprintf(error_file, "%d %d %.16lf %.16lf\n", j, simple_n_maxs[i],
                                                  fabs(U->vec[U->len-1].x - SolAnalythiqueSimpleRessort(U->t[U->len-1], simple_vars)),
                                                  // fabs(gU->gvecs[gU->len-1]->comp[0] - SolAnalythiqueSimpleRessort(gU->t[gU->len-1], simple_vars)),
                                                  t_max/simple_n_maxs[i]);
            free(write_file);
            // FreeGvecFuncVals(gU);
            FreeVecFuncVals(U);
        }
    }
#endif
#ifdef GIANT_VEC_SOLVE
    for(int i = 0; i < 7; i++){
        for(int j = 0; j < 4; j++){
            gU = InitGvecFuncValsInterval(simple_n_maxs[i], 0, t_max);
            SolGiantVecFuncVal(gU, &gVecRessort1Func, simple_vars, *gvec_cond_init, giant_solvers[j]);
            // U = InitVecFuncValsWithInterval(simple_n_maxs[i], 0, t_max);
            // SolVecFuncVal(U, &F_tab, simple_vars, vec_simple_cond_init, vec_solvers[j]);
            sprintf(s_int,"_%d", simple_n_maxs[i]);
            write_file = MakeFilename(simple_filenames[j], s_int);
            // EcritureVecFuncVals(*U, write_file);
            EcritureGiantVecFuncVals(*gU, write_file);
            // printf("%.16lf %d\n", fabs(gU->gvecs[gU->len-1]->comp[0] - SolAnalythiqueSimpleRessort(gU->t[gU->len-1], simple_vars)), simple_n_maxs[i]);
            fprintf(error_file, "%d %d %.16lf %.16lf\n", j, simple_n_maxs[i],
                                                  fabs(U->vec[U->len-1].x - SolAnalythiqueSimpleRessort(U->t[U->len-1], simple_vars)),
                                                  // fabs(gU->gvecs[gU->len-1]->comp[0] - SolAnalythiqueSimpleRessort(gU->t[gU->len-1], simple_vars)),
                                                  t_max/simple_n_maxs[i]);
            free(write_file);
            FreeGvecFuncVals(gU);
            // FreeVecFuncVals(U);
        }
    }
#endif
    fclose(error_file);
}

void SolveTaipei(){
        // TAIPEI
    double m1 = 264*pow(10, 6), k1 = 225*pow(10, 6), 
           m2 = 660*pow(10, 3), k2 = 510*pow(10, 3),
           taipei_c  =  52*pow(10, 3);
    tab_t* tapei_vars = InitTab(5);
    tapei_vars->vals[0] = m1;
    tapei_vars->vals[1] = k1;
    tapei_vars->vals[2] = m2;
    tapei_vars->vals[3] = k2;
    tapei_vars->vals[4] = taipei_c;

    double x0 = 3, x10 = x0, x20 = 0, v10 = 0, v20 = 0;
    giant_vec_t* taipei_cond_init = InitGVec(4);
    taipei_cond_init->comp[0] = x10;
    taipei_cond_init->comp[1] = x20;
    taipei_cond_init->comp[2] = v10;
    taipei_cond_init->comp[3] = v20;
    
    double taipei_t_max = 240;
    int taipei_n_max = 10000;
    char*  taipei_filenames[4] = {"Donnees/Taipei/euler", "Donnees/Taipei/RK2euler", "Donnees/Taipei/RK2heun", "Donnees/Taipei/RK4"};

    // pour x0 = 3
    giant_vec_func_vals_t* G; 
    taipei_n_max = 2000;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 4; j++){
            G = InitGvecFuncValsInterval(taipei_n_max, 0, taipei_t_max);
            SolGiantVecFuncVal(G, &TAIPEI, tapei_vars, *taipei_cond_init, giant_solvers[j]);
            sprintf(s_int, "x3_%d", taipei_n_max);
            write_file = MakeFilename(taipei_filenames[j], s_int);
            EcritureGiantVecFuncVals(*G, write_file);
            FreeGvecFuncVals(G);  
            free(write_file);
        }
        taipei_n_max *= 10;
    }

    // pour x0 = 0.25
    taipei_n_max = 2000;
    taipei_cond_init->comp[0] = x10;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 4; j++){
            G = InitGvecFuncValsInterval(taipei_n_max, 0, taipei_t_max);
            SolGiantVecFuncVal(G, &TAIPEI, tapei_vars, *taipei_cond_init, giant_solvers[j]);
            sprintf(s_int, "x025_%d", taipei_n_max);
            write_file = MakeFilename(taipei_filenames[j], s_int);
            EcritureGiantVecFuncVals(*G, write_file);
            FreeGvecFuncVals(G);  
            free(write_file);
        }
        taipei_n_max *= 10;
    }
}