#include <stdio.h>
#include <stdlib.h>

#include "functions.h"
#include "vectors.h"




void FreeFuncVals(func_vals_t* func){
    free(func->x);
    free(func->f);
    free(func);
}

void FreeFuncVals2(func_vals2_t* func){
    free(func->x);
    free(func->y);
    for(int i = 0; i < func->leny; i++){
        free(func->f[i]);
    }
    free(func->f);
    free(func);
}

void FreeVecFuncVals(vec_func_vals_t* func){
    free(func->t);
    free(func->vec);
    free(func);
}


void FreeTab(tab_t* tab){
    free(tab->vals);
    free(tab);
}

void PrintTab(tab_t tab){
    for(int i = 0; i < tab.len; i++)
        printf("%lf ", tab.vals[i]);
    puts("");
}

tab_t* InitTab(int len_tab){
    tab_t* new_tab = malloc(sizeof(tab_t));
    new_tab->len = len_tab;
    new_tab->vals = malloc(sizeof(double)*new_tab->len);
    return new_tab;
}

func_vals_t* InitFuncVals(int len){
    func_vals_t* new_func = malloc(sizeof(func_vals_t));
    new_func->len=len;
    new_func->x = malloc(sizeof(double)*new_func->len);
    new_func->f = malloc(sizeof(double)*new_func->len);
    return new_func;
}

func_vals_t* InitFuncValsWithInterval(int nb_intervals, double min, double max){
    // the +1s are nescescary to include the max 
    func_vals_t* new_func = InitFuncVals(nb_intervals+1);
    double h = (max - min)/nb_intervals;
    for(int i = 0; i < nb_intervals+1; i++){
        new_func->x[i] = min + h*i;
    }
    return new_func;
}

vec_func_vals_t* InitVecFuncVals(int len){
    vec_func_vals_t* new_func = malloc(sizeof(vec_func_vals_t));
    new_func->len = len;
    new_func->t   = malloc(sizeof(double)*new_func->len);
    new_func->vec = malloc(sizeof(vec2D_t)*new_func->len);
    return new_func;
}

vec_func_vals_t* InitVecFuncValsWithInterval(int nb_intervals, double min, double max){
    // the +1s are nescescary to include the max 
    vec_func_vals_t* new_func = InitVecFuncVals(nb_intervals+1);
    double h = (max - min)/nb_intervals;
    for(int i = 0; i < nb_intervals+1; i++){
        new_func->t[i] = min + h*i;
    }
    return new_func;
}

func_vals2_t* InitFuncVals2(int lenx, int leny){
    func_vals2_t* new_func = malloc(sizeof(func_vals_t));
    new_func->lenx=lenx;
    new_func->leny=leny;
    new_func->x = malloc(sizeof(double)*new_func->lenx);
    new_func->y = malloc(sizeof(double)*new_func->leny);
    new_func->f = malloc(sizeof(double*)*new_func->leny);
    for(int i = 0; i < new_func->leny; i++){
        new_func->f[i] = malloc(sizeof(double)*new_func->lenx);
    }
    return new_func;
}

tab_t* ReadTab(int len_tab, char* nom_fichier){
    tab_t* new_tab = InitTab(len_tab);

    FILE* fichier = fopen(nom_fichier, "r");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier);
        return NULL;
    }
    // maintenant que notre fichier est ouvert on peux le lire
    for(int i = 0; i < len_tab; i++){
        fscanf(fichier, "%lf", new_tab->vals + i);
    }
    fclose(fichier);
    return new_tab;
}

func_vals_t* ReadFuncAndValsXandF_InRows(char* nom_fichier){
    FILE* fichier = fopen(nom_fichier, "r");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier);
        return NULL;
    }

    func_vals_t* func = malloc(sizeof(func_vals_t));
    fscanf(fichier,"%d", &(func->len));
    func->x = malloc(sizeof(double)*func->len);
    func->f = malloc(sizeof(double)*func->len);
    for (int i = 0; i < func->len; i++){
        fscanf(fichier, "%lf ", func->x + i);
    }
    for (int i = 0; i < func->len; i++){
        fscanf(fichier, "%lf ", func->f + i);
    }
    
    fclose(fichier);
    return func;
}

func_vals_t* ReadFuncAndValsXandF_InColumns(char* nom_fichier){
    FILE* fichier = fopen(nom_fichier, "r");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier);
        return NULL;
    }

    func_vals_t* func = malloc(sizeof(func_vals_t));
    fscanf(fichier,"%d", &(func->len));
    func->x = malloc(sizeof(double)*func->len);
    func->f = malloc(sizeof(double)*func->len);
    for (int i = 0; i < func->len; i++){
        fscanf(fichier, "%lf ", func->x + i);
        fscanf(fichier, "%lf ", func->f + i);
    }
    
    fclose(fichier);
    return func;
}

void PrintFunc(func_vals_t func){
    for (int i = 0; i < func.len; i++){
        printf("%.16lf %.16lf\n", func.x[i], func.f[i]);
    }
}

void PrintVecFuncVals(vec_func_vals_t func){
    for (int i = 0; i < func.len; i++){
        // printf("%.16lf %.16lf %.16lf\n",func.t[i], func.x[i], func.y[i]);
        printf("%.16lf %.16lf %.16lf\n", func.t[i], func.vec[i].x, func.vec[i].y);
    }
}


void EcritureFunc_Vals(func_vals_t func_to_print, char* nom_fichier_de_print){
    FILE* fichier = fopen(nom_fichier_de_print, "w");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier_de_print);
    }
    for(int i = 0; i < func_to_print.len; i++){
        fprintf(fichier, "%.16lf %.16lf\n", func_to_print.x[i], func_to_print.f[i]);
    }
    fclose(fichier);
}


void EcritureVecFuncVals(vec_func_vals_t func_to_print, char* nom_fichier_de_print){
    FILE* fichier = fopen(nom_fichier_de_print, "w");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier_de_print);
    }
    for(int i = 0; i < func_to_print.len; i++){
        fprintf(fichier, "%.16lf %.16lf %.16lf\n", func_to_print.t[i], func_to_print.vec[i].x, func_to_print.vec[i].y);
    }
    fclose(fichier);
}

func_vals_t* GvecFuncTuFuncVals(int gvec_index, giant_vec_func_vals_t gvec_func_to_equalise){
    func_vals_t* new_func_vals = InitFuncVals(gvec_func_to_equalise.len);
    for(int i = 0; i < new_func_vals->len; i++){
        new_func_vals->x[i] = gvec_func_to_equalise.t[i];
        new_func_vals->f[i] = gvec_func_to_equalise.gvecs[i]->comp[gvec_index];
    }
    return new_func_vals;
}