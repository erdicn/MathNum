#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h> 
#include <math.h>
typedef struct Tableau{
    int len;
    double* vals;
} tab_t;

typedef struct Func_Vals{
    int len;
    double* x;
    double* f;
}func_vals_t;

void FreeFuncVals(func_vals_t* func){
    free(func->x);
    free(func->f);
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
        printf("%lf %lf\n", func.x[i], func.f[i]);
    }
}

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

void EcritureFunc_Vals(func_vals_t func_to_print, char* nom_fichier_de_print){
    FILE* fichier = fopen(nom_fichier_de_print, "w");
    if (fichier == NULL){
        printf("Unable to open %s file\n", nom_fichier_de_print);
    }
    for(int i = 0; i < func_to_print.len; i++){
        fprintf(fichier, "%lf %lf\n", func_to_print.x[i], func_to_print.f[i]);
    }
}

typedef struct Value_Error{
    double val;
    double err;
}val_err;

#define EPS 3e-12
// precondition les x de func sont ecarte equalement x1-x0 = x2-x1
val_err CompositeTrapezes(double lower_bound, double upper_bound, func_vals_t func){
    if (func.len < 2){
        printf("La fonction rentre a moins que 2 points besoin au moins de 2 points\n");
    }
    double h = func.x[1] - func.x[0];
    int counter = 0;
    double first = 0, last = 0, tot = 0;
    for(int i = 0; i < func.len; i ++){
        if (func.x[i] >= lower_bound - EPS && func.x[i] <= upper_bound + EPS){
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


#define REAL_VAL 448.0/9.0

int main(){
    char* nom_fichier_donnees_init = "donnees.dat";
    func_vals_t* func_TP1 = ReadFuncAndValsXandF_InColumns(nom_fichier_donnees_init);//ReadFuncAndValsXandF_InRows(nom_fichier_donnees_init);
    PrintFunc(*func_TP1);
    printf("%lf\n", Calcul_Pn(1, *func_TP1));
    
    FILE* fichier_erreurs_trap = fopen("trap_err.dat", "w");
    func_vals_t* points_calcule = Calcul_P_PointParPointAvecM(*func_TP1, 10);
    EcritureFunc_Vals(*points_calcule, "lagrange.dat");
    val_err result;

    for(int m = 10; m <= 1e7; m = m*10){
        points_calcule = Calcul_P_PointParPointAvecM(*func_TP1, m);
        result =  CompositeTrapezes(-2, 6, *points_calcule);
        printf("comp_trap = %.16lf +- %.16lf %.16lf %d\n",result.val, result.err, fabs(result.val-(REAL_VAL)), m);
        fprintf(fichier_erreurs_trap, "%.16lf %.16lf\n", points_calcule->x[1]-points_calcule->x[0], fabs(result.val-(REAL_VAL)));//, result.val, result.err, m);%.16lf %.16lf %d
        FreeFuncVals(points_calcule);
    }
    fclose(fichier_erreurs_trap);
    FILE* fichier_erreurs_simp = fopen("simp_err.dat", "w");

    for(int m = 10; m <= 1e7; m*=10){
        points_calcule = Calcul_P_PointParPointAvecM(*func_TP1, m);
        result =  CompositeSimposon(-2, 18, *points_calcule);
        printf("comp_simp = %.16lf +- %.16lf %.16lf %d\n", result.val, result.err, fabs(result.val-(REAL_VAL)), m);
        fprintf(fichier_erreurs_simp, "%.16lf %.16lf\n", points_calcule->x[1]-points_calcule->x[0],  fabs(result.val-(REAL_VAL)));//, result.val, result.err, m);%.16lf %.16lf %d
        FreeFuncVals(points_calcule);
    }
    fclose(fichier_erreurs_simp);

    for(double h= 1; h > 1e-7; h = h/10.){
        points_calcule = Calcul_P_PointParPointAvecH(*func_TP1, h);
        result =  CompositeTrapezes(-2, 6, *points_calcule);
        printf("comp_trap = %.16lf +- %.16lf %.16lf %.16lf\n",result.val, result.err, fabs(result.val-(REAL_VAL)), h);
        //fprintf(fichier_erreurs_trap, "%.16lf %.16lf\n", points_calcule->x[1]-points_calcule->x[0], fabs(result.val-(REAL_VAL)));//, result.val, result.err, m);%.16lf %.16lf %d
        FreeFuncVals(points_calcule);
    }

    FreeFuncVals(func_TP1);
    printf("vraie val = %.16lf\n", REAL_VAL);


    return 0;
}