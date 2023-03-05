#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../Include/functions.h"
#include "../Include/integration.h"
#include "../Include/interpolation.h"


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