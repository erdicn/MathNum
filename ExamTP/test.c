#include <stdio.h>
#include "../Include/edo_solver.h"


double FuncOldExam(double x, double y){
    return x+y;
}


int main(){
    printf("Result: %lf %lf\n", RK4Method(0, 0, &FuncOldExam, 1./2.), 41./(32.*12.));
    return 0;
}
