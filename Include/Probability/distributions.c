#define __USE_MISC
#include "distributions.h"
#include <math.h> // need tgamma and only from C11 and up
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>


// GAMMA 
// https://en.wikipedia.org/wiki/Gamma_distribution
myfloat pdfGAMMA(myfloat x, myfloat a, myfloat theta, 
                             myfloat* gamma_alpha_theta_a, bool recalculate){
    assert(a>0);
    assert(gamma_alpha_theta_a);
    myfloat g_a_t_a = recalculate ? (1./(MYGAMMA(a) * MYPOW(theta, a))) : 
                                    *gamma_alpha_theta_a;   
    return g_a_t_a*MYPOW(x, a-1)*MYEXP(-x/theta);
}  

myfloat pdfNormal(myfloat x, myfloat sig, myfloat mu){
    return 1./(sqrt(2*M_PI)*sig) * MYEXP(-pow(x-mu, 2)/(2*sig*sig));
} // TODO optimised no sqrt only precomputed 

// myfloat cdfGAMMA(myfloat x, myfloat a, myfloat theta, 
//                              myfloat* gamma_alpha_theta_a, bool recalculate){
//     assert(a>0);
//     assert(gamma_alpha_theta_a);
//     myfloat g_a_t_a = recalculate ? (1./(MYGAMMA(a) * MYPOW(theta, a))) : 
//                                     *gamma_alpha_theta_a;   
//     return g_a_t_a*MYPOW(x, a-1)*MYEXP(-x/theta);
// }  