#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "../myfloat.h"
#include <stdbool.h>

myfloat pdfGAMMA(myfloat x, myfloat a, myfloat theta, 
                 myfloat* gamma_alpha_theta_a, bool recalculate);
myfloat pdfNormal(myfloat x, myfloat sig, myfloat mu);
#endif /* DISTRIBUTIONS_H */