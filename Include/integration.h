#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "functions.h"

val_err CompositeTrapezes(double lower_bound, double upper_bound, func_vals_t func);
val_err CompositeSimposon(double lower_bound, double upper_bound, func_vals_t func);

#endif // INTEGRATION_H