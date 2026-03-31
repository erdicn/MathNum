#ifndef MYFLOAT_H
#define MYFLOAT_H

// typedef float myfloat;
// typedef double myfloat;
#ifndef MYFLOAT_DEFINED
// typedef float myfloat;
typedef double myfloat;
#endif

#include <math.h>

#define MYGAMMA(x) ( \
    sizeof(myfloat) == sizeof(float)       ? tgammaf(x) : \
    sizeof(myfloat) == sizeof(double)      ? tgamma(x)  : \
                                             tgammal(x)   \
)

#define MYPOW(x, p) ( \
    sizeof(myfloat) == sizeof(float)       ? powf(x, p) : \
    sizeof(myfloat) == sizeof(double)      ? pow(x, p)  : \
                                             powl(x, p)   \
)

#define MYCBRT(x) ( \
    sizeof(myfloat) == sizeof(float)       ? cbrtf(x) : \
    sizeof(myfloat) == sizeof(double)      ? cbrt(x )  : \
                                             cbrtl(x)   \
)

#define MYABS(x) ( \
    sizeof(myfloat) == sizeof(float)       ? fabsf(x) : \
    sizeof(myfloat) == sizeof(double)      ? fabs (x)  : \
                                             fabsl(x)   \
)


#define MYEXP(x) ( \
    sizeof(myfloat) == sizeof(float)       ? expf(x) : \
    sizeof(myfloat) == sizeof(double)      ? exp(x)  : \
                                             expl(x)   \
)

#endif /* MYFLOAT_H */
