#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#ifdef _OPENMP
#include <omp.h> // even thought we live in a 1D world we should enjoy using the whole computer
#include "OpenMP/test_openmp.h"
#endif

// #define USE_CUDA 1

#ifdef USE_CUDA
#define MYFLOAT_DEFINED 
typedef float myfloat;

#define HOST_DEVICE __host__ __device__
#include <cuda_runtime.h>
#else 
#define HOST_DEVICE 
#endif

#include "matrix.h" // my sparse matrix and vector header defines myfloat as double
#include "file_manips.h" // my sparse matrix and vector header defines myfloat as double


#define MAX(x, y) x > y ? x : y
#define MIN(x, y) x < y ? x : y

#define FOR_EACH_VEC(vec, i) for(int i = 0; i < vec->len; i++) vec->vals[i]
#define FOR_VEC_LEN(vec, i) for(int i = 0; i < vec->len; i++)
#define SWAP(x, y, T) do {T tmp = x; x = y; y = tmp;}while(0)

static inline void errorMessage(char* error, int line, const char* file){
    fprintf(stderr, "Error (%s) found in line %d in file %s\n", error, line, file);
    abort();
    return;
}

#define ERROR_MESSAGE(string_err_message)  errorMessage(string_err_message, __LINE__, __FILE__)
#define MY_INT_ASSERT(x, y)  if(x!=y) ERROR_MESSAGE("assert failed")
#define CHECK_ERROR(nb)   if(nb != 0) ERROR_MESSAGE("check_error")

#define MAX_FILENAME_LEN 128
#define G 9.80665
#define GAMMA 1

#ifdef EPS
#   undef EPS
#endif

#define EPS 1e-30
#define CFL 0.1
#define REFLECTIVE 1

fvec_t* getFaceCenteredPos(fvec_t* xp, fvec_t* faces, myfloat xmin, myfloat xmax){
    if(faces == NULL) faces = allocateFVec(faces, xp->len+1);

    faces->vals[0] = xmin;
    faces->vals[faces->len-1] = xmax;
    for(int i = 1; i < faces->len-1; i++){
        faces->vals[i] = 0.5*(xp->vals[i]-xp->vals[i-1]) + xp->vals[i-1];
    }
    return faces;
}

myfloat maxAbsEigenvalue(fvec_t* h, fvec_t* Q){
    assert(h->len == Q->len);
    myfloat max_val = -INFINITY;
    
    #if GAMMA!=1
        myfloat g2 = GAMMA*GAMMA;
    #endif
    #pragma omp parallel for
    for(int i = 0; i < h->len; i++){
        myfloat u = h->vals[i] > EPS ? Q->vals[i] / h->vals[i] : 0;
        #if GAMMA==1
            myfloat sq = sqrt(G*h->vals[i]);
            myfloat e1 = fabs(u + sq);
            myfloat e2 = fabs(u - sq);
        #else
            myfloat sq = sqrt((g2-GAMMA)*u*u + G*h->vals[i]);
            myfloat e1 = fabs(GAMMA*u + sq);
            myfloat e2 = fabs(GAMMA*u - sq);
        #endif
        myfloat maxe12 = MAX(e1, e2);
        max_val = MAX(max_val, maxe12);
    }
    return max_val;
}

// Foward first order derivative order 2 taken from https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
static inline myfloat forwardD1O2(myfloat dx, myfloat fx, myfloat fx_p1, myfloat fx_p2){
    return (-3*fx + 4*fx_p1 - fx_p2)/(2*dx);
}

static inline myfloat forwardD1O1(myfloat dx, myfloat fx, myfloat fx_p1){
    return (fx_p1 - fx)/dx;
}

static inline myfloat backwardD1O1(myfloat dx, myfloat fx, myfloat fx_m1){
    return (fx - fx_m1)/dx;
}

// Backwards first order derivative order 2 taken from https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
static inline myfloat backwardD1O2(myfloat dx, myfloat fx, myfloat fx_m1, myfloat fx_m2){
    // return (3*fx -4*fx_m1 + fx_m2)/(2*dx);
    return -forwardD1O2(dx, fx, fx_m1, fx_m2);
}

static inline myfloat funcKi(myfloat q, myfloat h){
    return h > 0 ? q*q/h + G*h*h*0.5 : 0;
}

typedef struct TimeSaver{
    uint64_t current_idx;
    uint64_t len;
    myfloat* to_save;
} TimeSaver_t;

int saveTimesteps(TimeSaver_t* times, myfloat current_time, fvec_t* vec_to_save, myfloat dt, char* filename_prefix){
    if(times->current_idx >= times->len){
        return -1;
    }
    myfloat current_time_to_save = times->to_save[times->current_idx];
    if(current_time_to_save-dt/2<= current_time &&  current_time <= current_time_to_save+dt/2 ){
        char full_filename[MAX_FILENAME_LEN+64];// the +64 is for time storing
        sprintf(full_filename, "%s%lf.dat", filename_prefix, current_time);
        saveVecToFile(vec_to_save, full_filename);
        times->current_idx++;
        return 1;
    }
    return 0;
}

void naiveUpwinding(fvec_t* h, fvec_t* q, myfloat dt, myfloat end_t, myfloat t, myfloat DX, TimeSaver_t* timesteps_to_save){
    fvec_t *h_tmp  = newCopyFVec(NULL, h);
    fvec_t *q_star = newCopyFVec(NULL, q);
    myfloat dq_dx, dki_dx, ki;
    int i;
    while(t <= end_t){
        // printf("t=%lf dt=%lf  DX=%lf end_t=%lf| ", t, dt, DX, end_t); printFVec(h);
        // calculating hn+1
        i = 1;
        dq_dx = q->vals[i] > 0 ? backwardD1O1(DX, q->vals[i], q->vals[i-1]) :
                                  forwardD1O2(DX, q->vals[i], q->vals[i+1], q->vals[i+2]) ;
        h_tmp->vals[i] = -dq_dx*dt + h->vals[i];
        h_tmp->vals[0] = h_tmp->vals[i];
        for(i = 2; i < h->len-2; i++) {
            // Upwinding
            dq_dx = q->vals[i] > 0 ? backwardD1O2(DX, q->vals[i], q->vals[i-1], q->vals[i-2]) :
                                      forwardD1O2(DX, q->vals[i], q->vals[i+1], q->vals[i+2]) ;
            h_tmp->vals[i] = -dq_dx*dt + h->vals[i];
        }
        dq_dx = q->vals[i] > 0 ? backwardD1O2(DX, q->vals[i], q->vals[i-1], q->vals[i-2]) :
                                  forwardD1O1(DX, q->vals[i], q->vals[i+1]) ;
        h_tmp->vals[i] = -dq_dx*dt + h->vals[i];
        h_tmp->vals[i+1] = h_tmp->vals[i];

        // Q*
        i = 1;
        ki = funcKi(q->vals[i], h_tmp->vals[i]);
        dki_dx = ki > 0 ? backwardD1O1(DX, funcKi( q->vals[i], h_tmp->vals[i]),  funcKi(q->vals[i-1], h_tmp->vals[i-1])) :
                           forwardD1O2(DX, funcKi( q->vals[i], h_tmp->vals[i]),  funcKi(q->vals[i+1], h_tmp->vals[i+1]),  funcKi(q->vals[i+2], h_tmp->vals[i+2])) ;
        q_star->vals[i] = -dki_dx*dt + q->vals[i];
        q_star->vals[0] = q_star->vals[i];
        for(i = 2; i < h->len-2; i++){
            // Upwinding
            ki = funcKi(q->vals[i], h_tmp->vals[i]);
            dki_dx = ki > 0 ? backwardD1O2(DX, funcKi( q->vals[i], h_tmp->vals[i]),  funcKi(q->vals[i-1], h_tmp->vals[i-1]),  funcKi(q->vals[i-2], h_tmp->vals[i-2])) :
                              forwardD1O2 (DX, funcKi( q->vals[i], h_tmp->vals[i]),  funcKi(q->vals[i+1], h_tmp->vals[i+1]),  funcKi(q->vals[i+2], h_tmp->vals[i+2])) ;
            q_star->vals[i] = -dki_dx*dt + q->vals[i];
        }
        ki = funcKi(q->vals[i], h_tmp->vals[i]);
        dki_dx = ki > 0 ? backwardD1O2(DX, funcKi( q->vals[i], h_tmp->vals[i]),  funcKi(q->vals[i-1], h_tmp->vals[i-1]), funcKi(q->vals[i-2], h_tmp->vals[i-2])) :
                           forwardD1O1(DX, funcKi( q->vals[i], h_tmp->vals[i]),  funcKi(q->vals[i+1], h_tmp->vals[i+1])) ;
        q_star->vals[i] = -dki_dx*dt + q->vals[i];
        q_star->vals[i+1] = q_star->vals[i];

        // Normal Q with the friction/slope or not
        for(i = 0; i < q->len; i++){
            q->vals[i] = q_star->vals[i];
        }

        char filename[64];
        sprintf(filename, "Out/");
        SWAP(h, h_tmp, fvec_t*);
        saveTimesteps(timesteps_to_save, t, h, dt, filename);
        // printf("%lf | \n", t);
        // printFVec(h);

        t += dt;
        // dt = CFL * DX / maxAbsEigenvalue(h, q);
        // printf("%lf\n", maxAbsEigenvalue(h, q));
    }
}

static inline myfloat computeWaveSpeed(myfloat u, myfloat h) {
    return fabs(u) + sqrt(G * h);
}

static inline myfloat rusanovMassFlux(myfloat ug, myfloat ud, myfloat hg, myfloat hd, myfloat c_max) {
    return 0.5 * (hg*ug + hd*ud) - 0.5 * c_max * (hd - hg);
}

static inline myfloat rusanovMomentumFlux(myfloat ug, myfloat ud, myfloat hg, myfloat hd, myfloat c_max) {
    myfloat mom_g = ug*ug*hg + 0.5*G*hg*hg;
    myfloat mom_d = ud*ud*hd + 0.5*G*hd*hd;
    return 0.5 * (mom_g + mom_d) - 0.5 * c_max * (hd*ud - hg*ug);
}

void rusanovMethod(fvec_t* h, fvec_t* q, fvec_t* Z, myfloat dt, myfloat end_t, myfloat t, myfloat DX, TimeSaver_t* timesteps_to_save){
    fvec_t *h_tmp  = newCopyFVec(NULL, h);
    fvec_t *q_star = newCopyFVec(NULL, q);
    myfloat dq_dx, dki_dx, ki;
    int i;
    
    fvec_t *f_mass = allocateFVec(NULL, q->len + 1);
    fvec_t *f_mom  = allocateFVec(NULL, q->len + 1);
    while(t <= end_t) {
        
        #pragma omp parallel for
        for(i = 1; i < h->len; i++) {
            myfloat hg = h->vals[i-1];
            myfloat hd = h->vals[i];
            
            myfloat ug = hg > EPS ? q->vals[i-1] / hg : 0.0;
            myfloat ud = hd > EPS ? q->vals[i] / hd : 0.0;

            myfloat cg = computeWaveSpeed(ug, hg);
            myfloat cd = computeWaveSpeed(ud, hd);
            myfloat c_max = MAX(cg, cd); 

            f_mass->vals[i] = rusanovMassFlux(ug, ud, hg, hd, c_max);
            f_mom->vals[i]  = rusanovMomentumFlux(ug, ud, hg, hd, c_max);
        }

        f_mass->vals[0] = f_mass->vals[1];
        f_mom->vals[0]  = f_mom->vals[1];
        f_mass->vals[h->len] = f_mass->vals[h->len-1];
        f_mom->vals[h->len]  = f_mom->vals[h->len-1];

        #pragma omp parallel for
        for(i = 0; i < h->len; i++) {
            h_tmp->vals[i] = h->vals[i] - (dt / DX) * (f_mass->vals[i+1] - f_mass->vals[i]);

            if (h_tmp->vals[i] > EPS) {
                myfloat dZ_dx = (Z->vals[i+1] - Z->vals[i-1]) / (2.0 * DX);
                myfloat source_term = -G * h->vals[i] * dZ_dx;
                q_star->vals[i] = q->vals[i] - (dt / DX) * (f_mom->vals[i+1] - f_mom->vals[i]) + (dt*source_term);
            } else {
                h_tmp->vals[i] = 0.0; 
                q_star->vals[i] = 0.0;
            }
        }

        h_tmp->vals[0] = h_tmp->vals[1];
        q_star->vals[0] = q_star->vals[1];
        h_tmp->vals[h->len-1] = h_tmp->vals[h->len-2];
        q_star->vals[h->len-1] = q_star->vals[h->len-2];

// Tried but not working 
// #       ifdef REFLECTIVE
//             h_tmp->vals[0] = h_tmp->vals[1];
//             q_star->vals[0] = -q_star->vals[1]; 
            
//             h_tmp->vals[h->len-1] = h_tmp->vals[h->len-2];
//             q_star->vals[h->len-1] = -q_star->vals[h->len-2];
// #       endif

        #pragma omp parallel for
        for(i = 0; i < q->len; i++) {
            q->vals[i] = q_star->vals[i];
        }

        char filename[64];
        sprintf(filename, "Out/");
        SWAP(h, h_tmp, fvec_t*);
        saveTimesteps(timesteps_to_save, t, h, dt, filename);

        t += dt;
        
        // (Optional: You can uncomment your dynamic dt calculation here)
        dt = CFL * DX / maxAbsEigenvalue(h, q);
    }
}

// I wax thinking of using cuda with an explicit scheme but since dx -> 0 for stability dt -> 0 too so abandonned that idea 
HOST_DEVICE static inline myfloat d2u_dx2(myfloat uj_m1, myfloat uj, myfloat uj_p1,
                              myfloat xj_m1, myfloat xj, myfloat xj_p1){
    return 2./(xj_p1-xj_m1) * ((uj_p1-uj)/(xj_p1-xj) - (uj-uj_m1)/(xj-xj_m1));
}

HOST_DEVICE static inline myfloat d2u_dx2Left(myfloat u1, myfloat u2, myfloat u3,
                                  myfloat x1, myfloat x2, myfloat x3){
    return 2*((x2-x1)*(u3-u1) - (x3-x1)*(u2-u1)) / ((x3-x1)*(x2-x1)*(x3-x2));
}

HOST_DEVICE static inline myfloat d2u_dx2Right(myfloat uN, myfloat uN1, myfloat uN2,
                                   myfloat xN, myfloat xN1, myfloat xN2){
    return 2*((xN1-xN)*(uN2-uN) - (xN2-xN)*(uN1-uN)) / ((xN2-xN)*(xN1-xN)*(xN2-xN1));
}

HOST_DEVICE static inline myfloat dxj_dt(myfloat rhoj_m1, myfloat rhoj, myfloat rhoj_p1,
                             myfloat xj_m1,   myfloat xj,   myfloat xj_p1,
               myfloat delta_xi, myfloat tau){ // bigger tau slower movement
    return 0.5*((rhoj_p1+rhoj)*(xj_p1-xj) - (rhoj+rhoj_m1)*(xj-xj_m1)) / (rhoj*tau*delta_xi*delta_xi);
}

HOST_DEVICE static inline myfloat dxBoundary_dt(){
    return 0;
}

HOST_DEVICE static inline myfloat rhoj(myfloat alha_h, myfloat d2u_dx2){
    return MYCBRT(1+ MYPOW(MYABS(d2u_dx2), 2));
}

HOST_DEVICE static inline myfloat alphaHSumElems(myfloat d2h_dx2_j, myfloat d2h_dx2_jm1, myfloat xj, myfloat xj_m1){
    return 0.5*(xj-xj_m1)*MYPOW(MYPOW(d2h_dx2_j, 2./3.)+MYPOW(d2h_dx2_jm1, 2./3.), 3);
}

// dx = xj+1-xj-1
HOST_DEVICE static inline myfloat RHS32(myfloat dxj_dt, myfloat dx, myfloat hj_p1, myfloat hj_m1){
    return (hj_p1-hj_m1)/dx * dxj_dt - (MYPOW(hj_p1, 3./2.)-MYPOW(hj_m1, 3./2.))/dx;
}

int max(int A[], int i, int j)
{
    int max_val = A[0];

    #pragma omp parallel for reduction(max:max_val) 
    for (int idx = i; idx < j; idx++)
       max_val = max_val > A[idx] ? max_val : A[idx];

    return max_val;
}

#define TAU 0.01
int mainAdaptiveh23Test(int argc, char* argv[]){
    assert(argc > 1);
    int N_grid = atoi(argv[1]);
    myfloat art_viscosity = 1e-5; // artificial viscosity
    testOpenMP(16);
    myfloat x_min = -2, x_max = 16, end_t = 1, t = 0;
    myfloat delta_xi = 1./(N_grid-1);

    int len_ts = 300;
    fvec_t* timesteps = linspaceFvecIncludeMax(NULL, 0, end_t, len_ts);
    printFVec(timesteps);
    TimeSaver_t* timesteps_to_save = &(TimeSaver_t){.current_idx = 0,
                                                    .len = len_ts,
                                                    .to_save = timesteps->vals};

    fvec_t *h         = allocateFVec(NULL, N_grid);
    fvec_t *rho       = allocateFVec(NULL, N_grid);
    fvec_t *rho_tmp   = allocateFVec(NULL, N_grid);
    fvec_t *d2h_dx2 = allocateFVec(NULL, N_grid);
    fvec_t *x  = linspaceFvecIncludeMax(NULL, x_min, x_max, N_grid);
    
    for(int i = 0; i < x->len; i++){
        h->vals[i] = 0.5 + exp(-pow(x->vals[i], 2));
    }
    fvec_t *h_tmp = newCopyFVec(NULL, h);

    fvec_t* advective_speed = newCopyFVec(NULL, h);
    FOR_EACH_VEC(advective_speed, i) =  3./2. * sqrt(advective_speed->vals[i]); // i became my own nightmare TODO better func 

    myfloat dt = 1;
    //  = CFL*/getMaxAbsFvec(advective_speed); 
    #pragma omp parallel for reduction(min:dt) 
    for (int i = 1; i<x->len-1; i++)
        dt = MIN(dt, (x->vals[i+1]-x->vals[i-1])*0.5/(MYABS(advective_speed->vals[i])));
    printf("DT=%lf\n", dt);

    {
        myfloat alpha_h = 0;
        for(int j = 1; j < h->len-1; j++){
            alpha_h += 0.5*(x->vals[j]-x->vals[j-1])*MYPOW(MYPOW(MYABS(d2h_dx2->vals[j]  ), 2./3.)+
            MYPOW(MYABS(d2h_dx2->vals[j-1]), 2./3.), 3);
        }
        alpha_h = MAX(1, alpha_h);
        #pragma omp parallel for
        for(int j = 0; j < h->len; j++){
            rho->vals[j] = rhoj(alpha_h, d2h_dx2->vals[j]);
        }
    }

    while(t<end_t){

        dt = 1.;
        for (int i = 1; i<x->len-1; i++)
            dt = MIN(dt, (x->vals[i+1]-x->vals[i-1])*0.5/(MYABS(advective_speed->vals[i])));
        printf("DT=%lf\n", dt);
        

        // solution
        d2h_dx2->vals[0] = d2u_dx2Left(h->vals[0], h->vals[1], h->vals[2],
                                       x->vals[0], x->vals[1], x->vals[2]);
        #pragma omp parallel for
        for(int j = 1; j < h->len-1; j++){
            d2h_dx2->vals[j] = d2u_dx2(h->vals[j-1], h->vals[j], h->vals[j+1],
                                         x->vals[j-1], x->vals[j], x->vals[j+1]);
            myfloat dx = x->vals[j+1] - x->vals[j-1];
            myfloat dxjdt = dxj_dt(rho->vals[j-1], rho->vals[j], rho->vals[j+1],
                                       x->vals[j-1], x->vals[j], x->vals[j+1], delta_xi, TAU);
            h_tmp->vals[j] = RHS32(dxjdt, dx, h->vals[j+1], h->vals[j-1])*dt + h->vals[j] + (art_viscosity * d2h_dx2->vals[j] * dt);
        }
        d2h_dx2->vals[N_grid-1] = d2u_dx2Right(h->vals[N_grid-1], h->vals[N_grid-2], h->vals[N_grid-3],
                                               x->vals[N_grid-1], x->vals[N_grid-2], x->vals[N_grid-3]);
    
        // new x 
        #pragma omp parallel for
        for(int j = 1; j < x->len-1; j++){
            // euler time stepping
            x->vals[j] += dxj_dt(rho->vals[j-1], rho->vals[j], rho->vals[j+1],
                                x->vals[j-1], x->vals[j], x->vals[j+1], delta_xi, TAU)*dt ;
        }

        // mesh
        myfloat alpha_h = 0;
        for(int j = 1; j < h->len-1; j++){
            alpha_h += 0.5*(x->vals[j]-x->vals[j-1])*MYPOW(MYPOW(MYABS(d2h_dx2->vals[j]  ), 2./3.)+
            MYPOW(MYABS(d2h_dx2->vals[j-1]), 2./3.), 3);
        }
        alpha_h = MAX(1, alpha_h);
        
        #pragma omp parallel for
        for(int j = 0; j < h->len; j++){
            rho->vals[j] = rhoj(alpha_h, d2h_dx2->vals[j]);
        }

        //making new rho smoothing 
        int nb_outer_rho_avg = 4;
        for(int i = 0; i < nb_outer_rho_avg; i++){
            rho_tmp->vals[0] = 0.5*(rho->vals[0] + rho->vals[1]); 
            #pragma omp parallel for
            for(int j = 1; j < rho_tmp->len-1; j++){
                rho_tmp->vals[j] = 0.25*(rho->vals[j-1] + rho->vals[j+1]) + 0.5*rho->vals[j];
            }
            rho_tmp->vals[N_grid-1] = 0.5*(rho->vals[N_grid-1] + rho->vals[N_grid-2]); 
            SWAP(rho, rho_tmp, fvec_t*);
        }

   

        char filename[64];
        sprintf(filename, "Out/");
        SWAP(h, h_tmp, fvec_t*);
        int saved = saveTimesteps(timesteps_to_save, t, h, dt, filename);
        if(saved){
            char full_filename[MAX_FILENAME_LEN+64];// the +64 is for time storing
            sprintf(full_filename, "%sx/%lf.dat", filename, t);
            saveVecToFile(x, full_filename);
        }
        t += dt;

       
    }

    // printFVec(h);
    return 0;
}


struct MeshParams {
    int N_grid;
    myfloat delta_xi;
    myfloat tau;
    myfloat art_viscosity;
};

int movingMeshRhs(double t, const double y[], double f[], void *params) {
    struct MeshParams *p = (struct MeshParams *)params;
    int N = p->N_grid;
    myfloat delta_xi = p->delta_xi;
    myfloat tau = p->tau;
    myfloat art_viscosity = p->art_viscosity;

    const double *h = y;
    const double *x = y + N;
    
    double *dh_dt = f;
    double *dx_dt = f + N;

    double d2h_dx2[N];
    double rho[N];
    double rho_tmp[N];

    d2h_dx2[0] = d2u_dx2Left(h[0], h[1], h[2], x[0], x[1], x[2]);
    #pragma omp parallel for
    for(int j = 1; j < N - 1; j++){
        d2h_dx2[j] = d2u_dx2(h[j-1], h[j], h[j+1], x[j-1], x[j], x[j+1]);
    }
    d2h_dx2[N-1] = d2u_dx2Right(h[N-1], h[N-2], h[N-3], x[N-1], x[N-2], x[N-3]);

    // rho
    myfloat alpha_h = 0;
    #pragma omp parallel for
    for(int j = 1; j < N - 1; j++){
        alpha_h += 0.5*(x[j]-x[j-1]) * MYPOW(MYPOW(MYABS(d2h_dx2[j]), 2./3.) + 
                                             MYPOW(MYABS(d2h_dx2[j-1]), 2./3.), 3);
    }
    alpha_h = MAX(1, alpha_h);

    #pragma omp parallel for
    for(int j = 0; j < N; j++){
        rho[j] = rhoj(alpha_h, d2h_dx2[j]);
    }

    // Smoothing rho
    int nb_outer_rho_avg = 4;
    for(int i = 0; i < nb_outer_rho_avg; i++){
        rho_tmp[0] = 0.5*(rho[0] + rho[1]); 
        #pragma omp parallel for
        for(int j = 1; j < N - 1; j++){
            rho_tmp[j] = 0.25*(rho[j-1] + rho[j+1]) + 0.5*rho[j];
        }
        rho_tmp[N-1] = 0.5*(rho[N-1] + rho[N-2]); 
        for(int j = 0; j < N; j++) rho[j] = rho_tmp[j];
    }

    
    // BC
    dx_dt[0] = 0.0;
    dh_dt[0] = 0.0;
    dx_dt[N-1] = 0.0;
    dh_dt[N-1] = 0.0;

    // Interior 
    #pragma omp parallel for
    for(int j = 1; j < N - 1; j++){
        myfloat dx_val = x[j+1] - x[j-1];
        myfloat dxjdt = dxj_dt(rho[j-1], rho[j], rho[j+1], 
                               x[j-1], x[j], x[j+1], delta_xi, tau);
        
        dx_dt[j] = dxjdt;
        dh_dt[j] = RHS32(dxjdt, dx_val, h[j+1], h[j-1]) + (art_viscosity * d2h_dx2[j]);
    }

    return GSL_SUCCESS;
}

// only used if  gsl_odeiv2_step_msbdf
int movingMeshJac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    struct MeshParams *p = (struct MeshParams *)params;
    int N_total = 2 * p->N_grid; 
    
    double f0[N_total];
    double f1[N_total];
    double y_pert[N_total];
    
    movingMeshRhs(t, y, f0, params);
    
    double eps = 1e-8; 
    
    for (int j = 0; j < N_total; j++) {
        for (int i = 0; i < N_total; i++) {
            y_pert[i] = y[i];
        }
        
        y_pert[j] += eps;
        
        movingMeshRhs(t, y_pert, f1, params);
        
        for (int i = 0; i < N_total; i++) {
            dfdy[i * N_total + j] = (f1[i] - f0[i]) / eps; // 1D row majro 
        }
    }
    // no explicit dependency in t     
    for (int i = 0; i < N_total; i++) {
        dfdt[i] = 0.0;
    }
    
    return GSL_SUCCESS;
}
int mainAdaptiveGSL(int argc, char* argv[]){
    assert(argc > 1);
    int N_grid = atoi(argv[1]);
    myfloat art_viscosity = 1e-5;
    myfloat x_min = -8, x_max = 16, end_t = 10;
    myfloat delta_xi = 1./(N_grid-1);

    int len_ts = 300;
    fvec_t* timesteps = linspaceFvecIncludeMax(NULL, 0, end_t, len_ts);
    
    fvec_t *h = allocateFVec(NULL, N_grid);
    fvec_t *x = linspaceFvecIncludeMax(NULL, x_min, x_max, N_grid);
    
    // IC
    #pragma omp parallel for
    for(int i = 0; i < x->len; i++){
        h->vals[i] = 0.5 + exp(-pow(x->vals[i], 2));
    }

    
    double *y = malloc(2 * N_grid * sizeof(double));
    for(int i = 0; i < N_grid; i++) {
        y[i] = h->vals[i];             
        y[N_grid + i] = x->vals[i];    
    }

    struct MeshParams params = {N_grid, delta_xi, TAU, art_viscosity};
    
    gsl_odeiv2_system sys = {movingMeshRhs, movingMeshJac, 2 * N_grid, &params};
    
    // automatically adjusts dt 
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf,//gsl_odeiv2_step_rk8pd,  // or gsl_odeiv2_step_msbdf
                                                          1e-6, 1e-6, 0.0);

    double t = 0.0;
    
    for (int i = 0; i < len_ts; i++) {
        double ti = timesteps->vals[i];
        
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf("GSL Error: Return value = %d\n", status);
            break;
        }

        #pragma omp parallel for
        for(int j = 0; j < N_grid; j++) {
            h->vals[j] = y[j];
            x->vals[j] = y[N_grid + j];
        }

        char filename[128];
        sprintf(filename, "Out/%lf.dat", t);
        saveVecToFile(h, filename);
        
        sprintf(filename, "Out/x/%lf.dat", t);
        saveVecToFile(x, filename);
        
        printf("Time t = %lf saved.\n", t);
    }

    gsl_odeiv2_driver_free(d);
    free(y);
    // TODO leaking my vecs but since afterwards the program closes its ok
    return 0;
}

int mainOld(int argc, char* argv[]);

int main(int argc, char* argv[]){
    mainAdaptiveGSL(argc, argv);
    // mainOld(argc, argv);
    return 0;
}

int mainOld(int argc, char* argv[]){
    assert(argc > 1);
    int N_grid = atoi(argv[1]);
    testOpenMP(16);
    myfloat x_min = -8, x_max = 32, end_t = 80, t = 0;
    fvec_t *h  = allocateFVec(NULL, N_grid);
    fvec_t *q  = allocateFVec(NULL, N_grid);
    fvec_t *Z  = allocateFVec(NULL, N_grid); // bottom
    fvec_t *x  = linspaceFvecIncludeMax(NULL, x_min, x_max, N_grid);

    fvec_t *dx = newCopyFVec(NULL, x);

    for(int i = 0; i < dx->len-1; i++){
        dx->vals[i] -= dx->vals[i+1];
        dx->vals[i] *= -1;
    } dx->vals[dx->len-1] = dx->vals[dx->len-2];

    myfloat DX = dx->vals[0];

    FOR_EACH_VEC(h, i) = (x->vals[i] < 0)*0.75 ;//+ (x->vals[i] < 4)*0.25;
    FOR_EACH_VEC(Z, i) = (x->vals[i] > 5)*((x->vals[i] - 5)*0.05);
    // for(int i = 0; i < x->len; i++){
    //     h->vals[i] = 0.5 + exp(-pow(x->vals[i], 2));
    // }

    
    myfloat dt = CFL * getMaxOfVec(dx) / maxAbsEigenvalue(h, q);
    printf("dx=%lf, dt=%lf\n", DX, dt);

    int len_ts = 300;
    fvec_t* timesteps = linspaceFvecIncludeMax(NULL, 0, end_t, len_ts);
    printFVec(timesteps);
    TimeSaver_t* timesteps_to_save = &(TimeSaver_t){.current_idx = 0,
                                                    .len = len_ts,
                                                    .to_save = timesteps->vals};
    printf("Timesteps to save (t_end = %lf) : ", end_t); 
    printFArrValues(timesteps_to_save->to_save, timesteps_to_save->len);

    saveVecToFile(Z, "z.dat");
    saveVecToFile(x, "x.dat");
    rusanovMethod(h, q, Z, dt, end_t, t, DX, timesteps_to_save);
    return 0;
}
