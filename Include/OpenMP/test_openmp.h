#ifndef TEST_OPENMP_H
#define TEST_OPENMP_H

#ifdef _OPENMP
#   include <omp.h>
#endif

void testOpenMP(int nb_cores){
    #pragma omp parallel for
    for(int i = 0; i < nb_cores; i++){
        printf("Thread %d of %d is processing index %d\n", omp_get_thread_num(), omp_get_num_threads(), i);
    }
}

#endif /* TEST_OPENMP_H */
