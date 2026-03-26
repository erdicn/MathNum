#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h> 

#include <math.h>
#define EPS  1e-7 // TODO de faire en tel sort que on mets ca dans chaque fonction par rapport a leur errors 

typedef struct Polynome{
    int len;
    double* func;
}poly;

void InitPolyToZero(poly* P){
    for(int i = 0; i < P->len; i++){
        P->func[i] = 0;
    }
}

poly* PolyFilledWithZeros(int len){
    poly* new_P = malloc(sizeof(poly));
    new_P->len = len;
    new_P->func = malloc(sizeof(double)*len);
    InitPolyToZero(new_P);
    return new_P;
}

poly* InitPoly(int len, double* func){
    poly* new_P = malloc(sizeof(poly));
    new_P->len = len;
    new_P->func = malloc(sizeof(double)*len); // We redo a malloc to not change the func data entered 
    for(int i = 0; i < len; i++){
        new_P->func[i] = func[i];
    }
    return new_P;
}

void RemoveZerosAtEnd(poly* P){
    double eps = EPS;
    int zeros = 0;
    for(int i = P->len - 1; i > 0; i--){
        if (P->func[i] > 0 - eps && P->func[i] < 0 + eps){
            zeros++;
        } else {
            break;
        }
    }
    if (zeros > 0){
        P->len = P->len - zeros;
        P->func = realloc( P->func, sizeof(double) * (P->len));
    }
}

int TestPolyEquality(poly *P1, poly *P2){
    if (P1->len != P2->len) {return false;}
    for(int i = 0; i < P1->len; i++){
        if (P1->func[i] != P2->func[i]){
            return false;
        }
    }
    return true;
}

void FreePoly(poly* P){
    free(P->func);
    free(P);
}

void PrintPoly(poly P){
    if (P.len > 0){
        for(int i = 0;  i < P.len-1; i++){
            printf("%.2lfx^%d + ", P.func[i], i);
        }
        printf("%.2lfx^%d", P.func[P.len-1], P.len-1);
        puts("");
    }
}

// void PolyInitToPoly(poly* P, poly* P_egalise){

// }

poly* PolyMultConst(double c, poly* P){
    for(int i = 0; i < P->len; i++){
        P->func[i] = P->func[i] * c;
    }
    return P;
}

poly* PolyDivConst(double c, poly* P){
    return PolyMultConst(1./c, P);    
}


poly* PolyMult(poly* P1, poly* P2){
    poly* new_P = malloc(sizeof(poly));
    new_P->len = P1->len + P2->len;
    new_P->func= malloc(sizeof(double)*(new_P->len));
    for(int i = 0; i < P1->len; i++){
        for(int j = 0; j < P2->len; j++){
            new_P->func[i+j] += P1->func[i] * P2->func[j];
        }
    }
    //RemoveZerosAtEnd(new_P); // TODO not sure if i need it or not
    return new_P;
}

poly* PolyAdd(poly* P1, poly* P2){
    poly* grand_P, * small_P; 
    if (P1->len > P2->len){
        grand_P = P1;
        small_P = P2;
    } else{
        grand_P = P2;
        small_P = P1;
    }
    poly* new_P = malloc(sizeof(poly));
    new_P->len = grand_P->len;
    new_P->func = malloc(sizeof(double)*new_P->len);
    for(int i = 0; i < new_P->len; i++){
        if (i < small_P->len){
            new_P->func[i] = small_P->func[i] + grand_P->func[i];
        } else {
            new_P->func[i] = grand_P->func[i];
        }
    }
    return new_P;
}


poly* LangrangePolynomial(int len_data, double* x_vals, double* f_vals){
    poly* P = PolyFilledWithZeros(len_data);
    poly* old_P, *old_L;
    double one[1] = {1};
    poly* L = InitPoly(1, one);
    double x_xj[2] = {0, 1}; // we initialise the 0 in the for
    poly* P_x_xj = InitPoly(2, x_xj);
    for(int i = 0; i < len_data; i++){
        L = InitPoly(1, one);
        for(int j = 0; j < len_data; j++){
            if (i == j) continue;
            old_L = L;
            P_x_xj->func[0] = -x_vals[j]; 
            L = PolyDivConst((x_vals[i] - x_vals[j]), PolyMult(L, P_x_xj));
            free(old_L);
        }
        old_P = P;
        P = PolyAdd(P, PolyMultConst(f_vals[i], L));
        free(old_P);
    }
    RemoveZerosAtEnd(P);
    return P;
}

typedef struct Matrix{
    int len_rows, len_cols;
    double** mat;
    int* before_pivots;
}matrix;

void FreeMatrix(matrix* M){
    // I am not free maybe the matrix can be, there is still hope for it unlike...
    free(M->before_pivots);
    for(int i = 0; i < M->len_rows; i++){
        free(M->mat[i]);
    }
    free(M->mat);
    free(M);
}

void SwapLines(matrix* M, int a, int b){
    double* tmp = M->mat[a];
    M->mat[a] = M->mat[b];
    M->mat[b] = tmp;
    int tmp_b_piv = M->before_pivots[b];
    M->before_pivots[b] = M->before_pivots[a];
    M->before_pivots[a] = tmp_b_piv;
}

// Uses buble sort to sort matrix acording to nb of zeros before pivots
void SortMatrixAcordingToPivots(matrix* M){
    int n =  M->len_rows; 
    for(int i = 0; i <n-1; i++){
        for(int j = 0; j < n-i-1; j++){
            if(M->before_pivots[j] > M->before_pivots[j+1]){
                SwapLines(M, j, j+1);
            }
        }
    }
}

void PrintMatrix(matrix M){
    for(int i = 0; i < M.len_rows; i++){
        for(int j = 0; j < M.len_cols; j++){
            printf("%.2lf ", M.mat[i][j]);
        }
        puts("");
    }
    if (M.before_pivots != NULL){
        for(int i = 0; i < M.len_rows; i++)
            printf("%d ", M.before_pivots[i]);
        puts("");
    }
}

void GaussLineSubstract(matrix* M, int substracting_line_index, int substractor_line_index){
    int zeros_before_pivot = 0;
    bool zeros_before = true;
    for(int i = 0; i < M->len_cols; i++){
        M->mat[substracting_line_index][i] = M->mat[substracting_line_index][i] - M->mat[substractor_line_index][i];
        if (M->mat[substracting_line_index][i] > 0 - EPS && 
                M->mat[substracting_line_index][i] < 0 + EPS &&
                    zeros_before){
            zeros_before_pivot++;
        } else {
            zeros_before = false;
        }
    }
    M->before_pivots[substracting_line_index] = zeros_before_pivot;
}

void GaussLineMultiply(matrix* M, int line_index, double multiplier){
    int zeros_before_pivot = 0;
    bool zeros_before = true;
    for(int i = 0; i < M->len_cols; i++){
        M->mat[line_index][i] = M->mat[line_index][i]*multiplier;
        if (M->mat[line_index][i] > 0 - EPS && 
                M->mat[line_index][i] < 0 + EPS &&
                    zeros_before){
            zeros_before_pivot++;
        } else {
            zeros_before = false;
        }
    }
    M->before_pivots[line_index] = zeros_before_pivot;
}

poly* ResolutionGauss(matrix M){
    // TODO maybe automatic matrix filler ?
    matrix* sol_M = malloc(sizeof(matrix));
    sol_M->len_rows = M.len_rows;
    sol_M->len_cols = M.len_cols;
    sol_M->before_pivots = malloc(sizeof(int)*(sol_M->len_rows));
    sol_M->mat = malloc(sizeof(double*)*(sol_M->len_rows));
    for(int i = 0; i < sol_M->len_rows; i++){sol_M->before_pivots[i] = 0;}
    for(int i = 0; i < sol_M->len_rows; i++){
        bool zero_counter = true;
        sol_M->mat[i] = malloc(sizeof(double)*(sol_M->len_cols));
        int j;
        for(j = 0; j < M.len_cols-1; j++){  // -1 car la case finale cest le resultat
            if (zero_counter && M.mat[i][j] == 0){
                sol_M->before_pivots[i] += 1;
            } else if (zero_counter){
                zero_counter = false;
            }
            sol_M->mat[i][j] = M.mat[i][j];
        }
        sol_M->mat[i][j] = M.mat[i][j];
    }
    SortMatrixAcordingToPivots(sol_M);
    // PrintMatrix(*sol_M);
    
    double multiplier;
    int pivot;
    for(int i = 0; i < sol_M->len_rows; i++){
        for(int j = 0; j < sol_M->len_rows; j++){
            if (sol_M->before_pivots[i] != sol_M->before_pivots[j] || j == i)
                continue;
            // PrintMatrix(*sol_M);
            pivot = sol_M->before_pivots[i];
            // printf("i=%d, j=%d, pivot=%d \n",i,j, pivot);
            multiplier = sol_M->mat[j][pivot] / sol_M->mat[i][pivot];
            GaussLineMultiply(sol_M, i, multiplier);
            GaussLineSubstract(sol_M, j, i);
        }
    }
    SortMatrixAcordingToPivots(sol_M);

    // PrintMatrix(*sol_M);
    // Error catching
    int before = -1;
    for(int i = 0; i < sol_M->len_rows; i++){
        if(sol_M->before_pivots[i] != before +1){
            printf("Error the gauss algorithm does not work properly or the matrix does not have a solurion probaply need a more thourough look at above\n");
        }
        before++;
    }
    if(sol_M->before_pivots[sol_M->len_rows-1] != sol_M->len_rows - 1){
        printf("solution does not exist\n");
        poly* error = malloc(sizeof(poly));
        error->len = 0;
        error->func = NULL;
        return error;
    }

    for(int i = sol_M->len_rows - 1; i > 0; i--){
        for(int j = i - 1; j >= 0; j--){
            multiplier = (double)sol_M->mat[j][i] / (double)sol_M->mat[i][i];
            GaussLineMultiply(sol_M, i, multiplier);
            GaussLineSubstract(sol_M, j, i);
        }
    }
    poly* solutions = malloc(sizeof(poly));
    solutions->len = sol_M->len_rows;
    solutions->func = malloc(sizeof(double)*solutions->len);
    for(int i = sol_M->len_rows - 1; i >= 0; i--){
        solutions->func[i] = sol_M->mat[i][sol_M->len_cols-1] / sol_M->mat[i][i];
    }
    // PrintMatrix(*sol_M);
    //FreeMatrix(sol_M);
    return solutions;
}

matrix* ReadMatrix(){
    FILE* file = fopen("matrix.txt", "r");
    if (file == NULL){
        printf("Errod matrix.txt is NULL\n");
        matrix* err = malloc(sizeof(matrix));
        err->len_rows = 0;
        err->len_cols = 0;
        return err;
    }
    matrix* M = malloc(sizeof(matrix));
    int nb_rows, nb_cols;
    fscanf(file, "%d %d", &nb_rows, &nb_cols);
    M->len_rows = nb_rows;
    M->len_cols = nb_cols;
    M->mat = malloc(sizeof(double*)*M->len_rows);
    // printf("%d %d\n", M->len_rows, M->len_cols);
    double readed = .13;
    bool zeros = true;
    int zero_counter = 0;
    M->before_pivots = malloc(sizeof(double)*M->len_rows);
    for(int i = 0; i < M->len_rows; i++){
        zero_counter = 0;
        zeros = true;
        M->mat[i] = malloc(sizeof(double)*M->len_cols);
        for(int j = 0; j < M->len_cols; j++){
            fscanf(file, "%lf", &readed);
            M->mat[i][j] = readed;
            // printf("%.0lf ", M->mat[i][j]);
            if(M->mat[i][j] > 0 - EPS && M->mat[i][j] < 0 + EPS && zeros){
                zero_counter++;
            } else{
                zeros = false;
            }
        }
        M->before_pivots[i] = zero_counter;
        // puts("");
    }
    fclose(file);
    return M;
}

poly* MoindresCarees(int len_data, double* x_vals, double* f_vals, int degree){
    if (len_data < degree + 2){
        printf("Pas assez de points pour faire le polynome degree voulu(%d) max degree possible(%d)\n",degree,len_data-2 );
        degree = len_data - 2;
    }
    matrix* Mat_to_solve = malloc(sizeof(matrix));
    Mat_to_solve->len_rows = degree + 1; // +1 car on veux inclure le degree
    Mat_to_solve->len_cols = degree + 2; // +1 pour le degree et +1 pour le resultat = +2
    Mat_to_solve->mat = malloc(sizeof(double*)*Mat_to_solve->len_rows);
    Mat_to_solve->before_pivots = NULL;
    // pour avoir les memes notations que le cours 
    int m = degree;
    //int n = len_data - 1;
    
    // On precalcule nos valeurs comme ca on a pas besoin de les recalculer
    double* list_with_precalculations = malloc(sizeof(double)*(2*m+1));
    for(int i = 0; i < (2*m+1); i++){
        list_with_precalculations[i] = 0;
        // printf("\n%d precalc %lf\n",i, list_with_precalculations[i]);
        for(int j = 0; j < len_data; j++){
            list_with_precalculations[i] += pow(x_vals[j], i);
            // printf("%lf ",pow(x_vals[i], i));
        }
        // printf("%d precalc %lf\n",i, list_with_precalculations[i]);
    }

    double fw_xw_multiplicatio_sum = 0;
    for(int i = 0; i < Mat_to_solve->len_rows; i++){
        Mat_to_solve->mat[i] = malloc(sizeof(double)*Mat_to_solve->len_cols);
        for(int j = 0; j < Mat_to_solve->len_rows; j++){
            Mat_to_solve->mat[i][j] = list_with_precalculations[i+j];
            // printf("%lf \n" ,list_with_precalculations[i+j]);
        }
        fw_xw_multiplicatio_sum = 0;
        for(int w = 0; w < len_data; w++){
            fw_xw_multiplicatio_sum += f_vals[w]*pow(x_vals[w], i);
        }
        // printf("%d %lf\n", i, fw_xw_multiplicatio_sum);
        Mat_to_solve->mat[i][Mat_to_solve->len_cols-1] = fw_xw_multiplicatio_sum;
    }
    // puts("");
    // PrintMatrix(*Mat_to_solve);
    // puts("");
    poly* solution = ResolutionGauss(*Mat_to_solve);
    // PrintMatrix(*Mat_to_solve);
    // FreeMatrix(Mat_to_solve);
    return solution;
}

// TODO faire lerreur

// TODO faire methode miniMax

// TODO all the frees 

void TestMoindresCarees();
void TestingLangrange();
void TestGauss();
 
int main(){
    printf("\nMoindre careees:\n"); TestMoindresCarees();
    printf("\nLangrange:\n");TestingLangrange();
    printf("\nGauss:\n");TestGauss();
    return 0;
}

void TestMoindresCarees(){
    double x_vals[] = {-2, 0, 4, 6};
    double f_vals[] = {3, 5, 8, 5};
    poly* moindre_0 = MoindresCarees(4, x_vals, f_vals, 0);
    PrintPoly(*moindre_0);
    poly* moindre_1 = MoindresCarees(4, x_vals, f_vals, 1); // TODO why does it gives -nan ???
    PrintPoly(*moindre_1);
    // TODO fix why there is a bug here double free or corruption (out) 
    poly* moindre_2 = MoindresCarees(4, x_vals, f_vals, 2);
    PrintPoly(*moindre_2);
    poly* moindre_3 = MoindresCarees(4, x_vals, f_vals, 3);
    PrintPoly(*moindre_3);
}

void TestGauss(){
    matrix* test_TD = ReadMatrix();
    // PrintMatrix(*test_TD);
    poly* sol = ResolutionGauss(*test_TD);
    PrintPoly(*sol);
}

void TestingLangrange(){
    poly testP;
    testP.len = 4;
    testP.func = malloc(sizeof(double)*testP.len);
    for(int i = 0; i < testP.len; i++){
        testP.func[i] = i+1;
    }
    poly* added2testP = PolyAdd(&testP, &testP);
    poly* multiplied2testP = PolyMultConst(2, &testP);
    PrintPoly(*added2testP);
    assert(TestPolyEquality(multiplied2testP, added2testP));
    FreePoly(added2testP);

    int f_len = 4;
    double f[4] = {1, 2, 0, 0}; 
    poly* P1 = InitPoly(f_len, f);
    printf("Before erasig the zeros ");PrintPoly(*P1);
    RemoveZerosAtEnd(P1);
    printf("After erasig the zeros  ");PrintPoly(*P1);
    double test_x[4] = {-2, 0, 4, 6};
    double test_f[4] = { 3, 5, 8, 5}; 
    PrintPoly(*LangrangePolynomial(4, test_x, test_f));

    double test_x2[24], test_f2[24];
    for(int i = 0; i < 24; i++){
        test_x2[i] = i;
        test_f2[i] = i*i;
    }
    printf("For y = x^2 ");
    PrintPoly(*LangrangePolynomial(24, test_x2, test_f2));
}