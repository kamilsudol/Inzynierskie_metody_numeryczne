#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

int l_function(int i, int j, int nx){
    return i+j*(nx+1);
}

double laplasjan(gsl_vector *T, int l, double delta, int nx){
    double lap = (gsl_vector_get(T, l+1) - 2*gsl_vector_get(T, l) + gsl_vector_get(T, l-1))/pow(delta, 2) + (gsl_vector_get(T, l+nx+1) - 2*gsl_vector_get(T, l) + gsl_vector_get(T, l-nx-1))/pow(delta, 2);
    return lap;
}

void wypelnij1(gsl_matrix *A, gsl_matrix *B, gsl_vector *c, double dt, double delta, double kb, double kd, double tb, double td, int nx, int ny, int n){
    for(int i = 1; i < nx; i++){
        for(int j = 1; j < ny; j++){
            int l = l_function(i, j, nx);
            gsl_matrix_set(A, l, l - nx - 1, dt/(2*pow(delta, 2)));
            gsl_matrix_set(A, l, l - 1, dt/(2*pow(delta, 2)));
            gsl_matrix_set(A, l, l + 1, dt/(2*pow(delta, 2)));
            gsl_matrix_set(A, l, l + nx + 1, dt/(2*pow(delta, 2)));
            gsl_matrix_set(A, l, l, -dt/(2*pow(delta, 2))-1);

            gsl_matrix_set(B, l, l - nx - 1, -dt/(2*pow(delta, 2)));
            gsl_matrix_set(B, l, l - 1, -dt/(2*pow(delta, 2)));
            gsl_matrix_set(B, l, l + 1, -dt/(2*pow(delta, 2)));
            gsl_matrix_set(B, l, l + nx + 1, -dt/(2*pow(delta, 2)));
            gsl_matrix_set(B, l, l, dt/(2*pow(delta, 2))-1);
        }
    }

    for(int j = 0; j <= ny; j++){
        int l = l_function(0, j, nx);
        gsl_matrix_set(A, l, l, 1);
        gsl_matrix_set(B, l, l, 1);
        gsl_vector_set(c, l, 0);
    }

    for(int j = 0; j <= ny; j++){
        int l = l_function(nx, j, nx);
        gsl_matrix_set(A, l, l, 1);
        gsl_matrix_set(B, l, l, 1);
        gsl_vector_set(c, l, 0);
    }

    for(int i = 1; i < nx; i++){
        int l = l_function(i, ny, nx);
        gsl_matrix_set(A, l, l - nx - 1, -1/(kb*delta));
        gsl_matrix_set(A, l, l, 1 + 1/(kb*delta));
        gsl_vector_set(c, l, tb);
    }

    for(int i = 1; i < nx; i++){
        int l = l_function(i, 0, nx);
        gsl_matrix_set(A, l, l + nx + 1, -1/(kd*delta));
        gsl_matrix_set(A, l, l, 1 + 1/(kd*delta));
        gsl_vector_set(c, l, td);
    }

    for(int k = 0; k < n; k++){
        for(int i = 1; i < nx; i++){
            int l = l_function(i, ny, nx);
            gsl_matrix_set(B, l, k, 0);
        }
    }

    for(int k = 0; k < n; k++){
        for(int i = 1; i < nx; i++){
            int l = l_function(i, 0, nx);
            gsl_matrix_set(B, l, k, 0);
        }
    }

}

void wypelnij2(gsl_vector *T, double ta, double tc, int nx, int ny){
    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            int l = l_function(i, j, nx);
            if(i == 0){
                gsl_vector_set(T, l, ta);
            }else if(i == nx){
                gsl_vector_set(T, l, tc);
            }else{
                gsl_vector_set(T, l, 0);
            }
        }
    }
}

void cn(gsl_matrix *A, gsl_matrix *B, gsl_vector *c, gsl_vector *d, gsl_vector *T, gsl_permutation *perm, int *sigma, int nx, int ny, double delta, int IT_MAX){
    FILE * file_ptr1 = fopen("temperatura1.out", "w");
    FILE * file_ptr2 = fopen("temperatura2.out", "w");
    FILE * file_ptr3 = fopen("temperatura3.out", "w");
    FILE * file_ptr4 = fopen("temperatura4.out", "w");
    FILE * file_ptr5 = fopen("temperatura5.out", "w");

    FILE * file_ptr6 = fopen("laplasjan1.out", "w");
    FILE * file_ptr7 = fopen("laplasjan2.out", "w");
    FILE * file_ptr8 = fopen("laplasjan3.out", "w");
    FILE * file_ptr9 = fopen("laplasjan4.out", "w");
    FILE * file_ptr10 = fopen("laplasjan5.out", "w");

    gsl_linalg_LU_decomp(A, perm, sigma);
    
    for(int it = 0; it < IT_MAX; it++){
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0.0, d);
        gsl_blas_daxpy(1.0, c, d);
        gsl_linalg_LU_solve(A, perm, d, T);
        if(it == 100){
            for(int i = 1; i < nx; i++){
                for(int j = 1; j < ny; j++){
                    fprintf(file_ptr1, "%d %d %f\n", i, j, gsl_vector_get(T, l_function(i, j, nx)));
                    fprintf(file_ptr6, "%d %d %f\n", i, j, laplasjan(T, l_function(i, j, nx), delta, nx));
                }
            }
        }else if(it == 200){
            for(int i = 1; i < nx; i++){
                for(int j = 1; j < ny; j++){
                    fprintf(file_ptr2, "%d %d %f\n", i, j, gsl_vector_get(T, l_function(i, j, nx)));
                    fprintf(file_ptr7, "%d %d %f\n", i, j, laplasjan(T, l_function(i, j, nx), delta, nx));
                }
            }
        }else if(it == 500){
            for(int i = 1; i < nx; i++){
                for(int j = 1; j < ny; j++){
                    fprintf(file_ptr3, "%d %d %f\n", i, j, gsl_vector_get(T, l_function(i, j, nx)));
                    fprintf(file_ptr8, "%d %d %f\n", i, j, laplasjan(T, l_function(i, j, nx), delta, nx));
                }
            }
        }else if(it == 1000){
            for(int i = 1; i < nx; i++){
                for(int j = 1; j < ny; j++){
                    fprintf(file_ptr4, "%d %d %f\n", i, j, gsl_vector_get(T, l_function(i, j, nx)));
                    fprintf(file_ptr9, "%d %d %f\n", i, j, laplasjan(T, l_function(i, j, nx), delta, nx));
                }
            }
        }else if(it == 2000){
            for(int i = 1; i < nx; i++){
                for(int j = 1; j < ny; j++){
                    fprintf(file_ptr5, "%d %d %f\n", i, j, gsl_vector_get(T, l_function(i, j, nx)));
                    fprintf(file_ptr10, "%d %d %f\n", i, j, laplasjan(T, l_function(i, j, nx), delta, nx));
                }
            }
        }
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
    fclose(file_ptr3);
    fclose(file_ptr4);
    fclose(file_ptr5);
    fclose(file_ptr6);
    fclose(file_ptr7);
    fclose(file_ptr8);
    fclose(file_ptr9);
    fclose(file_ptr10);
}

void dyfuzja(int nx, int ny, double delta, double dt, double ta, double tb, double tc, double td, double kb, double kd, int sigma, int IT_MAX){
    int N = (nx+1)*(ny+1);
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);
    gsl_permutation *perm = gsl_permutation_alloc(N);

    wypelnij1(A, B, c, dt, delta, kb, kd, tb, td, nx, ny, N);
    wypelnij2(T, ta, tc, nx, ny);

    cn(A, B, c, d, T, perm, &sigma, nx, ny, delta, IT_MAX);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(d);
    gsl_vector_free(T);
    gsl_permutation_free( perm );
}

void podpunkt1(){
    dyfuzja(40, 40, 1, 1, 40, 0, 30, 0, 0.1, 0.6, 0, 2000);
}

int main() {

    podpunkt1();

    return 0;
}
