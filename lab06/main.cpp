#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "mgmres.h"

/*
    IMN LAB06 KAMIL SUDOL
*/

double ro1(double x, double y, double sigma, double xmax, double ymax){
    return exp(-pow((x-0.25*xmax)/sigma, 2) - pow((y-0.5*ymax)/sigma, 2));
}

double ro2(double x, double y, double sigma, double xmax, double ymax){
    return -exp(-pow((x-0.75*xmax)/sigma, 2) - pow((y-0.5*ymax)/sigma, 2));
}

int j_function(int l, int nx){
    return floor(l/(nx+1));
}

int i_function(int l, int nx){
    return l-j_function(l, nx)*(nx+1);
}

double e_function(double e1, double e2, int l, double nx){
    if(i_function(l, nx) <= (nx/2.)){
        return e1;
    }else{
        return e2;
    }
}

int wypelnienie_macierzy(bool flag, bool file_flag, int N, double V1, double V2, double V3, double V4, int nx, int ny, double sigma, double delta, double xmax, double ymax, double e1, double e2, int *ia, int *ja, double *a, double *b, double *V){
    int k = -1;
    for(int l = 0; l < N; l++){
        int brzeg = 0;
        double vb = 0.0;

        if(i_function(l, nx) == 0){
            brzeg = 1;
            vb = V1;
        }

        if(j_function(l, nx) == ny){
            brzeg = 1;
            vb = V2;
        }

        if(i_function(l, nx) == nx){
            brzeg = 1;
            vb = V3;
        }

        if(j_function(l, nx) == 0){
            brzeg = 1;
            vb = V4;
        }

        if(flag){
            b[l] = 0.0;
        }else{
            b[l] = -(ro1(delta*i_function(l, nx), delta*j_function(l, nx), sigma, xmax, ymax)+ro2(delta*i_function(l, nx), delta*j_function(l, nx), sigma, xmax, ymax));
        }

        if(brzeg == 1){
            b[l] = vb;
        }

        ia[l] = -1;

        if(l - nx - 1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0) ia[l] = k;
            a[k] = e_function(e1, e2, l, nx)/pow(delta,2);
            ja[k] = l-nx-1;
        }

        if(l - 1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0) 
                ia[l] = k;
            a[k] = e_function(e1, e2, l, nx)/pow(delta,2);
            ja[k] = l-1;
        }

        k++;
        if(ia[l] < 0) 
            ia[l] = k;
        if(brzeg == 0){
            a[k] = -(2*e_function(e1, e2, l, nx)+e_function(e1, e2, l+1, nx)+e_function(e1, e2, l+nx+1, nx))/pow(delta,2);
        }else {
            a[k] = 1;
        }
        ja[k] = l;

        if(l < N && brzeg == 0){
            k++;
            a[k] = e_function(e1, e2, l+1, nx)/pow(delta,2);
            ja[k] = l+1;
        }

        if(l < N - nx - 1 && brzeg == 0){
            k++;
            a[k] = e_function(e1, e2, l+nx+1, nx)/pow(delta,2);
            ja[k] = l+nx+1;
        }
    }
    int nz_num = k + 1;
    ia[N] = nz_num;

    if(file_flag){
        FILE *file_ptr1 = fopen("wypelnienie_macierzy.out", "w");
        FILE *file_ptr2 = fopen("wypelnienie_wektora.out", "w");

        for(int l = 0; l < 5*N; l++){
            fprintf(file_ptr1,"%d %d %d %f\n", l, i_function(l, nx), j_function(l, nx), a[l]);
        }
        
        for(int l = 0; l < N; l++){
            fprintf(file_ptr2,"%d %d %d %f\n", l, i_function(l, nx), j_function(l, nx), b[l]);
        }

        fclose(file_ptr1);
        fclose(file_ptr2);
    }

    return nz_num;
}

void algeb(int nx, int ny, double V1, double V2, double V3, double V4, double e1, double e2, double delta, bool flag, bool file_flag, const char * plik){
    int N = (nx+1)*(ny+1);
    double xmax = delta*nx;
    double ymax = delta*ny;
    double sigma = xmax/10.;

    int ja[5*N], ia[N+1];
    double a[5*N], b[N], V[N];

    int itr_max = 500, mr = 500;
    double TOL_ABS = pow(10, -8), TOL_REL = pow(10,-8);

    int nz_num = wypelnienie_macierzy(flag, file_flag, N, V1, V2, V3, V4, nx, ny, sigma, delta, xmax, ymax, e1, e2, ia, ja, a, b, V);

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, TOL_ABS, TOL_REL);

    FILE * file_ptr = fopen(plik, "w");

    for(int i = 0; i < N; i++){
        fprintf(file_ptr, "%f %f %f\n", delta*i_function(i, nx), delta*j_function(i, nx), V[i]);
    }

    fclose(file_ptr);
}

int main() {
    double delta = 0.1;
    algeb(4, 4, 10, -10, 10, -10, 1, 1, delta, true, true, "mapa1.out");

    algeb(50, 50, 10, -10, 10, -10, 1, 1, delta, true, false, "mapa2.out");
    algeb(100, 100, 10, -10, 10, -10, 1, 1, delta, true, false, "mapa3.out");
    algeb(200, 200, 10, -10, 10, -10, 1, 1, delta, true, false, "mapa4.out");

    algeb(100, 100, 0, 0, 0, 0, 1, 1, delta, false, false, "mapa5.out");
    algeb(100, 100, 0, 0, 0, 0, 1, 2, delta, false, false, "mapa6.out");
    algeb(100, 100, 0, 0, 0, 0, 1, 10, delta, false, false, "mapa7.out");

    return 0;
}
