#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double function1(double y){
    double lambda = -1;
    return lambda*y;
}

double function2(double dt){
    double lambda = -1;
    return pow(M_E, lambda*dt);
}

double k_function(double lambda, double dt, double y, double k){
    return lambda*(y + dt*k);
}

void podpunkt1(double dt, const char * plik1, const char * plik2, const char * plik3){
    FILE *file_ptr1 = fopen(plik1, "w"), *file_ptr2 = fopen(plik2, "w"), *file_ptr3 = fopen(plik3, "w");

    int step = 5/dt;
    double y0 = 1, yn, t = 0;

    for(int i=0; i<step; i++){
        yn = y0 + dt*function1(y0);
        fprintf(file_ptr1,"%d %f\n", i, y0);
        fprintf(file_ptr2,"%d %f\n", i, function2(t));
        fprintf(file_ptr3,"%d %f\n", i, y0 - function2(t));
        t += dt;
        y0 = yn;
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
    fclose(file_ptr3);
}

void podpunkt2(double dt, const char * plik1, const char * plik2, const char * plik3){
    FILE *file_ptr1 = fopen(plik1, "w"), *file_ptr2 = fopen(plik2, "w"), *file_ptr3 = fopen(plik3, "w");

    int step = 5/dt;
    double y0 = 1, yn, k1, k2;
    double lambda = -1, t = 0;

    for(int i=0; i<step; i++){
        k1 = function1(y0);
        k2 = k_function(lambda, dt, y0, k1);
        yn = y0 + (dt/2.0)*(k1+k2);
        fprintf(file_ptr1,"%d %f\n", i, y0);
        fprintf(file_ptr2,"%d %f\n", i, function2(t));
        fprintf(file_ptr3,"%d %f\n", i, y0 - function2(t));
        t += dt;
        y0 = yn;
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
    fclose(file_ptr3);
}

void podpunkt3(double dt, const char * plik1, const char * plik2, const char * plik3){
    FILE *file_ptr1 = fopen(plik1, "w"), *file_ptr2 = fopen(plik2, "w"), *file_ptr3 = fopen(plik3, "w");

    int step = 5/dt;
    double y0 = 1, yn, k1, k2, k3, k4;
    double lambda = -1, t = 0;

    for(int i=0; i<step; i++){
        k1 = function1(y0);
        k2 = k_function(lambda, dt/2.0, y0, k1);
        k3 = k_function(lambda, dt/2.0, y0, k2);
        k4 = k_function(lambda, dt, y0, k3);
        yn = y0 + (dt/6.0)*(k1+2*k2+2*k3+k4);
        fprintf(file_ptr1,"%d %f\n", i, y0);
        fprintf(file_ptr2,"%d %f\n", i, function2(t));
        fprintf(file_ptr3,"%d %f\n", i, y0 - function2(t));
        y0 = yn;
        t += dt;
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
    fclose(file_ptr3);
}

void podpunkt4(){
    
}

int main(){

    podpunkt1(0.01, "euler1.txt", "eulera1.txt", "eulerb1.txt");
    podpunkt1(0.1, "euler2.txt", "eulera2.txt", "eulerb2.txt");
    podpunkt1(1.0, "euler3.txt", "eulera3.txt", "eulerb3.txt");

    podpunkt2(0.01, "RK21.txt", "RK2a1.txt", "RK2b1.txt");
    podpunkt2(0.1, "RK22.txt", "RK2a2.txt", "RK2b2.txt");
    podpunkt2(1.0, "RK23.txt", "RK2a3.txt", "RK2b3.txt");
    
    podpunkt3(0.01, "RK41.txt", "RK4a1.txt", "RK4b1.txt");
    podpunkt3(0.1, "RK42.txt", "RK4a2.txt", "RK4b2.txt");
    podpunkt3(1.0, "RK43.txt", "RK4a3.txt", "RK4b3.txt");
    
    return 0;
}