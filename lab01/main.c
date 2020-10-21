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

double potential(double omega, double dt){
    return 10*sin(omega*dt);
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

void podpunkt4(double omega, const char * plik1, const char * plik2){
    FILE *file_ptr1 = fopen(plik1, "w"), *file_ptr2 = fopen(plik2, "w");

    double kq1, kq2, kq3, kq4, ki1, ki2, ki3, ki4;

    double R = 100, L = 0.1, C = 0.001, dt = 0.0004;
    double w0 = 1/sqrt(L*C), T0 = (2*M_PI)/w0, t = 0;
    double Q0 = 0, Qn, I0 = 0, In;
    int i = 0;

    while(t <= 4*T0){
        kq1 = I0;
        ki1 = potential(omega*w0, t)/L -Q0/(L*C) - (R*I0)/L;
        kq2 = k_function(1, dt/2.0, I0, ki1);
        ki2 = potential(omega*w0, t + dt/2.0)/L - k_function(1/(L*C), dt/2.0, Q0, kq1)-k_function(R/L, dt/2.0, I0, ki1);
        kq3 = k_function(1, dt/2.0, I0, ki2);
        ki3 = potential(omega*w0, t + dt/2.0)/L - k_function(1/(L*C), dt/2.0, Q0, kq2)-k_function(R/L, dt/2.0, I0, ki2);
        kq4 = k_function(1, dt, I0, ki3);
        ki4 = potential(omega*w0, t + dt)/L - k_function(1/(L*C), dt, Q0, kq3)-k_function(R/L, dt, I0, ki3);

        Qn = Q0 + (dt/6.0)*(kq1+2*kq2+2*kq3+kq4);
        In = I0 + (dt/6.0)*(ki1+2*ki2+2*ki3+ki4);

        fprintf(file_ptr1, "%d %f\n", i, Q0);
        fprintf(file_ptr2, "%d %f\n", i, I0);

        i++;
        Q0 = Qn;
        I0 = In;
        t += dt;
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
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

    podpunkt4(0.5, "Q1.txt", "I1.txt");
    podpunkt4(0.8, "Q2.txt", "I2.txt");
    podpunkt4(1.0, "Q3.txt", "I3.txt");
    podpunkt4(1.2, "Q4.txt", "I4.txt");
    
    return 0;
}