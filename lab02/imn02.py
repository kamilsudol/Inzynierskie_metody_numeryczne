import numpy as np
import matplotlib.pyplot as plt

"""
    IMN LAB02 KAMIL SUDOL
"""

def function(beta, N, gamma, u):
    return (beta * N - gamma) * u - beta * (u ** 2)


def derivative_function(beta, N, gamma, u):
    return beta * N - gamma - 2 * beta * u


def function2(U1, U2, un, a1, a2, beta, N, gamma, dt):
    return U1 - un - dt*(a1 * function(beta, N, gamma, U1)+a2*function(beta, N, gamma, U2))


def matrix_function1(U, a, beta, N, gamma, dt):
    return (-1)*dt*a*derivative_function(beta, N, gamma, U)


def matrix_function2(U, a, beta, N, gamma, dt):
    return 1-dt*a*derivative_function(beta, N, gamma, U)


def dU_function(U1, U2, a11, a12, a21, a22, beta, N, gamma, dt, un):
    return (function2(U2, U1, un, a22, a21, beta, N, gamma, dt) * matrix_function1(U2, a12, beta, N, gamma, dt) - function2(U1, U2, un, a11, a12, beta, N, gamma, dt) * matrix_function2(U2, a22, beta, N, gamma, dt)) / (matrix_function2(U1, a11, beta, N, gamma, dt) * matrix_function2(U2, a22, beta, N, gamma, dt) - matrix_function1(U2, a12, beta, N, gamma, dt) * matrix_function1(U1, a21, beta, N, gamma, dt))


def ploter(x1, x2, title):
    plt.title(title)
    plt.plot(range(len(x1)), x1, 'k-', label='v(t)')
    plt.plot(range(len(x2)), x2, 'b-', label='u(t)')
    plt.legend(loc='upper right')
    plt.savefig(title+".png")
    #plt.show()
    plt.clf()


def podpunkt1():
    beta = 0.001
    N = 500
    gamma = 0.1
    tmax = 100
    dt = 0.1
    TOL = 10 ** (-6)

    step = tmax / dt
    u0 = 1.0
    un1 = u0

    v_t = []
    u_t = []

    for i in range(int(step)):
        u_t.append(un1)
        v_t.append(N - un1)
        iter = 0
        un = un1
        u_i = 0

        while abs(un1 - u_i) > TOL and iter <= 20:
            u_i = un1
            un1 = un + (dt / 2.0) * (function(beta, N, gamma, un) + function(beta, N, gamma, u_i))
            iter += 1

    ploter(v_t, u_t, "Metoda Pickarda")


def podpunkt2():
    beta = 0.001
    N = 500
    gamma = 0.1
    tmax = 100
    dt = 0.1
    TOL = 10 ** (-6)

    step = tmax / dt
    u0 = 1.0
    un1 = u0

    v_t = []
    u_t = []

    for i in range(int(step)):
        u_t.append(un1)
        v_t.append(N - un1)
        iter = 0
        un = un1
        u_i = 0

        while abs(un1 - u_i) > TOL and iter <= 20:
            u_i = un1
            un1 = u_i - (u_i - un - (dt / 2.0)*(function(beta, N, gamma, un)+function(beta, N, gamma, u_i))) / (1-(dt / 2.0)*derivative_function(beta, N, gamma, u_i))
            iter += 1

    ploter(v_t, u_t, "Iteracja Newtona")


def podpunkt3():
    beta = 0.001
    N = 500
    gamma = 0.1
    tmax = 100
    dt = 0.1
    TOL = 10 ** (-6)
    step = tmax / dt
    u0 = 1.0
    un1 = u0

    a11 = 0.25
    a12 = 0.25 - np.sqrt(3)/6
    a21 = 0.25 + np.sqrt(3)/6
    a22 = 0.25

    b1 = 0.5
    b2 = 0.5

    v_t = []
    u_t = []

    for i in range(int(step)):
        u_t.append(un1)
        v_t.append(N - un1)
        iter = 0
        U_i1 = 0
        U_i2 = 0
        un = un1
        U1 = un
        U2 = un

        while (abs(U1 - U_i1) > TOL or abs(U2 - U_i2) > TOL) and iter <= 20:
            U_i1 = U1
            U_i2 = U2
            dtU1 = dU_function(U1, U2, a11, a12, a21, a22, beta, N, gamma, dt, un)
            dtU2 = dU_function(U2, U1, a22, a21, a12, a22, beta, N, gamma, dt, un)
            U1 = U_i1+dtU1
            U2 = U_i2+dtU2
            iter += 1

        un1 = un + dt * (b1 * function(beta, N, gamma, U1) + b2 * function(beta, N, gamma,U2))

    ploter(v_t, u_t, "Niejawna RK2")


def main(args):
    podpunkt1()
    podpunkt2()
    podpunkt3()


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))