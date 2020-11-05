import matplotlib.pyplot as plt

"""
    IMN LAB03 KAMIL SUDOL
"""


def f_function(v):
    return v


def g_function(alfa, x, v):
    return alfa * (1 - x ** 2) * v - x


def schemat_trapez_1(x1, x, v1, v, dt, alfa):
    para = [x1-x-(dt/2.0)*(f_function(v)+f_function(v1)), v1-v-(dt/2.0)*(g_function(alfa, x, v)+g_function(alfa, x1, v1))]
    return para


def a11():
    return 1.0


def a12(dt):
    return -dt/2.0


def a21(dt, alfa, x, v):
    return (-dt/2.0)*((-2)*alfa*x*v - 1)


def a22(dt, alfa, x):
    return 1 - alfa*(dt/2.0)*(1-x**2)


def schemat_trapez_2(alfa, x, v, dt):
    sigma = 10**(-10)
    xkn = x
    vkn = v

    while True:
        k = schemat_trapez_1(xkn, x, vkn, v, dt, alfa)
        deltaX = (a22(dt, alfa, xkn)*(-k[0])-(-k[1])*a12(dt))/(a11()*a22(dt, alfa, xkn) - a12(dt)*a21(dt, alfa, xkn, vkn))
        deltaV = (a11()*(-k[1])-(-k[0])*a21(dt, alfa, xkn, vkn))/(a11()*a22(dt, alfa, xkn) - a12(dt)*a21(dt, alfa, xkn, vkn))

        xkn += deltaX
        vkn += deltaV

        if abs(deltaV) < sigma or abs(deltaX) < sigma:
            break

    para = [xkn, vkn]
    return para


def schemat_rk2_1(alfa, x, v):
    para = [f_function(v), g_function(alfa, x, v)]
    return para


def schemat_rk2_2(alfa, x, v, dt):
    k1 = schemat_rk2_1(alfa, x, v)
    para = [v + k1[1] * dt, alfa*(1-(x+dt*k1[0])**2)*(v+dt*k1[0])-(x+dt*k1[1])]
    return para


def schemat_rk2_3(alfa, x, v, dt):
    k1 = schemat_rk2_1(alfa, x, v)
    k2 = schemat_rk2_1(alfa, x, v)
    para = [x+(dt/2.0)*(k1[0]+k2[0]), v+(dt/2.0)*(k1[1]+k2[1])]
    return para


def kontrola_kroku_czasowego(funkcja, TOL):
    x0 = 0.01
    v0 = 0
    dt0 = 1
    S = 0.75
    p = 2
    tmax = 40
    alfa = 5
    t = 0
    dt = dt0
    xn = x0
    vn = v0

    czas = []
    delta_czas = []
    x_t = []
    v_t = []

    while True:
        k1 = eval(funkcja)(alfa, xn, vn, dt)
        k2 = eval(funkcja)(alfa, k1[0], k1[1], dt)
        k3 = eval(funkcja)(alfa, xn, vn, 2 * dt)

        Ex = (k2[0]-k3[0])/(2**p - 1)
        Ev = (k2[1]-k3[1])/(2**p - 1)

        if max(abs(Ex), abs(Ev)) < TOL:
            t += 2*dt
            xn = k2[0]
            vn = k2[1]
            czas.append(t)
            delta_czas.append(dt)
            v_t.append(vn)
            x_t.append(xn)
        else:
            dt = dt*((S*TOL)/(max(abs(Ex), abs(Ev))))**(1/(p+1))

        if t > tmax:
            break

    return [czas, delta_czas, x_t, v_t]


def podpunkt1():
    wynik1 = kontrola_kroku_czasowego("schemat_trapez_2", 10**(-2))
    wynik2 = kontrola_kroku_czasowego("schemat_trapez_2", 10**(-5))

    ploter(wynik1[0], wynik1[3], wynik2[0], wynik2[3], "v(t)", "tol = 10**(-2)", "tol = 10**(-5)", "trapez_v(t)")
    ploter(wynik1[0], wynik1[2], wynik2[0], wynik2[2], "x(t)", "tol = 10**(-2)", "tol = 10**(-5)", "trapez_x(t)")
    ploter(wynik1[0], wynik1[1], wynik2[0], wynik2[1], "dt(t)", "tol = 10**(-2)", "tol = 10**(-5)", "trapez_dt(t)")
    ploter(wynik1[2], wynik1[3], wynik2[2], wynik2[3], "v(x)", "tol = 10**(-2)", "tol = 10**(-5)", "trapez_v(x)")


def podpunkt2():
    wynik1 = kontrola_kroku_czasowego("schemat_rk2_3", 10**(-2))
    wynik2 = kontrola_kroku_czasowego("schemat_rk2_3", 10**(-5))

    ploter(wynik1[0], wynik1[3], wynik2[0], wynik2[3], "v(t)", "tol = 10**(-2)", "tol = 10**(-5)", "rk2_v(t)")
    ploter(wynik1[0], wynik1[2], wynik2[0], wynik2[2], "x(t)", "tol = 10**(-2)", "tol = 10**(-5)", "rk2_x(t)")
    ploter(wynik1[0], wynik1[1], wynik2[0], wynik2[1], "dt(t)", "tol = 10**(-2)", "tol = 10**(-5)", "rk2_dt(t)")
    ploter(wynik1[2], wynik1[3], wynik2[2], wynik2[3], "v(x)", "tol = 10**(-2)", "tol = 10**(-5)", "rk2_v(x)")


def ploter(x1, y1, x2, y2, title, label1, label2, file):
    plt.plot(x1, y1, 'r-', label=label1)
    plt.plot(x2, y2, 'b-', label=label2)
    plt.legend(loc='upper right')
    plt.title(title)
    plt.savefig(file + ".png")
    #plt.show()
    plt.clf()


def main(args):
    podpunkt1()
    podpunkt2()


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
