import matplotlib.pyplot as plt
import numpy as np

"""
    IMN LAB10 KAMIL SUDOL
"""


def x_i(i, delta):
    return delta*i


def t_i(i, dt):
    return dt*i


def kronecker(x, xf):
    if x == xf:
        return 1
    else:
        return 0


def a_f(i, xf, nt, dt, t):
    return np.cos(50*t/(nt*dt))*kronecker(i,xf)


def w_p(u, nx, xa, sigma, delta):
    for i in range(1, nx - 1):
        u[i] = np.exp((-(x_i(i, delta)-xa)**2)/(2*sigma**2))


def licz_e(delta, u,v, nx):
    s1 = (delta/4)*(((u[1]-u[0])/delta)**2 + ((u[nx]-u[nx-1])/delta)**2)
    s2 = 0
    for i in range(1,nx):
        s2+=v[i]**2 + ((u[i+1]-u[i-1])/(2*delta))**2
    return s1 + (delta/2)*s2


def wyznacz_a(a, u, u0, alfa, beta, delta, dt, nx, nt, t, xf):
    for i in range(1,nx):
        a[i]=(u[i+1]-2*u[i]+u[i-1])/(delta**2) - beta*(u[i]-u0[i])/dt + alfa*a_f(x_i(i, delta), xf, nt, dt, t)


def rownanie_falowe_alg(u, u0, v, vp, a, nx, nt, delta, dt, alfa, beta, xf):
    E = []
    U_wynik = []

    u0 = u[:]
    wyznacz_a(a, u, u0, alfa, beta, delta, dt, nx, nt, 0, xf)
    for it in range(1, nt+1):
        for i in range(1, nx):
            vp[i] = v[i] + dt * a[i] / 2
            u0[i] = u[i]
            u[i] = u[i] + dt * vp[i]

        wyznacz_a(a, u, u0, alfa, beta, delta, dt, nx, nt, t_i(it, dt), xf)
        for i in range(1, nx):
            v[i] = vp[i] + dt * a[i] / 2
        U_wynik.append(u[:])
        E.append(licz_e(delta, u, v, nx))

    return [E, U_wynik]


def rownanie_falowe(nx, nt, delta, dt, beta, alfa, xa, sigma, xf, flag):
    u = [0 for x in range(nx + 1)]
    u0 = [0 for x in range(nx + 1)]
    v = [0 for x in range(nx + 1)]
    vp = [0 for x in range(nx + 1)]
    a = [0 for x in range(nx + 1)]

    if flag:
        w_p(u, nx, xa, sigma, delta)

    wynik = rownanie_falowe_alg(u, u0, v, vp, a, nx, nt, delta, dt, alfa, beta, xf)
    ploter1(wynik[1], "beta=" + str(beta) + ", alfa=" + str(alfa), "struna("+str(alfa)+","+str(beta)+")")
    return wynik[0]


def podpunkt1():
    wynik = [rownanie_falowe(150, 1000, 0.1, 0.05, 0, 0, 7.5, 0.5, 2.5, True),
    rownanie_falowe(150, 1000, 0.1, 0.05, 0.1, 0, 7.5, 0.5, 2.5, True),
    rownanie_falowe(150, 1000, 0.1, 0.05, 1, 0, 7.5, 0.5, 2.5, True),
    rownanie_falowe(150, 1000, 0.1, 0.05, 1, 1, 7.5, 0.5, 2.5, False)]

    ploter2(wynik[0], wynik[1], wynik[2], "E(t)", "E(t)")
    ploter3(wynik[3], "E(t)", "E(t)pt4")


def ploter1(t, title, file):
    plt.imshow(t, cmap='gnuplot')
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.ylabel("t")
    plt.xlabel("x")
    plt.savefig(file + ".png")
    plt.show()
    plt.clf()


def ploter2(t1, t2, t3, title, file):
    plt.plot(t1,'-r', label = "b=0.0, a=0.0")
    plt.plot(t2,'-b', label = "b=0.1, a=0.0")
    plt.plot(t3,'-g', label = "b=1.0, a=0.0")
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel("x")
    plt.xlabel("t")
    plt.savefig(file + ".png")
    plt.show()
    plt.clf()

def ploter3(t1, title, file):
    plt.plot(t1,'-r', label = "b=0.0, a=0.0")
    plt.legend(loc='lower right')
    plt.title(title)
    plt.ylabel("x")
    plt.xlabel("t")
    plt.savefig(file + ".png")
    plt.show()
    plt.clf()


def main(args):
    podpunkt1()


if __name__ == '__main__':
    import sys

    sys.exit(main(sys.argv))
