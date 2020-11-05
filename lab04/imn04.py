import numpy as np
import matplotlib.pyplot as plt

"""
    IMN LAB04 KAMIL SUDOL
"""


def ro1(x, y, xmax, ymax, gammax, gammay):
    return np.e**(-((x-0.35*xmax)**2)/gammax**2-((y-0.5*ymax)**2)/gammay**2)


def ro2(x, y, xmax, ymax, gammax, gammay):
    return (-1)*np.e**(-((x - 0.65 * xmax) ** 2) / gammax ** 2 - ((y - 0.5 * ymax) ** 2) / gammay ** 2)


def ro(x, y, xmax, ymax, gammax, gammay):
    return ro1(x, y, xmax, ymax, gammax, gammay) + ro2(x, y, xmax, ymax, gammax, gammay)


def Gvij(vi_1j, v1_ij, vij_1, vi1_j, delta, epsilon, ro):
    return (1/4)*(vi_1j + v1_ij + vij_1 + vi1_j + (delta**2/epsilon)*ro)


def Lvij(omega, vij, vi_1j, v1_ij, vij_1, vi1_j, delta, epsilon, ro):
    return (1.0-omega)*vij+(omega/4)*(vi_1j + v1_ij + vij_1 + vi1_j + (delta**2/epsilon)*ro)


def stop(V, Ro, delta, sizex, sizey):
    suma = 0
    for i in range(sizex):
        for j in range(sizey):
            suma += (delta**2)*(0.5*((V[i+1][j]-V[i][j])/delta)**2 + 0.5*((V[i][j+1]-V[i][j])/delta)**2 - Ro[i][j]*V[i][j])
    return suma


def sigma(vi1j, vij, vi_1j, vij1, vij_1, delta, ro, epsilon):
    return (vi1j- 2*vij+ vi_1j)/delta**2 + (vij1 - 2*vij + vij_1)/delta**2 + ro/epsilon


def glob(omega):
    TOL = 10**(-8)
    nx = 150
    ny = 100
    delta = 0.1
    epsilon = 1.0
    V1 = 10.0
    Vstare = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    Vnowe = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    Sigma = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    Ro = [[0 for x in range(ny + 1)] for y in range(nx + 1)]

    for i in range(nx+1):
        Vstare[i][0]=V1
        Vnowe[i][0]=V1

    xmax = nx * delta
    ymax = ny * delta

    gammax = 0.1 * xmax
    gammay = 0.1 * ymax

    Sit = 1.0
    Suma = []

    for i in range(1, nx):
        for j in range(1, ny):
            Ro[i][j] = ro(i * delta, j * delta, xmax, ymax, gammax, gammay)

    while True:

        for i in range(1, nx):
            for j in range(1, ny):
                Vnowe[i][j] = Gvij(Vstare[i + 1][j], Vstare[i - 1][j], Vstare[i][j + 1], Vstare[i][j - 1], delta, epsilon, Ro[i][j])

        for j in range(1, ny+1):
            Vnowe[0][j] = Vnowe[1][j]
            Vnowe[nx][j] = Vnowe[nx-1][j]

        for i in range(nx + 1):
            for j in range(ny + 1):
                Vstare[i][j] = (1.0-omega)*Vstare[i][j]+omega*Vnowe[i][j]

        Sit_1 = Sit
        Sit = stop(Vnowe, Ro, delta, nx, ny)
        Suma.append(Sit)

        if abs((Sit - Sit_1)/Sit_1) < TOL:
            break

    for i in range(1, nx):
        for j in range(1, ny):
            Sigma[i][j]=sigma(Vnowe[i+1][j],Vnowe[i][j],Vnowe[i-1][j],Vnowe[i][j+1],Vnowe[i][j-1],delta, Ro[i][j], epsilon)

    ploter(Vnowe, "Wykres relaksacji (w="+str(omega)+")", "relaksacja_"+str(omega))
    ploter(Sigma, "Wykres relaksacji (w="+str(omega)+") - blad", "sigma_"+str(omega))
    return Suma


def lokal(omega):
    TOL = 10**(-8)
    nx = 150
    ny = 100
    delta = 0.1
    epsilon = 1.0
    V1 = 10.0

    V = [[0 for x in range(ny+1)] for y in range(nx+1)]
    Ro = [[0 for x in range(ny+1)] for y in range(nx+1)]

    for i in range(nx+1):
        V[i][0]=V1

    xmax = nx * delta
    ymax = ny * delta

    gammax = 0.1 * xmax
    gammay = 0.1 * ymax

    Sit = 1.0
    Suma = []

    for i in range(1, nx):
        for j in range(1, ny):
            Ro[i][j] = ro(i * delta, j * delta, xmax, ymax, gammax, gammay)

    while True:

        for i in range(1, nx):
            for j in range(1, ny):
                V[i][j] = Lvij(omega, V[i][j], V[i + 1][j], V[i - 1][j], V[i][j + 1], V[i][j - 1], delta, epsilon, Ro[i][j])

        for j in range(1, ny+1):
            V[0][j] = V[1][j]
            V[nx][j] = V[nx-1][j]

        Sit_1 = Sit
        Sit = stop(V, Ro, delta, nx, ny)
        Suma.append(Sit)
        if abs((Sit - Sit_1)/Sit_1) < TOL:
            break

    return Suma


def podpunkt1():
    title = "Relaksacja globalna"
    wynik = [glob(0.6), glob(1.0)]
    ploter1(wynik[0], wynik[1], 0.6, 1.0, title)


def podpunkt2():
    title = "Relaksacja lokalna"
    wynik = [lokal(1.0), lokal(1.4), lokal(1.8), lokal(1.9)]
    ploter2(wynik[0], wynik[1], wynik[2], wynik[3], 1.0, 1.4, 1.8, 1.9, title)


def ploter(t,title, file):
    plt.imshow(t, cmap='jet')
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.ylabel("x")
    plt.xlabel("y")
    plt.savefig(file+".png")
    plt.show()


def ploter1(s1, s2, o1, o2, title):
    plt.xscale("log")
    plt.plot(range(len(s1)), s1, 'r-', label="w="+str(o1))
    plt.plot(range(len(s2)), s2, 'b-', label="w="+str(o2))
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel("S")
    plt.xlabel("Iteracja")
    plt.savefig("s(it)_glob.png")
    plt.show()
    plt.clf()


def ploter2(s1, s2, s3, s4, o1, o2, o3, o4, title):
    plt.xscale("log")
    plt.plot(range(len(s1)), s1, 'r-', label="w="+str(o1))
    plt.plot(range(len(s2)), s2, 'b-', label="w="+str(o2))
    plt.plot(range(len(s3)), s3, 'k-', label="w="+str(o3))
    plt.plot(range(len(s4)), s4, 'g-', label="w="+str(o4))
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel("S")
    plt.xlabel("Iteracja")
    plt.savefig("s(it)_lok.png")
    plt.show()
    plt.clf()


def main(args):
    podpunkt1()


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
