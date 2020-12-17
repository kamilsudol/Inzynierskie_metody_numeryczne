import numpy as np
import matplotlib.pyplot as plt

"""
    IMN LAB05 KAMIL SUDOL
    Drobne uwagi:
        Jako, że na wykonanie programu trzeba poczekać te 15 minut, to w ramach oszczedności czasu
        załączyłem Panu od razu wyniki.
        Pozdrawiam cieplutko :)
"""


def vij(vi_1j, v1_ij, vij_1, vi1_j):
    return 0.25*(vi_1j + v1_ij + vij_1 + vi1_j)


def stop(V, k, delta, sizex, sizey):
    suma = 0
    for i in range(sizex-k):
        for j in range(sizey-k):
            suma += ((k*delta**2)/2)*(((V[i+k][j]-V[i][j])/(2*k*delta) + (V[i+k][j+k]-V[i][j+k])/(2*k*delta))**2 + ((V[i][j+k]-V[i][j])/(2*k*delta) + (V[i+k][j+k]-V[i+k][j])/(2*k*delta))**2)
    return suma


def zageszczenie1(vij, vi_kj, vij_k, vi_kj_k):
    return 0.25*(vij + vi_kj + vij_k + vi_kj_k)


def zageszczenie2(vi_kj, vi_kj_k):
    return 0.5*(vi_kj + vi_kj_k)


def siatka():
    TOL = 10**(-8)
    nx = 128
    ny = 128
    delta = 0.2
    V = [[0 for x in range(ny + 1)] for y in range(nx + 1)]

    xmax = nx * delta
    ymax = ny * delta

    for i in range(ny+1):
        V[0][i] = np.sin(np.pi*delta*(i/ymax))
        V[nx][i] = np.sin(np.pi*delta*(i/ymax))

    for i in range(nx+1):
        V[i][0] = np.sin(2*np.pi*delta*(i/xmax))
        V[i][ny] = (-1)*np.sin(2*np.pi*delta*(i/xmax))
    wynik = []
    k = 16
    while k > 0:
        Sit = 1.0
        Suma = []

        while True:
            for i in range(k, nx-k+1, k):
                for j in range(k, ny-k+1, k):
                    V[i][j] = vij(V[i + k][j], V[i - k][j], V[i][j + k], V[i][j - k])

            Sit_1 = Sit
            Sit = stop(V, k, delta, nx, ny)
            Suma.append(Sit)

            if abs((Sit - Sit_1)/Sit_1) < TOL:
                break

        for i in range(0, nx - k+1, k):
            for j in range(0, ny - k+1, k):
                V[int(i + k / 2.)][int(j + k / 2.)] = zageszczenie1(V[i][j], V[i + k][j], V[i][j + k], V[i + k][j + k])
                V[i + k][int(j + k / 2.)] = zageszczenie2(V[i + k][j], V[i + k][j + k])
                V[int(i + k / 2.)][j + k] = zageszczenie2(V[i][j + k], V[i + k][j + k])
                V[int(i + k / 2.)][j] = zageszczenie2(V[i][j], V[i + k][j])
                V[i][int(j + k / 2.)] = zageszczenie2(V[i][j], V[i][j + k])

        ploter1(V, "k = " + str(k), "mapa_"+str(k)+".png")
        wynik.append(Suma)
        k = int(k/2)
    ploter3(wynik)
    
def ploter3(wynik):
    ploter2(wynik[0], wynik[1], wynik[2], wynik[3], wynik[4], 16, 8, 4, 2, 1, "S(it)")


def podpunkt1():
    siatka()


def ploter1(t,title, file):
    plt.imshow(t, cmap='jet')
    plt.colorbar(orientation='vertical')
    plt.title(title)
    plt.ylabel("x")
    plt.xlabel("y")
    plt.savefig(file)
    #plt.show()
    plt.clf()


def ploter2(s1, s2, s3, s4, s5, k1, k2, k3, k4, k5, title):
    plt.plot(range(len(s1)), s1, 'r-', label="k="+str(k1))
    plt.plot(range(len(s1), len(s1) + len(s2)), s2, 'b-', label="k="+str(k2))
    plt.plot(range(len(s1) + len(s2), len(s1) + len(s2)+len(s3)), s3, 'k-', label="k="+str(k3))
    plt.plot(range(len(s1) + len(s2)+len(s3), len(s1) + len(s2)+len(s3)+len(s4)), s4, 'g-', label="k="+str(k4))
    plt.plot(range(len(s1) + len(s2)+len(s3)+len(s4), len(s1) + len(s2)+len(s3)+len(s4)+len(s5)), s5, 'y-', label="k="+str(k5))
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel("S")
    plt.xlabel("Iteracja")
    plt.savefig("s(it).png")
    #plt.show()
    plt.clf()


def main(args):
    podpunkt1()


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
