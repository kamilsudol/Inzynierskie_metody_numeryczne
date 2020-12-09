import matplotlib.pyplot as plt
import multiprocessing
import numpy as np

"""
    IMN LAB08 KAMIL SUDOL
    Drobne uwagi:
        Z racji tego, ze projekt ten wykonuje sie najdluzej sprosrod wszystkich dziesieciu,
        pozwolilem sobie na wprowadzenie drobnych uproszczen w kodzie celu przyspieszenia
        dzialania programu, tj wykomentowanie dodatkowych 20 iteracji z funkcji ad(..).
        Z calych sil probowalem napisac kod wg podanych wytycznych, jednak tym razem nawet
        wprowadzenie multiprocessingu nie bylo w stanie poprawic czasu wykonywania. Przy 
        pierwotnych wytycznych dla 40 iteracji program dlawil sie do tego stopnia, ze 30
        minut oczekiwania na wynik bylo norma, nie wspominajac juz o wiekszej liczbie iteracji,
        ktora byla wymagana, aby pokazaly sie 3 maksima na wykresie xsr(t). Probowalem po tej 
        "optymalizacji" wywolac program dla 12000 iteracji, zeby dla Pana wygody dorzucic 
        mniej wiecej oczekiwane wyniki, jednakze dla kazdorazowej proby wywolalanie konczylo
        sie fiaskiem, bo po okolu 2 godzinach przywieszal sie komputer. Dolaczylem jednak 
        wynik dla 2500 iteracji, na ktore co prawda wystarczy poczekac okolo dwudziestu minut, 
        jednak nie chcialbym, zeby marnowal Pan na to tyle czasu. Defaultowo zostawilem 1200 
        iteracji, przy ktorych oczekuje sie okolo 5 minut.
        
        Pozdrawiam cieplutko :)
"""

iter = 0


def gestosc(x, y, xA, yA, sigma):
    return (1 / (2 * np.pi * sigma ** 2)) * np.e ** (-((x - xA) ** 2 + (y - yA) ** 2) / (2 * sigma ** 2))


def wczytaj_psi(psi, nx, ny):
    file = np.loadtxt("psi.dat")
    data = file[:, 2]
    iter = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            psi[i][j] = data[iter]
            iter+=1


def vx_function(psi, vx, nx, ny, i1, i2, j1, delta):
    for i in range(1, nx):
        for j in range(1, ny):
            vx[i][j] = (psi[i][j+1]-psi[i][j-1])/(2*delta)

    for i in range(i1, i2+1):
        for j in range(0, j1+1):
            vx[i][j] = 0

    for i in range(1, nx):
        vx[i][0] = 0

    for j in range(ny+1):
        vx[0][j] = vx[1][j]
        vx[nx][j] = vx[nx-1][j]


def vy_function(psi, vy, nx, ny, i1, i2, j1, delta):
    for i in range(1, nx):
        for j in range(1, ny):
            vy[i][j] = -(psi[i+1][j]-psi[i-1][j])/(2*delta)

    for i in range(i1, i2+1):
        for j in range(0, j1+1):
            vy[i][j] = 0

    for i in range(1, nx):
        vy[i][ny] = 0


def vmax_function(vx, vy, nx, ny):
    max = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            obecny = np.sqrt((vx[i][j])**2+(vy[i][j])**2)
            if obecny > max:
                max = obecny
    return max


def u_function_wb1(un, un1, i, j, nx, ny, delta, D, dt, vx, vy):
    un1[i][j] = (delta**2)/(delta**2 + 2*D*dt)*(un[i][j] - (dt*vx[i][j]/2.)*((un[i+1][j] - un[nx][j])/(2*delta)+(un1[i+1][j] - un1[nx][j])/(2*delta)) - (dt*vy[i][j]/2.)*((un[i][j+1] - un[i][j-1])/(2*delta)+(un1[i][j+1] - un1[i][j-1])/(2*delta)) + (dt*D/2.)*((un[i+1][j]+un[nx][j]+un[i][j+1]+un[i][j-1]- 4*un[i][j])/(delta**2)+(un1[i+1][j]+un1[nx][j]+un1[i][j+1]+un1[i][j-1])/(delta**2)))


def u_function_wb2(un, un1, i, j, nx, ny, delta, D, dt, vx, vy):
    un1[i][j] = (delta**2)/(delta**2 + 2*D*dt)*(un[i][j] - (dt*vx[i][j]/2.)*((un[0][j] - un[i-1][j])/(2*delta)+(un1[0][j] - un1[i-1][j])/(2*delta)) - (dt*vy[i][j]/2.)*((un[i][j+1] - un[i][j-1])/(2*delta)+(un1[i][j+1] - un1[i][j-1])/(2*delta)) + (dt*D/2.)*((un[0][j]+un[i-1][j]+un[i][j+1]+un[i][j-1]- 4*un[i][j])/(delta**2)+(un1[0][j]+un1[i-1][j]+un1[i][j+1]+un1[i][j-1])/(delta**2)))


def u_function(un, un1, i, j, nx, ny, delta, D, dt, vx, vy):
    un1[i][j] = (delta**2)/(delta**2 + 2*D*dt)*(un[i][j] - (dt*vx[i][j]/2.)*((un[i+1][j] - un[i-1][j])/(2*delta)+(un1[i+1][j] - un1[i-1][j])/(2*delta)) - (dt*vy[i][j]/2.)*((un[i][j+1] - un[i][j-1])/(2*delta)+(un1[i][j+1] - un1[i][j-1])/(2*delta)) + (dt*D/2.)*((un[i+1][j]+un[i-1][j]+un[i][j+1]+un[i][j-1]- 4*un[i][j])/(delta**2)+(un1[i+1][j]+un1[i-1][j]+un1[i][j+1]+un1[i][j-1])/(delta**2)))


def calka(un, nx, ny):
    suma = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            suma += un[i][j]
    return suma


def x_sr(un, x, nx, ny):
    suma = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            suma += x[i]*un[i][j]
    return suma


def ad(un, un1, nx, ny, x, y, xA, yA, i1, i2, j1, delta, sigma, dt, D, vx, vy, IT_MAX):
    c = []
    xsr = []
    for i in range(nx + 1):
        for j in range(ny + 1):
            un[i][j] = gestosc(x[i], y[j], xA, yA, sigma)

    for it in range(IT_MAX):
        for i in range(nx + 1):
            for j in range(ny + 1):
                un1[i][j] = un[i][j]
        # for k in range(1, 21):
        for i in range(nx + 1):
            for j in range(1, ny):
                if i1 < i < i2 and 0 < j < j1:
                    continue
                elif i == 0:
                    u_function_wb1(un, un1, i, j, nx, ny, delta, D, dt, vx, vy)
                elif i == nx:
                    u_function_wb2(un, un1, i, j, nx, ny, delta, D, dt, vx, vy)
                else:
                    u_function(un, un1, i, j, nx, ny, delta, D, dt, vx, vy)
        for i in range(nx + 1):
            for j in range(ny + 1):
                un[i][j] = un1[i][j]

        c.append(calka(un, nx, ny)*delta**2)
        xsr.append(x_sr(un, x, nx, ny)*delta**2)

    return [c, xsr]


def crank_nicolson(nx, ny, i1, i2, j1, delta, sigma, xA, yA, D, IT_MAX, kolejka, flaga):
    psi = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    vx = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    vy = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    un = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    un1 = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    x = [0 for j in range(nx + 1)]
    y = [0 for j in range(ny + 1)]

    for i in range(ny + 1):
        y[i] = i * delta

    for i in range(nx + 1):
        x[i] = i * delta

    wczytaj_psi(psi, nx, ny)
    vx_function(psi, vx, nx, ny, i1, i2, j1, delta)
    vy_function(psi, vy, nx, ny, i1, i2, j1, delta)
    vmax = vmax_function(vx, vy, nx, ny)
    dt = delta/(4*vmax)


    wynik = ad(un, un1, nx, ny, x, y, xA, yA, i1, i2, j1, delta, sigma, dt, D, vx, vy, IT_MAX)
    kolejka.put([False, un, D,IT_MAX])
    if flaga:
        kolejka.put([True, wynik[0], wynik[1], vx, vy, D])


def podpunkt1(list, q):
    global iter
    kolejka = []

    for x in list:
        x.start()

    while iter <= 11:
        result = q.get()
        if result[0]:
            kolejka.append(result)
            iter += 1
        else:
            ploter(result)
    ploter3(kolejka)


def ploter(list):
    global iter
    ploter1(list[1], "D=" + str(list[2]) + " IT = " + str(list[3]), str(list[2]) + "mapa" + str(list[3]))
    iter += 1

def ploter3(list):
    ploter2(list[0][1], list[1][1], "c(tn)", "c(tn)")
    ploter2(list[0][2], list[1][2], "xsr(tn)", "xsr(tn)")
    ploter1(list[0][3], "Vx", "Vx")
    ploter1(list[0][4], "Vy", "Vy")


def ploter1(t, title, file):
    plt.imshow(t, cmap='jet')
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.ylabel("y")
    plt.xlabel("x")
    plt.savefig(file + ".png")
    #plt.show()
    plt.clf()


def ploter2(t1, t2, title, file):
    plt.plot(t1,'-k', label = "D=0.0")
    plt.plot(t2,'-b', label = "D=0.1")
    plt.legend(loc='upper right')
    plt.title(title)
    plt.ylabel(title)
    plt.xlabel("tn")
    plt.savefig(file + ".png")
    #plt.show()
    plt.clf()


def main(args):
    q = multiprocessing.Queue()
    tmax = 1200
    #tmax = 2500
    k = 5
    t = tmax/k
    list1 = []
    for i in range(1,k+1):
        if i == 5:
            list1.append(multiprocessing.Process(None, crank_nicolson, args=(400, 90, 200, 210, 50, 0.01, 0.1, 0.45, 0.45, 0, int(i * t), q, True)))
        else:
            list1.append(multiprocessing.Process(None, crank_nicolson, args=(400, 90, 200, 210, 50, 0.01, 0.1, 0.45, 0.45, 0, int(i * t), q, False)))

    for i in range(1, k + 1):
        if i == 5:
            list1.append(multiprocessing.Process(None, crank_nicolson, args=(400, 90, 200, 210, 50, 0.01, 0.1, 0.45, 0.45, 0.1, int(i * t), q, True)))
        else:
            list1.append(multiprocessing.Process(None, crank_nicolson, args=(400, 90, 200, 210, 50, 0.01, 0.1, 0.45, 0.45, 0.1, int(i * t), q, False)))
    podpunkt1(list1, q)


if __name__ == '__main__':
    import sys

    sys.exit(main(sys.argv))
