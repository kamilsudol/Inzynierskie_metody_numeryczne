import matplotlib.pyplot as plt
import multiprocessing

"""
    IMN LAB07 KAMIL SUDOL
    Drobne uwagi:
        Ze względu na to, że przy podstawowym zestawie parametrów zaraz po 2000 iteracji zaczęły się tworzyć 'NaN-y' oraz
        przez to, że program zwyczajnie wykonywał się nieakceptowalnie długo, postanowiłem nieco zmodyfikować program,
        aby zwrócił chociaż jakieś sensowne wyniki oraz przy okazji zaoszczędzić trochę Pańskiego czasu. W tym celu 
        stworzyłem trzy wersje wywołania w mainie:
        
        - podpunkt1(lista1) - przy obliczaniu relaksacji usunąłem część równania mnożonego przez omegę, która powodowała
        'NaN-y' oraz przy tej wersji program wykonuje się najszybciej, bo około 5 min.
        
        - podpunkt1(lista2) - w tej wersji część równania mnożonego przez omegę już występuje, jednakże w środku wykonywane
        są zaokrąglenia liczb w celu uniknięcia 'NaN-ów', przez co program wykonuje się najdłużej spośród wszystkich opcji, 
        bo około pół godziny.
        
        - podpunkt1(lista3) - w tej wersji nie modyfikowałem algorytmu z polecenia, jednakże odpowiednio zwiększyłem liczbę
        iteracji oraz próg zmiany omegi z zera na jeden w celu uniknięcia 'NaN-ów'. Wersja ta wykonuje się około 10 min.
        
        Domyślnie ustawiona jest wersja pierwsza, jednakże jeżeli miałby Pan czas oraz chęci, to w mainie na dole 
        może Pan odkomentować interesująca Pana wersję do wykonania.
        
        W przypadku uruchamiania programu na taurusie prosze o usunięcie całej sekcji tego komentarza, ponieważ pojawia
        się błąd spowodowany polskimi znakami.
        
        Pozdrawiam cieplutko :)
"""


iter = 0


def Qwy(Qwe, y, j1, ny):
    return Qwe*((y[ny] - y[j1])/y[ny])**3


def wb_psi(psi, Qwe, mi, i1, j1, y, nx, ny):
    for j in range(j1, ny+1):
        psi[0][j] = (Qwe/(2*mi))*((y[j]**3)/3-(y[j]**2)*(y[j1]+y[ny])/2 +y[j]*y[j1]*y[ny])

    for j in range(ny+1):
        psi[nx][j] = (Qwy(Qwe,y,j1, ny)/(2*mi))*((y[j]**3)/3-(y[j]**2)*(y[ny])/2) + Qwe*(y[j1]**2)*(3*y[ny] - y[j1])/(12*mi)

    for i in range(1, nx):
        psi[i][ny] = psi[0][ny]

    for i in range(i1, nx):
        psi[i][0] = psi[0][j1]

    for j in range(1, j1+1):
        psi[i1][j] = psi[0][j1]

    for i in range(1, i1+1):
        psi[i][j1] = psi[0][j1]


def wb_dzeta(dzeta, psi, Qwe, mi, i1, j1, y, nx, ny, delta):
    for j in range(j1, ny+1):
        dzeta[0][j] = (Qwe/(2*mi))*(2*y[j]-y[j1]-y[ny])

    for j in range(ny+1):
        dzeta[nx][j] = (Qwy(Qwe, y, j1, ny)/(2*mi))*(2*y[j]-y[ny])

    for i in range(1, nx):
        dzeta[i][ny] = (2/(delta**2))*(psi[i][ny-1] - psi[i][ny])

    for i in range(i1, nx):
        dzeta[i][0] = (2/(delta**2))*(psi[i][1] - psi[i][0])

    for j in range(1, j1):
        dzeta[i1][j] = (2/(delta**2))*(psi[i1+1][j] - psi[i1][j])

    for i in range(1, i1+1):
        dzeta[i][j1] = (2/(delta**2))*(psi[i][j1+1] - psi[i][j1])

    dzeta[i1][j1] = 0.5*(dzeta[i1-1][j1]-dzeta[i1][j1-1])


def tau(psi, dzeta, nx, delta, j2):
    suma = 0
    for i in range(1, nx):
        suma += psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4*psi[i][j2] - dzeta[i][j2]*delta**2
    return suma


def relaksacja1(psi, dzeta, mi, ro,  IT_MAX, y, nx, ny, delta, i1, j1, Qwe):
    wb_psi(psi, Qwe, mi, i1, j1, y, nx, ny)
    for it in range(1, IT_MAX + 1):
        for i in range(1, nx):
            for j in range(1, ny):
                if not (0 <= i <= i1  and 0 <= j <= j1) :
                    psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-(delta**2)*dzeta[i][j])
                    dzeta[i][j] = 0.25*(dzeta[i+1][j]+dzeta[i-1][j]+dzeta[i][j+1]+dzeta[i][j-1])
        wb_dzeta(dzeta, psi, Qwe, mi, i1, j1, y, nx, ny, delta)
        print("Q = ",Qwe," Iteracja: ", it, " Kontrola bledu: ", tau(psi, dzeta, nx, delta, j1+2))


def relaksacja2(psi, dzeta, mi, ro,  IT_MAX, y, nx, ny, delta, i1, j1, Qwe):
    wb_psi(psi, Qwe, mi, i1, j1, y, nx, ny)
    for it in range(1, IT_MAX + 1):
        if it < 2000:
            omega = 0
        else:
            omega = 1

        for i in range(1, nx):
            for j in range(1, ny):
                if not (0 <= i <= i1  and 0 <= j <= j1) :
                    psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-(delta**2)*dzeta[i][j])
                    dzeta[i][j] = 0.25 * (dzeta[i + 1][j] + dzeta[i - 1][j] + dzeta[i][j + 1] + dzeta[i][j - 1]) - (ro / (16 * mi)) * ((round(psi[i][j + 1], 4) - round(psi[i][j - 1], 4)) * (round(dzeta[i + 1][j], 4) - round(dzeta[i - 1][j], 4)) - (round(dzeta[i][j + 1], 4) - round(dzeta[i][j - 1], 4)) * (round(psi[i + 1][j], 4) - round(psi[i - 1][j],4)))
        wb_dzeta(dzeta, psi, Qwe, mi, i1, j1, y, nx, ny, delta)
        print("Q = ",Qwe," Iteracja: ", it, " Kontrola bledu: ", tau(psi, dzeta, nx, delta, j1+2))


def relaksacja3(psi, dzeta, mi, ro,  IT_MAX, y, nx, ny, delta, i1, j1, Qwe):
    wb_psi(psi, Qwe, mi, i1, j1, y, nx, ny)
    for it in range(1, IT_MAX + 1):
        if it < 20000:
            omega = 0
        else:
            omega = 1

        for i in range(1, nx):
            for j in range(1, ny):
                if not (0 <= i <= i1  and 0 <= j <= j1) :
                    psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-(delta**2)*dzeta[i][j])
                    dzeta[i][j] = 0.25*(dzeta[i+1][j]+dzeta[i-1][j]+dzeta[i][j+1]+dzeta[i][j-1]) - (omega*ro/(16*mi))*((psi[i][j+1]-psi[i][j-1])*(dzeta[i+1][j]-dzeta[i-1][j]) - (dzeta[i][j+1]-dzeta[i][j-1])*(psi[i+1][j]-psi[i-1][j]))
        wb_dzeta(dzeta, psi, Qwe, mi, i1, j1, y, nx, ny, delta)
        print("Q = ",Qwe," Iteracja: ", it, " Kontrola bledu: ", tau(psi, dzeta, nx, delta, j1+2))


def stokes(delta, ro, mi, nx, ny, i1, j1, IT_MAX, Qwe, kolejka, wersja):
    psi = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    dzeta = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    u = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    v = [[0 for x in range(ny + 1)] for y in range(nx + 1)]
    y = [0 for j in range(ny + 1)]

    for i in range(ny + 1):
        y[i] = i*delta

    if wersja == 1:
        relaksacja1(psi, dzeta, mi, ro, IT_MAX, y, nx, ny, delta, i1, j1, Qwe)
    elif wersja == 2:
        relaksacja2(psi, dzeta, mi, ro, IT_MAX, y, nx, ny, delta, i1, j1, Qwe)
    elif wersja == 3:
        relaksacja3(psi, dzeta, mi, ro, IT_MAX, y, nx, ny, delta, i1, j1, Qwe)

    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            if j1 < j < ny and 0 < i < i1:
                u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * delta)
                v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * delta)
            if i1 < i < nx and 0 < j < j1:
                u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * delta)
                v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * delta)
            if i1 <= i < nx and j1 <= j < ny:
                u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * delta)
                v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * delta)

    kolejka.put([psi, dzeta, u, v, Qwe])


def podpunkt1(list, q):
    global iter
    for x in list:
        x.start()
    while iter != 3:
        result = q.get()
        ploter(result)


def ploter(list):
    global iter
    ploter1(list[0], "Q="+str(list[4])+" psi(x,y)", str(list[4])+" psi(x,y)")
    ploter1(list[1], "Q="+str(list[4])+" dzeta(x,y)", str(list[4])+" dzeta(x,y)")
    ploter2(list[2], "Q="+str(list[4])+" u(x,y)", str(list[4])+" u(x,y)")
    ploter2(list[3], "Q="+str(list[4])+" v(x,y)", str(list[4])+" v(x,y)")
    iter += 1


def ploter1(t,title, file):
    plt.contour(t, cmap='gnuplot')
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.ylabel("y")
    plt.xlabel("x")
    plt.savefig(file+".png")
    # plt.show()
    plt.clf()


def ploter2(t,title, file):
    plt.imshow(t, cmap='jet')
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.ylabel("x")
    plt.xlabel("y")
    plt.savefig(file+".png")
    # plt.show()
    plt.clf()


def main(args):
    q = multiprocessing.Queue()
    list1 = [multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 20000, -1000, q, 1)),
            multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 20000, -4000, q, 1)),
            multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 20000, 4000, q, 1))]

    list2 = [multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 20000, -1000, q, 2)),
            multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 20000, -4000, q, 2)),
            multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 20000, 4000, q, 2))]

    list3 = [multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 25000, -1000, q, 3)),
             multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 25000, -4000, q, 3)),
             multiprocessing.Process(None, stokes, args=(0.01, 1, 1, 200, 90, 50, 55, 25000, 4000, q, 3))]

    podpunkt1(list1, q)
    # podpunkt1(list2, q)
    # podpunkt1(list3, q)



if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
