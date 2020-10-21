import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

for i in range(1,4):
    multiplier = 0
    if( i == 1):
        multiplier = 1
    elif i == 2:
        multiplier = 10
    elif i == 3:
        multiplier = 100
    file = "euler" + str(i) + ".txt"
    analit = np.loadtxt("eulera1.txt")
    data = np.loadtxt(file)

    x = data[:, 0]*multiplier
    y = data[:, 1]

    ax = analit[:, 0]
    ay = analit[:, 1]

    plt.plot(x, y, 'o', label = 'dt = '+ str(0.01 * multiplier))

plt.plot(ax, ay, 'b-', label = 'analitycznie')
plt.legend(loc = 'upper right')

plt.savefig("euler.png")
plt.clf()

for i in range(1,4):
    multiplier = 0
    if( i == 1):
        multiplier = 1
    elif i == 2:
        multiplier = 10
    elif i == 3:
        multiplier = 100
    file = "eulerb" + str(i) + ".txt"
    data = np.loadtxt(file)

    x = data[:, 0]*multiplier
    y = data[:, 1]

    plt.plot(x, y, '-', label = 'dt = '+ str(0.01 * multiplier))

plt.legend(loc = 'lower right')

plt.savefig("euler_blad.png")
plt.clf()

####################################################################

for i in range(1,4):
    multiplier = 0
    if( i == 1):
        multiplier = 1
    elif i == 2:
        multiplier = 10
    elif i == 3:
        multiplier = 100
    file = "RK2" + str(i) + ".txt"
    analit = np.loadtxt("RK2a1.txt")
    data = np.loadtxt(file)

    x = data[:, 0]*multiplier
    y = data[:, 1]

    ax = analit[:, 0]
    ay = analit[:, 1]

    plt.plot(x, y, 'o', label = 'dt = '+ str(0.01 * multiplier))

plt.plot(ax, ay, 'b-', label = 'analitycznie')
plt.legend(loc = 'upper right')

plt.savefig("RK2.png")
plt.clf()

for i in range(1,4):
    multiplier = 0
    if( i == 1):
        multiplier = 1
    elif i == 2:
        multiplier = 10
    elif i == 3:
        multiplier = 100
    file = "RK2b" + str(i) + ".txt"
    data = np.loadtxt(file)

    x = data[:, 0]*multiplier
    y = data[:, 1]

    plt.plot(x, y, '-', label = 'dt = '+ str(0.01 * multiplier))

plt.legend(loc = 'lower right')

plt.savefig("RK2_blad.png")
plt.clf()

####################################################################

for i in range(1,4):
    multiplier = 0
    if( i == 1):
        multiplier = 1
    elif i == 2:
        multiplier = 10
    elif i == 3:
        multiplier = 100
    file = "RK4" + str(i) + ".txt"
    analit = np.loadtxt("RK4a1.txt")
    data = np.loadtxt(file)

    x = data[:, 0]*multiplier
    y = data[:, 1]

    ax = analit[:, 0]
    ay = analit[:, 1]

    plt.plot(x, y, 'o', label = 'dt = '+ str(0.01 * multiplier))

plt.plot(ax, ay, 'b-', label = 'analitycznie')
plt.legend(loc = 'upper right')

plt.savefig("RK4.png")
plt.clf()

for i in range(1,4):
    multiplier = 0
    if( i == 1):
        multiplier = 1
    elif i == 2:
        multiplier = 10
    elif i == 3:
        multiplier = 100
    file = "RK4b" + str(i) + ".txt"
    data = np.loadtxt(file)

    x = data[:, 0]*multiplier
    y = data[:, 1]

    plt.plot(x, y, '-', label = 'dt = '+ str(0.01 * multiplier))

plt.legend(loc = 'lower right')

plt.savefig("RK4_blad.png")
plt.clf()