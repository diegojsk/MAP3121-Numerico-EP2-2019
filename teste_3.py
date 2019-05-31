from main import *
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    n = 10
    A = np.zeros((n, n))

    for i in range(n):
        A[i, i] = -2

    for i in range(n-1):
        A[i, i+1] = 1
        A[i+1, i] = 1

    print(A)

    def y(i):
        return i/(n + 1)

    X_0 = np.array([ np.sin(np.pi*y(i)) + np.sin(n*np.pi*y(i)) for i in range(1, n+1) ])

    lambda_1 = 2*(1 - np.cos(np.pi/(n+1)))

    lambda_2 = 2*(1 - np.cos(n*np.pi/(n+1)))

    def X_gab(t):
        return [np.exp(-lambda_1*t)*np.sin(np.pi*y(i)) +
                np.exp(-lambda_2*t)*np.sin(n*np.pi*y(i)) for i in range(1, n+1) ]

    def F(t, x):
        return np.matmul(A, x)

    output, ts = runge_kutta(F, X_0, 0.01, 0, 2)

    fig1, ax1 = plt.subplots()
    ax1.plot(ts, output)
    ax1.set_title("Resolução pelo método de Runge-Kutta")

    output2 = [X_gab(t) for t in ts]
    fig2, ax2 = plt.subplots()
    ax2.plot(ts, output2)
    ax2.set_title("Aplicação da solução analítica")

    plt.show()
