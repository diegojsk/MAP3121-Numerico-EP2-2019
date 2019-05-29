import numpy as np
from main import *
import matplotlib.pyplot as plt


def X_gab(t):
        X = np.array([
            np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t),
            np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t),
            -np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t),
            -np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t)])

        return X


if __name__ == "__main__":

    np.set_printoptions(precision=3, suppress=True)
    A = np.array([[-2, -1, -1, -2],
                  [1, -2,  2, -1],
                  [-1, -2, -2, -1],
                  [2, -1,  1, -2]]).astype(np.double)

    X_0 = np.array([1, 1, 1, -1])

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

    print(output)
