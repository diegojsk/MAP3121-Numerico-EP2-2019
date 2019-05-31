from main import *
import numpy as np
from math import sin, pi
import matplotlib.pyplot as plt



if __name__ == "__main__":

    g = 9.8
    l = 0.8
    w = np.sqrt(g/l)
    u = 2*w - 0.5

    Emax = 10**(-10)

    def F(t,y):
        """
        Sistema de equações do pêndulo
        """
        return np.array([ y[1], -(np.power(w,2))*sin(y[0]) - u*y[1]])

    X_0 = np.array([(5*pi)/2, 3])


    output, ts = runge_kutta(F, X_0, 0.01 , 0, 5)

    fig1, ax1 = plt.subplots()
    ax1.plot(ts, output)
    ax1.set_title("Resolução pelo método de Runge-Kutta")
    plt.show()

    # Agora com u = 0

    u = 0

    output, ts = runge_kutta(F, X_0, 0.01 , 0, 5)

    fig2, ax1 = plt.subplots()
    ax1.plot(ts, output)
    ax1.set_title("Resolução pelo método de Runge-Kutta")
    plt.show()

    output, ts = runge_kutta_automatico_pendulo(F, X_0, 0.01 , 0, 5,Emax, w)

    fig3, ax1 = plt.subplots()
    ax1.plot(ts, output)
    ax1.set_title("Resolução pelo método de Runge-Kutta")
    plt.show()







