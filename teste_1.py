from main import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    def F(t, x):
        return 1 + np.power((x - t), 2)

    x_0 = -18.95
    t_0 = 1.05
    t_f = 3
    output = runge_kutta(F, x_0, 0.01, t_0, t_f)
    plt.plot(output[:, 1], output[:, 0], 'b.')
    plt.show()
