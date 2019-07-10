import numpy as np


def runge_kutta(f_linha, f_0, h, t_0, t_f):
    """
    Aplica método de Runge-Kutta de ordem 4 para resolução de EDO
        :param f_linha: Função na forma f'(t, x), a ser resolvida numericamente
        :param f_0: Valor inicial da função
        :param h: Tamanho do passo de integração
        :param t_0: Instante inicial
        :param t_f: Instante final
    """

    num_passos = int((t_f - t_0)/h)
    t = t_0
    x = f_0

    xs = np.empty((num_passos, x.shape[0]))
    ts = np.empty((num_passos))

    i = 0

    while num_passos > i:

        k1 = h*f_linha(t, x)
        k2 = h*f_linha(t + 0.5*h, x + 0.5*k1)
        k3 = h*f_linha(t + 0.5*h, x + 0.5*k2)
        k4 = h*f_linha(t + h, x + k3)

        x = x + (k1 + 2*k2 + 2*k3 + k4)/6
        t = t + h

        xs[i, :] = x
        ts[i] = t

        i += 1

    return xs, ts

def depuracao(f_linha, f_0, h, t_0, t_f, sol):
    """
    Realiza a depuração indicada na seção 2.3
        :param f_linha: Função na forma f'(t, x), a ser resolvida numericamente
        :param f_0: Valor inicial da função
        :param h: Tamanho do passo de integração
        :param t_0: Instante inicial
        :param t_f: Instante final
        :param sol: Solução exata do equação

    """
    print("Depurando a solução, analisando a convergencia de p")
    h_0 = h
    for i in range(10):
        h = h_0/np.power(2, i)
        x1, t1 = runge_kutta(f_linha, f_0, h, t_0, t_f)
        x2, t2 = runge_kutta(f_linha, f_0, h/2, t_0, t_f)
        A = (x1.transpose()-sol(t1))
        B = (x2.transpose()-sol(t2))
        p = np.log2(abs(np.mean(A)/np.mean(B)))
        print("p = "+"{0:.3f}".format(p))

    return None


def calc_E_pendulo(theta, theta_ponto, omega):
    """
    Calcula o parâmetro energia definido na seção 3.2 para o pêndulo
        :param x: Valor da função x(t) em um dado instante
        :param y: Valor da função y(t) em um dado instante
        :param theta theta_ponto omega: Parâmetros do caso analisado
    """
    return np.power(theta_ponto, 2)/2 - np.cos(theta)*np.power(omega, 2)

def runge_kutta_automatico_pendulo(f_linha, f_0, h, t_0, t_f, E, w):
    """
    Aplica método de Runge-Kutta de ordem 4 para resolução de EDO com controle do passo no pendulo
        :param f_linha: Função na forma f'(t, x), a ser resolvida numericamente
        :param f_0: Valor inicial da função
        :param h: Tamanho do passo de integração
        :param t_0: Instante inicial
        :param t_f: Instante final
    """

    num_passos = int((t_f - t_0)/h)
    t = t_0
    x = f_0

    xs = np.empty((num_passos, x.shape[0]))
    ts = np.empty((num_passos))

    i = 0

    passo_adequado = True

    while num_passos > i and passo_adequado:

        k1 = h*f_linha(t, x)
        k2 = h*f_linha(t + 0.5*h, x + 0.5*k1)
        k3 = h*f_linha(t + 0.5*h, x + 0.5*k2)
        k4 = h*f_linha(t + h, x + k3)

        x = x + (k1 + 2*k2 + 2*k3 + k4)/6
        t = t + h

        xs[i, :] = x
        ts[i] = t

        i += 1

        erro = calc_E_pendulo(x[0],x[1],w)

        if erro < E :
            passo_adequado = False
            print("Passo não adequado !!!")

    print("Passo adequado !")

    return xs, ts

