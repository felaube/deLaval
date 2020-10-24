import math
import numpy as np
import matplotlib.pyplot as plt

"""
Definição inicial de constantes

GAMMA [-]: Razão de calores específicos do gás
M_TEST [-]: Número de Mach na seção de testes
Y_0 [-]: Raio na garganta do bocal, normalizado pelo comprimento do bocal
"""

GAMMA = 1.4
M_TEST = 4
Y_0 = 5.07
Y_1 = 10.14


def main():
    v_test_rad, v_test_deg = v(M_TEST)

    A_test_ratio = A_ratio(M_TEST)

    theta_1_max_rad = v_test_rad/2
    theta_1_max_deg = v_test_deg/2

    theta_1_rad, theta_1_deg = empirical_theta_1(v_test_rad, A_test_ratio)

    print(theta_1_max_deg)
    print(theta_1_deg)

    x_initial_curve, y_initial_curve = initial_curve(Y_0, Y_1, theta_1_rad)
    y_initial_curve_mirror = -y_initial_curve
    separating_dashed_line = np.zeros(x_initial_curve.shape)

    fig, ax = plt.subplots()
    ax.plot(x_initial_curve, y_initial_curve)
    ax.plot(x_initial_curve, y_initial_curve_mirror)
    ax.plot(x_initial_curve, separating_dashed_line, 'k--')
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")

    # plt.grid()
    plt.xlim(-25, 100)
    plt.xticks(np.arange(-25, 101, 5))
    plt.yticks(np.arange(-30, 31, 5))
    plt.show()


def empirical_theta_1(v_t, A_ratio=0, A_star=0, A_test=0):
    """
    Retorna o valor do parâmetro theta_1,
    calculado por meio de uma relação empírica,
    que utiliza a área na garganta (A_star),
    a área na seção de teste (A_test) e o ângulo
    de expansão na seção de teste (v_t).

    A_ratio espera a relação A/A*
    """

    if A_ratio == 0:
        A_ratio = A_test/A_star

    theta_1 = (1/A_ratio) ** (2/9) * (v_t/2)

    return (theta_1, math.degrees(theta_1))


def v(M, gamma=1.4):
    """
    Função de Prandtl-Meyer.

    Dado o número de Mach do escoamento, retorna o ângulo de expansão
    como uma tupla,em que o primeiro elemento é dado em radianos
    e o segundo em graus.
    """
    part_1 = math.sqrt((gamma + 1)/(gamma - 1))
    part_2 = math.atan(math.sqrt((gamma - 1)/(gamma + 1)*(M**2 - 1)))
    part_3 = math.atan(math.sqrt(M**2 - 1))

    number = part_1*part_2 - part_3

    return (number, math.degrees(number))


def A_ratio(M, gamma=1.4):
    """
    Retorna a relação de áreas A/A*
    """
    part_1 = 1/M
    part_2 = math.sqrt(((2 + (gamma - 1)*M ** 2)/(gamma + 1)))
    part_2 = part_2 ** ((gamma + 1)/(gamma - 1))

    return part_1*part_2


def initial_curve(y_0, y_1, theta_1):
    """
    Retorna uma lista com os pontos da parte do bocal
    da parte sônica ao ponto de inflexão,
    respeitando os requisitos de razão de área,
    inclinação nula na garganta, e valores de y, inclinação
    e curvatura iguais às da curva terminal (parte divergente)
    no ponto de inflexão.

    theta_1 [rad]
    """
    x_1 = 3/(2*math.tan(theta_1))*(y_1 - y_0)
    x = np.arange(-25, x_1, 0.01)
    y = y_0 + math.tan(theta_1)/x_1*(x ** 2)*(1 - x/(3*x_1))

    return (x, y)


if __name__ == "__main__":
    main()
