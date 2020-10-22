import math

"""
Definição inicial de constantes

GAMMA [-]: Razão de calores específicos do gás
M_TEST [-]: Número de Mach na seção de testes
Y_0 [-]: Raio na garganta do bocal, normalizado pelo comprimento do bocal
"""

GAMMA = 1.4
M_TEST = 4
Y_0 = 5


def main():
    theta_1_max = v(M_TEST)[1]/2

    print(theta_1_max)


def theta_1(A_star, A_test, v_t):
    """
    Retorna o valor do parâmetro theta_1,
    calculado por meio de uma relação empírica,
    que utiliza a área na garganta (A_star),
    a área na seção de teste (A_test) e o ângulo
    de expansão na seção de teste (v_t)
    """

    return (A_star/A_test) ** (2/9) * (v_t/2)


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


if __name__ == "__main__":
    main()
