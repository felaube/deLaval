import math
import decimal
import numpy as np
import matplotlib.pyplot as plt
from kaleido.scopes.plotly import PlotlyScope
import plotly.graph_objects as go
from plotly.subplots import make_subplots
"""
Definição inicial de constantes

GAMMA [-]: Razão de calores específicos do gás
M_TEST [-]: Número de Mach na seção de testes
Y_0 [-]: Raio na garganta do bocal, normalizado pelo comprimento do bocal
"""

GAMMA = 1.4
M_TEST = 4
Y_0 = 10


def main():
    v_test_rad, v_test_deg = v(M_TEST)

    A_test_ratio = A_ratio(M_TEST)

    # theta_1_max_rad = v_test_rad/2
    # theta_1_max_deg = v_test_deg/2

    theta_1_rad, theta_1_deg = empirical_theta_1(v_test_rad, A_test_ratio)

    v_1_rad = v_test_rad - theta_1_rad
    v_1_deg = v_test_deg - theta_1_deg

    M_1 = M_from_v(v_1_rad)

    x_final_curve, y_final_curve = foelsch(Y_0, theta_1_rad, v_1_rad,
                                           v_test_rad, M_1)

    y_1 = y_final_curve[0]

    x_initial_curve, y_initial_curve = initial_curve(Y_0, y_1, theta_1_rad)

    fig, ax = plt.subplots()
    ax.plot(x_initial_curve, y_initial_curve, 'r-')
    ax.plot(x_initial_curve, -y_initial_curve, 'r-')
    ax.plot(x_final_curve, y_final_curve, 'r-')
    ax.plot(x_final_curve, -y_final_curve, 'r-')
    ax.plot(x_initial_curve, np.zeros(x_initial_curve.shape), 'k--')
    ax.plot(x_final_curve, np.zeros(x_final_curve.shape), 'k--')
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")

    # plt.grid()
    plt.xlim(x_initial_curve[0], x_final_curve[-1]*1.05)
    plt.ylim(-y_final_curve[-1]*2, y_final_curve[-1]*2)
    # plt.xticks(np.arange(x_initial_curve[0], x_final_curve[-1]*1.05, 5))
    # plt.yticks(np.arange(-y_final_curve[-1]*1.25, y_final_curve[-1]*1.25, 5))
    plt.show()

    # Quasi Unidimensional
    graph_y_points = list(float_range(0, max(y_final_curve),max(y_final_curve)/len(y_final_curve)))
    mach_matriz = []
    for i in range(len(graph_y_points)):
        mach_matriz.append(add_line_on_quasi_graph( graph_y_points[i], y_final_curve, len(x_final_curve)))

    scope = PlotlyScope()
    fig_quasi_unidimensional = make_subplots(rows=1, cols=1)
    fig_quasi_unidimensional.add_trace(go.Heatmap(z=mach_matriz, colorbar=dict(title='Mach', titleside='right'), connectgaps=True, zsmooth='best'), 1, 1)
    fig_quasi_unidimensional.update
    fig_quasi_unidimensional.update_yaxes(title_text="Y axis (mm)", row=1, col=1)
    fig_quasi_unidimensional.update_xaxes(title_text="X axis (mm)", row=1, col=1)
    with open("quasi_unidimensional_graph.png", "wb") as f:
        f.write(scope.transform(fig_quasi_unidimensional, format="png"))


def empirical_theta_1(v_t, A_test_ratio=0, A_star=0, A_test=0):
    """
    Retorna o valor do parâmetro theta_1,
    calculado por meio de uma relação empírica,
    que utiliza a área na garganta (A_star),
    a área na seção de teste (A_test) e o ângulo
    de expansão na seção de teste (v_t).

    A_test_ratio - relação A/A*.
    v_t - [rad]
    """

    if A_test_ratio == 0:
        A_test_ratio = A_test/A_star

    theta_1 = ((1/A_test_ratio) ** (2/9)) * (v_t/2)

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
    da parte subsônica ao ponto de inflexão,
    respeitando os requisitos de razão de área,
    inclinação nula na garganta, e valores de y, inclinação
    e curvatura iguais às da curva terminal (parte divergente)
    no ponto de inflexão.

    theta_1 [rad]
    """
    x_1 = 3/(2*math.tan(theta_1))*(y_1 - y_0)
    x = np.arange(-50, x_1, 0.01)
    y = y_0 + math.tan(theta_1)/x_1*(x ** 2)*(1 - x/(3*x_1))

    return (x, y)


def foelsch(y_0, theta_1, v_1, v_test, M_1):
    """
    Calcula os pontos da curva final (após o ponto de inflexão do bocal),
    utilizando o Método de Foelsch.

    theta_1 - [rad]
    v_1 - [rad]
    v_test - [rad]
    """

    r_0 = y_0/theta_1

    A_1_ratio = A_ratio(M_1)

    y_1 = r_0*A_1_ratio*math.sin(theta_1)
    x_1 = 3/(2*math.tan(theta_1))*(y_1 - y_0)

    r_1 = y_1/math.sin(theta_1)

    M_array = np.arange(M_1, M_TEST, 0.01)
    x_array = np.zeros(M_array.shape)
    y_array = np.zeros(M_array.shape)

    for i, M in enumerate(M_array):
        a = math.asin(1/M)

        A_ratio_local = A_ratio(M)
        r_local = r_0*A_ratio_local
        v_local, _ = v(M)
        l_local = M*r_local*(v_local - v_1)

        x_2 = (x_1 - r_1*math.cos(theta_1) +
               r_local*math.cos(v_test - v_local))
        y_2 = r_local*math.sin(v_test - v_local)

        x_array[i] = x_2 + l_local*math.cos(v_test - v_local + a)
        y_array[i] = y_2 + l_local*math.sin(v_test - v_local + a)

    return (x_array, y_array)

def M_from_y(y_1, y_0 = Y_0, gamma=1.4):
    """
    Calcula o número de Mach dado o x em relação ao ponto inicial em mm,
    de forma iterativa. Assume-se intervalos para o número de Mach,
    que são atualizados e reduzidos a cada vez que há uma troca de sinal
    na subtração entre o valor estimado e o valor 'A_ratio_target'.

    y_1 - [mm]
    """
    lower_boundary = 1
    upper_boundary = 1000

    tolerance = 1e-6
    A_ratio_target = (y_1/y_0) ** 2

    while True:
        interval = np.linspace(lower_boundary, upper_boundary, 10)

        previous_est = math.inf
        current_est = math.inf

        for i, M in enumerate(interval):
            current_A_ratio = A_ratio(M)

            current_est = A_ratio_target - current_A_ratio

            if abs(current_est) <= tolerance:
                return M
            elif i != 0 and current_est*previous_est < 0:
                # Atualiza fronteiras de busca
                # Regra do sanduíche (se houve mudança de sinal entre dois
                # pontos, o zero da função está entre esses dois pontos)
                lower_boundary = interval[i-1]
                upper_boundary = M
                break
            previous_est = current_est




def M_from_v(v_target, gamma=1.4):
    """
    Calcula o número de Mach dado o ângulo de expansão em radianos,
    de forma iterativa. Assume-se intervalos para o número de Mach,
    que são atualizados e reduzidos a cada vez que há uma troca de sinal
    na subtração entre o valor estimado e o valor 'v_target'.

    v_target - [rad]
    """
    lower_boundary = 1
    upper_boundary = 1000

    tolerance = 1e-6

    while True:
        interval = np.linspace(lower_boundary, upper_boundary, 10)

        previous_est = math.inf
        current_est = math.inf

        for i, M in enumerate(interval):
            current_v, _ = v(M)

            current_est = v_target - current_v

            if abs(current_est) <= tolerance:
                return M
            elif i != 0 and current_est*previous_est < 0:
                # Atualiza fronteiras de busca
                # Regra do sanduíche (se houve mudança de sinal entre dois
                # pontos, o zero da função está entre esses dois pontos)
                lower_boundary = interval[i-1]
                upper_boundary = M
                break
            previous_est = current_est

def add_line_on_quasi_graph(graph_y_point, y_curve, index):
    graph_line = []
    for _, y_point in enumerate(y_curve):
        if (y_point < graph_y_point):
            graph_line.append(0)
        else:
            graph_line.append(M_from_y(y_point))
    
    return graph_line



def float_range(start, stop, step):
  while start < stop:
    yield float(start)
    start += decimal.Decimal(step)

if __name__ == "__main__":
    main()
