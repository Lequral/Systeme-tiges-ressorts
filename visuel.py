import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import calcul_analytique
from scipy.optimize import least_squares
import sympy as sp


def calcul(val=None):
    ax1.clear()

    w1_val = slider_w1.val
    w2_val = slider_w2.val
    biais_val = slider_biais.val
    F1_val = slider_F1.val
    F2_val = slider_F2.val

    # Donnée géométrique
    weights = [w1_val, w2_val, biais_val]
    sol = calcul_analytique.trouver_param_real(weights)
    donnée_géométrie = sol.x

    r_11, r_12, r_21, r_3, k_1, k_2, k_3, b = donnée_géométrie
    r_22 = r_12

    dico = calcul_analytique.equations(False)

    r_11_ex, r_12_ex, r_21_ex, r_22_ex, m_1_ex, m_2_ex, m_3_ex, g_ex, a_ex, b_ex, r_3_ex, L_ex, theta10_ex, theta20_ex, gamma10_ex, theta1_ex, theta2_ex, gamma1_ex, X01_ex, Y01_ex, F1_ex, T1_ex, X02_ex, Y02_ex, F2_ex, T2_ex, X03_ex, Y03_ex, k1_ex, k2_ex, k3_ex = dico[
        "param"]

    equation = dico["equation"]

    equation = [e.subs({
        r_11_ex: r_11,
        r_12_ex: r_12,
        r_21_ex: r_21,
        r_22_ex: r_22,
        g_ex: 9.81,
        a_ex: r_12+r_3,
        b_ex: b,
        r_3_ex: r_3,
        L_ex: 1,
        theta10_ex: 0,
        theta20_ex: 0,
        gamma10_ex: 0,
        F1_ex: F1_val,
        F2_ex: F2_val,
        k1_ex: k_1,
        k2_ex: k_2,
        k3_ex: k_3,
        m_1_ex: 1,
        m_2_ex: 1,
        m_3_ex: 1
    }) for e in equation]

    for e in equation:
        print(e)

    param = (theta1_ex, theta2_ex, gamma1_ex, X01_ex, Y01_ex,
             T1_ex, X02_ex, Y02_ex, T2_ex, X03_ex, Y03_ex)

    func = sp.lambdify(param, equation, 'numpy')

    def numpy_equation(vars):
        theta1, theta2, gamma1, X01, Y01, T_1, X02, Y02, T_2, X03, Y03 = vars
        return func(theta1, theta2, gamma1, X01, Y01, T_1, X02, Y02, T_2, X03, Y03)

    guess = [1 for i in range(11)]
    min_var = [0, 0, 0] + [-100 for j in range(8)]
    max_var = [2*np.pi, 2*np.pi, 2*np.pi] + [100 for j in range(8)]

    solu_param = least_squares(numpy_equation, guess, bounds=(
        min_var, max_var), max_nfev=1000)

    theta1 = solu_param.x[0]
    theta2 = solu_param.x[1]
    gamma1 = solu_param.x[2]

    point = np.array([
        [0.0, 0.0],  # 01 et G1
        [-r_11*np.cos(theta1), -r_11*np.sin(theta1)],  # A1
        [r_12*np.cos(theta1), r_12*np.sin(theta1)],  # B1
        [0.0, -2*b],  # 02 et G2
        [-r_21*np.cos(theta2), -2*b - r_21*np.sin(theta2)],  # A2
        [r_22*np.cos(theta2), -2*b + r_22*np.sin(theta2)],  # B2
        [abs(r_12 + r_3), -b],  # D
        [abs(r_12+r_3)-r_3*np.cos(gamma1), -b-r_3*np.sin(gamma1)]  # C
    ])

    barre = [
        (0, 1),  # O1A1
        (0, 2),  # O1B1
        (3, 5),  # O2B2
        (3, 4),  # O2A2
        (6, 7),  # DC
        (2, 7),  # B1C
        (5, 7)  # B2C

    ]

    def draw_torsion_spring(ax1, center, radius=0.1, turns=4, n=200, color="tab:green"):
        cx, cy = center
        theta = np.linspace(0, 2 * np.pi * turns, n)
        r = radius * (1 + 0.1 * np.sin(6 * theta))
        x = cx + r * np.cos(theta)
        y = cy + r * np.sin(theta)
        ax1.plot(x, y, color=color, linewidth=2)

    # barres
    for i, j in barre:
        ax1.plot(
            [point[i, 0], point[j, 0]],
            [point[i, 1], point[j, 1]],
            color="tab:blue",
            linewidth=3,
        )

    colors = ["blue", "blue", "green",
              "green", "purple", "yellow", "yellow"]

    for (i, j), c in zip(barre, colors):
        ax1.plot(
            [point[i, 0], point[j, 0]],
            [point[i, 1], point[j, 1]],
            color=c,
            linewidth=3,
        )

    # ressort
    draw_torsion_spring(ax1, (0, 0))
    draw_torsion_spring(ax1, (0.0, -2*b))
    draw_torsion_spring(ax1, (0.0, -2*b))

    # Noeuds
    ax1.scatter(point[:, 0], point[:, 1], color="black", zorder=3)

    # Index des noeuds (optionnel)
    for k, (x, y) in enumerate(point):
        node_labels = ["O1", "A1", "B1", "O2", "A2", "B2", "D", "C"]
        ax1.text(x + 0.02, y + 0.02, node_labels[k], fontsize=12)

    # Mise en forme
    ax1.set_aspect("equal")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_title("Structure 2D – tiges et ressorts (statique)")
    ax1.grid(True)
    ax1.set_xlim(-3, 3)
    ax1.set_ylim(-3, 3)

    fig1.canvas.draw_idle()
    print(
        f"theta 1 vaut {theta1}rad, theta 2 vaut {theta2}rad et l'output gamma 1 vaut {gamma1}")


# Choix poids et biais
w1_val = 0
ax_slider_w1 = plt.axes([0.2, 0.1, 0.6, 0.03])
slider_w1 = Slider(ax_slider_w1, 'w1', -1.0, 2.0, valinit=3.2234)


w2_val = 0
ax_slider_w2 = plt.axes([0.2, 0.2, 0.6, 0.03])
slider_w2 = Slider(ax_slider_w2, 'w2', -1.0, 2.0, valinit=3.0913)


biais_val = 0
ax_slider_biais = plt.axes([0.2, 0.3, 0.6, 0.03])
slider_biais = Slider(ax_slider_biais, 'biais', -1.0, 2.0, valinit=-0.0407)


F1_val = 0
ax_slider_F1 = plt.axes([0.2, 0.4, 0.6, 0.03])
slider_F1 = Slider(ax_slider_F1, 'F1', -1.0, 1.0, valinit=0)


F2_val = 0
ax_slider_F2 = plt.axes([0.2, 0.5, 0.6, 0.03])
slider_F2 = Slider(ax_slider_F2, 'F2', -1.0, 1.0, valinit=0)


# Tracer le graphique
fig1, ax1 = plt.subplots(figsize=(6, 6))

slider_w1.on_changed(calcul)
slider_w2.on_changed(calcul)
slider_biais.on_changed(calcul)
slider_F1.on_changed(calcul)
slider_F2.on_changed(calcul)

calcul()

plt.show()
