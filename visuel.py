import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import calcul_analytique

# ------------------------
# exemple
# ------------------------
nodes = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
])

rods = [
    (0, 1),
    (1, 2),
]

springs = [
    (0, 2),
]

# ------------------------
# Fonction de dessin d'un ressort
# ------------------------
def draw_spring(ax, p1, p2, n=20, amp=0.05):
    x1, y1 = p1
    x2, y2 = p2

    dx, dy = x2 - x1, y2 - y1
    L = np.hypot(dx, dy)

    t = np.linspace(0, 1, n)
    nx, ny = -dy / L, dx / L  # vecteur normal

    zigzag = amp * np.sin(2 * np.pi * t * (n // 2))

    x = x1 + dx * t + nx * zigzag
    y = y1 + dy * t + ny * zigzag

    ax.plot(x, y, color="tab:orange", linewidth=2)

# ------------------------
# Rendu
# ------------------------
fig, ax = plt.subplots(figsize=(6, 6))

# Tiges (segments droits)
for i, j in rods:
    ax.plot(
        [nodes[i, 0], nodes[j, 0]],
        [nodes[i, 1], nodes[j, 1]],
        color="tab:blue",
        linewidth=3,
    )

# Ressorts
for i, j in springs:
    draw_spring(ax, nodes[i], nodes[j])

# Noeuds
ax.scatter(nodes[:, 0], nodes[:, 1], color="black", zorder=3)

# Index des noeuds (optionnel)
for k, (x, y) in enumerate(nodes):
    ax.text(x + 0.02, y + 0.02, f"{k}", fontsize=12)

# Mise en forme
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Structure 2D – tiges et ressorts (statique)")
ax.grid(True)

plt.show()


# Notre modèle

# Choix poids et biais
w1_val = 0
ax_slider_w1 = plt.axes([0.2, 0.1, 0.6, 0.03])
slider_w1 = Slider(ax_slider_w1, 'w1', -1.0, 1.0, valinit=0)

def update_w1(val):
    global w1_val 
    w1_val = slider_w1.val

slider_w1.on_changed(update_w1)

w2_val = 0
ax_slider_w2 = plt.axes([0.2, 0.2, 0.6, 0.03])
slider_w2 = Slider(ax_slider_w2, 'w2',-1.0, 1.0, valinit=0)

def update_w2(val):
    global w2_val 
    w2_val = slider_w2.val

slider_w2.on_changed(update_w2)

b_val = 0
ax_slider_b = plt.axes([0.2, 0.3, 0.6, 0.03])
slider_b = Slider(ax_slider_b, 'biais', -1.0, 1.0, valinit=0)

def update_b(val):
    global b_val 
    b_val = slider_b.val

slider_b.on_changed(update_b)

#Donnée géométrique
weights = [w1_val, w2_val, b_val]
sol = calcul_analytique.trouver_param_real(weights)
donnée_géométrie = sol.x
print(donnée_géométrie)

r_11 = donnée_géométrie[0]
r_12 = donnée_géométrie[1]
r_21 = donnée_géométrie[2]
r_22 = donnée_géométrie[3]
r_3 = donnée_géométrie[4]
k_1 = donnée_géométrie[5]
k_2 = donnée_géométrie[6]
k_3 = donnée_géométrie[7]
b = donnée_géométrie[8]

dico = calcul_analytique.equations(True)
theta1 = dico["theta1_expr"] 
theta2 = dico["theta2_expr"] 
gamma1 = dico["gamma1_expr"]
O1O2_expr = dico["O1O2_expr"]

point = [
    [0.0, 0.0], #01 et G1
    [-r_11*np.cos(theta1), -r_11*np.sin(theta1)], #A1
    [r_12*np.cos(theta1), r_12*np.sin(theta1)], #B1
    [0.0, O1O2_expr], #02 et G2
    [-r_21*np.cos(theta2),O1O2_expr -r_21*np.sin(theta2) ], #A2
    [r_22*np.cos(theta2), O1O2_expr + r_22*np.sin(theta2)], #B2
    [], #D
    [] #C
]

barre = [
    (0, 1),
    (0, 2),
    (2, 5),
    (3, 5),
    (3, 4),
    (6, 7)

]

def draw_torsion_spring(ax, center, radius=0.1, turns=4, n=200, color="tab:green"):

    cx, cy = center
    theta = np.linspace(0, 2 * np.pi * turns, n)
    r = radius * (1 + 0.1 * np.sin(6 * theta))
    x = cx + r * np.cos(theta)
    y = cy + r * np.sin(theta)
    ax.plot(x, y, color=color, linewidth=2)

#Tracer le graphique 

fig1, ax1 = plt.subplots(figsize=(6, 6))

#barres
for i, j in barre:
    ax1.plot(
        [point[i, 0], point[j, 0]],
        [point[i, 1], point[j, 1]],
        color="tab:blue",
        linewidth=3,
    )

#ressort
draw_torsion_spring(ax1, (0,0))
draw_torsion_spring(ax1, (0.0, O1O2_expr))
draw_torsion_spring(ax1, ())

# Noeuds
ax1.scatter(point[:, 0], point[:, 1], color="black", zorder=3)

# Index des noeuds (optionnel)
for k, (x, y) in enumerate(point):
    ax1.text(x + 0.02, y + 0.02, f"{k}", fontsize=12)

# Mise en forme
ax1.set_aspect("equal")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("Structure 2D – tiges et ressorts (statique)")
ax1.grid(True)

plt.show()
