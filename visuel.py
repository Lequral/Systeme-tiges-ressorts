import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# Données du modèle
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
