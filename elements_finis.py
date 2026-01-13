import numpy as np
from scipy.optimize import root

r_11, r_12, r_21, r_22 = 1, 1, 1, 1
m_1, m_2, m_3, g = 1, 1, 1, 9.81
a, b, r_3, L = 1, 2, 1, 2
k_1, k_2, k_3 = 10, 10, 10

k_spring = 100  # raideur du ressort A2-C
L_rod = 2.0    # longueur A1-O1-A2
L_spring0 = 1.0

F_1 = np.array([0.0, 1.0])  # force sur A1 (verticale à changer)

O1 = np.array([0.0, 0.0])
C = np.array([1.0, -1.0])

x_A1, y_A1, x_A2, y_A2 = -1, 0, 1, 0
x = [x_A1, y_A1, x_A2, y_A2]


# def residual(x):
#     A1 = x[0:2]
#     A2 = x[2:4]

#     R1 = np.zeros(2)
#     R2 = np.zeros(2)

#     # --------------------
#     # Tige O1-A2 (ressort raide)
#     # --------------------
#     d = A2 - O1
#     L = np.linalg.norm(d)
#     u = d / L

#     F_rod = k_rod * (L - L_rod) * u

#     R2 -= F_rod  # force sur A2
#     # O1 est fixe → pas d'équation

#     # --------------------
#     # Contrainte A1 sur la tige
#     # --------------------
#     # Projection de A1 sur la droite O1-A2
#     t = np.dot(A1 - O1, u)
#     A1_proj = O1 + t * u

#     F_constraint = k_rod * (A1 - A1_proj)

#     R1 -= F_constraint
#     R2 += F_constraint * (t / L)

#     # --------------------
#     # Ressort A2-C
#     # --------------------
#     d = A2 - C
#     L = np.linalg.norm(d)
#     u = d / L

#     F_spring = k_spring * (L - L_spring0) * u
#     R2 -= F_spring

#     # --------------------
#     # Force externe sur A1
#     # --------------------
#     R1 += F_ext

#     return np.hstack([R1, R2])


# x0 = np.array([
#     -1.0, 0.0,   # A1 initial
#     1.0, 0.0    # A2 initial
# ])

# sol = root(residual, x0)

# A1, A2 = sol.x[0:2], sol.x[2:4]

# print("A1 =", A1)
# print("A2 =", A2)
