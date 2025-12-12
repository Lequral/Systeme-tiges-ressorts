import sympy as sp
from sympy.physics.mechanics import dynamicsymbols, Point, ReferenceFrame

sp.init_printing()
# ∑ = Barre1
# géométrie
r_11, r_12, m_1, theta10, g = sp.symbols(
    'r_11 r_12 m_1 theta10 g', real=True, positive=True, nonzero=True)
# paramètres de position
theta1 = sp.symbols('theta1', real=True)
# efforts
X01, Y01, F_1, T = sp.symbols('X01 Y01 F_1 T')
# réologie : facteur d'amortisement et raideur du ressort
k_1 = sp.symbols('k_1', real=True, positive=True, nonzero=True)

# Repères attachées au pièces en utilisant des RéférenceFrame plus adaptés à la mécanique des systèmes qu'un CoordSys3D
R0 = ReferenceFrame('R_0')
R1 = ReferenceFrame('R_1')

# positionnement des repères
# R1 se déduit de R0 par rotation d'angle alpha autour de z0=z1
R1.orient(R0, 'Axis', [theta1, R0.z])

# Points caractéristiques
O_1 = Point('O_1')  # origine
A_1 = Point('A_1')
B_1 = Point('B_1')
G_1 = Point('G_1')

# Et position associées
O_1.set_vel(R0, 0)  # O est un point de R0
O_1.set_vel(R1, 0)  # O est un point de R01

A_1.set_pos(O_1, -r_11 * R1.x)
A_1.set_vel(R1, 0)  # A est un point de R1

B_1.set_pos(O_1, r_12 * R1.x)  # position dans R0
B_1.set_vel(R1, 0)  # B est un point de R1

G_1.set_vel(R1, 0)
G_1.set_pos(O_1, 0*R1.x)

# Torseurs
R_01 = X01 * R0.x + Y01 * R0.y
M_01 = 0 * R0.z

R_P1 = -m_1*g*R1.y  # hypothèse la tige est à l'horizontale
M_P1 = G_1.pos_from(O_1).cross(R_P1)

R_Fr1 = 0*R0.x
M_Fr1 = k_1*(theta1-theta10)*R0.z

R_T = T*R1.y  # hypothèse la tension est perpendiculaire à la tige
M_T = B_1.pos_from(O_1).cross(R_T)

R_F1 = F_1*R1.y
M_F1 = A_1.pos_from(O_1).cross(R_F1)

# Théorème de la résultante statique
R_x = R_01.dot(R0.x) + R_P1.dot(R0.x) + R_Fr1.dot(R0.x) + \
    R_T.dot(R0.x) + R_F1.dot(R0.x)
R_y = R_01.dot(R0.y) + R_P1.dot(R0.y) + R_Fr1.dot(R0.y) + \
    R_T.dot(R0.y) + R_F1.dot(R0.y)
Eq_Rx_1 = sp.Eq(R_x, 0)
Eq_Ry_1 = sp.Eq(R_y, 0)
print('Théorème de la résultante')
print(Eq_Rx_1, Eq_Ry_1)
# display(Eq_Rx_1, Eq_Ry_1)

# Théorème du moment en O_1
M_O_1_z = M_01.dot(R0.z) + M_P1.dot(R0.z) + \
    M_Fr1.dot(R0.z) + M_T.dot(R0.z) + M_F1.dot(R0.z)
Eq_Mz_O_1 = sp.Eq(M_O_1_z, 0)
print('Théorème du moment, en O_1')
print(Eq_Mz_O_1)
