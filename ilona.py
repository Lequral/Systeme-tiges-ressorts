import sympy as sp
from sympy.physics.mechanics import dynamicsymbols, Point, ReferenceFrame

#paramètres
r_21 , r_22, r_3, theta_20, m2, g = sp.symbols('r_21 r__22 r_3 l theta_20 m2 g ', real = True, positive = True, nonzero = True)

# paramètres de position
theta_2, gamma_1 = sp.symbols('theta_2 gamma_1', real = True)


# efforts 
F_2, T_2, X02, Y02, N02, P2 = sp.symbols('F_2 T_2 X02 Y02 N02')

# raideur du ressort
k2 = sp.symbols('k2', real = True, positive = True, nonzero = True)

#repères
R0 = ReferenceFrame('R_0') 
R2 = ReferenceFrame('R_2')
R3 = ReferenceFrame('R_3')

#position des repères
R2.orient(R0,'Axis',[theta_2, R0.z])
R3.orient(R0,'Axis',[gamma_1, R0.z])

# Points caractéristiques 

A2 = Point('A2')
O2 = Point('O2')
B2 = Point('B2')
C = Point('C')
D = Point('D')

# position associées 
O2.set_vel(R0, 0) 
O2.set_vel(R2, 0) 
O2.set_vel(R3, 0) 

A2.set_pos(O2, -r_21 * R2.x )
A2.set_vel(R2, 0.0) # A est un point de R2 donc vitesse nulle par rapport à 2

#
B2.set_pos(O2, r_22 * R2.x ) 
B2.set_vel(R2, 0.0) 

#
C.set_pos(O2, r_22 * R2.x + l * R0.y)
C.set_vel(R3, 0.0) 

#
D.set_pos(O2, r_22 * R2.x + l * R0.y + r_3 * R3.x)
D.set_vel(R3, 0.0) 
D.set_vel(R0, 0.0) 

#BAME sur la barre 1

# Liaison 0/2
R_02 = X02 * R0.x + Y02 * R0.y # Pivot et analyse plane
M_02 = N02 * R0.z

# actions extérieure 
R_P2 = -m2*g * R0.y # en G2
F_2 = F_2 * R2.y
T_2 = T_2 * R2.y
MO2 = k2 *(theta_2 - theta_20)

#PFS 

R_x = R_02.dot(R0.x) + R_P2.dot(R0.x) + T_2.dot(R0.x) + F_2.dot(R0.x)
R_y = R_02.dot(R0.y) + R_P2.dot(R0.y) + T_2.dot(R0.y) + F_2.dot(R0.y)

O2A2 = A2.pos_from(O2)
M_O2 = MO2.dot(R0.z) + O2A2.cross(T_2)
Eq_Rx_1 = sp.Eq(R_x, 0)
Eq_Ry_1 = sp.Eq(R_y, 0)
Eq_M_O2 = sp.Eq(M_O2, 0)

print(Eq_Rx_1, Eq_Ry_1,Eq_M_O2)