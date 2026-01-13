import sympy as sp
from sympy.physics.mechanics import dynamicsymbols, Point, ReferenceFrame
from scipy.optimize import least_squares


def equations(est_simplifie: bool, debug=False) -> tuple:
    """
    Trouve les équations du problème avec ou sans hypothèse. Sans hypothèse, on ne peux obtenir une expression analytique de gamma1, etc...
    """
    equations = []

    # géométrie
    r_11, r_12, r_21, r_22 = sp.symbols(
        'r_11 r_12 r_21 r_22', real=True, positive=True, nonzero=True)
    m_1, m_2, m_3, g = sp.symbols(
        'm_1 m_2 m_3 g', real=True, positive=True, nonzero=True)
    a, b, r_3, L = sp.symbols('a b r_3 L', real=True,
                              positive=True, nonzero=True)
    W01, W02, B1, theta10, theta20, gamma10 = sp.symbols('w_01 w_02 b_1 theta_10 theta_20 gamma_10', real=True,
                                                         positive=False, nonzero=True)

    # paramètres de position
    theta1, theta2, gamma1 = sp.symbols('theta_1 theta_2 gamma_1', real=True)

    # efforts
    X01, Y01, F_1, T_1 = sp.symbols('X_01 Y_01 F_1 T_1')
    X02, Y02, F_2, T_2 = sp.symbols("X_02 Y_02 F_2 T_2")
    T_3, Tp_3 = sp.symbols("T_3 T'_3")
    X03, Y03 = sp.symbols("X_03 Y_03")

    # réologie : facteur d'amortisement et raideur du ressort
    k_1, k_2, k_3 = sp.symbols('k_1 k_2 k_3', real=True,
                               positive=True, nonzero=True)

    # Repères attachées au pièces en utilisant des RéférenceFrame plus adaptés à la mécanique des systèmes qu'un CoordSys3D
    R0 = ReferenceFrame('R_0')
    R1 = ReferenceFrame('R_1')
    R2 = ReferenceFrame('R_2')
    R3 = ReferenceFrame('R_3')

    # positionnement des repères
    # R1 se déduit de R0 par rotation d'angle theta1 autour de z0=z1
    R1.orient(R0, 'Axis', [theta1, R0.z])
    R2.orient(R0, 'Axis', [theta2, R0.z])
    R3.orient(R0, 'Axis', [gamma1, R0.z])

    # Points caractéristiques
    O_1 = Point('O_1')
    A_1 = Point('A_1')
    B_1 = Point('B_1')
    G_1 = Point('G_1')
    O_2 = Point('O_2')
    A_2 = Point('A_2')
    B_2 = Point('B_2')
    G_2 = Point('G_2')
    D = Point('D')
    C = Point('C')

    # Et position associées
    O_1.set_vel(R0, 0)  # O est un point de R0
    O_1.set_vel(R1, 0)

    A_1.set_pos(O_1, -r_11 * R1.x)
    A_1.set_vel(R1, 0)

    B_1.set_pos(O_1, r_12 * R1.x)
    B_1.set_vel(R1, 0)

    G_1.set_vel(R1, 0)
    G_1.set_pos(O_1, 0*R1.x)

    O_2.set_vel(R0, 0)
    O_2.set_vel(R2, 0)

    A_2.set_pos(O_2, -r_21 * R2.x)
    A_2.set_vel(R2, 0)

    B_2.set_pos(O_2, r_22 * R2.x)
    B_2.set_vel(R2, 0)

    G_2.set_vel(R2, 0)
    G_2.set_pos(O_2, 0*R2.x)

    D.set_vel(R0, 0)  # O est un point de R0
    D.set_vel(R3, 0)  # O est un point de R02
    D.set_pos(O_1, a*R0.x-b*R0.y)
    D.set_pos(O_2, a*R0.x+b*R0.y)

    C.set_vel(R3, 0)
    C.set_pos(D, -r_3*R3.x)

    O_1D = D.pos_from(O_1)
    O_2D = D.pos_from(O_2)
    DC = C.pos_from(D)
    O_1B_1 = B_1.pos_from(O_1)
    O_2B_2 = B_2.pos_from(O_2)

    # PFS ∑ = Barre1
    # Torseurs
    R_01 = X01 * R0.x + Y01 * R0.y
    M_01 = 0 * R0.z

    R_P1 = -m_1*g*R0.y
    M_P1 = G_1.pos_from(O_1).cross(R_P1)

    R_Fr1 = 0*R0.x
    M_Fr1 = k_1*(theta1-theta10)*R0.z

    R_T = T_1*R0.y
    if est_simplifie:
        R_T = T_1*R1.y  # hypothèse la tension est perpendiculaire à la tige
    M_T = B_1.pos_from(O_1).cross(R_T)

    R_F1 = F_1*R0.y
    if est_simplifie:
        R_F1 = F_1*R1.y  # hypothèse la tension est perpendiculaire à la tige
    M_F1 = A_1.pos_from(O_1).cross(R_F1)

    # Théorème de la résultante statique
    R_x = R_01.dot(R0.x) + R_P1.dot(R0.x) + R_Fr1.dot(R0.x) + \
        R_T.dot(R0.x) + R_F1.dot(R0.x)
    R_y = R_01.dot(R0.y) + R_P1.dot(R0.y) + R_Fr1.dot(R0.y) + \
        R_T.dot(R0.y) + R_F1.dot(R0.y)
    Eq_Rx_1 = sp.Eq(R_x, 0)
    Eq_Ry_1 = sp.Eq(R_y, 0)

    # Théorème du moment en O_1
    M_O_1_z = M_01.dot(R0.z) + M_P1.dot(R0.z) + \
        M_Fr1.dot(R0.z) + M_T.dot(R0.z) + M_F1.dot(R0.z)
    Eq_Mz_O_1 = sp.Eq(M_O_1_z, 0)

    # PFS ∑ = Barre2
    # Torseurs
    R_02 = X02 * R0.x + Y02 * R0.y
    M_02 = 0 * R0.z

    R_P2 = -m_2*g*R0.y
    M_P2 = G_2.pos_from(O_2).cross(R_P2)

    R_Fr2 = 0*R0.x
    M_Fr2 = k_2*(theta2-theta20)*R0.z

    R_Tp = T_2*R0.y
    if est_simplifie:
        R_Tp = T_2*R2.y  # hypothèse la tension est perpendiculaire à la tige
    M_Tp = B_2.pos_from(O_2).cross(R_Tp)

    R_F2 = F_2*R0.y
    if est_simplifie:
        R_F2 = F_2*R2.y  # hypothèse la tension est perpendiculaire à la tige
    M_F2 = A_2.pos_from(O_2).cross(R_F2)

    # Théorème de la résultante statique
    R_xp = R_02.dot(R0.x) + R_P2.dot(R0.x) + R_Fr2.dot(R0.x) + \
        R_Tp.dot(R0.x) + R_F2.dot(R0.x)
    R_yp = R_02.dot(R0.y) + R_P2.dot(R0.y) + R_Fr2.dot(R0.y) + \
        R_Tp.dot(R0.y) + R_F2.dot(R0.y)
    Eq_Rx_2 = sp.Eq(R_xp, 0)
    Eq_Ry_2 = sp.Eq(R_yp, 0)

    # Théorème du moment en O_2
    M_O_2_z = M_02.dot(R0.z) + M_P2.dot(R0.z) + \
        M_Fr2.dot(R0.z) + M_Tp.dot(R0.z) + M_F2.dot(R0.z)
    Eq_Mz_O_2 = sp.Eq(M_O_2_z, 0)

    # PFS ∑=fil1 puis ∑=fil2
    # Liaison bielle
    terme_gauche_bielle_1x = O_1D.dot(R0.x)+DC.dot(R0.x)-(O_1B_1.dot(R0.x))
    terme_gauche_bielle_1y = O_1D.dot(R0.y)+DC.dot(R0.y)-(O_1B_1.dot(R0.y))
    terme_gauche_bielle_2x = O_2D.dot(R0.x)+DC.dot(R0.x)-(O_2B_2.dot(R0.x))
    terme_gauche_bielle_2y = O_2D.dot(R0.y)+DC.dot(R0.y)-(O_2B_2.dot(R0.y))

    if est_simplifie:
        # Linéarisation
        terme_gauche_bielle_1y = sp.series(sp.series(
            terme_gauche_bielle_1y, theta1, 0, 3).removeO(), gamma1, 0, 3).removeO()
        terme_gauche_bielle_2y = sp.series(sp.series(
            terme_gauche_bielle_2y, theta2, 0, 3).removeO(), gamma1, 0, 3).removeO()

    Eq_bielle_1x = sp.Eq(terme_gauche_bielle_1x, L)
    Eq_bielle_1y = sp.Eq(terme_gauche_bielle_1y, 0)
    Eq_bielle_2x = sp.Eq(terme_gauche_bielle_2x, L)
    Eq_bielle_2y = sp.Eq(terme_gauche_bielle_2y, 0)

    # PFS ∑=barre3
    # Liaison 0/3
    R_03 = X03 * R0.x + Y03 * R0.y  # Pivot et analyse plane

    # actions extérieure
    R_P3 = -m_3*g * R0.y  # en G3 = D

    T_3 = T_1 * R3.y
    M_T_3 = C.pos_from(D).cross(T_3)

    Tp_3 = T_2 * R3.y
    M_Tp_3 = C.pos_from(D).cross(Tp_3)

    MO3 = k_3 * (gamma1 - gamma10)*R0.z

    R_x_3 = R_03.dot(R0.x) + R_P3.dot(R0.x) + T_3.dot(R0.x) + Tp_3.dot(R0.x)
    R_y_3 = R_03.dot(R0.y) + R_P3.dot(R0.y) + T_3.dot(R0.y)+Tp_3.dot(R0.y)
    Eq_Rx_3 = sp.Eq(R_x_3, 0)
    Eq_Ry_3 = sp.Eq(R_y_3, 0)

    M3 = MO3.dot(R0.z) + M_T_3.dot(R0.z) + M_Tp_3.dot(R0.z)

    Eq_M_O3 = sp.Eq(M3, 0)

    if not est_simplifie:
        if debug:
            try:
                display(Eq_Rx_1)
                display(Eq_Ry_1)
                display(Eq_Mz_O_1)
                display(Eq_Rx_2)
                display(Eq_Ry_2)
                display(Eq_Mz_O_2)
                display(Eq_bielle_1x)
                display(Eq_bielle_1y)
                display(Eq_bielle_2x)
                display(Eq_bielle_2y)
                display(Eq_Rx_3)
                display(Eq_Ry_3)
                display(Eq_M_O3)
            except:
                pass
        equations.append(Eq_Rx_1)
        equations.append(Eq_Ry_1)
        equations.append(Eq_Mz_O_1)
        equations.append(Eq_Rx_2)
        equations.append(Eq_Ry_2)
        equations.append(Eq_Mz_O_2)
        equations.append(Eq_bielle_1x)
        equations.append(Eq_bielle_1y)
        equations.append(Eq_bielle_2x)
        equations.append(Eq_bielle_2y)
        equations.append(Eq_Rx_3)
        equations.append(Eq_Ry_3)
        equations.append(Eq_M_O3)

        param = (r_11, r_12, r_21, r_22, m_1, m_2, m_3, g, a, b, r_3, L, theta10, theta20, gamma10,
                 theta1, theta2, gamma1, X01, Y01, F_1, T_1, X02, Y02, F_2, T_2, X03, Y03, k_1, k_2, k_3)
        func = sp.lambdify((
            param
        ), equations, 'numpy')

        def numpy_equations(vars):
            r_11, r_12, r_21, r_22, m_1, m_2, m_3, g, a, b, r_3, L, theta10, theta20, gamma10, theta1, theta2, gamma1, X01, Y01, F_1, T_1, X02, Y02, F_2, T_2, X03, Y03, k_1, k_2, k_3 = vars
            return func(r_11, r_12, r_21, r_22, m_1, m_2, m_3, g, a, b, r_3, L, theta10, theta20, gamma10, theta1, theta2, gamma1, X01, Y01, F_1, T_1, X02, Y02, F_2, T_2, X03, Y03, k_1, k_2, k_3)

        return {
            "numpy_equation": numpy_equations, "param": param, "equation": equations
        }
    else:
        # Variable à réutiliser
        expression_theta1 = sp.solve(Eq_bielle_1y, theta1)[0]
        expression_theta2 = sp.solve(Eq_bielle_2y, theta2)[0]

        Eq_Mz_O_1 = Eq_Mz_O_1.subs({theta1: expression_theta1})
        Eq_Mz_O_2 = Eq_Mz_O_2.subs({theta2: expression_theta2})

        expression_T_1 = sp.solve(Eq_Mz_O_1, T_1)[0]
        expression_T_2 = sp.solve(Eq_Mz_O_2, T_2)[0]

        Eq_M_O3 = Eq_M_O3.subs({T_1: expression_T_1, T_2: expression_T_2})

        # Equation finale
        expression_gamma1 = sp.solve(Eq_M_O3, gamma1)[0]

        expr_gamma1_simple = W01*F_1+W02*F_2+B1

        W01_value = -(r_11*r_12*r_22**2*r_3)/(k_1*r_22**2*r_3 **
                                              2 + k_2*r_12**2*r_3**2 - k_3*r_12**2*r_22**2)
        W02_value = -(r_12**2*r_21*r_22*r_3)/(k_1*r_22**2*r_3 **
                                              2 + k_2*r_12**2*r_3**2 - k_3*r_12**2*r_22**2)
        B1_value = (- b*k_1*r_22**2*r_3 + b*k_2*r_12**2*r_3 - gamma10*k_3*r_12**2*r_22**2 - k_1*r_12*r_22**2*r_3 *
                    theta10 - k_2*r_12**2*r_22*r_3*theta20)/(k_1*r_22**2*r_3**2 + k_2*r_12**2*r_3**2 - k_3*r_12**2*r_22**2)

        assert sp.simplify(expression_gamma1-expr_gamma1_simple.subs({
            W01: W01_value,
            W02: W02_value,
            B1: B1_value
        })) == 0, "Verification de l'égalité (entre l'équation simplifiée et la développée) :"

        expression_gamma1 = sp.simplify(
            expression_gamma1.subs({theta10: 0, theta20: 0, gamma10: 0}))
        B1_value = sp.simplify(B1_value.subs(
            {theta10: 0, theta20: 0, gamma10: 0}))

        if debug:
            try:
                # Si le fichier est executé dans un Jupyter Notebook
                display("Equations :")
                print('Théorème de la résultante ∑ = Barre1')
                display(Eq_Rx_1, Eq_Ry_1, Eq_Mz_O_1)
                print('Théorème de la résultante ∑ = Barre2')
                display(Eq_Rx_2, Eq_Ry_2, Eq_Mz_O_2)
                print("Equations bielle linéarisées")
                display(Eq_bielle_1x, Eq_bielle_1y, Eq_bielle_2x, Eq_bielle_2y)
                print('Théorème de la résultante à S3 isolé')
                display(Eq_Rx_3, Eq_Ry_3, Eq_M_O3)
                print("Variable à réutiliser grâce aux équations bielles")
                display(sp.Eq(theta1, expression_theta1))
                display(sp.Eq(theta2, expression_theta2))
                print("substitution dans les théorèmes des moments en O_1 et O_2")
                display(Eq_Mz_O_1, Eq_Mz_O_2)
                print("Expression de T et T'")
                display(sp.Eq(T_2, expression_T_2))
                display(sp.Eq(T_1, expression_T_1))
                print("substitution dans le théorème des moments en D")
                display(Eq_M_O3)
                display(sp.Eq(gamma1, expression_gamma1))
                display(sp.Eq(gamma1, expr_gamma1))

                display(sp.Eq(gamma1, expression_gamma1))
                display(sp.Eq(W01, W01_value))
                display(sp.Eq(W02, W02_value))
                display(sp.Eq(B1, B1_value))
            except:
                pass

        return {
            "weight_expression": [W01_value, W02_value, B1_value],
            "param": (r_11, r_12, r_21, r_22, r_3, k_1, k_2, k_3, b)
        }


def trouver_param_real(poids_reseaux: list, debug = False) -> list:
    """
    Trouve les valeurs du problèmes réels pour correspondre aux biais donné (poids_reseaux) en respectant les équations donnés (equationSimple)
    """
    dico = equations(True, debug)
    weight_expr, param = dico["weight_expression"], dico["param"]
    r_11, r_12, r_21, r_22, r_3, k_1, k_2, k_3, b =  param

    eq = [weight_expr[i] - poids_reseaux[i] for i in range(3)]

    func = sp.lambdify(param, eq, 'numpy')

    def eq(vars):
        r_11, r_12, r_21, r_22, r_3, k_1, k_2, k_3, b = vars
        return func(r_11, r_12, r_21, r_22, r_3, k_1, k_2, k_3, b)

    guess = [1 for i in range(9)]

    min_var = [-100, -100, -100, -100, -100, 1, 1, 1, -100]
    max_var = [100, 100, 100, 100, 100, 100, 100, 100, 100]

    solu = least_squares(eq, guess, bounds=(
        min_var, max_var), max_nfev=10000)

    return solu # A tester