# D'aprèsles sources suivantes, pour un modèle de réseau de neuronnes, la précision n'est pas très importante et on peut s'en sortir avec une précision de 10^-4
# https://docs.nvidia.com/deeplearning/performance/mixed-precision-training/index.html
# https://massedcompute.com/faq-answers/?question=Can%20you%20explain%20the%20impact%20of%20float16%20and%20float32%20on%20deep%20learning%20model%20accuracy?

import numpy as np
from scipy.optimize import fsolve

d = 4


def f(x):
    return np.abs(np.sin(x) - x) - 5*10**(-d)


# def g(x):
#     return np.abs(np.cos(x) - (1 - x**2/2)) - 5*10**(-d)


angle_max = fsolve(f, .1)[0]*180/np.pi

print(
    f"L'angle maximum pour conserver l'hypothèse des petits angles à une précision près de 10^{-d} est de {angle_max} °")
