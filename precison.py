# D'aprèsles sources suivantes, pour un modèle de réseau de neuronnes, la précision n'est pas très importante et on peut s'en sortir avec une précision de 10^-4
# https://docs.nvidia.com/deeplearning/performance/mixed-precision-training/index.html
# https://massedcompute.com/faq-answers/?question=Can%20you%20explain%20the%20impact%20of%20float16%20and%20float32%20on%20deep%20learning%20model%20accuracy?

import numpy as np
from scipy.optimize import fsolve

d=3

def f(x):
    return np.sin(x) - x + 10**(-d)

def g(x):
    return np.cos(x) - x + 10**(-d)

angle_max = min(fsolve(f, .1)[0],fsolve(g, .1)[0])

print(f"L'angle maximum pour conserver l'hypothèse des petits angles à une précision près de 10^{-d} est de {angle_max*180/np.pi} °")