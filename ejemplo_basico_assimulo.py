from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
import numpy as np
import matplotlib.pyplot as plt

# Define el sistema F(t, y, ydot)
def residual(t, y, ydot):
    y1, y2, y3 = y
    dy1, dy2, dy3 = ydot

    res_1 = dy1 + y2 - t
    res_2 = y1 + dy2 - np.sin(t)
    res_3 = dy3 + y3 - y1

    return np.array([res_1, res_2, res_3])

# Condiciones iniciales
y0 = np.array([0.0, 0.0, 1.0])         # y1, y2, y3
yd0 = np.array([0.0, 0.0, -1.0])       # valores estimados de dy/dt

# Crear el problema y el solver
problem = Implicit_Problem(residual, y0, yd0, 0.0)
solver = IDA(problem)

# Configurar tolerancias
solver.atol = [1e-6]*3
solver.rtol = 1e-6

# Simular
tfinal = 10.0
t, y, yd = solver.simulate(tfinal)

# Graficar resultados
plt.plot(t, y[:, 0], label='y1')
plt.plot(t, y[:, 1], label='y2')
plt.plot(t, y[:, 2], label='y3')
plt.xlabel('Time')
plt.ylabel('States')
plt.legend()
plt.grid(True)
plt.show()
