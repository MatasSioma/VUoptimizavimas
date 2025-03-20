import matplotlib.pyplot as plt
import numpy as np
import copy

while True:
    try:
        print("Pradzios taskas:")
        Xpr = np.array([int(input("X1: ")), int(input("X2: "))])
        epsilon = pow(10, int(input("ε = 10^")))
        break
    except:
        print("Įvestis turi būti skaičius.")

def tikslo_f(X):
    a, b = np.clip(X, -100, 100)
    return -1 / 8 * (a * b - a**2 * b - a * b**2)
 
def tikslo_fg(X):
    a, b = np.clip(X, -100, 100)
    dfda = -1 / 8 * (b - 2 * a * b - b**2)
    dfdb = -1 / 8 * (a - a**2 - 2 * a * b)
    return np.array([dfda, dfdb])

def gradiento(X0, epsilon, gamma = 0.1):
    Xi = np.array(copy.deepcopy(X0), dtype=float)
    i = 0
    kelias = [Xi.copy()]
    while True:
        gradientas = tikslo_fg(Xi)
        Xi = Xi - gamma * gradientas
        print(Xi)
        kelias.append(Xi)
        i += 1
        if abs(np.linalg.norm(tikslo_fg(Xi))) < epsilon:
            break
    
    return Xi, tikslo_f(Xi), i, np.array(kelias)

x1, x2 = np.meshgrid(np.arange(0, 0.6, 0.01), np.arange(0, 0.6, 0.01))
plt.contour(x1, x2, tikslo_f([x1, x2]), 30, cmap='RdGy')

X, Z, i, kelias = gradiento(Xpr, epsilon)
plt.plot(kelias[:, 0], kelias[:, 1])

plt.colorbar()
plt.show()