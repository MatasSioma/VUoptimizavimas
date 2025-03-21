import matplotlib.pyplot as plt
import numpy as np
import copy

epsilon = pow(10, -4)

def tikslo_f(X):
    a, b = np.clip(X, -100, 100)
    return -1 / 8 * (a * b - a**2 * b - a * b**2)
 
def tikslo_fg(X):
    a, b = np.clip(X, -100, 100)
    dfda = -1 / 8 * (b - 2 * a * b - b**2)
    dfdb = -1 / 8 * (a - a**2 - 2 * a * b)
    return np.array([dfda, dfdb])

def gradiento(X0, gamma):
    Xi = np.array(copy.deepcopy(X0), dtype=float)
    i = 0
    kelias = [Xi.copy()]
    gradientas = tikslo_fg(Xi)
    while True:
        Xi = Xi - gamma * gradientas
        i += 1
        kelias.append(Xi)
        gradientas = tikslo_fg(Xi)
        if abs(np.linalg.norm(gradientas)) < epsilon:
            break
    
    return Xi, tikslo_f(Xi), i, np.array(kelias)


def min_gamma(Xi, gradientas):
    t_f = lambda gamma: tikslo_f(Xi - gamma * gradientas)
    l = 0
    r = 25
    xm = (l + r) / 2
    ilgis = r - l
    fxm = t_f(xm)
    
    iteracijos = 0
    while True:
        iteracijos += 1
        
        x1 = l + ilgis / 4
        x2 = r - ilgis / 4
        fx1 = t_f(x1)
        fx2 = t_f(x2)
        
        if fx1 < fxm:
            r = xm
            xm = x1
            fxm = fx1
        elif fx2 < fxm:
            l = xm
            xm = x2
            fxm = fx2
        else:
            l = x1
            r = x2
        
        ilgis = r - l
        if ilgis < epsilon:
            break
    
    return xm #gamma kur f min

def greiciausias(X0):
    Xi = np.array(copy.deepcopy(X0), dtype=float)
    i = 0
    kelias = [Xi.copy()]
    gradientas = tikslo_fg(Xi)

    while True:
        gamma = min_gamma(Xi, gradientas)
        while True:
            Xi = Xi - gamma * gradientas
            if Xi[0] < 0 and Xi[1] < 0:
                Xi = Xi + gamma * gradientas
                gamma /= 2
                continue
            break
        i += 1
        kelias.append(Xi)
        gradientas = tikslo_fg(Xi)
        print(gamma, Xi)
        if abs(np.linalg.norm(gradientas)) < epsilon:
            break
    
    return Xi, tikslo_f(Xi), i, np.array(kelias)

x1, x2 = np.meshgrid(np.arange(0, 0.6, 0.01), np.arange(0, 0.6, 0.01))
plt.contour(x1, x2, tikslo_f([x1, x2]), 30, cmap='RdGy')

Xpr = np.array([1, 1])
X, Z, i, kelias = gradiento(Xpr, 0.3)
plt.plot(np.clip(kelias[:, 0], 0, 0.59), np.clip(kelias[:, 1], 0, 0.59), marker='o', label=f'Gradientinis nusileidimas γ: 0.3; X0: ({Xpr[0]}, {Xpr[1]})')

Xpr = np.array([0.1, 0.4])
X, Z, i, kelias = greiciausias(Xpr)
plt.plot(np.clip(kelias[:, 0], 0, 0.59), np.clip(kelias[:, 1], 0, 0.59), marker='o', label=f'Greiciausias nusileidimas X0: ({Xpr[0]}, {Xpr[1]})')

plt.legend()
plt.colorbar()
plt.show()