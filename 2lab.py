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
        if abs(np.linalg.norm(gradientas)) < epsilon:
            break
    
    return Xi, tikslo_f(Xi), i, np.array(kelias)

def simplex(X0, alpha=1, beta=0.5, gamma=2, niu=-0.5):
    n = len(X0)
    simplex = [X0]
    kelias = [X0.copy()]
    
    f_values = [tikslo_f(X0)]
    
    for i in range(n):
        X_new = copy.deepcopy(X0)
        X_new[i] += alpha
        simplex.append(X_new)
        f_values.append(tikslo_f(X_new))
    
    iterations = 0
    while True:
        iterations += 1
        sorted_indices = np.argsort(f_values)
        simplex = [simplex[i] for i in sorted_indices]
        f_values = [f_values[i] for i in sorted_indices]
        
        Xl, Xg, Xh = simplex[0], simplex[-2], simplex[-1]
        fl, fg, fh = f_values[0], f_values[-2], f_values[-1]
        
        Xc = np.mean(simplex[:-1], axis=0)
        kelias.append(Xl)
        
        Xn = Xc + gamma * (Xc - Xh)
        fn = tikslo_f(Xn)
        
        if fn < fg:
            Z, fZ = Xn, fn
        elif fn < fh:
            Z, fZ = Xn, fn
        else:
            Xn = Xc + beta * (Xh - Xc)
            fn = tikslo_f(Xn)
            if fn < fh:
                Z, fZ = Xn, fn
            else:
                for i in range(1, len(simplex)):
                    simplex[i] = Xl + niu * (simplex[i] - Xl)
                    f_values[i] = tikslo_f(simplex[i])
                continue
        
        simplex[-1] = Z
        f_values[-1] = fZ
        
        if np.linalg.norm(Xh - Xl) < epsilon:
            break
    
    return simplex[0], f_values[0], iterations, np.array(kelias)


x1, x2 = np.meshgrid(np.arange(0, 0.6, 0.01), np.arange(0, 0.6, 0.01))
plt.contour(x1, x2, tikslo_f([x1, x2]), 30, cmap='RdGy')

Xpr, g = np.array([1, 1]), 0.3
X, Z, i, kelias = gradiento(Xpr, g)
print(f'\nGradientinio nusileidimo metodu, γ: {g}; X0: ({Xpr[0]}, {Xpr[1]})\nAts.: {X} per {i} iteracijas')
plt.plot(np.clip(kelias[:, 0], 0, 0.59), np.clip(kelias[:, 1], 0, 0.59), marker='o', label=f'Gradientinis nusileidimas γ: {g}; X0: ({Xpr[0]}, {Xpr[1]})')

Xpr = np.array([0.1, 0.4])
X, Z, i, kelias = greiciausias(Xpr)
print(f'\nGreiciausio nusileidimo metodu, X0: ({Xpr[0]}, {Xpr[1]})\nAts.: {X} per {i} iteracijas')
plt.plot(np.clip(kelias[:, 0], 0, 0.59), np.clip(kelias[:, 1], 0, 0.59), marker='o', label=f'Greiciausias nusileidimas X0: ({Xpr[0]}, {Xpr[1]})')

Xpr = np.array([0.1, 0.1])
X, Z, i, kelias = simplex(Xpr)
print(f'\nDeformuojamas simpleksas X0: ({Xpr[0]}, {Xpr[1]})\nAts.: {X} per {i} iteracijas')
plt.plot(np.clip(kelias[:, 0], 0, 0.59), np.clip(kelias[:, 1], 0, 0.59), marker='o', label=f'Deformuojamas simpleksas X0: ({Xpr[0]}, {Xpr[1]})')

plt.legend()
plt.colorbar()
plt.show()