import numpy as np


def simplex(X0, tikslo_f, r, epsilon, alpha=0.1, beta=0.5, gamma=2.0, niu=-0.5):
    n = len(X0)
    simplex = [X0]
    kelias = [X0.copy()]

    delta1 = (np.sqrt(n + 1) + n - 1) / (n * np.sqrt(2)) * alpha
    delta2 = (np.sqrt(n + 1) - 1) / (n * np.sqrt(2)) * alpha

    f_kvietimai = 1
    def f(x):
        nonlocal f_kvietimai
        f_kvietimai += 1
        return tikslo_f(x, r)
    
    f_values = [f(X0)]
    
    for i in range(n):
        X_new = np.array(X0, dtype=np.float64)
        for j in range(n):
            if j == i:
                X_new[j] += delta2
            else:
                X_new[j] += delta1
        simplex.append(X_new)
        f_values.append(f(X_new))
    
    iterations = 0
    while True:
        iterations += 1
        indeksai = np.argsort(f_values)
        simplex = [simplex[i] for i in indeksai]
        f_values = [f_values[i] for i in indeksai]

        Xl, Xg, Xh = simplex[0], simplex[-2], simplex[-1]
        fl, fg, fh = f_values[0], f_values[-2], f_values[-1]
        
        Xc = np.mean(simplex[:-1], axis=0)
        
        Xn = Xc + gamma * (Xc - Xh)
        fn = f(Xn)
        
        if fn < fl:
            Z, fZ = Xn, fn
        elif fn < fg:
            Z, fZ = Xn, fn
        else:
            Xn = Xc + beta * (Xh - Xc)
            fn = f(Xn)
            if fn < fh:
                Z, fZ = Xn, fn
            else:
                for i in range(1, len(simplex)):
                    simplex[i] = Xl + niu * (simplex[i] - Xl)
                    f_values[i] = f(simplex[i])
                continue

        simplex[-1] = Z
        f_values[-1] = fZ
        
        kelias.append(Z)

        if np.linalg.norm(Xh - Xl) < epsilon:
            break
    
    return Xl, iterations