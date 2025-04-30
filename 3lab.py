import numpy as np
from simplex import simplex

def tikslo_f(x):
    a, b, c = x
    return -a * b * c

def g(x):
    a, b, c = x
    return 2 * a * b + 2 * a * c + 2 * b * c - 1

def h(x):
    return sum(max(0, -xi) ** 2 for xi in x)

def baudos_f(x, r):
    fx = tikslo_f(x)
    gx = g(x)
    hx = h(x)
    return fx + (1 / r) * (gx ** 2 + hx)

def main():
    # 2314004 - studid
    Xpr = [
        [0.0, 0.0, 0.0],  # X0
        [1.0, 1.0, 1.0],  # X1
        [0.0, 0.0, 0.4],  # Xm (a=0, b=0, c=4)
    ]

    alpha = 0.5
    epsilon = 1e-4

    for start in Xpr:
        X0 = start.copy()
        bauda_r = 20.0
        i = 0
        f_kvietimai = 0

        def bf(x, r):
            nonlocal f_kvietimai
            f_kvietimai += 1
            return baudos_f(x, r)

        while bauda_r > epsilon:
            X0, i_simplex = simplex(
                X0, bf, bauda_r, epsilon, alpha
            )
            i += i_simplex
            bauda_r /= 10.0

        print(f"Pradinis taškas: {start}")
        print(f"Rastas minimumo taškas: {X0}")
        print(f"Iteracijų kiekis: {i}")
        print(f"Tikslo funkcijos skaičiavimų kiekis: {f_kvietimai}")
        print(f"Tikslo funkcijos reikšmė minimumo taške: {baudos_f(X0, bauda_r)}\n")


if __name__ == "__main__":
    main()
