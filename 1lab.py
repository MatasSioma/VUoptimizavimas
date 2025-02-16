#2314004 - stud knygeles nr
import matplotlib.pyplot as pl
import numpy as np
import copy, math

# f(x) = (x^2−a)^2/b−1 -> (x^2)^2/4−1
def tiklsoF(x):
    # return pow(pow(x, 2), 2) / 3
    return pow(pow(x, 2) - 8, 2) / 7 #pvz su kitais a, b

while True:
    try:
        intervalas = (int(input("pradinis rėžis (l): ")), int(input("galinis rėžis (r): ")))
        epsilon = pow(10, int(input("ε = 10^")))
        if intervalas[0] > intervalas[1]:
            print("Klaida! 'l' yra mažiau negu 'r'.")
            continue
        break
    except ValueError:
        print("Įvestis turi būti skaičius.")

def dalijimoPusiau():
    global intervalas, epsilon
    l, r = copy.deepcopy(intervalas[0]), copy.deepcopy(intervalas[1])
    xm = (l + r) / 2
    ilgis = r - l

    iteracijos = 0
    while True:
        iteracijos += 1

        fxm = tiklsoF(xm)
        x1 = l + ilgis / 4
        x2 = r - ilgis / 4

        if tiklsoF(x1) < fxm:
            r = xm
            xm = x1
        elif tiklsoF(x2) < fxm:
            l = xm
            xm = x2
        else:
            l = x1
            r = x2

        ilgis = r - l
        if ilgis < epsilon:
            break

    print(f"Dalijimo pusiau metodo alogritmas rado atsakymą per {iteracijos} iteracijas.")
    print(f"su ε = {epsilon}, ats.: {tiklsoF(xm)}")
    pl.plot(xm, tiklsoF(xm), 'or')

def auksinioPjuvio():
    global intervalas, epsilon
    l, r = copy.deepcopy(intervalas[0]), copy.deepcopy(intervalas[1])

    phi = (math.sqrt(5) - 1) / 2

    ilgis = r - l
    x1 = r - phi*ilgis
    x2 = l + phi*ilgis

    fx1 = tiklsoF(x1)
    fx2 = tiklsoF(x2)

    iteracijos = 0
    while True:
        iteracijos += 1
        if fx1 > fx2:
            l = x1
            ilgis = r - l

            x1 = x2
            fx1 = fx2
            x2 = l + phi * ilgis
            fx2 = tiklsoF(x2)
        else:
            r = x2
            ilgis = r - l

            x2 = x1
            fx2 = fx1
            x1 = r - phi * ilgis
            fx1 = tiklsoF(x1)

        if ilgis < epsilon:
            break

    if fx1 < fx2: ats = (x1, fx1)
    else: ats = (x2, fx2)

    print(f"Auksinio pjūvio metodo alogritmas rado atsakymą per {iteracijos} iteracijas.")
    print(f"su ε = {epsilon}, ats.: {ats[1]}")
    pl.plot(ats[0], ats[1], 'og')


x = np.arange(intervalas[0], intervalas[1], 0.1)
pl.plot(x, tiklsoF(x))

dalijimoPusiau()
auksinioPjuvio()

pl.show()
