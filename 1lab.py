#2314004 - stud knygeles nr
import matplotlib.pyplot as pl
import numpy as np

# f(x) = (x^2−a)^2/b−1 -> (x^2)^2/4−1
def tiklsoF(x):
    return pow(pow(x, 2), 2) / 3
    # return pow(pow(x, 2) - 8, 2) / 7 #pvz su kitais a, b


while True:
    try:
        intervalas = (int(input("pradinis rėžis (l): ")), int(input("galinis rėžis (r): ")))
        epsilon = pow(10, int(input("ε = 10^")))
        l, r = intervalas[0], intervalas[1]
        if l > r:
            print("Klaida! 'l' yra mažiau negu 'r'.")
            continue
        break
    except ValueError:
        print("Įvestis turi būti skaičius.")

xm = (l + r) / 2
ilgis = r - l

interacijos = 0
while True:
    interacijos += 1

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

x = np.arange(0, 10, 0.1)
pl.plot(x, tiklsoF(x))
pl.plot(xm, tiklsoF(xm), 'or')
print(interacijos)
pl.show()
