#2314004 - stud knygeles nr
import matplotlib.pyplot as pl
import numpy as np

# f(x) = (x^2−a)^2/b−1 -> (x^2)^2/4−1
def tiklsoF(x):
    return pow(pow(x, 2), 2) / 3

while True:
    try:
        intervalas = (int(input("pradinis rėžis (l): ")), int(input("galinis rėžis (r): ")))
        l, r = intervalas[0], intervalas[1]
        if l > r:
            print("Klaida! 'l' yra mažiau negu 'r'.")
            continue
        break
    except ValueError:
        print("Įvestis turi būti skaičius.")

xm = (l + r) / 2
ilgis = r - l

fxm = tiklsoF(xm)
x1 = l + ilgis / 4
x2 = r - ilgis / 4

if tiklsoF(x1) < fxm:
    pass
elif tiklsoF(x2) < fxm:
    pass
else:
    pass

x = np.arange(0, 10, 0.1)
pl.plot(x, tiklsoF(x))
pl.show()
