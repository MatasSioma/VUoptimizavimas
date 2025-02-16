#2314004 - stud knygeles nr
import matplotlib.pyplot as pl
import numpy as np

# f(x) = (x^2−a)^2/b−1 -> (x^2)^2/4−1
def tiklsoF(x):
    return pow(pow(x, 2), 2) / 3

x = np.arange(0, 10, 0.1)

pl.plot(x, tiklsoF(x))
pl.show()
