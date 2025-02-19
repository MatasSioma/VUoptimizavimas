#2314004 - stud knygeles nr
import matplotlib.pyplot as pl
import numpy as np
import copy, math

# f(x) = (x^2−a)^2/b−1 -> (x^2)^2/4−1
def tiklsoF(x):
    return (x**4) / 4 - 1
    # return pow(pow(x, 2) - 8, 2) / 7 #pvz su kitais a, b

def tiklsoF_dx(x):
    return x**3

def tiklsoF_2dx(x):
    return 3 * x**2

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
    f_kvietimai = 0

    def f(x):
        nonlocal f_kvietimai
        f_kvietimai += 1
        return tiklsoF(x)
    
    l, r = copy.deepcopy(intervalas[0]), copy.deepcopy(intervalas[1])
    xm = (l + r) / 2
    ilgis = r - l
    fxm = f(xm)
    
    iteracijos = 0
    tarpiniai_t = []
    while True:
        iteracijos += 1
        
        x1 = l + ilgis / 4
        x2 = r - ilgis / 4
        fx1 = f(x1)
        fx2 = f(x2)
        
        tarpiniai_t.append((xm, fxm))
        
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

    print()
    print(f"Dalijimo pusiau metodo algoritmas rado atsakymą per {iteracijos} iteracijas.")
    print(f"su ε = {epsilon}, ats.: {fxm}")
    print(f"tiklsoF() buvo iškviestas {f_kvietimai} kartų.")
    
    for pt in tarpiniai_t:
        pl.plot(pt[0], pt[1], 'xc')

    pl.plot(xm, fxm, 'or')

def auksinioPjuvio():
    global intervalas, epsilon
    f_kvietimai = 0

    def f(x):
        nonlocal f_kvietimai
        f_kvietimai += 1
        return tiklsoF(x)
    
    l, r = copy.deepcopy(intervalas[0]), copy.deepcopy(intervalas[1])
    phi = (math.sqrt(5) - 1) / 2

    ilgis = r - l
    x1 = r - phi*ilgis
    x2 = l + phi*ilgis
    fx1 = f(x1)
    fx2 = f(x2)
    
    iteracijos = 0
    tarpiniai_t = []
    while True:
        iteracijos += 1
        
        tarpiniai_t.append((x1, fx1))
        tarpiniai_t.append((x2, fx2))
        
        if fx1 > fx2:
            l = x1
            ilgis = r - l
            x1 = x2
            fx1 = fx2
            x2 = l + phi * ilgis
            fx2 = f(x2)
        else:
            r = x2
            ilgis = r - l
            x2 = x1
            fx2 = fx1
            x1 = r - phi * ilgis
            fx1 = f(x1)
            
        if ilgis < epsilon:
            break
    
    if fx1 < fx2:
        ats = (x1, fx1)
    else:
        ats = (x2, fx2)
    
    print()
    print(f"Auksinio pjūvio metodo algoritmas rado atsakymą per {iteracijos} iteracijas.")
    print(f"su ε = {epsilon}, ats.: {ats[1]}")
    print(f"tiklsoF() buvo iškviestas {f_kvietimai} kartų.")
    
    for pt in tarpiniai_t:
        pl.plot(pt[0], pt[1], 'xm')
    pl.plot(ats[0], ats[1], 'og')

def niutonoMetodas():
    global intervalas, epsilon
    f_kvietimai = 0
    fdx_kvietimai = 0
    f2dx_kvietimai = 0
    
    def f(x):
        nonlocal f_kvietimai
        f_kvietimai += 1
        return tiklsoF(x)
    
    def f_dx(x):
        nonlocal fdx_kvietimai
        fdx_kvietimai += 1
        return tiklsoF_dx(x)
    
    def f_2dx(x):
        nonlocal f2dx_kvietimai
        f2dx_kvietimai += 1
        return tiklsoF_2dx(x)
    
    x = (copy.deepcopy(intervalas[0]) + copy.deepcopy(intervalas[1])) / 2
    iteracijos = 0
    tarpiniai_t = []
    while True:
        iteracijos += 1
        tarpiniai_t.append((x, f(x)))
        fprime = f_dx(x)
        fdouble = f_2dx(x)
        
        if abs(fprime) < epsilon:
            break
        
        x = x - fprime / fdouble
    
    print()
    print(f"Niutono metodo algoritmas rado atsakymą per {iteracijos} iteracijas.")
    print(f"su ε = {epsilon}, ats.: {tarpiniai_t[-1][1]}")
    print(f"tiklsoF() buvo iškviestas {f_kvietimai} kartų, tiklsoF_dx() {fdx_kvietimai} kartų, tiklsoF_2dx() {f2dx_kvietimai} kartų.")
    
    for pt in tarpiniai_t:
        pl.plot(pt[0], pt[1], 'xy')
    pl.plot(x, tarpiniai_t[-1][1], 'ob')

x = np.arange(intervalas[0], intervalas[1], 0.1)

pl.plot(x, tiklsoF(x), 'k-')
dalijimoPusiau()
pl.show()

pl.plot(x, tiklsoF(x), 'k-')
auksinioPjuvio()
pl.show()

pl.plot(x, tiklsoF(x), 'k-')
niutonoMetodas()
pl.show()