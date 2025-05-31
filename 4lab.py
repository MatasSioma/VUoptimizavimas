import numpy as np
from simpleksas import simplex
#2314004 - stud knygeles nr

A = np.array([
        [-1, 1, -1, -1],
        [2, 4, 0, 0],
        [0, 0, 1, 1]
    ])
c = np.array([2, -3, 0, -5])

b = np.array([8, 10, 3])
val_pr, point_pr, base_pr = simplex(A, b, c)
print("Pradinis uždavinys:")
print("Optimalus sprendinys:", point_pr)
print("Minimali tikslo funkcijos reikšmė:", val_pr)
print(f'Bazė:', " ".join(map(str, base_pr)))
print()

b = np.array([0, 0, 4])
val_ind, point_ind, base_ind = simplex(A, b, c)
print("Individualus uždavinys (...004):")
print("Optimalus sprendinys:", point_ind)
print("Minimali tikslo funkcijos reikšmė:", val_ind)
print(f'Bazė:', " ".join(map(str, base_ind)))

print()

print("Palyginimas:")
print(f"Pradinis:   min Z = {val_pr:.2f}, x = {point_pr}")
print(f"Individualus: min Z = {val_ind:.2f}, x = {point_ind}")