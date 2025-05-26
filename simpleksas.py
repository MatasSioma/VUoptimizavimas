import numpy as np

def simplex(a, b, c):
    # Sukūriame list'ą su bazinėmis reikšmėmis
    baze = np.arange(len(a)) + len(c)

    # Patikriname ar stulpelių kiekis a ir c kintamuosiuose sutampa su kintamųjų skaičiumi
    if len(a) + len(c) != len(a[0]):
        B = np.identity(len(a))
        a = np.hstack((a, B))

    # Sukūriame pradinę simplekso lentelę įrašant masyvus vertikaliai
    lentele = np.vstack((
        np.array([None] + [0] + list(c) + [0] * len(a)), # Pirmoji eilutė yra tikslinės funkcijos koeficientai
        np.hstack((np.transpose([baze]), np.transpose([b]), a)) # Likusios eilės yra apribojimai ir kintamieji
    ))


    iteracijos = 0
    while not np.all(lentele[0, 2:] >= 0):
        tikslo_funkcijos_koeficientai = lentele[0, 2:]
        ieinancio_kintamojo_indeksas = np.argmin(tikslo_funkcijos_koeficientai) + 2

        loginis_masyvas = lentele[:, ieinancio_kintamojo_indeksas] > 0
        pasirinktos_eilutes = lentele[loginis_masyvas]

        minimalus_indeksas = np.argmin(pasirinktos_eilutes[:, 1] / pasirinktos_eilutes[:, ieinancio_kintamojo_indeksas])
        iseinanti_eilute = np.nonzero(loginis_masyvas)[0][minimalus_indeksas]


        atraminis = lentele[iseinanti_eilute, ieinancio_kintamojo_indeksas]
        lentele[iseinanti_eilute, 1:] /= atraminis

        neiseinancios_eilutes = np.arange(len(lentele)) != iseinanti_eilute

        atraminis_faktorius = lentele[neiseinancios_eilutes, ieinancio_kintamojo_indeksas] / lentele[iseinanti_eilute, ieinancio_kintamojo_indeksas] / lentele[iseinanti_eilute, ieinancio_kintamojo_indeksas]
        atraminis_faktorius = atraminis_faktorius[:, np.newaxis]

        lentele[neiseinancios_eilutes, 1:] -= atraminis_faktorius * lentele[iseinanti_eilute, 1:]

        lentele[iseinanti_eilute, 0] = ieinancio_kintamojo_indeksas - 2

        iteracijos += 1

    optimalus_sprendinys = np.zeros(len(lentele[0, 2:]))
    indeksai = np.arange(1, len(lentele))

    geri_indeksai = np.logical_and(lentele[indeksai, 0] < len(lentele[0, 2:]), lentele[indeksai, 0] != None)
    sprendimo_kintamuju_indeksai = lentele[indeksai[geri_indeksai], 0].astype(int)
    optimalus_sprendinys[sprendimo_kintamuju_indeksai] = lentele[indeksai[geri_indeksai], 1]

    optimali_reiksme = -1 * lentele[0, 1]

    bazes_indeksai = np.where(optimalus_sprendinys != 0)[0] + 1

    return optimali_reiksme, optimalus_sprendinys, bazes_indeksai